module isotropicMG_class

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError, linFind, searchError
  use dictionary_class,      only : dictionary
  use IOdictionary_class,    only : IOdictionary
  use RNG_class,             only : RNG

  use perMaterialMgXs_inter, only : perMaterialMgXs
  use outscatterCDF_class,   only : outscatterCDF
  use xsMacroSet_class,      only : xsMacroSet
  use releaseMatrixMG_class, only : releaseMatrixMG

  implicit none
  private

  !!
  !! Small type to store decription of the material
  !!
  type, private :: materialData
    character(nameLen) :: matName
    logical(defBool)   :: isFissile = .false.
    logical(defBool)   :: isActive  = .true.
  end type materialData

  type, public, extends(perMaterialMgXs) :: isotropicMG
    private
    integer(shortInt)   :: nG
    integer(shortInt)   :: nMat

    type(xsMacroSet),dimension(:,:),pointer            :: XSs => null()  ! (energyGroup,matIdx)
    type(outscatterCDF), dimension(:,:), allocatable   :: transferMatrix ! (energyGroup,matIdx)
    type(outscatterCDF), dimension(:), allocatable     :: chiValues      ! (matIdx)
    type(releaseMatrixMG), dimension(:), allocatable   :: releaseData    ! (matIdx)
    type(materialData), dimension(:), allocatable      :: matData        ! (matIdx)
    real(defReal), dimension(:), allocatable           :: majorantXS     ! (energyGroup)

  contains

    !* INTERFACE PROCEDURES *!
    ! Initialisation procdures
    procedure :: init

    ! Procedures to obtain XSs
    procedure  :: getMatMacroXS
    procedure  :: getTransXS
    procedure  :: getMajorantXS
    procedure  :: getTotalMatXS

    ! Procedures to obtain emission data
    procedure  :: releaseAt
    procedure  :: sampleMuGout

    ! Procedures to access material information
    procedure  :: getMatIdx
    procedure  :: getMatName
    procedure  :: isFissileMat

    !* TYPE PROCEDURES *!
    procedure, private :: readMaterial
    procedure, private :: calculateMajorant

  end type isotropicMG

contains

  !!
  !! Load material data from a dictionary
  !!
  subroutine init(self,dict)
    class(isotropicMG), intent(inout) :: self
    type(dictionary), intent(in)      :: dict
    integer(shortInt)                 :: nG,nMat,matIdx

    ! Read number of energy groups and materials
    nG = dict % getInt('numberOfGroups')
    nMat = size( dict % keysDict() )

    self % nG   = nG
    self % nMat = nMat

    ! Allocate space fo XS data
    if (associated( self % XSs           )) deallocate (self % XSs            )
    if (allocated(  self % transferMatrix)) deallocate (self % transferMatrix )
    if (allocated(  self % chiValues     )) deallocate (self % chiValues      )
    if (allocated(  self % releaseData   )) deallocate (self % releaseData    )
    if (allocated(  self % matData       )) deallocate (self % matData        )

    allocate( self % XSs(nG,nMat)           )
    allocate( self % transferMatrix(ng,nMat))
    allocate( self % chiValues(nMat)        )
    allocate( self % releaseData(nMat)      )
    allocate( self % matData(nMat)          )

    ! Read material names
    self % matData(:) % matName = dict % keysDict()

    ! Read individual material data
    do matIdx=1,nMat
      call self % readMaterial(dict % getDict( self % matData(matIdx) % matName ),&
                               matIdx, &
                               self % nG)
    end do

    call self % calculateMajorant()

  end subroutine init

  !!
  !! Read material under index "idx" from dictionary "dict"
  !! For group transfer matrixes indexing convention is (G_out,G_in)
  !! When Matrixes are rank1 then the sequence follows SERPENT(2.1.30) output convention:
  !! G_in->1, G_in->2 ...
  !!
  subroutine readMaterial(self,dict,idx, nG)
    class(isotropicMG), intent(inout) :: self
    type(dictionary), intent(in)      :: dict
    integer(shortInt), intent(in)     :: idx
    integer(shortInt), intent(in)     :: nG
    real(defReal),dimension(nG)       :: tempXS
    real(defReal),dimension(nG*nG)    :: tempXSmatrix_rank1
    real(defReal),dimension(nG,nG)    :: tempXSmatrix
    type(IOdictionary)                :: xsDict
    integer(shortInt)                 :: i
    logical(defBool)                  :: isFissile
    character(100), parameter         :: Here='readMaterial (isotropicMG_class.f90)'

    ! Obtain path from dict and read into xsDict
    call xsDict % initFrom( dict % getChar('xsFile'))

    ! Check if material is fissile
    isFissile = xsDict % isPresent('fission')
    self % matData(idx) % isFissile = isFissile

    ! Verify size of stored data
    if (size(xsDict % getRealArray('capture')) /= nG) then
      call fatalError(Here,'capture xs are inconsistant with number of energy groups')

    elseif (size(xsDict % getRealArray('scattering_P0')) /= nG*nG) then
      call fatalError(Here,'scatter xs are inconsistant with number of energy groups')

    else if (size(xsDict % getRealArray('scatteringMultiplicity')) /= nG*nG) then
      call fatalError(Here,'scattering production data is inconsistant with number of energy groups')

    end if

    if (isFissile) then
      if (size(xsDict % getRealArray('fission')) /= nG) then
        call fatalError(Here,'fission xs are inconsistant with number of energy groups')

      else if (size(xsDict % getRealArray('chi')) /= nG) then
        call fatalError(Here,'chi data is inconsistant with number of energy groups')

      else if (size(xsDict % getRealArray('nu')) /= nG) then
        call fatalError(Here,'nu data is inconsistant with number of energy groups')

      end if
    end if

    ! Load and store capture XSs
    tempXS = xsDict % getRealArray('capture')
    if (any( tempXS < 0.0)) call fatalError(Here,'capture xss are -ve')
    self % XSs(:,idx) % captureXS = tempXS

    ! Load and store fission XSs
    if (isFissile) then
      tempXS = xsDict % getRealArray('fission')
    else
      tempXS = 0.0
    end if
    if (any( tempXS < 0.0)) call fatalError(Here,'fission xss are -ve')
    self % XSs(:,idx) % fissionXS = tempXS

    ! Load and store chi values
    if (isFissile) then
      tempXS = xsDict % getRealArray('chi')
    else
      tempXS = 0.0
      tempXS(1) = 1.0 ! Avoid Floating point exception
    end if
    if (any( tempXS < 0.0)) call fatalError(Here,'chi is -ve')
    call self % chiValues(idx) % init(tempXS)

    ! Load scattering matrix. Indexing convenction is (G_out,G_in)
    tempXSmatrix_rank1 = xsDict % getRealArray('scattering_P0')
    if (any( tempXSmatrix_rank1 < 0.0)) call fatalError(Here,'Scattering_P0 is -ve')
    tempXSmatrix = reshape(tempXSmatrix_rank1,[nG, nG])
    do i=1,nG
      call self % transferMatrix(i,idx) % init ( tempXSmatrix(:,i) )
    end do
    self % XSs(:,idx) % scatterXS = sum(tempXSmatrix,1)

    ! Load production matrix and nu. Indexing convenction is (G_out,G_in).
    tempXSmatrix_rank1 = xsDict % getRealArray('scatteringMultiplicity')
    if (any( tempXSmatrix_rank1 < 0.0)) call fatalError(Here,'scateringMultiplicity in -ve')
    tempXSmatrix = reshape(tempXSmatrix_rank1,[nG, nG])

    if (isFissile) then
      tempXS = xsDict % getRealArray('nu')
    else
      tempXS = 0.0
    end if

    call self % releaseData(idx) % init(tempXS, tempXSmatrix)

    ! Calculate total XS
    self % XSs(:,idx) % totalXS = self % XSs(:,idx) % scatterXS + &
                                  self % XSs(:,idx) % captureXS + &
                                  self % XSs(:,idx) % fissionXS

  end subroutine readMaterial

  !!
  !! Recalculates majorant using current active materials
  !!
  subroutine calculateMajorant(self)
    class(isotropicMG), intent(inout)           :: self
    integer(shortInt), dimension(:),allocatable :: activeIdx
    logical(defBool), dimension(:), allocatable :: mask
    integer(shortInt)                           :: i

    ! Find indexes of active materials
    mask = self % matData % isActive
    activeIdx = pack([(i,i=1,self % nG)], mask)

    ! Calculate majorantXS
    self % majorantXS = maxval( self % XSs(:,activeIdx) % totalXS, 2 )

  end subroutine calculateMajorant

  !!
  !! Attach pointer to approperiate XS data
  !!
  subroutine getMatMacroXS(self,macroXS,G,matIdx)
    class(isotropicMG), intent(inout)            :: self
    type(xsMacroSet),pointer,intent(inout)       :: macroXS
    integer(shortInt),intent(in)                 :: G
    integer(shortInt),intent(in)                 :: matIdx

    macroXS => self % XSs(G,matIdx)
    ! The weirdest error
    ! It seems that same call to a procedure is necessary to avoid "undefined reference error"
    ! in collision operator
    call macroXS % dummy()

  end subroutine getMatMacroXS

  !!
  !! Return transport XS (in general diffrent from total XS)
  !!
  function getTransXS(self,G,matIdx) result(xs)
    class(isotropicMG), intent(inout)  :: self
    integer(shortInt), intent(in)      :: G
    integer(shortInt), intent(in)      :: matIdx
    real(defReal)                      :: xs

    xs = self % XSs(G,matIdx) % totalXS

  end function getTransXS

  !!
  !! Return majorant XS (in general should be largest of TRANSPORT XSs)
  !!
  function getMajorantXS(self,G) result (xs)
    class(isotropicMG), intent(inout)  :: self
    integer(shortInt), intent(in)      :: G
    real(defReal)                      :: xs

    xs = self % majorantXS(G)

  end function getMajorantXS

  !!
  !! Return total XS of material
  !!
  function getTotalMatXS(self,G,matIdx) result (xs)
    class(isotropicMG), intent(inout)  :: self
    integer(shortInt), intent(in)      :: G
    integer(shortInt), intent(in)      :: matIdx
    real(defReal)                      :: xs

    xs = self % XSs(G,matIdx) % totalXS

  end function getTotalMatXS

  !!
  !! Returns neutron release at G_in for material matIdx and reaction MT
  !!
  function releaseAt(self,G_in,G_out,MT,matIdx) result(nu)
    class(isotropicMG), intent(inout)   :: self
    integer(shortInt), intent(in)       :: G_in
    integer(shortInt), intent(in)       :: G_out
    integer(shortInt), intent(in)       :: MT
    integer(shortInt), intent(in)       :: matIdx
    real(defReal)                       :: nu
    character(100), parameter           :: Here = 'releaseAt (isotropicMG_class.f90)'

    select case(MT)
      case(macroAllScatter)
        nu = self % releaseData(matIdx) % scatterRelease(G_in,G_out)

      case(macroFission)
        nu = self % releaseData(matIdx) % fissionRelease(G_in)

      case default
        call fatalError(Here,'Unrecoginsed MT number')

    end select

  end function releaseAt

  !!
  !! Samples deflection angle in LAB frame and post-colission energy group
  !! This implementation is CRAP. SHOULD BE BRANCHLESS. IMPROVE IT! **##~~##**
  !!
  subroutine sampleMuGout(self,mu,G_out,G_in,rand,MT,matIdx)
    class(isotropicMG), intent(inout)  :: self
    real(defReal), intent(out)         :: mu
    integer(shortInt), intent(out)     :: G_out
    integer(shortInt), intent(in)      :: G_in
    class(RNG), intent(inout)          :: rand
    integer(shortInt), intent(in)      :: MT
    integer(shortInt), intent(in)      :: matIdx
    real(defReal)                      :: r1
    character(100), parameter          :: Here ='sampleMuGout (isotropicMG_class.f90)'

    mu = TWO * rand % get() - ONE

    r1 = rand % get()

    select case(MT)
      case(macroAllScatter)
        G_out = self % transferMatrix(G_in, matIdx) % invert(r1)

      case(macroFission)
        G_out = self % chiValues(matIdx) % invert(r1)

      case default
        call fatalError(Here,'Unrecoginsed MT number')

    end select

  end subroutine sampleMuGout


  !!
  !! Return matIdx of material with matName
  !!
  function getMatIdx(self,matName) result(matIdx)
    class(isotropicMG), intent(in)      :: self
    character(*), intent(in)            :: matName
    integer(shortInt)                   :: matIdx
    character(100), parameter           :: Here ='getMatIdx (isotropicMG_class.f90)'

    matIdx = linFind(self % matData % matName, matName)
    call searchError(matIdx,Here)

  end function getMatIdx

  !!
  !! Returns matName of material with matIdx
  !!
  function getMatName(self,matIdx) result(matName)
    class(isotropicMG), intent(in)      :: self
    integer(shortInt), intent(in)       :: matIdx
    character(nameLen)                  :: matName

    matName = self % matData(matIdx) % matName

  end function getMatName

  !!
  !! Returns .true. if material is fissile
  !!
  function isFissileMat(self,matIdx) result(isIt)
    class(isotropicMG), intent(in)  :: self
    integer(shortInt), intent(in)   :: matIdx
    logical(defBool)                :: isIt

    isIt = self % matData(matIdx) % isFissile

  end function isFissileMat


end module isotropicMG_class
