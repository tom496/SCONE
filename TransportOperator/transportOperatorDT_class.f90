!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorDT
  contains
    procedure :: transit => deltaTracking
  end type transportOperatorDT

contains

  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorDT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance
    character(100), parameter :: Here = 'deltaTracking (transportOIperatorDT_class.f90)'

    ! Get majornat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getMajorantXS(p)

    DTLoop:do
      distance = -log( p% pRNG % get() ) * majorant_inv

      if (p % timeMax > ZERO .and. p % time + distance / p % getSpeed() > p % timeMax) then
        distance = distance * (p % timeMax - p % time)/(distance / p % getSpeed())
        p % fate = AGED_FATE
        !print *, '  -Particle time before aged fate:', numToChar(p%time)
        p % time = p % timeMax
        ! Move partice in the geometry
        call self % geom % teleport(p % coords, distance)
        return
      endif

      ! Move partice in the geometry
      call self % geom % teleport(p % coords, distance)

      ! Update time
      p % time = p % time + distance / p % getSpeed()

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

      ! Check for void
      if( p % matIdx() == VOID_MAT) cycle DTLoop

      ! Obtain the local cross-section
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Protect Against Sillines
      !if( sigmaT*majorant_inv < ZERO .or. ONE < sigmaT*majorant_inv) then
      !  call fatalError(Here, "TotalXS/MajorantXS is silly: "//numToChar(sigmaT*majorant_inv))
      !end if

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (p % pRNG % get() < sigmaT*majorant_inv) exit DTLoop

    end do DTLoop

    call tally % reportTrans(p)
  end subroutine deltaTracking


end module transportOperatorDT_class
