module particle_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use coord_class,       only : coordList
  use RNG_class,         only : RNG

  implicit none
  private

  !!
  !! Particle types paramethers
  !!
  integer(shortInt), parameter,public :: P_NEUTRON   = 1,&
                                         P_PHOTON    = 2,&
                                         P_PRECURSOR = 3

  !!
  !! Public particle type procedures
  !!
  public :: verifyType
  public :: printType

  !!
  !! Particle compressed for storage
  !!
  !! Public Members:
  !!   wgt      -> Weight of the particle
  !!   r        -> Global Position of the particle [cm]
  !!   dir      -> Direction vector of the particle (normalised to 1.0)
  !!   E        -> Energy of the particle [MeV]
  !!   G        -> Energy Group of the particle
  !!   isMG     -> True if particle uses MG data
  !!   type     -> Physical Type of the particle (NEUTRON, PHOTON etc.)
  !!   time     -> Position in time of the particle [s]
  !!   lambda_i -> Decay constants for precursor groups [/s]
  !!   fd_i     -> Fraction of precursors in each group (usually beta_i/beta)
  !!   E_out_i  -> Sampled energy out for each precursor's delayed neutron [MeV]
  !!   matIdx   -> material Index in which particle is present
  !!   cellIdx  -> Cell Index at the lowest level in which particle is present
  !!   uniqueID -> Unique ID of the cell at the lowest level in which particle is present
  !!
  !! Interface:
  !!   assignemnt(=)  -> Build particleState from particle
  !!   operator(.eq.) -> Return True if particle are exactly the same
  !!   display        -> Print debug information about the state to the console
  !!
  type, public :: particleState
    real(defReal)              :: wgt  = ZERO       ! Particle weight
    real(defReal),dimension(3) :: r    = ZERO       ! Global position
    real(defReal),dimension(3) :: dir  = ZERO       ! Global direction
    real(defReal)              :: E    = ZERO       ! Energy
    integer(shortInt)          :: G    = 0          ! Energy group
    logical(defBool)           :: isMG = .false.    ! Is neutron multi-group
    integer(shortInt)          :: type = P_NEUTRON  ! Particle physical type
    real(defReal)              :: time = ZERO       ! Particle time position
    real(defReal),dimension(precursorGroups):: lambda_i = -ONE   ! Precursor decay constants
    real(defReal),dimension(precursorGroups):: fd_i     = -ONE   ! Precursor group fractions
    real(defReal),dimension(precursorGroups):: E_out_i  = -ONE   ! Delayed neutron energies
    integer(shortInt)          :: matIdx   = -1     ! Material index where particle is
    integer(shortInt)          :: cellIdx  = -1     ! Cell idx at the lowest coord level
    integer(shortInt)          :: uniqueID = -1     ! Unique id at the lowest coord level
  contains
    generic    :: assignment(=)  => fromParticle
    generic    :: operator(.eq.) => equal_particleState
    procedure  :: display        => display_particleState
    procedure  :: fromParticle   => particleState_fromParticle
    procedure  :: kill           => kill_particleState

    ! Private procedures
    procedure,private :: equal_particleState
  end type particleState

!  !!
!  !! Archived state of the particle used for tallying transitions, fission matrixes etc.
!  !!
!  type, public,extends(phaseCoord) ::
!  contains
!    generic    :: operator(.eq.) => equal_particleState
!    procedure :: fromParticle    => particleState_fromParticle
!
!    ! Private procedures
!    procedure,private :: equal_particleState
!  end type

  !!
  !! This type represents particle
  !!
  !! In current form it was designed to support neutron and other neutral particles
  !! By extension(inheritance) support for photons or charged particles could be introduced
  !!
  !! Now also represents precursors.
  !!
  !! In addition to representing a particle by physical parameters, it containes additional
  !! data like RNG pointer. This is to enable access to objects associated with a particlular
  !! particle history from anywhere where a particle was passed. This it works kind of like a
  !! global scope for a specific particle. Such setup should enable particle to control how it
  !! is beeing processed during Monte Carlo transport.
  !!
  !! Particle also contains snapshots (states) at privious points in its history. This is used
  !! for tallying. It is necessary to compare current state with the earlier state to score e.g.
  !! collision probabilities.
  !!
  type, public :: particle
    ! Particle phase space data
    type(coordList)            :: coords
    real(defReal)              :: E         ! Particle Energy
    integer(shortInt)          :: G         ! Particle Energy Group
    real(defReal)              :: w         ! Particle Weight
    real(defReal)              :: time      ! Particle time point

    ! Precursor particle data
    real(defReal), dimension(precursorGroups) :: lambda_i  ! Precursor decay constants
    real(defReal), dimension(precursorGroups) :: fd_i      ! Precursor fractions
    real(defReal), dimension(precursorGroups) :: E_out_i   ! Delayed neutron energies

    ! Particle flags
    real(defReal)              :: w0             ! Particle initial weight (for implicit, variance reduction...)
    logical(defBool)           :: isDead
    logical(defBool)           :: isMG
    real(defReal)              :: timeMax = ZERO ! Maximum neutron time before cut-off
    integer(shortInt)          :: fate = 0       ! Neutron's fate after being subjected to an operator
    integer(shortInt)          :: type           ! Particle type

    ! Particle processing information
    class(RNG), pointer        :: pRNG  => null()  ! Pointer to RNG associated with the particle
    real(defReal)              :: k_eff            ! Value of default keff for implicit source generation
    integer(shortInt)          :: geomIdx          ! Index of the geometry used by the particle

    ! Archived snapshots of previous states
    type(particleState)        :: preHistory
    type(particleState)        :: preTransition
    type(particleState)        :: prePath
    type(particleState)        :: preCollision

  contains
     ! Build procedures
    generic              :: build => buildCE, buildMG
    generic              :: assignment(=) => particle_fromParticleState

    ! Inquiry about coordinates
    procedure                  :: rLocal
    procedure                  :: rGlobal
    procedure                  :: dirLocal
    procedure                  :: dirGlobal
    procedure                  :: nesting
    procedure                  :: getCellIdx
    procedure                  :: getUniIdx
    procedure                  :: matIdx
    procedure, non_overridable :: getType
    procedure                  :: getSpeed
    
    ! Inquiry about precursor parameters
    procedure               :: getPrecWeight
    procedure               :: getExpPrecWeight
    procedure               :: getTimedWeight
    procedure               :: getDelayedEnergy

    ! Operations on coordinates
    procedure            :: moveGlobal
    procedure            :: moveLocal
    procedure            :: rotate
    procedure            :: teleport
    procedure            :: point
    procedure            :: takeAboveGeom
    procedure            :: setMatIdx

    ! Save particle state information
    procedure, non_overridable  :: savePreHistory
    procedure, non_overridable  :: savePreTransition
    procedure, non_overridable  :: savePrePath
    procedure, non_overridable  :: savePreCollision

    ! Debug procedures
    procedure            :: display => display_particle
    procedure            :: typeToChar

    !! Private - Implementation specific procedures
    procedure,private                   :: buildCE
    procedure,private                   :: buildMG
    procedure,non_overridable,private   :: particle_fromParticleState

  end type particle

contains

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle build and assignment procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise CE particle
  !! Necessary arguments:
  !!   r   -> Global Position
  !!   dir -> Global direction
  !!   E   -> Energy [MeV]
  !!   w   -> particle weight
  !! Optional arguments:
  !!   t        -> particle time (default = 0.0)
  !!   lambda_i -> precursor decay constant (default = -1.0)
  !!   fd_i     -> precursor fraction (default = -1.0)
  !!   type     -> particle type (default = P_NEUTRON)
  !!
  pure subroutine buildCE(self, r, dir, E, w, t, type, lambda_i, fd_i, E_out_i)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r
    real(defReal),dimension(3),intent(in)   :: dir
    real(defReal),intent(in)                :: E
    real(defReal),intent(in)                :: w
    real(defReal),optional,intent(in)       :: t
    integer(shortInt),intent(in),optional   :: type
    real(defReal),dimension(precursorGroups),optional,intent(in) :: lambda_i
    real(defReal),dimension(precursorGroups),optional,intent(in) :: fd_i
    real(defReal),dimension(precursorGroups),optional,intent(in) :: E_out_i


    call self % coords % init(r, dir)
    self % E  = E
    self % w  = w
    self % w0 = w

    self % isDead = .false.
    self % isMG   = .false.

    if(present(t)) then
      self % time = t
    else
      self % time = ZERO
    end if
    
    if(present(lambda_i)) then
      self % lambda_i = lambda_i
    else
      self % lambda_i = -ONE
    end if
    
    if(present(fd_i)) then
      self % fd_i = fd_i
    else
      self % fd_i = -ONE
    end if
    
    if(present(E_out_i)) then
      self % E_out_i = E_out_i
    else
      self % E_out_i = -ONE
    end if

    if(present(type)) then
      self % type = type
    else
      self % type = P_NEUTRON
    end if

  end subroutine buildCE

  !!
  !! Initialise MG particle
  !! Necessary arguments:
  !!   r   -> Global Position
  !!   dir -> Global direction
  !!   G   -> Energy Group
  !!   w   -> particle weight
  !! Optional arguments:
  !!   t   -> particle time (default = 0.0)
  !!   lambda_i -> precursor decay constant (default = -1.0)
  !!   fd_i     -> precursor fraction (default = -1.0)
  !!   type-> particle type (default = P_NEUTRON)
  !!
  subroutine buildMG(self, r, dir, G, w, t, type, lambda_i, fd_i)
    class(particle), intent(inout)          :: self
    real(defReal),dimension(3),intent(in)   :: r
    real(defReal),dimension(3),intent(in)   :: dir
    real(defReal),intent(in)                :: w
    integer(shortInt),intent(in)            :: G
    real(defReal),intent(in),optional       :: t
    integer(shortInt),intent(in),optional   :: type
    real(defReal),dimension(precursorGroups),optional,intent(in) :: lambda_i
    real(defReal),dimension(precursorGroups),optional,intent(in) :: fd_i

    call self % coords % init(r, dir)
    self % G  = G
    self % w  = w
    self % w0 = w

    self % isDead = .false.
    self % isMG   = .true.

    if(present(t)) then
      self % time = t
    else
      self % time = ZERO
    end if
    
    if(present(lambda_i)) then
      self % lambda_i = lambda_i
    else
      self % lambda_i = -ONE
    end if
    
    if(present(fd_i)) then
      self % fd_i = fd_i
    else
      self % fd_i = -ONE
    end if

    if(present(type)) then
      self % type = type
    else
      self % type = P_NEUTRON
    end if

  end subroutine buildMG

  !!
  !! Copy phase coordinates into particle
  !!
  subroutine particle_fromParticleState(LHS,RHS)
    class(particle), intent(inout)   :: LHS
    type(particleState), intent(in)  :: RHS

    LHS % w                     = RHS % wgt
    LHS % w0                    = RHS % wgt
    call LHS % takeAboveGeom()
    LHS % coords % lvl(1) % r   = RHS % r
    LHS % coords % lvl(1) % dir = RHS % dir
    LHS % E                     = RHS % E
    LHS % G                     = RHS % G
    LHS % isMG                  = RHS % isMG
    LHS % type                  = RHS % type
    LHS % time                  = RHS % time
    LHS % lambda_i              = RHS % lambda_i
    LHS % fd_i                  = RHS % fd_i
    LHS % E_out_i               = RHS % E_out_i

  end subroutine particle_fromParticleState

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle coordinates inquiry procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Return the position either at the deepest nested level or a specified level
  !!
  function rLocal(self,n)result(r)
    class(particle), intent(in)             :: self
    integer(shortInt), intent(in), optional :: n
    real(defReal), dimension(3)             :: r
    integer(shortInt)                       :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    r = self % coords % lvl(n_loc) % r

  end function rLocal

  !!
  !! Return the position at the highest level
  !!
  pure function rGlobal(self)result(r)
    class(particle), intent(in) :: self
    real(defReal), dimension(3) :: r

    r = self % coords % lvl(1) % r

  end function rGlobal

  !!
  !! Return the direction either at the deepest nested level or at a specified level
  !!
  function dirLocal(self,n)result(dir)
    class(particle), intent(in) :: self
    integer(shortInt), optional :: n
    real(defReal), dimension(3) :: dir
    integer(shortInt)           :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    dir = self % coords % lvl(n_loc) % dir

  end function dirLocal

  !!
  !! Return the direction at the highest nesting level
  !!
  function dirGlobal(self)result(dir)
    class(particle), intent(in) :: self
    real(defReal), dimension(3) :: dir

    dir = self % coords % lvl(1) % dir

  end function dirGlobal

  !!
  !! Return the lowest nesting level of the particle
  !!
  function nesting(self) result(n)
    class(particle), intent(in) :: self
    integer(shortInt)           :: n

    n = self % coords % nesting

  end function nesting

  !!
  !! Return cell index at a given nesting level n
  !!  If no n is given return for lowest nesting level
  !!
  function getCellIdx(self,n) result(idx)
    class(particle), intent(in)             :: self
    integer(shortInt),optional, intent(in)  :: n
    integer(shortInt)                       :: idx
    integer(shortInt)                       :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    idx = self % coords % lvl(n_loc) % cellIdx

  end function getCellIdx

  !!
  !! Return universe index at a given nesting level n
  !!
  function getUniIdx(self,n) result(idx)
    class(particle), intent(in)             :: self
    integer(shortInt),optional, intent(in)  :: n
    integer(shortInt)                       :: idx
    integer(shortInt)                       :: n_loc

    if(present(n)) then
      n_loc = n
    else
      n_loc = self % coords % nesting
    end if

    idx = self % coords % lvl(n_loc) % uniIdx

  end function getUniIdx

  !!
  !! Return current material index
  !!
  pure function matIdx(self) result(Idx)
    class(particle), intent(in) :: self
    integer(shortInt)           :: idx

    idx = self % coords % matIdx

  end function matIdx

  !!
  !! Return one of the particle Tpes defined in universal variables
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   P_NEUTRON_CE, P_NEUTRON_MG
  !!
  !! Errors:
  !!   None
  !!
  pure function getType(self) result(type)
    class(particle), intent(in) :: self
    integer(shortInt)           :: type

    if (self % isMG) then
      type = P_NEUTRON_MG
    else
      type = P_NEUTRON_CE
    end if

  end function getType

  !!
  !! Return speed of the particle
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Speed of light for photon, speed from E for neutron
  !!   Made up number for MG (corresponds to 0.025 eV)
  !!
  !! Errors:
  !!   
  !!
  pure function getSpeed(self) result(speed)
    class(particle), intent(in) :: self
    real(defReal)               :: speed

    if (self % type == P_PHOTON) then
      speed = lightSpeed
    else
      if (self % isMG) then
        speed = 0.7294888_defReal
      else
        speed = sqrt(self % E * TWO/neutronMass) * lightSpeed
      end if
    end if

  end function getSpeed
  
  !!
  !! Return neutron weight for forced decay at time t
  !!
  !! Args:
  !!   t       -> time of forced decay
  !!   delta_T -> time step in which forced decay occurs
  !!
  !! Result:
  !!   Neutron weight for decay at time t
  !!   Returns ZERO if particle is not a precursor
  !!
  !! Errors:
  !!   
  !!
  function getPrecWeight(self, t, delta_T) result(w_d)
    class(particle), intent(in) :: self
    real(defReal), intent(in)   :: t
    real(defReal), intent(in)   :: delta_T
    integer(shortInt)           :: i
    real(defReal)               :: w_d
    
    w_d = ZERO
    if (self % type == P_PRECURSOR) then
        ! Loop over precursor groups
        do i=1, precursorGroups
            w_d = w_d + self % fd_i(i) * self % lambda_i(i) *&
                    exp(-self % lambda_i(i) * (t - self % time))
        end do
    end if
    w_d = w_d * self % w * delta_T
  end function getPrecWeight
  
  !!
  !! Return expected neutron weight for decay in next time step
  !!
  !! Args:
  !!   t1     -> start time of next time step
  !!   step_T -> length of next time step
  !!
  !! Result:
  !!   Expected neutron weight for decay in next time step
  !!   Returns ZERO if particle is not a precursor
  !!
  !! Errors:
  !!   
  !!
  function getExpPrecWeight(self, t1, step_T) result(w_d_av)
    class(particle), intent(in) :: self
    real(defReal), intent(in)   :: t1, step_T
    integer(shortInt)           :: i
    real(defReal)               :: w_d_av
    
    w_d_av = ZERO
    if (self % type == P_PRECURSOR) then
        ! Loop over precursor groups
        do i=1, precursorGroups
            w_d_av = w_d_av + self % fd_i(i) * &
                    (exp(-self % lambda_i(i) * (t1 - self % time)) - &
                    exp(-self % lambda_i(i) * (t1 + step_T - self % time)))
        end do
    end if
    w_d_av = w_d_av * self % w

  end function getExpPrecWeight
  
  !!
  !! Return current timed weight of particle
  !!
  !! Args:
  !!   t      -> current time [s]
  !!
  !! Result:
  !!   Current timed weight at time t
  !!   Returns ZERO if particle is not a precursor
  !!
  !! Errors:
  !!   
  !!
  function getTimedWeight(self, t) result(w_Timed)
    class(particle), intent(in) :: self
    real(defReal), intent(in)   :: t
    integer(shortInt)           :: i
    real(defReal)               :: w_Timed
    
    w_Timed = ZERO
    if (self % type == P_PRECURSOR) then
        ! Loop over precursor groups
        do i=1, precursorGroups
            w_Timed = w_Timed + self % fd_i(i) * exp(-self % lambda_i(i) * (t - self % time))
        end do 
        !print *, 'lambda_i:       ', numToChar(self % lambda_i)
        !print *
    end if
    ! Multiply by original weight
    w_Timed = w_Timed * self % w
    !print *, 'Self time:      ', numToChar(self % time)
    !print *, 'Regular weight: ', numToChar(self % w)
    !print *, 'Current time:   ', numToChar(t)
    !print *, 'Timed weight:   ', numToChar(w_Timed)
  end function getTimedWeight
  
  
  !!
  !! Sample and return precursor energy at t
  !!
  !! Args:
  !!   t      -> current time [s]
  !!
  !! Result:
  !!   Samples from array of energies (see sampleDelayed in fissionCE)
  !!   Should automatically select the only energy if particle represents
  !!   only 1 precursor group
  !!
  !! Errors:
  !!   
  !!
  function getDelayedEnergy(self, t) result(E_out)
    class(particle), intent(in)               :: self
    real(defReal), intent(in)                 :: t
    real(defReal)                             :: r1, w_Prob
    real(defReal), dimension(precursorGroups) :: probArray
    integer(shortInt)                         :: i
    real(defReal)                             :: E_out
    
    if (self % type == P_PRECURSOR) then
        r1 = self % pRNG % get()
        
        w_Prob = ZERO
        
        ! Create w_Prob by summing over all groups
        !print *, "E_out_i"
        do i=1, precursorGroups
            !print *, numToChar(self % E_out_i(i))
            w_Prob = w_Prob + self % fd_i(i) * self % lambda_i(i) *&
                                    exp(-self % lambda_i(i) * (t - self % time))
        end do
        
        ! Create array of probabilities of selecting each group
        do i=1, precursorGroups
            probArray(i) = self % fd_i(i) * self % lambda_i(i) *&
                                    exp(-self % lambda_i(i) * (t - self % time)) / w_Prob
        end do
        
        ! Select a precursor energy
        do i=1, precursorGroups
            r1 = r1 - probArray(i)
            if(r1 < ZERO) then
                E_out = self % E_out_i(i)
                return
            end if
        end do
        
        ! Sampling fails
        print *, 'Precursor energy sampling failed in particle object'
        E_out = self % E_out_i(precursorGroups)
    else
        E_out = ZERO
    end if
  end function getDelayedEnergy

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle operations on coordinates procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Move the particle above the geometry
  !! NOTE: regionID & matIdx will be reset!!!
  !!
  subroutine moveGlobal(self,distance)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: distance

    call self % coords % moveGlobal(distance)

  end subroutine moveGlobal

  !!
  !! Move particle in local co-ordinates down to nesting level n
  !!
  subroutine moveLocal(self,distance,n)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: distance
    integer(shortInt), intent(in)  :: n

    call self % coords % moveLocal(distance,n)

  end subroutine moveLocal

  !!
  !! Rotate particle
  !!  mu  -> cosine of deflection from current direction
  !!  phi -> azimuthal angle of rotation
  !!
  subroutine rotate(self,mu,phi)
    class(particle), intent(inout) :: self
    real(defReal), intent(in)      :: mu
    real(defReal), intent(in)      :: phi

    call self % coords % rotate(mu,phi)

  end subroutine rotate

  !!
  !! Place particle at an arbitrary point in above the geometry
  !!
  subroutine teleport(self, r)
    class(particle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: r

    call self % coords % assignPosition(r)

  end subroutine teleport

  !!
  !! Point particle in direction dir in highest nesting level
  !! Propagates new direction to lower levels
  !!
  subroutine point(self, dir)
    class(particle), intent(inout)          :: self
    real(defReal), dimension(3), intent(in) :: dir

    call self % coords % assignDirection(dir)

  end subroutine point

  !!
  !! Resets the particle's nesting level
  !!
  subroutine takeAboveGeom(self)
    class(particle), intent(inout) :: self

    call self % coords % takeAboveGeom()

  end subroutine takeAboveGeom

  !!
  !! Set Material index for testing purposes
  !!
  pure subroutine setMatIdx(self,matIdx)
    class(particle), intent(inout) :: self
    integer(shortInt), intent(in)  :: matIdx

    self % coords % matIdx = matIdx

  end subroutine setMatIdx

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle save state procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePreHistory(self)
    class(particle), intent(inout) :: self

    self % preHistory = self

  end subroutine savePreHistory

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePreTransition(self)
    class(particle), intent(inout) :: self

    self % preTransition = self

  end subroutine savePreTransition

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePrePath(self)
    class(particle), intent(inout) :: self

    self % prePath = self

  end subroutine savePrePath

  !!
  !! Save state of the particle at the beginning of history
  !!
  subroutine savePreCollision(self)
    class(particle), intent(inout) :: self

    self % preCollision = self

  end subroutine savePreCollision

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle debug procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Display state of a particle
  !!
  subroutine display_particle(self)
    class(particle), intent(in) :: self
    type(particleState)         :: state

    state = self
    call state % display()
    print *, self % coords % matIdx

  end subroutine display_particle

  !!
  !! Return character that describes the type of particle
  !!
  function typeToChar(self) result(c)
    class(particle), intent(in) :: self
    character(:), allocatable   :: c
    character(2)                :: eType

    if( self % isMG) then
      eType = 'MG'
    else
      eType = 'CE'
    end if

    c = eType // ' ' // trim(printType( self % type))

  end function typeToChar


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Particle state and phaseCoord procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Copy particle into phase coordinates
  !!
  subroutine particleState_fromParticle(LHS,RHS)
    class(particleState), intent(out)  :: LHS
    class(particle), intent(in)        :: RHS

    LHS % wgt  = RHS % w
    LHS % r    = RHS % rGlobal()
    LHS % dir  = RHS % dirGlobal()
    LHS % E    = RHS % E
    LHS % G    = RHS % G
    LHS % isMG = RHS % isMG
    LHS % type = RHS % type
    LHS % time = RHS % time

    ! Save all indexes
    LHS % matIdx   = RHS % coords % matIdx
    LHS % uniqueID = RHS % coords % uniqueId
    LHS % cellIdx  = RHS % coords % cell()

  end subroutine particleState_fromParticle

  !!
  !! Define equal operation on phase coordinates
  !!  Phase coords are equal if all their components are the same
  !!
  function equal_particleState(LHS,RHS) result(isEqual)
    class(particleState), intent(in) :: LHS
    class(particleState), intent(in) :: RHS
    logical(defBool)              :: isEqual

    isEqual = .true.
    isEqual = isEqual .and. LHS % wgt == RHS % wgt
    isEqual = isEqual .and. all(LHS % r   == RHS % r)
    isEqual = isEqual .and. all(LHS % dir == RHS % dir)
    isEqual = isEqual .and. LHS % time == RHS % time
    isEqual = isEqual .and. LHS % isMG .eqv. RHS % isMG
    isEqual = isEqual .and. LHS % type == RHS % type
    isEqual = isEqual .and. LHS % matIdx   == RHS % matIdx
    isEqual = isEqual .and. LHS % cellIdx  == RHS % cellIdx
    isEqual = isEqual .and. LHS % uniqueID == RHS % uniqueID

    if( LHS % isMG ) then
      isEqual = isEqual .and. LHS % G == RHS % G
    else
      isEqual = isEqual .and. LHS % E == RHS % E
    end if
  end function equal_particleState

!  !!
!  !! Copy particle state into archive object
!  !!
!  subroutine particleState_fromParticle(LHS,RHS)
!    class(particleState), intent(out) :: LHS
!    class(particle), intent(in)       :: RHS
!
!    ! Call superclass procedure
!    call phaseCoord_fromParticle(LHS,RHS)
!
!  end subroutine particleState_fromParticle

!  !!
!  !! Extend equality definition to particle state
!  !!
!  function equal_particleState(LHS,RHS) result(isEqual)
!    class(particleState), intent(in) :: LHS
!    type(particleState), intent(in)  :: RHS
!    logical(defBool)                 :: isEqual
!
!    ! Call superclass procedure
!    isEqual = equal_phaseCoord(LHS,RHS)
!
!  end function equal_particleState

  !!
  !! Prints state of the phaseCoord
  !!
  subroutine display_particleState(self)
    class(particleState), intent(in) :: self

    print *, self % r, self % dir, self % E, self % G, self % isMG, self % wgt, self % time

  end subroutine display_particleState

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill_particleState(self)
    class(particleState), intent(inout) :: self

    self % wgt  = ZERO
    self % r    = ZERO
    self % dir  = ZERO
    self % E    = ZERO
    self % G    = 0
    self % isMG = .false.
    self % type = P_NEUTRON
    self % lambda_i = -ONE
    self % fd_i     = -ONE
    self % time = ZERO
    self % matIdx   = -1
    self % cellIdx  = -1
    self % uniqueID = -1

  end subroutine kill_particleState


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Misc Procedures
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Returns .true. if the integer is valid particle type code
  !! Returns .false. otherwise
  !!
  elemental function verifyType(type) result(isValid)
    integer(shortInt), intent(in) :: type
    logical(defBool)              :: isValid

    isValid = .false.

    ! Check against particles types
    isValid = isValid .or. type == P_NEUTRON
    isValid = isValid .or. type == P_PHOTON

  end function verifyType

  !!
  !! Returns character with a description of particle type
  !!
  pure function printType(type) result(name)
    integer(shortInt), intent(in) :: type
    character(:),allocatable      :: name

    select case(type)
      case(P_NEUTRON)
        name = 'Neutron'

      case(P_PHOTON)
        name = 'Photon'

      case default
        name = 'INVALID PARTICLE TYPE'

    end select
  end function printType

end module particle_class
