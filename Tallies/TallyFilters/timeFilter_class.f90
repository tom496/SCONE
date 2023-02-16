module timeFilter_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use particle_class,    only : particleState
  use dictionary_class,  only : dictionary
  use tallyFilter_inter, only : tallyFilter

  implicit none
  private

  !!
  !! Filter that tests if the time of a particle is in an interval.
  !!
  !! Private members:
  !!   Tmin -> minimum value of time
  !!   Tmax -> maximum value of time
  !!
  !! Interface:
  !!   tallyFilter Interface
  !!   build -> build filter from components
  !!
  !! NOTE: Time is not enforced to be +ve. Could be useful in debugging
  !!
  !! Sample Dictionary Input:
  !!
  !!   CEfilter {
  !!     type timeFilter;
  !!     Tmin 1.0;
  !!     Tmax 2.0;
  !!   }
  type, public,extends(tallyFilter) :: timeFilter
    private
    real(defReal)     :: Tmin
    real(defReal)     :: Tmax
  contains
    procedure :: init
    procedure :: isPass

    !! Instance specific procedures
    generic :: build => build_CE
    procedure :: build_CE

  end type timeFilter

contains

  !!
  !! Initialise timeFilter from dictionary
  !!
  subroutine init(self,dict)
    class(timeFilter), intent(inout)   :: self
    class(dictionary), intent(in)      :: dict
    real(defReal)                      :: T1, T2

    ! CE Case
    call dict % get(T1,'Tmin')
    call dict % get(T2,'Tmax')
    call self % build(T1, T2)

  end subroutine init

  !!
  !! Returns true if time value is between specified bounds
  !!
  elemental function isPass(self,state) result(passed)
    class(timeFilter), intent(in)    :: self
    class(particleState), intent(in) :: state
    logical(defBool)                 :: passed
    real(defReal)                    :: T

    ! CE paricle
    T = state % time
    passed = (self % Tmin <= T) .and. (T <= self % Tmax)

  end function isPass

  !!
  !! Build timeFilter for CE particles only from components
  !!
  !! Args:
  !!   Tmin [in] -> minimum time [s]
  !!   Tmax [in] -> maximum time [s]
  !!
  !! Errors:
  !!   fatalError if Tmin > Tmax
  !!
  subroutine build_CE(self, Tmin, Tmax)
    class(timeFilter), intent(inout)   :: self
    real(defReal), intent(in)          :: Tmin
    real(defReal), intent(in)          :: Tmax
    character(100), parameter :: Here = 'build_CE (timeFilter_class.f90)'

    self % Tmin = Tmin
    self % Tmax = Tmax

    ! Verify bounds
    if( self % Tmax <= self % Tmin) then
      call fatalError(Here,'Tmin='// numToChar(self % Tmin) //' is larger or equal to Tmax=' // numToChar(self % Tmax))
    end if

  end subroutine build_CE


end module timeFilter_class
