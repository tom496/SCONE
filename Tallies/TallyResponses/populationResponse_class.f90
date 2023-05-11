module populationResponse_class

  use numPrecision
  use dictionary_class,    only : dictionary
  use particle_class,      only : particle
  use tallyResponse_inter, only : tallyResponse

  ! Nuclear Data interface
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! tallyResponse to score population
  !! needs to be divided by time step in post-processing
  !!
  !! Returns 1/particle speed
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!
  type, public,extends(tallyResponse) :: populationResponse
    private
  contains
    procedure :: init
    procedure :: get
    procedure :: kill
  end type populationResponse

contains

  !!
  !! Initialise Response from dictionary
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(populationResponse), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    ! Do nothing

  end subroutine init

  !!
  !! Get 1/particle speed
  !!
  !! See tallyResponse_inter for details
  !!
  function get(self, p, xsData) result(val)
    class(populationResponse), intent(in)       :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: val

    val = ONE/p % getSpeed()

  end function get

  !!
  !! Return to uninitialised State
  !!
  elemental subroutine kill(self)
    class(populationResponse), intent(inout) :: self

    ! Do nothing for nothing can be done

  end subroutine kill

end module populationResponse_class
