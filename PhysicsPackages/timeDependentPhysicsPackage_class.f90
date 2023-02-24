module timeDependentPhysicsPackage_class

  use numPrecision
  use universalVariables
  use endfConstants
  use genericProcedures,              only : fatalError, printFishLineR, numToChar, rotateVector
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Particle classes and Random number generator
  use particle_class,                 only : particle, P_NEUTRON
  use particleDungeon_class,          only : particleDungeon
  use source_inter,                   only : source
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                             gr_geomIdx  => geomIdx

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_activate    => activate ,&
                                             ndReg_display     => display, &
                                             ndReg_kill        => kill, &
                                             ndReg_get         => get ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase
  use neutronMaterial_inter,          only : neutronMaterial, neutronMaterial_CptrCast
  use ceNeutronMaterial_class,        only : ceNeutronMaterial
  use mgNeutronMaterial_inter,        only : mgNeutronMaterial
  use fissionCE_class,                only : fissionCE, fissionCE_TptrCast
  use fissionMG_class,                only : fissionMG, fissionMG_TptrCast
  use ceNeutronDatabase_inter,        only : ceNeutronDatabase, ceNeutronDatabase_CptrCast

  ! Operators
  use collisionOperator_class,        only : collisionOperator
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyAdmin_class,               only : tallyAdmin

  ! Factories
  use transportOperatorFactory_func,  only : new_transportOperator
  use sourceFactory_func,             only : new_source

  implicit none
  private

  !!
  !! Physics Package for time dependent calculations
  !!
  type, public,extends(physicsPackage) :: timeDependentPhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData => null()
    class(geometry), pointer               :: geom    => null()
    integer(shortInt)                      :: geomIdx = 0
    type(collisionOperator)                :: collOp
    class(transportOperator), allocatable  :: transOp
    class(RNG), pointer                    :: pRNG    => null()
    type(tallyAdmin),pointer               :: tally   => null()

    ! Settings
    integer(shortInt)  :: N_timeSteps
    integer(shortInt)  :: pop
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType
    integer(shortInt)  :: bufferSize
    real(defReal)      :: step_T
    integer(shortInt)  :: useCombing

    ! Calculation components
    type(particleDungeon), pointer :: thisTimeStep       => null()
    type(particleDungeon), pointer :: nextTimeStep       => null()
    type(particleDungeon), pointer :: temp_dungeon       => null()
    class(source), allocatable     :: timeDependent

    ! Timer bins
    integer(shortInt) :: timerMain
    real (defReal)     :: CPU_time_start
    real (defReal)     :: CPU_time_end

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: timeSteps
    procedure :: collectResults
    procedure :: run
    procedure :: kill

  end type timeDependentPhysicsPackage

contains

  subroutine run(self)
    class(timeDependentPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ TIME DEPENDENT CALCULATION /\/\"

    call self % timeSteps(self % tally, self % N_timeSteps, self % step_T)
    call self % collectResults()

    print *
    print *, "\/\/ END OF TIME DEPENDENT CALCULATION \/\/"
    print *
  end subroutine

  !!
  !!
  !!
  subroutine timeSteps(self, tally, N_timeSteps, step_T)
    class(timeDependentPhysicsPackage), intent(inout) :: self
    type(tallyAdmin), pointer,intent(inout)         :: tally
    integer(shortInt), intent(in)                   :: N_timeSteps
    integer(shortInt)                               :: i, n, nParticles
    type(particle), save                            :: p
    type(particleDungeon), save                     :: buffer
    type(collisionOperator), save                   :: collOp
    class(transportOperator), allocatable, save     :: transOp
    type(RNG), target, save                         :: pRNG
    real(defReal)                                   :: elapsed_T, end_T, T_toEnd, newTotalWeight
    real(defReal), intent(in)                       :: step_T
    integer(shortInt), dimension(N_timeSteps)       :: stepPopArray
    real(defReal), dimension(N_timeSteps)           :: stepWeightArray
    character(100),parameter :: Here ='time steps (timeDependentPhysicsPackage_class.f90)'
    !$omp threadprivate(p, buffer, collOp, transOp, pRNG)

    !$omp parallel
    ! Create particle buffer
    call buffer % init(self % bufferSize)

    ! Initialise neutron
    p % geomIdx = self % geomIdx
    p % k_eff = ONE

    ! Create a collision + transport operator which can be made thread private
    collOp = self % collOp
    transOp = self % transOp
    !$omp end parallel

    nParticles = self % pop

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    call self % timeDependent % generate(self % thisTimeStep, nParticles, self % pRNG)

    do i=1,N_timeSteps

      ! Send start of time step report
      if(self % printSource == 1) then
        call self % thisTimeStep % printToFile(trim(self % outputFile)//'_source'//numToChar(i))
      end if

      call tally % reportCycleStart(self % thisTimeStep)

      nParticles = self % thisTimeStep % popSize()
      
      ! previously set timeMax here but this caused omp issues
      

      !$omp parallel do schedule(dynamic)
      gen: do n = 1, nParticles

        ! TODO: Further work to ensure reproducibility!
        ! Create RNG which can be thread private
        pRNG = self % pRNG
        p % pRNG => pRNG
        call p % pRNG % stride(n)

        ! Obtain particle from dungeon
        call self % thisTimeStep % copy(p, n)

        !print *, 'previous particle timeMax: ', numToChar(p % timeMax)
        !print *, 'step_T:', numToChar(step_T)
        
        !print *, 'next particle timeMax: ', numToChar(p % timeMax)

        bufferLoop: do

          call self % geom % placeCoord(p % coords)
          
          ! Change particle timeMax
          p % timeMax = step_T * i
          
          !print *, '  NEXT PARTICLE'
          !print *, '    particle time: ', numToChar(p % time)
          !print *, '    particle timeMax:', numToChar(p % timeMax)
          !print *, '    pop in nextTimeStep:', numToChar(self % nextTimeStep % popSize())
          !print *, '    particles in buffer: ', numToChar(buffer % popSize())

          ! Save state
          call p % savePreHistory()

          p % fate = 0
          
          ! Transport particle untill its death
          history: do
            call transOp % transport(p, tally, buffer, buffer)

            ! Put particle in dungeon if past max time
            if(p % fate == AGED_FATE) then
              call self % nextTimeStep % detain(p)
              exit history
            endif

            if(p % isDead) exit history

            call collOp % collide(p, tally, buffer, buffer)
            if(p % isDead) exit history
          end do history

          ! Clear out buffer
          if (buffer % isEmpty()) then
            exit bufferLoop
          else
            call buffer % release(p)
          end if

        end do bufferLoop

      end do gen
      !$omp end parallel do
      
      ! Clear out this time step
      call self % thisTimeStep % cleanPop()

      ! Update RNG
      call self % pRNG % stride(self % pop)

      ! Send end of cycle report
      call tally % reportCycleEnd(self % thisTimeStep)

      if (self % useCombing == 0) then
          ! Temporary variable to hold total particle weight before population normalisation
          newTotalWeight = self % nextTimeStep % popWeight()

          print *, 'START Next time pop:', numToChar(self % nextTimeStep % popSize())
          print *, '      Next time weight:', numToChar(self % nextTimeStep % popWeight())

          ! Normalise population
          call self % nextTimeStep % normSize(self % pop, pRNG)
      
          print *, 'POP NORM Next time pop:', numToChar(self % nextTimeStep % popSize())
          print *, '         Next time weight:', numToChar(self % nextTimeStep % popWeight())
      
          ! Rescale the weight of each particle to match total weight before normalisation
          call self % nextTimeStep % normWeight(newTotalWeight)
          
          print *, 'WEIGHT NORM Next time pop:', numToChar(self % nextTimeStep % popSize())
          print *, '            Next time weight:', numToChar(self % nextTimeStep % popWeight())
      else
          print *, 'START Next time pop:', numToChar(self % nextTimeStep % popSize())
          print *, '      Next time weight:', numToChar(self % nextTimeStep % popWeight())
          
          call self % nextTimeStep % normCombing(self % pop, pRNG)
          
          print *, 'COMBING Next time pop:', numToChar(self % nextTimeStep % popSize())
          print *, '            Next time weight:', numToChar(self % nextTimeStep % popWeight())
      endif

      ! Add to array of weight
      stepWeightArray(i) = self % nextTimeStep % popWeight()
      
      ! Add to arrays of pop
      stepPopArray(i) = self % nextTimeStep % popSize()

      ! Flip timeStep dungeons
      self % temp_dungeon => self % nextTimeStep
      self % nextTimeStep    => self % thisTimeStep
      self % thisTimeStep    => self % temp_dungeon

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_timeSteps,defReal) * elapsed_T / i
      T_toEnd = max(ZERO, end_T - elapsed_T)
      
      nParticles = self % thisTimeStep % popSize()

      ! Display progress
      call printFishLineR(i)
      print *
      print *, 'Source batch:    ', numToChar(i), ' of ', numToChar(N_timeSteps)
      print *, 'Pop:             ', numToChar(self % thisTimeStep % popSize())
      print *, 'Next time pop:   ', numToChar(self % thisTimeStep % popSize())
      print *, 'Total weight:    ', numToChar(self % thisTimeStep % popWeight())
      print *, 'simulation time: ', numToChar(i * step_T)
      print *, 'Elapsed time:    ', trim(secToChar(elapsed_T))
      print *, 'End time:        ', trim(secToChar(end_T))
      print *, 'Time to end:     ', trim(secToChar(T_toEnd))
      call tally % display()
    end do
    print *
    print *, '-----------------------------------------------------------------'
    do i = 1, N_timeSteps
        print *, 'Time step ', numToChar(i), ' end pop/weight: ',&
                numToChar(stepPopArray(i)), ' / ', numToChar(stepWeightArray(i))
    end do
    print *, '-----------------------------------------------------------------'
    
  end subroutine timeSteps

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(timeDependentPhysicsPackage), intent(inout) :: self
    type(outputFile)                                :: out
    character(nameLen)                              :: name

    call out % init(self % outputFormat)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Source_batches'
    call out % printValue(self % N_timeSteps,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Transport_time'
    call out % printValue(timerTime(self % timerMain),name)

    ! Print tally
    call self % tally % print(out)

    call out % writeToFile(self % outputFile)

  end subroutine collectResults


  !!
  !! Initialise from individual components and dictionaries for source and tally
  !!
  subroutine init(self, dict)
    class(timeDependentPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)                :: dict
    class(dictionary),pointer                       :: tempDict
    integer(shortInt)                               :: seed_temp
    integer(longInt)                                :: seed
    character(10)                                   :: time
    character(8)                                    :: date
    character(:),allocatable                        :: string
    character(nameLen)                              :: nucData, energy, geomName
    type(outputFile)                                :: test_out
    character(100), parameter :: Here ='init (timeDependentPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_timeSteps,'timeSteps')
    call dict % get( nucData, 'XSdata')
    call dict % get( energy, 'dataType')
    call dict % get( self % step_T, 'timeStepLength')

    ! Process type of data
    select case(energy)
      case('mg')
        self % particleType = P_NEUTRON_MG
      case('ce')
        self % particleType = P_NEUTRON_CE
      case default
        call fatalError(Here,"dataType must be 'mg' or 'ce'.")
    end select

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Parallel buffer size
    call dict % getOrDefault( self % bufferSize, 'buffer', 10)
    
    ! Whether to use combing (default = no)
    call dict % getOrDefault(self % useCombing, 'combing', 0)

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Initialise RNG
    allocate(self % pRNG)

    ! *** It is a bit silly but dictionary cannot store longInt for now
    !     so seeds are limited to 32 bits (can be -ve)
    if( dict % isPresent('seed')) then
      call dict % get(seed_temp,'seed')

    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date // time
      call FNV_1(string,seed_temp)

    end if
    seed = seed_temp
    call self % pRNG % init(seed)

    ! Read whether to print particle source per cycle
    call dict % getOrDefault(self % printSource, 'printSource', 0)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'timeDependentGeom'
    call gr_addGeom(geomName, tempDict)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, self % geom % activeMats())
    self % nucData => ndReg_get(self % particleType)

    ! Read particle source definition
    tempDict => dict % getDictPtr('source')
    call new_source(self % timeDependent, tempDict, self % geom)

    ! Build collision operator
    tempDict => dict % getDictPtr('collisionOperator')
    call self % collOp % init(tempDict)

    ! Build transport operator
    tempDict => dict % getDictPtr('transportOperator')
    call new_transportOperator(self % transOp, tempDict)

    ! Initialise tally Admin
    tempDict => dict % getDictPtr('tally')
    allocate(self % tally)
    call self % tally % init(tempDict)

    ! Size particle dungeon
    allocate(self % thisTimeStep)
    call self % thisTimeStep % init(10 * self % pop)
    
    ! Size particle dungeon
    allocate(self % nextTimeStep)
    call self % nextTimeStep % init(10 * self % pop)

    call self % printSettings()

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(timeDependentPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(timeDependentPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ TIME DEPENDENT CALCULATION /\/\"
    print *, "Source batches:       ", numToChar(self % N_timeSteps)
    print *, "Population per batch: ", numToChar(self % pop)
    print *, "Initial RNG Seed:     ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings

end module timeDependentPhysicsPackage_class
