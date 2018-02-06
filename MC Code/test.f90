program test

  use numPrecision
  use RNG_class
  use genericProcedures
  use ByIsoNoMT_Data_class, only : byIsoNoMT_Data
  use aceNoMT_class

  use endfTable_class, only : endfTable

  use releaseLawENDF_class
  use constantRelease_class, only : constantRelease
  use polynomialRelease_class, only : polynomialRelease
  use tabularRelease_class, only : tabularRelease

  use muEndfPdf_class,   only : muEndfPdf, muEndfPdf_ptr
  use isotropicmu_class, only : isotropicmu
  use equiBin32mu_class, only : equiBin32mu
  use tabularmu_class,   only : tabularmu
  use tabularPdf_class,    only : tabularPdf
  use angleLawENDF_class,   only : angleLawENDF
  use tabularAngle_class, only : tabularAngle
  use noAngle_class,       only: noAngle

  use tabularEnergy_class, only: tabularEnergy
  use contTabularEnergy_class, only : contTabularEnergy
  use energyLawENDF_class,     only : energyLawENDF
  use noEnergy_class,          only : noEnergy

  use dictionary_class , only: dictionary, dictContent
  use IOdictionary_class, only : IOdictionary


  implicit none
  integer(kind=shortInt) :: i
  integer(longInt)       :: longI, longI2, rate
  integer(kind=shortInt),dimension(99),target :: A
  integer(kind=shortInt),dimension(:),pointer :: B
  INTEGER(SHORTiNT),DIMENSION(:), allocatable :: C
  type(ByIsoNoMT_Data)  :: CEdata
  character(len=pathLen)      :: matInput="./testInputLarge"
  character(len=pathLen)      :: isoInput="/home/mak60/myACE/JEF311.aceXS"
  character(len=99)      :: format
  character(len=99),dimension(2) :: Ach
  character(len=pathLen)    :: acePath = "/home/mak60/myACE/acedata/33075JEF311.ace"
  character(pathLen)       :: testDictFile = "./testDictInput"
  integer(shortInt)         :: firstLine = 1170286
  type(aceNoMT)             :: isotope
  real(defReal) :: kl, eps
  real(defReal),dimension(:),allocatable :: x, pdf, x2, pdf2,R

  type(RNG) :: random


  class(releaseLawENDF),pointer  :: release
  class(endfTable),pointer       :: eTable

  real(defReal), dimension(20) :: energy, second
  integer(shortInt), dimension(4) :: bounds, ENDF

  real, pointer :: p1,p2,p3

  type(muEndfPdf_ptr) :: myPtr, myPtr2
  type(tabularEnergy), pointer :: tabPtr
  class(angleLawENDF),pointer :: angle
  class(energyLawENDF), pointer :: energyT
  type(tabularEnergy),dimension(:),allocatable   :: tables
  type(contTabularEnergy) :: yjfj

  class(*),pointer   :: GP(:), GP2(:)
  real,pointer       :: a1(:),a2(:)
  real,pointer       :: rp1(:),rp2(:)
  character(:),pointer :: char_ptr

  type(dictionary) :: testDict
  type(dictionary) :: testDict2
  type(dictionary) :: testDict3
  type(dictionary) :: testDict4

  character(nameLen) :: abc = "aaa", bc
  character(20),dimension(:),allocatable  :: charT
  character(20),dimension(:),pointer      :: cA_ptr
  class(*),dimension(:),pointer           :: Uptr
  class(*),dimension(:),pointer           :: Uptr2
  character(20),dimension(:),pointer      :: cA_ptr2
  character(20),dimension(:),pointer      :: cA_ptr3
  character(20),dimension(:),pointer      :: cA_ptr4
  character(20),dimension(:),pointer      :: localPointer
  character(20),dimension(:), pointer     :: newMemory

  character(20),dimension(:),allocatable :: cA_alloc1
  character(20),dimension(:),allocatable :: cA_alloc2
  type(dictContent)     :: dictNode
  type(IOdictionary)    :: IOdictTest
  integer(shortInt) :: err
  character(100) :: msg


  charT = ['jgfjfjfj' ,&
           'kjgkgk  ' ,&
           '  ugkkik' ,&
           '  jgkj  ']



!  bc= charT(2)
!  print *, bc
!  stop
 !print *, index(trim(abc)," ")
!  abc = '7'
!  read (unit=abc,fmt="(ES30.30)",iostat=err,iomsg=msg) kl
!  print *, kl
!  print *, err
!  print *, msg
!  stop



  call IOdictTest % initFrom(testDictFile)

  print *, "TOP DICTIONARY"
  print *, IOdictTest % getRealArray('list1')
  print *, IOdictTest % getReal('keyword1')
  print *, IOdictTest % getInt('keyword2')
  print *, IOdictTest % getChar('keyword3')
  print *, IOdictTest % getInt('keyword4')
  print *, IOdictTest % getInt('keyword7')

  testDict = IOdictTest % getDict('dictionary1')

  print *, "NESTED DICTIONARY "
  print *, testDict % getInt('keyword')
  print *, testDict % getReal('keyword2')
  print *, testDict % getChar('keyword3')

  testDict2 = testDict % getDict('dictionary3')

  print *, "NESTED DICTIONARY2 "
  print *, testDict2 % getInt('key')
  print *, testDict2 % getInt('key2')

  testDict3 = testDict2 % getDict('dictionary4')

  print *, "NESTED DICTIONARY3"
  print *, testDict3 % getInt('key')
  print *, testDict3 % getReal('key2')

  testDict4 = testDict2 % getDict('dictionary5')

  print *, "NESTED DICTONARY 3.1"
  print *, testDict4 % getReal('key')
  print *, testDict4 % getInt('key2')


  stop
!  allocate(cA_ptr(size(charT)))
!  cA_ptr = charT
!
!  Uptr => cA_ptr
!
!  cA_ptr => null()
!
!  select type(Uptr)
!   type is (character(*))
!
!     localPointer(1:size(Uptr)) => Uptr(1:size(Uptr))
!     cA_alloc2 = localPointer
!     print *, cA_alloc2
!
!     allocate(newMemory(size(Uptr)) )
!     newMemory(1:size(Uptr)) =''
!
!     newMemory = localPointer
!     print *, localPointer
!     print *, newMemory
!
!  end select
!     Uptr2 => newMemory
!
!  select type(Uptr)
!    type is (character(*))
!      print *, size(Uptr)
!      cA_ptr2(1:size(Uptr)) => Uptr(1:size(Uptr))
!  end select
!
!  select type(Uptr2)
!    type is (character(*))
!        print *, size(Uptr2)
!      cA_ptr3(1:size(Uptr2)) => Uptr2(1:size(Uptr2))
!        print *, cA_ptr3
!  end select
!
!  newMemory(3) = 'HaHA'
!
!
!  print *, cA_ptr2
!  print *, cA_ptr3
!
!
!  stop
  call testDict % init(40,2)
  call testDict % store("aaa",68.9_8)
  call testDict % store("bbbb",76.8_8)
  call testDict % store("ccc",8.9_8)
  call testDict % store("ddd",8.6_8)
  call testDict % store("eee",1.2_8)
  call testDict % store("fff",[2, 4,5,6,7,8,9,1,2,3,4,5,6])
  call testDict % store("gggg",3.2_8)
  call testDict % store("hh",4.2_8)
  call testDict % store("i","jklJKL" )
  call testDict % store("jKL",charT)
  call testDict % store("OPR",7)
  call testDict % store("KLM",[3.4_8,4.5_8,1.2_8,4.3_8])


  call testDict2 % init(20)
  call testDict2 % store("aaa",testDict)

  call testDict % kill()

  testDict3 = testDict2 % getDict("aaa")

  print *, testDict3 % getReal("ddd")



  !call deepCopy(dictNode,testDict%entries(1))
!  call dictNode % deepCopy(testDict % entries(10))
!  call testDict % kill()
!  call testDict % init(40,2)
!  testDict % keywords(1) = 'a'
!  testDict % entries(1) = dictNode
!  call testDict%entries(2) % shallowCopy(dictNode)
!  !call shallowCopy(testDict%entries(2),dictNode)
!
!  print *, testDict % getCharArray('a')

  !print *, testDict % getCharArray('jKL')

!  allocate(a1(2))
!  allocate(a2(2))
!
!  a1 = [6.0 , 7.9]
!
!  GP => a1
!  GP2 => GP
!  rp1 => a1
!  a1 => null()
!  select type (GP2)
!    type is(real)
!
!      rp1 => GP2
!    type is(integer)
!     print *, 'kekekekekekekekeke'
!  end select
!
!  !deallocate(GP)
!
!  print*, rp1
!  i =10
!  allocate(character(i) :: char_ptr)
!
!  char_ptr = '12345678910'
!
!  print *, char_ptr
!
!  a1 => null()
!  rp1 => null()
!
!  print *, associated(rp1)


  !C=[1,2,3,4,5,6,7,8,9,10]
  !B => C(1:8)
  !print *, B(3:5)
  stop

  call isotope % init(acePath,1)
  !stop
  call CEdata % readFrom(matInput,isoInput)
  call CEdata % print()

  !call isotope % init(acePath,firstLine)


  !R= [1.1_8, 2.3_8, 3.6_8, 9.6_8, 11.9_8]

  !print *, binaryFloorIdxC(R,1.1_8)

  !C = [(2*i,i=1,10)]



  energy = [(real(i),i=1,20)]
  second = [(i*5.0,i=1,20)]
  bounds = [ 5, 10, 15, 20]
  ENDF = [5,1,3,4]
   !bounds = [20]
   !ENDF = [1]


!  second = [1.0_8,1.5_8,2.0_8,1.5_8,1.7_8,2.5_8,2.3_8,2.0_8,2.1_8,1.5_8]
!  energy = [(real(i),i=1,10)]
!  bounds = [3,7,10]
!  ENDF = [1,2,1]
!  !energy(10)= energy(9)
!  energy(4) = energy(3)
!  energy(6) = energy(5)
!

!  release => tabularRelease(energy,second,bounds,ENDF)!
!
!  eTable  => endfTable(energy,second,bounds,ENDF)
!
!  do i=1,100
!    kl = (19.99-1.0)/100.0 * i + 1.0
!    print *, kl, release % releaseAt(kl), eTable % at(kl), release % releaseAt(kl)-eTable % at(kl)
!  end do


!  x = [(-1.0+2.0/10*i,i=0,10)]
!  pdf = abs(x)
!
!  x2   = [ -1.0_8, 1.0_8]
!  pdf2 = [ 0.5_8, 0.5_8]
!
!  x = [ -1.0_8, 0.0_8, 1.0_8]
!  pdf = [ 1.0_8/3 ,2.0_8/3, 76876.0_8]

  !call table % init(x,pdf,0)


!
!  R(1) = -1.0
!  myPtr = equiBin32mu(R)
!  myPtr  = tabularmu(x2,pdf2,1)
!  myPtr2 = tabularmu(x,pdf,1)
!
!  allocate(tables(2))
!  call tables(1) % init (x2,pdf2,1)
!  call tables(2) % init (x,pdf,1)
!
!
!  angle  => tabularAngle([0.0_8, 1.0_8],[myPtr, myPtr2])
!  energyT => contTabularEnergy([0.0_8, 1.0_8],tables)
!
!  deallocate(energyT)
!
!  energyT => noEnergy()
!  angle => noAngle()
! ! myPtr = tabularmu(x,pdf,0)
!
!  do i=0,1000
!    kl = 2.0/1000 * i - 1.0
!     !eps = random % get()
!
!  !   print *, kl, angle % probabilityOf(kl,1.0_8), energyT % probabilityOf(kl,0.5_8)
!
!  end do



  !print *, binarySearch(energy,1.0_8)

  !print *, linearFloorIdx(C,21)
!  call system_clock(count_rate=rate)
!  call system_clock(count=longI)
!  do i=1,1
!    kl= interpolate(0.00001_8, 1.0_8, 0.01_8, 10.0_8, 0.5_8)
!  end do
!  call system_clock(count=longI2)
!  print *, 'end', (longI2-longI)/real(rate)
end program test

