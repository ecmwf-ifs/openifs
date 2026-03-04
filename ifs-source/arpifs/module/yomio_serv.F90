! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMIO_SERV

!**** *YOMIO_SERV*  - Defines io_serv main & ancillary derived types

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO-FRANCE*
!      Original : 01-01-2011

!      Modifications :
!      P.Marguinaud : 11-09-2012 : Cleaning
!      P.Marguinaud : 10-10-2013 : Remove fifo stuff, add support for
!                                  distribution over a fraction of the IO server
!                                  tasks
!      P.Marguinaud : 10-10-2014 : Cleaning
!      R. El Khatib : 09-03-2015 : flexible directory for output xml files
!      P.Marguinaud : 04-10-2016 : Port to single precision

USE PARKIND1, ONLY : JPIM, JPRB, JPIB, JPRD

USE YOMIO_SERV_REQ, ONLY : IO_SERV_REQ_PTR
USE YOMFP_SERV_DINF, ONLY : FP_SERV_DINF
USE MPL_MPIF, ONLY : MPI_COMM_NULL

USE PARFPOS, ONLY : JPOSLEN

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: NIO_SERV_SEND_BLOCKING_STD     = 1, &
                               & NIO_SERV_SEND_NON_BLOCKING_STD = 2

! processing levels
INTEGER(KIND=JPIM), PARAMETER :: NIO_SERV_PROCESS_FRST = 0, &
                               & NIO_SERV_PROCESS_NONE = 0, & ! IO server starts, but does nothing
                               & NIO_SERV_PROCESS_RECV = 1, & ! IO server receives
                               & NIO_SERV_PROCESS_COMP = 3, & ! IO server compresses data
                               & NIO_SERV_PROCESS_WRIT = 5, & ! IO server writes data
                               & NIO_SERV_PROCESS_LAST = 5

! processing levels ids
CHARACTER(LEN=10), PARAMETER :: CIO_SERV_PROCESS_LABL(NIO_SERV_PROCESS_FRST:NIO_SERV_PROCESS_LAST) &
                             & = (/ 'NONE      ', 'RECV      ', 'RECV--COMP', &
                             &      'COMP      ', 'COMP--WRIT', 'WRIT      ' /)

! Fortran unit used 
INTEGER (KIND=JPIM), PARAMETER :: NUNIT = 9_JPIM

! Base MPI tag
INTEGER(KIND=JPIM), PARAMETER :: NIO_SERV_TAG_ZERO = 0_JPIM

! Buffer allocated & managed by io_serv
TYPE IO_SERV_MEM_BLOCK
  REAL(KIND=JPRB), POINTER :: R(:,:) => NULL ()
! MPI request id
  INTEGER(KIND=JPIM)       :: NREQID = 0
! MPI send initiated ?
  LOGICAL                  :: LREQID = .FALSE.
END TYPE

! FA defaults (see arguments of FAVORI & FAGIOT)
TYPE IO_SERV_FAGIOT_ARGS
  INTEGER(KIND=JPIM) :: INGRIB = -99_JPIM             !   4
  INTEGER(KIND=JPIM) :: INBPDG = -99_JPIM             !   8
  INTEGER(KIND=JPIM) :: INBCSP = -99_JPIM             !  12
  INTEGER(KIND=JPIM) :: ISTRON = -99_JPIM             !  16
  INTEGER(KIND=JPIM) :: IPUILA = -99_JPIM             !  20
  INTEGER(KIND=JPIM) :: IDMOPL = -99_JPIM             !  24
END TYPE IO_SERV_FAGIOT_ARGS

TYPE IO_SERV_FAFRAME
! All parameters necessary for defining a FA frame
  INTEGER (KIND=JPIM)          :: NTYPTR = -99_JPIM
  INTEGER (KIND=JPIM)          :: NTRONC = -99_JPIM
  INTEGER (KIND=JPIM)          :: NNLATI = -99_JPIM
  INTEGER (KIND=JPIM)          :: NNXLON = -99_JPIM
  INTEGER (KIND=JPIM)          :: NNIVER = -99_JPIM
  INTEGER (KIND=JPIM), POINTER :: NNLOPA (:) => NULL ()
  INTEGER (KIND=JPIM), POINTER :: NNOZPA (:) => NULL ()
  REAL (KIND=JPRB)             :: RSLAPO = -99._JPRB
  REAL (KIND=JPRB)             :: RCLOPO = -99._JPRB
  REAL (KIND=JPRB)             :: RSLOPO = -99._JPRB
  REAL (KIND=JPRB)             :: RCODIL = -99._JPRB
  REAL (KIND=JPRB)             :: RREFER = -99._JPRB
  REAL (KIND=JPRB),    POINTER :: RSINLA (:) => NULL ()
  REAL (KIND=JPRB),    POINTER :: RAHYBR (:) => NULL ()
  REAL (KIND=JPRB),    POINTER :: RBHYBR (:) => NULL ()
  CHARACTER(LEN=32)            :: CNC = ""   ! frame name
END TYPE

TYPE IO_SERV_ECGRIB
! All parameters necessary for EC GRIB IO
  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_SH = -99_JPIM
  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_SH_ML = -99_JPIM
  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_GG = -99_JPIM
  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_GG2 = -99_JPIM
  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_GG_ML = -99_JPIM
  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_BUDG = -99_JPIM
  
  ! From YOMGRIB
  INTEGER (KIND=JPIM), POINTER :: NSFLEVS(:,:)=>NULL()
  INTEGER (KIND=JPIM)          :: NLOCGRB
  INTEGER (KIND=JPIM)          :: NSTREAM
  INTEGER (KIND=JPIM)          :: NLEG
  INTEGER (KIND=JPIM)          :: NREFERENCE
  INTEGER (KIND=JPIM)          :: NWINOFF_4V
  ! From YOMVAREPS
  INTEGER (KIND=JPIM)          :: NFCHO_TRUNC_INI
  INTEGER (KIND=JPIM)          :: NFCLENGTH_INI
  LOGICAL                      :: LVAREPS
  ! From YOMCT0
  INTEGER (KIND=JPIM)          :: NSTEPINI
  LOGICAL                      :: LFDBOP
  LOGICAL                      :: LSMSSIG
  CHARACTER (LEN = 25)         :: CMETER
  CHARACTER (LEN = 2)          :: CTYPE
  ! From YOMRIP
  REAL(KIND=JPRB)              :: TSTEP
  ! From YOMDYNCORE
  LOGICAL                      :: LPPSTEPS
  ! From YOMSATSIM
  INTEGER(KIND=JPIM), POINTER  :: MSERIES(:)=>NULL()
  INTEGER(KIND=JPIM), POINTER  :: MSATID(:)=>NULL() 
  INTEGER(KIND=JPIM), POINTER  :: MINST(:)=>NULL() 
  INTEGER(KIND=JPIM), POINTER  :: MCHAN(:)=>NULL() 
  REAL(KIND=JPRB),    POINTER  :: RCWN(:)=>NULL() 
  END TYPE IO_SERV_ECGRIB
  
  TYPE IO_SERV_DMPARAM
  INTEGER (KIND=JPIM) :: ISZGPG = 0_JPIM  ! Size of a grid-point field
  INTEGER (KIND=JPIM) :: ISZSPG = 0_JPIM  ! Size of a spectral field  
  INTEGER (KIND=JPIM) :: IADDPK = 0_JPIM  ! Extra 8-byte words for packing field
  ! (the maximum size of a packed field
  ! is ISZSPG+IADDPK, or ISZGPG+IADDPK)
  CHARACTER (LEN=JPOSLEN) :: CFPDOM  = '' ! YOMFPC
END TYPE
  
TYPE IO_SERV_WAVEMODEL
  LOGICAL                     :: LWACTIVE            ! True if wave model is active
  INTEGER(KIND=JPIM)          :: NTOTMX              ! Max grid points per proc
  INTEGER(KIND=JPIM)          :: NTOTG               ! Total grid points
  INTEGER(KIND=JPIM)          :: NGX                 ! YOWPARAM
  INTEGER(KIND=JPIM)          :: NGY                 ! YOWPARAM
  CHARACTER(LEN=1)            :: CLDOMAIN            ! YOWPARAM
  INTEGER(KIND=JPIM)          :: IU06                ! YOWTEST
  INTEGER(KIND=JPIM)          :: ITEST               ! YOWTEST
  LOGICAL                     :: LGRHDIFS            ! YOWGRIBHD
  REAL(KIND=JPRB)             :: PPMISS              ! YOWGRIBHD
  REAL(KIND=JPRB)             :: PPEPS               ! YOWGRIBHD
  REAL(KIND=JPRB)             :: PPREC               ! YOWGRIBHD
  REAL(KIND=JPRB)             :: PPRESOL             ! YOWGRIBHD
  REAL(KIND=JPRB)             :: PPMIN_RESET         ! YOWGRIBHD
  INTEGER(KIND=JPIM)          :: NTENCODE            ! YOWGRIBHD
  INTEGER(KIND=JPIM)          :: NGRBRESS            ! YOWGRIBHD
  LOGICAL                     :: LNEWLVTP            ! YOWGRIBHD
  LOGICAL                     :: LPADPOLES           ! YOWGRIBHD
  INTEGER(KIND=JPIM), ALLOCATABLE :: NLONRGG(:)      ! YOWGRID
  INTEGER(KIND=JPIM)          :: IRGG                ! YOWMAP
  REAL(KIND=JPRB)             :: AMONOP              ! YOWMAP
  REAL(KIND=JPRB)             :: AMOSOP              ! YOWMAP
  REAL(KIND=JPRB)             :: XDELLA              ! YOWMAP
  REAL(KIND=JPRB)             :: ZMISS               ! YOWPCONS
  INTEGER(KIND=JPIM)          :: NSTPW               ! YOEWCOU
  CHARACTER(LEN=256)          :: CFDB2DSP            ! YOWSTAT
  INTEGER(KIND=JPIM)          :: NPRECI              ! YOWMPP

  INTEGER(KIND=JPIM)          :: IMDLGRBID_G
  INTEGER(KIND=JPIM)          :: IMDLGRBID_M
  
  INTEGER(KIND=JPIM), POINTER :: ISORTL2G (:,:) => NULL () !  NPROC x      1 -> NGPTOT
  INTEGER(KIND=JPIM), POINTER :: IRANKSET (:,:) => NULL () !  NPROC x      1 -> MYPROC
  INTEGER(KIND=JPIM), POINTER :: ISORTCNT (:)   => NULL () !  NPROC
  INTEGER(KIND=JPIM), POINTER :: ISORTOFF (:)   => NULL () !  NPROC

  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_WAM_S = -99_JPIM
  INTEGER (KIND=JPIM)          :: NGRIB_HANDLE_WAM_I = -99_JPIM
  
END TYPE IO_SERV_WAVEMODEL
  
! Model parameters passed to io_serv (FA frames + a few module variables)
TYPE IO_SERV_MODELPAR
  INTEGER(KIND=JPIM)          :: NTIMEFMT       = 0       ! YOMOPH
  INTEGER(KIND=JPIM)          :: NSPEC2G        = 0       ! YOMDIM
  INTEGER(KIND=JPIM)          :: NSPEC2MX       = 0       ! YOMDIM
  INTEGER(KIND=JPIM)          :: NGPTOTG        = 0       ! YOMGEM
  INTEGER(KIND=JPIM)          :: NGPTOTMX       = 0       ! YOMGEM
  LOGICAL                     :: LECMWF         = .FALSE. ! YOMCT0
  LOGICAL                     :: LARPEGEF       = .TRUE.  ! YOMCT0
  LOGICAL                     :: LELAM          = .FALSE. ! YOMCT0
  CHARACTER(LEN=16)           :: CNMEXP         = ''      ! YOMCT0
  INTEGER(KIND=JPIM)          :: NPRTRW         = 0       ! YOMCT0
  INTEGER(KIND=JPIM)          :: NPRTRV         = 0       ! YOMCT0
  LOGICAL                     :: LWRSPECA_GP    = .FALSE. ! YOMCT0
  LOGICAL                     :: LSUSPECA_GP    = .FALSE. ! YOMCT0
  LOGICAL                     :: LINC           = .FALSE. ! YOMOP
  INTEGER(KIND=JPIM)          :: NFPDOM         = 0       ! YOMFPC
  INTEGER(KIND=JPIM)          :: NFPRGPLX       = 0       ! YOMFPG
  INTEGER(KIND=JPIM)          :: MBX_SIZE       = 0       ! YOMMP0
  INTEGER(KIND=JPIM)          :: NSMAX                    ! YOMDIM

! PROGRID parameters
  CHARACTER(LEN=64)  :: CMODEL   = ''
  INTEGER(KIND=JPIM) :: NIDCEN   = 0
  LOGICAL            :: LEXTERN  = .FALSE.
! FA compression defaults
  TYPE (IO_SERV_FAGIOT_ARGS) :: YFAGIOT_ARGS
! FA parameters; model
  TYPE (IO_SERV_FAFRAME) :: YFRAME
  TYPE (IO_SERV_DMPARAM) :: YPARAM
! FA parameters; fullpos
  TYPE (IO_SERV_FAFRAME), POINTER :: YFRAMEFP (:) => NULL ()
  TYPE (IO_SERV_DMPARAM), POINTER :: YPARAMFP (:) => NULL ()
! FA limits
  INTEGER (KIND=JPIM) :: JPXTRO = 0
  INTEGER (KIND=JPIM) :: JPXLAT = 0
  INTEGER (KIND=JPIM) :: JPXNIV = 0
  INTEGER (KIND=JPIM) :: JPNXFA = 0
  INTEGER (KIND=JPIM) :: JPNXCA = 0
! EC_GRIBIO
  TYPE (IO_SERV_ECGRIB) :: YECGRIB
! Wave model parameters
  TYPE (IO_SERV_WAVEMODEL) :: YWAM
! Model distribution
  INTEGER(KIND=JPIM), POINTER :: ISPSORTL2G (:,:) => NULL () ! Spectral
  INTEGER(KIND=JPIM), POINTER :: ISPRANKSET (:,:) => NULL () ! NPRTRW x NPRTRV -> MYPROC
  INTEGER(KIND=JPIM), POINTER :: ISPSORTCNT (:)   => NULL () ! NPRTRW 
  INTEGER(KIND=JPIM), POINTER :: ISPSORTOFF (:)   => NULL () ! NPRTRW 

  INTEGER(KIND=JPIM), POINTER :: IGPSORTL2G (:,:) => NULL () ! Grid-point
  INTEGER(KIND=JPIM), POINTER :: IGPRANKSET (:,:) => NULL () !  NPROC x      1 -> MYPROC
  INTEGER(KIND=JPIM), POINTER :: IGPSORTCNT (:)   => NULL () !  NPROC
  INTEGER(KIND=JPIM), POINTER :: IGPSORTOFF (:)   => NULL () !  NPROC

  INTEGER(KIND=JPIM), POINTER :: IFPSORTL2G (:,:) => NULL () ! Fullpos
  INTEGER(KIND=JPIM), POINTER :: IFPRANKSET (:,:) => NULL () !  NPROC x      1 -> MYPROC
  INTEGER(KIND=JPIM), POINTER :: IFPSORTCNT (:)   => NULL () !  NPROC 
  INTEGER(KIND=JPIM), POINTER :: IFPSORTOFF (:)   => NULL () !  NPROC 
END TYPE

TYPE IO_SERV

! IO server has sent some data
  LOGICAL :: LSENTFLD = .FALSE.

! common

! log levels
  LOGICAL :: LWARN = .FALSE.
  LOGICAL :: LINFO = .FALSE.
  LOGICAL :: LDBUG = .FALSE.

! log open ?
  LOGICAL :: LLOG_OPEN = .FALSE.

! log base time
  REAL(KIND=JPRB) :: ZBTIME

  INTEGER(KIND=JPIM) :: NMSG_LEVEL        = 1
  LOGICAL            :: LMSG_FLUSH        = .TRUE.

  INTEGER(KIND=JPIM) :: NMSG_LEVEL_SERVER = 0
  INTEGER(KIND=JPIM) :: NMSG_LEVEL_CLIENT = 0
  LOGICAL            :: LMSG_FLUSH_SERVER = .TRUE.
  LOGICAL            :: LMSG_FLUSH_CLIENT = .FALSE.

! Time after which a warning will be printed by compute tasks
! reporting time in IO server routines (in seconds)
  REAL(KIND=JPRD)    :: WAIT_TIME_THRESHOLD = 0.5_JPRD

! processing level (see enumeration above NIO_SERV_PROCESS_NONE ...)
  INTEGER(KIND=JPIM) :: NPROCESS_LEVEL = NIO_SERV_PROCESS_LAST

! Model MPI #1 + IO tasks communicator
  INTEGER(KIND=JPIM) :: NCOMM_W1IO  = MPI_COMM_NULL
! Model communicator
  INTEGER(KIND=JPIM) :: NCOMM_WR    = MPI_COMM_NULL
! IO server communicator
  INTEGER(KIND=JPIM) :: NCOMM_IO    = MPI_COMM_NULL
! MPI_COMM_WORLD communicator
  INTEGER(KIND=JPIM) :: NCOMM_WRIO  = MPI_COMM_NULL
! Rank of current task in MPI_COMM_WORLD
  INTEGER(KIND=JPIM) :: MYPROC_WRIO = -1
! Rank of current task in the IO group
  INTEGER(KIND=JPIM) :: MYPROC_IO   = -1
! Rank of current task in the working group (==YOMMP%MYPROC)
  INTEGER(KIND=JPIM) :: MYPROC_WR   = -1
! Number of tasks in COMM WORLD
  INTEGER(KIND=JPIM) :: NPROC_WRIO  = -1
! TRUE if the current task is running an IO server
  LOGICAL            :: LIO_SERVER  = .FALSE.
! TRUE is the current task is running an IO client
  LOGICAL            :: LIO_CLIENT  = .FALSE.
! True if gathered arrays are (NGPTOTG|NSPEC2G,NFLDS)
  LOGICAL            :: LARRAY2D    = .FALSE.
  
  INTEGER(KIND=JPIM) :: NPROC_IO = 0
! Ranks of all tasks running an IO server in MPI_COMM_WORLD`
  INTEGER(KIND=JPIM), POINTER :: MYPROCS_IO(:) => NULL ()
  
  INTEGER(KIND=JPIM) :: NPROC_WR = 0
! Ranks of all working tasks (running the forecast) in MPI_COMM_WORLD
  INTEGER(KIND=JPIM), POINTER :: MYPROCS_WR(:) => NULL ()

! Number of threads to be used by each IO server
  INTEGER(KIND=JPIM) :: NTHREAD_IO = 1

! Send specific

! Tag number (incremented automatically)
  INTEGER(KIND=JPIM) :: NIO_SERV_TAG         = NIO_SERV_TAG_ZERO
! NIO_SERV_METHOD : Blocking standard     = 1
!                   Non-blocking standard = 2
  INTEGER(KIND=JPIM) :: NIO_SERV_METHOD      = 2

! Max number of blocks that can be simultaneously allocated
  INTEGER(KIND=JPIM) :: NIO_SERV_BUF_MAXSIZE = 10 ! blocks
! Current number of allocated blocks
  INTEGER(KIND=JPIM) :: NIO_SERV_BUF_CURSIZE = 0

! List of currently allocated blocks (1:NIO_SERV_BUF_CURSIZE)
! size of IO_SERV_BUF is NIO_SERV_BUF_MAXSIZE
  TYPE(IO_SERV_MEM_BLOCK), POINTER :: IO_SERV_BUF (:) => NULL ()

! Number of send messages
  INTEGER(KIND=JPIM) :: ISEND = 0
! Number of 8-byte words send
  INTEGER(KIND=JPIB) :: ISENDSIZE = 0

! Recv specific

! Number of created IO_SERV_REQ since the beginning
  INTEGER(KIND=JPIM) :: IREQ_CREAT = 0
! Number of completed IO_SERV_REQ since the beginning
  INTEGER(KIND=JPIM) :: IREQ_COMPL = 0

! Number of received messages
  INTEGER(KIND=JPIM) :: IRECV = 0
! Number of 8-byte words received
  INTEGER(KIND=JPIB) :: IRECVSIZE = 0
! Send single field descriptors
  LOGICAL :: LFLDDESC_UNIQ    = .TRUE.
! Check consistency of multiple field descriptors
  LOGICAL :: LFLDDESC_CHECK   = .TRUE.
! Pending request list
  TYPE (IO_SERV_REQ_PTR), POINTER :: YREQ_PTR_LIST (:) => NULL ()

! Model parameters
  TYPE (IO_SERV_MODELPAR) :: MODELPAR
! Fraction of IO procs dedicated to model fields; the rest is for Fullpos
  REAL (KIND=JPRB) :: PIOPROCR_MDL = -1._JPRB
! IO ranks range dedicated to model
  REAL (KIND=JPRB) :: PIOPROC1_MDL = -1._JPRB
  REAL (KIND=JPRB) :: PIOPROC2_MDL = -1._JPRB
! IO ranks range dedicated to fullpos
  REAL (KIND=JPRB) :: PIOPROC1_FLP = -1._JPRB
  REAL (KIND=JPRB) :: PIOPROC2_FLP = -1._JPRB

! Date being processed

  INTEGER(KIND=JPIM) :: IDATEF (22)

! Input files directory

  CHARACTER (LEN=256) :: CIFDIR = ''

! Output files directory

  CHARACTER (LEN=256) :: COFDIR = '.'


! Output xml files directory

  CHARACTER (LEN=256) :: COXMLDIR = '.'

! IO server will write

  LOGICAL :: LIO_SERV_WR = .FALSE.

! IO server will read

  LOGICAL :: LIO_SERV_RD = .FALSE.

! Use V-set in output file names

  LOGICAL :: LIO_SERV_UVSIF = .TRUE.

! For FP_SERV

  TYPE (FP_SERV_DINF) :: YGPFSDINF

! Time offset for output files 
! For instance, ZTDEC=3600 implies that the file written at t=+1h
! will be named ICMSHFCST+0000

  REAL (KIND=JPRB) :: ZTDEC = 0._JPRB

END TYPE IO_SERV


! IO server used for forecast
TYPE(IO_SERV) :: IO_SERV_C001

SAVE

END MODULE YOMIO_SERV
