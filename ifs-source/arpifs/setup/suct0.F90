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

#ifdef VPP
!OCL SCALAR
#endif
SUBROUTINE SUCT0(KULOUT)

!**** *SUCT0*   - Routine to initialize level 0 control common

!     Purpose.
!     --------
!           Initialize level 0 control commons
!**   Interface.
!     ----------
!        *CALL* *SUCT0(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by SU0YOMA.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   23-Oct-2008 R. El Khatib CFPATH now a directory ending with a '/'
!   K. Yessad (Sep 2008): prune conf 951.
!   K. Yessad (Sep 2008): LREGETA -> LREGETA+LVFE_REGETA.
!   07-Aug-2009 A. Alias NSFXHISTS/NFRSFXHIS replace NSHISTS/NFSRHIS
!                        and LCALLSFX introduced
!   K. Yessad (Aug 2009): prune conf 912, externalise conf 911.
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Aug 2009): remove LSITRIC option
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!   A. Alias  (Feb 2011): LSFXLSM added
!   R. El Khatib 10-Aug-2011 LIOLEVG management
!   D. Degrauwe  (Feb 2012): LARPEGEF_WRGP_HIST added
!   K. Yessad (Feb 2012): add LR3D, LR2D, LRSHW, LRVEQ.
!   R. El Khatib : 01-Mar-2012 (LFPOS,LFPSPEC) => NFPOS ; LECFPOS added
!   R. El Khatib : 26-Jul-2012 Do not overwrite NFPOS if (LFPART2)
!   P. Marguinaud: 11-Sep-2012 More namelist parameters
!   P. Marguinaud: 15-May-2013 Remove LARPEGEF_WRGP_HIST
!   P. Marguinaud: 10-Oct-2013 Initialize namelist parameters NHISTSMIN, NSFXHISTSMIN, NPOSTSMIN
!   F. Vana  09-Jan-2014 LSLPHY available for confs 131 and 501
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!   K. Yessad (Oct 2013): NESCV for vertical displacement.
!   R. El Khatib : 17-Apr-2014 Remove LFPART2
!   K. Yessad (July 2014): Various setup and module refactoring.
!   R. El Khatib 02-Oct-2014 NOPT_MEMORY
!   P. Marguinaud : 10-Oct-2014 Add LGRIB_API
!   R. El Khatib : 03-Dec-2014 skeleton of the configuration 903
!   R. El Khatib : 06-Mar-2015 script name for pp server and real time LAM coupling in namelist
!   R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!   R. El Khatib 07-Mar-2016 Pruning of ISP
!   R. El Khatib : 23-Aug-2016 more interoperability
!   O. Marsden   :    Oct 2016 added L_OOPS(YOMCT0) setup
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   S. Saarinen : 26-Oct-2017 NOPT_MEMORY=2 to call better performing HEAP2-routines
!   P. Lopez     : 21-May-2018  Version for running multiple adjoint tests in one go.
!   R. El Khatib  03-Sep-2018 new configuration 904 which is a test program for change of resolution of an object FIELDS
!   H. Petithomme Sept 2019: add nml key lcorwat for water conservation
!   I. Polichtchouk (Jul 2021): Introduce LSACC model.
! End Modifications
!      ----------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

USE YOMARG   , ONLY : NUCONF=>NCONF      ,NUECMWF=>NECMWF, &
 &                    CUNMEXP=>CNMEXP    ,LUELAM=>LELAM  , &
 &                    UGEMU              ,NSUPERSEDE, NGRIBFILE, NFPSERVER
USE YOMCT0   , ONLY : JPNPST   ,NPOSTS   ,NHISTS   ,NGDITS   , &
 &                    NSDITS   ,NDHFGTS  ,NDHFZTS  ,NDHFDTS  ,NDHPTS   , &
 &                    NMASSCONS,CNMEXP   ,CFPNCF      , &
 &                    NOPT_MEMORY, &
 &                    CNPPATH  ,CFDIRLST ,CNDISPP  ,NCONF    ,NTASKS_CANARI, &
 &                    NQUAD    ,N2DINI   , &
 &                    NSPPR    ,NFRPOS   ,NFRHIS   ,N3DINI   , &
 &                    NFRGDI   ,NFRSDI   ,NCNTVAR  ,NFRDHFG  ,NFRDHFZ  , &
 &                    NFRDHFD  ,NFRDHP   ,NFRCO    , &
 &                    NFRMASSCON,N6BINS  ,LNF      ,LFDBOP   , &
 &                    LARPEGEF ,LSMSSIG  ,CMETER   ,CEVENT   , &
 &                    LARPEGEF_TRAJHR    ,LARPEGEF_TRAJBG,&
 &                    NFPOS    ,LOPDIS   , &
 &                    LCANARI  ,LOLDPP   ,LGUESS   ,LOBS     , &
 &                    LOBSC1   ,LOBSREF  ,LSIMOB   ,LELAM    ,LRPLANE  , &
 &                    LFBDAP   ,LBACKG   ,LMINIM   ,LSCREEN  ,LREFOUT  , &
 &                    LMONITORING, CSCRIPT_PPSERVER,CSCRIPT_LAMRTC     , &
 &                    L_SCREEN_CALL, L_SPLIT_SCREEN,LINFLAT  ,LINFLP9  , &
 &                    LINFLF1  ,LALLOPR  ,LREFGEN  , &
 &                    LIFSTRAJ ,LIFSMIN  , &
 &                    LAROME   ,LECFPOS  , &
 &                    NINTERPTRAJ, NINTERPINCR     ,LSCREEN_OPENMP,&
 &                    NFRCORM  ,LSFORC   ,LSFORCS  ,LSPBSBAL , &
 &                    NSFXHISTS,NFRSFXHIS, &
 &                    LARPEGEF_RDGP_INIT,LARPEGEF_RDGP_TRAJHR,LARPEGEF_RDGP_TRAJBG,&
 &                    LCALLSFX,LSFXLSM,LR3D,LR2D,LRSHW,LRVEQ,LIOLEVG,LECMWF,&
 &                    LWRSPECA_GP,LSUSPECA_GP,LWRSPECA_GP_UV,LSUSPECA_GP_UV,&
 &                    NHISTSMIN, NSFXHISTSMIN, NPOSTSMIN, LGRIB_API,&
 &                    NUNDEFLD, LCOUPLO4, L4DVAR, L_OOPS, &
 &                    NITER_ADTEST,LCONSERV,LCORWAT
USE YOMIOPNH , ONLY : LTRAJNH
USE YOMLUN   , ONLY : NULNAM   
USE YOMMP0   , ONLY : NPROC, NPRGPEW, NPRINTLEV
USE ALGORITHM_STATE_MOD  , ONLY : SET_OBS_IN_FC

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: J,ICONF  !,INSTOP !! suct0 can't know about time yet
INTEGER(KIND=JPIM) :: IERR
LOGICAL :: LL_CONF_EXISTS,LL_CONF_ONEPROCONLY
!!REAL(KIND=JPRB) :: ZUNIT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namct0.nam.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUCT0',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

!     1.1 General 0-level control variables, control of outputs, assimilation

IF (NCONF < 1) NCONF=1

! called from OOPS
L_OOPS=.FALSE.

! Forcing
LSFORC=.FALSE.
LSFORCS=.FALSE.

! Model
LELAM=LUELAM
IF (NSUPERSEDE==1) THEN
  IF (LELAM) THEN
    LRPLANE=UGEMU(2) >= 0.0_JPRB
  ELSE
    LRPLANE=.FALSE.
  ENDIF
ELSE
  LRPLANE=.FALSE.
ENDIF

! Version and name of experiment 
IF(NUECMWF == 1) THEN
  LECMWF=.TRUE.
ELSEIF (NUECMWF == 2) THEN
  LECMWF=.FALSE.
ENDIF
CNMEXP=CUNMEXP

! AROME key
LAROME=.FALSE.

! SURFEX key
LCALLSFX=.TRUE.
LSFXLSM=.FALSE.

! Quadrature
NQUAD =1

! Interpolation methods applied to trajectory and increments
NINTERPTRAJ =1
NINTERPINCR =1

! Configuration
IF (NUCONF > 0) THEN
  NCONF=NUCONF
ENDIF

IF(NUCONF == 302) THEN
  LIFSTRAJ=.TRUE.
  LOBSC1=.TRUE.
ENDIF

! Number of time steps
!IF(NSUPERSEDE > 0 .AND. NUSTOP >= 0) THEN
!  INSTOP=NUSTOP
!ELSE
!  INSTOP=2
!ENDIF


! Ocean coupling
LCOUPLO4=.FALSE.

! Other variables may depend of the previous ones and the command line arg.

LALLOPR  =.FALSE.
NFRGDI=1
NFRHIS=1
NFRPOS=1
NFRSDI=1
!!NFRMASSCON=INSTOP+999999  !! this is a silly initialisation, right? olivier
NFRMASSCON=999999
NMASSCONS(:)=0

LSCREEN =.FALSE.
LSCREEN_OPENMP =.TRUE.
L_SPLIT_SCREEN=.FALSE.
L_SCREEN_CALL=.TRUE.
LMONITORING=.FALSE.
LREFOUT =.FALSE.
LREFGEN =.FALSE.
LSPBSBAL=.FALSE.
LOLDPP  =.FALSE.
LOPDIS=.TRUE.
NTASKS_CANARI= 1
#ifdef __INTEL_COMPILER
NOPT_MEMORY=2
#elif defined(__GFORTRAN__)
NOPT_MEMORY=2
#else
NOPT_MEMORY=2
#endif
LINFLAT=.FALSE.
LINFLP9=.FALSE.
LINFLF1=.FALSE.
LIOLEVG=.TRUE.

NHISTSMIN=0
NSFXHISTSMIN=0
NPOSTSMIN=0

NITER_ADTEST=1

IF (LECMWF) THEN
  LOLDPP  =.TRUE.
  LTRAJNH=.FALSE.
  LOBSC1=.FALSE.
  LOBS  =.FALSE.
  NFRSFXHIS=1
  NSFXHISTS(:)=0
  IF(NUCONF == 2) THEN
    LOBSC1=.TRUE.
    NUCONF=1
    NCONF=1
    LIFSTRAJ=.TRUE.
  ELSE
    LIFSTRAJ=.FALSE.
  ENDIF
  IF((NCONF == 1.OR.NCONF == 302).AND.LOBSC1) THEN
    LOBS=.TRUE.
    NFRHIS=1000
    NHISTS(0)=1
    NHISTS(1)=1
  ENDIF
  IF(NCONF/100 == 1) THEN
    LOBS=.TRUE.
  ENDIF
  IF(NCONF /= 1 .AND. NCONF /= 131  .AND. NCONF /= 302 .AND. NCONF /= 401) THEN
    LOPDIS=.FALSE.
  ENDIF
  LSMSSIG =.FALSE.
  CMETER='smsmeter'
  CEVENT='smsevent'
  LNF   =.TRUE.
  N2DINI=1
  N3DINI=0
  NSPPR =0
!!  NFRPOS=INSTOP
!!  NFRHIS=INSTOP   !! olivier
  NFRPOS = -1       !! initialise to -1; if not written over by namelist, correct default value will be provided in suarg_datetime
  NFRCO =0
  NFRGDI=1
  NFRSDI=24
  NFRDHFG=4
  NFRDHFZ=4
  NFRDHFD=4
  NFRDHP=48
  N6BINS=6
  DO J=0,JPNPST
    NPOSTS(J)=0
    NHISTS(J)=0
    NGDITS(J)=0
    NSDITS(J)=0
    NDHFGTS(J)=0
    NDHFZTS(J)=0
    NDHFDTS(J)=0
    NDHPTS(J)=0
  ENDDO
  LCANARI=.FALSE.
  IF (NCONF==131) THEN
    LGUESS=.TRUE.
  ELSE
    LGUESS=.FALSE.
  ENDIF
  IF(NCONF/100 == 1)THEN
    LBACKG=.TRUE.
  ELSE
    LBACKG=.FALSE.
  ENDIF
  IF ((NCONF/100 == 1).OR.(NCONF/100 == 8).OR.(NCONF == 923))THEN
    LMINIM=.TRUE.
  ELSE
    LMINIM=.FALSE.
  ENDIF
  LOBSREF=.FALSE.
  LSIMOB=.FALSE.
  NCNTVAR = 2
ELSE
  ! ky: set provisional value to avoid unitialised values;
  !     some of them will be modified below according to conf.
  DO J=0,JPNPST
    NPOSTS(J)=0
    NHISTS(J)=0
    NGDITS(J)=0
    NSDITS(J)=0
    NDHFGTS(J)=0
    NDHFZTS(J)=0
    NDHFDTS(J)=0
    NDHPTS(J)=0
    NSFXHISTS(J)=0
  ENDDO

  IF(NUCONF == 2) THEN
    LOBSC1=.TRUE.
    NUCONF=1
    NCONF=1
    LOBS=.TRUE.
    NFRHIS=1000
    NHISTS(0)=1
    NHISTS(1)=1
  ELSE
    LOBSC1=.FALSE.
    LOBS=.FALSE.
  ENDIF
  LTRAJNH=.FALSE.
  NFRSFXHIS=1
  NSFXHISTS(:)=0
  IF (NSUPERSEDE==1) THEN
    IF (NCONF == 1 .OR. NCONF == 302) THEN
      NFRHIS=1
      NHISTS(0)=-10
      NHISTS(1)=-0
      NHISTS(2)=-6
      NHISTS(3)=-12
      NHISTS(4)=-18
      NHISTS(5)=-24
      NHISTS(6)=-30
      NHISTS(7)=-36
      NHISTS(8)=-48
      NHISTS(9)=-60
      NHISTS(10)=-72
      NFRPOS=1
      NPOSTS(0)=-10
      NPOSTS(1)=-0
      NPOSTS(2)=-6
      NPOSTS(3)=-12
      NPOSTS(4)=-18
      NPOSTS(5)=-24
      NPOSTS(6)=-30
      NPOSTS(7)=-36
      NPOSTS(8)=-48
      NPOSTS(9)=-60
      NPOSTS(10)=-72
      NFRSDI=20
      DO J=11,JPNPST
        NPOSTS(J)=0
        NHISTS(J)=0
      ENDDO
      DO J=0,JPNPST
        NGDITS(J)=0
        NSDITS(J)=0
        NDHFGTS(J)=0
        NDHFZTS(J)=0
        NDHFDTS(J)=0
        NDHPTS(J)=0
      ENDDO
    ELSEIF (NCONF == 701) THEN
      NFRHIS=1
      NFRPOS=4
      NFRSDI=1
      NHISTS(0)=-1
      NHISTS(1)=-0
      DO J=0,JPNPST
        NPOSTS(J)=0
        NGDITS(J)=0
        NSDITS(J)=0
        NDHFGTS(J)=0
        NDHFZTS(J)=0
        NDHFDTS(J)=0
        NDHPTS(J)=0
      ENDDO
      DO J=2,JPNPST
        NHISTS(J)=0
      ENDDO
    ELSEIF (NCONF == 131) THEN
      NFRHIS=1
      NFRPOS=4
      NFRSDI=1
      DO J=0,JPNPST
        NPOSTS(J)=0
        NHISTS(J)=0
        NGDITS(J)=0
        NSDITS(J)=0
        NDHFGTS(J)=0
        NDHFZTS(J)=0
        NDHFDTS(J)=0
        NDHPTS(J)=0
      ENDDO
    ENDIF
  ELSE
    NFRHIS=1
    NFRPOS=4
    NFRSDI=1
    DO J=0,JPNPST
      NPOSTS(J)=0
      NHISTS(J)=0
      NGDITS(J)=0
      NSDITS(J)=0
      NDHFGTS(J)=0
      NDHFZTS(J)=0
      NDHFDTS(J)=0
      NDHPTS(J)=0
    ENDDO
  ENDIF
  N6BINS=0
  IF (NCONF == 701) THEN
    LCANARI=.TRUE.
    LSIMOB=.FALSE.
    LOBS=.TRUE.
  ELSE
    LCANARI=.FALSE.
    LSIMOB=.TRUE.
    LOBS=.FALSE.
  ENDIF
  LSMSSIG =.FALSE.
  LNF   =.TRUE.
  N2DINI=1
  N3DINI=0
  NSPPR =0
  NFRCO =0
  NFRGDI=1
  NFRDHFG=4
  NFRDHFZ=4
  NFRDHFD=4
  NFRDHP=48
  LGUESS=.TRUE.
  IF(NCONF/100 == 1)THEN
    LBACKG=.TRUE.
  ELSE
    LBACKG=.FALSE.
  ENDIF
  IF ((NCONF/100 == 1).OR.(NCONF/100 == 8).OR.(NCONF == 923))THEN
    LMINIM=.TRUE.
  ELSE
    LMINIM=.FALSE.
  ENDIF
  LOBSREF=.FALSE.
  NCNTVAR = 2
ENDIF

NFRCORM=0
LCONSERV=.FALSE.
LCORWAT=.FALSE.
L4DVAR=.FALSE.

LGRIB_API = .NOT. LELAM

! Set default for index to be used for unused/undefined fields
NUNDEFLD=-99999999

!     1.2 File kind or content

IF (LECMWF) THEN
  LFDBOP=.TRUE.
  LFBDAP=.FALSE.
  LARPEGEF=.FALSE.
  LARPEGEF_TRAJHR=.FALSE.
  LARPEGEF_TRAJBG=.FALSE.
  LARPEGEF_RDGP_INIT=.FALSE.
  LARPEGEF_RDGP_TRAJHR=.FALSE.
  LARPEGEF_RDGP_TRAJBG=.FALSE.
ELSE
  LFDBOP=.FALSE.
  LFBDAP=.TRUE.
  LARPEGEF=.TRUE.
  LARPEGEF_TRAJHR=.TRUE.
  LARPEGEF_TRAJBG=.TRUE.
  LARPEGEF_RDGP_INIT=.FALSE.
  LARPEGEF_RDGP_TRAJHR=.FALSE.
  LARPEGEF_RDGP_TRAJBG=.FALSE.
ENDIF
CNDISPP=' '
LWRSPECA_GP=.FALSE.
LSUSPECA_GP=.FALSE.
LWRSPECA_GP_UV=.FALSE.
LSUSPECA_GP_UV=.FALSE.

!     1.3 General FULL-POS keys

NFPOS=0
CFPNCF='ncf927'
CFDIRLST=' '
CNPPATH=' '

!     1.4 Script files

IF (NFPSERVER == 0) THEN
  CSCRIPT_PPSERVER='cnt3_wait'
ELSE
  CSCRIPT_PPSERVER=' '
ENDIF
CSCRIPT_LAMRTC='atcp'

!      ----------------------------------------------------------------

!*       2.    READ NAMELIST.
!              --------------

CALL POSNAM(NULNAM,'NAMCT0')
READ(NULNAM,NAMCT0)

!      ----------------------------------------------------------------

!*       3.    RESET VALUES AND TEST.
!              ----------------------


IERR=0

!     3.1 General 0-level control variables, control of outputs, assimilation

!!**IF (NSUPERSEDE > 0) THEN
!!**  ZUNIT=3600._JPRB
!!**  IF ((TSTEP > 0.0_JPRB).AND.(NFRPOS < 0)) THEN
!!**    NFRPOS=NINT((REAL(-NFRPOS,JPRB)*ZUNIT)/TSTEP)
!!**  ENDIF
!!**  IF ((TSTEP > 0.0_JPRB).AND.(NFRHIS < 0)) THEN
!!**    NFRHIS=NINT((REAL(-NFRHIS,JPRB)*ZUNIT)/TSTEP)
!!**  ENDIF
!!**  IF ((TSTEP > 0.0_JPRB).AND.(NFRSFXHIS < 0)) THEN
!!**    NFRSFXHIS=NINT((REAL(-NFRSFXHIS,JPRB)*ZUNIT)/TSTEP)
!!**  ENDIF
!!**  IF ((TSTEP > 0.0_JPRB).AND.(NFRSDI < 0)) THEN
!!**    NFRSDI=NINT((REAL(-NFRSDI,JPRB)*ZUNIT)/TSTEP)
!!**  ENDIF
!!**ENDIF

ICONF=NCONF/100

IF ( NCONF /= 131 .AND. NCONF /= 701 ) THEN
  LSIMOB=.TRUE.
ENDIF

IF (NCONF /= 1 .AND. NCONF /= 302) LOBSC1=.FALSE.
IF (LOBSC1) LSIMOB=.FALSE.
IF (LSIMOB) NCNTVAR=1
IF (NCONF == 801 .AND. LBACKG) NCNTVAR=2
IF (NPRINTLEV > 0) LOPDIS = .TRUE.

IF (NCONF == 131) THEN
  LOBSREF = .TRUE.
ELSE
  LOBSREF = .FALSE.
ENDIF

IF (NCONF == 903) THEN
  NFPOS=MAX(NFPOS,1)
ENDIF

IF (NCONF == 201.OR.NCONF == 421.OR.NCONF == 521) THEN
  LRSHW=.TRUE.
  LRVEQ=.FALSE.
  LR3D=.FALSE.
ELSEIF (NCONF == 202.OR.NCONF == 422.OR.NCONF == 522) THEN
  LRSHW=.FALSE.
  LRVEQ=.TRUE.
  LR3D=.FALSE.
ELSEIF (NCONF==701 .OR. NCONF==901 .OR. NCONF==923 .OR. NCONF==931 .OR. NCONF==932 .OR. NCONF==933) THEN
  LRSHW=.FALSE.
  LRVEQ=.FALSE.
  LR3D=.FALSE.
ELSE
  LRSHW=.FALSE.
  LRVEQ=.FALSE.
  LR3D=.TRUE.
ENDIF
LR2D=LRSHW.OR.LRVEQ

LL_CONF_EXISTS=(NCONF == 1) &
 & .OR.(NCONF == 131) &
 & .OR.(NCONF == 201).OR.(NCONF == 202) &
 & .OR.(NCONF == 302) &
 & .OR.(NCONF == 401).OR.(NCONF == 421).OR.(NCONF == 422) &
 & .OR.(NCONF == 501).OR.(NCONF == 521).OR.(NCONF == 522) &
 & .OR.(NCONF == 601) &
 & .OR.(NCONF == 701) &
 & .OR.(NCONF == 801) &
 & .OR.(NCONF == 901) &
 & .OR.(NCONF == 903) &
 & .OR.(NCONF == 904) &
 & .OR.(NCONF == 923) &
 & .OR.(NCONF == 931).OR.(NCONF == 932).OR.(NCONF == 933)
IF (.NOT.LL_CONF_EXISTS) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUCT0: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*) ' NCONF=',NCONF,' does not exist.'
ENDIF



! * Ask for several processors for configurations coded for one processor only?
!   (quid about conf 931,932,933 ??)
LL_CONF_ONEPROCONLY=(NCONF == 901).OR.(NCONF == 923)
IF (NPROC /= 1 .AND. (LL_CONF_ONEPROCONLY)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUCT0: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*)' Must run conf number',NCONF,' with one proc!'
ENDIF

! * Ask for E-W distribution for configurations which do not support it?
IF (LR2D.AND.NPRGPEW > 1) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUCT0: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*) ' NPRGPEW MUST BE 1 FOR 2D-MODEL CONFIGURATIONS!'
ENDIF


ICONF=NCONF/100
LIFSMIN=ANY(SPREAD(ICONF,1,5) ==(/1,4,5,6,8/))
!LIFSMIN=ANY(SPREAD(ICONF,1,6) ==(/1,3,4,5,6,8/))

IF (.NOT.LECMWF) THEN
  ! remark KY+GD: this seems necessary for use at METEO-FRANCE.
  LIFSTRAJ=(NCONF==1)
  WRITE (KULOUT,*)' SUCT0 : LIFSTRAJ=(NCONF==1) ',LIFSTRAJ
ENDIF


! * Reset L4DVAR:
IF (NCONF /= 131) L4DVAR=.FALSE.


!     3.2 File kind or content

IF (NGRIBFILE == 1) THEN
  LARPEGEF_TRAJHR=.FALSE.
  LARPEGEF_TRAJBG=.FALSE.
  LARPEGEF_RDGP_INIT=.FALSE.
  LARPEGEF_RDGP_TRAJHR=.FALSE.
  LARPEGEF_RDGP_TRAJBG=.FALSE.
ENDIF  

! * Use LARPEGEF_RDGP_TRAJHR=T when LARPEGEF_TRAJHR=F?
IF (LARPEGEF_RDGP_TRAJHR.AND.(.NOT.LARPEGEF_TRAJHR)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUCT0: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*) ' If LARPEGEF_RDGP_TRAJHR=.TRUE., LARPEGEF_TRAJHR must be set to .TRUE.'
ENDIF

! * Use LARPEGEF_RDGP_TRAJBG=T when LARPEGEF_TRAJBG=F?
IF (LARPEGEF_RDGP_TRAJBG.AND.(.NOT.LARPEGEF_TRAJBG)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUCT0: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*) ' If LARPEGEF_RDGP_TRAJBG=.TRUE., LARPEGEF_TRAJBG must be set to .TRUE.'
ENDIF

! * Restricted usage of Fullpos by ECMWF :
LECFPOS=LECMWF.AND.(NFPOS /= 0).AND..NOT.LARPEGEF

IF (IERR >= 1) THEN
  CALL FLUSH(KULOUT)
  CALL ABOR1(' SUCT0: ABOR1 CALLED')
ENDIF

CALL SET_OBS_IN_FC(LOBSC1) ! Make sure L_OBS_IN_FC_CONFIG is consistent with LOBSC1 (OOPS)

!      -----------------------------------------------------------

!*       4.    PRINT FINAL VALUES.
!              ------------------

WRITE(UNIT=KULOUT,FMT='('' Printings for module YOMCT0 '')')

WRITE(KULOUT,*) '* SUCT0: General 0-level control variables.'
WRITE(UNIT=KULOUT,FMT='('' LECMWF = '',L2)') LECMWF
WRITE(UNIT=KULOUT,FMT='('' LELAM  = '',L2,'' LRPLANE  = '',L2)') LELAM,LRPLANE  
WRITE(UNIT=KULOUT,FMT='('' NCONF  = '',I6)') NCONF
WRITE(UNIT=KULOUT,FMT='('' CNMEXP = '',A16)') CNMEXP
WRITE(UNIT=KULOUT,FMT='('' LNF    = '',L2)') LNF
WRITE(UNIT=KULOUT,FMT='('' LTRAJNH = '',L2)') LTRAJNH
WRITE(UNIT=KULOUT,FMT='('' LSFORC = '',L2)') LSFORC
WRITE(UNIT=KULOUT,FMT='('' NQUAD  = '',I6)') NQUAD
WRITE(UNIT=KULOUT,FMT='('' N2DINI = '',I6)') N2DINI
WRITE(UNIT=KULOUT,FMT='('' N3DINI = '',I6)') N3DINI
WRITE(UNIT=KULOUT,FMT='('' NUNDEFLD = '',I10)') NUNDEFLD
WRITE(UNIT=KULOUT,FMT='('' NOPT_MEMORY = '',I2)') NOPT_MEMORY
WRITE(UNIT=KULOUT,FMT='('' NITER_ADTEST = '',I3)') NITER_ADTEST

WRITE(UNIT=KULOUT,FMT='('' LIOLEVG = '',L2)') LIOLEVG
WRITE(UNIT=KULOUT,FMT='('' LALLOPR = '',L2)') LALLOPR
WRITE(UNIT=KULOUT,FMT='('' LOLDPP ='',L2)') LOLDPP
WRITE(UNIT=KULOUT,FMT='('' LAROME ='',L2)') LAROME
WRITE(UNIT=KULOUT,FMT='('' LCALLSFX ='',L2,'' LSFXLSM ='',L2)') LCALLSFX,LSFXLSM  
WRITE(UNIT=KULOUT,FMT='('' LSFORC ='',L2,'' LSFORCS ='',L2)') LSFORC,LSFORCS
WRITE(UNIT=KULOUT,FMT='('' LCOUPLO4 ='',L2)') LCOUPLO4

WRITE(KULOUT,*) '* SUCT0: Assimilation general keys.'
WRITE(UNIT=KULOUT,FMT='('' L4DVAR = '',L2)') L4DVAR
WRITE(UNIT=KULOUT,FMT='('' NTASKS_CANARI = '',I6)') NTASKS_CANARI
WRITE(UNIT=KULOUT,FMT='('' NINTERPTRAJ = '',I2, &
 & '' NINTERPINCR = '',I2)') NINTERPTRAJ, NINTERPINCR
WRITE(UNIT=KULOUT,FMT='('' LCANARI= '',L2 &
 & ,'' LGUESS = '',L2,'' LOBSC1 = '',L2 &
 & ,'' LOBSREF= '',L2,'' LOBS = '',L2 &
 & ,'' LSIMOB = '',L2 )')& 
 & LCANARI,LGUESS,LOBSC1,LOBSREF,LOBS,LSIMOB
WRITE(UNIT=KULOUT,FMT='('' LSCREEN = '',L2,'' L_SCREEN_CALL = '',L2, &
 & '' L_SPLIT_SCREEN = '',L2,'' LBACKG = '',L2,'' LMINIM = '',L2)') &
 & LSCREEN,L_SCREEN_CALL,L_SPLIT_SCREEN,LBACKG,LMINIM
WRITE(UNIT=KULOUT,FMT='('' LSCREEN_OPENMP = '',L2)') LSCREEN_OPENMP
IF (NCONF/100 == 1.OR.NCONF == 801) THEN
  IF (NCNTVAR == 1) THEN
    WRITE(KULOUT,'('' NCNTVAR = 1: CONTROL VARIABLES ARE'',&
     & '' MODEL VARIABLES'')')
  ELSEIF (NCNTVAR == 2) THEN
    WRITE(KULOUT,'('' NCNTVAR = 2: CONTROL VARIABLES ARE'',&
     & '' DEPARTURES FROM THE FIRST GUESS '',/,&
     & '' NORMALIZED BY FORECAST ERROR '',&
     & '' STANDARD DEVIATIONS'',/,&
     & '' PROJECTED ONTO THE EIGENVALUES OF THE '',&
     & '' VERTICAL PREDICTION ERROR CORRELATION'',&
     & '' MATRICES'')')
  ELSE
    WRITE(KULOUT,'('' NCNTVAR = '',I5)') NCNTVAR
  ENDIF
ENDIF
WRITE(UNIT=KULOUT,FMT='('' LIFSTRAJ ='',L2)') LIFSTRAJ
WRITE(UNIT=KULOUT,FMT='('' LIFSMIN ='',L2)') LIFSMIN
WRITE(UNIT=KULOUT,FMT='('' LMONITORING = '',L2)') LMONITORING
WRITE(UNIT=KULOUT,FMT='('' LINFLAT ='',L2)') LINFLAT
WRITE(UNIT=KULOUT,FMT='('' LINFLP9 ='',L2)') LINFLP9
WRITE(UNIT=KULOUT,FMT='('' LINFLF1 ='',L2)') LINFLF1
WRITE(UNIT=KULOUT,FMT='('' LSPBSBAL ='',L2)') LSPBSBAL

WRITE(KULOUT,*) '* SUCT0: Control of outputs.'
WRITE(UNIT=KULOUT,FMT='('' NSPPR  = '',I6)') NSPPR  
WRITE(UNIT=KULOUT,FMT='('' NFRPOS = '',I6,'' NFRSFXHIS = '',I6 &
 & ,'' NFRHIS = '',I6,'' NFRGDI = '',I6,'' NFRSDI = '',I6)')&
 & NFRPOS,NFRSFXHIS,NFRHIS,NFRGDI,NFRSDI  
WRITE(UNIT=KULOUT,FMT='('' NFRDHFG = '',I6,'' NFRDHFZ = '',I6,&
 & '' NFRDHFD = '',I6,'' NFRDHP = '',I6,'' N6BINS = '',I6)')&
 & NFRDHFG,NFRDHFZ,NFRDHFD,NFRDHP,N6BINS  
WRITE(UNIT=KULOUT,FMT='('' NFRCO = '',I6)') NFRCO
WRITE(UNIT=KULOUT,FMT='('' NFRCORM = '',I6)') NFRCORM
WRITE(UNIT=KULOUT,FMT="(' LCONSERV/LCORWAT: ',2L2)") LCONSERV,LCORWAT
WRITE(KULOUT,*) ' NPOSTS =  ',NPOSTS(0),(NPOSTS(J),J=1,ABS(NPOSTS(0)))
WRITE(KULOUT,*) ' NPOSTSMIN =  ',NPOSTSMIN(0),(NPOSTSMIN(J),J=1,ABS(NPOSTSMIN(0)))
WRITE(KULOUT,*) ' NHISTS =  ',NHISTS(0),(NHISTS(J),J=1,ABS(NHISTS(0)))
WRITE(KULOUT,*) ' NHISTSMIN =  ',NHISTSMIN(0),(NHISTSMIN(J),J=1,ABS(NHISTSMIN(0)))
WRITE(KULOUT,*) ' NSFXHISTS = ',NSFXHISTS(0),(NSFXHISTS(J),J=1,ABS(NSFXHISTS(0)))
WRITE(KULOUT,*) ' NSFXHISTSMIN =  ',NSFXHISTSMIN(0),(NSFXHISTSMIN(J),J=1,ABS(NSFXHISTSMIN(0)))
WRITE(KULOUT,*) ' NGDITS =  ',NGDITS(0),(NGDITS(J),J=1,ABS(NGDITS(0)))
WRITE(KULOUT,*) ' NSDITS =  ',NSDITS(0),(NSDITS(J),J=1,ABS(NSDITS(0)))
WRITE(KULOUT,*) ' NDHPTS =  ',NDHPTS(0),(NDHPTS(J),J=1,ABS(NDHPTS(0)))
WRITE(KULOUT,*) ' NDHFGTS = ',NDHFGTS(0),(NDHFGTS(J),J=1,ABS(NDHFGTS(0)))
WRITE(KULOUT,*) ' NDHFZTS = ',NDHFZTS(0),(NDHFZTS(J),J=1,ABS(NDHFZTS(0)))
WRITE(KULOUT,*) ' NDHFDTS = ',NDHFDTS(0),(NDHFDTS(J),J=1,ABS(NDHFDTS(0)))
WRITE(KULOUT,*) ' NMASSCONS = ',NMASSCONS(0),(NMASSCONS(J),J=1,ABS(NMASSCONS(0)))  

WRITE(KULOUT,*) '* SUCT0: File kind or content.'
WRITE(UNIT=KULOUT,FMT='('' LFDBOP = '',L2,'' LFBDAP = '',L2)') LFDBOP,LFBDAP  
WRITE(UNIT=KULOUT,FMT='('' LGRIB_API = '',L2)') LGRIB_API
WRITE(UNIT=KULOUT,FMT='('' LARPEGEF = '',L2)') LARPEGEF
WRITE(UNIT=KULOUT,FMT='('' LARPEGEF_TRAJHR = '',L2)') LARPEGEF_TRAJHR
WRITE(UNIT=KULOUT,FMT='('' LARPEGEF_TRAJBG = '',L2)') LARPEGEF_TRAJBG
WRITE(UNIT=KULOUT,FMT='('' LARPEGEF_RDGP_INIT = '',L2)') LARPEGEF_RDGP_INIT
WRITE(UNIT=KULOUT,FMT='('' LARPEGEF_RDGP_TRAJHR = '',L2)') LARPEGEF_RDGP_TRAJHR
WRITE(UNIT=KULOUT,FMT='('' LARPEGEF_RDGP_TRAJBG = '',L2)') LARPEGEF_RDGP_TRAJBG
WRITE(UNIT=KULOUT,FMT='('' LWRSPECA_GP = '',L2)') LWRSPECA_GP
WRITE(UNIT=KULOUT,FMT='('' LSUSPECA_GP = '',L2)') LSUSPECA_GP
WRITE(UNIT=KULOUT,FMT='('' LWRSPECA_GP_UV = '',L2)') LWRSPECA_GP_UV
WRITE(UNIT=KULOUT,FMT='('' LSUSPECA_GP_UV = '',L2)') LSUSPECA_GP_UV

WRITE(KULOUT,*) '* SUCT0: General FULL-POS keys.'
WRITE(UNIT=KULOUT,FMT='('' NFPOS = '',I3, '' LECFPOS = '',L2)') NFPOS, LECFPOS

WRITE(KULOUT,*) '* SUCT0: Script filenames.'
WRITE(UNIT=KULOUT,FMT='('' CSCRIPT_LAMRTC  = '',A,'' CSCRIPT_PPSERVER = '',A)') TRIM(CSCRIPT_LAMRTC), TRIM(CSCRIPT_PPSERVER)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUCT0',1,ZHOOK_HANDLE)
END SUBROUTINE SUCT0
