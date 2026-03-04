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

SUBROUTINE SUARG (LDFAINIT)

!**** *SUARG*   - Routine to initialize common containing command line argument

!     Purpose.
!     --------
!           Initialize common YOMARG
!**   Interface.
!     ----------
!        *CALL* *SUARG

!        Explicit arguments :
!        --------------------
!        None

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
!      Ryad El Khatib *METEO-FRANCE*
!      Original : 93-05-06

!     Modifications.
!     --------------
!      Modified : 01-02-06 by S.Martinez : move GETOPT to GET_OPT (!!!)
!      Modified : 01-03-20 Ryad El Khatib : "927 without command line options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      G.Mozdzynski  03-12-12 Broadcast commmand line arguments
!      Modified : 05-03-01 Ryad El Khatib : Cleanups
!      M.Hamrud      01-May-2006 Generalized IO scheme
!      R. El Ouaraini & JM. Audoin 03-Oct-06: NVGRIB default value is that of
!                                                            the initial file
!      R. El Khatib  27-Apr-2007 Re-enable the possibility to switch off MPI
!      J.Woyciech/R. El Khatib Dec 09 Bf: call to fanion replaced by call to faveur 
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!        M. Mile  : 20-Okt-2010 GRIB_API for setup args from grib
!      P. Moll      26-Jui-2011 BF for CALL INI_IOSTREAM
!      R. El Khatib 22-Mar-2012 Fix uninitialized variables
!      R. El Khatib 30-Mar-2012 fanmsg to re-direct FA software verbosity towards nulout
!      R. El Khatib 13-Aug-2013 NUCMAX=NSMAX
!      T. Wilhelmsson 1-July-2013 Use GRIB_API_INTERFACE
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Various setup and module refactoring.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El Khatib 03-Dec-2014 NGRIBFILE=1 if input files in GRIB ; =0 if input file in FA
!      R. El Khatib 18-mar-2016 NFPSERVER
!      O. Marsden      May-2016 Replace call to CHIEN by call to RIEN
!      O. Marsden      Nov 2017 Removal of NUSTOP and UTSTEP from YOMARG, replaced by CSTOP and TSTEP in YOMRIP
!     ------------------------------------------------------------------

USE PARKIND1          , ONLY : JPIM     ,JPRB
USE YOMHOOK           , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARDIM            , ONLY : JPMXLE    ,JPMXGL
USE YOMARG            , ONLY : NULOEN   ,NUMEN    ,NULIM    ,NCONF    ,NFPSERVER,&
 &                             NECMWF   ,CNMEXP   ,LELAM    ,UVALH    ,UVBH     ,NGRIBFILE,&
 &                             UGEMU    ,NUCMAX   ,NUDGL    ,NUDLON   ,NUFLEV   ,NUSMAX   ,NUHTYP   ,&
 &                             NUSTTYP  ,NUDATE   ,NUSSSS   ,USTRET   ,UMUCEN   ,ULOCEN   ,NSUPERSEDE,&
 &                             LECMWF
USE YOMLUN            , ONLY : NULOUT   ,NINISH   ,NULNAM
USE YOMCT0            , ONLY : NSTEPINI
USE YOMOPH0           , ONLY : CFNISH, CFNIGG
USE YOMMP0            , ONLY : NPROC, MYPROC
USE IOSTREAM_MIX      , ONLY : INI_IOSTREAM, SETUP_IOSTREAM, SETUP_IOREQUEST,  &
 &                             CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST, IO_INQUIRE, &
 &                             CLOSE_IOREQUEST
USE MPL_MODULE        , ONLY : MPL_ALLGATHERV, MPL_BROADCAST
USE FA_MOD            , ONLY : JD_YEA, JD_MON, JD_DAY, JD_SEM, JD_SET
USE GRIB_API_INTERFACE, ONLY : IGRIB_OPEN_FILE, IGRIB_NEW_FROM_FILE, &
 &                             IGRIB_GET_VALUE, IGRIB_RELEASE, IGRIB_CLOSE_FILE

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, OPTIONAL, INTENT (IN) :: LDFAINIT

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IDATEF(22), ILMO(12)
INTEGER(KIND=JPIM) :: IULIM(JPMXGL/2), IULOEN(JPMXGL),IUMEN(JPMXGL)

CHARACTER (LEN = 50) :: CLFNISH
CHARACTER (LEN = 16) :: CLLEC
CHARACTER (LEN = 14) :: CLFNIGG
CHARACTER (LEN =4)   :: CLSTRING
!!not used!!  CHARACTER (LEN =4 ), PARAMETER :: CDPREF='SURF'
!!not used!!  CHARACTER (LEN =16), PARAMETER :: CDSUFF='PRESSION        '

REAL(KIND=JPRB) :: ZGEMU(JPMXGL), ZVALH(0:JPMXLE), ZVBH (0:JPMXLE)
REAL(KIND=JPRB) :: ZPV(0:JPMXLE*2+1)

INTEGER(KIND=JPIM) :: IABORT, IDGNH, ILEN1, &
 & ILEN2, ILEN4, ILENR2, ILEVMN, ILEVMX, INBARI, INC, IREP, &
 & ISDAY, ISHOUR, ISSSS, ITAG, IUDATE, IUDGL, IUDLON, IUFLEV, &
 & IUNTIN, IUQUAD, IUSMAX, IUSSSS, IUSTTYP, JGL, JLEV,ISTEP  ,&
 & IKNGRIB, INGRIB, IKNBPDG, IKNBCSP, IKSTRON, ISTRON, IKPUILA, IKDMOPL

INTEGER(KIND=JPIM) :: IVERBOSE ! Verbosity of FA software
INTEGER(KIND=JPIM) :: IEXIST,IIEXIST(NPROC)

LOGICAL :: LLEXIST, LLFICP, LLGGEXIST, LLMAP, LLSPEXIST, LLFAINIT

REAL(KIND=JPRB) :: Z_RINC, ZEPS, ZLOCEN, ZMUCEN, ZPREF, ZSTRET
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST

INTEGER(KIND=JPIM) :: I, JGRIB, IFILE, ITYPEOFGRID, IPLPRESENT, IDZONL, ISTROE

!     ------------------------------------------------------------------

#include "rien.h"
#include "erien.h"
#include "sufainit.intfb.h"

#include "abor1.intfb.h"
#include "mod_ini.intfb.h"
#include "updcal.intfb.h"
#include "namarg.nam.h"
#include "posnam.intfb.h"
#include "suoph0.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUARG',0,ZHOOK_HANDLE)

ASSOCIATE(LUELAM => LELAM, NUCONF => NCONF, LUECMWF => LECMWF, &
 &        NUECMWF => NECMWF) 

LLFAINIT = .TRUE.
IF (PRESENT (LDFAINIT)) LLFAINIT = LDFAINIT

!     ------------------------------------------------------------------

!   c: configuration (default=1)
!   v: version of model (no default)
!   e: experiment identifier (no default)
!   s: time step (default -999., it will be set according to resolution
!         and advection scheme)
!   f: forecast span (dxx: days; hxxx: hours; txxx: time-steps)
!   a: advection scheme (eul: Eulerian; sli: interpolating semi-Lag.
!   m: model (arpifs: ARPEGE/IFS ; aladin: ALADIN ; default=arpifs)
!      ----------------------------------------------------------------

!*       1.    DEFAULT VALUES FOR PARAMETERS ARE IN DATA
!              -----------------------------------------


ZEPS=EPSILON(1.0_JPRB)*10000._JPRB
ISHOUR=3600
ISDAY=3600*24
ILEN1=125
ILEN2=22+JPMXGL
ILENR2=10+2*(JPMXLE+1)
ILEN4=1
ILEVMN=0
ILEVMX=JPMXLE
NSTEPINI=0

!        2     Set default values and read namelist
!              --------------------------------

LUELAM=.FALSE.
NUCONF=1
LUECMWF=.TRUE.
CNMEXP='0123'
NSUPERSEDE=1
NFPSERVER=0


CALL POSNAM(NULNAM,'NAMARG')
READ(NULNAM,NAMARG)

!        2.1   Information from NAMARG

IF( LUECMWF )THEN
  NUECMWF=1
ELSE
  NUECMWF=2
ENDIF

!*    Initialize file names
WRITE(NULOUT,*) '---- Set up files names -----'
CALL SUOPH0(CNMEXP)

! Guess what kind of file in input
IF(NSUPERSEDE==1) THEN
  IF (MYPROC==1) THEN
    IUNTIN=NINISH
    CLFNISH=CFNISH
    INQUIRE(FILE=CLFNISH,EXIST=LLSPEXIST)
    IF (LLSPEXIST) THEN
      OPEN(FILE=CLFNISH,UNIT=IUNTIN,FORM='FORMATTED')
      REWIND(IUNTIN)
      READ(UNIT=IUNTIN,FMT='(A4)') CLSTRING
      CLOSE(IUNTIN)
      IF (CLSTRING=='GRIB') THEN
        NGRIBFILE=1
      ELSE
        NGRIBFILE=0
      ENDIF
    ELSE
      CLFNIGG=CFNIGG
      INQUIRE(FILE=CLFNIGG,EXIST=LLGGEXIST)
      IF (LLGGEXIST) THEN
        OPEN(FILE=CLFNIGG,UNIT=IUNTIN,FORM='FORMATTED')
        REWIND(IUNTIN)
        READ(UNIT=IUNTIN,FMT='(A4)') CLSTRING
        CLOSE(IUNTIN)
        IF (CLSTRING=='GRIB') THEN
          NGRIBFILE=1
        ELSE
          NGRIBFILE=0
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  CALL MPL_BROADCAST(NGRIBFILE,KTAG=0,KROOT=1,CDSTRING='SUARG:')
ELSE
  NGRIBFILE=2-NUECMWF
ENDIF

! Defines parameters for FA (unconditional, for interoperability purpose) :
IF (LLFAINIT) CALL SUFAINIT

!        2.2   Information from initial conditions file
IF(NSUPERSEDE==1) THEN
  IF (NUCONF /= 901 ) THEN
    CALL INI_IOSTREAM(KULOUT=NULOUT,LDIFS_INTERNAL=.TRUE.)
  ENDIF

  IF (NGRIBFILE==0 .AND. NUECMWF/=1) THEN

!   Let proc #1 only to be verbose
    IF (MYPROC==1) THEN
      IVERBOSE=1
    ELSE
      IVERBOSE=0
    ENDIF
    CALL FANMSG(IVERBOSE,NULOUT)
    CLFNISH=CFNISH
    IUNTIN=NINISH

    IF ( NPROC > 1 ) THEN

!      Wait until all nodes gives the same answer

       CHECK_EXIST : DO
         IEXIST = 0
         INQUIRE(FILE=CLFNISH,EXIST=LLEXIST)
         IF ( LLEXIST ) IEXIST=1
         CALL MPL_ALLGATHERV(IEXIST,IIEXIST,CDSTRING='SUARG:')
         IF ( ALL(IIEXIST==1) .OR. ALL(IIEXIST==0) ) THEN
            LLEXIST = ( IIEXIST(1) == 1 )
            EXIT CHECK_EXIST
         ENDIF
       ENDDO CHECK_EXIST
    ELSE
       INQUIRE(FILE=CLFNISH,EXIST=LLEXIST)
    ENDIF

!     WAIT UNTIL THE FILE CLFNISH EXIST. This part allows to execute
!     different applications in the same time without waiting that
!     the file CLFNISH exists.

    IF (LLEXIST) THEN
      CLLEC='CADRE LECTURE   '
      INBARI=0
      CALL FAITOU(IREP,IUNTIN,.TRUE.,CLFNISH,'OLD',.TRUE.,&
       & .TRUE.,1,1,INBARI,CLLEC)
      CALL LFIMST(IREP,IUNTIN,.FALSE.)
      NUFLEV=JPMXLE
      IF (.NOT.LUELAM) THEN
        CALL RIEN(CLLEC,NUSTTYP,UMUCEN,ULOCEN,USTRET,&
         & NUSMAX,NUDGL,NUDLON,IULOEN,IUMEN,ZGEMU,&
         & NUHTYP,NUFLEV,ZPREF,UVALH,UVBH,&
         & IUQUAD,1,JPMXGL,ZEPS,LLFICP,&
         & NULOUT)  
        DO JGL=1,(NUDGL+1)/2
          NULOEN(JGL)=IULOEN(JGL)
          NUMEN(JGL) =IUMEN(JGL)
        ENDDO
      ELSE
        CALL ERIEN(CLLEC,NUSTTYP,LLMAP,NUSMAX,NUDGL,&
         & NUDLON,IULIM,UGEMU,NUFLEV,ZPREF,UVALH,UVBH,&
         & ZEPS,NULOUT)  
        UMUCEN=1._JPRB
        ULOCEN=0._JPRB
        IF (LLMAP) THEN
          USTRET=1._JPRB
        ELSE
          USTRET=-1.0_JPRB
        ENDIF
        NUHTYP=0
        DO JGL=1,8
          NULIM(JGL)=IULIM(JGL)
        ENDDO
        NULOEN(:)=NUDLON
      ENDIF
      CALL FADIEX(IREP,IUNTIN,IDATEF)
      CALL FAVORI(IKNGRIB,IKNBPDG,IKNBCSP,IKSTRON,IKPUILA,IKDMOPL)
      CALL FAVEUR(IREP,IUNTIN,INGRIB,IKNBPDG,IKNBCSP,ISTRON,IKPUILA,IKDMOPL)
      IF (IKNGRIB /= INGRIB) THEN
        CALL FAGIOT(INGRIB,IKNBPDG,IKNBCSP,IKSTRON,IKPUILA,IKDMOPL)
        WRITE(NULOUT,*) 'The default value of NVGRIB is not that given by FAVORI:',&
         & IKNGRIB,'but that of the initial file :',INGRIB
      ENDIF
      CALL FAIRME(IREP,IUNTIN,'UNKNOWN')
      ISSSS = IDATEF (JD_SEM)+IDATEF (JD_SET)
      NUSSSS=MOD(ISSSS,ISDAY)
      Z_RINC=REAL(ISSSS,JPRB)/REAL(ISDAY,JPRB)
      INC=INT(Z_RINC)
      IF (INC /= 0) THEN
        CALL UPDCAL(IDATEF(JD_DAY),IDATEF(JD_MON),IDATEF(JD_YEA),INC,&
         & IDATEF(JD_DAY),IDATEF(JD_MON),IDATEF(JD_YEA),ILMO,NULOUT)  
      ENDIF
      NUDATE=IDATEF(JD_YEA)*10000+IDATEF(JD_MON)*100+IDATEF(JD_DAY)

      DO JLEV=0,NUFLEV
        UVALH(JLEV)=UVALH(JLEV)*ZPREF
      ENDDO
      IF (MYPROC/=1) THEN
!       Return to FA default verbosity
        IVERBOSE=1
        CALL FANMSG(IVERBOSE,NULOUT)
      ENDIF
    ELSE
      WRITE(NULOUT,*) ' INPUT FILE NOT FOUND '
      CALL ABOR1 ('INPUT FILE NOT FOUND')
    ENDIF

  ELSE

!*       2.4   Information from GRIB initial conditions files
!              LECMWF 

!        2.4.1 Information on spectral data          

    CLFNISH=CFNISH
    ITAG = 88888
    IF (NUCONF == 901) THEN
      INQUIRE(FILE=CLFNISH,EXIST=LLSPEXIST)
    ELSE
      LLSPEXIST = .TRUE.
    ENDIF
    IF (LLSPEXIST) THEN
      IF (NUCONF == 901 ) THEN
        CALL IGRIB_OPEN_FILE(IFILE, CLFNISH,'R')
        CALL IGRIB_NEW_FROM_FILE(IFILE,JGRIB)
        CALL IGRIB_GET_VALUE(JGRIB,'gridDefinitionTemplateNumber',ITYPEOFGRID)
        WRITE(NULOUT,*) 'typeofgrid=',ITYPEOFGRID
        CALL IGRIB_GET_VALUE(JGRIB,'dataDate',IUDATE)
        CALL IGRIB_GET_VALUE(JGRIB,'dataTime',IUSSSS)
        CALL IGRIB_GET_VALUE(JGRIB,'step',ISTEP)
        CALL IGRIB_GET_VALUE(JGRIB,'pentagonalResolutionParameterJ',NUSMAX)
        CALL IGRIB_GET_VALUE(JGRIB,'numberOfVerticalCoordinateValues',IUFLEV)
        IUFLEV=(IUFLEV/2)-1
        CALL IGRIB_GET_VALUE(JGRIB,'pv' ,ZPV)
        ZVALH(0:IUFLEV)=ZPV(0:IUFLEV)
        ZVBH(0:IUFLEV)=ZPV(IUFLEV+1:IUFLEV*2+1)
        CALL IGRIB_RELEASE(JGRIB)
        CALL IGRIB_CLOSE_FILE(IFILE)
      ELSE
        CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CLFNISH),CDMODE='r',KIOMASTER=1)
        CALL SETUP_IOREQUEST(YL_IOREQUEST,'SPECTRAL_FIELDS',LDGRIB=.TRUE.,&
         & LDINTONLY=.TRUE.,KCHUNKSIZE=1)
        CALL IO_INQUIRE(YL_IOSTREAM,YL_IOREQUEST,&
         & KDATE=IUDATE,KTIME=IUSSSS,KSTEP=ISTEP,&
         & KSMAX=NUSMAX,KFLEV=IUFLEV,PVALH=ZVALH,PVBH=ZVBH)
        CALL CLOSE_IOREQUEST(YL_IOREQUEST)
        CALL CLOSE_IOSTREAM(YL_IOSTREAM)
        ENDIF     
      NUSTTYP=1
      USTRET=1.0_JPRB
      UMUCEN=1.0_JPRB
      ULOCEN=0.0_JPRB
      CALL MOD_INI(IUDATE,IUSSSS,ISTEP,NUCONF,-1,& 
       & NSTEPINI,-1.0_JPRB,NULOUT) 
      IF (LUELAM) THEN
        ! Needs remake to new LAM domain definition parameters
        CALL ABOR1('SUARG: REMAKE ETGBSEC2')
!       CALL ETGBSEC2(ISEC2,ILEN2,ZSEC2,ILENR2,IUDGL,IUDLON,&
!        &UBETA,IUSMAX,IUSTTYP,NULIM(1),JPMXGL,&
!        &UGEMU,JPMXGL,IUFLEV,ZVALH,ZVBH,&
!        &ILEVMN,ILEVMX,NULOUT)
!                 Truncations
        NUSMAX=IUSMAX
        NUSTTYP=IUSTTYP
!                 Horizontal geometry
        USTRET=1._JPRB
        UMUCEN=1._JPRB
        ULOCEN=0._JPRB
      ENDIF

!              Vertical coordinate
      IF    (NUCONF == 201.OR.NUCONF == 202 &
         & .OR.NUCONF == 421.OR.NUCONF == 422 &
         & .OR.NUCONF == 521.OR.NUCONF == 522 &
         & )THEN  
        NUFLEV=1
        UVALH(0)=0._JPRB
        UVALH(1)=0._JPRB
        UVBH (0)=0._JPRB
        UVBH (1)=1._JPRB
      ELSE
        NUFLEV=IUFLEV
        DO JLEV=0,NUFLEV
          UVALH(JLEV)=ZVALH(JLEV)
          UVBH (JLEV)=ZVBH(JLEV)
        ENDDO
      ENDIF
!              Date
      NUDATE=IUDATE
      NUSSSS=IUSSSS

    ELSE
      NUSMAX=1
    ENDIF

!        2.4.2 Information on gridpoint data         

    CLFNIGG=CFNIGG
    IF (NUCONF == 901) THEN
      INQUIRE(FILE=CLFNIGG,EXIST=LLGGEXIST)
    ELSE
      LLGGEXIST = .TRUE.
    ENDIF
    IF (LLGGEXIST) THEN
      IF (NUCONF == 901 ) THEN
        CALL IGRIB_OPEN_FILE(IFILE, CLFNIGG,'r')
        CALL IGRIB_NEW_FROM_FILE(IFILE,JGRIB)
        CALL IGRIB_GET_VALUE(JGRIB,'gridDefinitionTemplateNumber',ITYPEOFGRID)
        WRITE(NULOUT,*) 'typeofgrid=',ITYPEOFGRID
        CALL IGRIB_GET_VALUE(JGRIB,'dataDate',IUDATE)
        CALL IGRIB_GET_VALUE(JGRIB,'dataTime',IUSSSS)
        CALL IGRIB_GET_VALUE(JGRIB,'step',ISTEP)
        CALL IGRIB_GET_VALUE(JGRIB,'NV',IUFLEV)
        CALL IGRIB_GET_VALUE(JGRIB,'Nj' ,IUDGL)
        WRITE(NULOUT,*) "Nj =",IUDGL
      
        USTRET=1._JPRB
        UMUCEN=1._JPRB
        ULOCEN=0._JPRB

        IUSTTYP=1
        NUSTTYP=IUSTTYP

        CALL IGRIB_GET_VALUE(JGRIB,'PLPresent' ,IPLPRESENT)
        WRITE(NULOUT,*) ' PLPresent=',IPLPRESENT
        IF (IPLPRESENT/=0) THEN
          CALL IGRIB_GET_VALUE(JGRIB,'pl',IULOEN)
        ELSE
          DO I=1,IUDGL
            IULOEN(I)=IUDGL*2
          ENDDO
        ENDIF
        CALL IGRIB_CLOSE_FILE(IFILE)
      ELSE
        CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CLFNIGG),CDMODE='r',KIOMASTER=1)
        CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,&
          & LDINTONLY=.TRUE.,KCHUNKSIZE=1)
        CALL IO_INQUIRE(YL_IOSTREAM,YL_IOREQUEST,&
          & KDATE=IUDATE,KTIME=IUSSSS,KSTEP=ISTEP,&
  !       & KFLEV=IUFLEV,PVALH=ZVALH,PVBH=ZVBH,&
          & KLATS=IUDGL,KLONS=IULOEN)
        CALL CLOSE_IOREQUEST(YL_IOREQUEST)
        CALL CLOSE_IOSTREAM(YL_IOSTREAM)
       ENDIF

      CALL MOD_INI(IUDATE,IUSSSS,ISTEP,NUCONF,-1,&  
       & NSTEPINI,-1.0_JPRB,NULOUT)  
      IF (.NOT.LUELAM) THEN

!                 Grid definition
        NUHTYP=0
        NUDGL =IUDGL
!        NUDLON=IUDGL*2
        NUDLON=MAXVAL(IULOEN(1:(NUDGL+1)/2))
        DO JGL=1,(NUDGL+1)/2
          NULOEN(JGL)=IULOEN(JGL)
          IF(NULOEN(JGL) /= NULOEN(1)) NUHTYP=2
        ENDDO
      ELSE
        ! Needs remake to new LAM domain definition parameters
        CALL ABOR1('SUARG: REMAKE ETGBSEC2')
!       CALL ETGBSEC2(ISEC2,ILEN2,ZSEC2,ILENR2,IUDGL,IUDLON,&
!        &UBETA,IUSMAX,IUSTTYP,NULIM,JPMXGL,&
!        &UGEMU,JPMXGL,IUFLEV,ZVALH,ZVBH,&
!        &ILEVMN,ILEVMX,NULOUT)

!                 Grid definition
        NUHTYP=0
        NUDGL =IUDGL
        NUDLON=IUDLON
        NULOEN(:)=NUDLON
      ENDIF
      IF (LLSPEXIST) THEN
        IABORT=0
        IF (NUDATE /= IUDATE) THEN
          WRITE(NULOUT,*) ' NOT THE SAME DATE !'
          WRITE(NULOUT,*) ' IN SP FILE : ',NUDATE
          WRITE(NULOUT,*) ' IN GG FILE : ',IUDATE
!                    IABORT=1 
        ELSE
          IF (NUSSSS /= IUSSSS) THEN
            WRITE(NULOUT,*) ' NOT THE SAME TIME !'
            WRITE(NULOUT,*) ' IN SP FILE : ', NUSSSS
            WRITE(NULOUT,*) ' IN GG FILE : ', IUSSSS
!                      IABORT=1
          ENDIF
        ENDIF
        IF (IABORT == 1) THEN
          CALL ABOR1('SUARG: ABOR1 CALLED')
        ENDIF
      ELSE
!                Geometry
        IF (.NOT.LUELAM) THEN
          NUSTTYP=IUSTTYP
          USTRET=ZSTRET
          UMUCEN=ZMUCEN
          ULOCEN=ZLOCEN
        ELSE
          USTRET=1._JPRB
          UMUCEN=1._JPRB
          ULOCEN=0._JPRB
        ENDIF
        IUFLEV=NUFLEV
        DO JLEV=0,NUFLEV
          UVALH(JLEV)=ZVALH(JLEV)
          UVBH (JLEV)=ZVBH(JLEV)
        ENDDO
!                Date
!                 NUDATE=IUDATE
!                 NUSSSS=IUSSSS
      ENDIF

    ELSE
      WRITE(NULOUT,*) ' GRIDPOINT FILE SEARCHED:',CLFNIGG
      WRITE(NULOUT,*) ' GRIDPOINT INPUT FILE NOT FOUND'
      IF (LLSPEXIST) THEN
        NUHTYP=0
        NUDGL =1
        NUDLON=2
      ELSE
        WRITE(NULOUT,*) ' >>> NO INPUT FILE AT ALL !'
        CALL ABOR1('SUARG: ABOR1 CALLED')
      ENDIF
    ENDIF

  ENDIF

  NUCMAX=NUSMAX

ENDIF

!*       3.    PRINT FINAL VALUES
!              ------------------

WRITE(UNIT=NULOUT,FMT='('' MODULE YOMARG: '')')
WRITE(UNIT=NULOUT,FMT='('' * NAMARG '')')
IF(NUCONF /= -1) WRITE(UNIT=NULOUT,FMT='(''   NUCONF = '',I3)') NUCONF
IF(NUECMWF == 1) WRITE(UNIT=NULOUT,FMT='(''   ECMWF '')')
IF(NUECMWF == 2) WRITE(UNIT=NULOUT,FMT='(''   METEO-FRANCE '')')
WRITE(UNIT=NULOUT,FMT='(''   CUNMEXP = '',A4)') CNMEXP
!IF(NUSLAG == 1) WRITE(UNIT=NULOUT,FMT='(''   EULERIAN ADVECTION '')')
!IF(NUSLAG == 3) WRITE(UNIT=NULOUT,FMT='(''   SEMI-LAGRANGIAN ADVECTION '')')
IF (.NOT.LUELAM) THEN
  WRITE(UNIT=NULOUT,FMT='(''   GLOBAL ARPEGE/IFS MODEL '')')
ELSE
  WRITE(UNIT=NULOUT,FMT='(''   LIMITED AREA MODEL '')')
ENDIF
WRITE(UNIT=NULOUT,FMT='('' NSUPERSEDE = '',I3,'' NFPSERVER = '',I3)') NSUPERSEDE, NFPSERVER
WRITE(UNIT=NULOUT,FMT='('' NGRIBFILE = '',I2)') NGRIBFILE

IF(NSUPERSEDE==1) THEN
  WRITE(UNIT=NULOUT,FMT='('' * Other variables '')')
  WRITE(UNIT=NULOUT,FMT='('' NUCMAX = '',I5)') NUCMAX  
  WRITE(UNIT=NULOUT,FMT='('' NUDATE = '',I8,'' NUSSSS = '',I6)') NUDATE, NUSSSS  
  WRITE(UNIT=NULOUT,FMT='('' NSTEPINI = '',I5)') NSTEPINI
  WRITE(UNIT=NULOUT,FMT='('' NUSMAX = '',I5,'' NUDGL = ''&
   & ,I4,'' NUDLON = '',I5,'' NUFLEV = '',I3)')&
   & NUSMAX,NUDGL,NUDLON,NUFLEV  

  DO JLEV=0,NUFLEV
    WRITE(NULOUT,FMT='('' UVALH('',I3,'') = '',F20.10,&
     & ''     UVBH('',I3,'') = '',F20.10)')&
     & JLEV,UVALH(JLEV),JLEV,UVBH(JLEV)  
  ENDDO

  IF(.NOT.LUELAM) THEN
    WRITE(UNIT=NULOUT,FMT='('' NUSTTYP = '',I1,'' USTRET = ''&
     & ,F9.6,'' UMUCEN = '',F9.6,'' ULOCEN = '',F9.6)')&
     & NUSTTYP,USTRET,UMUCEN,ULOCEN  
    WRITE(UNIT=NULOUT,FMT='('' NUHTYP = '',I1)') NUHTYP
    IDGNH=(NUDGL+1)/2
    DO JGL=1,IDGNH
      WRITE(NULOUT,&
       & FMT='(2('' NULOEN ('',I4,'') = '',I4,&
       & '' NUMEN   ('',I4,'') = '',I4,3X))')&
       & JGL,NULOEN(JGL),JGL,NUMEN(JGL),&
       & JGL+IDGNH,NULOEN(IDGNH+1-JGL),JGL+IDGNH,NUMEN(IDGNH+1-JGL)  
    ENDDO
  ELSE
    WRITE(UNIT=NULOUT,FMT='('' NUSTTYP = '',I3)') NUSTTYP
    WRITE(UNIT=NULOUT,FMT='(&
     & '' ISTROE = '',I3, '' IDOM   = '',I3, &
     & '' IDLUN  = '',I3, '' IDLUX  = '',I3)')&
     & NULIM(1), NULIM(2), NULIM(3), NULIM(4)   
    WRITE(UNIT=NULOUT,FMT='(&
     & '' IDGUN  = '',I3, '' IDGUX  = '',I3,&
     & '' IDZONL = '',I3, '' IDZONG = '',I3)')&
     & NULIM(5), NULIM(6), NULIM(7), NULIM(8)   
    WRITE(UNIT=NULOUT,FMT='('' ZRPK = '',E13.3)') UGEMU(2)
    WRITE(UNIT=NULOUT,FMT='(&
     & '' ZLON0 = '',E13.3, '' ZLAT0 = '',E13.3, &
     & '' ZLONC = '',E13.3, '' ZLATC = '',E13.3)') &
     & UGEMU(3), UGEMU(4), UGEMU(5), UGEMU(6)    
    WRITE(UNIT=NULOUT,FMT='(&
     & '' ZDELX = '',E13.3, '' ZDELY = '',E13.3)') &
     & UGEMU(7), UGEMU(8)  
  ENDIF
ENDIF


END ASSOCIATE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUARG',1,ZHOOK_HANDLE)
END SUBROUTINE SUARG
