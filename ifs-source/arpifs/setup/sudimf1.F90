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

SUBROUTINE SUDIMF1(YDMODEL)

!     ----------------------------------------------------------------------
!**** *SUDIMF1 * - Set up number of fields dimensions

!     Purpose.
!     --------
!           Initialization of YOMDIMF and YOMDPHY with some prints

!**   Interface.
!     ----------
!        *CALL* *SUDIMF1*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Tomas Wilhelmsson *ECMWF* 
!      Original based on SUDIM1 : 2013-08-01 

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      E.Dutra/G.Arduini (Jan 2018): added NCSNEC snow vertical dimension
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!     ----------------------------------------------------------------------

USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULOUT, NULNAM, NULERR
USE YOMCT0    , ONLY : LR3D, LR2D, NCONF, LECMWF, LELAM,L_OOPS
USE YOMVAR    , ONLY : LJC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL), INTENT(INOUT), TARGET :: YDMODEL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), POINTER :: NVXP, NVXP2, NVEXTR,  NVEXTRRAD,  NVEXTRDYN,  &
 & NVXTR2, NCEXTR,  NVCLIS,  NTOZ1D,  NTOZ2D,  NTOZ3D,  NTSSG, NVECOUT, NCSNEC
LOGICAL, POINTER :: LTPROF, LDIRCLSMOD
LOGICAL :: LEURBAN

#include "namdphy.nam.h"

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDIMF1',0,ZHOOK_HANDLE)
ASSOCIATE(YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDMCC=>YDMODEL%YRML_AOC%YRMCC,YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, &
 & YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF)
ASSOCIATE(LADER=>YDDIMF%LADER, LSPT=>YDDIMF%LSPT, LUVDER=>YDDIMF%LUVDER, &
 & LVOR=>YDDIMF%LVOR, NFC2D=>YDDIMF%NFC2D,  &
 & NCOM=>YDDPHY%NCOM, NCSI=>YDDPHY%NCSI,    &
 & NCSS=>YDDPHY%NCSS, NCXP=>YDDPHY%NCXP, NTILES=>YDDPHY%NTILES, &
 & LSLHD=>YDDYNA%LSLHD,LSLHD_STATIC=>YDDYNA%LSLHD_STATIC,L3DTURB=>YDDYNA%L3DTURB,LCOMAD=>YDDYNA%LCOMAD, &
 & LEFLAKE=>YDEPHY%LEFLAKE, LEURBAN=>YDEPHY%LEURBAN, LEO3CH=>YDEPHY%LEO3CH, LEOCML=>YDEPHY%LEOCML, &
 & LO3CH_BMS=>YDEPHY%LO3CH_BMS, LO3CH_HLO=>YDEPHY%LO3CH_HLO, LOCMLTKE=>YDEPHY%LOCMLTKE, &
 & LMCC01=>YDMCC%LMCC01, LDIAGTURB_EC=>YDEPHY%LDIAGTURB_EC, LDIAGTURBGRAD_EC=>YDEPHY%LDIAGTURBGRAD_EC,&
 & LOZONE=>YDPHY%LOZONE, LMSE=>YDARPHY%LMSE, &
 & JPCXP=>YDPHYDS%JPCXP, JPVXP=>YDPHYDS%JPVXP)
! Associate pointers for variables in namelist
NVXP       => YDDPHY%NVXP
NVXP2      => YDDPHY%NVXP2
NVEXTR     => YDDPHY%NVEXTR
NVEXTRRAD  => YDDPHY%NVEXTRRAD
NVEXTRDYN  => YDDPHY%NVEXTRDYN
NVXTR2     => YDDPHY%NVXTR2
NVECOUT    => YDDPHY%NVECOUT
NCEXTR     => YDDPHY%NCEXTR
NVCLIS     => YDDPHY%NVCLIS
NTOZ1D     => YDDPHY%NTOZ1D
NTOZ2D     => YDDPHY%NTOZ2D
NTOZ3D     => YDDPHY%NTOZ3D
NTSSG      => YDDPHY%NTSSG
LTPROF     => YDDPHY%LTPROF
LDIRCLSMOD => YDDPHY%LDIRCLSMOD
NCSNEC     => YDDPHY%NCSNEC

!     ------------------------------------------------------------------
!*       1. SET DEFAULT VALUES FOR NAMDPHY VARIABLES, AND READ NAMDPHY.
!           -----------------------------------------------------------

NVXP=0
NCXP=0 ! not in NAMDPHY for the time being, but could be in the future

NVXP2=0
IF (LECMWF) THEN
  NVECOUT=1
ELSE
  NVECOUT=0
ENDIF

NVEXTR=0
NVEXTRRAD=0
NVEXTRDYN=0
NVXTR2=0
NCEXTR=0
NTSSG=0
LTPROF=.FALSE.
LDIRCLSMOD=.FALSE.

! - set number of reservoirs: NCSS
IF(.NOT. LR2D) THEN
  IF (LECMWF) THEN
    NCSS=4
  ELSE
    NCSS=1
  ENDIF
ELSE
  NCSS=0
ENDIF

IF (LECMWF) THEN
  IF ((NCONF == 1.OR.NCONF == 302.OR.L_OOPS).AND.LEO3CH) THEN
    NTOZ3D=0
    NTOZ2D=1
    NTOZ1D=0
    IF (LO3CH_BMS.OR.LO3CH_HLO) THEN
      NVCLIS=7
    ELSE
      NVCLIS=8
    ENDIF
  ELSE
    NTOZ3D=0
    NTOZ2D=0
    NTOZ1D=0
    NVCLIS=0
  ENDIF
ELSE
  NTOZ3D=0
  NTOZ2D=0
  NTOZ1D=0
  NVCLIS=0
ENDIF

IF (ANY(SPREAD(NCONF,1,11)==(/1,131,302,401,501,601,701,801,901,903,904/))) THEN
  NCSNEC=1
ELSE
  NCSNEC=0
ENDIF

CALL POSNAM(NULNAM,'NAMDPHY')
READ(NULNAM,NAMDPHY)

!     ------------------------------------------------------------------
!*       2. CHECKINGS AND RESETTINGS FOR YOMDPHY VARIABLES.
!           -----------------------------------------------

!        2.1    Checking and resettings on NAMDPHY variables.

! Reset NVXP2.
IF (NCONF == 701) NVXP2=11

! Reset NCSS
IF (LMSE) THEN
  LTPROF=.FALSE.
  WRITE(NULOUT,*) ' INFO ---- LTPROF reset to .F. when using surfex '
ENDIF
IF (.NOT.LECMWF .AND. LTPROF) NCSS=3

! Wrong value for NTSSG?
IF (NTSSG == 1) THEN
  WRITE(UNIT=NULERR,FMT='('' SUDIMF1: NTSSG == 1 IS NOT IMPLEMENTED'')')
  CALL ABOR1('SUDIMF1: NTSSG == 1 IS NOT IMPLEMENTED ')
ENDIF
IF (NTSSG > 0.AND..NOT.LMCC01)THEN
  WRITE(UNIT=NULERR,FMT='('' SUDIMF1: NTSSG MUST BE 0 IF LMCC01=F '')')
  CALL ABOR1('SUDIMF1: NTSSG MUST BE 0 IF LMCC01=F ')
ENDIF

! Reset NTSSG to 0 for some configurations.
IF (ANY(SPREAD(NCONF,1,4) == (/923,931,932,933/))) NTSSG=0

! Inconsistent values for NVCLIS, NTOZ1D, NTOZ2D, NTOZ3D?
IF (LOZONE) THEN
  IF (NVCLIS*(NTOZ1D+NTOZ2D+NTOZ3D) < 4) THEN
    WRITE(UNIT=NULERR,FMT='('' SUDIMF1: LOZONE INCONSISTENT WITH YOMDPHY'')')
    WRITE(UNIT=NULERR,FMT='('' SUDIMF1: LOZONE=T => NVCLIS*(NTOZ1D+NTOZ2D+NTOZ3D) must be >= 4 '')')
    CALL ABOR1('SUDIMF1: LOZONE INCONSISTENT WITH YOMDPHY ')
  ENDIF
ENDIF
IF (LEO3CH) THEN
! nvclis 7 for BMS and HLO O3 schemes
 IF (LO3CH_BMS.OR.LO3CH_HLO) THEN
  IF (NVCLIS*(NTOZ1D+NTOZ2D+NTOZ3D) < 7) THEN
    WRITE(UNIT=NULERR,FMT='('' SUDIMF1: LEO3CH INCONSISTENT WITH YOMDPHY'')')
    WRITE(UNIT=NULERR,FMT='('' SUDIMF1: LEO3CH=T => NVCLIS*(NTOZ1D+NTOZ2D+NTOZ3D) must be >= 7 '')')
    CALL ABOR1('SUDIMF1: LEO3CH INCONSISTENT WITH YOMDPHY ')
  ENDIF
 ELSE
  IF (NVCLIS*(NTOZ1D+NTOZ2D+NTOZ3D) < 8) THEN
    WRITE(UNIT=NULERR,FMT='('' SUDIMF1: LEO3CH INCONSISTENT WITH YOMDPHY'')')
    WRITE(UNIT=NULERR,FMT='('' SUDIMF1: LEO3CH=T => NVCLIS*(NTOZ1D+NTOZ2D+NTOZ3D) must be >= 8 '')')
    CALL ABOR1('SUDIMF1: LEO3CH INCONSISTENT WITH YOMDPHY ')
  ENDIF
 ENDIF ! end of YOEPHY%LO3CH_BMS
ENDIF ! end of LEO3CH

! Reset NVCLIS to 0 for some configurations.
IF (ANY(SPREAD(NCONF,1,4) == (/923,931,932,933/))) NVCLIS=0

! Inconsistent values for NVEXTR and NVEXTRDYN?
IF(NVEXTR < NVEXTRDYN) THEN
  WRITE(NULERR,*) ' SUDIMF1: NVEXTR=',NVEXTR,' NVEXTRDYN=',NVEXTRDYN
  CALL ABOR1('SUDIMF1: NVEXTR cannot be < NVEXTRDYN')
ENDIF

! Wrong values for NVXP and NCXP?
IF(NVXP > JPVXP) THEN
  WRITE(NULERR,'(''NVXP MUST BE SMALLER OR EQ. TO JPVXP'')')
  CALL ABOR1(' SUDIMF1: PROBLEM WITH NVXP ')
ENDIF
IF(NCXP > JPCXP) THEN
  WRITE(NULERR,'(''NCXP MUST BE SMALLER OR EQ. TO JPCXP'')')
  CALL ABOR1(' SUDIMF1: PROBLEM WITH NCXP ')
ENDIF

!        2.2    Other YOMDPHY variables: NCSI, NCSNEC, NTILES, NCOM

NCSI=0
! NCSNEC=0
NTILES=0
NCOM=0
IF (LECMWF .AND. LR3D) THEN
  NCSI=4
! DEFAULT VALUE
!   NCSNEC=1
  NTILES=8
  IF (LEFLAKE) NTILES=9 !FLAKE
  IF (LEURBAN) NTILES=10 !URBAN
  IF (LEOCML) NCOM=35 !KPP
  IF (LOCMLTKE) NCOM=9 !TKE
ENDIF
!     ------------------------------------------------------------------
!*       3. CHECKINGS AND RESETTINGS FOR YOMDIMF VARIABLES.
!           -----------------------------------------------

!!OM  NPPM is now (CY45) declared as a parameter in PARDIM
!!NPPM=4

IF(LR3D .OR. ANY(SPREAD(NCONF,1,2) == (/701,901/))) THEN  
  NFC2D=1
  LSPT=.TRUE.
  IF(YDDYNA%LSLAG) THEN
    IF (YDDYNA%LNHDYN.OR.(LSLHD.AND..NOT.LSLHD_STATIC)&
      & .OR.(NCONF == 131.AND.LJC).OR.L3DTURB.OR.LCOMAD.OR.(LDIAGTURB_EC.AND.LDIAGTURBGRAD_EC)) THEN
      LUVDER=.TRUE.
    ELSE
      LUVDER=.FALSE.
    ENDIF
  ELSE
    LUVDER=.TRUE.
  ENDIF
  IF(NCONF == 701)THEN
    LADER=.FALSE.
  ELSE
    LADER=.TRUE.
  ENDIF
  IF(NCONF==1.AND..NOT.(LELAM.OR.YDDYNA%LNHDYN)) THEN
    LVOR=((.NOT.YDDYNA%LSLAG).OR.(YDDYNA%LSLAG.AND.LSLHD)).AND.LADER
  ELSE
    LVOR=LADER
  ENDIF
ELSEIF(LR2D) THEN  
  NFC2D=1
  LSPT=.FALSE.
  LUVDER=.NOT.YDDYNA%LSLAG
  LADER=.TRUE.
  LVOR=LADER
ELSEIF (ANY(SPREAD(NCONF,1,4) == (/923,931,932,933/))) THEN
  NFC2D=0
  LSPT=.FALSE.
  LUVDER=.FALSE.
  LADER=.FALSE.
  LVOR=LADER
ELSE
  WRITE(UNIT=NULERR,FMT='('' SUDIMF1: CONFIGURATION NOT EXPECTED'')')
  CALL ABOR1('SUDIMF1: part 3 CONFIGURATION NOT EXPECTED')
ENDIF

!     ------------------------------------------------------------------
!*       4. PRINTINGS.
!           ----------

WRITE(NULOUT,*) ' PRINTINGS IN SUDIMF1: PHYSICS DIMENSIONING (YOMDPHY) '
WRITE(NULOUT,*) ' -- PHYSICS DIMENSIONING (YOMDPHY) '
WRITE(NULOUT,*) ' NCSS = ',NCSS
WRITE(NULOUT,*) ' NCSNEC = ',NCSNEC
WRITE(NULOUT,*) ' NVXP = ',NVXP,' NVXP2 = ',NVXP2
WRITE(NULOUT,*) ' NCXP = ',NCXP,' NTILES = ',NTILES  
WRITE(NULOUT,*)&
 & ' NVEXTR = ',NVEXTR,' NVEXTRRAD = ',NVEXTRRAD,' NVXTR2 = ',NVXTR2,&
 & ' NVEXTRDYN = ',NVEXTRDYN,' NVECOUT = ',NVECOUT
WRITE(NULOUT,*) ' NCEXTR= ',NCEXTR
WRITE(NULOUT,*)' NVCLIS = ',NVCLIS,' NTSSG = ',NTSSG,&
 & ' NTOZ1D = ',NTOZ1D,' NTOZ2D = ',NTOZ2D,' NTOZ3D = ',NTOZ3D,&
 & ' LTPROF = ',LTPROF
WRITE(NULOUT,*) ' LDIRCLSMOD = ',LDIRCLSMOD

! ky: for the time being, final values of LVOR, LUVDER and LADER
!     are provided in SUDIMF2, and printed in SUDIMF2.
!     Moving SUDIMF2 code into SUDIMF1 must be studied.

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDIMF1',1,ZHOOK_HANDLE)
END SUBROUTINE SUDIMF1
