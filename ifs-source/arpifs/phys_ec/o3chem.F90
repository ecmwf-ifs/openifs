! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE O3CHEM ( YDRIP,YDEPHY,YDOZO,YDPHY2,KIDIA,KFDIA,KLON,KTDIA,KLEV,KVCLIS,PGEMU,PMU0,&
 ! - INPUT  2D .
 & PAPRS,PAPRSF,PKOZO,PDELP,PT,PO3,&
 ! - UPDATE 2D .
 & PTENO3)  

!**** *O3CHEM* PHOTOCHEMICAL FLUX AND TENDENCY

!     PURPOSE.
!     -------
!     - CALCULATION OF OZONE FLUX AND TENDENCY
!       DUE TO PHOTOCHEMICAL SOURCES AND SINKS
!       (MIXING RATIO RELATIVE TO MASS)

!**   Interface.
!     ----------
!        *CALL* *O3CHEM*

!-----------------------------------------------------------------------

! -   INPUT ARGUMENTS.
!     -------------------

! - DIMENSIONS ETC.

!       KIDIA     : START OF HORIZONTAL LOOP
!       KFDIA     : END   OF HORIZONTAL LOOP
!       KLON      : HORIZONTAL DIMENSION
!       KTDIA     : START OF THE VERTICAL LOOP IN THE PHYSICS
!                   (IF SOME LEVELS ARE SKIPPED AT THE TOP OF THE MODEL).
!       KLEV      : END OF VERTICAL LOOP AND VERTICAL DIMENSION
!       KVCLIS    : NUMBER OF PHOTOCHEMICAL COEFFICIENTS
!       PGEMU     : SINE OF LATITUDE
!       PMU0      : LOCAL COSINE OF INSTANTANEOUS SOLAR ZENITH ANGLE

! - 2D (0:KLEV) .

!       PAPRS     : PRESSURE ON HALF-LEVELS.

! - 2D (1:KLEV) .

!       PAPRSF    : PRESSURE ON FULL LEVELS. 
!       PDELP     : LAYER THICKNESS IN PRESSURE UNITS.
!       PKOZO     : PHOTOCHEMICAL COEFFICIENTS COMPUTED FROM A
!                   2D PHOTOCHEMICAL MODEL (KVCLIS=8)
!       PT        : TEMPERATURE
!       PO3       : MASS MIXING RATIO OF OZONE

! -   UPDATED ARGUMENTS.
!     -----------------
! - 2D (1:KLEV) TENDENCY .
!       PTENO3     :  PHOTOCHEMICAL TENDENCY FOR OZONE

!-----------------------------------------------------------------------

!     Externals.
!     ---------

!     Method.
!     -------

!      D(PO3) / DT = PKOZO2 + PKOZO3 * ( PO3  - PKOZO1 )
!                           + PKOZO5 * ( PT   - PKOZO4 )
!                           + PKOZO7 * ( ZCO3 - PKOZO6 )
!                           - PKOZO8 * PPEECL * PPEECL * PO3
!      where PPEECL is the Equivalent Chlorine content of the stratosphere
!      for the actual year. The chemistry in only on in daylight, and the
!      heterogeneous chemistry is only turned on below a threshold
!      temperature (ZTEMP = 195 K).

!     Reference.
!    -----------
!         CARIOLLE ET DEQUE, JGR,91,10.825,1986

!     Author
!    --------
!      96-12-05, A. Untch (based on ACOZONE by A. Lasserre-Bigorry)

!     Modifications :
!    ----------------
!      2002-01-17, A. Dethof: Introduce check for negative ozone:
!                             ZANEX(JLON,JLEV)=MAX(1.E-9_JPRB,ZANEX(JLON,JLEV)).
!                             Use chemistry all the time, apart from
!                             heterogeneous chemistry which is still only
!                             used in sunlight.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      2006-09-08, A. Untch: Changed solar elevation angle from 90deg to 87deg
!                             in heterogeneous chemistry (recommended by D. Cariolle)
!                             Scaling of het. chem. term by p/p0 (recommended 
!                             by D. Cariolle)
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      A. Untch (Oct 2008): Mods for chemistry version 2.9
!      B. Monge-Sanz (Feb 2015): Switching to BMS coefficients for new O3 scheme, based 
!                                on method from Monge-Sanz et al. (2011, ACP)
!      J. Flemming (Dec 2017): LO3CHEM_SAFE switch to clip positive PKOZO and very large T increments   
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG,RPI
USE YOMPHY2  , ONLY : TPHY2
USE YOMRIP   , ONLY : TRIP
USE YOMOZO   , ONLY : TOZO
USE YOEPHY   , ONLY : TEPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(TOZO)        ,INTENT(INOUT) :: YDOZO
TYPE(TPHY2)       ,INTENT(INOUT) :: YDPHY2
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVCLIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(KLON,KLEV,KVCLIS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENO3(KLON,KLEV) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZCO3(KLON,KLEV),ZANEX(KLON,KLEV),ZDELPS(KLON,KLEV)
REAL(KIND=JPRB) :: ZKOZO(KLON,KLEV,KVCLIS)

INTEGER(KIND=JPIM) :: JLEV, JLON

REAL(KIND=JPRB) :: ZEPS, ZRAPP,ZLIGHT(KLON), ZTEMP, ZTANLAT, ZRPREF, ZCOSSOEL
REAL(KIND=JPRB) ::  ZPRSAFE, ZDTSAFE, ZT, ZSCALE, ZOZO3SAFE, ZOZO3MAX, ZRLXTIME, ZA1, ZA2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('O3CHEM',0,ZHOOK_HANDLE)
ASSOCIATE(LO3CH_BMS=>YDEPHY%LO3CH_BMS, LO3CH_HLO=>YDEPHY%LO3CH_HLO, LO3CH_OLDVER=>YDEPHY%LO3CH_OLDVER, &
 & TEECL2=>YDOZO%TEECL2,LO3CH_SAFE=>YDEPHY%LO3CH_SAFE,  &
 & RDECLI=>YDRIP%RDECLI, TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

!         PRELIMINARY COMPUTATIONS AND VARIOUS INITIALIZATIONS

ZEPS=1.E-11_JPRB
! ZCOSSOEL= cos(87deg) (recommended by D Cariolle)
ZCOSSOEL=0.052336_JPRB

ZOZO3MAX=-1.0E-10_JPRB
IF (LO3CH_SAFE) THEN
  ! MAX allowed T difference  
  ZDTSAFE=15.0_JPRB   
  ! above more drastic relaxation is enforced (1000 Pa )   
  ZPRSAFE=1000.0_JPRB  
! 24  hour time scale 
  ZRLXTIME=1.0_JPRB*24.0_JPRB
! 1 hour time scale 
!  ZRLXTIME=1.0_JPRB
  ZOZO3SAFE=-1.0_JPRB/(ZRLXTIME*3600.0_JPRB)
ENDIF 

!         CHEMISTRY PARAMETERS
!         1. Use chemistry in daylight only and there scale relaxation 
!            rate by the fractional length of the day (limited < 2)
!         2. Turn heterogeneous destruction on below 195K,( in presence of
!            clouds) and at high latidudes, above 50N/S

! local copy 
ZKOZO(KIDIA:KFDIA,1:KLEV,1:KVCLIS)=PKOZO(KIDIA:KFDIA,1:KLEV,1:KVCLIS)

DO JLON=KIDIA,KFDIA
  ZTANLAT=TAN(RDECLI)*PGEMU(JLON)/SQRT(ZEPS+1.0_JPRB-PGEMU(JLON)**2)
  IF(ABS(ZTANLAT) >= 1.0_JPRB)THEN
    ZLIGHT(JLON)=0.5_JPRB+SIGN(0.5_JPRB,PMU0(JLON)-ZCOSSOEL)
  ELSE
    ZLIGHT(JLON)=(0.5_JPRB+SIGN(0.5_JPRB,PMU0(JLON)-ZCOSSOEL))&
     & /MAX((1.0_JPRB-ACOS(ZTANLAT)/RPI),0.5_JPRB)  
  ENDIF
ENDDO
ZTEMP=195._JPRB
!    COMPUTATION OF UPPER HALF-LAYER DEPTH FOR THE WHOLE ATMOSPHERIC COLUMN

!       ZDELPS   : UPPER HALF-LAYER DEPTH

DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    ZDELPS(JLON,JLEV)=PAPRSF(JLON,JLEV)-PAPRS (JLON,JLEV-1)
  ENDDO
ENDDO

IF(LO3CH_OLDVER) THEN
! reference pressure 50hPa
!(Scaling with p/50hPa only needed for chemistry version 2.3. 
! For version 2.9 this is already incorporated in the coefficients.)
  ZRPREF=1.0_JPRB/5000.0_JPRB
ENDIF

!    INITIALIZATION TO ZERO
ZCO3(:,1)=0.0_JPRB

!*   COMPUTATION OF MIXING RATIO AT TIME T+DT

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA

  ZT=PT(JLON,JLEV) 

  IF (LO3CH_SAFE) THEN
! T increments smaller than ZDTSAFE 
    IF ( ZT - ZKOZO(JLON,JLEV,4) > ZDTSAFE) ZT = ZKOZO(JLON,JLEV,4) + ZDTSAFE        
    IF ( ZT - ZKOZO(JLON,JLEV,4) < -1.0_JPRB * ZDTSAFE) ZT = ZKOZO(JLON,JLEV,4) - ZDTSAFE        

    IF(.NOT.LO3CH_HLO) THEN
! clip positive values 
      ZKOZO(JLON,JLEV,3)  =MIN(ZOZO3MAX,  ZKOZO(JLON,JLEV,3) ) 
      IF ( PAPRS(JLON,JLEV) <= ZPRSAFE  ) THEN    
!    weighted mean of ZKOZO(JLON,JLEV,3) and ZOZO3SAFE depending on model pressure   
!       ZA1 =MIN(1.0_JPRB,MAX(0.0_JPRB, (LOG(PAPRS(JLON,JLEV))-LOG(PAPRS(JLON,1)))/(LOG(ZPRSAFE)-LOG(PAPRS(JLON,1)))))
        ZA1 =MIN(1.0_JPRB,MAX(0.0_JPRB, ((PAPRS(JLON,JLEV))-(PAPRS(JLON,1)))/((ZPRSAFE)-(PAPRS(JLON,1)))))
        ZA2 = 1.0_JPRB - ZA1   
        ZKOZO(JLON,JLEV,3)  =MIN(ZKOZO(JLON,JLEV,3) , ZA2 * ZOZO3SAFE + ZA1 * ZKOZO(JLON,JLEV,3)) 
        ZSCALE = ABS(PKOZO(JLON,JLEV,3)) / ABS(ZKOZO(JLON,JLEV,3)) 
! reduce T and TCO3 parameter by the same fraction as O3 term is increased  
!       ZKOZO(JLON,JLEV,5) =  ZSCALE * ZKOZO(JLON,JLEV,5)
!       ZKOZO(JLON,JLEV,7) =  ZSCALE * ZKOZO(JLON,JLEV,7)
      ENDIF
    ENDIF
  ENDIF 

    ZANEX(JLON,JLEV)=PO3(JLON,JLEV)&
     & +TSPHY&
     & *( ZKOZO(JLON,JLEV,2)&
     & +ZKOZO(JLON,JLEV,5)*(ZT          -ZKOZO(JLON,JLEV,4))&
     & +ZKOZO(JLON,JLEV,7)*(ZCO3(JLON,JLEV)-ZKOZO(JLON,JLEV,6))&
     & -ZKOZO(JLON,JLEV,1)*                 ZKOZO(JLON,JLEV,3) )  
    IF(LO3CH_OLDVER) THEN
!    Scaling of heterogeneous chemistry term with p/50hPa 
!    (*PAPRSF(JLON,JLEV)*ZRPREF) only needed for version 2.3
      ZRAPP=1.0_JPRB&
       & -TSPHY&
       & *( ZKOZO(JLON,JLEV,3)&
       & +ZKOZO(JLON,JLEV,7)*ZDELPS(JLON,JLEV)/RG&
       & - ZLIGHT(JLON)*(0.5_JPRB+SIGN(0.5_JPRB,ZTEMP-PT(JLON,JLEV)))&
       & *PAPRSF(JLON,JLEV)*ZRPREF&
       & *TEECL2*ZKOZO(JLON,JLEV,8) )
    ELSE
!    When using BMS or HLO coeffs (7 coeffs with het. chem. embedded)
      IF (LO3CH_BMS.OR.LO3CH_HLO) THEN
        ZRAPP=1.0_JPRB&
       & -TSPHY&
       & *( ZKOZO(JLON,JLEV,3)&
       & +ZKOZO(JLON,JLEV,7)*ZDELPS(JLON,JLEV)/RG)
      ELSE
!     When using Cariolle's coeffs --> use additional heterogeneous term (8th coeff)
      ZRAPP=1.0_JPRB&
       & -TSPHY&
       & *( ZKOZO(JLON,JLEV,3)&
       & +ZKOZO(JLON,JLEV,7)*ZDELPS(JLON,JLEV)/RG&
       & - ZLIGHT(JLON)*(0.5_JPRB+SIGN(0.5_JPRB,ZTEMP-PT(JLON,JLEV)))&
       & *TEECL2*ZKOZO(JLON,JLEV,8) ) 
      ENDIF
    ENDIF

    ZANEX(JLON,JLEV)=ZANEX(JLON,JLEV)/ZRAPP
    ZANEX(JLON,JLEV)=MAX(1.E-9_JPRB,ZANEX(JLON,JLEV))
    IF(JLEV < KLEV)&
     & ZCO3(JLON,JLEV+1)=ZCO3(JLON,JLEV)&
     & +PDELP(JLON,JLEV)*ZANEX(JLON,JLEV)/RG  
  ENDDO
ENDDO

!*   COMPUTATION OF THE ACTUAL TREND

DO JLEV=KTDIA,KLEV
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    PTENO3(JLON,JLEV)=(ZANEX(JLON,JLEV)-PO3(JLON,JLEV))/TSPHY
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('O3CHEM',1,ZHOOK_HANDLE)
END SUBROUTINE O3CHEM
