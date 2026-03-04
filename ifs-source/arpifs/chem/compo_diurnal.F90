! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COMPO_DIURNAL &
  &( YDRIP,KIDIA, KFDIA, KLON, CDTYPE, PGELAM, PGELAT, PSCALE, PAMPLITUDE, PHOURPEAK, PLSM, PTANDEC, PBASELINE )

!*** * COMPO_DIURNAL* - COMPUTE SCALING FACTOR BASED ON DIURNAL CYCLE

!**   INTERFACE.
!     ----------
!          *COMPO_DIURNAL* IS CALLED FROM *AER_SRC*, *AER_DRYDEP*, *AER_SO2SO4*
!          and
!          *CHEM_INITFLUX*

!* INPUTS:
!  -------
! PGELAM     - Geographic longitude in radians
! PGELAT     - Geographic latidude in radians
! CDTYPE     - type of function
! PAMPLITUDE - Amplitude of the scaling factor
! PHOURPEAK  - hour of peak (local solar time)
! PLSM       - Land Sea Mask - if present, only apply diurnal cycle where >=0.1
! PTANDEC    - Tangent of solar declination 
! PBASELINE  - Baseline scale factor (e.g. night-time)

!* OUTPUTS:
!  --------
! PSCALE         DIURNAL CYCLE SCALING FACTOR

!     AUTHOR.
!     -------
!        SAMUEL REMY, from chemical sources
!        Johannes Flemming: general for all BB emission tracers, PHOURPEAK   


!     MODIFICATIONS.
!     --------------
!        J. Flemming 20180517:  GFAS diurnal function and KTYP 
!        J. Flemming 20190710:  VOC diurnal cycle including day length and anthropogenic diurnal cycle  
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    ,ONLY : RPI
USE YOMLUN    ,ONLY : NULOUT
USE YOMRIP    ,ONLY : TRIP

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP)        ,INTENT(IN)  :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON, KIDIA, KFDIA
CHARACTER(LEN=*)  ,INTENT(IN)  :: CDTYPE
REAL(KIND=JPRB)   ,INTENT(IN)  :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PGELAT(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSCALE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PAMPLITUDE
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PHOURPEAK 
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PTANDEC  
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PBASELINE
!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZLTRAD(KLON), ZLTH(KLON), ZFAC(KLON), ZSIG_VOC(KLON), ZARG(KLON), ZDAYLENGTH(KLON)
REAL(KIND=JPRB) :: ZFAC_GFAS
REAL(KIND=JPRB) :: ZTEST, ZLTH_TEST
REAL(KIND=JPRB) :: ZSCALETEST(24)
INTEGER(KIND=JPIM) :: IHOUR(KLON)
INTEGER(KIND=JPIM) :: JL, JHOUR
! widths of day time peak / sigma in hours   
REAL(KIND=JPRB) , PARAMETER :: ZSIG_GFAS=2.0_JPRB 

! Legacy explicit anthropogenic diurnal cycles
! ZANTEMI removed because it's unused
REAL(KIND=JPRB) , DIMENSION(24), PARAMETER :: & 
   & ZANTEMI_NOX = (/0.25,0.17,0.15,0.14,0.17,0.30,0.96,1.75,1.77,1.46,1.27,1.19,1.25,1.33,1.32,1.41,1.72,1.79, & 
   & 1.48,1.20,0.94,0.84,0.73,0.43/)
REAL(KIND=JPRB) , DIMENSION(24), PARAMETER :: & 
   & ZANTEMI_CO = (/0.32,0.28,0.27,0.27,0.29,0.42,1.09,1.62,1.66,1.52,1.32,1.17,1.15,1.17,1.14,1.16,1.30,1.41, & 
   & 1.44,1.38,1.20,1.13,0.88,0.43 /)
REAL(KIND=JPRB) , DIMENSION(24), PARAMETER :: & 
   & ZANTEMI_SOG = (/0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.9,2.7,3.05,3.05,3.05,3.05,2.7,1.3,0.4, & 
   & 0.2,0.2,0.2,0.2,0.2,0.2 /)
! ZANTEMI_NH3 removed because it's the same as ZCG_AGSL
! ZANTEMI_SO2 removed because it's the same as ZCG_ENE

! CAMS-GLOB explicit diurnal cycles as provided by Marc Guevara under CAMS_81
!   Note: SHP, FEF, TNR and SWD are all uniform so omitted.
!   Note: AGS and AGL are identical, so a single AGSL is included here.
REAL(KIND=JPRB), DIMENSION(24), PARAMETER :: &
 & ZCG_ENE = (/ 0.79,   0.72,   0.72,   0.71,   0.74,   0.80,   0.92,   1.08,   1.19,   1.22,   1.21,   1.21,     &
 &              1.17,   1.15,   1.14,   1.13,   1.10,   1.07,   1.04,   1.02,   1.02,   1.01,   0.96,   0.88  /)
REAL(KIND=JPRB), DIMENSION(24), PARAMETER :: &
 & ZCG_TRO = (/ 0.19,   0.09,   0.06,   0.05,   0.09,   0.22,   0.86,   1.84,   1.86,   1.41,   1.24,   1.20,     &
 &              1.32,   1.44,   1.45,   1.59,   2.03,   2.08,   1.51,   1.06,   0.74,   0.62,   0.61,   0.44  /)
REAL(KIND=JPRB), DIMENSION(24), PARAMETER :: &
 & ZCG_RES = (/ 0.38,   0.36,   0.36,   0.36,   0.37,   0.5 ,   1.19,   1.53,   1.57,   1.56,   1.35,   1.16,     &
 &              1.07,   1.06,   1.00,   0.98,   0.99,   1.12,   1.41,   1.52,   1.39,   1.35,   1.00,   0.42  /)
REAL(KIND=JPRB), DIMENSION(24), PARAMETER :: &
 & ZCG_IND = (/ 0.75,   0.75,   0.78,   0.82,   0.88,   0.95,   1.02,   1.09,   1.16,   1.22,   1.28,   1.30,     &
 &              1.22,   1.24,   1.25,   1.16,   1.08,   1.01,   0.95,   0.90,   0.85,   0.81,   0.78,   0.75  /)
REAL(KIND=JPRB), DIMENSION(24), PARAMETER :: &
 & ZCG_SLV = (/ 0.50,   0.35,   0.20,   0.10,   0.10,   0.20,   0.75,   1.25,   1.40,   1.50,   1.50,   1.50,     &
 &              1.50,   1.50,   1.50,   1.50,   1.50,   1.40,   1.25,   1.10,   1.00,   0.90,   0.80,   0.70  /)
REAL(KIND=JPRB), DIMENSION(24), PARAMETER :: &
 & ZCG_AGSL= (/ 0.60,   0.60,   0.60,   0.60,   0.60,   0.65,   0.75,   0.90,   1.10,   1.35,   1.45,   1.60,     &
 &              1.65,   1.75,   1.70,   1.55,   1.35,   1.1 ,   0.90,   0.75,   0.65,   0.60,   0.60,   0.60  /)
REAL(KIND=JPRB), DIMENSION(24), PARAMETER :: &
 & ZCG_AWB = (/ 0.06,   0.06,   0.06,   0.07,   0.07,   0.07,   0.20,   0.20,   0.20,   1.82,   1.82,   1.82,     &
 &              3.39,   3.39,   3.39,   1.68,   1.68,   1.68,   0.56,   0.56,   0.56,   0.22,   0.22,   0.22  /)

LOGICAL :: LLTEST 



!-----------------------------------------------------------------------
#include "abor1.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COMPO_DIURNAL',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
LLTEST=.FALSE.

! The integral of the function is suposed to be 1 * 2RPI, i.e. the daily average is not chnaged  
PSCALE(KIDIA:KFDIA) = 1.0_JPRB 

ZLTRAD(KIDIA:KFDIA) = YDRIP%RWSOVR + PGELAM(KIDIA:KFDIA)
! hour 
ZLTH(KIDIA:KFDIA) = MODULO(ZLTRAD(KIDIA:KFDIA)*12.0_JPRB/RPI, 24.0_JPRB)

WHERE(ZLTH(KIDIA:KFDIA) == 24.0_JPRB)
  ! Edge case: in single precision at least, MODULO(..., 24.0) can
  ! occasionally return exactly 24.0, rather than a value which is
  ! strictly less than this as expected. Ensure that we don't end
  ! up setting IHOUR=25 in this case!
  IHOUR(KIDIA:KFDIA) = 24_JPIM
ELSEWHERE
  IHOUR(KIDIA:KFDIA) = INT(ZLTH(KIDIA:KFDIA), KIND=JPIM) + 1_JPIM
ENDWHERE

IF (ANY(IHOUR(KIDIA:KFDIA) < 1)) THEN
  WRITE(NULOUT,*) "ZLTRAD", ZLTRAD
  WRITE(NULOUT,*) "ZLTH", ZLTH
  WRITE(NULOUT,*) "IHOUR", IHOUR
  CALL ABOR1("COMPO_DIURNAL: IHOUR HAS VALUES BELOW 1") 
ENDIF
IF (ANY(IHOUR(KIDIA:KFDIA) > 24)) THEN
  WRITE(NULOUT,*) "ZLTRAD", ZLTRAD
  WRITE(NULOUT,*) "ZLTH", ZLTH
  WRITE(NULOUT,*) "IHOUR", IHOUR
  CALL ABOR1("COMPO_DIURNAL: IHOUR HAS VALUES ABOVE 24") 
ENDIF

SELECT CASE (CDTYPE) 
  CASE ('Uniform', 'CAMS_GLOB_SHP', 'CAMS_GLOB_FEF', 'CAMS_GLOB_TNR', 'CAMS_GLOB_SWD')
    PSCALE(KIDIA:KFDIA) = 1.0_JPRB

  ! Hard-coded explicit cycles defined as hourly arrays
  CASE('AnthroNOx')     ; PSCALE(KIDIA:KFDIA) = ZANTEMI_NOX(IHOUR(KIDIA:KFDIA))
  CASE('AnthroCO')      ; PSCALE(KIDIA:KFDIA) = ZANTEMI_CO(IHOUR(KIDIA:KFDIA))
  CASE('AnthroSOG')     ; PSCALE(KIDIA:KFDIA) = ZANTEMI_SOG(IHOUR(KIDIA:KFDIA))
  CASE('CAMS_GLOB_TRO') ; PSCALE(KIDIA:KFDIA) = ZCG_TRO(IHOUR(KIDIA:KFDIA))
  CASE('CAMS_GLOB_RES') ; PSCALE(KIDIA:KFDIA) = ZCG_RES(IHOUR(KIDIA:KFDIA))
  CASE('CAMS_GLOB_IND') ; PSCALE(KIDIA:KFDIA) = ZCG_IND(IHOUR(KIDIA:KFDIA))
  CASE('CAMS_GLOB_SLV') ; PSCALE(KIDIA:KFDIA) = ZCG_SLV(IHOUR(KIDIA:KFDIA))
  CASE('CAMS_GLOB_AWB') ; PSCALE(KIDIA:KFDIA) = ZCG_AWB(IHOUR(KIDIA:KFDIA))

  CASE('CAMS_GLOB_ENE') ; PSCALE(KIDIA:KFDIA) = ZCG_ENE(IHOUR(KIDIA:KFDIA))

  CASE('CAMS_GLOB_AGL', 'CAMS_GLOB_AGS')
                          PSCALE(KIDIA:KFDIA) = ZCG_AGSL(IHOUR(KIDIA:KFDIA))

  ! Functional/parametric cycles below

  CASE ('Sine') ! symmentrical cosine   
    IF (.NOT. (PRESENT(PAMPLITUDE) .AND. PRESENT(PHOURPEAK))) THEN
      CALL ABOR1("COMPO_DIURNAL: CDTYPE='Sine' needs PAMPLITUDE and PHOURPEAK") 
    ENDIF 
    PSCALE(KIDIA:KFDIA) = 1.0_JPRB + COS(ZLTRAD(KIDIA:KFDIA) - PHOURPEAK/12.0_JPRB * RPI) * PAMPLITUDE

  CASE('ZeroNight') ! "zero night", 18 to 6 LT 
    IF (.NOT. (PRESENT(PHOURPEAK))) THEN
      CALL ABOR1("COMPO_DIURNAL: CDTYPE='ZeroNight' needs PHOURPEAK") 
    ENDIF 
    PSCALE(KIDIA:KFDIA) = RPI * MAX(0._JPRB, COS(ZLTRAD(KIDIA:KFDIA) - PHOURPEAK/12.0_JPRB * RPI) )

  CASE('GFAS') ! GFAS diurnal profile, Kaiser et al. 2009 TM596, Fig 1
    IF (.NOT. (PRESENT(PHOURPEAK) .AND. PRESENT(PBASELINE))) THEN
      CALL ABOR1("COMPO_DIURNAL: CDTYPE='GFAS' needs PHOURPEAK and PBASELINE") 
    ENDIF 
    ZFAC_GFAS = (1.0_JPRB - PBASELINE) * (24.0_JPRB / (ZSIG_GFAS * SQRT( 2.0_JPRB * RPI)))   
    PSCALE(KIDIA:KFDIA) = PBASELINE + ZFAC_GFAS * EXP(-0.5_JPRB * ((ZLTH(KIDIA:KFDIA) - PHOURPEAK)/ZSIG_GFAS)**2.0_JPRB)  

  CASE('VOC') ! VOC profile  
    IF (.NOT. (PRESENT(PBASELINE) .AND. PRESENT(PHOURPEAK) .AND. PRESENT(PTANDEC))) THEN
      CALL ABOR1("COMPO_DIURNAL: CDTYPE='VOC' needs PHOURPEAK, PBASELINE and PTANDEC") 
    ENDIF
    ! calculate half day length 
    ZARG(KIDIA:KFDIA) = -PTANDEC * ATAN(PGELAT(KIDIA:KFDIA))   
    WHERE ( ZARG(KIDIA:KFDIA) < -0.98_JPRB )
      ZDAYLENGTH(KIDIA:KFDIA) = 24.0_JPRB  
    ELSEWHERE ( ZARG(KIDIA:KFDIA) > 0.98_JPRB )
      ZDAYLENGTH(KIDIA:KFDIA) = 0.0_JPRB  
    ELSEWHERE
      ZDAYLENGTH(KIDIA:KFDIA) = 24.0_JPRB * ACOS(-PTANDEC * ATAN(PGELAT(KIDIA:KFDIA))) / RPI 
    ENDWHERE

    WHERE ( ZDAYLENGTH(KIDIA:KFDIA) < 24.0_JPRB .AND. ZDAYLENGTH(KIDIA:KFDIA) > 0.0_JPRB )
      ZSIG_VOC(KIDIA:KFDIA) = 0.25_JPRB * ZDAYLENGTH(KIDIA:KFDIA)   
      ZFAC(KIDIA:KFDIA)  = (1.0_JPRB - PBASELINE) * (24.0_JPRB / (ZSIG_VOC(KIDIA:KFDIA) * SQRT( 2.0_JPRB * RPI)))
      PSCALE(KIDIA:KFDIA) = PBASELINE &
                          & + ZFAC(KIDIA:KFDIA) &
                          &   * EXP(-0.5_JPRB * ((ZLTH(KIDIA:KFDIA) - PHOURPEAK)/ZSIG_VOC(KIDIA:KFDIA))**2.0_JPRB)
    ELSEWHERE
      PSCALE(KIDIA:KFDIA) = 1.0_JPRB 
    ENDWHERE

    IF (LLTEST) THEN 
      IF ( ZDAYLENGTH(KFDIA) < 24.0_JPRB .AND. ZDAYLENGTH(KFDIA) > 0.0_JPRB ) THEN   
        DO JHOUR=1,24   
          ! hour
          ZLTH_TEST = REAL(JHOUR-1, KIND=JPRB) 
          ZSCALETEST(JHOUR) = PBASELINE + ZFAC(KFDIA) * EXP(-0.5_JPRB * ((ZLTH_TEST - PHOURPEAK)/ZSIG_VOC(KFDIA))**2.0_JPRB)
        ENDDO 
        ZTEST=SUM(ZSCALETEST)/24._JPRB
        IF (ZTEST > 1.01_JPRB .OR. ZTEST < 0.99_JPRB) THEN
          CALL ABOR1(" SCALES NOT 1 COMPO_DIURNAL  ") 
        ENDIF 
       ENDIF
     ENDIF 

  CASE DEFAULT 
    CALL ABOR1("COMPO_DIURNAL: CDTYPE='"//TRIM(CDTYPE)//"' not recognised") 
END SELECT 

IF (PRESENT(PLSM) ) THEN 
  WHERE (PLSM(KIDIA:KFDIA) < 0.1_JPRB) PSCALE(KIDIA:KFDIA) = 1.0_JPRB
ENDIF 

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COMPO_DIURNAL',1,ZHOOK_HANDLE)
END SUBROUTINE COMPO_DIURNAL
