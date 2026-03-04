! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CULIGHT &
 & ( YDEPHY,  YGFL  ,  YDECUMF, &
 &   KIDIA ,  KFDIA ,  KLON,    KLEV,  PGAW, PGELAT, &
 &   PAP,     PAPH  ,  PAPHI,   PAPHIF,  LDLAND,&
 &   PT    ,  PLU   ,  PMFU,    PCAPE, &
 &   PFPLCL,  PFPLCN,  PQPFROZ, &
 &   LDCUM ,  KCBOT ,  KCTOP,&
! Outputs
 &   LDLINOX, PLIGH_TOT,  PLIGH_CTG,  PCTOPH, &
 &   PPRECMX, PICE,   PCDEPTH,  PWMFU,  PCHARGE) 

!    THIS ROUTINE CALCULATES LIGHTNING FLASH RATES.

!    P. Lopez     ECMWF   (03/2005)
!
!
!    PURPOSE.
!    --------
!    TO CALCULATE LIGHTNING FLASH RATES (FLASHES/KM2/DAY)
!    FROM INFORMATION COMING FROM THE CONVECTION SCHEME.
!    TOTAL AND CLOUD-TO-GROUND FLASH RATES ARE COMPUTED.
!    Partly inspired from Kurz and Grewe (2002).
! 
!    The following lightning parameterizations are available:
!
!      NLIMODE = 1  -> Lopez, 2005, version 1.
!                2  -> Lopez, 2005, version 2.
!                3  -> Lopez, 2005, version 3.
!                4  -> TM5 Meijer, 2001.
!                5  -> Price and Rind, 1993.
!                6  -> Lopez, 2015, latest version.
!                7  -> Grewe et al., 2001.
!                8  -> Allen and Prickering, 2002.

!    INTERFACE
!    ---------
!    THIS ROUTINE IS CALLED FROM *LIGHTNING_LAYER*.

!    METHOD.
!    --------  

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

!    *PGAW*         GAUSSIAN WEIGHTS Reduced ~ grid box area
!    *PGELAT*       Latitude (radians) 
!    *PAP*          PRESSURE ON FULL LEVELS                PA
!    *PAPH*         PRESSURE ON HALF LEVELS                PA
!    *PAPHI*        GEOPOTENTIAL ON HALF LEVELS            M2/S2
!    *PAPHIF*       GEOPOTENTIAL ON FULL LEVELS            M2/S2
!    *PT*           TEMPERATURE                            K
!    *PLU*          CONVECTIVE CLOUD CONDENSATE            KG/KG
!    *PMFU*         CONVECTIVE UPDRAUGHT MASS FLUX         KG/(SM2)
!    *PFPLCL*       CONVECTIVE RAIN FLUX                   KG/(SM2)
!    *PFPLCN*       CONVECTIVE SNOW FLUX                   KG/(SM2)
!    *PQPFROZ*      CONVECTIVE FROZEN PRECIP CONTENT       KG/KG

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDLAND*       LAND-SEA MASK 

!    INPUT PARAMETERS (INTEGER):

!    *KCBOT*       CONVECTIVE CLOUD BASE LEVEL     
!    *KCTOP*       CONVECTIVE CLOUD TOP LEVEL

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDLINOX*     GRID-POINT FLAG: .TRUE. FOR LIGHTNING NOx COMPUTATIONS IN CULINOX

!    OUTPUT PARAMETERS (REAL):

!    *PLIGH_TOT*   TOTAL LIGHTNING FLASH RATES               FL/KM2/DAY
!    *PLIGH_CTG*   CLOUD-TO-GROUND LIGHTNING FLASH RATES     FL/KM2/DAY
!    *PCTOPH*      CONVECTIVE CLOUD TOP HEIGHT               KM
!    *PPRECMX*     MAXIMUM CONVECTIVE PRECIPITATION IN THE   KG/(SM2)
!                  VERTICAL.
!    *PICE*        TOTAL CONVECTIVE CLOUD ICE CONTENT        KG/M2 
!    *PCDEPTH*     DEPTH OF CLOUD ABOVE FREEZING LEVEL       KM 
!    *PWMFU*       MEAN UPDRAUGHT MASS FLUX                  M/S
!    *PCHARGE*     ELECTRIC CHARGE PRODUCTION RATE             

!    EXTERNALS
!    ---------

!    MODIFICATIONS
!    -------------
!    J Flemming   ECMWF   (04/2010) ! Prince & Rind 1993 and Meijer 2001 added
!    P Bechtold 18/05/2012   Use RDAYI for 86400
!    J  Flemming 03/07/2013  Remove difference in CTG and IC lighting , introdruce new profile according to Ott et al., 2010
!    N. Semane+P.Bechtold 04/10/2012   add RPLRG for small planet 
!    P. Lopez 24/07/2015     New ECMWF lightning parameterization plus added Grewe (2001) and Allen-Pickering (2002).
!    P. Lopez 08/10/2015     Moved lightning NOx computations to CULINOX.
!    P. Lopez 25/03/2021     Extended usage to the case of explicit convection (i.e. at very high resolutions).
!
!----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    , ONLY : RG       ,RLVTT    ,RLSTT    ,RTT      ,&
                     & RD       ,RPI      ,RA       ,RDAYI
USE YOETHF    , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
                     & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
                     & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
                     & RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
USE YOMDYNCORE, ONLY : RPLRG
USE YOEPHY    , ONLY : TEPHY
USE YOM_YGFL  , ONLY : TYPE_GFLD
USE YOECUMF   , ONLY : TECUMF
! USE YOMLUN    , ONLY : NULOUT


IMPLICIT NONE

TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
TYPE(TECUMF)      ,INTENT(INOUT) :: YDECUMF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAPE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCN(KLON,0:KLEV)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQPFROZ(KLON,KLEV)  
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON)  
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON)  
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDLINOX(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLIGH_TOT(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLIGH_CTG(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTOPH(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPRECMX(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PICE(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCDEPTH(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWMFU(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCHARGE(KLON)

!             LOCAL STORAGE
!             ----- -------

INTEGER(KIND=JPIM) :: JK, JL, IFREEZ, IFREEZ25, IFREEZ15
INTEGER(KIND=JPIM) :: ILEVTOPG, ILEVBOTG

REAL(KIND=JPRB) :: ZGEO2KM, ZBASHG
REAL(KIND=JPRB) :: ZPCTG, Z1G, ZDEPTH, ZRHO, ZAREA_REF
REAL(KIND=JPRB) :: ZAREA(KLON)
REAL(KIND=JPRB) :: ZCDEPTH(KLON), ZPREC_CV(KLON), ZCBASEH(KLON)
REAL(KIND=JPRB) :: ZQGRAUP, ZQSNOW, ZBETA, ZPFROZEN, ZCOLDDPT, ZHBMAX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE   


#include "fcttre.func.h"

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CULIGHT',0,ZHOOK_HANDLE)

ASSOCIATE(NLIMODE=>YDEPHY%NLIMODE, NCHEM=>YGFL%NCHEM, LMFPEN=>YDECUMF%LMFPEN)

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!----------------------------------------------------------------------

ZHBMAX=1.8_JPRB
Z1G=1.0_JPRB/RG
ZGEO2KM=0.001_JPRB*Z1G

! Initializations
PLIGH_TOT(:)=0.0_JPRB
PLIGH_CTG(:)=0.0_JPRB
PCTOPH(:)=-99.0_JPRB
PPRECMX(:)=0.0_JPRB
PICE(:)=0.0_JPRB
PCDEPTH(:)=0.0_JPRB
PWMFU(:)=-99.0_JPRB
ZCDEPTH(:)=0.0_JPRB
PCHARGE(:)=0.0_JPRB
ZPREC_CV(:)=0.0_JPRB
ZCBASEH(:)=-99.0_JPRB

LDLINOX(:)=.FALSE.

! grid box size in km**2 
ZAREA(KIDIA:KFDIA)=1.0E-6_JPRB*4.0_JPRB*RPI*RA*RA*PGAW(KIDIA:KFDIA)
! see Allen and Pickering JGR, 2002, 2.0 x 2.5 reference grid box centred at 30 deg latitude.
!LOP:Incorrect ZAREA_REF=(20000.0/FLOAT(90))*(20000.0/FLOAT(72)) 
ZAREA_REF=(20015.1_JPRB/90.0_JPRB)*(17333.6_JPRB/72.0_JPRB) 

!----------------------------------------------------------------------
!     1.           CALCULATE LIGHTNING FLASH RATES
!----------------------------------------------------------------------

IF (NCHEM /= 0) THEN
  NLIMODE = 4
!  WRITE(NULOUT,*) 'CAMS-SPECIFIC CONFIGURATION CHANGE: NLIMODE HAS BEEN SET TO 4'
ENDIF

DO JL=KIDIA,KFDIA

  IF (LDCUM(JL) .AND. KCTOP(JL) >= 1 .AND. KCTOP(JL) <= KLEV) THEN

! Compute lightning occurrences only if cloud base height is below 
! 4 km and if convective precipitation is present at least at one level
! between cloud base and cloud top.

    ZBASHG =(PAPHIF(JL,KCBOT(JL))-PAPHIF(JL,KLEV))*ZGEO2KM
    PPRECMX(JL)=MAXVAL(PFPLCL(JL,KCTOP(JL):KCBOT(JL))+PFPLCN(JL,KCTOP(JL):KCBOT(JL)))

! naj is this cloud base <  4 km OK for TM5's and MOZART's parameterisations ?
    IF (ZBASHG < (4.0_JPRB/RPLRG) .AND. PPRECMX(JL) > 1.E-10_JPRB) THEN

      ! Set flag for grid-points where NOx computations should be performed in CULINOX.
      LDLINOX(JL)=.TRUE.

      IFREEZ=-1
      IFREEZ25=-1
      IFREEZ15=-1
      DO JK=KLEV,1,-1
        IF (PT(JL,JK) <= RTT .AND. IFREEZ == -1) IFREEZ=JK  ! below 0C
        IF (PT(JL,JK) <= RTT - 25.0_JPRB  .AND. IFREEZ25 == -1) IFREEZ25=JK ! below -25C
        IF (PT(JL,JK) <= RTT - 15.0_JPRB  .AND. IFREEZ15 == -1) IFREEZ15=JK ! below -15C
      ENDDO

      ! Total cloud depth
      PCTOPH(JL)=MAX((PAPHI(JL,KCTOP(JL))-PAPHI(JL,KCBOT(JL)))*ZGEO2KM,0.0_JPRB)

      ! Cold cloud depth
      PCDEPTH(JL)=MAX((PAPHI(JL,KCTOP(JL))-PAPHI(JL,IFREEZ))*ZGEO2KM,0.0_JPRB)
      PCDEPTH(JL)=MIN(PCDEPTH(JL),PCTOPH(JL))

      ! Diagnose updraft velocity for lightning formulae.

      IF ( NLIMODE <= 3 .OR. NLIMODE == 7 ) THEN
        PWMFU(JL)=0.0_JPRB
        ZDEPTH=0.0_JPRB
        DO JK=KCTOP(JL),KCBOT(JL)
          IF (PMFU(JL,JK) > 0.0_JPRB) THEN
            ZRHO=PAP(JL,JK)/(RD*PT(JL,JK))
            ZDEPTH=ZDEPTH+(PAPHI(JL,JK-1)-PAPHI(JL,JK))
            PWMFU(JL)=PWMFU(JL)+PMFU(JL,JK)*(PAPHI(JL,JK-1)-PAPHI(JL,JK))/ZRHO
          ENDIF
        ENDDO
        IF (ZDEPTH > 0.0_JPRB) THEN
          PWMFU(JL)=PWMFU(JL)/ZDEPTH
        ELSE
          PWMFU(JL)=-99.0_JPRB
        ENDIF
      ENDIF

      IF ( NLIMODE <= 3 ) THEN
        ! Total convective cloud ice.
        DO JK=1,KLEV
          PICE(JL) = PICE(JL) &
           & + MAX(0.0_JPRB,PLU(JL,JK)*(1.0_JPRB-FOEALFCU(PT(JL,JK)))) &
           & * (PAPH(JL,JK)-PAPH(JL,JK-1))*Z1G
        ENDDO
      ENDIF

      IF ( NLIMODE == 6 ) THEN
        ILEVTOPG = MAX(IFREEZ25,KCTOP(JL))
        ILEVBOTG = MIN(IFREEZ,KCBOT(JL))

        ! Convective cloud base height.
        ZCBASEH(JL)=MAX(0.0_JPRB,(PAPHI(JL,KCBOT(JL))-PAPHI(JL,KLEV))*ZGEO2KM)

        ! Convective cloud depth.
        ZCDEPTH(JL)=MAX(0.0_JPRB,(PAPHI(JL,KCTOP(JL))-PAPHI(JL,KCBOT(JL)))*ZGEO2KM)

        ! Land/sea-dependent factor for the partitioning of frozen precipitation into graupel and snow.
        IF (LDLAND(JL)) THEN
          ZBETA = 0.70_JPRB
        ELSE
          ZBETA = 0.45_JPRB
        ENDIF

        ! Electric charge production rate.
        DO JK=ILEVTOPG,ILEVBOTG
          ! Fraction of frozen precipitation flux in the form of graupel inside this layer.
          IF (LMFPEN) THEN
            ! With parameterized convection, first convert input frozen precipitation flux into content
            ! and then partition into graupel and snow.
            ZRHO = PAP(JL,JK)/(RD*PT(JL,JK))
            ZPFROZEN = MAX(0.0_JPRB, PFPLCN(JL,JK))
            ZQGRAUP = ZBETA * ZPFROZEN / (ZRHO * 3.0_JPRB)
            ZQSNOW = (1.0_JPRB - ZBETA) * ZPFROZEN / (ZRHO * 0.5_JPRB)
          ELSE
            ! With explicit convection, partition input frozen precipitation content into 
            ! graupel and snow.
            ZQGRAUP = ZBETA * PQPFROZ(JL,JK)
            ZQSNOW = (1.0_JPRB - ZBETA) * PQPFROZ(JL,JK)  
          ENDIF
          PCHARGE(JL) = PCHARGE(JL) &
           & + ZQGRAUP * ( MAX(0.0_JPRB, PLU(JL,JK)) + ZQSNOW ) &
           & * (PAPH(JL,JK)-PAPH(JL,JK-1)) * Z1G
        ENDDO

      ENDIF

      IF ( NLIMODE == 8 ) THEN
        ! For Allen and Prickering (2002) param.
        ! Total convective surface precipitation (in mm/day; capped at 90 mm/day).
        ZPREC_CV(JL) = MIN((PFPLCL(JL,KLEV) + PFPLCN(JL,KLEV)) * 86400.0_JPRB, 90.0_JPRB) 
      ENDIF

! Lightning flash rates.

      IF ( NLIMODE == 1 ) THEN
        ! PL seems TechMemo code.
        IF (PCDEPTH(JL) > 1.0_JPRB .AND. (PCTOPH(JL)-PCDEPTH(JL)) > 0.5_JPRB .AND. &
          & PWMFU(JL) > 0.0_JPRB) THEN
          IF (LDLAND(JL)) THEN
            PLIGH_TOT(JL)=3.1E-01_JPRB*PICE(JL)*PWMFU(JL)*PCDEPTH(JL)
          ELSE
            PLIGH_TOT(JL)=1.3E-02_JPRB*PICE(JL)*PWMFU(JL)*PCDEPTH(JL)
          ENDIF
        ENDIF

      ELSEIF ( NLIMODE == 2 ) THEN
        ! PL was commented out. 
        IF (PCDEPTH(JL) > 1.0_JPRB .AND. (PCTOPH(JL)-PCDEPTH(JL)) > 0.7_JPRB .AND. &
          & PCTOPH(JL) > 5.0_JPRB .AND. PWMFU(JL) > 0.05_JPRB) THEN
          IF (LDLAND(JL)) THEN
            PLIGH_TOT(JL)=6.5E-01_JPRB*PICE(JL)*PWMFU(JL)*PCDEPTH(JL)
          ELSE
            PLIGH_TOT(JL)=8.0E-02_JPRB*PICE(JL)*PWMFU(JL)*PCDEPTH(JL)
          ENDIF
        ENDIF

      ELSEIF ( NLIMODE == 3 ) THEN
        ! PL version in code as is. 
        IF (PCDEPTH(JL) > 1.0_JPRB .AND. (PCTOPH(JL)-PCDEPTH(JL)) > 0.7_JPRB .AND. &
          & PCTOPH(JL) > 5.0_JPRB .AND. PWMFU(JL) > 0.05_JPRB) THEN
          IF (LDLAND(JL)) THEN
            PLIGH_TOT(JL)=4.5E-03_JPRB*(PICE(JL)*PCDEPTH(JL))**2.0_JPRB*PWMFU(JL)**0.5_JPRB
          ELSE
            PLIGH_TOT(JL)=6.6E-04_JPRB*(PICE(JL)*PCDEPTH(JL))**2.0_JPRB*PWMFU(JL)**0.5_JPRB
          ENDIF
        ENDIF

      ELSEIF ( NLIMODE == 4 ) THEN
        ! TM5 Meijer 2001.
        IF ( PCTOPH(JL) > 5.0_JPRB .AND. IFREEZ25 > 0 ) THEN
          IF (LDLAND(JL)) THEN
            PLIGH_TOT(JL)= 5.0_JPRB * 4.0E6_JPRB * PFPLCL(JL,KLEV) * 1.0E-9_JPRB * RDAYI 
          ELSE
            PLIGH_TOT(JL)= 0.1_JPRB * (5.0_JPRB * 4.0E6_JPRB * PFPLCL(JL,KLEV) * 1.0E-9_JPRB * RDAYI) 
          ENDIF
        ENDIF

      ELSEIF ( NLIMODE == 5 ) THEN
        ! Price and Rind 1993 as coded in mozart/ Base/mo_hook.F90 
        ! Flash rate in flashes / min -> convert to /(day*km2)
        IF ( IFREEZ > 0 ) THEN
          IF (LDLAND(JL)) THEN
             PLIGH_TOT(JL) = 3.44E-5_JPRB * PCTOPH(JL)**4.9_JPRB * (24.0_JPRB * 60.0_JPRB / ZAREA(JL)) 
          ELSE
             PLIGH_TOT(JL) = 6.4E-4_JPRB * PCTOPH(JL)**1.73_JPRB * (24.0_JPRB * 60.0_JPRB / ZAREA(JL))  
          ENDIF
        ENDIF  

      ELSEIF ( NLIMODE == 6 ) THEN
        ! New lightning parameterization (Ph. Lopez, 2015; flash density in flashes/km2/day).
        IF (ZCDEPTH(JL) > 4.0_JPRB .AND. PCAPE(JL) > 175._JPRB .AND. ZCBASEH(JL) > 0.1_JPRB .AND. &
          & KCTOP(JL) < IFREEZ15 .AND. PCDEPTH(JL) > 1.0_JPRB) THEN
          PLIGH_TOT(JL) = 36.0_JPRB * MIN(1.0_JPRB, (SQRT(ZAREA(JL)) / 30.0_JPRB)**0.15_JPRB) &
                      & * PCHARGE(JL) * SQRT(PCAPE(JL)) * MIN(ZCBASEH(JL),ZHBMAX)**2
        ENDIF

      ELSEIF ( NLIMODE == 7 ) THEN
        ! Grewe et al. (2001). 
        ! Flash rate in flashes/min  -> convert to flashes/(day*km2).
        IF (PWMFU(JL) > 0.0_JPRB) THEN
          PLIGH_TOT(JL)= 1.54E-5_JPRB * (PWMFU(JL) * SQRT(PCTOPH(JL) * 1000.0_JPRB))**4.9_JPRB &
                     & * (1440.0_JPRB / ZAREA(JL)) 
        ENDIF

      ELSEIF ( NLIMODE == 8 ) THEN
        ! Allen and Pickering (2002). 
        ! Note: This parameterization actually computes cloud-to-ground lightning.
        IF ( ZPREC_CV(JL) > 7.0_JPRB ) THEN
          IF (LDLAND(JL)) THEN
            PLIGH_CTG(JL) = 3.75E-2_JPRB - 4.76E-2_JPRB * ZPREC_CV(JL) + 5.41E-3_JPRB * ZPREC_CV(JL)**2 &
                        & + 3.21E-4_JPRB * ZPREC_CV(JL)**3 - 2.93E-6_JPRB * ZPREC_CV(JL)**4 
          ELSE
            PLIGH_CTG(JL) = 5.23E-2_JPRB - 4.80E-2_JPRB * ZPREC_CV(JL) + 5.45E-3_JPRB * ZPREC_CV(JL)**2 &
                        & + 3.68E-5_JPRB * ZPREC_CV(JL)**3 - 2.42E-7_JPRB * ZPREC_CV(JL)**4 
          ENDIF
          PLIGH_CTG(JL) = MAX(PLIGH_CTG(JL),0.0_JPRB)
          ! Cloud-to-ground flash rate in flashes/min  -> convert to flashes/(day*km2).
          PLIGH_CTG(JL) = PLIGH_CTG(JL) * (1440.0_JPRB / ZAREA_REF)
          ! Convert to total lightning using Price and Rind (1994).
          ZCOLDDPT = MAX(PCDEPTH(JL),5.5_JPRB)
          IF ( ZCOLDDPT < 14.0_JPRB ) THEN
            ZPCTG = 1.0_JPRB / &
             &   (0.021_JPRB*(ZCOLDDPT**4) & 
             &  - 0.648_JPRB*(ZCOLDDPT**3) &
             &  + 7.49_JPRB *(ZCOLDDPT**2) &
             &  - 36.54_JPRB*ZCOLDDPT &
             &  + 64.09_JPRB)
            ZPCTG = MAX(ZPCTG,0.0_JPRB)
          ELSE
            ZPCTG = 0.02_JPRB 
          ENDIF
          ! Total flash rate in flashes/(day*km2).
          PLIGH_TOT(JL) = PLIGH_CTG(JL) / ZPCTG
        ENDIF
      ENDIF 

! Cloud-to-ground lightning frequency as a function of depth of cloud
! cold sector (if higher than 5.5 km and less than 14 km).
! Note: if using Allen and Pickering (2002), CTG should already be available.

      IF ( NLIMODE /= 8 ) THEN
        IF ( PCDEPTH(JL) > 5.5_JPRB .AND. PCDEPTH(JL) < 14.0_JPRB ) THEN
          ZPCTG = 1.0_JPRB / &
             &   (0.021_JPRB*(PCDEPTH(JL)**4.0_JPRB) & 
             &  - 0.648_JPRB*(PCDEPTH(JL)**3.0_JPRB) &
             &  + 7.49_JPRB*(PCDEPTH(JL)**2.0_JPRB) &
             &  - 36.54_JPRB*PCDEPTH(JL) &
             &  + 64.09_JPRB)
          PLIGH_CTG(JL) = ZPCTG * PLIGH_TOT(JL)
        ELSE
          PLIGH_CTG(JL) = PLIGH_TOT(JL) 
        ENDIF
      ENDIF

    ENDIF
  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CULIGHT',1,ZHOOK_HANDLE)

END SUBROUTINE CULIGHT
