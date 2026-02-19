! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CULINOX &
 & ( YDCHEM   , KIDIA    ,  KFDIA ,  KLON,    KLEV, &
 &   PGELAT   , PAPH     ,  PAPHI ,  PAPHIF, &
 &   LDLAND   , LDLINOX  , &
 &   PT       , KCTOP    , &
 &   PLIGH_TOT, PLIGH_CTG, &
! Outputs.
 &   PNOEMI, PNOEMI2D) 

!    THIS ROUTINE COMPUTES NOx PRODUCTION BY LIGHTNING.

!    J Flemming   ECMWF   (04/2010) ! NO emissions  
!
!
!    PURPOSE.
!    --------
!    TO CALCULATE NOx PRODUCTION BY LIGHTNING.
! 

!    INTERFACE
!    ---------
!    THIS ROUTINE IS CALLED FROM *LIGHTNING_LAYER*.

!    METHOD.
!    -------
!    

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

!    *PGELAT*       Latitude (radians) 
!    *PAPH*         PRESSURE ON HALF LEVELS                PA
!    *PAPHI*        GEOPOTENTIAL ON HALF LEVELS            M2/S2
!    *PAPHIF*       GEOPOTENTIAL ON FULL LEVELS            M2/S2
!    *PT*           TEMPERATURE                            K
!    *PLIGH_TOT*    TOTAL LIGHTNING FLASH RATES            FL/KM2/DAY
!    *PLIGH_CTG*    CLOUD-TO-GROUND LIGHTNING FLASH RATES  FL/KM2/DAY

!    INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND-SEA MASK 
!    *LDLINOX*      GRID-POINT FLAG: .TRUE. FOR LIGHTNING NOx COMPUTATIONS 

!    INPUT PARAMETERS (INTEGER):

!    *KCTOP*       CONVECTIVE CLOUD TOP LEVEL

!    OUTPUT PARAMETERS (REAL):

!    *PNOEMI*       3D NO Emissions in                      kg/m*2/s
!    *PNOEMI2D*     2D total NO Emissions in                kg/m*2/s

!    EXTERNALS
!    ---------

!    MODIFICATIONS
!    -------------
!    P. Lopez 08/10/2015 : Separate routine from lightning computations.
!
!----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM,   JPRB
USE YOMHOOK   , ONLY : LHOOK,  DR_HOOK, JPHOOK

USE YOMCST    , ONLY : RG,  RTT,  RPI,  RNAVO,  RDAYI
USE YOMCHEM   , ONLY : TCHEM


IMPLICIT NONE

TYPE(TCHEM)       ,INTENT(INOUT):: YDCHEM
INTEGER(KIND=JPIM),INTENT(IN)   :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)   :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)   :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)   :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PGELAT(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PAPH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PAPHI(KLON,0:KLEV) 
LOGICAL           ,INTENT(IN)   :: LDLAND(KLON) 
LOGICAL           ,INTENT(IN)   :: LDLINOX(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)   :: KCTOP(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PLIGH_TOT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)   :: PLIGH_CTG(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)  :: PNOEMI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)  :: PNOEMI2D(KLON)
 
!             LOCAL STORAGE
!             ----- -------

INTEGER(KIND=JPIM) :: JK, JL, IFREEZ15, IX(1), IL1, IL2, IHKM, ISURF 

REAL(KIND=JPRB) :: Z1G, ZNO, ZNO_CG, ZNO_IC, ZDELPC, ZEN_IC, &
                 & ZEN_CG, ZTEST, ZDELH, ZASL, ZSCALE 
REAL(KIND=JPRB) :: ZLAT  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE 

LOGICAL         :: LLTEST, LLCSHAPE

! Percent of NO emissions in levels 0-17 km for SupTrop(1), Mid-Lat(2), tropics continental(3), tropics marine(4) 
! Ott et al, table2
REAL(KIND=JPRB),DIMENSION(4,17) , PARAMETER :: ZNOPROF=RESHAPE( & 
&  (/ 1.0_JPRB,  2.4_JPRB,   0.2_JPRB,  0.6_JPRB , & 
&     2.1_JPRB,  5.0_JPRB,   0.5_JPRB,  1.5_JPRB , &
&     3.9_JPRB,  7.4_JPRB,   0.6_JPRB,  2.9_JPRB , &
&     5.8_JPRB,  9.3_JPRB,   1.4_JPRB,  4.3_JPRB , &
&     7.7_JPRB, 10.6_JPRB,   2.7_JPRB,  5.4_JPRB , &
&     9.3_JPRB, 11.4_JPRB,   4.0_JPRB,  6.7_JPRB , &
&    10.5_JPRB, 11.5_JPRB,   5.0_JPRB,  7.7_JPRB , &
&    11.0_JPRB, 11.0_JPRB,   6.2_JPRB,  8.5_JPRB , &
&    11.0_JPRB,  9.9_JPRB,   8.6_JPRB,  9.6_JPRB , &
&    10.4_JPRB,  8.3_JPRB,  10.3_JPRB, 10.2_JPRB , &
&     9.2_JPRB,  6.3_JPRB,  11.6_JPRB, 10.5_JPRB , &
&     7.5_JPRB,  4.2_JPRB,  12.4_JPRB, 10.2_JPRB , &
&     5.5_JPRB,  2.2_JPRB,  12.7_JPRB,  8.2_JPRB , &
&     3.4_JPRB,  0.5_JPRB,  12.4_JPRB,  6.5_JPRB , &
&     1.5_JPRB,  0.0_JPRB,   7.6_JPRB,  4.5_JPRB , &
&     0.2_JPRB , 0.0_JPRB,   3.0_JPRB,  2.2_JPRB , &
&     0.0_JPRB , 0.0_JPRB,   0.8_JPRB,  0.5_JPRB   /),(/4,17/))

!#include "abor1.intfb.h"

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CULINOX',0,ZHOOK_HANDLE)

ASSOCIATE(LCHEM_CSHAPE=>YDCHEM%LCHEM_CSHAPE, RCHEM_LINOX_SCALING=>YDCHEM%RCHEM_LINOX_SCALING)
LLTEST=.FALSE.
LLCSHAPE=LCHEM_CSHAPE
!LLCSHAPE=.TRUE.

! ZLI_SCAL now moved to RCHEM_LINOX_SCALING in YOMCHEM to allow
! local adjustment via namelist.

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!----------------------------------------------------------------------

Z1G=1.0_JPRB/RG

! Initializations
PNOEMI(:,:)=0.0_JPRB
PNOEMI2D(:)=0.0_JPRB

!----------------------------------------------------------------------
!     1.           CALCULATE LIGHTNING FLASH RATES
!----------------------------------------------------------------------

DO JL=KIDIA,KFDIA

  IF (LDLINOX(JL)) THEN

! naj is this cloud base <  4 km OK for TM5's and MOZART's parameterisations ?
    IFREEZ15=-1
    DO JK=KLEV,1,-1
      IF (PT(JL,JK) <= RTT - 15.0_JPRB  .AND. IFREEZ15 == -1) IFREEZ15=JK ! below -15C
    ENDDO

! NO emissions for TM5 scheme 
    IF ( PLIGH_CTG(JL) > 0.0_JPRB ) THEN
!     Flash energy - intra cloud 6.7e8 J, cloud to ground 6.7e8 J J/km2/day 
      ZEN_CG =  PLIGH_CTG(JL) * 6.7E8_JPRB       
      ZEN_IC = (PLIGH_TOT(JL) - PLIGH_CTG(JL) ) * 6.7E8_JPRB      

!     Flash NO production is 10e16 molecules NO/J 
!     #NO/km**2/day  
!     is this factor correct ??, Meijer says 1.0e16
      ZNO_CG = 1.0E17_JPRB * ZEN_CG 
      ZNO_IC = 1.0E17_JPRB * ZEN_IC 
!     Convert to kg / m2 / s 
      ZNO_CG = RCHEM_LINOX_SCALING* ( 0.0300061_JPRB  * ZNO_CG / RNAVO ) /( RDAYI * 1.0E6_JPRB ) 
      ZNO_IC = RCHEM_LINOX_SCALING* ( 0.0300061_JPRB  * ZNO_IC / RNAVO ) /( RDAYI * 1.0E6_JPRB )
      PNOEMI2D(JL) = ZNO_IC + ZNO_CG 
!     Distribute vertically    
!     DISTRIBUTION of LNOx over the COLUMNi - according to TM5
!       assume all IC-LNO and 70% of CG-LNOx betweem t=-15 and cloudtop;
!       assume                 10% of CG-LNOx between EARTH SURFACE and t=-15
!       assume                 20% of CG-LNOx in BOUNDARY LAYER
!         
! - LNO within one of these three regions is distributed proportional to the mass of each layer.
!                      
      IF ( LLCSHAPE ) THEN 
!       Distribute all IC-LNO and 70% of CG-LNOx betweem t=-15 and cloudtop;
 
        ZNO = ZNO_IC + 0.7_JPRB * ZNO_CG
!       some checks for IFREEZ15 == -1 or above KCTOP 
        IF ( IFREEZ15 < 0 ) IFREEZ15 = KCTOP(JL)
        IF ( IFREEZ15 < KCTOP(JL) ) IFREEZ15 = KCTOP(JL)
        ZDELPC=PAPH(JL,IFREEZ15 ) - PAPH(JL,KCTOP(JL) - 1) 
        DO JK= IFREEZ15, KCTOP(JL), -1 
          PNOEMI(JL, JK) = PNOEMI(JL, JK) + ZNO * (PAPH(JL,JK) - PAPH(JL,JK-1)) / ZDELPC
        ENDDO
!       distributing 10% of CG LNOx between EARTH SURFACE and t=-15
        ZNO = 0.1_JPRB * ZNO_CG
        ZDELPC = PAPH(JL,KLEV) - PAPH(JL,IFREEZ15-1) 
        DO JK= KLEV, IFREEZ15, -1 
          PNOEMI(JL, JK) =  PNOEMI(JL, JK) + ZNO * (PAPH(JL,JK ) - PAPH(JL,JK - 1)) / ZDELPC
        ENDDO
!       distributing 20% of CG LNOx in lowest 10 levels 
        ZNO = 0.2_JPRB * ZNO_CG
        ZDELPC = PAPH(JL,KLEV) - PAPH(JL,KLEV - 10 - 1) 
        DO JK = KLEV, KLEV-10, -1 
          PNOEMI(JL, JK) =  PNOEMI(JL, JK) + ZNO * (PAPH(JL,JK ) - PAPH(JL,JK - 1)) / ZDELPC
        ENDDO
      ELSE

! use Ott et al , 2010, distribution
! profile type according to latitude
        ZLAT=(180.0_JPRB/RPI)*PGELAT(JL)
! mid-lats
        IF ( ABS(ZLAT) >= 35.0_JPRB ) ISURF = 2_JPIM
! subtropics 
        IF ( ABS(ZLAT) > 20.0_JPRB .AND. ABS(ZLAT) < 35.0_JPRB ) ISURF = 1_JPIM
! tropics 
        IF ( ABS(ZLAT) <= 20.0_JPRB ) THEN
          ISURF = 4_JPIM 
          IF (LDLAND(JL)) ISURF = 3_JPIM
        ENDIF

! loop over 1km layers to 17 km
        IL1=KLEV
        ZSCALE=0.0_JPRB
!        ZTEST=0.0_JPRB  
        DO IHKM = 1, 17 
          ZASL = FLOAT(IHKM)*1000.0_JPRB + (PAPHI(JL,KLEV))*Z1G ! in m
          IX=MINLOC( ABS( PAPHIF(JL,1:KLEV)*Z1G - ZASL ))
          IL2=IX(1)
          ZDELH=PAPHI(JL,IL2-1) - PAPHI(JL,IL1)
          ZNO=ZNOPROF(ISURF, IHKM) *  PNOEMI2D(JL)*0.01_JPRB
          DO JK = IL2, IL1 
            ZSCALE = ZSCALE + ZNOPROF(ISURF, IHKM)*0.01_JPRB*(PAPHI(JL,JK-1) - PAPHI(JL,JK)) / ZDELH
            PNOEMI(JL, JK) = ZNO * (PAPHI(JL,JK-1 ) - PAPHI(JL,JK)) / ZDELH
!            ZTEST = ZTEST+ PNOEMI(JL, JK) 
!            print*,IL1, IL2, ZNOPROF(ISURF, IHKM)*0.01_JPRB, ZSCALE, PNOEMI2D(JL), PNOEMI(JL, JK), ZTEST   
          ENDDO
          IL1 = IL2 - 1_JPIM
        ENDDO
        PNOEMI(JL, 1:KLEV) = PNOEMI(JL, 1:KLEV) * ZSCALE  
      ENDIF
 
! test vertical sum
      IF (LLTEST) THEN
        ZTEST=0.0_JPRB  
        DO JK=1,KLEV
          ZTEST=ZTEST+ PNOEMI(JL, JK) 
        ENDDO
        IF (ABS( 2.0_JPRB *(ZTEST - PNOEMI2D(JL))/(ZTEST + PNOEMI2D(JL))) > 0.001_JPRB ) THEN
          WRITE(*,*) ' ERROR IN VERTICAL DISTRIBUTION OF LIGHTNING ',ZTEST, PNOEMI2D(JL)
!          CALL ABOR1(' ERROR IN VERTICAL DISTRIBUTION OF LIGHTNING ')
         ENDIF 
      ENDIF 

    ENDIF

  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CULINOX',1,ZHOOK_HANDLE)

END SUBROUTINE CULINOX
