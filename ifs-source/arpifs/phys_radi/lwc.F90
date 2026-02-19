! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LWC &
 & ( YDML_PHY_RAD,YDEPHLI,KIDIA , KFDIA, KLON  , KLEV, KTOPC,&
 & PBINT , PBSUI, PCLFR, PCLDLD, PCLDLU,&
 & PCNTRB, PEMIT, PFLUC,&
 & PFLUX                              &
 & )  

!**** *LWC*   - LONGWAVE RADIATION, CLOUD EFFECTS

!     PURPOSE.
!     --------
!           INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
!           RADIANCES

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PBINT  : (KLON,KLEV+1)       ; HALF LEVEL PLANCK FUNCTION
! PBSUI  : (KLON)              ; SURFACE PLANCK FUNCTION
! PCLFR  : (KLON,KLEV)         ; CLOUD FRACTIONAL COVER
! PCLDLD : (KLON,KLEV)         ; DOWNWARD EFFECTIVE CLOUD FRACTION
! PCLDLU : (KLON,KLEV)         ; UPWARD EFFECTIVE CLOUD FRACTION
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PEMIT  : (KLON)              ; SURFACE TOTAL LW EMISSIVITY
! PFLUC  : (KLON,2,KLEV+1)     ; CLEAR-SKY LW RADIATIVE FLUXES
!     ==== OUTPUTS ===
! PFLUX  : (KLON,2,KLEV+1)     ; TOTAL SKY LW RADIATIVE FLUXES :
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
!          2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
!          3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
!     CLOUDS

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        JJ Morcrette 97-04-18   Cleaning
!        JJMorcrette 01-02-16 Hogan & Illingworth (2001)'s mixed overlap
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Janiskova 24-Jan-2006 Raisanen overlap assumption 
!                               & optimization for cloudy part
!        M.Janiskova 02-Mar-2012 bug correction for IKL
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEPHLI   , ONLY : TEPHLI

IMPLICIT NONE

TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
TYPE(MODEL_PHYSICS_RADIATION_TYPE),INTENT(INOUT):: YDML_PHY_RAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOPC
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBINT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBSUI(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLFR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDLD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCNTRB(KLON,KLEV+1,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIT(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFLUC(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUX(KLON,2,KLEV+1) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZCLEAR(KLON)            , ZCLOUD(KLON)&
 & ,  ZCLM(KLON,KLEV+1,KLEV+1), ZDNF(KLON,KLEV+1,KLEV+1)&
 & ,  ZFD(KLON)               , ZFU(KLON)&
 & ,  ZUPF(KLON,KLEV+1,KLEV+1)  
REAL(KIND=JPRB) :: ZCLFR(KLON,KLEV)
REAL(KIND=JPRB) :: ZACLFR(KLON,KLEV),ZCLDLD(KLON,KLEV),ZCLDLU(KLON,KLEV)
REAL(KIND=JPRB) :: ZEMIS(KLON,KLEV)
REAL(KIND=JPRB) :: ZSUM(KLON)
REAL(KIND=JPRB) :: ZTMP1(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZTMP3(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZTMP4(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZTMP5(KLON,KLEV+1)

INTEGER(KIND=JPIM) :: IKCP1, IKM1, IKP1, IMAXC, JCLOUD,&
 & JK, JK1, JK2, JKJ, JL, IKL, IJK 

REAL(KIND=JPRB) :: ZALPHA1, ZCFRAC, ZNUM, ZDEN
REAL(KIND=JPRB) :: ZTR1, ZTR2, ZCADJ
REAL(KIND=JPRB) :: ZDIV1, ZDIV2, ZDIV3, ZDIV4, ZDIV5, ZDIV6
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.     INITIALIZATION
!                 --------------

!100  CONTINUE

!      print *,' Enter LWC '
IF (LHOOK) CALL DR_HOOK('LWC',0,ZHOOK_HANDLE)
ASSOCIATE(RA1OVLP=>YDML_PHY_RAD%YREOVLP%RA1OVLP, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, LRAISANEN=>YDEPHLI%LRAISANEN, &
 & NOVLP=>YDML_PHY_RAD%YRERAD%NOVLP, &
 & REPCLC=>YDML_PHY_RAD%YRERDI%REPCLC)
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFLUX(JL,1,JK) = PFLUC(JL,1,JK)
    PFLUX(JL,2,JK) = PFLUC(JL,2,JK)
  ENDDO
ENDDO

!GM*******
IMAXC=KLEV
!GM*******

!Experimentation with reduced ITOPC 
!ITOPC=KLEV-10
!ITOPC=KLEV-20
!ITOPC=KLEV

!     ------------------------------------------------------------------

!*         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
!                  ---------------------------------------

!*         2.0     INITIALIZE TO CLEAR-SKY FLUXES
!                  ------------------------------

!200  CONTINUE

DO JK = 1,KTOPC+1
  DO JL = KIDIA,KFDIA
    ZDNF(JL,1,JK)=PFLUC(JL,2,JK)
  ENDDO
ENDDO

DO JK = 1,KLEV+1
  DO JL = KIDIA,KFDIA
    ZUPF(JL,1,JK) =PFLUC(JL,1,JK)
  ENDDO
ENDDO

!      print *,' LWC after Initialisation to clear-sky fluxes'

!*         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
!                  ----------------------------------------------

!210  CONTINUE

DO JCLOUD = 1 , KTOPC
  IKCP1=JCLOUD+1

!*         2.1.1   ABOVE THE CLOUD
!                  ---------------

!2110 CONTINUE

  DO JK=IKCP1,KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZFU(JL)=0.0_JPRB
    ENDDO

    IF (JK  >  IKCP1) THEN
      DO JKJ=IKCP1,IKM1
        DO JL = KIDIA,KFDIA
          ZFU(JL) = ZFU(JL) + PCNTRB(JL,JKJ,JK)
        ENDDO
      ENDDO
    ENDIF

    DO JL = KIDIA,KFDIA
      ZUPF(JL,IKCP1,JK)=PBINT(JL,JK)-ZFU(JL)
    ENDDO
  ENDDO

!*         2.1.2   BELOW THE CLOUD
!                  ---------------

!2120 CONTINUE

  DO JK=1,JCLOUD
    IKP1=JK+1
    DO JL = KIDIA,KFDIA
      ZFD(JL)=0.0_JPRB
    ENDDO

    IF (JK  <  JCLOUD) THEN
      DO JKJ=IKP1,JCLOUD
        DO JL = KIDIA,KFDIA
          ZFD(JL) = ZFD(JL) + PCNTRB(JL,JKJ,JK)
        ENDDO
      ENDDO
    ENDIF

    DO JL = KIDIA,KFDIA
      ZDNF(JL,IKCP1,JK)=-PBINT(JL,JK)-ZFD(JL)
    ENDDO
  ENDDO

ENDDO
!      print *,' LWC after 213: Fluxes for unity emissivity'

!*         2.2     CLOUD COVER MATRIX
!                  ------------------

!*    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
!     HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1

!220  CONTINUE

! matrix notation for ZCLM is not standard, but:
!    z11  z21  z31  ...
!    z12  z22  z32  ...
!    z13  z23  z33  ...
!     :    :    :    :

!DO JK = KTOPC,KLEV+1
!  DO JL = KIDIA,KFDIA
!    ZCLM(JL,KTOPC,JK) = 0.0_JPRB
!    ZCLM(JL,KTOPC+1,JK) = 0.0_JPRB
!  ENDDO
!ENDDO

DO JK2 = KTOPC , KLEV+1
  DO JK1 = KTOPC , KLEV+1
    DO JL = KIDIA,KFDIA
      ZCLM(JL,JK1,JK2) = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

IF (.NOT.LRAISANEN) THEN
  DO JK = KTOPC+2,KLEV
    DO JL = KIDIA,KFDIA
      ZCLM(JL,JK,JK+1) = 0.0_JPRB
    ENDDO
  ENDDO
ENDIF

! Securities

!DO JK = 1,KLEV
!  IKL = KLEV + 1 - JK
!  DO JL = KIDIA,KFDIA
!    IF (PCLFR(JL,IKL) <= REPCLC) THEN
!      ZCLDLD(JL,JK) = 0.0_JPRB
!      ZCLDLU(JL,JK) = 0.0_JPRB
!    ELSE
!      ZCLDLD(JL,JK) = PCLDLD(JL,JK)
!      ZCLDLU(JL,JK) = PCLDLU(JL,JK)
!    ENDIF
!  ENDDO
!ENDDO

IF (KTOPC /= KLEV) THEN
  DO JK = 1,KTOPC
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      ZACLFR(JL,IKL) = PCLFR(JL,IKL)

      ZCLDLD(JL,JK) = PCLDLD(JL,JK)
      ZCLDLU(JL,JK) = PCLDLU(JL,JK)
    ENDDO
  ENDDO
  DO JK = KTOPC+1,KLEV
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      ZACLFR(JL,IKL) = REPCLC

      ZCLDLD(JL,JK) = 0.0_JPRB
      ZCLDLU(JL,JK) = 0.0_JPRB
    ENDDO
  ENDDO
ELSE
  DO JK = 1,KLEV
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      ZACLFR(JL,IKL) = PCLFR(JL,IKL)

      ZCLDLD(JL,JK) = PCLDLD(JL,JK)
      ZCLDLU(JL,JK) = PCLDLU(JL,JK)
    ENDDO
  ENDDO
ENDIF

IF (LRAISANEN) THEN
  DO JK = 1,KLEV
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      IF (PCLDLU(JL,JK) <= REPCLC .OR. PCLDLD(JL,JK) <= REPCLC) THEN
        ZACLFR(JL,IKL) = REPCLC
      ELSE
        ZACLFR(JL,IKL) = PCLFR(JL,IKL)
      ENDIF
    ENDDO
  ENDDO
ENDIF

!      print *,' LWC after Initialisation CC matrix'

!*         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
!                  ------------------------------------------

!240  CONTINUE

IF (LRAISANEN) THEN
!* maximum-random (Raisanen 1998)

  ZEMIS(:,:)  = 0.0_JPRB
  DO JK = 1,KLEV
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      ZCLFR (JL,IKL) = ZACLFR (JL,IKL)
      IF (ZCLFR(JL,IKL) <= REPCLC) THEN
        ZCLFR(JL,IKL) = REPCLC
      ENDIF
      IF (ZCLFR(JL,IKL) >= 1.0_JPRB-REPCLC) THEN
        ZCLFR(JL,IKL) = 1.0_JPRB-REPCLC
      ENDIF
      ZDIV1 = 1.0_JPRB / ZCLFR(JL,IKL)
      ZEMIS(JL,JK) = ZCLDLU(JL,JK) * ZDIV1
    ENDDO
  ENDDO

  DO JK = 1,KTOPC-1
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      IF (ZCLDLU(JL,JK+1) <= REPCLC) THEN
        ZTMP1(JL,JK) = 1.0_JPRB/REPCLC
      ELSE
        ZTMP1(JL,JK) = 1.0_JPRB/ZCLDLU(JL,JK+1)
      ENDIF
      IF (ZCLFR(JL,IKL) <= ZCLFR(JL,IKL-1)) THEN
        ZCADJ  = ZCLFR(JL,IKL)
      ELSE
        ZCADJ  = ZCLFR(JL,IKL-1)
      ENDIF      
!* Modified for regularization 
      IF (LPHYLIN .AND. ZCLFR(JL,IKL-1) > 0.999_JPRB) THEN
        ZCLFR(JL,IKL-1) = 0.999_JPRB
      ENDIF
!* End of modification
      ZTMP3(JL,JK) = 1.0_JPRB/(1.0_JPRB-ZCLFR(JL,IKL-1))
      ZTMP4(JL,JK) = ZCADJ*(1.0_JPRB-ZEMIS(JL,JK+1))
      ZTMP5 (JL,JK) = ZCLFR(JL,IKL)-ZCADJ
    ENDDO
  ENDDO

  DO JK1 = 2 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZCLOUD (JL) = ZCLDLU(JL,JK1-1)
      ZSUM (JL)  = ZCLOUD(JL)
      ZCLM (JL,JK1-1,JK1) = ZSUM(JL)
    ENDDO
    IJK=MIN(KTOPC+1,JK1)

    DO JK = IJK - 2 , 1 , -1
      IKL = KLEV+1 - JK
      DO JL = KIDIA,KFDIA
        ZTR1 = ZCLOUD(JL)*ZTMP1(JL,JK)
        ZTR2 = (1.0_JPRB-(ZSUM(JL)-ZCLOUD(JL))-ZTR1*ZCLFR(JL,IKL-1))&
         & * ZTMP3(JL,JK)  
        ZCLOUD(JL) = ZEMIS(JL,JK)&
         & * (ZTMP4(JL,JK)*ZTR1 + ZTMP5(JL,JK)*ZTR2)  
        ZSUM(JL)  = ZSUM(JL) + ZCLOUD(JL)
        ZCLM (JL,JK,JK1) = ZSUM(JL)
      ENDDO
    ENDDO
  ENDDO

ELSE

  DO JK1 = 2 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZCLEAR(JL)=1.0_JPRB
      ZCLOUD(JL)=0.0_JPRB
    ENDDO

    IJK=MIN(KTOPC,JK1)
    DO JK = IJK - 1 , 1, -1
      ZALPHA1=RA1OVLP(KLEV+1-JK)

      DO JL = KIDIA,KFDIA
        IF ((NOVLP==1).OR.(NOVLP==6).OR.(NOVLP==8)) THEN
!* maximum-random    Geleyn and Hollingsworth 1979
          IF (ZCLDLU(JL,JK) >= ZCLOUD(JL)) THEN
            ZNUM = ZCLDLU (JL,JK)
          ELSE
            ZNUM = ZCLOUD(JL)
          ENDIF
          IF (ZCLOUD(JL) <= 1.0_JPRB-REPCLC) THEN
            ZDEN = ZCLOUD(JL)
          ELSE
            ZDEN= 1.0_JPRB-REPCLC
          ENDIF
          ZDIV2 = 1.0_JPRB/(1.0_JPRB-ZDEN)
          ZCLEAR(JL)=ZCLEAR(JL)*(1.0_JPRB-ZNUM)*ZDIV2
          ZCLM(JL,JK,JK1) = 1.0_JPRB - ZCLEAR(JL)
          ZCLOUD (JL) = ZCLDLU(JL,JK)
        ELSEIF ((NOVLP==2).OR.(NOVLP==7)) THEN
!* maximum
          IF (ZCLDLU(JL,JK) >= ZCLOUD(JL)) THEN
            ZCLOUD (JL) = ZCLDLU(JL,JK)
          ENDIF
          ZCLM (JL,JK,JK1) = ZCLOUD(JL)
        ELSEIF ((NOVLP == 3).OR.(NOVLP==5)) THEN
!* random
          ZCLEAR(JL) = ZCLEAR(JL)*(1.0_JPRB - ZCLDLU(JL,JK))
          ZCLOUD(JL) = 1.0_JPRB - ZCLEAR(JL)
          ZCLM (JL,JK,JK1) = ZCLOUD(JL)
        ELSEIF (NOVLP == 4) THEN
!** Hogan & Illingworth (2001)
          IF (ZCLDLU(JL,JK) >= ZCLOUD(JL)) THEN
            ZNUM = ZCLDLU(JL,JK)
          ELSE
            ZNUM = ZCLOUD(JL)
          ENDIF
          IF (ZCLOUD(JL) <= 1.0_JPRB-REPCLC) THEN
            ZDEN = ZCLOUD(JL)
          ELSE
            ZDEN= 1.0_JPRB-REPCLC
          ENDIF
          ZDIV3 = 1.0_JPRB/(1.0_JPRB-ZDEN)
          ZCLEAR(JL) = ZCLEAR(JL)&
           & * (ZALPHA1*(1.0_JPRB-ZNUM)*ZDIV3&
           & + (1.0_JPRB-ZALPHA1)*(1.0_JPRB-ZCLDLU(JL,JK)) ) 
          ZCLM(JL,JK,JK1) = 1.0_JPRB - ZCLEAR(JL)
          ZCLOUD(JL) = ZCLDLU(JL,JK)
        ENDIF
      ENDDO
    ENDDO

  ENDDO

ENDIF

!      print *,' LWC after 244: CC below level of calculation'

!*         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
!                  ------------------------------------------

!250  CONTINUE

IF (LRAISANEN) THEN
!* maximum-random (Raisanen 1998)

  ZEMIS(:,:)  = 0.0_JPRB

  DO JK = 1,KLEV
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      ZDIV4 = 1.0_JPRB / ZCLFR(JL,IKL)
      ZEMIS(JL,JK) = ZCLDLD(JL,JK) * ZDIV4
    ENDDO
  ENDDO

  DO JK = 2,KTOPC
    IKL = KLEV + 1 - JK
    DO JL = KIDIA,KFDIA
      IF (ZCLDLD(JL,JK-1) <= REPCLC) THEN
        ZTMP1(JL,JK) = 1.0_JPRB/REPCLC
      ELSE
        ZTMP1(JL,JK) = 1.0_JPRB/ZCLDLD(JL,JK-1)
      ENDIF
      IF (ZCLFR(JL,IKL) <= ZCLFR(JL,IKL+1)) THEN
        ZCADJ  = ZCLFR(JL,IKL)
      ELSE
        ZCADJ  = ZCLFR(JL,IKL+1)
      ENDIF     
!* Modified for regularization
      IF (LPHYLIN .AND. ZCLFR(JL,IKL+1) > 0.999_JPRB) THEN
        ZCLFR(JL,IKL+1) = 0.999_JPRB
      ENDIF
!* End of modification
      ZTMP3(JL,JK) = 1.0_JPRB/(1.0_JPRB-ZCLFR(JL,IKL+1))
      ZTMP4(JL,JK) = ZCADJ*(1.0_JPRB-ZEMIS(JL,JK-1))
      ZTMP5 (JL,JK) = ZCLFR(JL,IKL)-ZCADJ
    ENDDO
  ENDDO

  DO JK1 = 1 ,KTOPC 
    DO JL = KIDIA,KFDIA
      ZCLOUD (JL) = ZCLDLD(JL,JK1)
      ZSUM (JL)  = ZCLOUD(JL)
      ZCLM (JL,JK1,JK1) = ZSUM(JL)
    ENDDO

    DO JK = JK1+1 , KTOPC
      IKL = KLEV + 1 - JK
      DO JL = KIDIA,KFDIA
        ZTR1 = ZCLOUD(JL)*ZTMP1(JL,JK)
        ZTR2 = (1.0_JPRB-(ZSUM(JL)-ZCLOUD(JL))-ZTR1*ZCLFR(JL,IKL+1))&
         & * ZTMP3(JL,JK)  
        ZCLOUD(JL) = ZEMIS(JL,JK)&
         & * (ZTMP4(JL,JK)*ZTR1 + ZTMP5(JL,JK)*ZTR2)  
        ZSUM(JL)  = ZSUM(JL) + ZCLOUD(JL)
        ZCLM (JL,JK,JK1)  = ZSUM(JL)
      ENDDO
    ENDDO
  ENDDO

ELSE

  DO JK1 = 1 , KTOPC
    DO JL = KIDIA,KFDIA
      ZCLEAR(JL)=1.0_JPRB
      ZCLOUD(JL)=0.0_JPRB
    ENDDO

    DO JK = JK1 , KTOPC
      ZALPHA1=RA1OVLP(KLEV+1-JK)

      DO JL = KIDIA,KFDIA
        IF ((NOVLP == 1).OR.(NOVLP==6).OR.(NOVLP==8)) THEN
!* maximum-random
          IF (ZCLDLD(JL,JK) >= ZCLOUD(JL)) THEN
            ZNUM = ZCLDLD (JL,JK)
          ELSE
            ZNUM = ZCLOUD(JL)
          ENDIF
          IF (ZCLOUD(JL) <= 1.0_JPRB-REPCLC) THEN
            ZDEN = ZCLOUD(JL)
          ELSE
            ZDEN= 1.0_JPRB-REPCLC
          ENDIF
          ZDIV5 = 1.0_JPRB/(1.0_JPRB-ZDEN)
          ZCLEAR(JL)=ZCLEAR(JL)*(1.0_JPRB-ZNUM)*ZDIV5
          ZCLM(JL,JK,JK1) = 1.0_JPRB - ZCLEAR(JL)
          ZCLOUD (JL) = ZCLDLD(JL,JK)
        ELSEIF ((NOVLP == 2).OR.(NOVLP==7)) THEN
!* maximum
          IF (ZCLDLD(JL,JK) >= ZCLOUD(JL)) THEN
            ZCLOUD (JL) = ZCLDLD(JL,JK)
          ENDIF
          ZCLM (JL,JK,JK1) = ZCLOUD(JL)
        ELSEIF ((NOVLP == 3).OR.(NOVLP==5)) THEN
!* random
          ZCLEAR(JL) = ZCLEAR(JL)*(1.0_JPRB - ZCLDLD(JL,JK))
          ZCLOUD(JL) = 1.0_JPRB - ZCLEAR(JL)
          ZCLM (JL,JK,JK1) = ZCLOUD(JL)
        ELSEIF (NOVLP == 4) THEN
!** Hogan & Illingworth (2001)      
          IF (ZCLDLD(JL,JK) >= ZCLOUD(JL)) THEN
            ZNUM = ZCLDLD(JL,JK)
          ELSE
            ZNUM = ZCLOUD(JL)
          ENDIF
          IF (ZCLOUD(JL) <= 1.0_JPRB-REPCLC) THEN
            ZDEN = ZCLOUD(JL)
          ELSE
            ZDEN= 1.0_JPRB-REPCLC
          ENDIF
          ZDIV6 = 1.0_JPRB/(1.0_JPRB-ZDEN)
          ZCLEAR(JL)=ZCLEAR(JL)&
           & * (ZALPHA1*(1.0_JPRB-ZNUM)*ZDIV6&
           & + (1.0_JPRB-ZALPHA1)*(1.0_JPRB-ZCLDLD(JL,JK)) ) 
          ZCLM(JL,JK,JK1) = 1.0_JPRB - ZCLEAR(JL)
          ZCLOUD(JL) = ZCLDLD(JL,JK)
        ENDIF
      ENDDO
    ENDDO
  ENDDO

ENDIF

!      print *,' LWC after 254: CC above level of calculation'

!*         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
!                  ----------------------------------------------

!300  CONTINUE

!*         3.1     DOWNWARD FLUXES
!                  ---------------

!310  CONTINUE

DO JL = KIDIA,KFDIA
  PFLUX(JL,2,KLEV+1) = 0.0_JPRB
ENDDO

DO JK1 = 1, KTOPC

!*                 CONTRIBUTION FROM CLEAR-SKY FRACTION

  DO JL = KIDIA,KFDIA
    ZFD (JL) = (1.0_JPRB - ZCLM(JL,KTOPC,JK1)) * ZDNF(JL,1,JK1)

!*                 CONTRIBUTION FROM ADJACENT CLOUD

    ZFD(JL) = ZFD(JL) + ZCLM(JL,JK1,JK1) * ZDNF(JL,JK1+1,JK1)
  ENDDO

!*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

  DO JK = JK1, KTOPC-1 
    DO JL = KIDIA,KFDIA
      ZCFRAC = ZCLM(JL,JK+1,JK1) - ZCLM(JL,JK,JK1)
      ZFD(JL) =  ZFD(JL) + ZCFRAC * ZDNF(JL,JK+2,JK1)
    ENDDO
  ENDDO

  DO JL = KIDIA,KFDIA
    PFLUX(JL,2,JK1) = ZFD (JL)
  ENDDO

ENDDO
!      print *,' LWC after 317: Downward fluxes'

!*         3.2     UPWARD FLUX AT THE SURFACE
!                  --------------------------

!320  CONTINUE

DO JL = KIDIA,KFDIA
  PFLUX(JL,1,1) = PEMIT(JL)*PBSUI(JL)-(1.0_JPRB-PEMIT(JL))*PFLUX(JL,2,1)
ENDDO

!*         3.3     UPWARD FLUXES
!                  -------------

!330  CONTINUE

DO JK1 = 2 , KLEV+1

!*                 CONTRIBUTION FROM CLEAR-SKY FRACTION

  DO JL = KIDIA,KFDIA
    ZFU (JL) = (1.0_JPRB - ZCLM(JL,1,JK1)) * ZUPF(JL,1,JK1)

!!*                 CONTRIBUTION FROM ADJACENT CLOUD
!    ZFU(JL) =  ZFU(JL) + ZCLM(JL,JK1,JK1-1) * ZUPF(JL,JK1,JK1)

  ENDDO

!*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

  IJK=MIN(KTOPC+1,JK1)
  DO JK = 2 , IJK-1
    DO JL = KIDIA,KFDIA
      ZCFRAC = ZCLM(JL,JK-1,JK1) - ZCLM(JL,JK,JK1)
      ZFU(JL) =  ZFU(JL) + ZCFRAC * ZUPF(JL,JK  ,JK1)
    ENDDO
  ENDDO

!*                 CONTRIBUTION FROM ADJACENT CLOUD

  DO JL = KIDIA,KFDIA
    ZFU(JL) =  ZFU(JL) + ZCLM(JL,JK1-1,IJK) * ZUPF(JL,IJK,IJK)
  ENDDO

  DO JL = KIDIA,KFDIA
    PFLUX(JL,1,JK1) = ZFU (JL)
  ENDDO

ENDDO
!      print *,' LWC after 337: Upward fluxes'

!-----------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LWC',1,ZHOOK_HANDLE)
END SUBROUTINE LWC
