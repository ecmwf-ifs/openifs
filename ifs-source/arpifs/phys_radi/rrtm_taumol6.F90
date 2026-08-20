! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
!----------------------------------------------------------------------------
SUBROUTINE RRTM_TAUMOL6 (KIDIA,KFDIA,KLEV,P_TAU,P_WX,&
 & P_TAUAERL,P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,P_FORFRAC,K_INDFOR,K_JP,K_JT,K_JT1,&
 & P_COLH2O,P_COLCO2,P_COLDRY,K_LAYTROP,P_SELFFAC,P_SELFFRAC,K_INDSELF,PFRAC,PMINORFRAC,KINDMINOR)  

!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 201306 updated to rrtmg v4.85
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE PARRRTM  , ONLY : JPBAND ,JPXSEC
USE YOERRTM  , ONLY : JPGPT  ,NG6   ,NGS5
USE YOERRTWN , ONLY : NSPA   
USE YOERRTA6 , ONLY : ABSA   ,KA_MCO2 ,CFC11ADJ , CFC12  ,&
 & FRACREFA,SELFREF,FORREF  
USE YOERRTRF, ONLY : CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAU(KIDIA:KFDIA,JPGPT,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_WX(KIDIA:KFDIA,JPXSEC,KLEV) ! Amount of trace gases
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLCO2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLDRY(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRAC(KIDIA:KFDIA,JPGPT,KLEV) 

INTEGER(KIND=JPIM),INTENT(IN)   :: K_INDFOR(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: P_FORFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: P_FORFRAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PMINORFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)   :: KINDMINOR(KIDIA:KFDIA,KLEV)

! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IND0(KLEV),IND1(KLEV),INDS(KLEV),INDF(KLEV),INDM(KLEV)

INTEGER(KIND=JPIM) :: IG, JLAY
INTEGER(KIND=JPIM) :: JLON

REAL(KIND=JPRB) :: ZADJFAC,ZADJCOLCO2(KIDIA:KFDIA,KLEV),ZRATCO2,ZCHI_CO2
REAL(KIND=JPRB) :: ZTAUFOR,ZTAUSELF,ZABSCO2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12


!     Compute the optical depth by interpolating in ln(pressure) and
!     temperature. The water vapor self- and foreign- continuum is interpolated
!     (in temperature) separately.  

ASSOCIATE(NFLEVG=>KLEV)
IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL6',0,ZHOOK_HANDLE)

DO JLAY = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (JLAY <= K_LAYTROP(JLON)) THEN
! In atmospheres where the amount of CO2 is too great to be considered
! a minor species, adjust the column amount of CO2 by an empirical factor 
! to obtain the proper contribution.
      ZCHI_CO2 = P_COLCO2(JLON,JLAY)/P_COLDRY(JLON,JLAY)
      ZRATCO2 = 1.E20_JPRB*ZCHI_CO2/CHI_MLS(2,K_JP(JLON,JLAY)+1)
      IF (ZRATCO2 > 3.0_JPRB) THEN
         ZADJFAC = 2.0_JPRB+(ZRATCO2-2.0_JPRB)**0.77_JPRB
         ZADJCOLCO2(JLON,JLAY) = ZADJFAC*CHI_MLS(2,K_JP(JLON,JLAY)+1)*P_COLDRY(JLON,JLAY)*1.E-20_JPRB
      ELSE
         ZADJCOLCO2(JLON,JLAY) = P_COLCO2(JLON,JLAY)
      ENDIF

      IND0(JLAY) = ((K_JP(JLON,JLAY)-1)*5+(K_JT(JLON,JLAY)-1))*NSPA(6) + 1
      IND1(JLAY) = (K_JP(JLON,JLAY)*5+(K_JT1(JLON,JLAY)-1))*NSPA(6) + 1
      INDS(JLAY) = K_INDSELF(JLON,JLAY)
      INDF(JLAY) = K_INDFOR(JLON,JLAY)
      INDM(JLAY) = KINDMINOR(JLON,JLAY)

!-- DS_000515  
!CDIR UNROLL=NG6
      DO IG = 1, NG6
!-- DS_000515  
        ZTAUSELF = P_SELFFAC(JLON,JLAY)* (SELFREF(INDS(JLAY),IG) + P_SELFFRAC(JLON,JLAY) * &
            &     (SELFREF(INDS(JLAY)+1,IG) - SELFREF(INDS(JLAY),IG)))
        ZTAUFOR = P_FORFAC(JLON,JLAY) * (FORREF(INDF(JLAY),IG) + P_FORFRAC(JLON,JLAY) * &
            &     (FORREF(INDF(JLAY)+1,IG) - FORREF(INDF(JLAY),IG))) 
        ZABSCO2 = KA_MCO2(INDM(JLAY),IG) + PMINORFRAC(JLON,JLAY) * &
            &     (KA_MCO2(INDM(JLAY)+1,IG) - KA_MCO2(INDM(JLAY),IG))

        P_TAU(JLON,NGS5+IG,JLAY) = P_COLH2O(JLON,JLAY) *&
         & (P_FAC00(JLON,JLAY) * ABSA(IND0(JLAY)  ,IG) +&
         & P_FAC10(JLON,JLAY) * ABSA(IND0(JLAY)+1,IG) +&
         & P_FAC01(JLON,JLAY) * ABSA(IND1(JLAY)  ,IG) +&
         & P_FAC11(JLON,JLAY) * ABSA(IND1(JLAY)+1,IG)) +&
         & ZTAUSELF + ZTAUFOR &
         & + P_WX(JLON,2,JLAY) * CFC11ADJ(IG)&
         & + P_WX(JLON,3,JLAY) * CFC12(IG)&
         & + ZADJCOLCO2(JLON,JLAY) * ZABSCO2 &
         & + P_TAUAERL(JLON,JLAY,6)  
        PFRAC(JLON,NGS5+IG,JLAY) = FRACREFA(IG)
      ENDDO
    ENDIF

!     Nothing important goes on above LAYTROP in this band.
!-- JJM_000517
    IF (JLAY > K_LAYTROP(JLON)) THEN
!CDIR UNROLL=NG6
      DO IG = 1, NG6
!-- JJM_000517
        P_TAU(JLON,NGS5+IG,JLAY) = 0.0_JPRB &
         & + P_WX(JLON,2,JLAY) * CFC11ADJ(IG)&
         & + P_WX(JLON,3,JLAY) * CFC12(IG)&
         & + P_TAUAERL(JLON,JLAY,6)  
        PFRAC(JLON,NGS5+IG,JLAY) = FRACREFA(IG)
      ENDDO
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL6',1,ZHOOK_HANDLE)

END ASSOCIATE
END SUBROUTINE RRTM_TAUMOL6
