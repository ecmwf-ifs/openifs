! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
!******************************************************************************
SUBROUTINE RRTM_TAUMOL14 (KIDIA,KFDIA,KLEV,P_TAU,&
 & P_TAUAERL,P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,P_FORFRAC,K_INDFOR,K_JP,K_JT,K_JT1,&
 & P_COLCO2,K_LAYTROP,P_SELFFAC,P_SELFFRAC,K_INDSELF,PFRAC)  

!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)

!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF

!     MODIFICATIONS.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        NEC           25-Oct-2007 Optimisations
!        JJMorcrette 20110613 flexible number of g-points
!        ABozzo 201306 updated to rrtmg v4.85
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NGS13  ,NG14
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA14, ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,SELFREF,FORREF

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAU(KIDIA:KFDIA,JPGPT,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,KLEV,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLCO2(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRAC(KIDIA:KFDIA,JPGPT,KLEV) 

INTEGER(KIND=JPIM),INTENT(IN)   :: K_INDFOR(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)   :: P_FORFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: P_FORFRAC(KIDIA:KFDIA,KLEV) 
! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS, INDF, JLAY
INTEGER(KIND=JPIM) :: JLON
REAL(KIND=JPRB) :: ZTAUFOR,ZTAUSELF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum 
! and foreign continuum is interpolated (in temperature) separately.  

ASSOCIATE(NFLEVG=>KLEV)
IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL14',0,ZHOOK_HANDLE)

DO JLAY = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (JLAY <= K_LAYTROP(JLON)) THEN
      IND0 = ((K_JP(JLON,JLAY)-1)*5+(K_JT(JLON,JLAY)-1))*NSPA(14) + 1
      IND1 = (K_JP(JLON,JLAY)*5+(K_JT1(JLON,JLAY)-1))*NSPA(14) + 1
      INDS = K_INDSELF(JLON,JLAY)
      INDF = K_INDFOR(JLON,JLAY)

!-- DS_990714  
!-- jjm20110728 re-establishing the loop instead of specified IG to allow a flexible number of NG14
      DO IG = 1, NG14
!      IG=1
        ZTAUSELF = P_SELFFAC(JLON,JLAY)* (SELFREF(INDS,IG) + P_SELFFRAC(JLON,JLAY) * &
          &       (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
        ZTAUFOR = P_FORFAC(JLON,JLAY) * (FORREF(INDF,IG) + P_FORFRAC(JLON,JLAY) * &
          &       (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 

        P_TAU(JLON,NGS13+IG,JLAY) = P_COLCO2(JLON,JLAY) *&
         & (P_FAC00(JLON,JLAY) * ABSA(IND0  ,IG) +&
         & P_FAC10(JLON,JLAY) * ABSA(IND0+1,IG) +&
         & P_FAC01(JLON,JLAY) * ABSA(IND1  ,IG) +&
         & P_FAC11(JLON,JLAY) * ABSA(IND1+1,IG)) &
         & + ZTAUSELF + ZTAUFOR &
         & + P_TAUAERL(JLON,JLAY,14)  
        PFRAC(JLON,NGS13+IG,JLAY) = FRACREFA(IG)
      ENDDO
!-- jjm20110728
!-- DS_990714  
    ENDIF

    IF (JLAY > K_LAYTROP(JLON)) THEN
      IND0 = ((K_JP(JLON,JLAY)-13)*5+(K_JT(JLON,JLAY)-1))*NSPB(14) + 1
      IND1 = ((K_JP(JLON,JLAY)-12)*5+(K_JT1(JLON,JLAY)-1))*NSPB(14) + 1
!-- DS_990714  
!-- jjm20110728 re-establishing the loop instead of specified IG to allow a flexible number of NG14
      DO IG = 1, NG14
!      IG=1
        P_TAU(JLON,NGS13+IG,JLAY) = P_COLCO2(JLON,JLAY) *&
         & (P_FAC00(JLON,JLAY) * ABSB(IND0  ,IG) +&
         & P_FAC10(JLON,JLAY) * ABSB(IND0+1,IG) +&
         & P_FAC01(JLON,JLAY) * ABSB(IND1  ,IG) +&
         & P_FAC11(JLON,JLAY) * ABSB(IND1+1,IG)) &
         & + P_TAUAERL(JLON,JLAY,14)  
        PFRAC(JLON,NGS13+IG,JLAY) = FRACREFB(IG)
      ENDDO
!-- jjm20110728
!-- DS_990714  
    ENDIF
  ENDDO
ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL14',1,ZHOOK_HANDLE)

END ASSOCIATE
END SUBROUTINE RRTM_TAUMOL14
