! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUCCDIA &
 & (  YDERAD,YDEPHLI,YDEPHY,KIDIA,    KFDIA,    KLON,    KLEV,&
 & KSTEP,    KCBOT,    KCTOP,&
 & LDCUM,    PQU,      PLU,      PMFU,     PRAIN,&
 & PARPRC,   KTOPC,    KBASEC                   )  

!**** *CUCCDIA*- UPDATES PRECIPITAION, CLOUD BASE AND CLOUD TOP
!                FOR DIAGNOSTIC SCHEME FOR CONVECTIVE CLOUDS

!          M.TIEDTKE         E.C.M.W.F.    12/89

!**   INTERFACE.
!     ----------

!          *CUCCDIA* IS CALLED FROM *CUCALL*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL

!     INPUT PARAMETERS (LOGICAL)

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

!     INPUT PARAMETERS (REAL)

!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    *PARPRC*       ACCUMULATED PRECIPITATION AMMOUNT             KG/(M2*S)
!                   FOR RADIATION CALCULATION
!    *KTOPC*        CONVECTIVE CLOUD TOP LEVEL FOR RADIATION
!    *KBASEC*       CONVECTIVE CLOUD BASE LEVEL FOR RADIATION

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEPHY   , ONLY : TEPHY
USE YOERAD   , ONLY : TERAD
USE YOEPHLI  , ONLY : TEPHLI

IMPLICIT NONE

TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAIN(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PARPRC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTOPC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KBASEC(KLON) 
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IKB, IKT, JL

REAL(KIND=JPRB) :: ZDMFQ, ZNORMR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

!*    1.0          STORE CLOUD PARAMETERS FOR RADIATION CALCULATION
!                  -----------------------------------------------

!********************************
IF (LHOOK) CALL DR_HOOK('CUCCDIA',0,ZHOOK_HANDLE)
ASSOCIATE(LPHYLIN=>YDEPHLI%LPHYLIN, &
 & LECUMF=>YDEPHY%LECUMF, LERADI=>YDEPHY%LERADI, &
 & NRADFR=>YDERAD%NRADFR)
IF (.NOT.LPHYLIN.AND.NRADFR /= 1) THEN
!********************************      
  IF(LECUMF.AND.LERADI.AND.MOD(KSTEP+1,NRADFR) /= 0) THEN
    ZNORMR=1.0_JPRB/REAL(MAX(NRADFR-1,1),JPRB)
!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        KBASEC(JL)=MAX(KBASEC(JL),KCBOT(JL))
        KTOPC(JL) =MIN(KTOPC(JL) ,KCTOP(JL))
        IF(KTOPC(JL) == 1) KTOPC(JL)=KCTOP(JL)
        IKB=KCBOT(JL)
        IKT=KCTOP(JL)
        ZDMFQ=PMFU(JL,IKB)*MAX(PQU(JL,IKB)+PLU(JL,IKB)-PQU(JL,&
         & IKT),&
         & PQU(JL,IKB)*0.05_JPRB)  
        PARPRC(JL)=PARPRC(JL)+MAX(ZDMFQ,PRAIN(JL))*ZNORMR
      ENDIF
    ENDDO
  ENDIF
!********************************
ELSE
!********************************      
  IF(LECUMF) THEN
!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        KBASEC(JL)=KCBOT(JL)
        KTOPC(JL) =KCTOP(JL)
        IF(KTOPC(JL) == 1) KTOPC(JL)=KCTOP(JL)
        IKB=KCBOT(JL)
        IKT=KCTOP(JL)
        ZDMFQ=PMFU(JL,IKB)*MAX(PQU(JL,IKB)+PLU(JL,IKB)-PQU(JL,&
         & IKT),&
         & PQU(JL,IKB)*0.05_JPRB)  
        PARPRC(JL)=MAX(ZDMFQ,PRAIN(JL))
      ENDIF
    ENDDO
  ENDIF
!********************************        
ENDIF
!********************************      

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUCCDIA',1,ZHOOK_HANDLE)
END SUBROUTINE CUCCDIA
