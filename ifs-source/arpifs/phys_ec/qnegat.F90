! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE QNEGAT &
 & ( YDECND, KIDIA, KFDIA, KLON, KLEV,&
 & PTSPHY,&
 & PQM,  PQE,    PRESH,&
 & PFCQNG, PITEND        )  

!**** *QNEGAT* - AVOIDS NEGATIVE VALUES FOR SPECIFIC HUMIDITY.

!     B.RITTER     E.C.M.W.F.     29/04/83.

!     PURPOSE.
!     --------

!          THIS ROUTINE RESETS TO A VERY SMALL MINIMUM VALUE ALL THE
!     SPECIFIC HUMIDITIES SMALLER THAN THIS (IN PARTICULAR ALL THE
!     NEGATIVE VALUES THAT MIGHT HAVE APPEARED IN THE COURSE OF THE
!     TIME STEP).

!**   INTERFACE.
!     ----------

!      *QNEGAT* IS CALLED FROM *CALLPAR*.
!       THE ROUTINE TAKES ITS INPUT FROM Q0 OR Q0 AND ITS
!       TENDENCY (MAKING IT EFFECTIVELY TO BE Q1).
!       APPROPRIATE PRESSURE FIELD HAS TO BE PROVIDED AT THE
!       SAME TIME LEVEL.
!      IT RETURNS LOCAL TENDENCY CORRECTING NEGATIVE MOISTURE
!      AND ITS PROJECTION TO THE FLUX.
!      THE TEMPERATURE WILL NOT BE ADJUSTED AND SURFACE
!      FLUXES (EVAPORATION) WILL NOT BE MODIFIED.

!     METHOD.
!     -------

!          FOR EACH LAYER (FROM TOP TO BOTTOM) THE SPECIFIC HUMIDITY IS
!     COMPARED WITH ITS CRITICAL VALUE. IF NECESSARY THE ADJUSTMENT IS
!     PERFORMED AND THE NECESSARY MOISTURE IS TAKEN FROM THE LAYER
!     UNDERNEATH.

!     EXTERNALS.
!     ----------

!          NONE.

!     REFERENCE.
!     ----------

!          SEE RELEVANT PART OF THE DOCUMENTATION.

!     MODIFICATIONS.
!     --------------

!     F. Vana  19-Nov-2016   Total rewrite and bug fix.

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG
USE YOECND   , ONLY : TECND

IMPLICIT NONE

TYPE(TECND)       ,INTENT(INOUT) :: YDECND
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY                ! Time-step
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM(KLON,KLEV)        ! Moisture field 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQE(KLON,KLEV)        ! Tendency increment by q_neg corr
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESH(KLON,0:KLEV)    ! Half level pressure at appr. time
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCQNG(KLON,0:KLEV)   ! Pseudo-flux to correct for q<0.
REAL(KIND=JPRB), OPTIONAL  ,INTENT(IN) :: PITEND(KLON,KLEV)  ! Input tendency of q

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZDP(KLON,KLEV)
REAL(KIND=JPRB) :: ZFCQNG(KLON,0:KLEV) 
REAL(KIND=JPRB) :: ZRTSPHY, ZDQ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('QNEGAT',0,ZHOOK_HANDLE)

ASSOCIATE(REPQMI=>YDECND%REPQMI)
!     ------------------------------------------------------------------

!*       1.     ADJUSTMENT FOR ALL LAYERS.
!               ---------- --- --- -------

ZRTSPHY=1.0_JPRB/PTSPHY

!*       1.1      PRELIMINARY COMPUTATIONS FOR THE PRESSURE DIFFERENCE

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZDP(JL,JK)=PRESH(JL,JK)-PRESH(JL,JK-1)
  ENDDO
ENDDO


!*        1.2     INITIALIZATION OF THE PSEUDO FLUX COMPENSATING NEGATIVE HUMIDITY

IF (.NOT. PRESENT(PITEND)) THEN
  PQE(KIDIA:KFDIA,1:KLEV)=0._JPRB
ELSE
  PQE(KIDIA:KFDIA,1:KLEV)=PITEND(KIDIA:KFDIA,1:KLEV)
ENDIF


!*        1.3     VERTICAL LOOP EXCLUDING THE BOTTOM LAYER.
                  
DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZDQ=ZRTSPHY*MAX(0.0_JPRB,REPQMI-(PQM(JL,JK)+PTSPHY*PQE(JL,JK)))
    PQE(JL,JK)=PQE(JL,JK)+ZDQ
    ZDQ=ZDQ*ZDP(JL,JK)/ZDP(JL,JK+1)
    PQE(JL,JK+1)=PQE(JL,JK+1)-ZDQ
  ENDDO
ENDDO


!*        1.4      BOTTOM LAYER
!   Note: This last drop adds moisture which can't be compensated
!          from the level bellow

DO JL=KIDIA,KFDIA
  ZDQ=ZRTSPHY*MAX(0.0_JPRB,REPQMI-(PQM(JL,KLEV)+PTSPHY*(PQE(JL,KLEV))))
  PQE(JL,KLEV)=PQE(JL,KLEV)+ZDQ
ENDDO


!*        1.5      EXTRACT THE TENDENCY FROM NEGATIVE Q CORRECTION

IF (PRESENT(PITEND)) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PQE(JL,JK)=PQE(JL,JK)-PITEND(JL,JK)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*         2.       FLUX COMPUTATION

DO JL=KIDIA,KFDIA
  ZFCQNG(JL,0)=0.0_JPRB
ENDDO
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZFCQNG(JL,JK)=-PQE(JL,JK)*ZDP(JL,JK)/RG + ZFCQNG(JL,JK-1)
    PFCQNG(JL,JK)=PFCQNG(JL,JK)+ZFCQNG(JL,JK)
  ENDDO
ENDDO

END ASSOCIATE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('QNEGAT',1,ZHOOK_HANDLE)
END SUBROUTINE QNEGAT
