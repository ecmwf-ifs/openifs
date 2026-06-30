MODULE SRFWNG_MOD
CONTAINS
SUBROUTINE SRFWNG(KIDIA,KFDIA,KLEVS,PTMST,KSOTY,&
 & PWL,PWLMX,PWSA,&
 & YDSOIL,&
 & PROS,PROD,PWFSD,&
 & LDLAND,PDHWLS)  
 
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RHOH2O
USE YOS_SOIL  , ONLY : TSOIL

! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFWNG* -  COMPUTES CORRECTIONS TO AVOID NEGATIVE SOIL MOISTURE.

!     PURPOSE.
!     --------
!          THIS ROUTINE COMPUTES CORRECTIONS IN THE SOIL MOISTURE OF
!     THE SKIN, SURFACE AND DEEP RESERVOIR AND RUN-OFF TO AVOID NEGATIVE
!     SOIL MOISTURE VALUES.

!**   INTERFACE.
!     ----------
!          *SRFWNG* IS CALLED FROM *SURF*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLEVS*      NUMBER OF SOIL LAYERS
!    *KSOTY*      SOIL TYPE                                     (1-7)

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP FOR THE PHYSICS                       S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1 (ACCUMULATED,REAL):

!     INPUT/OUTPUT PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PWL*        SKIN RESERVOIR WATER CONTENT                   kg/m**2
!    *PWLMX*      MAXIMUM SKIN RESERVOIR CONTENT                 kg/m**2
!    *PWSA*       SOIL MOISTURE                                M**3/M**3

!     INPUT/OUTPUT PARAMETERS (REAL):
!    *PROS*       RUN-OFF FOR THE SURFACE LAYER                 kg/m**2
!    *PROD*       RUN-OFF FOR DEEPER LAYERS (EXCEPT DRAINAGE)   kg/m**2
!    *PWFSD*      WATER FLUX BETWEEN SURFACE AND DEEP LAYER     kg/m**2/s

!     OUTPUT PARAMETERS (DIAGNOSTIC):
!    *PDHWLS*     Diagnostic array for soil water (see module yomcdh)

!     METHOD.
!     -------
!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION.

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     Original    P.VITERBO      E.C.M.W.F.     16/01/89
!     Modified    P.VITERBO    99-03-26     Tiling of the land surface
!     Modified    J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!     P. Viterbo    24-05-2004      Change surface units
!     G. Balsamo    03-07-2006      Add soil type
!     G. Balsamo    18-08-2015      Rewritten for soil multi-layer 
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS
INTEGER(KIND=JPIM), INTENT(IN)   :: KSOTY(:)

REAL(KIND=JPRB),    INTENT(IN)   :: PTMST
REAL(KIND=JPRB),    INTENT(IN)   :: PWLMX(:)
LOGICAL,            INTENT(IN)   :: LDLAND(:)
TYPE(TSOIL),        INTENT(IN)   :: YDSOIL

REAL(KIND=JPRB),    INTENT(INOUT):: PWL(:)
REAL(KIND=JPRB),    INTENT(INOUT):: PROS(:)
REAL(KIND=JPRB),    INTENT(INOUT):: PWFSD(:)
REAL(KIND=JPRB),    INTENT(INOUT):: PWSA(:,:)
REAL(KIND=JPRB),    INTENT(INOUT):: PDHWLS(:,:,:)

REAL(KIND=JPRB),    INTENT(OUT)  :: PROD(:)

!*         0.2    DECLARATION OF LOCAL VARIABLES.
!                 ----------- -- ----- ----------

INTEGER(KIND=JPIM) :: JL, JK, JS

REAL(KIND=JPRB) :: ZDI(KLEVS), ZIDI(KLEVS), ZRRRI(KLEVS)
REAL(KIND=JPRB) :: ZDWD, ZDZWL,ZWSAT,ZTMST,ZHOH2O
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
!*         1.    SET UP SOME CONSTANTS.
!                --- -- ---- ----------

!*               PHYSICAL CONSTANTS.
!                -------- ----------

IF (LHOOK) CALL DR_HOOK('SRFWNG_MOD:SRFWNG',0,ZHOOK_HANDLE)
ASSOCIATE(LEVGEN=>YDSOIL%LEVGEN, RDAW=>YDSOIL%RDAW, RWSAT=>YDSOIL%RWSAT, &
 & RWSATM=>YDSOIL%RWSATM)

ZTMST=1.0_JPRB/PTMST
ZHOH2O=1.0_JPRB/RHOH2O
DO JK=1,KLEVS
  ZDI(JK)= RDAW(JK)
  ZIDI(JK)= 1.0_JPRB/RDAW(JK)
ENDDO

!     ------------------------------------------------------------------
!*         2. WATER CORRECTIONS AND VEGETATION RATIO UPDATE.
!             ----- ----------- --- ---------- ----- -------
DO JL=KIDIA,KFDIA
  
  IF (LDLAND(JL)) THEN

    IF (LEVGEN) THEN
      JS=KSOTY(JL)
      ZWSAT=RWSATM(JS)
    ELSE
      ZWSAT=RWSAT
    ENDIF
    PROD(JL) = 0._JPRB ! initialize sub-surface runoff to zero 
!          LIMIT PW1,...,PW(n)       0. < PWD < WSAT
    DO JK=1,KLEVS
      IF ( JK == 1_JPIM ) THEN
!          LIMIT PWL        PWL < ZWLMX
        ZDZWL=MAX(0.0_JPRB,PWL(JL)-PWLMX(JL))
        PWL(JL)=PWL(JL)-ZDZWL
        PWSA(JL,JK)=PWSA(JL,JK)+ZDZWL*ZHOH2O*ZIDI(JK)
      ENDIF
!          LIMIT PW(i)       0. < PW(i) < WSAT
      ZDWD=MIN(PWSA(JL,JK),0.0_JPRB)
      PWSA(JL,JK)=PWSA(JL,JK)-ZDWD
      IF (JK < KLEVS) THEN
!          ACCOUNT for water correction on next level
        PWSA(JL,JK+1)=PWSA(JL,JK+1)+ZDWD*ZDI(JK)*ZIDI(JK+1)
      ENDIF
      ZRRRI(JK)=MAX(0.0_JPRB,PWSA(JL,JK)-ZWSAT)
      
      IF ( JK == 1_JPIM ) THEN
!       TOP LAYER RUNOFF AND WATER TRANSFER TO DEEPER LAYER
        PROS(JL)=PROS(JL)+ZDI(JK)*RHOH2O*ZRRRI(JK)
        PWFSD(JL)=PWFSD(JL)+ZDWD*RHOH2O*ZDI(JK)*ZTMST
      ELSE
        ! Sub-surface runoff 
        PROD(JL)=PROD(JL)+ZDI(JK)*RHOH2O*ZRRRI(JK)
      ENDIF  
      
      PWSA(JL,JK)=PWSA(JL,JK)-ZRRRI(JK)
!     DDH diagnostics (runoff)
      IF (SIZE(PDHWLS) > 0) THEN
!      (negative values mean water lost by the layer)
        PDHWLS(JL,JK,6)=PDHWLS(JL,JK,6)-(RHOH2O*ZTMST)*ZDI(JK)*(ZDWD+ZRRRI(JK))
        IF (JK < KLEVS) THEN
          PDHWLS(JL,JK+1,6)=PDHWLS(JL,JK+1,6)+(RHOH2O*ZTMST)*ZDI(JK)*(ZDWD)
        ENDIF
      ENDIF
    ENDDO
  ELSE
!          SEA POINTS
    PROD(JL)=0.0_JPRB
    PWFSD(JL)=0.0_JPRB
  ENDIF
ENDDO
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFWNG_MOD:SRFWNG',1,ZHOOK_HANDLE)

END SUBROUTINE SRFWNG
END MODULE SRFWNG_MOD
