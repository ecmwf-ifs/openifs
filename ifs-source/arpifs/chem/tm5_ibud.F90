! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_IBUD(YGFL,KIDIA,KFDIA,KLON,PRR,PRJ,PY,PCR2,PCR3)

!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!--------------------------------------------------------------------------
!
!*** CALCULATION OF REACTION BUDGETS due to gas-phase chemistry and photolysis
!
!--------------------------------------------------------------------------
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5_incbud* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! PRR   :  gas phase reaction rates at current level
! PRJ   :  photolysis rates at current level
! PY(KLON,NCHEM)        : concentrations after reaction  (molec/cm3?)
!
!
! OUTPUTS:
! -------
! PCR2 (KLON,NPHOTO)    : budget contribution due to photolysis           
! PCR3 (KLON,NCHEM)     : budget contribution due to chem          
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!        TM5-community    
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2010-10-08



USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE TM5_CHEM_MODULE , ONLY : NRJ, NRR, NREAC
USE TM5_PHOTOLYSIS , ONLY : NPHOTO

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON
REAL(KIND=JPRB),INTENT(IN)    :: PY(KLON,YGFL%NCHEM+3)  ! concentrations
REAL(KIND=JPRB),INTENT(IN)    :: PRJ(KLON,NPHOTO)  ! photolysis reaction rates at specific level
REAL(KIND=JPRB),INTENT(IN)    :: PRR(KLON,NREAC)   ! chemistry reaction rates at specific level
REAL(KIND=JPRB),INTENT(INOUT) :: PCR2(KLON,NPHOTO) ! photolysis budget accumulators
REAL(KIND=JPRB),INTENT(INOUT) :: PCR3(KLON,NREAC)  ! chemistry  budget accumulators
! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JR, J1, J2


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_IBUD',0,ZHOOK_HANDLE )

 DO JR=1_JPIM,NPHOTO !photolysis rates
    J1=NRJ(JR)
    DO JL=KIDIA,KFDIA
      IF(J1 > 0_JPIM) PCR2(JL,JR)=PCR2(JL,JR)+PRJ(JL,JR)*PY(JL,J1)
    ENDDO
 ENDDO!JR=1,nj


 DO JR=1_JPIM,NREAC !reactions
    J1=NRR(JR,1) !make sure J1 > 0
    J2=NRR(JR,2)
    IF (J2 > 0_JPIM) THEN
      DO JL=KIDIA,KFDIA
        PCR3(JL,JR)= PCR3(JL,JR)+PY(JL,J1)*PY(JL,J2)*PRR(JL,JR)
      ENDDO
    ELSE
      DO JL=KIDIA,KFDIA
        PCR3(JL,JR)= PCR3(JL,JR)+PY(JL,J1)*PRR(JL,JR)
      ENDDO
    ENDIF
 ENDDO  !JR=1,nreac



IF (LHOOK) CALL DR_HOOK('TM5_IBUD',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_IBUD
