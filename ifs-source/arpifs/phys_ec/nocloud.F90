! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE NOCLOUD(YDECLDP,KDIM, &
 ! Output quantities
 &  FLUX, PDIAG, TENDENCY_LOC, TENDENCY_SATADJ)

!**** *NOCLOUD* - Routine called when cloud schemes are by-passed

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! KEYMASK  : Mask of GFL fields usage

!     ==== Input/Output ====
! FLUX     : Derived variable for fluxes
! PDIAG    : Derived variable for diagnostics

!     ==== Output ====
! TENDENCY_LOC    :  local tendencies from cloudsc
! TENDENCY_SATADJ :  local tendencies from cloud_satadj

!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-28  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      F. Vana  08-Apr-2016

!-----------------------------------------------------------------------

USE YOECLDP  , ONLY : TECLDP, NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE, FLUX_TYPE, STATE_TYPE, MASK_GFL_TYPE, &
 &                    AUX_DIAG_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TECLDP)                   , INTENT(INOUT) :: YDECLDP
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_SATADJ

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NOCLOUD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*         1.     NECESSARY COMPUTATIONS IF CLOUD SCHEMES ARE BY-PASSSED

DO JK=0,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    FLUX%PFHPSL(JL,JK)=0.0_JPRB
    FLUX%PFHPSN(JL,JK)=0.0_JPRB
    FLUX%PFPLSL(JL,JK)=0.0_JPRB
    FLUX%PFPLSN(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

PDIAG%PPRECTYPE(KDIM%KIDIA:KDIM%KFDIA) = 0._JPRB
PDIAG%PFZRA    (KDIM%KIDIA:KDIM%KFDIA) = 0._JPRB

! zeroing process tendencies
DO JK=1,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    TENDENCY_LOC%T(JL,JK)=0.0_JPRB
    TENDENCY_LOC%Q(JL,JK)=0.0_JPRB
    IF (TENDENCY_LOC%KEYMASK%A ) TENDENCY_LOC%A(JL,JK)         =0.0_JPRB
    IF (TENDENCY_LOC%KEYMASK%QL) TENDENCY_LOC%CLD(JL,JK,NCLDQL)=0.0_JPRB
    IF (TENDENCY_LOC%KEYMASK%QI) TENDENCY_LOC%CLD(JL,JK,NCLDQI)=0.0_JPRB
    IF (TENDENCY_LOC%KEYMASK%QR) TENDENCY_LOC%CLD(JL,JK,NCLDQR)=0.0_JPRB
    IF (TENDENCY_LOC%KEYMASK%QS) TENDENCY_LOC%CLD(JL,JK,NCLDQS)=0.0_JPRB
    TENDENCY_SATADJ%T(JL,JK)=0.0_JPRB
    TENDENCY_SATADJ%Q(JL,JK)=0.0_JPRB
    IF (TENDENCY_SATADJ%KEYMASK%A ) TENDENCY_SATADJ%A(JL,JK)         =0.0_JPRB
    IF (TENDENCY_SATADJ%KEYMASK%QL) TENDENCY_SATADJ%CLD(JL,JK,NCLDQL)=0.0_JPRB
    IF (TENDENCY_SATADJ%KEYMASK%QI) TENDENCY_SATADJ%CLD(JL,JK,NCLDQI)=0.0_JPRB
    IF (TENDENCY_SATADJ%KEYMASK%QR) TENDENCY_SATADJ%CLD(JL,JK,NCLDQR)=0.0_JPRB
    IF (TENDENCY_SATADJ%KEYMASK%QS) TENDENCY_SATADJ%CLD(JL,JK,NCLDQS)=0.0_JPRB
  ENDDO
ENDDO


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('NOCLOUD',1,ZHOOK_HANDLE)
END SUBROUTINE NOCLOUD
