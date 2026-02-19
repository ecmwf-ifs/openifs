! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE NOCONVECTION(YDEPHY,KDIM, &
 ! Output quantities
 &  FLUX, LLKEYS, PDIAG)

!**** *NOCONVECTION* - Routine called when convection is bypassed

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions

!     ==== Input/Output ====
! FLUX         : Derived variable for fluxes
! LLKEYS       : Derived variable with keys
! PDIAG        : Derived variable for diagnostic quantities

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
!      Original : 2012-11-22  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM ,   JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER, ONLY : DIMENSION_TYPE, AUX_DIAG_TYPE, &
  &                   KEYS_LOCAL_TYPE, FLUX_TYPE
USE YOEPHY   , ONLY : TEPHY

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TEPHY)                     ,INTENT(INOUT) :: YDEPHY
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (KEYS_LOCAL_TYPE)         , INTENT(INOUT) :: LLKEYS
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NOCONVECTION',0,ZHOOK_HANDLE)
ASSOCIATE(LEPCLD=>YDEPHY%LEPCLD)
!     ------------------------------------------------------------------

!*         1.     NECESSARY COMPUTATIONS IF CONVECTION IS BY-PASSSED


DO JK=0,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    FLUX%PDIFCS(JL,JK)=0.0_JPRB
    FLUX%PDIFCQ(JL,JK)=0.0_JPRB
    FLUX%PFHPCL(JL,JK)=0.0_JPRB
    FLUX%PFHPCN(JL,JK)=0.0_JPRB
    FLUX%PFPLCL(JL,JK)=0.0_JPRB
    FLUX%PFPLCN(JL,JK)=0.0_JPRB
    FLUX%PSTRCU(JL,JK)=0.0_JPRB
    FLUX%PSTRCV(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
IF(LEPCLD) THEN
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    LLKEYS%LLCUM(JL)=.FALSE.
    LLKEYS%LLSC (JL)=.FALSE.
    PDIAG%PCAPE(JL)=0.0_JPRB
    PDIAG%ITYPE(JL)=0.0_JPRB
    PDIAG%PVDISCU(JL)=0.0_JPRB
  ENDDO
  DO JK=1,KDIM%KLEV
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PDIAG%PMFD (JL,JK)=0.0_JPRB
      PDIAG%PMFU (JL,JK)=0.0_JPRB
      PDIAG%ZLU  (JL,JK)=0.0_JPRB
      PDIAG%ZLUDE(JL,JK)=0.0_JPRB
      PDIAG%ZLUDELI(JL,JK,1)=0.0_JPRB
      PDIAG%ZLUDELI(JL,JK,2)=0.0_JPRB
      PDIAG%ZLUDELI(JL,JK,3)=0.0_JPRB
      PDIAG%ZLUDELI(JL,JK,4)=0.0_JPRB
      PDIAG%ZSNDE(JL,JK,1)=0.0_JPRB
      PDIAG%ZSNDE(JL,JK,2)=0.0_JPRB
      PDIAG%PRSUD(JL,JK,1)=0.0_JPRB
      PDIAG%PRSUD(JL,JK,2)=0.0_JPRB
    ENDDO
  ENDDO
  DO JK=0,KDIM%KLEV
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FLUX%PFCCQL(JL,JK)=0.0_JPRB
      FLUX%PFCCQN(JL,JK)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('NOCONVECTION',1,ZHOOK_HANDLE)
END SUBROUTINE NOCONVECTION
