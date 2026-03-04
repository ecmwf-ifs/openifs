! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE NOSURFTSTP(KDIM, &
 ! Output quantities
 &  FLUX, PDDHS)

!**** *NOSURFTSTP* - Routine called when new surface variable computation is by-passed

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
! PDIAG        : Derived variable for diagnostic quantities
! AUXL         : Derived variables for local auxiliary quantities

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

USE PARKIND1  ,ONLY : JPIM ,   JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, FLUX_TYPE, DDH_SURF_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (DDH_SURF_TYPE)           , INTENT(INOUT) :: PDDHS

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NOSURFTSTP',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*         1.     NECESSARY COMPUTATIONS IF SURFACE TIME STEPPING IS BY-PASSSED

DO JL=KDIM%KIDIA,KDIM%KFDIA
  FLUX%PFTG12(JL)=0.0_JPRB
  FLUX%PFWROD(JL)=0.0_JPRB
  FLUX%PFWRO1(JL)=0.0_JPRB
  FLUX%PFWMLT(JL)=0.0_JPRB
  FLUX%PFWG12(JL)=0.0_JPRB
  FLUX%PFWEV(JL)=0.0_JPRB
ENDDO
! DDH
PDDHS%PDHTSS(KDIM%KIDIA:KDIM%KFDIA,:,1:5)=0.0_JPRB
PDDHS%PDHTSS(KDIM%KIDIA:KDIM%KFDIA,:,11:KDIM%KDHVTSS+KDIM%KDHFTSS)=0.0_JPRB
PDDHS%PDHTTS(KDIM%KIDIA:KDIM%KFDIA,:,1:KDIM%KDHVTTS)=0.0_JPRB
PDDHS%PDHTTS(KDIM%KIDIA:KDIM%KFDIA,:,9:KDIM%KDHVTTS+KDIM%KDHFTTS)=0.0_JPRB
PDDHS%PDHTIS(KDIM%KIDIA:KDIM%KFDIA,:,1:KDIM%KDHVTIS)=0.0_JPRB
PDDHS%PDHTIS(KDIM%KIDIA:KDIM%KFDIA,:,9:KDIM%KDHVTIS+KDIM%KDHFTIS)=0.0_JPRB
PDDHS%PDHSSS(KDIM%KIDIA:KDIM%KFDIA,:,:)=0.0_JPRB
PDDHS%PDHIIS(KDIM%KIDIA:KDIM%KFDIA,:)=0.0_JPRB
PDDHS%PDHWLS(KDIM%KIDIA:KDIM%KFDIA,:,:)=0.0_JPRB

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('NOSURFTSTP',1,ZHOOK_HANDLE)
END SUBROUTINE NOSURFTSTP
