! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE NOTURBULENCE(YDSURF,KDIM, &
 ! Output quantities
  & PDIAG,FLUX,PSURF,SURFL,AUXL,PDDHS,TENDENCY_VDF)

!**** *NOTURBULENCE* - Routine called when turbulence scheme is bypassed

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions

!     ==== Input/Output ====
! PDIAG    : Derived variable for diagnostic quantities
! FLUX     : Derived variable for fluxes
! PSURF    : Derived variable for general surface quantities
! SURFL    : Derived variable for local surface quantities
! AUXL     : Derived variables for local auxiliary quantities
! PDDHS    : Derived variable for surface DDH quantities

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
! B. Ingleby      2019-01-17  Change PQCFL to Y2M

!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM ,   JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, &
   & FLUX_TYPE,SURF_AND_MORE_TYPE,SURF_AND_MORE_LOCAL_TYPE,DDH_SURF_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (DDH_SURF_TYPE)           , INTENT(INOUT) :: PDDHS
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_VDF
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JK, JL, JTILE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NOTURBULENCE',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VD=>YDSURF%YSD_VD, YSP_RR=>YDSURF%YSP_RR)

!     ------------------------------------------------------------------

!*         1.     NECESSARY COMPUTATIONS IF TURBULENT SCJEME IS BY-PASSSED

DO JL=KDIM%KIDIA,KDIM%KFDIA
  PDIAG%PVDIS(JL)=0.0_JPRB
  PDIAG%PUSTRG(JL)=0.0_JPRB
  PDIAG%PVSTRG(JL)=0.0_JPRB
  PDIAG%PVDISG(JL)=0.0_JPRB
  FLUX%PFTLHEV(JL)=0.0_JPRB
  FLUX%PFTLHSB(JL)=0.0_JPRB
  PSURF%PSD_VD(JL,YSD_VD%Y10U%MP)=0.0_JPRB
  PSURF%PSD_VD(JL,YSD_VD%Y10V%MP)=0.0_JPRB
  PSURF%PSD_VD(JL,YSD_VD%Y2T%MP) =0.0_JPRB
  PSURF%PSD_VD(JL,YSD_VD%Y2D%MP) =0.0_JPRB
  PSURF%PSD_VD(JL,YSD_VD%Y2SH%MP)=0.0_JPRB
  PSURF%PEVAPMU(JL) = 0.0_JPRB
  PSURF%PSD_VD(JL,YSD_VD%YBLH%MP)=0.0_JPRB
  FLUX%PFWSB(JL)=0.0_JPRB
  PDIAG%PI10FG(JL)=0.0_JPRB
  PDIAG%PUSTRG(JL)=0.0_JPRB
  PDIAG%PVSTRG(JL)=0.0_JPRB
  PDIAG%PVDISG(JL)=0.0_JPRB
ENDDO

SURFL%ZEVAPSNW(KDIM%KIDIA:KDIM%KFDIA)=0.0_JPRB
! DDH
PDDHS%PDHTLS(KDIM%KIDIA:KDIM%KFDIA,:,:)=0.0_JPRB
PDDHS%PDHTSS(KDIM%KIDIA:KDIM%KFDIA,:,6)=0.0_JPRB
PDDHS%PDHTSS(KDIM%KIDIA:KDIM%KFDIA,:,KDIM%KDHVTSS+1:KDIM%KDHVTSS+4)=0.0_JPRB
PDDHS%PDHTTS(KDIM%KIDIA:KDIM%KFDIA,:,KDIM%KDHVTTS+1:KDIM%KDHVTTS+4)=0.0_JPRB
PDDHS%PDHTIS(KDIM%KIDIA:KDIM%KFDIA,:,KDIM%KDHVTIS+1:KDIM%KDHVTIS+4)=0.0_JPRB
DO JK=0,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    FLUX%PDIFTQ(JL,JK)=0.0_JPRB
    FLUX%PDIFTS(JL,JK)=0.0_JPRB
    FLUX%PDIFTL(JL,JK)=0.0_JPRB
    FLUX%PDIFTI(JL,JK)=0.0_JPRB
    FLUX%PSTRTU(JL,JK)=0.0_JPRB
    FLUX%PSTRTV(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JTILE=1,KDIM%KTILES
  SURFL%ZFRSOTI(KDIM%KIDIA:KDIM%KFDIA,JTILE)=FLUX%PFRSO(KDIM%KIDIA:KDIM%KFDIA,KDIM%KLEV)
  PSURF%PUSTRTI(KDIM%KIDIA:KDIM%KFDIA,JTILE)=0.0_JPRB
  PSURF%PVSTRTI(KDIM%KIDIA:KDIM%KFDIA,JTILE)=0.0_JPRB
  PSURF%PAHFSTI(KDIM%KIDIA:KDIM%KFDIA,JTILE)=0.0_JPRB
  PSURF%PEVAPTI(KDIM%KIDIA:KDIM%KFDIA,JTILE)=0.0_JPRB
  PSURF%PTSKTI(KDIM%KIDIA:KDIM%KFDIA,JTILE)=PSURF%PSP_RR(KDIM%KIDIA:KDIM%KFDIA,YSP_RR%YT%MP9)
ENDDO
TENDENCY_VDF%U(:,:)=0.0_JPRB
TENDENCY_VDF%V(:,:)=0.0_JPRB
TENDENCY_VDF%T(:,:)=0.0_JPRB
TENDENCY_VDF%Q(:,:)=0.0_JPRB
TENDENCY_VDF%A(:,:)=0.0_JPRB
TENDENCY_VDF%O3(:,:)=0.0_JPRB
TENDENCY_VDF%CLD(:,:,:)=0.0_JPRB


!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('NOTURBULENCE',1,ZHOOK_HANDLE)
END SUBROUTINE NOTURBULENCE
