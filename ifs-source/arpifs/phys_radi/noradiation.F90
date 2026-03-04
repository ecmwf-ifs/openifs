! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE NORADIATION(YDSURF,KDIM, &
 ! Output quantities
 & PRAD, FLUX, AUXL, PSURF, SURFL) 

!**** *NORADIATION* - Routine called when radiation scheme is by-passed

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions

!     ==== Input/Output ====
! PRAD     : Derived variable for radiative quantites
! FLUX     : Derived variable for fluxes
! AUXL     : Derived variables for local auxiliary quantities
! PSURF    : Derived variables for general surface quantities
! SURFL    : Derived variables for local surface quantities

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
!      Original : 2012-12-04  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM ,   JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, FLUX_TYPE, AUX_RAD_TYPE, AUX_DIAG_LOCAL_TYPE, &
   & SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_RAD_TYPE)            , INTENT(INOUT) :: PRAD
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL, JK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NORADIATION',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VF=>YDSURF%YSD_VF, YSP_RR=>YDSURF%YSP_RR)

!     ------------------------------------------------------------------

!*         1.     NECESSARY COMPUTATIONS IF RADIATION IS BY-PASSSED

DO JK=1,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PRAD%PNEB(JL,JK)=0.0_JPRB
    PRAD%PQLI(JL,JK)=0.0_JPRB
    PRAD%PQICE(JL,JK)=0.0_JPRB
    PRAD%PHRSW(JL,JK)=0.0_JPRB
    PRAD%PHRLW(JL,JK)=0.0_JPRB
    PRAD%PHRSC(JL,JK)=0.0_JPRB
    PRAD%PHRLC(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JL=KDIM%KIDIA,KDIM%KFDIA
  FLUX%PFRSOC(JL,0)=0.0_JPRB
  FLUX%PFRTHC(JL,0)=0.0_JPRB
  FLUX%PFRSOC(JL,1)=0.0_JPRB
  FLUX%PFRTHC(JL,1)=0.0_JPRB
  FLUX%PFRSOD(JL)=0.0_JPRB
  FLUX%PFRSODC(JL)=0.0_JPRB 
  FLUX%PFRTHD(JL)=0.0_JPRB
  FLUX%PFRTHDC(JL)=0.0_JPRB 
  PRAD%PISUND(JL)=0.0_JPRB
  PRAD%PDSRP (JL)=0.0_JPRB
  FLUX%PUVDF (JL)=0.0_JPRB
  FLUX%PPARF (JL)=0.0_JPRB
  FLUX%PPARCF(JL)=0.0_JPRB
  PRAD%PTINCF(JL)=0.0_JPRB
  PRAD%PFDIR (JL)=0.0_JPRB
  PRAD%PCDIR (JL)=0.0_JPRB
ENDDO
DO JK=0,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    FLUX%PFRSO(JL,JK)=0.0_JPRB
    FLUX%PFRTH(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
! Set surface solar flux to a minimum value for stability in vertical diffusion
DO JL=KDIM%KIDIA,KDIM%KFDIA
  FLUX%PFRSO(JL,KDIM%KLEV)=1.0E-15_JPRB
ENDDO
DO JK=1,KDIM%KTILES
  SURFL%ZALBTI(KDIM%KIDIA:KDIM%KFDIA,JK)=PSURF%PSD_VF(KDIM%KIDIA:KDIM%KFDIA,YSD_VF%YALBF%MP)
ENDDO
PSURF%PEMIS(KDIM%KIDIA:KDIM%KFDIA)=PSURF%PSD_VF(KDIM%KIDIA:KDIM%KFDIA,YSD_VF%YEMISF%MP)
PRAD%PISUND(KDIM%KIDIA:KDIM%KFDIA)=0.0_JPRB
AUXL%ZTSKRAD(KDIM%KIDIA:KDIM%KFDIA)=PSURF%PSP_RR(KDIM%KIDIA:KDIM%KFDIA,YSP_RR%YT%MP9)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('NORADIATION',1,ZHOOK_HANDLE)
END SUBROUTINE NORADIATION
