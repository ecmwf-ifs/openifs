! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE PHYS_ARRAYS_INI(YDSURF,YDEPHY,YGFL,DIMS,PAUX,ZDDHS,ZDIAG,PSURF,PRAD,FLUX)

!----compiled for Cray with -h contiguous----

!**** *Initialization of derived variables local to physics sub-space

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------

! Derived arrays           Reserved space
! -----------------------------------------------
! PAUX                     PAUX_ARR
! ZDDHS                    PDDHS_ARR
! ZDIAG                    PDIAG_ARR
! PSURF                    PSURF_ARR
!                          KSURF_ARR
! PRAD                     PRAD_ARR
! FLUX                     PFLUX_ARR

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
!      Original : 31-Jan-2013  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     JJMorcrette 20131001 additional volcanic aerosol diagnostics 
!     P. Lopez, ECMWF, July 2015 Added lightning fields in ZDIAG.
!     F. Vana 19-Nov-2016  Negative q pseudo-flux must be initialized to 0 here.
!     E. Dutra/G.Arduini, Jan 2018: snow multi-layer and new snow liquid water content 
!     B. Ingleby,      2019-01-17 Remove PQCFL
!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER          , ONLY : AUX_RAD_TYPE, SURF_AND_MORE_TYPE, AUX_TYPE, &
   &                            AUX_DIAG_TYPE, FLUX_TYPE, DDH_SURF_TYPE, DIMENSION_TYPE
USE YOM_YGFL           , ONLY : TYPE_GFLD
USE YOEPHY             , ONLY : TEPHY
USE YOE_AERODIAG       , ONLY : NPAERODIAG

!     -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)    ,INTENT(INOUT) :: YDSURF
TYPE(TEPHY)    ,INTENT(INOUT) :: YDEPHY
TYPE(TYPE_GFLD),INTENT(INOUT) :: YGFL
TYPE (DIMENSION_TYPE)    ,INTENT(IN)    :: DIMS
TYPE (AUX_TYPE)          ,INTENT(INOUT) :: PAUX
TYPE (DDH_SURF_TYPE)     ,INTENT(INOUT) :: ZDDHS
TYPE (AUX_DIAG_TYPE)     ,INTENT(INOUT) :: ZDIAG
TYPE (SURF_AND_MORE_TYPE),INTENT(INOUT) :: PSURF
TYPE (AUX_RAD_TYPE)      ,INTENT(INOUT) :: PRAD
TYPE (FLUX_TYPE)         ,INTENT(INOUT) :: FLUX
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('PHYS_ARRAYS_INI',0,ZHOOK_HANDLE)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, &
 & LECURR=>YDEPHY%LECURR, &
 & YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSP_OMD=>YDSURF%YSP_OMD)
!-----------------------------------------------------------------------

!   initialize surface DDH fields
ZDDHS%PDHTLS(:,:,:) = 0.0_JPRB
ZDDHS%PDHTSS(:,:,:) = 0.0_JPRB
ZDDHS%PDHTTS(:,:,:) = 0.0_JPRB
ZDDHS%PDHTIS(:,:,:) = 0.0_JPRB
ZDDHS%PDHSSS(:,:,:) = 0.0_JPRB
ZDDHS%PDHIIS(:,:)   = 0.0_JPRB
ZDDHS%PDHWLS(:,:,:) = 0.0_JPRB


!   initialize the diag fields
ZDIAG%PMFUDE_RATE(:,:) = 0.0_JPRB
ZDIAG%PMFDDE_RATE(:,:) = 0.0_JPRB
ZDIAG%PMFU(:,:)        = 0.0_JPRB
ZDIAG%PMFD(:,:)        = 0.0_JPRB
ZDIAG%ZLUDE(:,:)       = 0.0_JPRB
ZDIAG%ZLUDELI(:,:,:)   = 0.0_JPRB
ZDIAG%ZSNDE(:,:,:)     = 0.0_JPRB
ZDIAG%PRSUD(:,:,:)     = 0.0_JPRB
ZDIAG%ZLU(:,:)         = 0.0_JPRB
ZDIAG%ZLRAIN(:,:)      = 0.0_JPRB
ZDIAG%PCOVPTOT(:,:)    = 0.0_JPRB
ZDIAG%ZCAPE(:,:)       = 0.0_JPRB
ZDIAG%PVISIH(:)        = 0.0_JPRB
ZDIAG%PLIGH_TOT(:)     = 0.0_JPRB
ZDIAG%PLIGH_CTG(:)     = 0.0_JPRB
ZDIAG%PCHARGE_LIG(:)   = 0.0_JPRB

ZDIAG%ITYPE(:)         =0
ZDIAG%IPBLTYPE(:)      =0
ZDIAG%ICBOT(:)         =0
ZDIAG%ICTOP(:)         =0
ZDIAG%ICBOT_LIG(:)     =0
ZDIAG%ICTOP_LIG(:)     =0

!*        SET SURFACE TENDENCIES TO ZERO
PSURF%PTSAE1(:,:)=0.0_JPRB
PSURF%PWSAE1(:,:)=0.0_JPRB
PSURF%PTIAE1(:,:)=0.0_JPRB
PSURF%PTLE1(:)=0.0_JPRB
PSURF%PWLE1(:)=0.0_JPRB
PSURF%PSNSE1(:,:)=0.0_JPRB       ! Tendency of snow mass
PSURF%PASNE1(:,:)=0.0_JPRB       ! Tendency of snow albedo (only 1st layer used)
PSURF%PRSNE1(:,:)=0.0_JPRB       ! Tendency of snow density
PSURF%PTSNE1(:,:)=0.0_JPRB       ! Tendency of snow temperature
PSURF%PWSNE1(:,:)=0.0_JPRB       ! Tendency of snow liquid water content
PSURF%PTLICEE1(:)=0.0_JPRB       ! tendency of lake ice temperature
PSURF%PTLMNWE1(:)=0.0_JPRB       ! tendency of lake totat layer temperature
PSURF%PTLWMLE1(:)=0.0_JPRB       ! tendency of lake mixed layer temperature
PSURF%PTLBOTE1(:)=0.0_JPRB       ! tendency of lake bottom layer temperature
PSURF%PTLSFE1(:)=0.0_JPRB        ! tendency of lake shape factor - 
PSURF%PHLICEE1(:)=0.0_JPRB       ! tendency of lake ice depth m
PSURF%PHLMLE1(:)=0.0_JPRB        ! tendency of lake mixed layer depth m/s

IF (.NOT. LECURR) THEN
  PSURF%PSD_VF(:,YSD_VF%YUCUR%MP)=0.0_JPRB
  PSURF%PSD_VF(:,YSD_VF%YVCUR%MP)=0.0_JPRB
ENDIF

!- set all aerosol optical depths to zero (except accumulated)
PSURF%PSD_VD(:,YSD_VD%YODSS%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODDU%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODOM%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODBC%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODSU%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODNI%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODAM%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODVFA%MP) =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YODVSU%MP) =0._JPRB
!- set all PMs to zero
PSURF%PSD_VD(:,YSD_VD%YAEPM1%MP) =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YAEPM25%MP)=0._JPRB
PSURF%PSD_VD(:,YSD_VD%YAEPM10%MP)=0._JPRB
!- set all UVs to zero
PSURF%PSD_VD(:,YSD_VD%YUVBED%MP)  =0._JPRB
PSURF%PSD_VD(:,YSD_VD%YUVBEDCS%MP)=0._JPRB


!*      SET FLUXES TO ZERO
FLUX%PDIFTI(:,:) = 0.0_JPRB
FLUX%PDIFTL(:,:) = 0.0_JPRB
FLUX%PFCQNNG(:,:) = 0.0_JPRB
FLUX%PFCQLNG(:,:) = 0.0_JPRB
FLUX%PFCSQN(:,:) = 0.0_JPRB
FLUX%PFCSQL(:,:) = 0.0_JPRB
FLUX%PFSQRF(:,:) = 0.0_JPRB
FLUX%PFSQSF(:,:) = 0.0_JPRB
FLUX%PFCQRNG(:,:) = 0.0_JPRB
FLUX%PFCQSNG(:,:) = 0.0_JPRB
FLUX%PFSQLTUR(:,:) = 0.0_JPRB
FLUX%PFSQITUR(:,:) = 0.0_JPRB
FLUX%PFCCQL(:,:) = 0.0_JPRB
FLUX%PFCCQN(:,:) = 0.0_JPRB
FLUX%PFPLCL(:,:) = 0.0_JPRB
FLUX%PFPLCN(:,:) = 0.0_JPRB
FLUX%PFPLSL(:,:) = 0.0_JPRB
FLUX%PFPLSN(:,:) = 0.0_JPRB
FLUX%PFCQNG(:,:) = 0.0_JPRB
!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PHYS_ARRAYS_INI',1,ZHOOK_HANDLE)
END SUBROUTINE PHYS_ARRAYS_INI



