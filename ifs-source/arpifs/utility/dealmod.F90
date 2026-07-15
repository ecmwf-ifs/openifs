! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE DEALMOD(YDRADF,YDSLPHY,YDSPPT,YDSPPT_CONFIG)

!**** *DEALMOD* - Deallocate unnecessary model fields for screening

!     Purpose.
!     -------
!           Deallocate unnecessary model fields for screening

!**   Interface.
!**   ---------
!**         DEALMOD is called from routine SCREEN

!     Author.
!     -------
!        H. JARVINEN *ECMWF*

!     Modifications.
!     --------------
!        Original : 96-12-10
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!        JJMorcrette 20070321 Prognostic aerosols for rad and clouds
!        JJMorcrette 20091201 Total and clear-sky direct SW radiation flux at surface
!        M.Leutbecher  04-Nov-2009  Fields for stochastic tendency pertns
!        R. El Khatib  20-Aug-2012  GAUXBUF removed and replaced by HFPBUF
!        M Ahlgrimm    31 Oct 2011  Surface Downward clear-sky LW and SW fluxes
!        M. Fisher     7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!        R. Hogan      June 2014    Added DerivativeLw
!        SJ Lock       Jan-2016     Tidying SPPT fields
!------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMSLPHY , ONLY : TSLPHY
USE YOMRADF  , ONLY : TRADF
USE YOMSPSDT , ONLY : TSPPT_CONFIG, TSPPT_DATA
USE SPECTRAL_ARP_MOD,     ONLY : DEALLOCATE_ARP
USE GRIDPOINT_FIELDS_MIX, ONLY : DEALLOCATE_GRID
USE DEALLOCATE_IF_ASSOCIATED_MOD, ONLY : DEALLOCATE_IF_ASSOCIATED
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=), DEALLOCATE_SPEC

!------------------------------------------------------

IMPLICIT NONE

TYPE(TRADF),        INTENT(INOUT) :: YDRADF
TYPE(TSLPHY),       INTENT(INOUT) :: YDSLPHY
TYPE(TSPPT_DATA),   INTENT(INOUT) :: YDSPPT
TYPE(TSPPT_CONFIG), INTENT(INOUT) :: YDSPPT_CONFIG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
TYPE(SPECTRAL_FIELD) :: YDSPEC
!------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DEALMOD',0,ZHOOK_HANDLE)
!------------------------------------------------------

!*        1.    DEALLOCATE
!               ----------

!*           1.7   DECLARED IN YOMSP

CALL DEALLOCATE_SPEC(YDSPEC)

IF(ALLOCATED (YDSLPHY%SAVTEND)  )  DEALLOCATE (YDSLPHY%SAVTEND)

IF (ALLOCATED(YDRADF%EMTD      )) DEALLOCATE        (YDRADF%EMTD      )
IF (ALLOCATED(YDRADF%TRSW      )) DEALLOCATE        (YDRADF%TRSW      )
IF (ALLOCATED(YDRADF%EMTC      )) DEALLOCATE        (YDRADF%EMTC      )
IF (ALLOCATED(YDRADF%TRSC      )) DEALLOCATE        (YDRADF%TRSC      )
IF (ALLOCATED(YDRADF%TAUAER    )) DEALLOCATE        (YDRADF%TAUAER    )
IF (ALLOCATED(YDRADF%SRSWD     )) DEALLOCATE        (YDRADF%SRSWD     )
IF (ALLOCATED(YDRADF%SRSWDC    )) DEALLOCATE        (YDRADF%SRSWDC    )
IF (ALLOCATED(YDRADF%SRLWD     )) DEALLOCATE        (YDRADF%SRLWD     )
IF (ALLOCATED(YDRADF%SRLWDC    )) DEALLOCATE        (YDRADF%SRLWDC    )
IF (ALLOCATED(YDRADF%SRSWDCS   )) DEALLOCATE        (YDRADF%SRSWDCS   )
IF (ALLOCATED(YDRADF%SRLWDCS   )) DEALLOCATE        (YDRADF%SRLWDCS   )
IF (ALLOCATED(YDRADF%SRSWDV    )) DEALLOCATE        (YDRADF%SRSWDV    )
IF (ALLOCATED(YDRADF%SRSWDUV   )) DEALLOCATE        (YDRADF%SRSWDUV   )
IF (ALLOCATED(YDRADF%EDRO      )) DEALLOCATE        (YDRADF%EDRO      )
IF (ALLOCATED(YDRADF%SRSWPAR   )) DEALLOCATE        (YDRADF%SRSWPAR   )
IF (ALLOCATED(YDRADF%SRSWUVB   )) DEALLOCATE        (YDRADF%SRSWUVB   )
IF (ALLOCATED(YDRADF%SRSWPARC  )) DEALLOCATE        (YDRADF%SRSWPARC  )
IF (ALLOCATED(YDRADF%SRSWTINC  )) DEALLOCATE        (YDRADF%SRSWTINC  )
IF (ALLOCATED(YDRADF%SRFDIR    )) DEALLOCATE        (YDRADF%SRFDIR    )
IF (ALLOCATED(YDRADF%SRCDIR    )) DEALLOCATE        (YDRADF%SRCDIR    )
CALL DEALLOCATE_IF_ASSOCIATED(YDRADF%DERIVATIVELW)

!*           1.8   DECLARED IN YOMSPSDT

IF (YDSPPT_CONFIG%LSPSDT) THEN
  YDSPPT%YSPSDT_AR1            => NULL()
  YDSPPT%YGPSDT                => NULL()
  YDSPPT%NSEED_SDT             => NULL()
  YDSPPT_CONFIG%CIPATINIT_SDT  => NULL()
  YDSPPT_CONFIG%COPATTRUN_SDT  => NULL()
  YDSPPT_CONFIG%COPATSP_SDT    => NULL()
  YDSPPT_CONFIG%COPATGP_SDT    => NULL()
ENDIF
!------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DEALMOD',1,ZHOOK_HANDLE)
END SUBROUTINE DEALMOD
