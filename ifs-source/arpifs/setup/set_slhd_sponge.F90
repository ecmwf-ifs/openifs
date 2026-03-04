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

SUBROUTINE SET_SLHD_SPONGE(YDGEOMETRY,YDDYN)

!**** *SET_SLHD_SPONGE*  - Prescribe the activity for SLHD

!     Purpose.
!     --------
!        Initialize masking functions for areas where SLHD
!        is not needed. 

!        ECMWF uses SLHD as grid-point storm killer for high 
!        atmosphere. This is esential to get rid of spurious
!        noise resulting in SL departure point perturbation 
!        in TL/AD code to be 40-70 stencils away from the 
!        original interpolation stencil.

!**   Interface.
!     ----------
!        *CALL* *SET_SLHD_SPONGE

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Filip Vana  *ECMWF*
!      Original : 14-Jun-2018

!     Modifications.
!     --------------
!     F. Vana 21-May-2019  Setup for NLEV_SPONGE
!     F. Vana 11-Jul-2019: Better control for the defaults

!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : LOUTPUT, NPRINTLEV
USE YOMDYN   , ONLY : TDYN
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN

INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SET_SLHD_SPONGE',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, STPRE=>YDGEOMETRY%YRSTA%STPRE, &
 &  LSLHDHEAT=>YDDYN%LSLHDHEAT, NLEV_SPONGE=>YDDYN%NLEV_SPONGE, &
 &  SLHD_P_LOW=>YDDYN%SLHD_P_LOW, SLHD_P_HIGH=>YDDYN%SLHD_P_HIGH, &
 &  SLHD_MASK_U=>YDDYN%SLHD_MASK_U, SLHD_MASK_T=>YDDYN%SLHD_MASK_T)

!     ------------------------------------------------------------------

!*    1.  Set masking functions
!          ------------------

DO JLEV=1,NFLEVG
  SLHD_MASK_U(JLEV)=MIN(1._JPRB,MAX(0._JPRB,(SLHD_P_LOW-STPRE(JLEV))/(SLHD_P_LOW-SLHD_P_HIGH)))
  IF (SLHD_MASK_U(JLEV) == 0._JPRB) NLEV_SPONGE=MIN(NLEV_SPONGE,JLEV-1)
ENDDO

IF (LSLHDHEAT) SLHD_MASK_T(1:NFLEVG)=SLHD_MASK_U(1:NFLEVG)

IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  WRITE(UNIT=NULOUT,FMT='('' SLHD MASKING FUNCTION '')')
  DO JLEV=1,NFLEVG
    WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,F20.10)')&
     & JLEV,SLHD_MASK_U(JLEV)
  ENDDO
  WRITE(UNIT=NULOUT,FMT='('' NLEV_SPONGE =  '',I6)') NLEV_SPONGE
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SET_SLHD_SPONGE',1,ZHOOK_HANDLE)
END SUBROUTINE SET_SLHD_SPONGE
