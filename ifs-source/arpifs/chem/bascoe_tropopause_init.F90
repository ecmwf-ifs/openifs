! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_TROPOPAUSE_INIT
!**   DESCRIPTION 
!     ----------
!
!   Part of BASCOE / TM5 routines for IFS chemistry: 
!     AUTHOR.
!     -------
!        Coded in C-IFS by VINCENT HUIJNEN    *KNMI*
!        Original code from BASCOE_CTM v4s09, simonc@oma.be, June 2008
!
!-----------------------------------------------------------------------
! Calculate jtropop, the pressure index of tropopause
!                                                simonc, v2s35, Jul 2003
! VH - Here setup the fall-back option
!-----------------------------------------------------------------------
USE BASCOE_MODULE      , ONLY : PTROP_SOC_B,NLAT_PTROPO, LATS_PTROPO
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK


IMPLICIT NONE
!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------

! SOCRATES tropopause vertical index in SOCRATES grid
INTEGER(KIND=JPIM),PARAMETER   :: JNLAT_SOC=35
INTEGER(KIND=JPIM),DIMENSION(35),PARAMETER :: IZM=(/ &
!VH     &    2*9, 2*10, 2*11, 12, 13, 14, 15, 2*16, 11*17, 2*16, &
!VH     &           15, 14, 13, 12, 2*11, 2*10, 2*9 /)
     &    9,9, 10,10, 11,11, 12, 13, 14, 15, 16,16, &
     &   17,17,17,17,17,17,17,17,17,17,17, 16,16, &
     &           15, 14, 13, 12, 11,11, 10,10, 9,9 /)
REAL(KIND=JPRB), DIMENSION(35) :: ZLAT_SOC, ZPTROP_SOC
REAL(KIND=JPRB), PARAMETER     :: ZPSURF_STD = 101325. ! std p at surf (Pa)

INTEGER(KIND=JPIM)             :: JLAT
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
LOGICAL                        :: LL_OK
!----------------------------------------------------------------------------
#include "abor1.intfb.h"
#include "bascoe_interp8.intfb.h"
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_TROPOPAUSE_INIT',0,ZHOOK_HANDLE )


 ! interpolate to 1-deg resolution
 LATS_PTROPO(1:NLAT_PTROPO) = (/ ( ( 180._JPRB/NLAT_PTROPO ) * REAL(JLAT) - 90._JPRB, JLAT=1, NLAT_PTROPO) /)

 ! SOCRATES Climatology
 ZLAT_SOC(1:JNLAT_SOC)  = (/ ( 5. * REAL(JLAT) - 90., JLAT=1, JNLAT_SOC) /)
 ZPTROP_SOC(1:JNLAT_SOC) = ZPSURF_STD * EXP( - REAL(IZM(1:JNLAT_SOC)-1) / 7. )
 CALL BASCOE_INTERP8( NLAT_PTROPO, LATS_PTROPO, PTROP_SOC_B, JNLAT_SOC, ZLAT_SOC, ZPTROP_SOC, &
     &                LL_OK, 'asbndry')
 IF( .NOT. LL_OK ) THEN
    CALL ABOR1("BASCOE_TROPOPAUSE_INIT: INTERP of ZPTROP_SOC failed")
 ENDIF

IF (LHOOK) CALL DR_HOOK('BASCOE_TROPOPAUSE_INIT',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_TROPOPAUSE_INIT
