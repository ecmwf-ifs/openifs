! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_GS_LIQ(KSTEP, KMONTH, KIDIA, KFDIA, KLON, KLEV, KTROPOP, PRSF1, PLAT, PTP,PAER,PSA_SIZEDIST,PAER_INFO)


!**   DESCRIPTION 
! ----------------------------------------------------------------------
USE PARKIND1  ,    ONLY : JPIM,   JPRB
USE YOMHOOK   ,    ONLY : LHOOK,  DR_HOOK, JPHOOK
USE BASCOE_MODULE, ONLY :  NBINS, NAER, IAER_NTOT, IAER_SAD, &
! J. Debosscher: related to SAD climatology:
  & NLAT_CLIM, NLEV_CLIM, SAD_CLIM,P_CLIM,LAT_CLIM
!-----------------------------------------------------------------------
IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
INTEGER(KIND=JPIM),INTENT(IN)      :: KSTEP, KMONTH, KIDIA, KFDIA, KLON, KLEV
INTEGER(KIND=JPIM),INTENT(IN)      :: KTROPOP(KLON)
REAL   (KIND=JPRB),INTENT(IN)      :: PTP(KLON,KLEV),PRSF1(KLON,KLEV)
REAL   (KIND=JPRB),INTENT(IN)      :: PLAT(KLON)
REAL   (KIND=JPRB),INTENT(OUT)     :: PAER(KLON,KLEV,NAER)
REAL   (KIND=JPRB),INTENT(OUT)     :: PSA_SIZEDIST(KLON,KLEV,NBINS)
REAL   (KIND=JPRB),INTENT(INOUT)   :: PAER_INFO(KLON,KLEV)


!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM)        :: JTROPOP(KLON)
REAL(KIND=JPRB),PARAMETER :: ZXR = 0.07E-6, ZXMED = 1.76  ! def. of aerosol distribution
REAL(KIND=JPRB)           :: ZXN, ZTEMP, ZND(NBINS)
! LOGICAL                   :: LL_PRESENCE_STS

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JL,JK,JN
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM)         :: JLAT,JLEV
!-----------------------------------------------------------------------
!-------------------------------------------------------------------
#include "bascoe_initsp.intfb.h"
#include "bascoe_sage_number_density.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_GS_LIQ',0,ZHOOK_HANDLE )

! Only read in climatological field once
IF (KSTEP == 0_JPIM) THEN
  DO JL=KIDIA,KFDIA
    DO JK = 1,KLEV
      CALL BASCOE_SAGE_NUMBER_DENSITY( PRSF1(JL,JK), ABS(PLAT(JL)), ZXR, ZXMED, ZXN )

      ! Scale to something more realistic (Post-Pinatubo)
      ZXN = ZXN * 0.2
      
      ! fill field for later use
      PAER_INFO(JL,JK) = ZXN
    ENDDO
  ENDDO
  ! initialize down to surface level onwards
  JTROPOP(KIDIA:KFDIA) = KLEV
ELSE
  ! consequtive time steps: only constrain down to tropopause level 
  JTROPOP(KIDIA:KFDIA) = KTROPOP(KIDIA:KFDIA)+6
ENDIF

! Initialize...
PAER(KIDIA:KFDIA,1:KLEV,1:NAER) = 0._JPRB
PSA_SIZEDIST(KIDIA:KFDIA,1:KLEV,1:NBINS) = 0._JPRB

DO JL=KIDIA,KFDIA
  ! Only fix in stratosphere...
  DO JK = 1,JTROPOP(JL)
    !VH CALL BASCOE_SAGE_NUMBER_DENSITY( PRSF1(JL,JK), ABS(PLAT(JL)), ZXR, ZXMED, ZXN )

    ! Scale to something more realistic (Post-Pinatubo)
    !VH ZXN = ZXN * 0.2
    ZXN = PAER_INFO(JL,JK)

    ZTEMP = MAX( PTP(JL,JK), 200._JPRB )
    CALL BASCOE_INITSP( ZXN*1.0E-6, ZXR, ZXMED, ZTEMP, PRSF1(JL,JK), ZND )
!VH Try simple climatology
!VH    CALL BASCOE_INITSP( 0._JPRB, ZXR, ZXMED, ZTEMP, PRSF1(JL,JK), ZND )
    DO JN = 1, NBINS
      PSA_SIZEDIST(JL,JK,JN) = ZND(JN)  ! particles / (kg of air)
    ENDDO
    PAER(JL,JK,IAER_NTOT) = SUM(ZND)
     

! ----------------------------------------------------------------------
!  Get basic air parameters
! ----------------------------------------------------------------------

!    LL_PRESENCE_STS = ANY( ZND(:) > 1.e-2_JPRB )

! ----------------------------------------------------------------------
!  Where particles ares present, calculate surface area densities gs_sts
!  conv = 1e-2*rho_air (density, kg/m3), PRSF1 in Pa.
!   1e-2 for m2/m3 -> cm2/cm3. For volume density, use 1e2*conv
! ----------------------------------------------------------------------
!    IF( LL_PRESENCE_STS ) THEN
!      ZCONV = 1.E-2_JPRB * PRSF1(JL,JK) * 28.96_JPRB / ( PTP(JL,JK)*8.3E3_JPRB )
!      !VH put code in line ... CALL SURFCE( NBINS, ZND(1), PTSIZE, ZX3 )
!      !     This subroutine calculates the total surface area density of an
!      !     ensemble of particles of a given type, i.e. particle surface area
!      !     per kg of air.
!      ZX3=0.0_JPRB
!      DO JB=1,NBINS
!        ZX3=ZX3+PTSIZE(JB,2)*ZND(JB)
!      ENDDO
!
!      PAER(JL,JK,IAER_SAD) = ZX3*ZCONV
!
!    ELSE
!      PAER(JL,JK,IAER_SAD) = 0._JPRB
!    ENDIF
!---------------------------------------------------------
! J. Debosscher: use SAD climatology from file in stead
    JLAT=MINLOC(ABS(LAT_CLIM(1:NLAT_CLIM)-PLAT(JL)),1)
    JLEV=MINLOC(ABS(P_CLIM(KMONTH,JLAT,1:NLEV_CLIM)-PRSF1(JL,JK)),1)

    !write(NULOUT,*)'JLAT',JLAT
    !write(NULOUT,*)'KMONTH',KMONTH


    PAER(JL,JK,IAER_SAD)=1.E-2_JPRB*SAD_CLIM(KMONTH,JLAT,JLEV)
    !write(NULOUT,*)'SAD CLIM OK'

!---------------------------------------------------------
  ENDDO
ENDDO



IF (LHOOK) CALL DR_HOOK('BASCOE_GS_LIQ',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_GS_LIQ

