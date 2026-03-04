! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_J_CALC( KLEV, PZKM, PSZA, PTP, PDENSA , PALB, PFBEAMR, PAJVAL)

!**   DESCRIPTION
!     ----------
!
!   Computes photodissociation rates for the stratospheric chemistry
!       on the set of KLEV vertical model levels, using an adapted
!       version of TUV
!   The vertical domain for TUV is made of the model levels for which
!       J rates are calculated, and extra levels on top from the
!       absorption model defined in BASCOE_J_MODULE.
!   To speed-up processing, TUV is computed on a subset of these extra
!       levels, specified by 'INJLEV'
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_J_CALC* IS CALLED FROM *CHEM_BASCOE[TM5]*.
!
! INPUTS:
! -------
! KLEV          : lowest model level where to compute J rates
! PZKM(KLEV)    : model height (km)
! PSZA          : Solar zenith angle
! PTP(KLEV)     : model temperature
! PDENSA(KLEV,NABSPEC)  : densities of absorbing species
! PALB          : surface albedo
! PFBEAMR       : solar flux
!
! OUTPUTS:
! -------
! PAJVAL(KLEV,NDISS)    : photolysis rates
!
! LOCAL:
! -------
! INJLEV   : step to select levels from  BASCOE_J_MODULE above model top
!               if 0, select only top level
! ZTP
! ZS
! ZCP
! ZHDENS
! ZJDENS
! ZJRATES
!


USE PARKIND1 ,          ONLY : JPIM     ,JPRB
USE YOMHOOK  ,          ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE , ONLY : nabspec, abs_molmass, mxwvn, &
            & air, O2abs, O3abs, NOabs, CO2abs, NO2abs
USE BASCOE_J_MODULE ,   ONLY : NDISS
USE BASCOE_J_EXT_MODULE,ONLY : NABSLAYER, LTHICK, TEMPER, DENS, HDENS
USE YOMCST,             ONLY : RG

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM),INTENT(IN)     :: KLEV         ! Number of model levels
  REAL(KIND=JPRB), INTENT(IN)       :: PZKM(KLEV)   ! Altitude of model levels (km)
  REAL(KIND=JPRB), INTENT(IN)       :: PSZA         ! Solar Zenith Angle (degree)
  REAL(KIND=JPRB), INTENT(IN)       :: PTP(KLEV)    ! Temperature (K)
  REAL(KIND=JPRB), INTENT(IN)       :: PDENSA(KLEV,NABSPEC)    ! densities
  REAL(KIND=JPRB), INTENT(IN)       :: PALB         ! Surface albedo
  REAL(KIND=JPRB), INTENT(IN)       :: PFBEAMR(mxwvn)   ! solar flux
  REAL(KIND=JPRB), DIMENSION(KLEV,NDISS), INTENT(OUT)  :: PAJVAL


!-----------------------------------------------------------------------
!*      0.2 Local Parameters
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM), PARAMETER :: INJLEV   = 5     ! selection of level above top
  REAL(KIND=JPRB), PARAMETER    :: ZR0            = 6371.0E0    ! Effective earth radius (km)
  REAL(KIND=JPRB), PARAMETER    :: ZRgas         = 8.3144621   ! Perfect gas constant (J/K/mol)

!-----------------------------------------------------------------------
!*      0.3 Local variables
!---------------------------------------------------------------------
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  ! variables related to TUV levels
  INTEGER(KIND=JPIM)    :: IABOVE       ! 0-based index of lowest level above model
  INTEGER(KIND=JPIM)    :: ILEV, JLEV   ! 0-based number of levels, and counter

  REAL(KIND=JPRB)   :: ZS0              ! top height (km)
  REAL(KIND=JPRB)   :: ZABOVE           ! thickness above model (km)
  REAL(KIND=JPRB)   :: ZG               ! G at level (cm/s2)
  REAL(KIND=JPRB)   :: ZH               ! 
  REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE :: ZS          ! level heights (km)
  REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE :: ZTP         ! temperature (K)
  REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE :: ZCP       ! specific heat of air
  REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE :: ZDENS     ! densities (of abs. species)
  REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE :: ZHDENS    ! scaled height (of abs. species)
  REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE :: ZJRATES   ! computed photolysis rates



!-----------------------------------------------------------------------
#include "bascoe_tuv.intfb.h"

  IF (LHOOK) CALL DR_HOOK('BASCOE_J_CALC',0,ZHOOK_HANDLE )

  ! ----------------------------------------------------------------------------------
  ! Initializations for BASCOE_TUV
  !     ...define levels
  ! ----------------------------------------------------------------------------------
  ZS0 = NABSLAYER * LTHICK        ! height of absorption top level
  ZABOVE = ZS0 - PZKM(1)    ! thickness above model
  IF (ZABOVE < LTHICK) THEN
    ! Model top is too high to use climatological data defined in BASCOE_J_MODULE
    !   Use only the model levels to compute photolysis rates using TUV
    IABOVE = -1
  ELSEIF (INJLEV /= 0) THEN
    ! Group absorption layers above model top by number 'INJLEV'
    IABOVE = INT( ZABOVE / REAL(INJLEV, KIND=JPRB) )
    IF (ZABOVE - REAL(IABOVE*INJLEV, KIND=JPRB) < LTHICK) THEN
        ! Layer above is thick enough
        IABOVE = IABOVE -1
    ENDIF
  ELSE
    ! use a single absorption layer above model top
    IABOVE = 0
  ENDIF

  ILEV = IABOVE + KLEV

  ALLOCATE( ZS(0:ILEV),             &
        &   ZTP(0:ILEV),            &
        &   ZCP(0:ILEV,NABSPEC),    &
        &   ZDENS(0:ILEV,NABSPEC),  &
        &   ZHDENS(0:ILEV,NABSPEC), &
        &   ZJRATES(NDISS,0:ILEV) )
  ZCP(0:ILEV,1:NABSPEC) = 1.005_JPRB

  ! ----------------------------------------------------------------------------------
  !     ...levels above model top
  ! ----------------------------------------------------------------------------------
  IF (IABOVE >= 0) THEN
    ! there is at least one level above model top:
    !   fill with climatological data defined in BASCOE_J_MODULE
    ZS(0)       = ZS0
    ZTP(0)      = TEMPER(0)
    ZDENS(0,1:NABSPEC)    = DENS(0,1:NABSPEC)
    ZHDENS(0,1:NABSPEC)   = HDENS(0,1:NABSPEC)
    IF (INJLEV /= 0) THEN
      ! other levels above model top
      DO JLEV = 1,IABOVE
        ZS(JLEV)       = ZS(JLEV-1) - REAL(INJLEV, KIND=JPRB)
        ZTP(JLEV)      = TEMPER(JLEV*INJLEV)
        ZDENS(JLEV,1:NABSPEC)    = DENS(JLEV*INJLEV,1:NABSPEC)
        ZHDENS(JLEV,1:NABSPEC)   = HDENS(JLEV*INJLEV,1:NABSPEC)
      ENDDO
    ENDIF
  ENDIF

  ! ----------------------------------------------------------------------------------
  !     ...model levels
  ! ----------------------------------------------------------------------------------
  DO JLEV = IABOVE+1, ILEV
    ZS(JLEV)    = PZKM(JLEV-IABOVE)
    ZTP(JLEV)   = PTP(JLEV-IABOVE)
    ZDENS(JLEV,:) = PDENSA(JLEV-IABOVE,:)
    ZG = 1.e2_JPRB * RG * ( ZR0 / (ZR0+zs(JLEV)) )**2.         ! cm/s2
    ZH = 1.e7_JPRB * ZRgas * ZTP(JLEV) / ZG                !
    ZHDENS(JLEV,O2abs)  = ZH / abs_molmass(O2abs)
    ZHDENS(JLEV,O3abs)  = ZH / abs_molmass(O3abs)
    ZHDENS(JLEV,NOabs)  = ZH / abs_molmass(NOabs)
    ZHDENS(JLEV,CO2abs) = ZH / abs_molmass(CO2abs)
    ZHDENS(JLEV,NO2abs) = ZH / abs_molmass(NO2abs)
  ENDDO

  ! ----------------------------------------------------------------------------------
  !     ...obtain J rates using TUV over levels above model top + model levels
  ! ----------------------------------------------------------------------------------
  ! YC: ignore level of turbopause as in bascoe j online ..??
  CALL BASCOE_TUV( ILEV, PALB, PFBEAMR=PFBEAMR, &
   &               P_HDENS=ZHDENS, PDENS=ZDENS, PTEMPER=ZTP, PSZA=PSZA, PZS=ZS, P_CP=ZCP, &
   &               PDRAT=ZJRATES )

  ! ----------------------------------------------------------------------------------
  !     ...transfer J rates on model levels to output
  ! ----------------------------------------------------------------------------------
  DO JLEV = 1, KLEV
    PAJVAL(JLEV,1:NDISS) = ZJRATES(1:NDISS,JLEV + IABOVE)
  ENDDO

  DEALLOCATE(ZS,ZTP,ZCP,ZDENS,ZHDENS,ZJRATES)


  IF (LHOOK) CALL DR_HOOK('BASCOE_J_CALC',1,ZHOOK_HANDLE )

END SUBROUTINE BASCOE_J_CALC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
