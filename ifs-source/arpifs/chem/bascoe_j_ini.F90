! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_J_INI
!**   DESCRIPTION
!     ----------
!
!   initialize auxiliary variables (defined in module BASCOE_J_MODULE)
!       for IFS(BASCOE) stratospheric chemistry, using TUV
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_J_INI* IS CALLED FROM *BASCOE[TM5]_CHEM_INI*.

!
!     AUTHOR.
!     -------
!        Yves Christophe     *BIRA*

!-----------------------------------------------------------------------

USE YOMLUN              , ONLY : NULOUT
USE PARKIND1            , ONLY : JPIM     ,JPRB
USE YOMHOOK             , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_EXT_MODULE , ONLY : NABSLAYER, LTHICK, WMOLE, TEMPER, DENS, HDENS, &
            & DENS_AIR, DENS_O2
USE BASCOE_TUV_MODULE   , ONLY : nabspec, ABS_MOLMASS, abs_name, &
            & air, O2abs, O3abs, NOabs, CO2abs, NO2abs
USE YOMCST,             ONLY : RG

IMPLICIT NONE

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  CHARACTER(LEN=*), PARAMETER   :: CL_MY_NAME   = 'BASCOE_J_INI'

  REAL(KIND=JPRB), PARAMETER    :: ZR0            = 6371.0d0    ! Effective earth radius (km)
  REAL(KIND=JPRB), PARAMETER    :: ZRGAS          = 8.3144621   ! Perfect gas constant (J/K/mol)
  REAL(KIND=JPRB), PARAMETER    :: ZKM_TURBO     = 97.         ! altitude of turbopause

  INTEGER(KIND=JPIM)            :: ILEVTURBO
  INTEGER(KIND=JPIM)            :: ILEV, IABS

  REAL(KIND=JPRB), DIMENSION(0:NABSLAYER) :: ZS, ZG, ZHM
  REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE

  CHARACTER(LEN=128)            :: CL_FILENM

#include "bascoe_tuv_ini.intfb.h"
#include "bascoe_read_coldat.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_J_INI',0,ZHOOK_HANDLE )

    !-----------------------------------------------------------------------
    !   ... Run BASCOE_TUV_INI i.e. initialize the photolysis code
    !-----------------------------------------------------------------------
    CALL BASCOE_TUV_INI( )

    !-----------------------------------------------------------------------
    !   ... Vertical profiles zs, zg; dens(O2), dens(air) (from MSIS)
    !-----------------------------------------------------------------------

    ILEVTURBO = 0
    ZS(0) = NABSLAYER * LTHICK
    DO ILEV = 1, NABSLAYER
       ZS(ILEV) = ZS(ILEV-1) - LTHICK
       IF( ZS(ILEV) > ZKM_TURBO ) ILEVTURBO = ILEV
    ENDDO
    WRITE(nulout,*) CL_MY_NAME // ': altitude grid: ZS(0)= ',ZS(0) &
   &   ,' to ZS(NABSLAYER)= ',ZS(NABSLAYER),' by ',LTHICK,' km'
    WRITE(nulout,*) CL_MY_NAME // ': alt of turbopause(km) : ',ZS(ILEVTURBO)
    ZG(:) = 1.e2_JPRB * RG * ( ZR0 / (ZR0+ZS(:)) )**2                ! cm/s2
    ZHM(:) = 1.e7_JPRB * ZRGAS * temper(:) / ZG(:)                  !
    DENS(:,O2abs) = DENS_O2
    DENS(:,air) = DENS_AIR

    !-----------------------------------------------------------------------
    !   ... read nb density profiles of other absorbing species
    !-----------------------------------------------------------------------
    DO IABS = O3abs, nabspec
       CL_FILENM = TRIM(ADJUSTL( abs_name(IABS) )) // '_vertprof.dat'
       WRITE(NULOUT,*) CL_MY_NAME // ': reading >' // TRIM( CL_FILENM ) // '<'
       CALL BASCOE_READ_COLDAT( CL_FILENM, NABSLAYER+1 &
   &                   , ZS(NABSLAYER:0:-1) &
   &                   , DENS(NABSLAYER:0:-1,IABS), 'asbndry' )
    ENDDO

    !-----------------------------------------------------------------------
    !       ... Scale heights
    !-----------------------------------------------------------------------
    Hdens(:,air) = ZHM(:) / WMOLE(:)
    WRITE(NULOUT,*) CL_MY_NAME // ': air scale height at surface (cm): ', &
   &                                                  Hdens(NABSLAYER,air)
    Hdens(:,O2abs)  = ZHM(:) / ABS_MOLMASS(O2abs)
    Hdens(:,O3abs)  = ZHM(:) / ABS_MOLMASS(O3abs)
    Hdens(:,NOabs)  = ZHM(:) / ABS_MOLMASS(NOabs)
    Hdens(:,CO2abs) = ZHM(:) / ABS_MOLMASS(CO2abs)
    Hdens(:,NO2abs) = ZHM(:) / ABS_MOLMASS(NO2abs)


IF (LHOOK) CALL DR_HOOK('BASCOE_J_INI',1,ZHOOK_HANDLE )

END SUBROUTINE BASCOE_J_INI

