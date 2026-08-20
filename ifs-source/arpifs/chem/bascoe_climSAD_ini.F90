! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_CLIMSAD_INI
!**   DESCRIPTION
!     ----------
!
!   Initialize and read from an ASCII input file the Surface Aerosol Densities (SAD) 
!       of stratospheric liquid sulfate aerosols for (BASCOE) stratospheric chemistry
!
!   Two input files are vailable and chosen (for now) through integer parameter IMODE_SAD below:
!
!      'ICBG_SAD_clim.txt' : a year-independent climatology averaged over 1997-2005 
!                              from a run by IFS-CB05-BASCOE-GLOMAP (ICBG) for CAMS_43
!                            This was the only option at the end of CAMS_42 (June 2021)
!
!      'CMIP6_strato-aerosol-SAD_1990-2018_v4.txt' : dataset developed and used for CMIP6, 
!                        from output of GCCM (SOCOL?) scaled to in-site obs (TBC) downloaded on 2021-11-23
!                         as CMIP_1850_2018_sad_v4.nc and converted to ASCII for year range 1990-2018.
!                   
!
!**   AUTHOR.
!     -------
!        Jonas Debosscher, 2021-11-30
!
!-----------------------------------------------------------------------

  USE PARKIND1           , ONLY : JPIM     ,JPRB
  USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
  USE YOMLUN             , ONLY : NULOUT, NULERR

  USE BASCOE_MODULE, ONLY : NSYEAR_CLIM, NEYEAR_CLIM, NYEAR_CLIM, NMONTH_CLIM, NLAT_CLIM,&
  &  NLEV_CLIM,SAD_CLIM,P_CLIM,LAT_CLIM

  IMPLICIT NONE

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  CHARACTER(LEN=*),   PARAMETER :: SAD_FLNM_ICBG='ICBG_SAD_clim.txt'
  CHARACTER(LEN=*),   PARAMETER :: SAD_FLNM_CMIP6='CMIP6_strato-aerosol-SAD_1990-2018_v4.txt'
  INTEGER(KIND=JPIM), PARAMETER :: IMODE_SAD = 2    ! 1: use ICBG (BE.a051); 2: use CMIP6 (BE.a050,BE.a052)
  INTEGER(KIND=JPIM), PARAMETER :: IUTMP=77
  REAL(KIND=JPHOOK)             :: ZHOOK_HANDLE
  INTEGER(KIND=JPIM)            :: IOS, JM, JL, JK, JY, IYEAR, IMONTH

  LOGICAL :: LLFILE_EXISTS


#include "abor1.intfb.h"
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_CLIMSAD_INI',0,ZHOOK_HANDLE)

IF( IMODE_SAD == 1_JPIM ) THEN
    NYEAR_CLIM=1
    OPEN(IUTMP,file=SAD_FLNM_ICBG, action='read', iostat=IOS)
    IF( IOS /= 0 ) THEN
        WRITE(NULOUT,'(a,i3,a)') 'BASCOE_CLIMSAD_INI:  failed (ios=',ios,') to open ' // SAD_FLNM_ICBG
        WRITE(NULERR,'(a,i3,a)') 'BASCOE_CLIMSAD_INI:  failed (ios=',ios,') to open ' // SAD_FLNM_ICBG
        CALL ABOR1(' error: unable to open ICBG SAD Climatology')
    ENDIF
    DO JK=1,10
        READ(IUTMP,*)
    ENDDO
    READ(IUTMP,*)NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM
    ALLOCATE(SAD_CLIM(1,NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM),&
            &  P_CLIM(1,NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM),LAT_CLIM(NLAT_CLIM))
    DO JM=1,NMONTH_CLIM
        DO JL=1,NLAT_CLIM
            IF(JM==1_JPIM)THEN
                READ(IUTMP,*)IMONTH,LAT_CLIM(JL)
            ELSE
                READ(IUTMP,*)
            ENDIF
            DO JK=1,NLEV_CLIM
                READ(IUTMP,*)P_CLIM(1,JM,JL,JK),SAD_CLIM(1,JM,JL,JK)
            ENDDO
        ENDDO
    ENDDO
    CLOSE(IUTMP)

    !Change units. P_CLIM is assumed to be read in hPa, but application in IFS is in Pa..
    P_CLIM(1,:,:,:) = P_CLIM(1,:,:,:) * 100_JPRB

ELSEIF( IMODE_SAD == 2_JPIM ) THEN
    OPEN(IUTMP,file=SAD_FLNM_CMIP6, action='read', iostat=IOS)
    IF( IOS /= 0 ) THEN
       WRITE(NULOUT,'(a,i3,a)') 'BASCOE_CLIMSAD_INI:  failed (ios=',ios,') to open ' // SAD_FLNM_CMIP6
       WRITE(NULERR,'(a,i3,a)') 'BASCOE_CLIMSAD_INI:  failed (ios=',ios,') to open ' // SAD_FLNM_CMIP6
       CALL ABOR1(' error: unable to open CMIP6 SAD dataset')
    ENDIF
    DO JK=1,12
        READ(IUTMP,*)
    ENDDO
    READ(IUTMP,*)NSYEAR_CLIM,NEYEAR_CLIM,NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM
    NYEAR_CLIM=NEYEAR_CLIM-NSYEAR_CLIM+1
    ALLOCATE( SAD_CLIM(NYEAR_CLIM,NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM), &
            & P_CLIM(NYEAR_CLIM,NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM),LAT_CLIM(NLAT_CLIM) )
    DO JY=1,NYEAR_CLIM
        DO JM=1,NMONTH_CLIM
            DO JL=1,NLAT_CLIM
                IF((JY==1_JPIM).and.(JM==1_JPIM))THEN
                    READ(IUTMP,*)IYEAR,IMONTH,LAT_CLIM(JL)
                ELSE
                    READ(IUTMP,*)
                ENDIF

                DO JK=1,NLEV_CLIM
                    READ(IUTMP,*)P_CLIM(JY,JM,JL,JK),SAD_CLIM(JY,JM,JL,JK)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    CLOSE(IUTMP)
    
    !Change units. P_CLIM is assumed to be read in hPa, but application in IFS is in Pa..
    P_CLIM(:,:,:,:) = P_CLIM(:,:,:,:) * 100_JPRB

ELSE
   WRITE(NULERR,'(a,i3,a)') 'BASCOE_CLIMSAD_INI: error: unknown switch value: ',IMODE_SAD
   CALL ABOR1(' error: unknown switch value')
ENDIF

IF (LHOOK) CALL DR_HOOK('BASCOE_CLIMSAD_INI',1,ZHOOK_HANDLE)

END SUBROUTINE BASCOE_CLIMSAD_INI
