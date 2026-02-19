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
!   Initialize and read SAD climatology file 
!   for use in (BASCOE) stratospheric chemistry

!**   AUTHOR.
!     -------
!        Jonas Debosscher, 2021-05-20
!
!-----------------------------------------------------------------------

  USE PARKIND1           , ONLY : JPIM     ,JPRB
  USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
  USE YOMLUN             , ONLY : NULOUT
  
  USE BASCOE_MODULE, ONLY : NMONTH_CLIM, NLAT_CLIM, NLEV_CLIM, MONTHNUM,&
  & SAD_CLIM,P_CLIM,LAT_CLIM
  
  IMPLICIT NONE
  
!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  CHARACTER(LEN=*),   PARAMETER :: CL_SADCLIMFILE='ICBG_SAD_clim.txt'
  INTEGER(KIND=JPIM), PARAMETER :: IUTMP=77
  REAL(KIND=JPHOOK)             :: ZHOOK_HANDLE
  INTEGER(KIND=JPIM)            :: JM, JL, JK
  
  LOGICAL :: LLFILE_EXISTS
  
  
#include "abor1.intfb.h"
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_CLIMSAD_INI',0,ZHOOK_HANDLE)
  
INQUIRE(FILE=CL_SADCLIMFILE, EXIST=LLFILE_EXISTS)
write(NULOUT,*) 'BASCOE_CLIMSAD_INI: SAD clim file found:',CL_SADCLIMFILE, LLFILE_EXISTS  
! Read SAD climatology file if it exists
IF(LLFILE_EXISTS)THEN
   OPEN(IUTMP,file=CL_SADCLIMFILE, action='read')
   DO JK=1,10
       READ(IUTMP,*)
   ENDDO
   READ(IUTMP,*)NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM
   ALLOCATE(SAD_CLIM(NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM),&
 &     P_CLIM(NMONTH_CLIM,NLAT_CLIM,NLEV_CLIM),LAT_CLIM(NLAT_CLIM))
   DO JM=1,NMONTH_CLIM
           DO JL=1,NLAT_CLIM
                   IF(JM==1_JPIM)THEN
                       READ(IUTMP,*)MONTHNUM,LAT_CLIM(JL)
                   ELSE
                       READ(IUTMP,*)
                   ENDIF
                   DO JK=1,NLEV_CLIM
                       READ(IUTMP,*)P_CLIM(JM,JL,JK),SAD_CLIM(JM,JL,JK)
                   ENDDO
           ENDDO
   ENDDO
   CLOSE(IUTMP)

   !Change units. P_CLIM is assumed to be read in hPa, but application in IFS is in Pa..
   P_CLIM(:,:,:) = P_CLIM(:,:,:) * 100_JPRB 
ELSE
      CALL ABOR1(' error: SAD Climatology file not found.')
ENDIF

  
IF (LHOOK) CALL DR_HOOK('BASCOE_CLIMSAD_INI',1,ZHOOK_HANDLE)

END SUBROUTINE BASCOE_CLIMSAD_INI
