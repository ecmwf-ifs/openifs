! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU1C( YDMODEL )

USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  , ONLY : JPIM, JPRB, JPRM
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0    , ONLY : LSFCFLX, LROUGH, REXTZ0M, REXTZ0H
USE YOMLUN    , ONLY : NULNAM
USE YOMLOG1C

#ifdef DOC

!**** *SU1C*  - Initialize common YOMLOG1C controlling SCM

!     Purpose.
!     --------
!        Initialize YOMLOG1C, the common that includes 
!        the basic switches for the single column model.

!**   Interface.
!     ----------
!        *CALL* *SU1C from SU0YOM1C

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the SCM

!     Author.
!     -------
!        Joao Teixeira   *ECMWF*

!     Modifications.
!     --------------
!        Original      94-04-27
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM) :: INCID, VARID, ISTATUS, UNLIMDIMID, DIMLEN
REAL(KIND=JPRM)    :: TIME(2)
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "netcdf.inc"
#include "nam1c.nam.h"      

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU1C',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDMODEL%YRML_GCONF%YRRIP%TSTEP)

!*       1.    Set default values.
!              -------------------


!*       1.1   Set default values for the one-column model switches.
!              -----------------------------------------------------


LDYNFOR = .TRUE.
LUGVG   = .TRUE.
LVERVEL = .TRUE.
LETADOT = .FALSE.
LUPWIND = .TRUE.
LWADV   = .TRUE.
!LWADVCLD= .FALSE.
LWADVCLD= .TRUE.
LTADV   = .FALSE.
LQADV   = .FALSE.
LUVADV  = .FALSE.
LVARFOR = .FALSE.
LVARSST = .FALSE.
LRELAX  = .FALSE.
LUVREL  = .TRUE.
LTQREL  = .TRUE.
RUVREL_PLEV = 2000.E2_JPRB ! 2000hPa ensures all levels relaxed by default
RTQREL_PLEV = 2000.E2_JPRB ! 2000hPa ensures all levels relaxed by default
LSFCFLX = .FALSE.
LROUGH  = .FALSE.
REXTZ0M = 0.0_JPRB 
REXTZ0H = 0.0_JPRB

NPOSASC = 50

OUTFORM = 'netcdf'
CMODID  = 'model code identification'
CSIMID  = 'simulation identification'

RDTRELAX= 10800._JPRB


! --- forcing frequency diagnosed from input file ---
istatus = NF_OPEN            ('scm_in.nc', nf_nowrite, incid)
call handle_err_nc(istatus)

istatus = NF_INQ_UNLIMDIM    (incid, unlimdimid)
call handle_err_nc(istatus)
istatus = NF_INQ_DIMLEN      (incid, unlimdimid, dimlen)
call handle_err_nc(istatus)

if ( dimlen > 1 ) then

  istatus = NF_INQ_VARID     (incid, 'time', varid)
  call handle_err_nc(istatus)
  istatus = NF_GET_VARA_REAL (incid, varid, 1, 2, time)
  call handle_err_nc(istatus)

  if ( mod ( time(2) - time(1) , real(tstep) ) == 0 ) then
    NFRFOR = ( time(2) - time(1) ) / tstep
  else
    write(*,*) 'forcing time interval is NO multiple of time-step!'
!   stop
  endif

else

  NFRFOR = 1

endif

istatus = NF_CLOSE           (incid)
call handle_err_nc(istatus)

NFRSST  = NFRFOR
NFROBS  = NFRFOR
NSTRTINI= 1


!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------


REWIND(NULNAM)
READ(NULNAM,NAM1C)


!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SU1C',1,ZHOOK_HANDLE)
END SUBROUTINE SU1C
