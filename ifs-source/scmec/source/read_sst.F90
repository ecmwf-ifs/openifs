! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

subroutine read_sst (nt, zlat, zlon, zsst, zspl)
!     ------------------------------------------------------------------
!     Purpose.
!     --------
!        Reads latitude (zlat), longitude (zlon), SST (zsst) and 
!        log(surface pressure) (zspl) from scm_in.nc at time step nt.
!
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPRM

IMPLICIT NONE

real(KIND=JPRB)    :: zlat, zlon, zsst, zspl
integer(KIND=JPIM) :: nt, incid, varid, istatus, start1, count1
real(KIND=JPRM)    :: temp0

#include "netcdf.inc"

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
!     ------------------------------------------------------------------


!        1.    OPEN INPUT FILE.

istatus = NF_OPEN ('scm_in.nc', nf_nowrite, incid)
call handle_err_nc(istatus)


!        2.    READ TIME STEP.

start1 = nt
count1 = 1

istatus = NF_INQ_VARID     (incid, 'lat', varid)
call handle_err_nc(istatus)
istatus = NF_GET_VARA_REAL (incid, varid, start1, count1, temp0)
call handle_err_nc(istatus)
zlat = temp0

istatus = NF_INQ_VARID     (incid, 'lon', varid)
call handle_err_nc(istatus)
istatus = NF_GET_VARA_REAL (incid, varid, start1, count1, temp0)
call handle_err_nc(istatus)
zlon = temp0

istatus = NF_INQ_VARID     (incid, 'open_sst', varid)
call handle_err_nc(istatus)
istatus = NF_GET_VARA_REAL (incid, varid, start1, count1, temp0)
call handle_err_nc(istatus)
zsst = temp0

istatus = NF_INQ_VARID     (incid, 'ps', varid)
call handle_err_nc(istatus)
istatus = NF_GET_VARA_REAL (incid, varid, start1, count1, temp0)
call handle_err_nc(istatus)
zspl = log(temp0)


!        3.    CLOSE INPUT FILE.

istatus = NF_CLOSE (incid)
call handle_err_nc(istatus)

end subroutine read_sst
