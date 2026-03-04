! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

subroutine varwrite1c_nc (kncid, nlev, kstart, kcount, varname, pvar)
!-----------------------------------------------------------------------
! Purpose: Write 0-dim or 1-dim variable in NetCDF file.
!
! Variables: kncid   - NetCDF file ID
!            nlev    - size of written array (1 for 0-dim)
!            kstart  - starting index(ices) of NetCDF file variable
!            kcount  - # of values written
!            varname - variable name
!            pvar    - variable
!
! Attention: in SCM pvar is defined as REAL(KIND=JPRB), but netcdf
!            uses different REAL!
!
!        M. Ko"hler     11-2000
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPRM

implicit none

!---global variables
character (len = *) :: varname
integer(KIND=JPIM)  :: kncid, nlev, kstart(min(nlev,2)), kcount(min(nlev,2))
REAL(KIND=JPRB)     :: pvar(nlev)

!---local variables
real(KIND=JPRM)     :: zvar1d_temp(nlev) !must be single precision (no f90 -r8)
integer(KIND=JPIM)  :: ivarid, istatus

#include "netcdf.inc"

!-----------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
!-----------------------------------------------------------------------

!write(*,*) 'Writing to NetCDF (', nlev, ' levels):  ', varname, &
!           '         mn ', sum(pvar)/nlev
!write(*,*) pvar

zvar1d_temp = pvar             !convert REAL*8 to NetCDF REAL format

istatus = NF_INQ_VARID    (kncid, varname, ivarid)
call handle_err_nc(istatus)
istatus = NF_PUT_VARA_REAL(kncid,ivarid,kstart,kcount, zvar1d_temp)     
call handle_err_nc(istatus)

end
