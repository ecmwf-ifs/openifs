! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

subroutine handle_err_nc(kstatus)

!**** *handle_err_nc*  - Handle NetCDF errors

!     Purpose.
!     --------
!     Handle NetCDF errors

!**   Interface.
!     ----------
!        *CALL* *handle_err_nc

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        Taken and adopted from NetCDF documentation Version 3.

!     Author.
!     -------
!        Martin Koehler

!     Modifications.
!     --------------
!        Original    00-09-12
!        M. Ko"hler  6-6-2006  Single Column Model integration within IFS 

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KSTATUS

#include "netcdf.inc"

!     ------------------------------------------------------------------

IF (KSTATUS .NE. NF_NOERR) THEN
  PRINT *, 'NETCDF ERROR:  ', NF_STRERROR(KSTATUS)
  KSTATUS = KSTATUS/0.0              !OPTIONAL CORE FILE PRODUCTION
  STOP 'STOPPED WITH NETCDF ERROR'
ENDIF

END
