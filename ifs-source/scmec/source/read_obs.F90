! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE READ_OBS(YDDIMV,NT, U, V, T, Q)

!     ------------------------------------------------------------------
!     Purpose.
!     --------
!        Reads u-winv, v-wind, temperature, specific humidity at step nt.
!
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!     ------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPRM

IMPLICIT NONE

TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
REAL(KIND=JPRB)    :: U(YDDIMV%NFLEVG), V(YDDIMV%NFLEVG), T(YDDIMV%NFLEVG), Q(YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: NT, INCID, VARID, ISTATUS, START2(2), COUNT2(2)
REAL(KIND=JPRM)    :: TEMP1A(YDDIMV%NFLEVG)

#include "netcdf.inc"

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
!     ------------------------------------------------------------------


!        1.    OPEN INPUT FILE.

ISTATUS = NF_OPEN ('scm_in.nc', NF_NOWRITE, INCID)
CALL HANDLE_ERR_NC(ISTATUS)


!        2.    READ TIME STEP.

START2 = (/ 1     , NT /)
COUNT2 = (/ YDDIMV%NFLEVG, 1  /)

ISTATUS = NF_INQ_VARID     (INCID, 'u', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
CALL HANDLE_ERR_NC(ISTATUS)
U(1:YDDIMV%NFLEVG) = TEMP1A

ISTATUS = NF_INQ_VARID     (INCID, 'v', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
CALL HANDLE_ERR_NC(ISTATUS)
V(1:YDDIMV%NFLEVG) = TEMP1A

ISTATUS = NF_INQ_VARID     (INCID, 't', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
CALL HANDLE_ERR_NC(ISTATUS)
T(1:YDDIMV%NFLEVG) = TEMP1A

ISTATUS = NF_INQ_VARID     (INCID, 'q', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
CALL HANDLE_ERR_NC(ISTATUS)
Q(1:YDDIMV%NFLEVG) = TEMP1A


!        3.    CLOSE INPUT FILE.

ISTATUS = NF_CLOSE (INCID)
CALL HANDLE_ERR_NC(ISTATUS)

END SUBROUTINE READ_OBS
