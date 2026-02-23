! (C) Copyright 1995- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE CNTEND

!**** *CNTEND*  - Closes netCDF datafiles.

!     Purpose.
!     --------
!          Data file closure after finishing the run.

!***  Interface.
!     ----------
!        *CALL* *CNTEND

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd Hurk *KNMI*

!     Modifications.
!     --------------
!        Original      : 00/07/14
!        Martin Koehler: 00/09/10 adoption to one column model
!        M. Ko"hler      6-6-2006 Single Column Model integration within IFS 

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMLOG1C , ONLY : OUTFORM  ,NPOSPRG, NPOSDIA, NPOSDIA2

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: JPNCDF=3
INTEGER(KIND=JPIM) :: IPOS(JPNCDF), NPOS, J, ISTATUS

#include "netcdf.inc"

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
!     ------------------------------------------------------------------

IF (OUTFORM .EQ. 'netcdf') THEN

  IPOS(1) = NPOSPRG
  IPOS(2) = NPOSDIA
  IPOS(3) = NPOSDIA2
  DO J = 1,JPNCDF
    NPOS    = IPOS(J)
    istatus = NF_CLOSE(npos)
    call handle_err_nc(istatus)
  ENDDO

ENDIF

WRITE(*,*) 'Single Column Model Completed Successfully.'

END SUBROUTINE CNTEND
