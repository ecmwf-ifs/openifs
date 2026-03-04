! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SURIP1C(KULOUT)

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMCST , ONLY : RDAY
USE YOMRIP0  , ONLY : NSSSSS, NINDAT, RTIMST
USE YOMLOG1C
USE YOMCT0   , ONLY : LSFCFLX

#ifdef DOC

!**** *SURIP1C * - Routine to initialize NSSSSS & NINDAT in SCM

!     Purpose.
!     --------
!           Initialize the common YOMRIP

!**   Interface.
!     ----------
!        *CALL* *SURIP1C(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit of the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMRIP
!        COMMON YOMARG

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
!        Martin Koehler  *ECMWF*

!     Modifications.
!     --------------
!        Original    2006-04-27
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: INCID, VARID, ISTATUS, ITEMP, NT
INTEGER(KIND=JPIM) :: IA, ID, IM

#include "netcdf.inc"
#include "handle_err_nc.intfb.h"
#include "nam1c.nam.h"      
#include "fcttim.func.h"

!     ------------------------------------------------------------------

!*       1.    Initialize YOMRIP.
!              ------------------

IF (LHOOK) CALL DR_HOOK('SURIP1C',0,ZHOOK_HANDLE)

!        1.1 Initial time from scm input file

REWIND(NULNAM)
READ(NULNAM,NAM1C)   ! read nstrtini temporarily

NT = NSTRTINI
WRITE(*,*) 'NSTRTINI',NSTRTINI

ISTATUS = NF_OPEN           ('scm_in.nc', NF_NOWRITE, INCID)
CALL HANDLE_ERR_NC(ISTATUS)

ISTATUS = NF_INQ_VARID      (INCID, 'date', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_INT   (INCID, VARID, NT, 1, ITEMP)
CALL HANDLE_ERR_NC(ISTATUS)
NINDAT = ITEMP

ISTATUS = NF_INQ_VARID      (INCID, 'second', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_INT   (INCID, VARID, NT, 1, ITEMP)
CALL HANDLE_ERR_NC(ISTATUS)
NSSSSS = ITEMP

ISTATUS = NF_CLOSE          (INCID)
CALL HANDLE_ERR_NC(ISTATUS)

ID=NDD(NINDAT)
IM=NMM(NINDAT)
IA=NCCAA(NINDAT)
RTIMST=RTIME(IA,IM,ID,NSSSSS,RDAY)

!      ----------------------------------------------------------------

!*       2.    PRINT PART OF YOMRIP.
!              ---------------------

WRITE(UNIT=KULOUT,FMT='('' SCM: NINDAT='',I10,'' NSSSSS='',I6)')NINDAT,NSSSSS

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURIP1C',1,ZHOOK_HANDLE)
END SUBROUTINE SURIP1C
