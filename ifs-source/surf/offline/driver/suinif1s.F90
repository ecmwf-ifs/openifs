SUBROUTINE SUINIF1S(LOPLEFT)
USE PARKIND1  ,ONLY : JPIM     ,JPRB,  JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMCT01S , ONLY : NSTART
USE YOMDPHY  , ONLY : NPOI     ,NTRAC
USE YOMFORC1S, ONLY : GFOR     ,JPSTPFC  ,JPTYFC   ,UFI     ,&
                      &VFI      ,TFI      ,QFI      ,PSFI    ,&
                      &SRFFI    ,TRFFI    ,R30FI    ,S30FI   ,&
                      &R30FI_C  ,S30FI_C, CO2FI
USE YOMLUN1S , ONLY : NULOUT
USE YOMLOG1S , ONLY : CFFORC   ,CFOUT,CFINIT,CFSURF
USE YOMGF1S  , ONLY : UNLEV0   ,VNLEV0   ,TNLEV0   ,QNLEV0   ,CNLEV0 ,&
                      &PNLP0    ,UNLEV1   ,VNLEV1   ,TNLEV1   ,&
                      &QNLEV1   ,CNLEV1   ,FSSRD    ,FSTRD    ,FLSRF    ,&
                      &FCRF     ,FLSSF    ,FCSF , PNLP1
USE MPL_MODULE

! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUINIF1S*  - Initialize the fields of the one column surface model.

!     Purpose.
!     --------
!           Initialize the fields of the one column surface model.

!***  Interface.
!     ----------
!        *CALL* *SUINIF1S*

!        Explicit arguments :
!        --------------------
!         LOPLEFT: FLAG INDICATING WHETHER MORE POINTS NEED TO BE DONE

!        Implicit arguments :
!        --------------------
!         NONE

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by CNT31S

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the one column model


!     Author.
!     -------
!        Jean-Francois Mahfouf and Pedro Viterbo   *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-03-21
!        Bart vd HURK (KNMI) Multi-column use of netCDF and restart utility
!        Anna Agusti-Panareda 2020-11-17 Atmospheric CO2 forcing

IMPLICIT NONE

INTEGER :: NCID,IERR ! for netcdf interface

INTEGER(KIND=JPIM) :: MUFI,MVFI,MTFI,MQFI,MPSFI,MSRFFI,MTRFFI,MR30FI,MS30FI, &
   &   MR30CFI,MS30CFI,MCO2FI

LOGICAL LOPLEFT ! FLAG FOR RUNNING MORE GRID POINTS
INTEGER(KIND=JPIM) :: MYPROC, NPROC
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "sudyn1s.intfb.h"
#include "suct01s.intfb.h"
#include "rdres.intfb.h"
#include "sugc1s.intfb.h"
#include "sugpd1s.intfb.h"
#include "sugp1s.intfb.h"
#include "sufcdf.intfb.h"
#include "sugdi1s.intfb.h"
#include "sucdh1s.intfb.h"

      
#include "netcdf.inc"

IF (LHOOK) CALL DR_HOOK('SUINIF1S',0,ZHOOK_HANDLE)


MYPROC = MPL_MYRANK()
NPROC  = MPL_NPROC()

!*       SETUP OF TIME VARIABLES.

CALL SUDYN1S (NULOUT)

!*       RUN AND OUTPUT CONTROLS.

CALL SUCT01S (NULOUT)

!*       READ RESTART

IF(NSTART /= 0)THEN

  IF(CFOUT == 'netcdf')THEN
    CALL RDRES(LOPLEFT)
  ELSE
    WRITE(NULOUT,*)'RESTART FILE REQUIRES NETCDF CONFIGURATION'
    CALL ABOR1('RESTART FILE REQUIRES NETCDF CONFIGURATION')
  ENDIF

ELSE

!*       SETUP OF GEOGRAPHICAL CONSTANTS.

!* Open surface climate field

  IF (CFSURF == 'netcdf')THEN
    IF( MYPROC == 1 ) NCID = NCOPN('surfclim', NCNOWRIT, IERR)
    WRITE(NULOUT,*)'NETCDF-FILE surfclim OPENED ON UNIT ',NCID
  ENDIF
  CALL SUGC1S (NCID,LOPLEFT)

!*       SETUP OF SOIL/VEGETATION CONSTANTS.

  CALL SUGPD1S (NCID)
  IF (CFSURF == 'netcdf')THEN
    IF( MYPROC == 1 ) CALL NCCLOS(NCID,IERR)
  ENDIF

!*       INITIALISATION OF SOIL PROGNOSTIC VARIABLES.

!* Open surface prognostics field

  IF (CFINIT == 'netcdf')THEN
    IF( MYPROC == 1 ) NCID = NCOPN('soilinit', NCNOWRIT, IERR)
    WRITE(NULOUT,*)'NETCDF-FILE soilinit OPENED ON UNIT ',NCID
  ENDIF
  CALL SUGP1S (NCID)
  IF (CFINIT == 'netcdf')THEN
    IF( MYPROC == 1 ) CALL NCCLOS(NCID,IERR)
  ENDIF

ENDIF

!*       INITIALISATION OF THE ATMOSPHERIC FORCING.

!*       ALLOCATION OF INPUT FORCING
ALLOCATE (GFOR(NPOI,JPSTPFC,JPTYFC))

MUFI=1
MVFI=MUFI+1
MTFI=MVFI+1
MQFI=MTFI+1
MPSFI=MQFI+1
MSRFFI=MPSFI+1
MTRFFI=MSRFFI+1
MR30FI=MTRFFI+1
MS30FI=MR30FI+1
MR30CFI=MS30FI+1
MS30CFI=MR30CFI+1
MCO2FI=MS30CFI+1

IF(MCO2FI /= JPTYFC)THEN
  WRITE(*,*)'NUMBER OF FORCING VARIABLES NOT COMPATIBLE: ',MCO2FI,JPTYFC
  CALL ABORT
ENDIF

UFI => GFOR(:,:,MUFI)
VFI => GFOR(:,:,MVFI)
TFI => GFOR(:,:,MTFI)
QFI => GFOR(:,:,MQFI)
PSFI => GFOR(:,:,MPSFI)
SRFFI => GFOR(:,:,MSRFFI)
TRFFI => GFOR(:,:,MTRFFI)
R30FI => GFOR(:,:,MR30FI)
S30FI => GFOR(:,:,MS30FI)
R30FI_C => GFOR(:,:,MR30CFI)
S30FI_C => GFOR(:,:,MS30CFI)
CO2FI => GFOR(:,:,MCO2FI)

IF (CFFORC == 'netcdf')THEN
  CALL SUFCDF
ELSE
  WRITE (NULOUT,"(' Forcing data set ',a,' not known ')")&
   &CFFORC(1:len_trim(CFFORC))
  CALL ABORT
ENDIF

ALLOCATE (UNLEV0(NPOI))
ALLOCATE (VNLEV0(NPOI))
ALLOCATE (TNLEV0(NPOI))
ALLOCATE (QNLEV0(NPOI))
ALLOCATE (CNLEV0(NPOI,NTRAC)) ! check what happens if ntrac = 0
ALLOCATE (PNLP0(NPOI))
ALLOCATE (PNLP1(NPOI))
ALLOCATE (UNLEV1(NPOI))
ALLOCATE (VNLEV1(NPOI))
ALLOCATE (TNLEV1(NPOI))
ALLOCATE (QNLEV1(NPOI))
ALLOCATE (CNLEV1(NPOI,NTRAC))
ALLOCATE (FSSRD(NPOI))
ALLOCATE (FSTRD(NPOI))
ALLOCATE (FLSRF(NPOI))
ALLOCATE (FCRF(NPOI))
ALLOCATE (FLSSF(NPOI))
ALLOCATE (FCSF(NPOI))

CNLEV0(:,:) = 0._JPRB
CNLEV1(:,:) = 0._JPRB 
!*       ALLOCATES THE DIAGNOSTIC VARIABLES.

CALL SUGDI1S
CALL SUCDH1S

CALL MPL_BARRIER()

IF (LHOOK) CALL DR_HOOK('SUINIF1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SUINIF1S
