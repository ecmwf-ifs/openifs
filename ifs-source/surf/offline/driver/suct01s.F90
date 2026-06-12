SUBROUTINE SUCT01S(KULOUT)
USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN1S , ONLY : NULNAM
USE YOMDYN1S , ONLY : TSTEP
USE YOMCT01S , ONLY : JPNPST   ,NPOSTS   ,NHISTS   ,&
            &CNMEXP   ,NFRPLT   ,NCYCLE   ,NSTART   ,NSTOP    ,&
            &NFRPOS   ,NFRHIS   ,LNF      ,LMPLOT   ,NFRRES, LSCMEC


#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUCT01S*   - Routine to initialize level 0 control common

!     Purpose.
!     --------
!           Initialize level 0 control commons
!***  Interface.
!     ----------
!        *CALL* *SUCT01S(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMCT01S

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by SUINIF1S.

!     Reference.
!     ----------
!        ECMWF Research Department documentation 
!        of the one column surface model

!     Author.
!     -------
!        Jean-Francois Mahfouf and Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-03-21
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE
INTEGER(KIND=JPIM) :: KULOUT

INTEGER(KIND=JPIM) :: J
REAL(KIND=JPRB) :: ZUNIT

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "namct01s.h"

IF (LHOOK) CALL DR_HOOK('SUCT01S',0,ZHOOK_HANDLE)

!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

!     1.1 Set implicit default values
!     Cycle and  name of experiment 

NCYCLE=16
CNMEXP='0123'

!     Plotting outputs

LNF   =.TRUE.
LMPLOT=.FALSE.
NFRPLT=1

NSTART=0
NSTOP=480

!    Frequency of post-processing

NFRPOS=NSTOP
NFRHIS=NSTOP
NFRRES=0
DO J=0,JPNPST
  NPOSTS(J)=0
  NHISTS(J)=0
ENDDO

! ECMWF Single Column Model off by default
LSCMEC=.FALSE.



!      ----------------------------------------------------------------

!*       2.    Modifies default values.
!              ------------------------

!        2.1   Read namelist

REWIND(NULNAM)
READ(NULNAM,NAMCT01S)

!        2.3   Reset values and test

ZUNIT=3600.
IF ((TSTEP > 0.).AND.(NFRPOS < 0)) THEN
  NFRPOS=NINT((REAL(-NFRPOS,KIND=JPRD)*ZUNIT)/TSTEP)
ENDIF
IF ((TSTEP > 0.).AND.(NFRRES < 0)) THEN
  NFRRES=NINT((REAL(-NFRRES,KIND=JPRD)*ZUNIT)/TSTEP)
ENDIF
IF ((TSTEP > 0.).AND.(NFRHIS < 0)) THEN
  NFRHIS=NINT((REAL(-NFRHIS,KIND=JPRD)*ZUNIT)/TSTEP)
ENDIF

!      -----------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMCT01S '')')
WRITE(UNIT=KULOUT,FMT='('' LNF    = '',L2,'' LMPLOT = '',L2)') LNF, LMPLOT
WRITE(UNIT=KULOUT,FMT='('' NFRPLT = '',I6,'' NCYCLE  = '',I6)') NFRPLT,NCYCLE
WRITE(UNIT=KULOUT,FMT='('' CNMEXP = '',A16)')CNMEXP
WRITE(UNIT=KULOUT,FMT='('' NSTART = '',I6,'' NSTOP  = '',I6 &
     &,'' NFRPOS = '',I6 &
     &,'' NFRRES = '',I6 &
     &,'' NFRHIS = '',I6)')&
 &NSTART,NSTOP ,NFRPOS,NFRRES,NFRHIS

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUCT01S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SUCT01S
