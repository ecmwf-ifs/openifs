SUBROUTINE SUDYN1S(KULOUT)
USE PARKIND1  ,ONLY : JPIM     ,JPRB          ,JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMDYN1S , ONLY : NSTEP    ,TSTEP    ,TDT, NACCTYPE,LPREINT,LSWINT,LFLXINT,TCOUPFREQ
USE YOMLUN1S , ONLY : NULNAM
USE YOMRIP   , ONLY : RTIMTR   ,RTIMST

#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUDYN1S*   - Initialize constants and control for temporal evolution
!                   of the one-column surface model

!     Purpose.
!     --------
!           Initialize YOMDYN1S

!**   Interface.
!     ----------
!        *CALL* *SUDYN1S(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMDYN1S
!        COMMON YOMLUN1S

!     Method.
!     -------
!        See documentation

!     Reference.
!     ----------
!        ECMWF Research Department documentation 
!        of the one-column surface model

!     Author.
!     -------
!        Jean-Francois Mahfouf and Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        95-03-21
!      03/12/2013 : E. DUTRA ADD LOGICAL CONTROLING FORCING INTERPOLATION 


#endif
IMPLICIT NONE
INTEGER(KIND=JPIM) :: KULOUT

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "namdyn1s.h"

IF (LHOOK) CALL DR_HOOK('SUDYN1S',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

TSTEP=1800.
TDT  =TSTEP
NSTEP=0
NACCTYPE=0 ! Forcing fluxes assument centered on timestamp (linear interp.)
!            NACCTYPE=1 ! forward accumulation (starting at timestamp)
!            NACCTYPE=2 ! backward accumulation (ending at timestamp)

LPREINT=.FALSE. ! PRECIPITATION DISTRIBUTION OF FORCING 
LSWINT=.FALSE. ! SOLAR ANGLE INTERPOLATION 
LFLXINT=.FALSE. ! Linear interpolation of solar and thermal rad 
TCOUPFREQ=86400._JPRB ! Coupling frequency with cama-flood 
!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

!        2.1   Read namelist.
REWIND(NULNAM)
READ(NULNAM,NAMDYN1S)

!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

RTIMTR=RTIMST+0.5_JPRD*TSTEP
WRITE(UNIT=KULOUT,FMT='('' COMMON YOMDYN1S '')')
WRITE(UNIT=KULOUT,FMT='('' TSTEP  = '',E14.8 &
     &,'' TDT    = '',E14.8  )') TSTEP,TDT
WRITE(UNIT=KULOUT,FMT='('' NACCTYPE  = '',I5)') NACCTYPE
WRITE(UNIT=KULOUT,FMT='('' LPREINT  = '',L5)') LPREINT 
WRITE(UNIT=KULOUT,FMT='('' LSWINT  = '',L5)') LSWINT
WRITE(UNIT=KULOUT,FMT='('' LSWINT  = '',L5)') LFLXINT
WRITE(UNIT=KULOUT,FMT='('' TCOUPFREQ  = '',E14.8 )') TCOUPFREQ
!     ------------------------------------------------------------------

IF ( LSWINT .AND. NACCTYPE .NE. 2 ) THEN
    WRITE(UNIT=KULOUT,FMT='(A)') "ERROR! Solar angle interpolation only possible with NACCTYPE==2"
    STOP
ENDIF 
IF ( LPREINT.AND. NACCTYPE .NE. 2 ) THEN
    WRITE(UNIT=KULOUT,FMT='(A)') "ERROR! Precipitation distribution (LPREINT) only possible with NACCTYPE==2"
    STOP
ENDIF 
IF ( NACCTYPE .EQ. 1 ) THEN
    WRITE(UNIT=KULOUT,FMT='(A)') "ERROR/Warning NACCTYPE==1 has not been fully tested"
    WRITE(UNIT=KULOUT,FMT='(A)') "Remove stop from sudyn1s and check dtforc!"
    STOP
ENDIF 

IF (LHOOK) CALL DR_HOOK('SUDYN1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SUDYN1S
