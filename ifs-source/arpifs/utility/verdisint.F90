! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE VERDISINT(YDVFE,YDCVER,CDOPER,CDBC,KPROMA,KSTART,KPROF,KFLEV,PIN,POUT,KCHUNK)

!**** *VERDISINT*   VERtical DIScretization -
!                INTerface for finite element type vertical operations:
!                derivative or integral

!     Purpose.
!     --------
!          This subroutine prepares an interface to VERINT
!          computing vertical integral with respect to eta
!          and VERDER computing vertical derivative with 
!          respect to eta of a function given at full or 
!          half model levels using a general scheme.

!**   Interface.
!     ----------
!        *CALL* *VERDISINT(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          CDOPER    - type of integral or derivative applied:
!                      'ITOP' - integral from top
!                      'IBOT' - integral from bottom
!                      'INGW' - invertible integral operator
!                      'HDER' - first derivative at half levels
!                      'FDER' - first derivative on full levels
!                      'DDER' - second derivative on full levels
!                      'DEGW' - invertible derivative operator
!          CDBC      - boundary conditions used ('00','01','10','11', 
!                      first digit for top, second digit for bottom,
!                      0 for value prescribed, 1 for derivative prescribed)
!          KPROMA    - horizontal dimension.
!          KSTART    - first element of work.
!          KPROF     - depth of work.
!          KFLEV     - vertical dimension for array PIN.
!          PIN       - input field
!        OPTIONAL INPUT:
!          KCHUNK    - Use NPROMA as blocking factor when called outside 
!                      OpenMP threaded region

!        OUTPUT: 
!          POUT      - integral or derivative of PIN according to CDOPER

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        none

!     Reference.
!     ----------

!     Author.
!     -------
!        P. Smolikova (CHMI/LACE/ALADIN)

!     Modifications.
!     --------------
!        Original : Sep 2017
!        P.Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE PARKIND1,ONLY : JPIM, JPRB, JPRD
USE YOMCVER ,ONLY : TCVER
USE YOMHOOK ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  ,ONLY : NULERR
USE YOMVERT ,ONLY : TVFE

!     ---------------------------------------------------------------------------

IMPLICIT NONE
TYPE(TVFE),TARGET ,INTENT(IN) :: YDVFE
TYPE(TCVER)       ,INTENT(IN)    :: YDCVER
CHARACTER(LEN=*)  ,INTENT(IN) :: CDOPER, CDBC
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA, KSTART, KPROF, KFLEV
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KCHUNK
REAL(KIND=JPRB)   ,INTENT(IN) :: PIN(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT):: POUT(KPROMA,KFLEV+1) 

!     ------------------------------------------------------------------

CHARACTER(LEN=2)   :: CLBC
INTEGER(KIND=JPIM) :: ILEVIN, ILEVOUT, IND, IEND, ITYPE, JCHUNK
REAL(KIND=JPRD),CONTIGUOUS,POINTER :: ZOPER(:,:)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verder.intfb.h"
#include "verint.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VERDISINT',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*   1. Initialization and control
!    -----------------------------

IF (YDCVER%LVFE_ECMWF) THEN
  IF (CDOPER /= 'FDER'.OR.YDCVER%LVFE_NOBC) THEN
    CLBC='XX' ! no BC applied for derivatives
  ELSE
    CLBC='BC' ! BC according to Hortal applied on derivatives
  ENDIF
ELSE
  CLBC=CDBC
ENDIF


!*   2. Set operator size according to boundary conditions applied
!    --------------------------------------------------------------

IF (CDOPER=='FDER'.OR.CDOPER=='DDER') THEN
  ILEVOUT = KFLEV
ELSE
  ILEVOUT = KFLEV+1
ENDIF

IF (CDOPER=='IBOT') THEN
  ITYPE=1
ELSE
  ITYPE=0
ENDIF

IF (CDOPER=='INGW'.OR.CDOPER=='DEGW') THEN
  IND=0
  ILEVIN = KFLEV+1
ELSEIF (CLBC=='XX') THEN
  IND=1
  ILEVIN = KFLEV
ELSE
  IND=0
  ILEVIN = KFLEV+2
ENDIF
IEND=ILEVIN-1+IND

!*   3. Set operator according to boundary conditions applied
!    --------------------------------------------------------

!*   3.1 Vertical integral
!    ---------------------
IF (CDOPER=='INGW') THEN
  ZOPER => YDVFE%RINTGW(:,:)
ELSEIF (ANY(CDOPER == (/'ITOP','IBOT'/))) THEN
  IF (CLBC=='XX') THEN
    ZOPER => YDVFE%RINTE(:,:)
  ELSEIF (CLBC=='00') THEN
    ZOPER => YDVFE%RINTBF00(:,:)
  ELSEIF (CLBC=='11') THEN
    ZOPER => YDVFE%RINTBF11(:,:)
  ELSE
    WRITE(NULERR,*) 'VERDISINT: ITOP/IBOT NOT IMPLEMENTED FOR CDBC=',CLBC
    CALL ABOR1(' VERDISINT: ABOR1 CALLED')
  ENDIF
ENDIF

!*   3.2 Vertical derivative
!    -----------------------

IF (CDOPER=='DEGW') THEN
  ZOPER => YDVFE%RDERGW(:,:)
ELSEIF (CDOPER=='HDER') THEN
  IF (CLBC=='00') THEN
    ZOPER => YDVFE%RDERBH00(:,:)
  ELSEIF (CLBC=='01') THEN
    ZOPER => YDVFE%RDERBH01(:,:)
  ELSE
    WRITE(NULERR,*) 'VERDISINT: HDER NOT IMPLEMENTED FOR CDBC=',CLBC
    CALL ABOR1(' VERDISINT: ABOR1 CALLED')
  ENDIF
ELSEIF (CDOPER=='FDER') THEN
  IF (CLBC=='XX') THEN
    ZOPER => YDVFE%RDERI(:,:)
  ELSEIF (CLBC=='BC') THEN
    ZOPER => YDVFE%RDERB(:,:)
  ELSEIF (CLBC=='00') THEN
    ZOPER => YDVFE%RDERBF00(:,:)
  ELSEIF (CLBC=='01') THEN
    ZOPER => YDVFE%RDERBF01(:,:)
  ELSEIF (CLBC=='10') THEN
    ZOPER => YDVFE%RDERBF10(:,:)
  ELSEIF (CLBC=='11') THEN
    ZOPER => YDVFE%RDERBF11(:,:)
  ELSE
    WRITE(NULERR,*) 'VERDISINT: FDER NOT IMPLEMENTED FOR CDBC=',CLBC
    CALL ABOR1(' VERDISINT: ABOR1 CALLED')
  ENDIF
ELSEIF (CDOPER=='DDER') THEN
  IF (CLBC=='XX') THEN
    ZOPER => YDVFE%RDDERI(:,:)
  ELSEIF (CLBC=='01') THEN
    ZOPER => YDVFE%RDDERBF01(:,:)
  ELSE
    WRITE(NULERR,*) 'VERDISINT: DDER NOT IMPLEMENTED FOR CDBC=',CLBC
    CALL ABOR1(' VERDISINT: ABOR1 CALLED')
  ENDIF
ENDIF

!*   4. Apply the required operation
!    --------------------------------

IF (ANY(CDOPER == (/'INGW','ITOP','IBOT'/))) THEN
  IF(PRESENT(KCHUNK)) THEN
    JCHUNK=KCHUNK
  ELSE
    JCHUNK=1
  ENDIF
  CALL VERINT(KPROMA,KSTART,KPROF,ILEVIN,ILEVOUT,ZOPER,PIN(:,IND:IEND),&
   & POUT,ITYPE,KCHUNK=JCHUNK)
ELSEIF (ANY(CDOPER == (/'FDER','HDER','DEGW','DDER'/))) THEN
  CALL VERDER(KPROMA,KSTART,KPROF,ILEVIN,ILEVOUT,ZOPER,PIN(:,IND:IEND),POUT)
ELSE 
  WRITE(NULERR,*) 'VERDISINT: UNKNOWN CDOPER=', CDOPER
  CALL ABOR1(' VERDISINT: ABOR1 CALLED')
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VERDISINT',1,ZHOOK_HANDLE)
END SUBROUTINE VERDISINT
