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

SUBROUTINE VERDISINT(YDVFE,YDCVER,CDOPER,CDBC,KPROMA,KSTART,KPROF,KFLEV,PIN,POUT,KBC,KCHUNK)

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
!        OPTIONAL:
!          KBC       - determine the definition of boundary conditions used;
!                      KBC=NVFE_DERBC/INTBC if not present
!          KCHUNK    - OMP chunking stride

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
!        P. Dueben & M. Diamantakis: Dec 2019 Improvements for mass conservation 
!                                    in single precision (equivalent changes with CY45R1)
!        F. Vana  14-Jan-2020: Simplifying and unifying VFE execution to be always double
!        P. Gillies & F. Vana 22-Jan-2020: Bit reproducibility for OMP
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB    ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCVER   ,ONLY : TCVER
USE YOMLUN    ,ONLY : NULERR
USE YOMVERT   ,ONLY : TVFE

!     ---------------------------------------------------------------------------

IMPLICIT NONE
TYPE(TVFE)        ,INTENT(IN)    :: YDVFE
TYPE(TCVER)       ,INTENT(IN)    :: YDCVER
CHARACTER(LEN=4)  ,INTENT(IN)    :: CDOPER
CHARACTER(LEN=2)  ,INTENT(IN)    :: CDBC
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(KPROMA,KFLEV+1) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KBC
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KCHUNK
!     ------------------------------------------------------------------

CHARACTER(LEN=15)  :: CLOPER
INTEGER(KIND=JPIM) :: IFLEV_IN, IFLEV_OUT, IND, ITYPE, IBC, ICHUNK
REAL(KIND=JPRD)    :: ZIN(KPROMA,0:KFLEV+1)
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! This needs to be double precision
REAL(KIND=JPRD),ALLOCATABLE :: ZOPER(:,:) 

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verder.intfb.h"
#include "verint.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VERDISINT',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF(.NOT.YDCVER%LVERTFE) CALL ABOR1(' ! LVERTFE=.F. IN VERDISINT !!!')

IF (PRESENT(KBC)) THEN
  IBC=KBC
ELSEIF (CDOPER=='INGW'.OR.CDOPER=='ITOP'.OR.CDOPER=='IBOT') THEN
  IBC=YDCVER%NVFE_INTBC
ELSEIF (CDOPER=='DEGW'.OR.CDOPER=='FDER'.OR.CDOPER=='HDER'.OR.CDOPER=='DDER') THEN
  IBC=YDCVER%NVFE_DERBC
ENDIF

IF (IBC > 3) THEN
  WRITE(NULERR,*) 'VERDISINT: IBC=',IBC,' MAY NOT BE USED'
  CALL ABOR1(' VERDISINT: ABOR1 CALLED')
ENDIF

!*   1. Set operator size according to boundary conditions applied
!    --------------------------------------------------------------

!*   1.1 Vertical integral
!    ---------------------
IF (CDOPER=='INGW') THEN
  IND=0
  IFLEV_IN = KFLEV+1
  IFLEV_OUT = KFLEV+1
  IF (IBC==0) THEN
    ITYPE=1
  ELSE
    ITYPE=0
  ENDIF
ELSEIF (CDOPER=='ITOP') THEN
  IND=1
  IFLEV_IN = KFLEV
  IFLEV_OUT = KFLEV+1
  ITYPE=0
ELSEIF (CDOPER=='IBOT') THEN
  IND=1
  IFLEV_IN = KFLEV
  IFLEV_OUT = KFLEV+1
  ITYPE=1
ENDIF
IF (CDOPER=='INGW'.OR.CDOPER=='ITOP'.OR.CDOPER=='IBOT') THEN
  IF (IBC==1.OR.IBC==3) THEN
    IND=0
    IFLEV_IN = KFLEV+2
  ENDIF
ENDIF

!*   1.2 Vertical derivative
!    -----------------------
IF (CDOPER=='DEGW') THEN
  IF (IBC==1 .OR. IBC==3) THEN
    IFLEV_IN = KFLEV+2
  ELSE
    IFLEV_IN = KFLEV+1
  ENDIF
  IND=0
  IFLEV_OUT = KFLEV+1
ELSEIF (CDOPER=='HDER') THEN
  IF (IBC==1 .OR. IBC==3) THEN
    IND=0
    IFLEV_IN = KFLEV+2
  ELSE
    IND=1
    IFLEV_IN = KFLEV
  ENDIF
  IFLEV_OUT = KFLEV+1
ELSEIF (CDOPER=='FDER'.OR.CDOPER=='DDER') THEN
  IF (IBC==1 .OR. IBC==3) THEN
    IND=0
    IFLEV_IN = KFLEV+2
  ELSE
    IND=1
    IFLEV_IN = KFLEV
  ENDIF
  IFLEV_OUT = KFLEV
ENDIF

     
!*   1.3 Rewrite for special cases
!    -----------------------------
IF (CDBC=='00') THEN
  IND=0
  IFLEV_IN = KFLEV+2
ENDIF

! Allocate working array
ALLOCATE(ZOPER(IFLEV_OUT,IND:IFLEV_IN))

!*   2. Set operator according to boundary conditions applied
!    --------------------------------------------------------

!*   2.1 Vertical integral
!    ---------------------
IF (CDOPER=='INGW') THEN

  IF (IBC==0) THEN
    IF(ALLOCATED(YDVFE%RINTE))THEN
      ZOPER = YDVFE%RINTE 
    ELSE
      CALL ABOR1("(E) IN VERDISINT RINTE NOT ALLOCATED")
    ENDIF
    CLOPER='RINTE'
  ELSE
    IF(ALLOCATED(YDVFE%RINTGW))THEN
      ZOPER = YDVFE%RINTGW
    ELSE
      CALL ABOR1("(E) IN VERDISINT RINTGW NOT ALLOCATED")
    ENDIF
    CLOPER='RINTGW'
  ENDIF

ELSEIF (CDOPER=='ITOP'.OR.CDOPER=='IBOT') THEN

  IF (IBC==0) THEN
    IF(ALLOCATED(YDVFE%RINTE))THEN
      ZOPER = YDVFE%RINTE 
    ELSE
      CALL ABOR1("(E) IN VERDISINT RINTE NOT ALLOCATED")
    ENDIF
    CLOPER='RINTE'
  ELSEIF (IBC==1.OR.IBC==3) THEN
    IF (CDBC=='00') THEN
      IF(ALLOCATED(YDVFE%RINTBF00))THEN
        ZOPER = YDVFE%RINTBF00 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RINTBF00 NOT ALLOCATED")
      ENDIF
      CLOPER='RINTBF00'
    ELSE
      IF(ALLOCATED(YDVFE%RINTBF11))THEN
        ZOPER = YDVFE%RINTBF11 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RINTBF11 NOT ALLOCATED")
      ENDIF
      CLOPER='RINTBF11'
      ZIN(KSTART:KPROF,0)=0.0_JPRD
      ZIN(KSTART:KPROF,KFLEV+1)=0.0_JPRD
    ENDIF
  ELSEIF (IBC==2) THEN
    IF (CDBC=='00') THEN
      IF(ALLOCATED(YDVFE%RINTBF00))THEN
        ZOPER = YDVFE%RINTBF00 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RINTBF00 NOT ALLOCATED")
      ENDIF
      CLOPER='RINTBF00'
    ELSE
      IF(ALLOCATED(YDVFE%RINTBF11_IMPL))THEN
        ZOPER = YDVFE%RINTBF11_IMPL 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RINTBF11_IMPL NOT ALLOCATED")
      ENDIF
      CLOPER='RINTBF11_IMPL'
      ZIN(KSTART:KPROF,0)=0.0_JPRD
      ZIN(KSTART:KPROF,KFLEV+1)=0.0_JPRD
    ENDIF
  ENDIF

ENDIF

ZIN(KSTART:KPROF,1:KFLEV)=PIN(KSTART:KPROF,1:KFLEV)

!*   2.2 Vertical derivative
!    -----------------------

IF (CDOPER=='DEGW') THEN

  IF(ALLOCATED(YDVFE%RDERGW))THEN
    ZOPER = YDVFE%RDERGW 
  ELSE
    CALL ABOR1("(E) IN VERDISINT RDERGW NOT ALLOCATED")
  ENDIF
  CLOPER='RDERGW'

ELSEIF (CDOPER=='HDER') THEN

  IF(ALLOCATED(YDVFE%RDERBH00))THEN
    ZOPER = YDVFE%RDERBH00 
  ELSE
    CALL ABOR1("(E) IN VERDISINT RDERBH00 NOT ALLOCATED")
  ENDIF
  CLOPER='RDERBH00'
  IF (CDBC=='01') THEN
    IF (IBC==2) THEN
      IF(ALLOCATED(YDVFE%RDERBH01_IMPL))THEN
        ZOPER = YDVFE%RDERBH01_IMPL 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBH01_IMPL NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBH01_IMPL'
    ELSEIF (IBC==3) THEN
      IF(ALLOCATED(YDVFE%RDERBH01))THEN
        ZOPER = YDVFE%RDERBH01 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBH01 NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBH01'
    ELSE
      WRITE(NULERR,*) 'VERDISINT: HDER NOT IMPLEMENTED &
       & FOR CDBC=',CDBC,' AND IBC=',IBC
      CALL ABOR1(' VERDISINT: ABOR1 CALLED')
    ENDIF
  ELSEIF (CDBC/='00') THEN
    WRITE(NULERR,*) 'VERDISINT: HDER NOT IMPLEMENTED &
     & FOR CDBC=', CDBC
    CALL ABOR1(' VERDISINT: ABOR1 CALLED')
  ENDIF

ELSEIF (CDOPER=='FDER') THEN

  IF (IBC==0) THEN
    IF(ALLOCATED(YDVFE%RDERI))THEN
      ZOPER = YDVFE%RDERI 
    ELSE
      CALL ABOR1("(E) IN VERDISINT RDERI NOT ALLOCATED")
    ENDIF
    CLOPER='RDERI'
  ELSEIF (IBC==1.OR.IBC==3) THEN
    IF (CDBC=='00') THEN
      IF(ALLOCATED(YDVFE%RDERBF00))THEN
        ZOPER = YDVFE%RDERBF00 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBF00 NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBF00'
    ELSEIF (CDBC=='01') THEN
      IF(ALLOCATED(YDVFE%RDERBF01))THEN
        ZOPER = YDVFE%RDERBF01 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBF01 NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBF01'
    ELSEIF (CDBC=='10') THEN
      IF(ALLOCATED(YDVFE%RDERBF10))THEN
        ZOPER = YDVFE%RDERBF10 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBF10 NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBF10'
    ELSEIF (CDBC=='11') THEN
      IF(ALLOCATED(YDVFE%RDERBF11))THEN
        ZOPER = YDVFE%RDERBF11 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBF11 NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBF11'
    ENDIF
  ELSEIF (IBC==2) THEN
    IF (CDBC=='00') THEN
      WRITE(NULERR,*) 'VERDISINT: RDERBF00_IMPL NOT DEFINED IN SUVERTFE!!'
      CALL ABOR1(' VERDISINT: ABOR1 CALLED')
    ELSEIF (CDBC=='01') THEN
      IF(ALLOCATED(YDVFE%RDERBF01_IMPL))THEN
        ZOPER = YDVFE%RDERBF01_IMPL 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBF01_IMPL NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBF01_IMPL'
    ELSEIF (CDBC=='10') THEN
      IF(ALLOCATED(YDVFE%RDERBF10_IMPL))THEN
        ZOPER = YDVFE%RDERBF10_IMPL 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBF10_IMPL NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBF10_IMPL'
    ELSEIF (CDBC=='11') THEN
      IF(ALLOCATED(YDVFE%RDERBF11_IMPL))THEN
        ZOPER = YDVFE%RDERBF11_IMPL 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDERBF11_IMPL NOT ALLOCATED")
      ENDIF
      CLOPER='RDERBF11_IMPL'
    ENDIF
  ENDIF

ELSEIF (CDOPER=='DDER') THEN

  IF (IBC==0) THEN
    IF(ALLOCATED(YDVFE%RDDERI))THEN
      ZOPER = YDVFE%RDDERI 
    ELSE
      CALL ABOR1("(E) IN VERDISINT RDDERI NOT ALLOCATED")
    ENDIF
    CLOPER='RDDERI'
  ELSEIF (IBC==1) THEN
    IF (CDBC=='00') THEN
      IF(ALLOCATED(YDVFE%RDDERI))THEN
        ZOPER = YDVFE%RDDERI 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDDERI NOT ALLOCATED")
      ENDIF
      CLOPER='RDDERI'
    ELSEIF (CDBC=='01') THEN
      IF(ALLOCATED(YDVFE%RDDERBF01))THEN
        ZOPER = YDVFE%RDDERBF01 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDDERBF01 NOT ALLOCATED")
      ENDIF
      CLOPER='RDDERBF01'
    ELSE
      WRITE(NULERR,*) 'VERDISINT: WRONG BOUNDARY CONDITIONS USED.'
      WRITE(NULERR,*) 'VERDISINT: ONLY RDDERI OR RDDERBF01 FOR IBC=',IBC
      CALL ABOR1(' VERDISINT: ABOR1 CALLED')
    ENDIF
  ELSEIF (IBC==2) THEN
    IF (CDBC=='11') THEN
      IF(ALLOCATED(YDVFE%RDDERBF11_IMPL))THEN
        ZOPER = YDVFE%RDDERBF11_IMPL 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDDERBF11_IMPL NOT ALLOCATED")
      ENDIF
      CLOPER='RDDERBF11_IMPL'
    ELSE
      WRITE(NULERR,*) 'VERDISINT: WRONG BOUNDARY CONDITIONS USED.'
      WRITE(NULERR,*) 'VERDISINT: ONLY RDDERBF11_IMPL FOR IBC == 2!!'
      CALL ABOR1(' VERDISINT: ABOR1 CALLED')
    ENDIF
  ELSEIF (IBC==3) THEN
    IF (CDBC=='11') THEN
      IF(ALLOCATED(YDVFE%RDDERBF11))THEN
        ZOPER = YDVFE%RDDERBF11 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDDERBF11 NOT ALLOCATED")
      ENDIF
      CLOPER='RDDERBF11'
    ELSEIF (CDBC=='01') THEN
      IF(ALLOCATED(YDVFE%RDDERBF01))THEN
        ZOPER = YDVFE%RDDERBF01 
      ELSE
        CALL ABOR1("(E) IN VERDISINT RDDERBF01 NOT ALLOCATED")
      ENDIF
      CLOPER='RDDERBF01'
    ELSE
      WRITE(NULERR,*) 'VERDISINT: WRONG BOUNDARY CONDITIONS USED.'
      WRITE(NULERR,*) 'VERDISINT: ONLY RDDERBF01 OR RDDERBF11 FOR &
       & IBC=',IBC
      CALL ABOR1(' VERDISINT: ABOR1 CALLED')
    ENDIF
  ENDIF

ENDIF

!*   3. Apply the required operation
!    --------------------------------

IF (CDOPER=='INGW'.OR.CDOPER=='ITOP'.OR.CDOPER=='IBOT') THEN
  IF(PRESENT(KCHUNK)) THEN
    ICHUNK=KCHUNK
  ELSE
    ICHUNK=1
  ENDIF
  CALL VERINT(KPROMA,KSTART,KPROF,IFLEV_IN,IFLEV_OUT,ZOPER,ZIN(1,IND),POUT,ITYPE,ICHUNK)
ELSEIF (CDOPER=='DEGW'.OR.CDOPER=='FDER'.OR.CDOPER=='HDER'.OR.CDOPER=='DDER') THEN
  CALL VERDER(KPROMA,KSTART,KPROF,IFLEV_IN,IFLEV_OUT,ZOPER,PIN(1,IND),POUT)
ELSE 
  WRITE(NULERR,*) 'VERDISINT: UNKNOWN CDOPER=', CDOPER
  CALL ABOR1(' VERDISINT: ABOR1 CALLED')
ENDIF

IF (ALLOCATED(ZOPER)) DEALLOCATE(ZOPER)
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VERDISINT',1,ZHOOK_HANDLE)
END SUBROUTINE VERDISINT
