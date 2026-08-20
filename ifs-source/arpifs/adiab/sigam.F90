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

SUBROUTINE SIGAM(YDGEOMETRY,YDDYN,KLEV,KLON,PD,PT,PSP,KNLON,KFLEVG)

!**** *SIGAM* - Solve hydrostatic operator in semi-implicit

!     Purpose.
!     --------
!           Operator gamma to compute p.

!**   Interface.
!     ----------
!        *CALL* *SIGAM(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PD    : DIVERGENCE       (output)
!        PT    : TEMPERATURE      (input)
!        PSP   : SURFACE PRESSURE (input)
!        KNLON : NUMBER OF VERTICAL COLUMNS TREATED
!        KFLEVG: NUMBER OF ELEMENTS IN A VERTICAL COLUMN

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      Modified : 09-Oct-2007 by K. YESSAD: possibility to have a specific
!                 value of LVERTFE in the SI linear model.
!      F. Vana + NEC 28-Apr-2009: OpenMP
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      G. Mozdzynski Oct 2012: OpenMP optimization
!      K. Yessad (Dec 2016): Prune obsolete options.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      R.Brozkova + NEC 03-Mar-2021: Optimization for vector (NEC)
!      R. El Khatib 28-Feb-2023 Bugfixes for open-mp
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : RD
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(1+(KFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(1+(KFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KNLON) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPHI(KNLON,0:KFLEVG+1)
REAL(KIND=JPRB) :: ZOUT(KNLON,0:KFLEVG)
REAL(KIND=JPRB) :: ZSPHIX(0:KFLEVG,KNLON)
INTEGER(KIND=JPIM) :: IDT, JLEV, JLON
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGAM',0,ZHOOK_HANDLE)
ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SILNPR=>YDDYN%SILNPR, SIRPRG=>YDDYN%SIRPRG)
!     ------------------------------------------------------------------

!*       1.    SUM GEOPOTENTIAL, COMPUTES P AND PUT IT IN PD.
!              ----------------------------------------------

IF(YDCVER%LVERTFE) THEN

!$OMP PARALLEL PRIVATE(JLEV,JLON,IDT,ZDETAH)
!$OMP DO SCHEDULE(STATIC) 
  DO JLEV=1,KFLEVG
    ZDETAH=YDVETA%VFE_RDETAH(JLEV)
    DO JLON=1,KNLON
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZSPHI(JLON,JLEV)=-RD*PT(IDT)*SILNPR(JLEV)*ZDETAH
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ZSPHI(1:KNLON,0)=0.0_JPRB
  ZSPHI(1:KNLON,KFLEVG+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'IBOT','11',KNLON,1,KNLON,KFLEVG,ZSPHI,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)

!$OMP PARALLEL PRIVATE(JLEV,JLON,IDT)
!$OMP DO SCHEDULE(STATIC) 
  DO JLEV=1,KFLEVG
    DO JLON=1,KNLON
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      PD(IDT)=ZOUT(JLON,JLEV-1)+PSP(JLON)*SIRPRG
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

ELSE

#ifndef __NEC__
!$OMP PARALLEL PRIVATE(JLON,JLEV,IDT)
!$OMP DO SCHEDULE(STATIC)
  DO JLON=1,KNLON
    ZSPHIX(KFLEVG,JLON)=0.0_JPRB
    DO JLEV=KFLEVG,1,-1
#else
!NEC_PL: move initialization out of OMP loop  
  DO JLON=1,KNLON
    ZSPHIX(KFLEVG,JLON)=0.0_JPRB
  ENDDO
! Warning : loop interchange to eliminate dependency on ZSPHIX is not compatible
! with open-mp parallelization of the outer loop. REK.
  DO JLEV=KFLEVG,1,-1
    DO JLON=1,KNLON
#endif
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZSPHIX(JLEV-1,JLON)=ZSPHIX(JLEV,JLON)+RD*PT(IDT)*SILNPR(JLEV)
      PD(IDT)=ZSPHIX(JLEV,JLON)+SIALPH(JLEV)*RD*PT(IDT)+PSP(JLON)*SIRPRG
    ENDDO
  ENDDO
#ifndef __NEC__
!$OMP END DO
!$OMP END PARALLEL
#endif

ENDIF

!      -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SIGAM',1,ZHOOK_HANDLE)
END SUBROUTINE SIGAM
