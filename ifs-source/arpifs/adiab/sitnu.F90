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

SUBROUTINE SITNU(YDGEOMETRY,YDDYN,KLEV,KLON,PD,PT,PSP,KNLON)

!**** *SITNU*   - Continuity equation for semi-implicit.

!     Purpose.
!     --------
!           Evaluate operators Tau and Nu in semi-implicit.

!**   Interface.
!     ----------
!        *CALL* *SITNU(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PD    : DIVERGENCE
!        PT    : TEMPERATURE
!        PSP   : SURFACE PRESSURE
!        KNLON : NUMBER OF VERTICAL COLUMNS TREATED

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None.
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
!                 value of LVERTFE in the SI NH linear model.
!      F. Vana + NEC 28-Apr-2009: OpenMP + optimization
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      G. Mozdzynski Oct 2012: OpenMP optimization
!      K. Yessad (Dec 2016): Prune obsolete options.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : RKAPPA
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP(KNLON) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIV(KNLON,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZOUT(KNLON,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSDIVX(0:YDGEOMETRY%YRDIMV%NFLEVG,KNLON)
INTEGER(KIND=JPIM) :: IDT, JLEV, JLON
REAL(KIND=JPRB) :: ZREC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SITNU',0,ZHOOK_HANDLE)
ASSOCIATE( YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,YDCVER=>YDGEOMETRY%YRCVER) 
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & SIALPH=>YDDYN%SIALPH, &
 & SIDELP=>YDDYN%SIDELP, SILNPR=>YDDYN%SILNPR, SIRDEL=>YDDYN%SIRDEL, &
 & SIRPRN=>YDDYN%SIRPRN, SITLAF=>YDDYN%SITLAF, SITR=>YDDYN%SITR, &
 & YDDIMV=>YDGEOMETRY%YRDIMV)
!     ------------------------------------------------------------------

!*       1.    SUM DIVERGENCE AND COMPUTES TEMPERATURE.
!              ----------------------------------------

IF(YDCVER%LVERTFE) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLON=1,KNLON
    DO JLEV=1,NFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZSDIV(JLON,JLEV)=PD(IDT)*SIDELP(JLEV)*YDVETA%VFE_RDETAH(JLEV)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  IF (KNLON>=1) THEN
!DEC$ IVDEP
    ZSDIV(1:KNLON,0)=0.0_JPRB
    ZSDIV(1:KNLON,NFLEVG+1)=0.0_JPRB
    CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',KNLON,1,KNLON,NFLEVG,ZSDIV,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)
  ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT,ZREC)
  DO JLON=1,KNLON
    DO JLEV=1,NFLEVG
      ZREC=1.0_JPRB/SITLAF(JLEV)
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      PT(IDT)=RKAPPA*SITR*ZOUT(JLON,JLEV-1)*ZREC
    ENDDO
    PSP(JLON)=ZOUT(JLON,NFLEVG)*SIRPRN
  ENDDO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLON=1,KNLON
    ZSDIVX(0,JLON)=0.0_JPRB
    DO JLEV=1,NFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZSDIVX(JLEV,JLON)=ZSDIVX(JLEV-1,JLON)+PD(IDT)*SIDELP(JLEV)
      PT(IDT)=RKAPPA*SITR*(SIRDEL(JLEV)*SILNPR(JLEV)*ZSDIVX(JLEV-1,&
       & JLON)&
       & +SIALPH(JLEV)*PD(IDT))  
    ENDDO
    PSP(JLON)=ZSDIVX(NFLEVG,JLON)*SIRPRN
  ENDDO
!$OMP END PARALLEL DO
ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SITNU',1,ZHOOK_HANDLE)
END SUBROUTINE SITNU
