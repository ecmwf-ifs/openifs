! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE EC_PHYS_LSLPHY(&
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDGEOMETRY,YDMODEL, &
 & KST,KEND,LDLFSTEP,PDTPHY,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PULP9,PVLP9,PTLP9,PGFLLP9,PSPPTCLEARLP9,&
 & PSAVTEND,PGFLSLP)

!**** *EC_PHYS_LSLPHY* - Grid point calculations: get EC physical tendencies.

!     Purpose.
!     --------
!           Grid point calculations: get EC physical tendencies
!           computed at the previous time step (lagged physics)
!           to simulate a non lagged physics.

!**   Interface.
!     ----------
!        *CALL* *EC_PHYS_LSLPHY(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        LDLSLPHY  : .T.: split physics activated.
!        LDLFSTEP  : .T.: first time-step?
!        PDTPHY    : timestep used in the physics.

!     INPUT/OUTPUT:
!     -------------
!        PULP9     : quantity in U-wind eqn for split physics to be
!                    interpolated at the origin point 
!        PVLP9     : quantity in V-wind eqn for split physics to be
!                    interpolated at the origin point 
!        PTLP9     : quantity in T eqn for split physics to be
!                    interpolated at the origin point 
!        PGFLLP9   : unified_treatment grid-point fields to be
!                    interpolated at the origin point for split physics

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!       K. YESSAD after old version of CPG. 

!     Modifications.
!     --------------
!   Original : 01-08-13
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   K. Yessad (July 2014): Move some variables.
!   F. Vana   11-Sep-2020 Cleaning
!   F. Vana and S. Lang January 2023: SPPT updated to new seq. physics order
!     ------------------------------------------------------------------

USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE TYPE_MODEL             , ONLY : MODEL
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDLFSTEP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTPHY 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PULP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVLP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTLP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLLP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN,&
 &                                          YDMODEL%YRML_GCONF%YGFL%NDIMSLP) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPPTCLEARLP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSAVTEND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRSLPHY%NVTEND) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFLSLP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMSLP) 

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gpaddslphy.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EC_PHYS_LSLPHY',0,ZHOOK_HANDLE)
ASSOCIATE(YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
  &   YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2)
ASSOCIATE(TDT=>YDRIP%TDT,TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

IF(PDTPHY /= TDT) THEN
  TSPHY = MAX(PDTPHY,1.0_JPRB)
ENDIF

CALL GPADDSLPHY(YDGEOMETRY,YDMODEL,YGFL,LDLFSTEP,KST,KEND,TSPHY,&
 & PSAVTEND,PGFLSLP,&
 & PULP9,PVLP9,PTLP9,PGFLLP9,PSPPTCLEARLP9)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EC_PHYS_LSLPHY',1,ZHOOK_HANDLE)
END SUBROUTINE EC_PHYS_LSLPHY
