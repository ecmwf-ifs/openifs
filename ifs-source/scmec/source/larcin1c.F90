! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LARCIN1C(YDGEOMETRY,YDMODEL,KFLEV,PLEV&
           &,PQSL,PTSL,PUSL,PVSL,PASL,PLSL,PISL&
           &,PQF,PTF,PUF,PVF,PAF,PLF,PIF,KWISA) !,PWRL0,PWF0)

#ifdef DOC
!**** *LARCIN1C  -  semi-LAgrangian scheme:
!                 Research of the Coordinates (of the medium 
!                 or origin point) and INterpolations.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LARCIN1C

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KFLEV   - vertical dimension.
!          KWISA   - switches.
!          PLEV    - vertical coordinate of the interpolation point.
!          PVETAF  - Values of ETA at full levels
!          PVCUICO - Contains denominators of the weights of cubic
!                    vertical interpolations.
!          PWRL0   - W-component of the wind.
!          PUSL    - Quantity to interpolate at the origin or medium point
!                    in the U-component of the wind equation.
!          PVSL    - Quantity to interpolate at the origin or medium point
!                    in the V-component of the wind equation.
!          PTSL    - Quantity to interpolate at the origin or medium point
!                    in the temperature equation.
!          PQSL    - Quantity to interpolate at the origin or medium point
!                    in the humidity equation.

!        OUTPUT:
!          PWF0    - Interpolated W-wind at the medium point.
!          PUF     - Interpolated quantity at the origin or medium point
!                    in the U-component of the wind equation.
!          PVF     - Interpolated quantity at the origin or medium point
!                    in the V-component of the wind equation.
!          PTF     - Interpolated quantity at the origin or medium point
!                    in the temperature equation.
!          PQF     - Interpolated quantity at the origin or medium point
!                    in the humidity equation.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!        Called by LARMES1C and LAINOR1C.
!        Calls  LAITRI1C, LAITLI1C AND LASCAW1C.

!     Reference.
!     ----------

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original    19-05-1994
!        J.Teixeira   Jun.-95:   quasi-monotone interpolations.
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 

!     ------------------------------------------------------------------
#endif

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1 , ONLY : JPIM, JPRB

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM) :: KFLEV
INTEGER(KIND=JPIM) :: KWISA

!     INPUT:
REAL(KIND=JPRB) :: PLEV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PVETAF(0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: PVCUICO(4,0:YDGEOMETRY%YRDIMV%NFLEVG-1)
!REAL(KIND=JPRB) :: PWRL0(YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PUSL (YDGEOMETRY%YRDIMV%NFLEVG) , PVSL (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTSL (YDGEOMETRY%YRDIMV%NFLEVG) , PQSL (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PASL (YDGEOMETRY%YRDIMV%NFLEVG) , PLSL (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PISL (YDGEOMETRY%YRDIMV%NFLEVG)
!     OUTPUT:
!REAL(KIND=JPRB) :: PWF0(YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PUF (YDGEOMETRY%YRDIMV%NFLEVG) , PVF (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTF (YDGEOMETRY%YRDIMV%NFLEVG) , PQF (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PAF (YDGEOMETRY%YRDIMV%NFLEVG) , PLF (YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: PIF (YDGEOMETRY%YRDIMV%NFLEVG)
!     LOCAL:
INTEGER(KIND=JPIM) :: ILEV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDVER(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVINTW(YDGEOMETRY%YRDIMV%NFLEVG,2:4)
LOGICAL :: LLQMQ, LLQMCLD

!     ------------------------------------------------------------------
#include "laitri1c.intfb.h"
#include "laitqm1c.intfb.h"
#include "laitli1c.intfb.h"
#include "lascaw1c.intfb.h"
!     ------------------------------------------------------------------
ASSOCIATE(YDDYN=>YDMODEL%YRML_DYN%YRDYN)

LLQMQ=.TRUE.    ! QM setting for moisture
LLQMCLD=.FALSE. ! QM setting for other cloud variables

!*       1.    COMPUTATION OF COORDINATES AND WEIGHTS 
!              OF MEDIUM AND ORIGIN POINTS. 
!              --------------------------------------

CALL LASCAW1C(YDGEOMETRY,YDMODEL%YRML_DYN%YRSLINT,KFLEV,PLEV , ZVINTW , ILEV , ZDVER , KWISA)


!     ------------------------------------------------------------------

!*       2.    MEDIUM POINT INTERPOLATIONS OF WIND
!              FOR COMPUTATION OF MEDIUM POINT.
!              -----------------------------------


IF (KWISA == 1) THEN

!*    LINEAR INTERPOLATIONS.

  CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,PUSL,PUF)
  CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,PVSL,PVF)
  CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,PTSL,PTF)
  CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,PQSL,PQF)
  CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,PASL,PAF)
  CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,PLSL,PLF)
  CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,PISL,PIF)

ENDIF


!     ------------------------------------------------------------------

!*       4.    ORIGIN POINT INTERPOLATIONS.
!              ----------------------------

IF (KWISA == 3) THEN

!     - for the U-wind equation quantity.
!     - for the V-wind equation quantity.
!     - for the temperature equation quantity.
!     - for the humidity equation quantity.
!     - for the cloud fraction equation quantity.
!     - for the liquid water equation quantity.
!     - for the ice equation quantity.

  IF (YDDYN%LQMW) THEN
    CALL LAITQM1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PUSL,PUF )
    CALL LAITQM1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PVSL,PVF )
  ELSE
    CALL LAITRI1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PUSL,PUF )
    CALL LAITRI1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PVSL,PVF )
  ENDIF

  
  IF (YDDYN%LQMT) THEN
    CALL LAITQM1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PTSL,PTF )
  ELSE
    CALL LAITRI1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PTSL,PTF )
  ENDIF

  IF (LLQMQ) THEN
    CALL LAITQM1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PQSL,PQF )
  ELSE
    CALL LAITRI1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PQSL,PQF )
  ENDIF

  IF (LLQMCLD) THEN
    CALL LAITQM1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PASL,PAF )
    CALL LAITQM1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PLSL,PLF )
    CALL LAITQM1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PISL,PIF )
  ELSE
    CALL LAITRI1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PASL,PAF )
    CALL LAITRI1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PLSL,PLF )
    CALL LAITRI1C(YDGEOMETRY%YRDIMV,KFLEV , ZVINTW ,ILEV ,PISL,PIF )
  ENDIF

ENDIF


END ASSOCIATE
END SUBROUTINE LARCIN1C
