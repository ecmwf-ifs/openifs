! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LARMES1C(YDGEOMETRY,YDMODEL,KFLEV, PDT, PWRL, PWRL9, PLEV)

#ifdef DOC
!**** *LARMES1C - semi-LAgrangian scheme:
!                 Research of the MEdium point.

!     Purpose.
!     --------
!      The computation of the location of the origine point O of
!      the lagrangian trajectory is performed by an iterative
!      method described by Robert. 

!**   Interface.
!     ----------
!        *CALL* *LARMES1C

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KFLEV   - vertical dimension.
!          PWRL    - W-component of the wind. (eta dot)
!          PWRL9   - W-component of the wind from previous time-step. 
!          PDT     - Time step
!        OUTPUT:
!          PLEV    - vertical coordinate of the interpolation point.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by LAPINE1C.

!     Reference.
!     ----------

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original    19-05-1994
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!        F. Vana     9-Jul-2014  rewritten to deliver the origin point only

!     ------------------------------------------------------------------
#endif

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM) :: KFLEV

!     DUMMY REAL SCALARS
REAL(KIND=JPRB) :: PDT

!     INPUT:
REAL(KIND=JPRB) :: PWRL(KFLEV)
REAL(KIND=JPRB) :: PWRL9(KFLEV)

!     OUTPUT:
REAL(KIND=JPRB) :: PLEV(KFLEV)

!     LOCAL SCALARS
INTEGER(KIND=JPIM) :: ISTESB, ISTEST, JITER, JLEV, KWISA
REAL(KIND=JPRB) :: ZLEVB, ZLEVO, ZLEVT

!     LOCAL ARRAYS:
INTEGER(KIND=JPIM) :: ILEV(KFLEV)
REAL(KIND=JPRB) :: ZWF(KFLEV),  ZWRL(KFLEV)
REAL(KIND=JPRB) :: ZDUM(KFLEV)
REAL(KIND=JPRB) :: ZDVER(KFLEV), ZVINTW(KFLEV,2:4)


!     ------------------------------------------------------------------
#include "lascaw1c.intfb.h"
#include "laitli1c.intfb.h"
!     ------------------------------------------------------------------
ASSOCIATE(YDDYN=>YDMODEL%YRML_DYN%YRDYN)

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

KWISA=1  ! Linear interpolation 

! Top and bottom limits for departure points 
ZLEVT=YDGEOMETRY%YRVETA%VETAF(0)
ZLEVB=YDGEOMETRY%YRVETA%VETAF(KFLEV+1)

! SETTLST scheme
ZWRL(1:KFLEV)=2._JPRB*PWRL(1:KFLEV)-PWRL9(1:KFLEV)
! Store W9 for next timestep
PWRL9(1:KFLEV)=PWRL(1:KFLEV)

! Initial value for medium point wind
ZWF(1:KFLEV)=ZWRL(1:KFLEV)


!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------

DO JITER=1,YDDYN%NITMP

!*       2.1   DETERMINATION OF THE DEPARTURE POINT 

  ISTEST=0
  ISTESB=0

  DO JLEV=1,KFLEV

    ZLEVO=YDGEOMETRY%YRVETA%VETAF(JLEV)-PDT*0.5_JPRB*(ZWF(JLEV)+PWRL(JLEV))
    !ZLEVO=MAX(YRVETA%VETAF(1),MIN(YRVETA%VETAF(KFLEV),ZLEVO))
    ISTEST=ISTEST-MIN(0,MAX(-1,NINT(ZLEVO-ZLEVT-0.5_JPRB)))
    ISTESB=ISTESB-MIN(0,MAX(-1,NINT(ZLEVB-ZLEVO-0.5_JPRB)))
    ZLEVO=MIN(ZLEVB,MAX(ZLEVT,ZLEVO))
    PLEV(JLEV)=ZLEVO

  ENDDO

  IF(ISTEST /= 0) THEN
    WRITE(NULOUT,'(A,I6,A)') ' SMILAG TRAJECTORY OUT OF ATM ',&
     &ISTEST,' TIMES.'
  ENDIF
  IF(ISTESB /= 0) THEN
    WRITE(NULOUT,'(A,I6,A)') ' SMILAG TRAJECTORY UNDERGROUND ',&
     &ISTESB,' TIMES.'
  ENDIF

!*       2.2   DETERMINATION OF THE WIND AT THE DEPARTURE POINT.

  IF(JITER /= YDDYN%NITMP) THEN

    CALL LASCAW1C(YDGEOMETRY,YDMODEL%YRML_DYN%YRSLINT,KFLEV,PLEV , ZVINTW , ILEV , ZDVER , KWISA)
    CALL LAITLI1C(YDGEOMETRY%YRDIMV,KFLEV,ZDVER,ILEV,ZWRL,ZWF)

  ENDIF

ENDDO

END ASSOCIATE
END SUBROUTINE LARMES1C
