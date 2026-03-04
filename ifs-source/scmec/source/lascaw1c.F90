! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LASCAW1C(YDGEOMETRY, YDSLINT, KFLEV, PLEV, PVINTW, KLEV, PDVER, KWIS )

#ifdef DOC
!**** *LASCAW  -  semi-LAgrangian scheme:
!                 Storage of Coordinates And Weights.

!     Purpose.
!     --------
!       Computation and storage of coordinates and weights.

!**   Interface.
!     ----------
!        *CALL* *LASCAW1C

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KFLEV   - vertical dimension.
!          PLEV    - vertical coordinate of the interpolation point.
!                    coefficients
!        OUTPUT:
!          KLEV    - lower level of the vertical interpolation
!                    grid needed for vertical interpolations (LAIVRE).
!          PDVER   - weights (distances) for vertical linear interpolation
!                    on a same vertical.
!          PVINTW  - Vertical cubic interpolation weights

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        No external.
!        Called by LARCIN1C.

!     Reference.
!     ----------

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original    20-05-1994
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 

!     ------------------------------------------------------------------
#endif

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMSLINT,     ONLY : TSLINT
USE PARKIND1,     ONLY : JPIM, JPRB

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TSLINT),   INTENT(IN) :: YDSLINT
INTEGER(KIND=JPIM) :: KFLEV
INTEGER(KIND=JPIM) :: KWIS
REAL(KIND=JPRB) :: PLEV (KFLEV)

!     OUTPUT:
INTEGER(KIND=JPIM) :: KLEV(KFLEV)
REAL(KIND=JPRB) :: PDVER(KFLEV)
REAL(KIND=JPRB) :: PVINTW(KFLEV,2:4)
REAL(KIND=JPRB) :: PXF(KFLEV)
REAL(KIND=JPRB) :: PXSL(KFLEV)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ILEV, ILEVV, JLEV, IFLVM2

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZD1, ZD2, ZD3, ZD4, ZFAC


!     ----------------------------------------------------------------
!
IFLVM2=KFLEV-2

!        1.  COORDINATES AND WEIGHTS FOR LINEAR INTERPOLATIONS.
!        ------------------------------------------------------    

IF (KWIS == 1) THEN

  ZFAC=YDSLINT%YRVSLETA%VRLEVX/(YDGEOMETRY%YRVETA%VETAF(KFLEV+1)-YDGEOMETRY%YRVETA%VETAF(0))

  DO JLEV=1,KFLEV

    ILEV  = YDSLINT%YRVSLETA%NVAUTF(INT(PLEV(JLEV)*ZFAC))-1
    IF(ILEV < IFLVM2.AND.&
       & (PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEV+2)) > 0.0_JPRB) ILEV=ILEV+1  
    KLEV(JLEV)=ILEV
    PDVER(JLEV)=(PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEV+1))/(YDGEOMETRY%YRVETA%VETAF(ILEV+2)-YDGEOMETRY%YRVETA%VETAF(ILEV+1))

  ENDDO

ENDIF


!        2.  COORDINATES AND WEIGHTS FOR VERTICAL CUBIC INT.
!        ---------------------------------------------------

IF (KWIS == 3) THEN

  ZFAC=YDSLINT%YRVSLETA%VRLEVX/(YDGEOMETRY%YRVETA%VETAF(KFLEV+1)-YDGEOMETRY%YRVETA%VETAF(0))

  DO JLEV=1,KFLEV

    ILEVV=YDSLINT%YRVSLETA%NVAUTF(INT(PLEV(JLEV)*ZFAC))-1
    IF(ILEVV < IFLVM2.AND.&
       & (PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEVV+2)) > 0.0_JPRB) ILEVV=ILEVV+1  
    KLEV(JLEV)=ILEVV
    PDVER(JLEV)=(PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEVV+1))/(YDGEOMETRY%YRVETA%VETAF(ILEVV+2)-YDGEOMETRY%YRVETA%VETAF(ILEVV+1))

    ILEVV=KLEV(JLEV)
    IF(ILEVV >= 1.AND.ILEVV <= KFLEV-3) THEN
      ZD1=PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEVV  )
      ZD2=PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEVV+1)
      ZD3=PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEVV+2)
      ZD4=PLEV(JLEV)-YDGEOMETRY%YRVETA%VETAF(ILEVV+3)
      PVINTW(JLEV,2)=ZD1    *ZD3*ZD4*YDSLINT%YRVSLETA%VCUICO(2,ILEVV)
      PVINTW(JLEV,3)=ZD1*ZD2    *ZD4*YDSLINT%YRVSLETA%VCUICO(3,ILEVV)
      PVINTW(JLEV,4)=ZD1*ZD2*ZD3    *YDSLINT%YRVSLETA%VCUICO(4,ILEVV)
    ELSE
      PVINTW(JLEV,2)=1.0_JPRB-PDVER(JLEV)
      PVINTW(JLEV,3)=PDVER(JLEV)
      PVINTW(JLEV,4)=0.0_JPRB
    ENDIF

  ENDDO

ENDIF


END SUBROUTINE LASCAW1C
