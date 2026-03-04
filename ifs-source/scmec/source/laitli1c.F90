! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LAITLI1C(YDDIMV, KFLEV, PDVER, KLEV, PXSL, PXF )

#ifdef DOC
!**** *LAITLI1C  -  semi-LAgrangian scheme:
!                   Linear interpolation for one variable.

!     Purpose.
!     --------
!       Performs linear interpolation for one variable.

!**   Interface.
!     ----------
!        *CALL* *LAITLI1C

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KFLEV   - vertical dimension.
!          PDVER   - weights (distances) for vertical linear interpolation
!                    on a same vertical.
!          PXSL    - semi-lagrangian variable.

!        OUTPUT:
!          PXF     - interpolated variable.

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
!        Original  20-05-1994
!        M. Ko"hler  6-6-2006  Single Column Model integration within IFS 

!     ------------------------------------------------------------------
#endif

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM     ,JPRB

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
INTEGER(KIND=JPIM) :: KFLEV
INTEGER(KIND=JPIM) :: KLEV(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PDVER(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PXSL(YDDIMV%NFLEVG)

!     OUTPUT:
REAL(KIND=JPRB) :: PXF(YDDIMV%NFLEVG)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ILEV, JLEV

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZDVER, ZINF, ZSUP


!     ------------------------------------------------------------------

!*       1.    INTERPOLATIONS.
!              ---------------

DO JLEV=1,KFLEV

!     Computation of coordinates and distances.

  ILEV=KLEV(JLEV)

  ZDVER = PDVER(JLEV)

!     Interpolation.

  ZSUP  = PXSL (ILEV+1)
  ZINF  = PXSL (ILEV+2)

  PXF(JLEV)= ZSUP + ZDVER*(ZINF-ZSUP)

ENDDO


END SUBROUTINE LAITLI1C
