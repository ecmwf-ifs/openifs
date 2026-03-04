! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LACDYN1C(YDVETA,KLEV,PVVEL,PETADOTDPDETA,PDELP,PWRL0)

#ifdef DOC
!**** *LACDYN1C*  Semi-Lagrangian scheme.

!     Purpose.
!     --------
!          This subroutine computes the vertical velocity
!          in eta coordinates, from dp/dt prescribed by the user.
!          It assumes implicitly that the tendency of ps is zero and 
!          that V.GRAD(p) is negligible. 

!**   Interface.
!     ----------
!        *CALL* *LACDYN1C

!        Explicit arguments :
!        --------------------

!          PDELP   - pressure depth of layers at t.           (input)
!          PVVEL   - full level omega                         (input)
!          PETADOTDPDETA - half level omega                   (input) 
!          PWRL0   - eta dot                                  (output)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           Called by CPG1C.

!     Reference.
!     ----------

!     Author.
!     -------
!        Joao Teixeira   *ECMWF*

!     Modifications.
!     --------------
!        Original    94-05-18
!        M. Ko"hler  6-6-2006  Single Column Model integration within IFS 
!        F. Vana     9-Jul-2014 phasing with CY40R1
!     ------------------------------------------------------------------
#endif

USE YOMVERT  , ONLY : TVAB, TVETA, TVFE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMLOG1C , ONLY : LETADOT

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
TYPE(TVETA), INTENT(INOUT) :: YDVETA
INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
REAL(KIND=JPRB), INTENT(IN)    :: PDELP(KLEV)
REAL(KIND=JPRB), INTENT(IN)    :: PETADOTDPDETA(0:KLEV)
REAL(KIND=JPRB), INTENT(IN)    :: PVVEL(KLEV)
REAL(KIND=JPRB), INTENT(OUT)   :: PWRL0(KLEV)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JLEV


!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF D(ETA)/DT.
!              -------------------------

DO JLEV=1,KLEV
  IF (LETADOT) THEN
    PWRL0(JLEV)=0.5_JPRB*( PETADOTDPDETA(JLEV)+PETADOTDPDETA(JLEV-1) )
  ELSE
    PWRL0(JLEV)=PVVEL(JLEV)
  ENDIF
  PWRL0(JLEV)= PWRL0(JLEV)*((YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))/PDELP(JLEV))
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE LACDYN1C
