! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GPTFS1C(KGPP ,PEPSP1 ,PGP9 ,PGP0 ,PGP1)

!**** *GPTFS1C* - Time filtering and swapping of surface variables

!     Purpose.
!     --------
!             Time filtering and swapping of the surface variables

!***  Interface.
!     ----------
!        *CALL* *GPTFS1C(...)

!        Explicit arguments :
!        --------------------
!                              KGPP   - dimension of GP*
!                              PEPSP1 - CONSTANT FOR THE TIME FILTER
!                              PGP9   - surf. variables at t-dt
!                              PGP0   - surf. variables at t
!                              PGP1   - surf. variables at t+dt

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the one column model

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original : 94-01-14

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM) :: KGPP

!     DUMMY REAL SCALARS
REAL(KIND=JPRB) :: PEPSP1
REAL(KIND=JPRB) :: PGP9(KGPP), PGP0(KGPP), PGP1(KGPP)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JSLEV

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZZPHY

!     ------------------------------------------------------------------

!*       1. PERFORM TIME FILTER AND SWAPPING OF VARIABLES.
!           ----------------------------------------------

ZZPHY=1.0_JPRB-PEPSP1

DO JSLEV=1,KGPP
  PGP9(JSLEV)=PEPSP1*PGP1(JSLEV)+ZZPHY*PGP0(JSLEV)
  PGP0(JSLEV)=PGP1(JSLEV)
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE GPTFS1C



