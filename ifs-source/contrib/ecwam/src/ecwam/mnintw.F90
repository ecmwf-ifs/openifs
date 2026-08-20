! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE  MNINTW (                                              &
!                     DIMENSIONS
     &                  NSPECW, MPARTSW,                                &
!                     INPUT
     &                  NPARTW, NIPARTINFW,                             &
!                     OUTPUT
     &                  NINTW)

!---------------------------------------------------------------------

!**** *MNINTW*  READS GRID OUTPUT FROM PREPROC.

!     R.BROKOPF MAY 94.

!*    PURPOSE.
!     --------
!       COMPUTES INDEX OF WIND-SEA PARTITIONING FOR EACH SPECTRUM.

!**   INTERFACE.
!     ----------
!       *CALL* *MNINTW(
!                     DIMENSIONS
!     &                  NSPECW , NSPECW , MPARTSW ,
!                     INPUT
!     &                  NPARTW, NIPARTINFW,
!                     OUTPUT
!     &                  NINTW)*

!     METHOD.
!     -------

!       NONE

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: NSPECW , MPARTSW

      INTEGER(KIND=JWIM) :: NPARTW(NSPECW) , NINTW(NSPECW),             &
     &                      NIPARTINFW(NSPECW,MPARTSW)

      INTEGER(KIND=JWIM) :: IPART, ISPEC

! ----------------------------------------------------------------------

      DO ISPEC = 1, NSPECW
         NINTW(ISPEC) = 0
      END DO

      DO ISPEC = 1, NSPECW
        DO IPART = 1,NPARTW(ISPEC)
          IF(NIPARTINFW(ISPEC,IPART).NE.0)                              &
     &       NINTW(ISPEC)=IPART
        ENDDO
      ENDDO

      END SUBROUTINE MNINTW
