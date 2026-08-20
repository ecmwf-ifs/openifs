! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE READFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,       &
     &                  FILENAME, IUNIT, LOUNIT, LCUNIT, LRSTPARAL)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997

!*    PURPOSE.
!     --------
!     READS ARRAY FL

!**   INTERFACE.
!     ----------
!     CALL *READFL*(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,
!                   FILENAME, IUNIT, LOUNIT, LCUNI, LRSTPARAL)
!     *FL*       ARRAY TO BE WRITTEN TO FILE
!     *IJINF*    FIRST LOWER DIMENSION OF FL
!     *IJSUP*    FIRST UPPER DIMENSION OF FL
!     *KINF*     SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*     SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*     THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*     THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE
!     *IUNIT*    FILE UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LOUNIT*  LOGICAL, TRUE IF FREE UNIT HAS TO BE FOUND
!     *LCUNIT* LOGICAL, TRUE IF UNIT HAs TO BE CLOSED 
!     *LRSTPARAL* LOGICAL, TRUE THEN READ IN PARALLEL


!     METHOD.
!     -------
!     READS ARRAY FL FROM FILE (FORTRAN READ)

!     IF PBIO IS USED, THEN PBREAD IS CALLED TO READ (MSUP-MINF+1) 
!     FREQUENCY CONTRIBUTION TO FL.


!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D     ,LLUNSTR
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "iwam_get_unit.intfb.h"
      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IUNIT

      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL

      CHARACTER(LEN=296), INTENT(IN) :: FILENAME

      LOGICAL, INTENT(IN) :: LOUNIT, LCUNIT, LRSTPARAL

      INTEGER(KIND=JWIM) :: LFILE, IJ, J2, J3

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP) :: FL_G

      LOGICAL :: LLEXIST

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READFL',0,ZHOOK_HANDLE)

      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)

      IF (LOUNIT) THEN
        LLEXIST=.FALSE.
        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF (.NOT. LLEXIST) THEN
          WRITE (IU06,*) '*************************************'
          WRITE (IU06,*) '*                                   *'
          WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
          WRITE (IU06,*) '*  IN READFL :                      *'
          WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILENAME
          WRITE (IU06,*) '*                                   *'
          WRITE (IU06,*) '*************************************'
          WRITE (*,*) '*************************************'
          WRITE (*,*) '*                                   *'
          WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
          WRITE (*,*) '*  IN READFL :                      *'
          WRITE (*,*) '*  COULD NOT FIND FILE ',FILENAME
          WRITE (*,*) '*                                   *'
          WRITE (*,*) '*************************************'
          CALL ABORT1
        ENDIF
        IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE), 'r', 'u',0,'READWRITE')
      ENDIF

      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        READ(IUNIT) (((FL_G(IJ,J2,J3),                                  &
     &                  IJ=IJINF,IJSUP),                                &
     &                  J2=KINF,KSUP),                                  &
     &                  J3=MINF,MSUP)

        CALL UNBLKRORD(-1,IJINF,IJSUP,KINF,KSUP,MINF,MSUP,              &
     &                 FL(IJINF:IJSUP,KINF:KSUP,MINF:MSUP),             &
     &               FL_G(IJINF:IJSUP,KINF:KSUP,MINF:MSUP))
#else
      CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
        READ(IUNIT) (((FL(IJ,J2,J3),                                    &
     &                  IJ=IJINF,IJSUP),                                &
     &                  J2=KINF,KSUP),                                  &
     &                  J3=MINF,MSUP)
      ELSE
!       WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!       BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
        READ(IUNIT) (((FL(IJ2NEWIJ(IJ),J2,J3),                          &
     &                  IJ=IJINF,IJSUP),                                &
     &                  J2=KINF,KSUP),                                  &
     &                  J3=MINF,MSUP)

      ENDIF
      IF (LCUNIT) CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('READFL',1,ZHOOK_HANDLE)

      END SUBROUTINE READFL
