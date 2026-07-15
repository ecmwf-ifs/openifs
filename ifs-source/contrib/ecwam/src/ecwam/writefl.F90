! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,      &
     &                   FILENAME, IUNIT, LOUNIT, LRSTPARAL)

! ---------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     WRITES ARRAY FL TO FILE.

!**   INTERFACE.
!     ----------
!     CALL *WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,
!    &              FILENAME,IUNIT,LOUNIT,LRSTPARAL)
!     *FL*        ARRAY TO BE WRITTEN TO FILE
!     *IJINF*     FIRST LOWER DIMENSION OF FL
!     *IJSUP*     FIRST UPPER DIMENSION OF FL
!     *KINF*      SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*      SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*      THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*      THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME*  FILENAME (INCLUDING PATH) OF TARGET FILE
!     *IUNIT*     PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LOUNIT*    LOGICAL, TRUE IF FREE UNIT HAS TO BE FOUND
!     *LRSTPARAL* LOGICAL, TRUE THEN WRITE IN PARALLEL

!     METHOD.
!     -------
!     WRITES ARRAY FL TO FILE (A FORTRAN WRITE)

!     FL IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
!     TO UNIT IUNIT.
!     IT IS IN A FORM THAT IS INDEPENDENT OF THE MODEL DECOMPOSITION.



!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D, LLUNSTR
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWABORT , ONLY : WAM_ABORT
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "iwam_get_unit.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IUNIT

      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL

      CHARACTER(LEN=296), INTENT(IN) :: FILENAME

      LOGICAL, INTENT(IN) :: LOUNIT, LRSTPARAL

      INTEGER(KIND=JWIM) :: IJ, K, M, J2, J3
      INTEGER(KIND=JWIM) :: KRET, KOUNT, LFILE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP) :: FL_G
      REAL(KIND=JWRB), ALLOCATABLE :: RFL(:,:,:)

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WRITEFL',0,ZHOOK_HANDLE)

      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
      IF (LOUNIT) THEN
        IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'w', 'u', 0, 'READWRITE')
      ELSE
        IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'a', 'u', 0, 'READWRITE')
      ENDIF

      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        FL_G(0,:,:)=0.0_JWRB
        CALL UNBLKRORD(1,IJINF,IJSUP,KINF,KSUP,MINF,MSUP,               &
     &                 FL(IJINF:IJSUP,KINF:KSUP,MINF:MSUP),             &
     &               FL_G(IJINF:IJSUP,KINF:KSUP,MINF:MSUP))

        WRITE(IUNIT) (((FL_G(IJ,J2,J3),IJ=IJINF,IJSUP),J2=KINF,KSUP),J3=MINF,MSUP)
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
        WRITE(IUNIT) (((FL(IJ,J2,J3),IJ=IJINF,IJSUP),J2=KINF,KSUP),J3=MINF,MSUP)
      ELSE
!     WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!     BUT THE SINGLE BINARY INPUT FILES SHOULD BE IN THE OLD MAPPING
        WRITE(IUNIT) (((FL(IJ2NEWIJ(IJ),J2,J3),IJ=IJINF,IJSUP),J2=KINF,KSUP),J3=MINF,MSUP)
      ENDIF

      CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('WRITEFL',1,ZHOOK_HANDLE)

      END SUBROUTINE WRITEFL
