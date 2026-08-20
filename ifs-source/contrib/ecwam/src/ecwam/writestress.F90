! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WRITESTRESS(IJINF, IJSUP, NREAL, RFIELD, FILENAME, LRSTPARAL) 

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      APRIL 1997

!*    PURPOSE.
!     --------
!     WRITES THE RESTART WIND AND STRESS FIELDS TO FILE.

!**   INTERFACE.
!     ----------

!     CALL *WRITESTRESS*(IJINF, IJSUP, NREAL, RFIELD,
!                        FILENAME,LRSTPARAL)
!     *IJINF*    FIRST LOWER DIMENSION OF RFIELD 
!     *IJSUP*    FIRST UPPER DIMENSION OF RFIELD 
!     *NREAL*    NUMBER OF FIELDS IN RFIELD
!     *RFIELD*   REAL RESTART FIELDS
!     *FILENAME* FILENAME (INCLUDING PATH) OF TARGET FILE
!     *LRSTPARAL* LOGICAL, TRUE THEN WRITE IN PARALLEL

!     METHOD.
!     -------
!     WRITES ARRAYS TO FILE (FORTRAN WRITES).

!     THE OUTPUT IS WRITTEN AS UNFORMATTED BINARY TO FILENAME
!     IT IS IN A FORM THAT IS INDEPENDENT OF THE MODEL DECOMPOSITION.
!     IF LRSTPARAL IS FALSE


!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D, LLUNSTR
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOWABORT , ONLY : WAM_ABORT
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "iwam_get_unit.intfb.h"
      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, NREAL
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP,NREAL),INTENT(INOUT) :: RFIELD
      CHARACTER(LEN=296), INTENT(IN) :: FILENAME
      LOGICAL, INTENT(IN) :: LRSTPARAL

      INTEGER(KIND=JWIM) :: IFLD, IJ, IG
      INTEGER(KIND=JWIM) :: LFILE, IUNIT

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP) :: RFIELD_G

      CHARACTER(LEN=1) :: MODE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WRITESTRESS',0,ZHOOK_HANDLE)

      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
      IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'w', 'u', 0,'READWRITE')

      WRITE(IUNIT) CDTPRO, CDATEWO, CDAWIFL, CDATEFL

      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        DO IFLD=1,NREAL
          CALL UNBLKRORD(1,IJINF,IJSUP,1,1,1,1,RFIELD(1,IFLD),RFIELD_G(1))
          WRITE(IUNIT)(RFIELD_G(IJ),IJ=IJINF,IJSUP)
        ENDDO
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
        DO IFLD=1,NREAL
          WRITE(IUNIT)(RFIELD(IJ,IFLD),IJ=IJINF,IJSUP)
        ENDDO
      ELSE
        DO IFLD=1,NREAL
          WRITE(IUNIT) (RFIELD(IJ2NEWIJ(IJ),IFLD),IJ=IJINF,IJSUP)
        ENDDO
      ENDIF

      CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('WRITESTRESS',1,ZHOOK_HANDLE)

      END SUBROUTINE WRITESTRESS
