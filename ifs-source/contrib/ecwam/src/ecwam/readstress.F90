! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE READSTRESS(IJINF, IJSUP, NREAL, RFIELD, FILENAME, LRSTPARAL)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     READS THE RESTART WIND AND STRESS FIELDS

!**   INTERFACE.
!     ----------
!     CALL *READSTRESS*(NREAL, RFIELD, FILENAME)  
!     *NREAL*    NUMBER OF FIELDS IN RFIELD
!     *RFIELD*     REAL RESTART FIELDS
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE

!     METHOD.
!     -------
!     READS ARRAYS FROM FILE. FORTRAN READS

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D     ,LLUNSTR
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOWABORT,  ONLY : WAM_ABORT
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "iwam_get_unit.intfb.h"
      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, NREAL
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP,NREAL),INTENT(OUT) :: RFIELD
      CHARACTER(LEN=296), INTENT(IN) :: FILENAME
      LOGICAL, INTENT(IN) :: LRSTPARAL


      INTEGER(KIND=JWIM) :: LFILE, IUNIT
      INTEGER(KIND=JWIM) :: IFLD, IJ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP) :: RFIELD_G 

      LOGICAL :: LLEXIST

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READSTRESS',0,ZHOOK_HANDLE)

      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
      LLEXIST=.FALSE.
      INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
      IF (.NOT. LLEXIST) THEN
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (IU06,*) '*  IN READSTRESS :                  *'
        WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        WRITE (*,*) '*************************************'
        WRITE (*,*) '*                                   *'
        WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (*,*) '*  IN READSTRESS :                  *'
        WRITE (*,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (*,*) '*                                   *'
        WRITE (*,*) '*************************************'
        CALL ABORT1
      ENDIF
      IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0, 'READWRITE')

      READ(IUNIT) CDTPRO,CDATEWO,CDAWIFL,CDATEFL


      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        DO IFLD=1,NREAL
          READ(IUNIT)(RFIELD_G(IJ),IJ=IJINF,IJSUP)
          CALL UNBLKRORD(-1,IJINF,IJSUP,1,1,1,1,RFIELD(1,IFLD),RFIELD_G(1))
        ENDDO
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
!     ALL PE'S READ OR 1D DECOMPOSITION 
        DO IFLD=1,NREAL
          READ(IUNIT)(RFIELD(IJ,IFLD),IJ=IJINF,IJSUP)
        ENDDO
      ELSE
!     WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE
!     RE-LABELLED BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
        DO IFLD=1,NREAL
          READ(IUNIT)(RFIELD(IJ2NEWIJ(IJ),IFLD),IJ=IJINF,IJSUP)
        ENDDO
      ENDIF

      CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('READSTRESS',1,ZHOOK_HANDLE)

      END SUBROUTINE READSTRESS
