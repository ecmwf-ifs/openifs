! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MBOUNC (IU09, IU19, IFORM)

! ----------------------------------------------------------------------

!**** *MBOUNC* - MAKE COARSE GRID BOUNDARY.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     -------

!       COMPUTE ALL INFORMATION FOR COARSE GRID BOUNDARY VALUE
!       OUTPUT (COMMON CBOUND).

!**   INTERFACE.
!     ----------

!       *CALL* *MBOUNC (IU09, IU19, IFORM)*
!          *IU09*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
!          *IU19*   - LOGICAL UNIT FOR    FORMATED WRITE.
!          *IFORM*  - FORMAT OPTION  = 1  UNFORMATED WRITE.
!                                    = 2  FORMATED WRITE.
!                                      OTHERWISE BOTH.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *FINDB*     - FIND BLOCK AND SEA POINT NUMBERS.
!       *MBOXB*     - FIND LAT AND LONG OF BOUNDARY POINTS.
!       *PACKI*     - PACKS AN INTEGER ARRAY.
!       *PACKC*     - PACKS A  REAL    ARRAY.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC    ,IGARC    ,DLAMAC   ,    &
     &            DPHIAC   ,AMOSOC   ,AMONOC   ,AMOEAC   ,AMOWEC   ,    &
     &            BLATC    ,BLNGC    ,GBOUNC   ,IPOGBO
      USE YOWMAP   , ONLY : NGX      ,NGY      ,XDELLA   ,XDELLO
      USE YOWPRPROC,ONLY : NBMAX
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "findb.intfb.h"
#include "mboxb.intfb.h"
#include "packi.intfb.h"
#include "packr.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU09, IU19, IFORM

      INTEGER(KIND=JWIM) :: I, II, IO, INI 
      INTEGER(KIND=JWIM) :: NBOUNEW
      INTEGER(KIND=JWIM) :: ListSTART(1), ListEND(1)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IPOGBN(:)

! ----------------------------------------------------------------------

  998 FORMAT(10I8)
  999 FORMAT(5E16.7)

! ----------------------------------------------------------------------

!*    1. INITIAL.
!        --------
      
      DPHIAC = XDELLA
      DLAMAC = XDELLO

      NBMAX = GBOUNC*((NGX+NGY)*2-4)

      ALLOCATE(IPOGBO(0:GBOUNC))
      ALLOCATE(IPOGBN(0:GBOUNC))
      IF (ALLOCATED(IJARC)) DEALLOCATE(IJARC)
      ALLOCATE(IJARC(NBMAX))
      IF (ALLOCATED(IGARC)) DEALLOCATE(IGARC)
      ALLOCATE(IGARC(NBMAX))
      IF (ALLOCATED(BLATC)) DEALLOCATE(BLATC)
      ALLOCATE(BLATC(NBMAX))
      IF (ALLOCATED(BLNGC)) DEALLOCATE(BLNGC)
      ALLOCATE(BLNGC(NBMAX))

      IPOGBN(0) = 0
      IPOGBO(0) = 0
      DO I = 1,NBMAX
        IJARC(I) = 0
        ! IGARC obsolete
        IGARC(I) = 1
        BLATC(I) = 0.0_JWRB
        BLNGC(I) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTED THE SQUARE BOX.
!        ------------------------
      DO II=1,GBOUNC

        CALL MBOXB (NBOUNC,AMOWEC(II),AMOSOC(II),AMOEAC(II),AMONOC(II), &
     &   BLATC(IPOGBN(II-1)+1), BLNGC(IPOGBN(II-1)+1))
       
        IPOGBN(II)=IPOGBN(II-1)+NBOUNC
! ----------------------------------------------------------------------

!*    3. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!        -----------------------------------------
        ListSTART(1)=1
        ListEND(1)=1
        CALL FINDB (NBMAX, NBOUNC, BLATC(IPOGBN(II-1)+1),               &
     &    BLNGC(IPOGBN(II-1)+1),                                        &
     &    IJARC(IPOGBN(II-1)+1),                                        &
     &     1, ListSTART,ListEND, 1, 1)

! ----------------------------------------------------------------------

!*    4. PACKED ALL ARRAYS.
!        -------------------

        CALL PACKR (NBOUNC, NBOUNEW, NBMAX, IGARC(IPOGBN(II-1)+1),      &
     &    BLATC(IPOGBN(II-1)+1))
        CALL PACKR (NBOUNC, NBOUNEW, NBMAX, IGARC(IPOGBN(II-1)+1),      &
     &    BLNGC(IPOGBN(II-1)+1))
        CALL PACKI (NBOUNC, NBOUNEW, NBMAX, IGARC(IPOGBN(II-1)+1),      &
     &   IJARC(IPOGBN(II-1)+1))
        CALL PACKI (NBOUNC, NBOUNEW, NBMAX, IGARC(IPOGBN(II-1)+1),      &
     &   IGARC(IPOGBN(II-1)+1))
        IPOGBO(II) = IPOGBO(II-1)+NBOUNEW
        IPOGBN(II) = IPOGBO(II)

! ----------------------------------------------------------------------

!*    5. PRINTER PROTOCOL.
!        -----------------

        WRITE (IU06,'(''1BOUNDARY OUTPUT POINTS OF THIS GRID:'')')
        WRITE (IU06,*) '===================================='
        WRITE (IU06,*) '  NUMBER OF BOUNDARY POINTS IS NBOUNC = ',      &
     &    NBOUNEW
        WRITE (IU06,'(4X,''     |-----INPUT-----|'',                    &
     &                 ''-POINT INDEX--|'')')
        WRITE (IU06,'(4X,''  NO.    LAT.   LONG.  BLOCK.'',             &
     &              ''  POINT.'')')
        DO IO=IPOGBO(II-1)+1,IPOGBO(II)
          WRITE (IU06,'(4X,I5,2F8.2,2I8)')                              &
     &     IO-IPOGBO(II-1), BLATC(IO), BLNGC(IO),                       &
     &     IGARC(IO), IJARC(IO)
        ENDDO

      ENDDO

! ----------------------------------------------------------------------

!*    6. WRITE COMMON CBOUND.
!        --------------------
      IF (IFORM.NE.2) THEN
        WRITE(IU09) GBOUNC
        WRITE(IU09) (IPOGBN(I),I=1,GBOUNC)
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE(IU19,998) GBOUNC
        WRITE(IU19,998) (IPOGBN(I),I=1,GBOUNC)
      ENDIF
 
      DO II=1,GBOUNC
        NBOUNC = IPOGBN(II)
        INI = IPOGBN(II-1) + 1
        IF (IFORM.NE.2) THEN
          WRITE(IU09)NBOUNC-INI+1
          WRITE(IU09)(IGARC(I),I=INI,NBOUNC)
          WRITE(IU09)(IJARC(I),I=INI,NBOUNC)
          WRITE(IU09)DLAMAC, DPHIAC, AMOSOC(II), AMONOC(II),            &
     &     AMOEAC(II), AMOWEC(II),                                      &
     &     (BLNGC(I),I=INI,NBOUNC), (BLATC(I),I=INI,NBOUNC)
        ENDIF
        IF (IFORM.NE.1) THEN
          WRITE(IU19,998)NBOUNC-INI+1
          WRITE(IU19,998)(IGARC(I),I=INI,NBOUNC)
          WRITE(IU19,998)(IJARC(I),I=INI,NBOUNC)
          WRITE(IU19,999)DLAMAC,DPHIAC, AMOSOC(II), AMONOC(II),         &
     &     AMOEAC(II), AMOWEC(II),                                      &
     &     (BLNGC(I),I=INI,NBOUNC), (BLATC(I),I=INI,NBOUNC)
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

      END SUBROUTINE MBOUNC
