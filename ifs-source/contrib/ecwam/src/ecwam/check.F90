! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE CHECK (IINPC)

! ----------------------------------------------------------------------

!**** *CHECK* - ROUTINE TO CHECK CONSISTENCY BETWEEN COMPUTED BLOCKS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       *CHECK* CHECKS CONSISTENCY BETWEEN BLOCK INDICES.

!**   INTERFACE.
!     ----------

!       *CALL* *CHECK (IINPC)*
!          *IINPC*   - NUMBER INPUT POINTS FROM A PREVIOUS COARSE GRID.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED
      USE YOWPCONS , ONLY : DEG
      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC   ,IJARC
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF   ,IJARF
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMAP   , ONLY : BLK2GLO  ,NGX      ,NGY      ,NIBLO,    &
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLO
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IINPC

      INTEGER(KIND=JWIM) :: IO
      INTEGER(KIND=JWIM) :: I, K, IJ, IERR
      INTEGER(KIND=JWIM) :: ILEN, IPAGE, LAST, L, IA, IE

      REAL(KIND=JWRB) :: BMOWEP, BMOEAP
      REAL(KIND=JWRB) :: GRID(NGX,NGY)

      CHARACTER(LEN=100) :: TITL
      CHARACTER(LEN=1) :: LST(NGX,NGY)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *LST*       CHARACTER  LAND SEA TABLE  L = LAND
!                                             S = SEA
!                                             + = SEA AND OUTPUT POINT.
!      *GRID*      REAL       ARRAY FOR GRIDDED PRINT OUTPUT.

! ----------------------------------------------------------------------


!*    2. GENERATE LAND SEA TABLE FROM INDEX ARRAYS.
!        ------------------------------------------

      DO K=1,NGY
        DO I=1,NGX
          LST(I,K) = 'L'
        ENDDO
      ENDDO

      IERR = 0
      DO IJ=IJS,IJL
      IF (BLK2GLO%IXLG(IJ) /= 0 .OR. BLK2GLO%KXLT(IJ) /= 0) LST(BLK2GLO%IXLG(IJ),BLK2GLO%KXLT(IJ)) = 'S'
      ENDDO


!*    2.2 INCLUDE COARSE GRID NEST OUTPUT POINTS.
!         ---------------------------------------

      IF (IBOUNC == 1) THEN
        DO IO=1,NBOUNC
          IJ = IJARC(IO)
          IF (IJ < IJS .OR. IJ > IJL) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * GRID POINT NUMBER OF OUTPUT POINT IS*'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * COARSE BOUNDARY POINT NUMBER IS IO= ',IO
            WRITE (IU06,*) ' * GRID POINT NUMBER IS   IJ = ', IJ
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR > 20) CALL ABORT1
          ENDIF
          IF (BLK2GLO%IXLG(IJ) /= 0 .OR. BLK2GLO%KXLT(IJ) /= 0)                     &
     &     LST(BLK2GLO%IXLG(IJ),BLK2GLO%KXLT(IJ)) = '/'
        ENDDO
      ENDIF

!*    2.3 INCLUDE FINE GRID NEST INPUT POINTS.
!         ------------------------------------

      IF (IBOUNF == 1) THEN
        DO IO=1,NBOUNF
          IJ = IJARF(IO)
          IF (IJ < IJS .OR. IJ > IJL) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * GRID POINT NUMBER OF OUTPUT POINT IS*'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * FINE BOUNDARY POINT NUMBER IS IO= ', IO
            WRITE (IU06,*) ' * GRID POINT NUMBER IS   IJ = ', IJ
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR > 20) CALL ABORT1
          ENDIF
          IF (BLK2GLO%IXLG(IJ) /= 0 .OR. BLK2GLO%KXLT(IJ) /= 0)         &
     &     LST(BLK2GLO%IXLG(IJ),BLK2GLO%KXLT(IJ)) = 'B'
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    5. OUTPUT OF OVERALL GRID INFORMATION.
!        -----------------------------------

 5000 CONTINUE
      WRITE (IU06,'(1H1,'' GRID SUMMERY:'')')

! ----------------------------------------------------------------------

!*    6. OUTPUT OF OPTIMAL DIMENSIONS.
!        -----------------------------

      WRITE (IU06,'(//,'' DIMENSIONS OF ARRAYS, WHICH ARE USED'',       &
     &             '' IN PRESET AND CHIEF '',/)')
      WRITE (IU06,'(''                                     DEFINED'',   &
     &           ''      USED'',''  REQUIRED'')')
      WRITE (IU06,'('' NUMBER OF DIRECTIONS        NANG '', 3I10)')     &
     &           NANG, NANG, NANG 
      WRITE (IU06,'('' NUMBER OF FREQUENCIES       NFRE '', 3I10)')     &
     &           NFRE, NFRE, NFRE 
      WRITE (IU06,'('' NUMBER OF FREQUENCIES   NFRE_RED '', 3I10)')     &
     &           NFRE_RED, NFRE_RED, NFRE_RED 
      WRITE (IU06,'('' NUMBER LONGITUDE GRID POINTS NGX '', 3I10)')     &
     &           NGX, NGX, NGX
      WRITE (IU06,'('' NUMBER LATITUDE GRID POINTS  NGY '', 3I10)')     &
     &           NGY, NGY, NGY
      WRITE (IU06,'('' MAXIMUM BLOCK LENGTH       NIBLO '', 3I10)')     &
     &           NIBLO

      WRITE (IU06,'(/,'' THE DIMENSIONS IN PRESET AND CHIEF HAVE TO '', &
     &             '' BE THE VALUES IN COLUMN - REQUIRED - '')')
      WRITE (IU06,'(  '' IF YOU WANT TO USE THE OPTIMAL DIMENSION'',    &
     &             '' LENGTH IN THE WAMODEL, THEN  '',/,                &
     &             '' RERUN PREPROC WITH THE DIMENSION'',               &
     &             '' GIVEN AS -USED-'')')


      END SUBROUTINE CHECK
