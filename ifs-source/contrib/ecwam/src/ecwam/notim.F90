! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE NOTIM (CDTWIS, CDTWIE,              &
 &                NXS, NXE, NYS, NYE, FIELDG,  &
 &                BLK2LOC, WVENVI, FF_NEXT,    &
 &                IREAD, LWCUR, NEMO2WAM)

! ----------------------------------------------------------------------

!**** *NOTIM* - STEERING MODULE IF NO TIME INTERPOLATION WANTED.

!*    PURPOSE.
!     --------

!       NOTIM NO TIME INTERPOLATION: PROCESS FORCING FIELDS.

!**   INTERFACE.
!     ----------  

!       *CALL* *NOTIM (CDTWIS, CDTWIE,
!                      NXS, NXE, NYS, NYE, FIELDG,  &
!                      BLK2LOC, WVENVI, FF_NEXT,  
!                      IREAD, LWCUR, NEMO2WAM)
!          *CDTWIS*   - DATE OF FIRST WIND FIELD.
!          *CDTWIE*   - DATE OF LAST FIRST WIND FIELD.
!          *NXS:NXE*  - FIRST DIMENSION OF FIELDG
!          *NYS:NYE*  - SECOND DIMENSION OF FIELDG
!          *FIELDG*   - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!          *BLK2LOC*  - POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *WVENVI*   - WAVE ENVIRONMENT.
!          *FF_NEXT*  - DATA STRUCTURE WITH THE NEXT FORCING FIELDS
!          *IREAD*    - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*    - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!          *NEMO2WAM* - FIELDS FRON OCEAN MODEL to WAM


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWPHYS  , ONLY : XNLEV
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWO
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDATEWL  ,CDTNEXT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "getwnd.intfb.h"
#include "incdate.intfb.h"

      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS, CDTWIE
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NEXT
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM


      INTEGER(KIND=JWIM) :: ICODE_WND

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDTWIH

      LOGICAL :: LWNDFILE, LCLOSEWND

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('NOTIM',0,ZHOOK_HANDLE)


      CDTWIH = CDTWIS

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF

! ----------------------------------------------------------------------

!*    3. LOOP OVER OUTPUT WIND TIMES.
!        ----------------------------

!*    3.1 READ ONE WIND FIELD AND TRANSFORM TO BLOCKS.
!         --------------------------------------------

!*    3.2 SAVE BLOCKED WIND FIELD.
!         ------------------------

      IF (IDELPRO > IDELWO) THEN
!*      WIND FIELD IS NOT CONSTANT FOR ONE PROPAGATION TIME STEP:
!       THIS OPTION IS NO LONGER AVAILABLE !!!
        WRITE (IU06,*) ' '
        WRITE (IU06,*) ' **********************************'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' *   FATAL ERROR IN SUB. NOTIM:   *'
        WRITE (IU06,*) ' *   ==========================   *'
        WRITE (IU06,*) ' *   THE OPTION IDELPRO > IDELWO  *'
        WRITE (IU06,*) ' *      IS LONGER AVAILABLE !!!!! *'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' *          PROGRAM ABORTS        *'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' **********************************'
        CALL ABORT1
      ELSE
!*      WIND FIELD IS CONSTANT FOR ONE PROPAGATION TIME STEP:
!       -----------------------------------------------------

        CDTNEXT=CDTWIH

        CALL GETWND (BLK2LOC,                               &
     &               NXS, NXE, NYS, NYE, FIELDG,            &
     &               WVENVI,                                &
     &               FF_NEXT,                               &
     &               CDTWIH, LWNDFILE, LCLOSEWND, IREAD,    &
     &               LWCUR, NEMO2WAM,                       &
     &               ICODE_WND)


!*      UPDATE WIND FIELD REQUEST TIME.
        CALL INCDATE (CDTWIH,IDELWO)

!*      IF TIME LEFT BRANCH NOT ALLOWED TO LOOP ANYMORE 
        IF (CDTWIH <= CDTWIE) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' ********************************************'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
          WRITE (IU06,*) ' *       ==========================         *'
          WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
          WRITE (IU06,*) ' * MORE THAN ONE LOOP NOT ALLOWED ANY LONGER*'
          WRITE (IU06,*) ' * WHEN USING ARRAY INSTEAD OF WRITE TO FILE*'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' ********************************************'
          CALL ABORT1
        ENDIF

      ENDIF

      IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)

END SUBROUTINE NOTIM
