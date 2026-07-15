! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE CHECKCFL (KIJS, KIJL, DEPTH, DTC,               &
     &                     CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2,       &
     &                     CFLSO2, CFLTP,CFLTM,CFLOP,CFLOM)

! ----------------------------------------------------------------------

!**** *CHECKCFL* - 

!     J. BIDLOT    ECMWF   2005

!*    PURPOSE.
!     --------

!     CHECK THAT THE CFL CRITERIA ARE NOT VIOLATED ANYWHERE.

!**   INTERFACE.
!     ----------

!       *CALL* *CHECKCFL(KIJS, KIJL, DEPTH, DTC
!    &                   CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2,CFLSO2, 
!    &                   CFLTP,CFLTM,CFLOP,CFLOM)*
!          *KIJS* - INDEX OF FIRST POINT
!          *KIJL* - INDEX OF LAST POINT
!          *DEPTH* - WATER DEPTH
!          AND ALL COEFFICIENTS OF THE UPWIND SCHEME.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!       NONE.

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NSUP 
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: DTC
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2 
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CFLSO2,CFLTP,CFLTM,CFLOP,CFLOM


      INTEGER(KIND=JWIM) :: IJ, IC, ICP, ICL
      LOGICAL :: LLABORTCFL

! ----------------------------------------------------------------------

      LLABORTCFL=.FALSE.

      DO IJ = KIJS,KIJL

        IF (CFLEA(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE EAST  *'
          WRITE(IU06,*) ' * CFLEA = ',CFLEA(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLWE(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *          *                   *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE WEST  *'
          WRITE(IU06,*) ' * CFLWE = ',CFLWE(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLNO(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE NORTH *'
          WRITE(IU06,*) ' * CFLNO = ',CFLNO(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLSO(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE SOUTH *'
          WRITE(IU06,*) ' * CFLSO = ',CFLSO(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLNO2(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE NORTH *'
          WRITE(IU06,*) ' * CFLNO2 = ',CFLNO2(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLSO2(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE SOUTH *'
          WRITE(IU06,*) ' * CFLSO2 = ',CFLSO2(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLTP(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR POSITIVE DIRECTIONS      *'
          WRITE(IU06,*) ' * CFLTP = ',CFLTP(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' * IJ = ',IJ, DEPTH(IJ)
          DO ICL=1,2
            DO IC=1,2
              ICP=KLAT(IJ,IC,ICL)
              IF (ICP /= NSUP+1) WRITE(IU06,*)' * KLAT_',IC,ICL, DEPTH(ICP)
            ENDDO
          ENDDO
          DO IC=1,2
            ICP=KLON(IJ,IC)
            IF (ICP /= NSUP+1) WRITE(IU06,*)' * KLON_',IC, DEPTH(ICP)
          ENDDO
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLTM(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR NEGATIVE DIRECTIONS      *'
          WRITE(IU06,*) ' * CFLTM = ',CFLTM(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' * IJ = ',IJ, DEPTH(IJ)
          DO ICL=1,2
            DO IC=1,2
              ICP=KLAT(IJ,IC,ICL)
              IF (ICP /= NSUP+1) WRITE(IU06,*)' * KLAT_',IC,ICL, DEPTH(ICP)
            ENDDO
          ENDDO
          DO IC=1,2
            ICP=KLON(IJ,IC)
            IF (ICP /= NSUP+1) WRITE(IU06,*)' * KLON_',IC, DEPTH(ICP)
          ENDDO
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLOP(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR FREQUENCY TO THE RIGHT   *'
          WRITE(IU06,*) ' * CFLOP = ',CFLOP(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (CFLOM(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' *  ==========================  *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR FREQUENCY TO THE LEFT    *'
          WRITE(IU06,*) ' * CFLOM = ',CFLOM(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF

        IF (DTC(IJ) > 1.0_JWRB) THEN
          WRITE(IU06,*) ' ************************************'
          WRITE(IU06,*) ' *                                  *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL     *'
          WRITE(IU06,*) ' * ============================     *'
          WRITE(IU06,*) ' * STABILITY CONDITION NOT FULFILLED*'
          WRITE(IU06,*) ' * DTC = ',DTC(IJ)
          WRITE(IU06,*) ' *                                  *'
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY       *'
          WRITE(IU06,*) ' ************************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF (LLABORTCFL) THEN
          WRITE(IU06,*) ' ************************************'
          WRITE(IU06,*) ' * ABORTING IN WAVE MODEL CHECKFL !!*'
          WRITE(IU06,*) ' ************************************'
          CALL FLUSH(IU06)
          CALL ABORT1
        ENDIF
      ENDDO

      END SUBROUTINE CHECKCFL
