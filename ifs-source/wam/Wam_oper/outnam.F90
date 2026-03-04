! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction


      SUBROUTINE OUTNAM (                                                &
     &  KNANG, KNFRE,                                                    &
     &  KNGX, KNGY, KNIBLO, KNOVER, KNGOUT, KNOUTT,                      &
     &  KFRH, MFRSTLW, MLSTHG,                                           &
     &  KNMAXC, KNMAXF, KNBINP, KNIBL1, KNIBLD, KNIBLC,                  &
     &  KITAUMAX, KJUMAX, KIUSTAR, KIALPHA, KNDEPTH, KIREFRA, iper)

! ----------------------------------------------------------------------

!*       PARAMETER NAMELIST AS IN COMMON PARAM OF THE WAM-MODEL.
!        -------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) ::                                  &
     &  KNANG, KNFRE,                                                    &
     &  KNGX, KNGY, KNIBLO, KNOVER, KNGOUT, KNOUTT,                      &
     &  KFRH, MFRSTLW, MLSTHG,                                           &
     &  KNMAXC, KNMAXF, KNBINP, KNIBL1, KNIBLD, KNIBLC,                  &
     &  KITAUMAX, KJUMAX, KIUSTAR, KIALPHA, KNDEPTH, KIREFRA, iper

      INTEGER(KIND=JWIM) ::                                              &
     &  iang, ifre,                                                      &
     &  igx, igy, iiblo, iover, ioutp, ioutt,                            &
     &  ifrh, ifrstlw, ilsthg,                                           &
     &  imaxc, imaxf, ibinp, iibl1, iibld, iiblc,                        &
     &  itaumax, iumax, iustar, ialpha, idepth, irefra

      INTEGER(KIND=JWIM) :: IU06

      LOGICAL :: llper          !  the land sea mask is periodic (T)

      NAMELIST /PARWAM/                                                 &
     &  iang, ifre,                                                     &
     &  igx, igy, iiblo, iover, ioutp, ioutt,                           &
     &  ifrh, ifrstlw, ilsthg,                                          &
     &  imaxc, imaxf, ibinp, iibl1, iibld, iiblc,                       &
     &  itaumax, iumax, iustar, ialpha, idepth, irefra, llper

      IU06=6

!***  1.CHECK ON PERIODICITY OF GRID.
!     -------------------------------

      IF (IPER.EQ.1) THEN
        llper=.TRUE.
      ELSEIF (IPER.EQ.0) THEN
        llper=.FALSE.
      ELSE
        WRITE(IU06,*)' ************************************************'
        WRITE(IU06,*)'                                                 '
        WRITE(IU06,*)' SUBROUTINE  OUTNAM:                             '
        WRITE(IU06,*)'                                                 '
        WRITE(IU06,*)' IPER SHOULD BE 1 = PERIODIC GRID or             '
        WRITE(IU06,*)'                0 = NONPERIODIC GRID             '
        WRITE(IU06,*)' BUT IS >',IPER,'<                               '
        WRITE(IU06,*)'                                                 '
        WRITE(IU06,*)' PROGRAM ABORTS NOW    PROGRAM ABORTS NOW        '
        WRITE(IU06,*)' ************************************************'
        CALL ABORT1
      ENDIF

!***  2. DETERMINE PARAMETERS NAMELIST.
!     ---------------------------------

      iang    = KNANG
      ifre    = KNFRE
      igx     = KNGX
      igy     = KNGY
      iiblo   = KNIBLO
      iover   = KNOVER
      ioutp   = KNGOUT
      ioutt   = KNOUTT
      ifrh    = KFRH
      ifrstlw = MFRSTLW
      ilsthg  = MLSTHG
      imaxc   = KNMAXC
      imaxf   = KNMAXF
      ibinp   = KNBINP
      iibl1   = KNIBL1
      iibld   = KNIBLD
      iiblc   = KNIBLC
      itaumax = KITAUMAX
      iumax   = KJUMAX
      iustar  = KIUSTAR
      ialpha  = KIALPHA
      idepth  = KNDEPTH
      irefra  = KIREFRA

!***  3. WRITE TO FILE.
!     -----------------

      OPEN(UNIT=91, FILE='PARWAM')
      WRITE (91,PARWAM)
      CLOSE (UNIT=91)

      END SUBROUTINE OUTNAM
