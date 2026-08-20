! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE HEADBC (NBOUNC, IDELPRO, TH0, FR1, IUOUT, IU06)

! ----------------------------------------------------------------------

!**** *HEADBC* - OUTPUT OF THE COARSE GRID BOUNDARY VALUE FILE HEADER.

!     R. PORTZ     MPI      JANUARY 1991
!     J. BIDLOT    ECMWF    JANUARY 1997 : TH0 IS PASSED AS AN ARGUMENT 

!*    PURPOSE.
!     --------

!       WRITE A HEADER TO THE BOUNDARY VALUE OUTPUT FILE.

!**   INTERFACE.
!     ----------

!       *CALL* *HEADBC (NBOUNC, IDELPRO, TH0, FR1, IUOUT, IU06)*
!          *NBOUNC*  - NUMBER OF OUTPUT POINTS.
!          *IDELPRO* - PROPAGATION = OUTPUT TIMESTEP (SECONDS).
!          *TH0*     - FIRST DIRECTION (IN RADIAN)
!          *FR1*     - FIRST FREQUENCY (HERTZ).
!          *IUOUT*   - OUTPUT UNIT OF BOUNDARY VALUES.
!          *IU06*    - PRINTER OUTPUT UNIT.

!     METHOD.
!     -------

!       SEQUENTIAL UNFORMATED WRITE TO UNIT.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NBOUNC, IDELPRO, IUOUT, IU06
      REAL(KIND=JWRB), INTENT(IN) :: TH0, FR1

      REAL(KIND=JWRB) :: XANG, XFRE, XBOU, XDEL

!*    1. FORMAT AND WRITE HEADER.
!        ------------------------

      XANG = REAL(NANG)
      XFRE = REAL(NFRE)
      XBOU = REAL(NBOUNC)
      XDEL = REAL(IDELPRO)
      WRITE(IUOUT,ERR=2000) XANG, XFRE, TH0, FR1, FRATIO, XBOU, XDEL

      RETURN

! ----------------------------------------------------------------------

!*    2. ERROR MESSAGE.
!        --------------
2000  CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*    FATAL ERROR IN SUB. HEADBC        *'
      WRITE(IU06,*) '*    ===========================       *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*    WRITE ERROR FROM UNIT : IUOUT = ', IUOUT
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        PROGRAM ABORTS                *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

      END SUBROUTINE HEADBC
