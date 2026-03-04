C Copyright 1981- ECMWF.
C
C This software is licensed under the terms of the Apache Licence
C Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
C
C In applying this licence, ECMWF does not waive the privileges and immunities
C granted to it by virtue of its status as an intergovernmental organisation
C nor does it submit to any jurisdiction.
C

      INTEGER FUNCTION EMOSNUM(KONOFF)
C
C---->
C**** EMOSNUM
C
C     Purpose
C     -------
C
C     Returns current EMOSLIB version number.
C
C
C     Interface
C     ---------
C
C     INUM = EMOSNUM(KONOFF)
C
C     Input
C     -----
C
C     KONOFF - switch for displayed message
C              = 0 if display is required on the first call
C              = non-zero if display is not required
C
C
C     Output
C     ------
C
C     Returns a 6-digit version number, aaabbc, where:
C       aaa = 3-digit major number
C       bb  = 2-digit minor number
C       c   = 1-digit spare number (normally 0)
C
C
C     Method
C     ------
C
C     Reads a 6-digit version number from the environment variable
C     EMOS_CYCLE. If this does not give a 6-digit number, an internal
C     hard-coded default value is used.
C
C     On the first call, the function (optionally) displays a message:
C
C       ***********************************
C       * EMOSLIB version number = aaabbc *
C       ***********************************
C
C
C     Externals
C     ---------
C
C     None.
C
C
C     Author
C     ------
C
C     J.D.Chambers     ECMWF     May 1998
C
C
C----<
C ------------------------------------------------------------------
C*    Section 0.   Variables.
C ------------------------------------------------------------------
C
C
      IMPLICIT NONE
C
#include "common/grprs.h"
C
C     Function arguments
C
      INTEGER KONOFF
C
C     Local variables
C
      INTEGER INUMBER, ICOUNT !, IOFFSET
      SAVE INUMBER, ICOUNT
      CHARACTER*42 CMESS
C     CHARACTER*20 YNUMBER
C
      DATA INUMBER/000467/, ICOUNT/0/
      DATA CMESS/'* BUFRDC/EMOSLIB version number = ****** *'/
C
C ------------------------------------------------------------------
C*    Section 1.   Initialise
C ------------------------------------------------------------------
C
      EMOSNUM = INUMBER
C
      RETURN
      END

