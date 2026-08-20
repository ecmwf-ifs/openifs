! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WAM_U2L1CR (CSTRNG)
 
!**** WAM_U2L1CR - CONVERT UPPER CASE LETTERS TO LOWER CASE.
 
!     PURPOSE.
!     --------
 
!           CONVERT CHARACTER STRING TO ALL LOWER CASE.
 
!**   INTERFACE.
!     ----------
 
!           CALL WAM_U2L1CR (CSTRNG)
 
!               INPUT PARAMETERS.
!               -----------------
 
!               CSTRNG     - STRING IN UPPER AND/OR LOWER CASE.
 
!               OUTPUT PARAMETERS.
!               ------------------
 
!               CSTRNG     - CHARACTER STRING IN ALL LOWER CASE.
 
!     METHOD.
!     -------
 
!           NUMERIC VALUES OF EACH BYTE EXAMINED AND ALTERED IF IN THE
!           RANGE OF LOWER CASE LETTERS.
 
!     EXTERNALS.
!     ----------
 
!           ICHAR
!           CHAR
 
!     AUTHOR.
!     -------
 
!           J. HENNESSY      ECMWF      09:08:90.
 
!     MODIFICATIONS
!     --------------
 
!           None.
!     ----------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
 
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: ILEN, J, ICH, IA, IZ, INULL
 
      CHARACTER(LEN=*) :: CSTRNG
 
!     ASCII Representation of upper case characters A and Z.
 
      DATA IA  /65/
      DATA IZ  /90/
 
!     ASCII REPRESENTATION OF NULL AND BLANK CHARACTERS.
 
      DATA INULL /0/
 
!     ----------------------------------------------------------------
 
      ILEN = LEN (CSTRNG)
 
      DO J=1,ILEN
        ICH = ICHAR (CSTRNG(J:J))
        IF ( (ICH.GE.IA).AND.(ICH.LE.IZ) ) CSTRNG (J:J) = CHAR (ICH+32)
        IF (ICH.EQ.INULL) CSTRNG (J:J) = CHAR (32)
      ENDDO
 
      END SUBROUTINE WAM_U2L1CR
