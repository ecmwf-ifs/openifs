! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUTPP (CDATE, IUOUT, ID1, ID2, NGX, NGY, TITL, CONST,  &
     &                  ARRAY, AMOWEP, AMOSOP, AMOEAP, AMONOP)

! ---------------------------------------------------------------

!**** *OUTPP* - FORMATED OUTPUT OF AN ARRAY.

!     H. GUNTHER       ECMWF    NOVEMBER 1989

!*    PURPOSE.
!     --------

!       FORMATED OUTPUT OF AN ARRAY.

!**   INTERFACE.
!     ---------

!       *CALL* *OUTPP (CDATE, IUOUT, ID1, ID2, NGX, NGY, TITL, CONST,
!                      ARRAY, AMOWEP, AMOSOP, AMOEAP, AMONOP)*
!          *CDATE*    INTEGER   DATE (YYYYMMDDHHMM).
!          *IUOUT*    INTEGER   OUTPUT UNIT.
!          *ID1*      INTEGER   FIRST DIMENSION.
!          *ID2*      INTEGER   SECOND DIMENSION.
!          *NGX*      INTEGER   FIRST DIMENSION  USED.
!          *NGY*      INTEGER   SECOND DIMENSION USED.
!          *TITL*     CHARACTER HEADER TO BE PRINTED.
!          *CONST*    REAL      SCALING FACTOR.
!          *ARRAY*    REAL      ARRAY TO BE PRINTED.
!          *AMOWEP*   REAL      MOST WESTERN LONGITUDE (DEGREE).
!          *AMOSOP*   REAL      MOST SOUTHERN LATITUDE (DEGREE).
!          *AMOEAP*   REAL      MOST EASTERN LONGITUDE (DEGREE).
!          *AMONOP*   REAL      MOST NORTHERN LATITUDE (DEGREE).

!     METHOD.
!     -------

!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF
!       30 COLUMNS PER PAGE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IUOUT, ID1, ID2, NGX, NGY

      REAL(KIND=JWRB), INTENT(IN) :: CONST
      REAL(KIND=JWRB), INTENT(IN) :: AMOWEP, AMOSOP, AMOEAP, AMONOP
      REAL(KIND=JWRB), DIMENSION(ID1,ID2), INTENT(IN) :: ARRAY

      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      CHARACTER(LEN=100), INTENT(IN) :: TITL


      INTEGER(KIND=JWIM) :: I, J, NP
      INTEGER(KIND=JWIM) :: NPTS, NPAGE, ISTART, IEND
      INTEGER(KIND=JWIM), DIMENSION(NGX) :: ILON
      INTEGER(KIND=JWIM), DIMENSION(ID1,ID2) :: IARRAY

      REAL(KIND=JWRB) :: DLAMA, DPHIA
      REAL(KIND=JWRB), DIMENSION(NGY) :: YLAT

! ----------------------------------------------------------------------


!*    1. INITIALIZATION
!     -----------------

      IF(NGX.NE.1) THEN
        DLAMA = (AMOEAP-AMOWEP)/REAL(NGX-1)
      ELSE
        DLAMA = 0
      ENDIF
      DO I=1,NGX
        ILON(I)=NINT(AMOWEP + (I-1)*DLAMA)
      ENDDO

      IF(NGY.NE.1) THEN
        DPHIA = (AMOSOP-AMONOP)/REAL(NGY-1)
      ELSE
        DPHIA = 0
      ENDIF
      DO J=1,NGY
        YLAT(J)=AMONOP + (J-1)*DPHIA
      ENDDO

      NPTS=30
      NPAGE=(NGX+NPTS-1)/NPTS

      DO J=1,NGY
        DO I=1,NGX
          IARRAY(I,J) = NINT(CONST*ARRAY(I,J))
        ENDDO
      ENDDO

      ISTART=-NPTS+1
      IEND=ISTART+NPTS-1

      DO NP=1,NPAGE
        WRITE (IUOUT,300) CDATE,TITL,NP
        ISTART=ISTART+NPTS
        IEND= MIN(IEND+NPTS,NGX)
        WRITE (IUOUT,302) (I,I=ISTART,IEND)
        WRITE (IUOUT,303) (ILON(I),I=ISTART,IEND)
        WRITE (IUOUT,304)
        DO J=1,NGY
          WRITE (IUOUT,305) J,YLAT(J),(IARRAY(I,J),I=ISTART,IEND)
        ENDDO
      ENDDO
      WRITE (IUOUT,*) ' '

 300  FORMAT('1',6X,A14,2X,A100,5X,'PAGE ',I2,/)
 302  FORMAT(7X,'I=',30I4)
 303  FORMAT(5X,'LON=',30I4)
 304  FORMAT('   J LAT',/)
 305  FORMAT(1X,I2,F5.1,1X,30I4)

      END SUBROUTINE OUTPP
