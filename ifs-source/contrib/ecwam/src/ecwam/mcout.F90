! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MCOUT

! ----------------------------------------------------------------------

!**** *MCOUT* - ROUTINE TO COMPUTE OUTPUT INDICES (COMMON COUT).

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       *MCOUT* COMPUTES THE INDICES OF SPECTRA OUTPUT POINTS.

!**   INTERFACE.
!     ----------

!       *CALL* *MCOUT*

!     METHOD.
!     -------

!       THE LATITUDE AND LOGITUDE ARE CONVERTED TO INDICES.

!     EXTERNALS.
!     ----------

!       *FINDB*     - FIND BLOCK AND SEA POINT NUMBERS.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NGOUT    ,IJAR     ,OUTLONG  ,OUTLAT 
      USE YOWMAP   , ONLY : BLK2GLO  ,AMOWEP   ,AMOSOP   ,XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : NPROC
      USE YOWSPEC,   ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "findb.intfb.h"

      INTEGER(KIND=JWIM) :: IO, KX, IX, NG
      INTEGER(KIND=JWIM) :: NGOUTNEW
      INTEGER(KIND=JWIM), ALLOCATABLE :: IDUM(:)

      REAL(KIND=JWRB) :: ALONG, ALAT



!*    1. NO OUTPUT POINTS SPECIFIED.
!        ---------------------------

      IF (NGOUT == 0) THEN
        WRITE(IU06,'(1H1,'' SPECIAL OUTPUT POINTS FOR SPECTRA:'')')
        WRITE(IU06,*) 'OUTPUT POINTS ARE NOT DEFINED IN USER INPUT'
        RETURN
      ELSE
        ALLOCATE(IJAR(NGOUT))
      ENDIF

! ----------------------------------------------------------------------

!*    2. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!        -----------------------------------------

      CALL FINDB (NGOUT, NGOUT, OUTLAT, OUTLONG, IJAR,            &
     &            NPROC, NSTART, NEND, 1, NPROC)

! ----------------------------------------------------------------------

!*    3. PRINTER PROTOCOL.
!        -----------------

      WRITE(IU06,'(1H1,'' SPECIAL OUTPUT POINTS FOR SPECTRA:'')')
      WRITE(IU06,'(''    NUMBER OF OUTPUT POINTS IS NGOUT = '',I4)') NGOUT
      WRITE(IU06,'(4X,''  NO.    LAT.   LONG.    LAT.   LONG.  INDEX.'')')
      DO IO=1,NGOUT
        IF (IJAR(IO) > 0) THEN
          IX  = BLK2GLO%IXLG(IJAR(IO))
          KX  = BLK2GLO%KXLT(IJAR(IO))
          ALONG = AMOWEP + (IX-1)*ZDELLO(KX)
          ALAT  = AMOSOP + (KX-1)*XDELLA
        ELSE
          ALONG = 9999999
          ALAT  = 9999999
        ENDIF
        WRITE(IU06,'(4X,I5,4F8.2,I8)') IO, OUTLAT(IO), OUTLONG(IO), ALAT, ALONG, IJAR(IO)
      ENDDO

! ----------------------------------------------------------------------

!*    4. REMOVE OUTPUT POINTS WHICH ARE NOT IN GRID.
!        -------------------------------------------

      NG = 0
      DO IO=1,NGOUT
        IF (NG > 0 .AND. IO-NG > 0) THEN
          IJAR(IO-NG) = IJAR(IO)
        ENDIF
        IF (IJAR(IO) == 0) NG = NG+1
      ENDDO
      NGOUTNEW = NGOUT-NG
      IF (NG > 0) THEN
        WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
        WRITE (IU06,*) ' +                                       +'
        WRITE (IU06,*) ' +     WARNING ERROR FROM SUB. MCOUT     +'
        WRITE (IU06,*) ' +     =============================     +'
        WRITE (IU06,*) ' +                                       +'
        WRITE (IU06,*) ' + NO SEAPOINT FOUND FOR NG = ',NG
        WRITE (IU06,*) ' + OUTPUT POINT REQUESTS (SEE ABOVE LIST)+'
        WRITE (IU06,*) ' + THESE POINTS WILL NOT BE TAKEN.       +'
        WRITE (IU06,*) ' + NUMBER OF OUTPUT POINTS IS NGOUT = ',NGOUTNEW
        WRITE (IU06,*) ' +                                       +'
        WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
        ALLOCATE(IDUM(NGOUTNEW))
        DO IO=1,NGOUTNEW
           IDUM(IO) = IJAR(IO)
        ENDDO
        DEALLOCATE(IJAR)
        ALLOCATE(IJAR(NGOUTNEW))
        DO IO=1,NGOUTNEW
           IJAR(IO) = IDUM(IO)
        ENDDO
        DEALLOCATE(IDUM)
        NGOUT = NGOUTNEW
      ENDIF

      END SUBROUTINE MCOUT
