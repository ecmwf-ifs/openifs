! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE READBOU (IU09, IU10, IU06)

! ----------------------------------------------------------------------

!**** *READBOU* - READ COARSE AND FINE GRID BOUNDARY INFORMATION FILE.

!     R. PORTZ     MPI          JANUARY 1991

!*    PURPOSE.
!     --------

!       READ COMMON CBOUND AND FBOUND AS WRITTEN BY PREPROC.

!**   INTERFACE.
!     ----------

!       *CALL* *READBOU (IU09, IU10, IU06)*
!          *IU09*    - INPUT  UNIT OF COMMON CBOUND.
!          *IU10*    - INPUT  UNIT OF COMMON FBOUND.
!          *IU06*    - PRINTER OUTPUT UNIT.

!     METHOD.
!     -------

!       SEQUENCIAL UNFORMATED WRITE TO UNIT.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC   ,IJARC    ,IGARC,        &
     &            GBOUNC, IPOGBO
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF   ,IJARF    ,IGARF    ,    &
     &            IBFL     ,IBFR     ,BFW
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D
      USE YOWSPEC  , ONLY : IJ2NEWIJ

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU09, IU10, IU06


      INTEGER(KIND=JWIM) :: I, II
      INTEGER(KIND=JWIM) :: NBOUNC_LOC

      REAL(KIND=JWRB) :: DLA, DPH, AMOS, AMON, AMOE, AMOW
      REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: BLNGC, BLATC

! ----------------------------------------------------------------------

      IF (IBOUNC.EQ.1) THEN
        READ(IU09, ERR=2000) GBOUNC
        ALLOCATE(IPOGBO(0:GBOUNC))
        READ(IU09, ERR=2000)(IPOGBO(I),I=1,GBOUNC)
        NBOUNC=IPOGBO(GBOUNC)
        IPOGBO(0)=0
        ALLOCATE(IGARC(NBOUNC))
        ALLOCATE(IJARC(NBOUNC))
        ALLOCATE(BLNGC(NBOUNC))
        ALLOCATE(BLATC(NBOUNC))
!       READ COMMON BOUNC
        DO II=1,GBOUNC
          READ(IU09, ERR=2000) NBOUNC_LOC
          READ (IU09,ERR=2000) (IGARC(IPOGBO(II-1)+I),I=1,NBOUNC_LOC)
          READ (IU09,ERR=2000) (IJARC(IPOGBO(II-1)+I),I=1,NBOUNC_LOC)
          READ(IU09,ERR=2000) DLA, DPH, AMOS, AMON, AMOE, AMOW,         &
     &     (BLNGC(I),I=1,NBOUNC_LOC), (BLATC(I),I=1,NBOUNC_LOC)
        ENDDO
      ENDIF

      IF (IBOUNF.EQ.1) THEN

!     READ COMMON BOUNF

        READ (IU10,ERR=2000) NBOUNF

        ALLOCATE(IGARF(NBOUNF))
        ALLOCATE(IJARF(NBOUNF))
        ALLOCATE(IBFL(NBOUNF))
        ALLOCATE(IBFR(NBOUNF))
        ALLOCATE(BFW(NBOUNF))

        READ (IU10,ERR=2000) (IGARF(I),I=1,NBOUNF),                     &
     &   (IJARF(I),I=1,NBOUNF),                                         &
     &   (IBFL(I),I=1,NBOUNF),                                          &
     &   (IBFR(I),I=1,NBOUNF),                                          &
     &   (BFW(I),I=1,NBOUNF)

!       CONVERT IJARF TO NEW ORDER IF 2D DECOMPOSITION
        IF(.NOT.LL1D .AND. NPROC.GT.1 ) THEN
          DO I=1,NBOUNF
            IJARF(I)=IJ2NEWIJ(IJARF(I))
          ENDDO
        ENDIF

      ENDIF
      RETURN
 2000 CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*    FATAL ERROR IN SUB. READBOU       *'
      WRITE(IU06,*) '*    ===========================       *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        READ ERROR                    *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        PROGRAM ABORTS                *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

      END SUBROUTINE READBOU
