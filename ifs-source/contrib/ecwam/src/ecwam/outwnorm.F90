! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUTWNORM(LDREPROD, BOUT)

! ----------------------------------------------------------------------

!**** *OUTWNORM* - PRINTS AVERAGE, MINIMUM AND MAXIMUM VALUES OF
!                  INTEGRATED PARAMETER FIELDS 

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!     LDREPROD -  Reproducibility flag for SUMmation-operator.
!                 .TRUE. ==> Use home-written binary tree, No MPI_ALLREDUCE used.
!                 .FALSE. ==> let MPI_ALLREDUCE do the summation.
!                 NOTE that it is overuled when global norms have been reuested
!                 see namelist LLNORMWAMOUT_GLOBAL
!    *BOUT*    - OUTPUT PARAMETERS BUFFER

!     METHOD.
!     -------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NFLAG, NFLAGALL, JPPFLAG,                   &
     &                      COUTDESCRIPTION, ITOBOUT ,NIPRMOUT,         &
     &                      COUTNAME
      USE YOWCOUP  , ONLY : LLNORMWAMOUT_GLOBAL, CNORMWAMOUT_FILE
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWMPP   , ONLY : IRANK   ,NPROC
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "mpminmaxavg.intfb.h"

      LOGICAL, INTENT(IN) :: LDREPROD
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT


      INTEGER(KIND=JWIM) :: ITG, IT, I, IRECV, INFO
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(4,NIPRMOUT) :: WNORM
      CHARACTER(LEN=5) :: CINFO
      CHARACTER(LEN=80),parameter :: LH=REPEAT(' ',80)
      INTEGER :: IUNORM
      LOGICAL :: LUNORM_OPENED

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTWNORM',0,ZHOOK_HANDLE)

      IF (NFLAGALL .AND. NIPRMOUT > 0) THEN

        IUNORM = -1
        LUNORM_OPENED = .FALSE.
        IF( LEN_TRIM(CNORMWAMOUT_FILE) > 0 ) THEN
          INQUIRE( FILE=TRIM(CNORMWAMOUT_FILE), OPENED=LUNORM_OPENED )
          IF( LUNORM_OPENED ) THEN
            INQUIRE( FILE=TRIM(CNORMWAMOUT_FILE), NUMBER=IUNORM )
          ELSE
            OPEN( NEWUNIT=IUNORM, FILE=TRIM(CNORMWAMOUT_FILE) )
            ! It is guaranteed by F2008 standard that IUNORM /= -1
          ENDIF
        ENDIF

        IRECV=1
        CALL MPMINMAXAVG(LLNORMWAMOUT_GLOBAL, IRECV, LDREPROD, BOUT, WNORM)

        WRITE(IU06,*) ' ' 
        WRITE(IU06,*) '  WAMNORM ON ',CDTPRO

        IF (LLNORMWAMOUT_GLOBAL) THEN
          CINFO='IRECV'
          INFO=IRECV
          IF (IRANK == IRECV) THEN
            WRITE(IU06,*) '  !!!!!!!!! REPRODUCIBLE NORMS !!!!!!'
            IF (IUNORM /= -1 .AND. .NOT. LUNORM_OPENED) THEN
              WRITE(IUNORM,'(A)') '# REPRODUCIBLE STATISTICS (LLNORMWAMOUT_GLOBAL=T)'
            ENDIF
          ELSE
            WRITE(IU06,*) '  !!!!!!!!! SEE LOG PE ', IRECV,' FOR NORMS'
          ENDIF
        ELSE
          CINFO='NPROC'
          INFO=NPROC
          IF (LDREPROD) THEN
            WRITE(IU06,*) '  !!!!!!!!! REPRODUCIBLE ONLY IF SAME NPROC!'
            IF (IUNORM /= -1 .AND. .NOT. LUNORM_OPENED) THEN
              WRITE(IUNORM,'(A,I0,A)') '# REPRODUCIBLE STATISTICS WHEN NPROC=',NPROC,' (LLNORMWAMOUT_GLOBAL=F)'
            ENDIF
          ELSE
            WRITE(IU06,*) '  !!!!!!!!! NON reproducible WAMNORM !!!'
            IF (IUNORM /= -1 .AND. .NOT. LUNORM_OPENED) THEN
              WRITE(IUNORM,'(A,I0,A)') '# NON REPRODUCIBLE STATISTICS (LLNORMWAMOUT_GLOBAL=F, LDREPROD=F)'
            ENDIF
          ENDIF
        ENDIF

        WRITE(IU06,'(1x,3A)') '          AVERAGE,              MINIMUM,  ',    &
     &   '           MAXIMUM,   NON MISSING POINTS, ',CINFO

        IF (IUNORM /= -1) THEN
          IF (.NOT. LUNORM_OPENED) THEN
            WRITE(IUNORM,'(A)') "#"
            DO ITG = 1, JPPFLAG
              IF (NFLAG(ITG)) THEN
                IT = ITOBOUT(ITG)
                WRITE(IUNORM,'(A,  T3,A,          T19,A, A )') &
                &             '#', COUTNAME(ITG), ": ",  TRIM(COUTDESCRIPTION(ITG))
              ENDIF
            ENDDO
            WRITE(IUNORM,'(A)') "#"
            WRITE(IUNORM,'(A,  T3,A,   T19,A, T23,A,  T29,A,     T74,A,     T119,A,    T164,A )') &
            &             '#', 'DATE', 'IDX', 'NAME', 'AVERAGE (DEC, HEX)', 'MINIMUM (DEC, HEX)', &
            & 'MAXIMUM (DEC, HEX)', 'NON MISSING POINTS'
          ENDIF
          WRITE(IUNORM,'(A)') ''
        ENDIF
!       WRITE NORM TO LOGFILE
        DO ITG = 1, JPPFLAG
          IF (NFLAG(ITG)) THEN
            IT = ITOBOUT(ITG)
            WRITE(IU06,*) '  WAMNORM FOR ',TRIM(COUTDESCRIPTION(ITG))
            WRITE(IU06,*) '  ',(WNORM(I,IT),I=1,3)
            WRITE(IU06,111) (WNORM(I,IT),I=1,3),INT(WNORM(4,IT)),INFO
            IF (IUNORM /= -1) THEN
              WRITE(IUNORM,112) CDTPRO, ITG, COUTNAME(ITG), (WNORM(I,IT),REAL(WNORM(I,IT),JWRU),I=1,3),INT(WNORM(4,IT))
            ENDIF
          ENDIF
        ENDDO
111     FORMAT(5x,'HEX: ',3(Z16.16,2x),1x,i8,1x,i6)
112     FORMAT(T3,A, T19,I0, T23,A, T28,3(E23.16,' 0x',Z16,3x),T164,I0)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTWNORM',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTWNORM
