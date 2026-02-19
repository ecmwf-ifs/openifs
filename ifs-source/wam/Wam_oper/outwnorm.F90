! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

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
     &                      COUTNAME, ITOBOUT ,NIPRMOUT
      USE YOWCOUP  , ONLY : LLNORMWAMOUT_GLOBAL
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

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTWNORM',0,ZHOOK_HANDLE)

      IF (NFLAGALL .AND. NIPRMOUT > 0) THEN

        IRECV=1
        CALL MPMINMAXAVG(LLNORMWAMOUT_GLOBAL, IRECV, LDREPROD, BOUT, WNORM)

        WRITE(IU06,*) ' ' 
        WRITE(IU06,*) '  WAMNORM ON ',CDTPRO

        IF (LLNORMWAMOUT_GLOBAL) THEN
          CINFO='IRECV'
          INFO=IRECV
          IF (IRANK == IRECV) THEN
            WRITE(IU06,*) '  !!!!!!!!! REPRODUCIBLE NORMS !!!!!!'
          ELSE
            WRITE(IU06,*) '  !!!!!!!!! SEE LOG PE ', IRECV,' FOR NORMS'
          ENDIF
        ELSE
          CINFO='NPROC'
          INFO=NPROC
          IF (LDREPROD) THEN
            WRITE(IU06,*) '  !!!!!!!!! REPRODUCIBLE ONLY IF SAME NPROC!'
          ELSE
             WRITE(IU06,*) '  !!!!!!!!! NON reproducible WAMNORM !!!'
          ENDIF
        ENDIF

        WRITE(IU06,*) '          AVERAGE,              MINIMUM,  ',    &
     &   '           MAXIMUM,   NON MISSING POINTS, ',CINFO

!       WRITE NORM TO LOGFILE
        DO ITG = 1, JPPFLAG
          IF (NFLAG(ITG)) THEN
            IT = ITOBOUT(ITG)
            WRITE(IU06,*) '  WAMNORM FOR ',COUTNAME(ITG)
            WRITE(IU06,*) '  ',(WNORM(I,IT),I=1,3),INT(WNORM(4,IT)),INFO
            WRITE(IU06,111) (WNORM(I,IT),I=1,3),INT(WNORM(4,IT)),INFO
          ENDIF
        ENDDO
111     FORMAT(6x,'HEX: ',3(Z16.16,2x),1x,i8,1x,i6)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTWNORM',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTWNORM
