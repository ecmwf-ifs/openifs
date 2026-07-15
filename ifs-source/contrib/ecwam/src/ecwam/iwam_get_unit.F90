! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTEGER(KIND=JWIM) FUNCTION IWAM_GET_UNIT (KUSO, CDNAME, CDACCESS, CDFORM, KRECL, CDACTION)

!     FIND A FREE UNIT TO ACCESS A FILE (CDNAME) AND OPEN IT.

!--------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWABORT, ONLY : WAM_ABORT
      USE EC_LUN   , ONLY : NULERR

!--------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KUSO
      CHARACTER(LEN=*),  INTENT(IN) :: CDNAME, CDACTION
      CHARACTER(LEN=1), INTENT(IN) :: CDACCESS, CDFORM
      INTEGER(KIND=JWIM), INTENT(IN) :: KRECL

!!      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT

      INTEGER(KIND=JWIM) :: MINUNIT, MAXUNIT, ILEN, JUME

      LOGICAL :: LLOPENED, LLUSED

      CHARACTER(LEN=20) :: CLACCESS, CLFORM, CLPOSITION

      CHARACTER(len=1024) :: CERRMSG

! ----------------------------------------------------------------------

      CLACCESS   = "SEQUENTIAL"
      CLFORM     = "FORMATTED"
      CLPOSITION = "ASIS"

      MINUNIT=140
      MAXUNIT=300

      ILEN=LEN(CDNAME)
      IF ( ILEN == 0) then
        WRITE(NULERR,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
        WRITE(NULERR,*) '!!!!!!!!! LENGTH OF FILENAME=0  !!!!!!'
        CALL WAM_ABORT("Filename must not be empty",__FILENAME__,__LINE__)
      ENDIF
      ILEN=LEN_TRIM(CDNAME)
      IF ( ILEN == 0) then
        WRITE(NULERR,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
        WRITE(NULERR,*) '!!!!!!!!! FILENAME ALL BLANKS   !!!!!!'
        CALL WAM_ABORT("Filename must not be empty",__FILENAME__,__LINE__)
      ENDIF
         
      INQUIRE ( FILE=CDNAME, OPENED=LLOPENED, NUMBER=JUME)

      IF (.NOT. LLOPENED) THEN
        IF ( CDACCESS == 'D' .OR. CDACCESS == 'd') CLACCESS="DIRECT"
        IF ( CDFORM   == 'U' .OR. CDFORM   == 'u') CLFORM="UNFORMATTED"
        IF ( CDACCESS == 'A' .OR. CDACCESS == 'a') CLPOSITION="APPEND" 
 
        DO JUME=MINUNIT,MAXUNIT,1
          INQUIRE ( UNIT=JUME, OPENED=LLUSED)
          IF ( .NOT. LLUSED) EXIT
        ENDDO
 
        IF (JUME > MAXUNIT) THEN
          WRITE(NULERR,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
          WRITE(NULERR,*) '!!!!!!!!! NO FREE UNIT FOUND BETWEEN !'
          WRITE(NULERR,*) '!!!!!!!!! ',MINUNIT,' AND ',MAXUNIT
          WRITE(NULERR,*) '!!!!!!!!! CALL ABORT                 !'
          WRITE(CERRMSG,'(A,I0,A,I0,A,A,A)') "No free unit found between ",MINUNIT," and ",&
     &                                       MAXUNIT," to open file '",TRIM(CDNAME),"'"
          CALL WAM_ABORT(CERRMSG,__FILENAME__,__LINE__)
        ENDIF
        IF (KUSO >= 0 ) WRITE(KUSO,2004)        &
     &   '  IWAM_GET_UNIT: U=',JUME,              &
     &   ' F=',CDNAME(1:LEN_TRIM(CDNAME)),        &
     &   ' F=',CLFORM(1:LEN_TRIM(CLFORM)),        &
     &   ' A=',CLACCESS(1:LEN_TRIM(CLACCESS)),    &
     &   ' P=',CLPOSITION(1:LEN_TRIM(CLPOSITION))
 2004   FORMAT(A,I3,A,A,A,A,A,A,A,A)

        IF ( CLACCESS(1:6) == "DIRECT" ) THEN
          IF ( CLFORM == "FORMATTED" )THEN
             OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, RECL=KRECL, ACTION=CDACTION)
          ELSE
             OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, RECL=KRECL, ACTION=CDACTION, CONVERT='BIG_ENDIAN')
          ENDIF
        ELSE
          IF ( CLFORM == "FORMATTED" )THEN
             OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, POSITION=CLPOSITION, ACTION=CDACTION)
          ELSE
             OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, POSITION=CLPOSITION, ACTION=CDACTION, CONVERT='BIG_ENDIAN')
          ENDIF
        ENDIF
      ENDIF

      IWAM_GET_UNIT=JUME

END FUNCTION IWAM_GET_UNIT
