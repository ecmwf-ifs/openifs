! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DR_HACK_MOD

USE EC_PARKIND  ,ONLY : JPIM

IMPLICIT NONE

PRIVATE

LOGICAL, SAVE,      PUBLIC :: LL_DRHACK = .FALSE. ! Will be set to .TRUE. if envvar DR_HACK=1
INTEGER(KIND=JPIM), PUBLIC :: NULDRHACK = 999 ! Output unit for drHack pseudo xml file 

PUBLIC :: DR_HACK_INIT
PUBLIC :: DR_HACK
PUBLIC :: DR_HACK_END

CONTAINS

SUBROUTINE DR_HACK_INIT()
  IMPLICIT NONE
  CHARACTER(LEN=512) :: CLENV
  CALL GET_ENVIRONMENT_VARIABLE('DR_HACK',CLENV)
  IF( CLENV == 'yes'  .OR. CLENV == 'YES'  .OR. &
    & CLENV == 'true' .OR. CLENV == 'TRUE' .OR. &
    & CLENV == 'on'   .OR. CLENV == 'ON'   .OR. &
    & CLENV == '1' ) THEN
    LL_DRHACK=.TRUE.
    IF( MYPROC() == 1 ) THEN
      OPEN (UNIT = NULDRHACK, file = "drhack.txt",position="append",action="write")
    ENDIF
  ENDIF
CONTAINS
  FUNCTION MYPROC()
    USE MPL_DATA_MODULE ,ONLY : MPL_NUMPROC
    USE MPL_MYRANK_MOD  ,ONLY : MPL_MYRANK
    INTEGER(KIND=JPIM) :: MYPROC
    IF( MPL_NUMPROC > 0 ) THEN
      MYPROC = MPL_MYRANK()
    ELSE
      MYPROC = 1
    ENDIF
  END FUNCTION MYPROC
END SUBROUTINE DR_HACK_INIT

SUBROUTINE DR_HACK_END()
  IMPLICIT NONE
  LOGICAL LOPENED 
  INQUIRE(UNIT=NULDRHACK, OPENED=LOPENED) 
  IF( LOPENED ) THEN
    CLOSE (NULDRHACK)
  ENDIF
END SUBROUTINE DR_HACK_END

SUBROUTINE DR_HACK(ROUTINE,START)
!
! Florian Suzat (METEO-FRANCE) Sept 2017 : add drHack functionality
!
! drHack documentation:
! ----------------------------------
! ARPIFS has become a huge and complicated program. Debugging it can be very
! painful especially for newbies. Documenting it is also is a huge and tedious
! job.
! The idea behind “drHack” is basically to hack drHook: using the calls 
! "IF (LHOOK) CALL DR_HOOK('XXX',I,ZHOOK_HANDLE)" 
! (where XXX is the name of a routine, and I is 0 at the beginning of the
! routine and 1 at
! the end) in order to build a big XML file describing the ARPIFS calling tree.
! At initialization, if both environmental variables DR_HOOK and DR_HACK are set
! equal to 1, 
! then the hack is activated, otherwise everything works as usual.

! IMPORTANT: for the moment, it does not work with openmp
! (need to run with openmp=1) 

! When active, we first open a file drhack.txt.
! Every time the program enters a routine, we append <ROUTINE_NAME> to the
! file, and every time the routine is left, we append </ROUTINE_NAME> (mind the
! “/” extra character).
! Then, at the end of the run, the (big!) file drhack.txt contains the calling
! tree of the MPI processor number 0 as an XML file:
! <MASTER>
! <STACK_MIX_INIT_STACK>
! <STACK_MIX_GETSTACKUSAGEX>
! </STACK_MIX_GETSTACKUSAGEX>
! </STACK_MIX_INIT_STACK>
! <CNT0>
! <GEOMETRY_MOD_GEOMETRY_SET>
! </GEOMETRY_MOD_GEOMETRY_SET>
! ....
! 
! The resulting files are not usable as is (because they are too big). But with
! a few
! lines of python, it is easy to produce a condensed version of the drhack.txt
! file
! (if you want an example script, you may ask florian.suzat@meteo.fr).
! Then, with html and javascript, these condensed files are read and a
! dynamic collapsible search tree is built.
! Illustrations of such pages can be seen at http://intra.cnrm.meteo.fr/drhack/ 
! (only from the MeteoFrance network... If you want an export, mail
! florian.suzat@meteo.fr)

! Hope this help...

! -----------------------------------------------------------------
! Different implementation of this have been tested, but this one, even if it is
! not elegant at all, is almost fast.... 

IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN) :: ROUTINE
INTEGER(KIND=JPIM),INTENT(IN) :: START
INTEGER(KIND=JPIM) :: i
CHARACTER(LEN(ROUTINE)) :: ROUTINE_CLEAN

! replace some special character
DO i = 1,LEN(ROUTINE)
  SELECT CASE (ROUTINE(i:i))
  CASE ("<")
    ROUTINE_CLEAN (i:i)="_"
  CASE (">")
    ROUTINE_CLEAN (i:i)="_"
  CASE (":")
    ROUTINE_CLEAN (i:i)="_"
  CASE (" ")
    ROUTINE_CLEAN (i:i)="_"
  CASE DEFAULT
    ROUTINE_CLEAN (i:i)=ROUTINE(i:i)
  END SELECT
END DO

IF (START==0) THEN
  WRITE(NULDRHACK,*) '<',ROUTINE_CLEAN,'>'
ELSE
  WRITE(NULDRHACK,*) '</',ROUTINE_CLEAN,'>'
  !CLOSE FILE IF LAST ROUTINE
  IF (ROUTINE_CLEAN .eq. 'MODEL_MOD_MODEL_DELETE') THEN
    CLOSE(NULDRHACK)
  ENDIF
ENDIF
END SUBROUTINE DR_HACK

END MODULE

