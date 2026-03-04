! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

FUNCTION GETOPT(Y_OPTSTR, Y_OPTARG)

USE EC_PARKIND, ONLY : JPIM
IMPLICIT NONE

INTEGER(KIND=JPIM) :: GETOPT
CHARACTER(LEN=*)   :: Y_OPTSTR, Y_OPTARG

INTEGER(KIND=JPIM) :: I,INITIAL,L_OPTSTR
INTEGER(KIND=JPIM) :: N_ARG,N_ARGS,NDX_ARG,NDX_OPT
CHARACTER(LEN=512) :: Y_ARG
CHARACTER(LEN=1)   :: Y_OPT
LOGICAL            :: LL_ENDOPTS

COMMON/GOPT_COM1/ INITIAL  ,N_ARGS  ,N_ARG   ,NDX_ARG ,LL_ENDOPTS
COMMON/GOPT_COM2/ Y_ARG

!  Initialise on 1st. call
!  -----------------------

IF (INITIAL .NE. 123456) THEN
   INITIAL    = 123456
   N_ARGS     = COMMAND_ARGUMENT_COUNT()
   N_ARG      = 0
   NDX_ARG    = 0
   LL_ENDOPTS = .FALSE.
ENDIF

Y_OPTARG = ' '

!  Get length of "y_optstr"
!  ------------------------

L_OPTSTR = 0
DO 10 I = LEN(Y_OPTSTR), 1, -1
IF (Y_OPTSTR(I : I) .NE. ' ') THEN
   L_OPTSTR = I
   GOTO 20
ENDIF
10 CONTINUE
20 CONTINUE

!  If already at end of options, return
!  ------------------------------------

IF ((N_ARG .GT. N_ARGS) .OR. LL_ENDOPTS) THEN
   LL_ENDOPTS = .TRUE.
   GETOPT = -1
   RETURN
ENDIF

!  If we need to get the next argument, do so. Check for end of options
!  --------------------------------------------------------------------

IF (NDX_ARG .EQ. 0) THEN
   N_ARG = N_ARG + 1

   IF (N_ARG .GT. N_ARGS) THEN
      LL_ENDOPTS = .TRUE.
      GETOPT     = -1
      RETURN
   ENDIF

   Y_ARG   = ' '
   CALL GETARG(N_ARG, Y_ARG)
   NDX_ARG = 1

   IF (Y_ARG(1 : 1) .NE. '-') THEN
      LL_ENDOPTS = .TRUE.
      GETOPT     = -1
      RETURN
   ELSE IF (Y_ARG .EQ. '--') THEN
      LL_ENDOPTS = .TRUE.
      N_ARG      = N_ARG + 1
      Y_ARG      = ' '
      IF (N_ARG .LE. N_ARGS) THEN
         CALL GETARG(N_ARG, Y_ARG)
      ENDIF
      GETOPT     = -1
      RETURN
   ENDIF
   NDX_ARG = 2
ENDIF

!  We have an option, now see if it is valid
!  -----------------------------------------

Y_OPT   = Y_ARG(NDX_ARG : NDX_ARG)
NDX_OPT = INDEX(Y_OPTSTR, Y_OPT)

IF (NDX_OPT .EQ. 0) THEN
   NDX_ARG = NDX_ARG + 1
   IF (Y_ARG(NDX_ARG : ) .EQ. ' ') THEN
      NDX_ARG = 0
   ENDIF
   GETOPT = ICHAR(Y_OPT)
   RETURN
ENDIF

!  We have a valid option, see if it should have an argument
!  ---------------------------------------------------------

NDX_ARG = NDX_ARG + 1
IF (Y_ARG(NDX_ARG :) .EQ. ' ') THEN
   NDX_ARG = 0
ENDIF

IF (NDX_OPT .EQ. L_OPTSTR) THEN
   GETOPT = ICHAR(Y_OPT)
   RETURN
ELSE IF (Y_OPTSTR(NDX_OPT + 1 : NDX_OPT + 1) .NE. ':') THEN
   GETOPT = ICHAR(Y_OPT)
   RETURN
ENDIF


!  A valid option with an argument
!  -------------------------------

IF (NDX_ARG .EQ. 0) THEN
   IF (N_ARG .EQ. N_ARGS) THEN
      LL_ENDOPTS = .TRUE.
      GETOPT     = 0
      RETURN
   ENDIF

   N_ARG = N_ARG + 1
   Y_ARG = ' '
   CALL GETARG(N_ARG, Y_ARG)
   NDX_ARG = 1
   LL_ENDOPTS = N_ARG .EQ. N_ARGS
ENDIF

Y_OPTARG = Y_ARG(NDX_ARG : )
NDX_ARG  = 0
GETOPT   = ICHAR(Y_OPT)

ENDFUNCTION GETOPT
