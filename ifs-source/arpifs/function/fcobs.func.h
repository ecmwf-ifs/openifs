! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 


!*    YOMOBFUN - STATEMENT FUNCTIONS FOR OBS. PROCESSING

!        D. VASILJEVIC    ECMWF    09/09/1991

!**   MODIFICATIONS
!     -------------
!     2007-May-30   D. Tan         Unreject functions for Aeolus DWL
!     2011-Dec-19   C. Payan       Clear functions for events, base 37t1_bf.04 and 38_main
!     ------------------------------------------------------------------

!     IBITS (INVAR,IPOS,ICOVR)=AND(SHIFTR(INVAR,IPOS),ICOVR)
!     STATEMENT FUNCTION TO UNPACK ANY BIT PATTERN

!     VARIABLE   / TYPE /                   MEANING                    /
!     --------     -----                    --------
!     IBITST         I        STATEMENT FUNCTION
!     INVAR          I        WORD TO BE UNPACKED
!     IPOS           I        BIT POSITION OF BIT STRUCT. TO BE UNPACKED
!     ICOVR          I        NO. OF BIT OCC. BY THE BIT STRUCTURE
!     not used

!     INSERT(INVAR,INVAL,IBIT)=INVAR.OR.(SHIFTL(INVAL,IBIT))
!     STATEMENT FUNCTION TO INSERT ANY BIT PATTERN

!     VARIABLE   / TYPE /                   MEANING                    /
!     --------     -----                    --------
!     INSERT         I        STATEMENT FUNCTION
!     INVAR          I        WORD TO APPLY FUN. ON
!     INVAL          I        BIT PATTERN TO BE INSERTED
!     IBIT           I        BIT POSITION

!     WE CHANGED SHIFTL T0 ISHFT AND .OR. TO IOR



!     ICHBPT(IVAR,IBP,IBO)    =IVAR.AND.
!    &                      (SHIFTL(177777777777777777777B,IBP+IBO).OR.
!                            SHIFTR(177777777777777777777B,64 -IBP)   )
!     STATEMENT FUNCTION TO CHANGE BIT PATTERN OF ONES TO ZEROS

!     VARIABLE   / TYPE /                   MEANING                    /
!     --------     -----                    --------
!     ICHBPT         I        STATEMENT FUNCTION
!     IVAR           I        INPUT WORD TO BE DEALT WITH
!     IBP            I        BIT POSITION OF BIT PATTERN
!     IBO            I        NUMBER OF BITS TAKEN BY BIT PATTERN
!     not used

!***********************************************************************
!     WIND COMPONENTS

!     UCOM(PDD,PFF)= -PFF * SIN(PDD*RADIANS)
!     VCOM(PDD,PFF)= -PFF * COS(PDD*RADIANS)

!     INSERT(INVAR,INVAL,INBIT)=INVAR.OR.(SHIFTL(INVAL,INBIT))

INTEGER(KIND=JPIM) :: INSERT,INVAR,INVAL,INBIT
INSERT(INVAR,INVAL,INBIT) = IOR(INVAR,ISHFT(INVAL,INBIT))
REAL(KIND=JPRD) :: ZINSERT, ZINVAR, ZBCLEAR
ZINSERT(ZINVAR,INVAL,INBIT) = IOR(INT(ZINVAR),ISHFT(INVAL,INBIT))
ZBCLEAR(ZINVAR,INBIT)     = IBCLR(INT(ZINVAR),INBIT)
INTEGER(KIND=JPIM) :: INTREA
REAL(KIND=JPRD) :: ZARG
INTREA(ZARG)=MIN(MAX(-2147483647,INT(ZARG)),2147483647)
REAL(KIND=JPRD) :: REAINT
INTEGER(KIND=JPIM) :: IARG
REAINT(IARG)=REAL(MIN(MAX(-2147483647,IARG),2147483647),JPRB)

INTEGER(KIND=JPIM) :: ICHSTAT_ACTIVE, ICHSTAT_PASSIVE, ICHSTAT_REJECT, ICHSTAT_BLACKLST
ICHSTAT_ACTIVE(IARG)   = INSERT(IARG, 1, 0)
ICHSTAT_PASSIVE(IARG)  = INSERT(IARG, 1, 1)
ICHSTAT_REJECT(IARG)   = INSERT(IARG, 1, 2)
ICHSTAT_BLACKLST(IARG) = INSERT(IARG, 1, 3)

REAL(KIND=JPRD) :: ZCHSTAT_ACTIVE, ZCHSTAT_PASSIVE, ZCHSTAT_REJECT, ZCHSTAT_BLACKLST, ZCHSTAT_USE_EMISKF_ONLY
ZCHSTAT_ACTIVE(ZARG)   = ZINSERT(ZARG, 1, 0)
ZCHSTAT_PASSIVE(ZARG)  = ZINSERT(ZARG, 1, 1)
ZCHSTAT_REJECT(ZARG)   = ZINSERT(ZARG, 1, 2)
ZCHSTAT_BLACKLST(ZARG) = ZINSERT(ZARG, 1, 3)
ZCHSTAT_USE_EMISKF_ONLY(ZARG) = ZINSERT(ZARG, 1, 4)

INTEGER(KIND=JPIM) :: ICHEV1, ICHEV2, IEVNT
ICHEV1(IARG, IEVNT) = INSERT(IARG, 1, IEVNT-1)
ICHEV2(IARG, IEVNT) = INSERT(IARG, 1, IEVNT-1)

REAL(KIND=JPRD) :: ZCHEV1, ZCHEV2
ZCHEV1(ZARG, IEVNT) = ZINSERT(ZARG, 1, IEVNT-1)
ZCHEV2(ZARG, IEVNT) = ZINSERT(ZARG, 1, IEVNT-1)

REAL(KIND=JPRD) :: ZCLEV1, ZCLEV2
ZCLEV1(ZARG, IEVNT) = ZBCLEAR(ZARG, IEVNT-1)
ZCLEV2(ZARG, IEVNT) = ZBCLEAR(ZARG, IEVNT-1)

INTEGER(KIND=JPIM), PARAMETER :: IREJECTED         = 4
INTEGER(KIND=JPIM), PARAMETER :: IREJECTED_REVERSE = IEOR(IREJECTED,-1)

INTEGER(KIND=JPIM) :: ICHSTAT_UNREJECT
ICHSTAT_UNREJECT(IARG) = IAND(IARG,IREJECTED_REVERSE)

REAL(KIND=JPRD) :: ZCHSTAT_UNREJECT
ZCHSTAT_UNREJECT(ZARG) = ICHSTAT_UNREJECT(INT(ZARG))
!     ------------------------------------------------------------------
