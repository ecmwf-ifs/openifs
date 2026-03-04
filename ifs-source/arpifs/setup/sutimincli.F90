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

SUBROUTINE SUTIMINCLI(KDAY,KMONTH,KYEAR,KMMCLI,KULOUT,LDFILE,PALFA)

!**** *SUTIMINCLI*  - SetUp TIMe INterpolation for CLImatology

!     Purpose.
!     --------
!           Initialize the variables needed for the time interpolation of
!           the monthly climatology data

!**   Interface.
!     ----------
!        *CALL* *SUTIMINCLI(...)*

!        Explicit arguments : 
!        ------------------
!          KDAY    - tarjet day of time interpolation
!          KMONTH  - tarjet month of time interpolation
!          KYEAR   - tarjet year of time interpolation
!          KMMCLI  - month of the climatology files
!                    KMMCLI should be set to a value out of [1,12] 
!                    if the file is missing
!          LDFILE  - .TRUE. if climatology file fits the model date
!          KULOUT  - output unit number
!          PALFA   - coefficient of time interpolation

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation
!        X = alfa . file1 + (1-alfa) . file2

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*

!     Modifications.
!     --------------
!        ORIGINAL : 00-02-15 from SUCACLIASM/DM
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDAY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMONTH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KYEAR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMMCLI(2) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
LOGICAL           ,INTENT(INOUT) :: LDFILE(2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALFA 
INTEGER(KIND=JPIM) :: IMONTHLEN(12) ! Number of days for each month of the year
INTEGER(KIND=JPIM) :: IMIDDLE       ! Middle of the target month
INTEGER(KIND=JPIM) :: IMONTH(2)     ! Value of the monthes : 1 current one ; 2 = nearest
INTEGER(KIND=JPIM) :: IFILE(2)      ! array telling which is file #1 and which is file #2
INTEGER(KIND=JPIM) :: ID, IL, JJ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*      1.  PREPARATIONS
!           ------------

!*      1.1   Compute the monthes length of the year

IF (LHOOK) CALL DR_HOOK('SUTIMINCLI',0,ZHOOK_HANDLE)
IMONTHLEN( 1)=31

IMONTHLEN( 2)=28
! Bissextile years :
IF (MOD(KYEAR,4) == 0) THEN
  IMONTHLEN(2)=29
ENDIF
! Year 2000 :
IF (MOD(KYEAR,100) == 0.AND.MOD(KYEAR,400) /= 0) THEN
  IMONTHLEN(2)=28
ENDIF

IMONTHLEN( 3)=31
IMONTHLEN( 4)=30
IMONTHLEN( 5)=31
IMONTHLEN( 6)=30
IMONTHLEN( 7)=31
IMONTHLEN( 8)=31
IMONTHLEN( 9)=30
IMONTHLEN(10)=31
IMONTHLEN(11)=30
IMONTHLEN(12)=31

!*      1.2 Compute the middle of the current month

IMIDDLE=(IMONTHLEN(KMONTH)+1)/2

!*      1.3 Compute which month is the nearest from the current one

IMONTH(1)=KMONTH
IF (KDAY < IMIDDLE) THEN
  IMONTH(2)=IMONTH(1)-1
  IF (IMONTH(2) == 0) THEN
    IMONTH(2)=12
  ENDIF
ELSE
  IMONTH(2)=IMONTH(1)+1
  IF (IMONTH(2) == 13) THEN
    IMONTH(2)=1
  ENDIF
ENDIF
WRITE(KULOUT,'('' CURRENT MONTH = '',I2.2,'' NEAREST MONTH = '',I2.2 )') &
 & IMONTH(1),IMONTH(2)   

!*      2.  TEST THE MONTH OF THE CLIMATOLOGY FILES
!           ---------------------------------------

!*      2.1 Reject the useless climatology files

DO JJ=1,2
  IF (KMMCLI(JJ) > 0 .AND. KMMCLI(JJ) < 13) THEN
    LDFILE(JJ) = (KMMCLI(JJ) == IMONTH(1) .OR. KMMCLI(JJ) == IMONTH(2))
  ELSE
    LDFILE(JJ)=.FALSE.
  ENDIF
ENDDO

IF (KMMCLI(2) == KMMCLI(1)) THEN
  LDFILE(2)=.FALSE.
ENDIF

!*      2.2 Analyse the climatology files at disposal

IF (.NOT.LDFILE(1).AND..NOT.LDFILE(2)) THEN
  WRITE(KULOUT,'('' NO CLIMATOLOGY AT DISPOSAL '')')
ELSEIF (LDFILE(1).AND..NOT.LDFILE(2)) THEN
  WRITE(KULOUT,'('' ONLY ONE CLIMATOLOGY FILE AT DISPOSAL (#1) '')')
ELSEIF (LDFILE(2).AND..NOT.LDFILE(1)) THEN
  WRITE(KULOUT,'('' ONLY ONE CLIMATOLOGY FILE  AT DISPOSAL (#2) '')')
ELSEIF (KMMCLI(1) == IMONTH(1) .AND. KMMCLI(2) == IMONTH(2)) THEN
  WRITE(KULOUT,'('' TIME INTERPOLATION OF CLIMATOLOGY FILES '')')
  IFILE(1)=1
  IFILE(2)=2
ELSEIF (KMMCLI(2) == IMONTH(1) .AND. KMMCLI(1) == IMONTH(2)) THEN
  WRITE(KULOUT,'('' TIME INTERPOLATION OF REVERSE CLIMATOLOGY FILES '')')
  IFILE(1)=2
  IFILE(2)=1
ENDIF

!*      3.  COMPUTE THE INTERPOLATION COEFFICIENT
!           -------------------------------------

IF (ALL(LDFILE)) THEN
  IF (KDAY > IMIDDLE) THEN
    IL=IMONTHLEN(KMMCLI(IFILE(1)))
  ELSE
    IL=IMONTHLEN(KMMCLI(IFILE(2)))
  ENDIF
  ID=ABS(IMIDDLE-KDAY)
  IF (IFILE(1) == 1 .AND. IFILE(2) == 2) THEN
    PALFA=REAL(IL-ID,JPRB)/REAL(IL,JPRB)
  ELSEIF(IFILE(1) == 2 .AND. IFILE(2) == 1) THEN 
!   Reverse files <=> exchange coefficient
    PALFA=REAL(ID,JPRB)/REAL(IL,JPRB)
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SUTIMINCLI',1,ZHOOK_HANDLE)
END SUBROUTINE SUTIMINCLI
