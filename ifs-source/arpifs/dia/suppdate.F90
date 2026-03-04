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

SUBROUTINE SUPPDATE(YDRIP,KSTEP,KDTIME,KEVENTTS,KEVENTTSMIN,KDIME,KFREVENT,K1EVENT,KPPDATE)

!**** *SUPPDATE*  - Setup post-processind date

!     Purpose.
!     --------
!           To initialize the date to write out on an output file
!           (historical file or post-processing file)

!**   Interface.
!     ----------
!        *CALL* *SUPPDATE(...)*

!        Explicit arguments :
!        ------------------
!         KSTEP    - model time step
!         KDTIME   - kind of date for the file : 
!                    0 <=> date of the model
!                    1 <=> date of the starting file + model forecast range
!                   (resulting date will differ if the starting file is already
!                    a forecast)
!                    2 <=> date of the input file (unchanged)
!         KEVENTTS - array of time step events
!         KDIME    - dimension of KEVENTTS array
!         KFREVENT - frequency events
!         K1EVENT  - overriding switch for events
!         KPPDATE  - date in output file

!        Implicit arguments :
!        --------------------
!        See 'USE MODULE' above.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        ABOR1

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 99-07-20 from NEWFA

!     Modifications.
!     --------------
!      R. El Khatib : 01-03-16 : KPPDATE(10)=time of the previous event if any ;
!                    =0 else
!      R. El Khatib : 03-09-12 : fix for negative NSTEP
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 05-02-21 : bugfixes on this mess ...
!      F. Taillefer : 21-03-2013 : add case 2 for KDTIME
!      P. Marguinaud: 10-10-2013 : Avoid reading again input file
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NINISH, NULOUT
USE YOMCT0   , ONLY : NCONF
USE YOMRIP0  , ONLY : NINDAT, NSSSSS
USE YOMRIP   , ONLY : TRIP
USE YOMINI   , ONLY : LDFI
USE YOMOPH0  , ONLY : CFNISH, NTIMEFMT 
USE FA_MOD   , ONLY : &
                      & JD_YEA, JD_MON, JD_DAY, &  
                      & JD_HOU, JD_MIN, JD_TUN, &  
                      & JD_THO, JD_GR8, JD_IAN, &  
                      & JD_CU1, JD_CU2,         &  
                      & JD_DEX,         JD_SEM, &  
                      & JD_SET, JD_CE1, JD_CE2, &  
                      & JD_TST, JD_FMT,         &  
                      & JD_SIZ

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP)        ,INTENT(IN)           :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)           :: KSTEP
INTEGER(KIND=JPIM),INTENT(IN)           :: KDIME 
INTEGER(KIND=JPIM),INTENT(IN)           :: KDTIME 
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KEVENTTS(0:KDIME) 
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KEVENTTSMIN(0:KDIME) 
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KFREVENT
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: K1EVENT
INTEGER(KIND=JPIM),INTENT(OUT)          :: KPPDATE(:) 

!     ------------------------------------------------------------------
#include "abor1.intfb.h"
#include "monio_t.intfb.h"
#include "fcttim.func.h"

INTEGER(KIND=JPIM), ALLOCATABLE :: IEVENTS(:)

INTEGER(KIND=JPIM) :: IDATEM (JD_SIZ) ! Model
INTEGER(KIND=JPIM) :: IDATEF (JD_SIZ) ! File
INTEGER(KIND=JPIM) :: IDATEP (JD_SIZ) ! Result

INTEGER(KIND=JPIM) :: ISS0, J, II

REAL (KIND=JPRB) :: ZTT_PRE
INTEGER (KIND=JPIM) :: ISS_PRE, IHH_PRE, ISTEP_PRE
REAL (KIND=JPRB) :: ZTT_NOW
INTEGER (KIND=JPIM) :: ISS_NOW, IHH_NOW, ISTEP_NOW

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUPPDATE',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IDATEP = 0

IDATEM = 0
IDATEM (JD_YEA) = NCCAA(NINDAT)
IDATEM (JD_MON) = NMM(NINDAT)
IDATEM (JD_DAY) = NDD(NINDAT)
IDATEM (JD_HOU) = NSSSSS / 3600
IF (NTIMEFMT > 0) THEN
  IDATEM (JD_MIN) = (NSSSSS - IDATEM (JD_HOU) * 3600) / 60
ELSE
  IDATEM (JD_MIN) = 0
ENDIF
IDATEM (JD_DEX) = 1
IDATEM (JD_SEM) = NSSSSS
IDATEM (JD_SET) = INT (YDRIP%RSTATI, JPIM)
IDATEM (JD_CE1) = 0
IDATEM (JD_CE2) = 0
IDATEM (JD_TST) = YDRIP%TSTEP
IDATEM (JD_FMT) = NTIMEFMT

ISTEP_NOW = KSTEP

!*    1.    DECIDE SUITABLE DATE
!           --------------------

SELECT CASE (KDTIME) 

  CASE (0)

    IDATEP = IDATEM 
    ISS0 = 0

  CASE (1) 

    CALL SUPPDATE_READ (CFNISH, NINISH, IDATEF)
    
    IDATEP (JD_YEA:JD_MIN) = IDATEF (JD_YEA:JD_MIN)
    IDATEP (JD_DEX:JD_SIZ) = IDATEF (JD_DEX:JD_SIZ)

    ISS0 = IDATEF (JD_SET)

  CASE (2) 

    CALL SUPPDATE_READ (CFNISH, NINISH, IDATEF)
    IDATEP = IDATEF

  CASE DEFAULT

    CALL ABOR1('SUPPDATE: UNKNOWN VALUE KDTIME')

END SELECT

IF (KDTIME < 2) THEN

  ZTT_NOW = REAL (ISTEP_NOW, JPRB) * YDRIP%TSTEP + REAL (ISS0, JPRB)
  ISS_NOW = NINT (ZTT_NOW)
  IHH_NOW = NINT (ZTT_NOW / 3600._JPRB)

  IF (ISS_NOW <= 0) THEN

    IF (KDTIME == 0) THEN
      IDATEP (JD_TUN:JD_GR8) = (/ 1, 0, 0 /)
      IDATEP (JD_SET) = 0
    ELSE
      IDATEP (JD_TUN:JD_GR8) = IDATEF (JD_TUN:JD_GR8)  
      IDATEP (JD_SET) = IDATEF (JD_SET)
    ENDIF

    IF (NCONF == 1 .OR. NCONF == 302) THEN
!     Time for 4D-Var 
      IDATEP (JD_IAN) = 1
    ELSE
!     Old times of analysis ...
      IF (NCONF == 701 .OR. NCONF == 131) THEN
!       Unitialised analysis : 
        IDATEP (JD_IAN) = 0
      ELSEIF (.NOT. LDFI) THEN
!       If the initial file is initialized, then the output file is 
!       necessarily initialized (we cannot remove the initialization of 
!       initialized data !
        IF (KDTIME == 0) THEN
          CALL SUPPDATE_READ (CFNISH, NINISH, IDATEF)
        ENDIF
        IDATEP (JD_IAN) = IDATEF (JD_IAN)
      ELSE
!       Initialised analysis in any case : 
        IDATEP (JD_IAN) = 1
      ENDIF
    ENDIF

    IF (KDTIME == 0) THEN
      IDATEP (JD_CU1) = 0
      IDATEP (JD_CU2) = 0
      IDATEP (JD_CE1) = 0
      IDATEP (JD_CE2) = 0
    ELSE
      IDATEP (JD_CU1) = IDATEF (JD_CU1)
      IDATEP (JD_CU2) = IDATEF (JD_CU2)
      IDATEP (JD_CE1) = IDATEF (JD_CE1)
      IDATEP (JD_CE2) = IDATEF (JD_CE2)
    ENDIF

  ELSE

!   Forecast : 
    IF (IHH_NOW < 65000) THEN
!     In hours : 
      IDATEP (JD_TUN) = 1
      IDATEP (JD_THO) = IHH_NOW
    ELSE
!     In days :
      IDATEP (JD_TUN) = 2
      IDATEP (JD_THO) = IHH_NOW / 24
    ENDIF

    IDATEP (JD_GR8) = 0
    IDATEP (JD_IAN) = 10

    IDATEP (JD_SET) = ISS_NOW

    IF (ISTEP_NOW <= YDRIP%NSTART) THEN

      IF (KDTIME == 1) THEN
        IDATEP (JD_CU1) = IDATEF (JD_CU1)
        IDATEP (JD_CU2) = IDATEF (JD_CU2)
        IDATEP (JD_CE1) = IDATEF (JD_CE1)
        IDATEP (JD_CE2) = IDATEF (JD_CE2)
      ELSE
        IDATEP (JD_CU1) = 0
        IDATEP (JD_CU2) = 0
        IDATEP (JD_CE1) = 0
        IDATEP (JD_CE2) = 0
      ENDIF

    ELSEIF (PRESENT(KEVENTTS) .AND. PRESENT(K1EVENT) .AND. PRESENT(KFREVENT)) THEN

!     Compute previous time of output since model start
      ALLOCATE (IEVENTS (0:YDRIP%NSTOP/KFREVENT))
      IEVENTS = 0

      CALL MONIO_T(0,YDRIP,IEVENTS,K1EVENT,KFREVENT,KEVENTTS,&
                  & KN___TSMIN=KEVENTTSMIN)

      II=0
      DO J = ISTEP_NOW/KFREVENT-1, YDRIP%NSTART/KFREVENT, -1
        II=1
        IF (IEVENTS (J) /= 0) EXIT
      ENDDO
      ISTEP_PRE = J * KFREVENT * II

      DEALLOCATE (IEVENTS)

      ZTT_PRE = REAL (ISTEP_PRE, JPRB) * YDRIP%TSTEP + REAL (ISS0, JPRB)
      ISS_PRE = NINT (ZTT_PRE)
      IHH_PRE = NINT (ZTT_PRE / 3600._JPRB)

      IF (IDATEP (JD_TUN) == 1) THEN
!       In hours :
        IDATEP (JD_CU1) = IHH_PRE
      ELSEIF(IDATEP (JD_TUN) == 2) THEN
!       In days :
        IDATEP (JD_CU1) = IHH_PRE / 24
      ENDIF

      IDATEP (JD_CE1) = ISS_PRE

      IDATEP (JD_CU2) = 0
      IDATEP (JD_CE2) = 0

    ELSE
!     No past event, put 0.
      IDATEP (JD_CU1) = 0
      IDATEP (JD_CE1) = 0
      IDATEP (JD_CU2) = 0
      IDATEP (JD_CE2) = 0
    ENDIF

  ENDIF

ENDIF

!*    2.    PRINT OUT THE DATE
!           ------------------

WRITE(NULOUT,*) 'SUPPDATE WRITES OUT IDATEF = ',IDATEP, size(IDATEP),size(KPPDATE)
call flush(nulout)

KPPDATE = IDATEP (1:SIZE (KPPDATE))

IF (LHOOK) CALL DR_HOOK('SUPPDATE',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE SUPPDATE_READ (CDFILE, KUNIT, KDATEF)

USE SUPPDATE_MOD, ONLY : IDATEF_CFNISH
USE YOMFA, ONLY : LSUPPDATE

CHARACTER (LEN=*),   INTENT (IN)  :: CDFILE
INTEGER (KIND=JPIM), INTENT (IN)  :: KUNIT
INTEGER (KIND=JPIM), INTENT (OUT) :: KDATEF (JD_SIZ)

INTEGER (KIND=JPIM) :: IREP, INBARI
INTEGER (KIND=JPIM), POINTER :: IDATEF (:)
CHARACTER (LEN=16), PARAMETER ::  CLLEC='CADRE LECTURE   '

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUPPDATE:SUPPDATE_READ',0,ZHOOK_HANDLE)

IF (CDFILE == CFNISH) THEN
  IDATEF => IDATEF_CFNISH
ELSE
  CALL ABOR1 ('SUPPDATE:SUPPDATE_READ: UNKNOWN FILE')
ENDIF

IF (LSUPPDATE) THEN
  IDATEF (1:JD_SIZ) = 0
ENDIF

IF (ALL (IDATEF == 0)) THEN

  INBARI=0
  CALL FAITOU(IREP,KUNIT,.TRUE.,CDFILE,'OLD',.TRUE.,.TRUE.,1,1,INBARI,CLLEC)
  CALL LFIMST(IREP,KUNIT,.FALSE.)
  CALL FADIEX(IREP,KUNIT,IDATEF)
  CALL FAIRME(IREP,KUNIT,'UNKNOWN')

ENDIF

KDATEF = IDATEF 

IF (LHOOK) CALL DR_HOOK('SUPPDATE:SUPPDATE_READ',1,ZHOOK_HANDLE)

END SUBROUTINE SUPPDATE_READ

END SUBROUTINE SUPPDATE

