! (C) Copyright 1989- Meteo-France.

SUBROUTINE RDCLIMO(YDGFP_PHYDS,KDOM,CDFILE,KMONTH,YDFPUSERGEO,PFIELD,KGP,KFIELDS,KPTR,KCOD,KFPCLI)

!**** *RDCLIMO*  - READ OUTPUT CLIMATOLOGY

!     PURPOSE.
!     --------
!        Open output climatology files and read the needed fields
!        In case of a file ARPEGE, poles values are not read.

!**   INTERFACE.
!     ----------
!       *CALL* *RDCLIMO*

!        EXPLICIT ARGUMENTS
!        --------------------
!        KDOM   : domain number
!        CDFILE : file name
!        KMONTH : month
!        PFIELD : OUTPUT FIELDS ARRAY
!        KGP    : NUMBER OF OUTPUT POINTS
!        KFIELDS: NUMBER OF OUTPUT FIELDS
!        KPTR   : pointer to the subdomain
!        KCOD   : FIELDS CODES
!        KFPCLI : climatology usage
!                 1 => constant climatology values
!                 2 or 3 => monthly climatology values

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-03-16 IDATE(10)=time of the previous event ; LRFILAF
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 05-02-22 Change interface to SUPPDATE
!      R. El Khatib : 05-03-15 Cleanings
!      R. El Khatib 31-Aug-2012 new E-zone management
!      P.Marguinaud : Refactor using FPSELEZO
!      R. El Khatib 01-Aug-2013 remove LDMASK
!      P.Marguinaud : Change INI3WRFP arguments
!      R. El Khatib 25-Feb-2016 getfplun + local computation of filename
!      R. El Khatib 22-mar-2016 flexible climatology filename with SUFPCLIFNAME
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE FA_MOD   , ONLY : JD_SIZ   ,JD_MON
USE YOMCT1   , ONLY : LRFILAF
USE TYPE_FPDSPHYS, ONLY : FPDSPHY
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE PARFPOS  , ONLY : JPOSPHY

IMPLICIT NONE

TYPE(FPDSPHY),     INTENT(IN)    :: YDGFP_PHYDS(JPOSPHY)
INTEGER(KIND=JPIM),INTENT(IN)    :: KDOM
CHARACTER(LEN=*),  INTENT(IN)    :: CDFILE
INTEGER(KIND=JPIM),INTENT(IN)    :: KMONTH
TYPE (TFPUSERGEO) ,INTENT(IN)    :: YDFPUSERGEO
INTEGER(KIND=JPIM),INTENT(IN)    :: KGP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTR
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFIELD(KGP,KFIELDS) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOD(KFIELDS) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPCLI
! ZREAL    : field read on a file ARPEGE(/ALADIN), 
!            including poles(/possible extension zone)

INTEGER(KIND=JPIM) :: IDATEFA(JD_SIZ)

REAL(KIND=JPRB), ALLOCATABLE :: ZREAL(:)

CHARACTER(LEN=16) :: CLNAME
CHARACTER(LEN=16) :: CLLEC

INTEGER(KIND=JPIM) :: IREP, JFLD, IULFP, IFILE
INTEGER(KIND=JPIM) :: IFXLAT, IFXLON ! number of latitudes & longitudes in file

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "getfplun.intfb.h"
#include "openfpfa.intfb.h"
#include "fpselezo.intfb.h"

!     ------------------------------------------------------------------

!*       1. PREPARATIONS
!           ------------

IF (LHOOK) CALL DR_HOOK('RDCLIMO',0,ZHOOK_HANDLE)

IF (KFPCLI <= 1) THEN
  ! constant values => do not try to test the date
  IFILE=1
ELSE
  ! monthly values => test the date
  IFILE=2
ENDIF

CALL GETFPLUN(KDOM,IULFP)
IF (YDFPUSERGEO%CFPGRID=='GAUSS') THEN
  WRITE(UNIT=CLLEC,FMT='(''FULPOS.ARP.'',I1.1,I2.2,I2.2)') IFILE, KDOM, IULFP
ELSE
  WRITE(UNIT=CLLEC,FMT='(''FULPOS.ALD.'',I1.1,I2.2,I2.2)') IFILE, KDOM, IULFP
ENDIF

!*       2.1 OPEN AND CHECK FILE

IDATEFA(:)=HUGE(1_JPIM)
IDATEFA(JD_MON) = KMONTH

CALL OPENFPFA(YDFPUSERGEO,IULFP,TRIM(CDFILE),CLLEC,IFILE,IDATEFA,LRFILAF,IFXLAT,IFXLON)

!*       2.4 READ FIELDS

IF (YDFPUSERGEO%CFPGRID=='LELAM') THEN
  ALLOCATE(ZREAL(IFXLAT*IFXLON))
ENDIF
DO JFLD=1,KFIELDS
  CLNAME=YDGFP_PHYDS(KCOD(JFLD))%CLNAME
  WRITE(NULOUT,FMT='('' OUTPUT '',A16, '' READ FROM ARPEGE/ALADIN CLIMATE FILE '',A)') CLNAME, TRIM(CDFILE)
  IF (YDFPUSERGEO%CFPGRID=='LELAM') THEN
    CALL FACILE(IREP,IULFP,CLNAME(1:4),1,CLNAME(5:16),ZREAL,.FALSE.)
    CALL FPSELEZO (YDFPUSERGEO, KPTR, KGP, IFXLAT, IFXLON, ZREAL, PFIELD (:,JFLD))
  ELSE
    CALL FACILE(IREP,IULFP,CLNAME(1:4),1,CLNAME(5:16),PFIELD(KPTR,JFLD),.FALSE.)
  ENDIF
ENDDO

!*       2.5 CLOSE FILE

CALL FAIRME(IREP,IULFP,'UNKNOWN')

IF (LHOOK) CALL DR_HOOK('RDCLIMO',1,ZHOOK_HANDLE)
END SUBROUTINE RDCLIMO

