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

SUBROUTINE SUMPIOH(KPROC,KPRIO,KFLG,KFLD,KFLDOFF,LDISIO,KMYPROC,KAVOID,KSHIFT)

!**** *SUMPIOH*  - SET-UP MESSAGE PASSING FOR I/O HANDLINGS

!     PURPOSE.
!     --------
!        To setup the distribution of fields among I/O processors

!**   INTERFACE.
!     ----------
!       *CALL* *SUMPIOH()*

!        EXPLICIT ARGUMENTS
!        ------------------
!           KFLG     : number of fields to distribute
!           KPROC    : total number of processors
!           KPRIO    : number of processors for I/O
!           KFLD     : number of fields for each processor
!           KFLDOFF  : fields offset after distribution among KPRIO procs.
!           LDISIO   : on output, true if KMYPROC is involved in IO
!           KMYPROC  : processor polled for IO
!           KAVOID   : processor to avoid for I/O (optional)

!        IMPLICIT ARGUMENTS
!        ------------------
!         See modules above

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
!        RYAD EL KHATIB *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 99-02-02
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Marguinaud  28-May-2010 Add processor IO probe
!        O.Vignes      27-May-2010 Added optional argument KAVOID
!        P.Marguinaud  11-Sep-2012 Add KSHIFT argument
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)            :: KPROC 
INTEGER(KIND=JPIM),INTENT(IN)            :: KPRIO 
INTEGER(KIND=JPIM),INTENT(IN),  OPTIONAL :: KFLG 
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KFLD(KPROC) 
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KFLDOFF(KPROC) 
LOGICAL,           INTENT(OUT), OPTIONAL :: LDISIO
INTEGER(KIND=JPIM),INTENT(IN),  OPTIONAL :: KMYPROC
INTEGER(KIND=JPIM),INTENT(IN),  OPTIONAL :: KAVOID
INTEGER(KIND=JPIM),INTENT(IN),  OPTIONAL :: KSHIFT

!     IPROC  : processor number each field belongs

#include "abor1.intfb.h"

INTEGER(KIND=JPIM), ALLOCATABLE :: IPROC(:)

INTEGER(KIND=JPIM) :: JF, JROC

INTEGER(KIND=JPIM) :: KCALCPROC
INTEGER(KIND=JPIM) :: KJF
INTEGER(KIND=JPIM) :: ISHIFT


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*      1. SETUP DISTRIBUTION AMONG "I/O" PROCESSORS

KCALCPROC(KJF) = 1+((KPROC*MODULO(KJF-1+ISHIFT,KPRIO))/KPRIO)
!KCALCPROC(KJF) = MOD(KJF-1,KPRIO)+1

IF (LHOOK) CALL DR_HOOK('SUMPIOH',0,ZHOOK_HANDLE)

ISHIFT = 0
IF (PRESENT (KSHIFT)) ISHIFT = KSHIFT

IF (PRESENT(LDISIO) .AND. PRESENT(KMYPROC)) THEN

  LDISIO = .FALSE.
  DO JF=1,KPRIO
    IF(KMYPROC==KCALCPROC(JF))THEN
      LDISIO = .TRUE.
      EXIT
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('SUMPIOH',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IF (.NOT.(PRESENT(KFLG).AND.PRESENT(KFLD).AND.PRESENT(KFLDOFF))) THEN
  CALL ABOR1('SUMPIOH: INCORRECT ARGUMENT SET')
ENDIF

ALLOCATE(IPROC(KFLG))

DO JF=1,KFLG
  IPROC(JF)=KCALCPROC(JF)
ENDDO


DO JROC=1,KPROC
  KFLD(JROC)=COUNT(IPROC == JROC)
ENDDO
KFLDOFF(1)=0
DO JROC=2,KPROC
  KFLDOFF(JROC)=KFLDOFF(JROC-1)+KFLD(JROC-1)
ENDDO

IF (ALLOCATED(IPROC)) DEALLOCATE(IPROC)

IF (LHOOK) CALL DR_HOOK('SUMPIOH',1,ZHOOK_HANDLE)
END SUBROUTINE SUMPIOH


