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

SUBROUTINE NEWFA(KNBARP,CDNMCA,KUNTIN,CDFN,LDIMST,KNIMES,CDIDEN,KDATE)

!**** *NEWFA*  - Open a new file Arpege

!     Purpose.
!     --------
!           To open a new file Arpege/Aladin

!**   Interface.
!     ----------
!        *CALL* *NEWFA(...)*

!        Explicit arguments :
!        ------------------
!        INPUT:
!         KNBARP  - Number of fields expected to be written
!         CDNMCA  - file frame
!         KUNTIN  - file unit number
!         CDFN    - file name
!         LDIMST  - .TRUE. to print out figures at close. 
!         KNIMES  - Message level
!         CDIDEN  - file identificator
!         KDATE  - date to write

!        Implicit arguments :
!        --------------------
!        See 'USE MODULE' above.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Calls 'FA' and LFI routines.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-02-12

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KNBARP 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDNMCA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KUNTIN 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDFN 
LOGICAL           ,INTENT(IN)    :: LDIMST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNIMES 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDIDEN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDATE(:) 
CHARACTER (LEN = 7) ::  CLSTTU

INTEGER(KIND=JPIM) :: INBARI, IREP

LOGICAL :: LLERFA, LLNOM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('NEWFA',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*    1.    OPEN A NEW FILE.
!           ---------------

INBARI=0
LLERFA=.TRUE.
LLNOM=.TRUE.
CLSTTU='UNKNOWN'
CALL FAITOU(IREP,KUNTIN,LLNOM,CDFN,CLSTTU,LLERFA,LDIMST,&
 & KNIMES,KNBARP,INBARI,CDNMCA)  
CALL LFIMST(IREP,KUNTIN,.FALSE.)

!*    2.    SET FILE IDENTIFICATOR
!           ----------------------

IF (CDIDEN /= ' ') THEN
  CALL FAUTIF(IREP,KUNTIN,CDIDEN)
ENDIF

!*    3.    WRITE THE DATE
!           --------------

CALL FANDAX(IREP,KUNTIN,KDATE)

IF (LHOOK) CALL DR_HOOK('NEWFA',1,ZHOOK_HANDLE)
END SUBROUTINE NEWFA
