! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPDOM(CDOM,KXDOM,KFPDOM,CDSTRIN,KDIM,KSTRIN,KNDOM,KDOM)

!**** *SUFPDOM*  - INITIALIZE LIST OF POST-PROCESSING SUBDOMAINS

!     PURPOSE.
!     --------
!        To read the list of subdomains as a character string and convert it
!        into an integer array. In the string, the subdomain names are
!        separated with a ":"
!        In the string is blank, no subdomains are required
!
!       WARNING: The default behaviour has been modified: if the string is blank
!                no domains are required (2010-08-11)

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPDOM*

!        EXPLICIT ARGUMENTS
!        --------------------
!                INPUT :
!           CDOM    : possible names of the subdomains
!           KXDOM   : dimension of CDOM
!           KFPDOM  : maximum number of subdomains
!           CDSTRIN : character array containing the list of subdomains
!           KDIM    : maximum number of character strings
!           KSTRIN  : number of character strings
!                 OUTPUT :
!           KNDOM   : number of subdomains
!           KDOM    : subdomains indexes for each string

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-04-08
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F.Voitus      11-Aug-2010 if the string is blank no domains are required
!                                  (change of default behaviour)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KXDOM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPDOM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTRIN 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDOM(KXDOM) 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDSTRIN(KDIM) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNDOM(KSTRIN) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KDOM(KFPDOM,KSTRIN) 
CHARACTER (LEN = 1) ::  CLSEP
INTEGER(KIND=JPIM) :: ID, IEND, ILOC, IST, JD, JS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ILEND   : string length of CDOM
!     ILEN    : string length of CDSTRIN

!*       1. LOOP ON STRINGS
!           ---------------

IF (LHOOK) CALL DR_HOOK('SUFPDOM',0,ZHOOK_HANDLE)
DO JS=1,KSTRIN
  IF (CDSTRIN(JS) == ' ') THEN
    KNDOM(JS)=0
    DO JD=1,KFPDOM
      KDOM(JD,JS)=1
    ENDDO
  ELSE
    KNDOM(JS)=0
    CLSEP=':'
    IEND=-1
!        loop on substrings
    112 CONTINUE

    IST=IEND+2
    ILOC=INDEX(CDSTRIN(JS)(IST:),CLSEP)
    IF (ILOC == 0) THEN
      IEND=LEN(CDSTRIN(JS))
    ELSE
      IEND=(ILOC-1)+(IST-1)
    ENDIF
!        treatment
    ID=0
    DO JD=1,KFPDOM
      IF (CDSTRIN(JS)(IST:IEND) == CDOM(JD)) ID=JD
    ENDDO
    IF (ID == 0) THEN
      CALL ABOR1(' UNKNOWN SUBDOMAIN ')
    ELSE
      KNDOM(JS)=KNDOM(JS)+1
      KDOM(KNDOM(JS),JS)=ID
    ENDIF

    IF (IEND < LEN(CDSTRIN(JS))) GO TO 112
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('SUFPDOM',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPDOM

