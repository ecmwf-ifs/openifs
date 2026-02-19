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

SUBROUTINE WRPGR(CDFNAM,KSTEP)

!**** *WRPGR *  - DDH output : code and write out pseudogrib

!     Purpose.
!     --------
!     Prepare for pseudo-grib coding of DDH data by creating pseudo-grib
!     header and integer section from index and integer data. Code by
!     calling ecCodes. Write out coded data.

!**   Interface.
!     ----------
!        *CALL* *WRPGR(..)*

!        Explicit arguments :     CDFNAM - file name
!        --------------------     KSTEP  - time step

!        Implicit arguments :      None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-03-01
!        Modified : 02-06-10 by D.Dent - improved DDH restart handling
!        M.Hamrud   03-10-01 CY28 Cleaning
!        T.Wilhelmssson 19-11-27 Use ecCodes
!     ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN,   ONLY : NULERR
USE YOMPGO,   ONLY : NPTCH, NPTIN, NPTRE, NPTIND, NLENCHA, NLENIND, NLENINT  
USE YOMPGOM,  ONLY : CHAPG, REAPG, NINDPG, NINTPG
USE IOSTREAM_MIX, ONLY : YGBH
USE GRIB_API_INTERFACE, ONLY : IGRIB_CLONE, IGRIB_SET_VALUE, IGRIB_OPEN_FILE, &
 & IGRIB_WRITE, IGRIB_RELEASE
IMPLICIT NONE

CHARACTER(LEN=14) ,INTENT(IN)    :: CDFNAM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 

INTEGER(KIND=JPIM) :: INTPG(NLENIND+NLENINT), ICHAPG(NLENCHA)
CHARACTER(LEN=9) ,DIMENSION(3)  :: CLFNAM=(/'         ','         ','         '/)
INTEGER(KIND=JPIM), SAVE :: IUDDH(3)
INTEGER(KIND=JPIM) :: IERR, IFILE, ILENIN, J, IHANDLE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRPGR',0,ZHOOK_HANDLE)

!*       1.  PREPARE FOR WRITING OUT PSEUDO-GRIB.
!            ------------------------------------

!*       1.1  GRIB HEADER

CALL IGRIB_CLONE(YGBH%NGRIB_HANDLE_DIAG,IHANDLE)

CALL IGRIB_SET_VALUE(IHANDLE,'indicatorOfUnitOfTimeRange',1)
CALL IGRIB_SET_VALUE(IHANDLE,'startStep',KSTEP)

!*       1.2  ADJUST POINTERS TO POINT TO END OF DATA

NPTCH=NPTCH-1
NPTIN=NPTIN-1
NPTRE=NPTRE-1

!*       1.3  CREATE INTEGER SECTION FROM INDEX AND INTEGER DATA

DO J=1,NPTIND
  INTPG(J)=NINDPG(J)
  IF(INTPG(J) < 0) THEN
    INTPG(J)=INTPG(J)-NPTIND
  ENDIF
ENDDO
DO J=1,NPTIN
  INTPG(J+NPTIND)=NINTPG(J)
ENDDO
ILENIN=NPTIND+NPTIN

!*       1.4  CREATE PSEUDO-GRIB

CALL IGRIB_SET_VALUE(IHANDLE,'numberOfCharacters',NPTCH)
CALL IGRIB_SET_VALUE(IHANDLE,'numberOfFloats',NPTRE)
CALL IGRIB_SET_VALUE(IHANDLE,'numberOfIntegers',ILENIN)
CALL IGRIB_SET_VALUE(IHANDLE,'numberOfLogicals',0)
CALL IGRIB_SET_VALUE(IHANDLE,'floatValues',REAPG(1:NPTRE))
CALL IGRIB_SET_VALUE(IHANDLE,'integerValues',INTPG(1:ILENIN))
! DIAG charValues encoded as one byte integers in ecCodes 
DO J=1,NPTCH
  ICHAPG(J)=IACHAR(CHAPG(J:J))
ENDDO
CALL IGRIB_SET_VALUE(IHANDLE,'charValues',ICHAPG(1:NPTCH))

!     ------------------------------------------------------------------

!*       2.  WRITE OUT PSEUDO-GRIB.
!            ----------------------

!*       2.2  OPEN, WRITE AND CLOSE FILE

DO J=1,3
  IFILE=J
  IF(CDFNAM(1:9) == CLFNAM(J)) THEN
    EXIT
  ELSEIF(CLFNAM(J) == '         ') THEN
    CLFNAM(J)=CDFNAM(1:9)
!     N.B.  use append mode to allow for and existing file when restarting
    CALL IGRIB_OPEN_FILE(IUDDH(IFILE),CLFNAM(IFILE),'a')
    EXIT
  ENDIF
ENDDO

CALL IGRIB_WRITE(IUDDH(IFILE),IHANDLE)
CALL IGRIB_RELEASE(IHANDLE)
! CALL IGRIB_CLOSE(IUDDH)

!     ------------------------------------------------------------------

!*       3.  RELEASE SPACE.
!            --------------

IF (ALLOCATED(NINDPG)) DEALLOCATE(NINDPG)
IF (ALLOCATED(NINTPG)) DEALLOCATE(NINTPG)
IF (ALLOCATED(REAPG )) DEALLOCATE(REAPG )

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRPGR',1,ZHOOK_HANDLE)
END SUBROUTINE WRPGR

