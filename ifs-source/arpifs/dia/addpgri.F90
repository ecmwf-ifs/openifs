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

SUBROUTINE ADDPGRI(CDRECNM,KBLOCK,KDATA,KLEN)

!**** *ADDPGRI* - DDH output : add data to integer block for pseudogrib

!     Purpose.
!     --------
!     To prepare for output of DDH in pseudo grib format. Adds data
!     to integer part and a record name to the character part.

!**   Interface.
!     ----------
!        *CALL* *ADDPGRI(..)*

!        Explicit arguments :     CDRECNM - record name
!        --------------------     KBLOCK  - block number
!                                 KDATA   - data
!                                 KLEN    - number of data words

!        Implicit arguments :      None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-03-01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMPGO   , ONLY : NPTCH    ,NPTIN    ,NPTIND   ,NLENNAM  ,&
 &                    NLENIND  ,NLENINT  ,NLENCHA  
USE YOMPGOM  , ONLY : CHAPG    ,NINDPG   ,NINTPG

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEN 
CHARACTER         ,INTENT(IN)    :: CDRECNM*(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBLOCK 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDATA(KLEN) 
INTEGER(KIND=JPIM) :: ILEN, IPTCHB, IPTCHE, J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1.  ADD RECORD NAME AND CORRESPONDING ENTRY IN INDEX.
!            -------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ADDPGRI',0,ZHOOK_HANDLE)
IF(KBLOCK == 1) THEN
  IPTCHB=NPTCH
  IPTCHE=NPTCH+NLENNAM-1
  ILEN=MIN(LEN(CDRECNM),NLENNAM)
  IF(IPTCHE <= NLENCHA) THEN
    CHAPG(IPTCHB:IPTCHE)=CDRECNM(1:ILEN)
    NPTCH=IPTCHE+1
  ELSE
    WRITE(NULOUT,*)' CHARACTER SECTION FOR PSEUDO-GRIB TO SHORT'
    CALL ABOR1('ADDPGRI')
  ENDIF

  IF(NPTIND+2 <= NLENIND) THEN
    NPTIND=NPTIND+1
    NINDPG(NPTIND)=-NPTIN
    NPTIND=NPTIND+1
    NINDPG(NPTIND)=KLEN
  ELSE
    WRITE(NULOUT,*)' INDEX SECTION FOR PSEUDO-GRIB TO SHORT'
    CALL ABOR1('ADDPGRI')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.  ADD DATA INTO INTEGER PART, UPDATE POINTER.
!            -------------------------------------------

IF(NPTIN+KLEN-1 <= NLENINT) THEN
  DO J=1,KLEN
    NINTPG(NPTIN+J-1)=KDATA(J)
  ENDDO
  NPTIN=NPTIN+KLEN
ELSE
  WRITE(NULOUT,*) ' INTEGER SECTION FOR PSEUDO-GRIB TO SHORT'
  CALL ABOR1('ADDPGRI')
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ADDPGRI',1,ZHOOK_HANDLE)
END SUBROUTINE ADDPGRI

