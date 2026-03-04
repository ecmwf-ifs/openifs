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

SUBROUTINE ADDPGR(CDRECNM,KBLOCK,PDATA,KLEN)

!**** *ADDPGR*  - DDH output : add data to real block for pseudogrib

!     Purpose.
!     --------
!     To prepare for output of DDH in pseudo grib format. Adds data
!     to real part and a record name to the character part.

!**   Interface.
!     ----------
!        *CALL* *ADDPGR(..)*

!        Explicit arguments :     CDRECNM - record name
!        --------------------     KBLOCK  - block number
!                                 PDATA   - data
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

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMPGO   , ONLY : NPTCH    ,NPTRE    ,NPTIND   ,NLENNAM  ,&
 &                    NLENRBL2 ,NLENREA  ,NLENIND  ,NLENCHA  
USE YOMPGOM  , ONLY : CHAPG    ,REAPG    ,NINDPG

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEN 
CHARACTER         ,INTENT(IN)    :: CDRECNM*(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBLOCK 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDATA(KLEN) 
INTEGER(KIND=JPIM) :: ILEN, IPTCHB, IPTCHE, J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1.  ADD RECORD NAME AND CORRESPONDING ENTRY IN INDEX.
!            -------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ADDPGR',0,ZHOOK_HANDLE)
IF(KBLOCK == 1) THEN
  IPTCHB=NPTCH
  IPTCHE=NPTCH+NLENNAM-1
  ILEN=MIN(LEN(CDRECNM),NLENNAM)
  IF(IPTCHE <= NLENCHA) THEN
    CHAPG(IPTCHB:IPTCHE)=CDRECNM(1:ILEN)
    NPTCH=IPTCHE+1
  ELSE
    WRITE(NULOUT,*)' CHARACTER SECTION FOR PSEUDO-GRIB TO SHORT'
    CALL ABOR1('ADDPGR')
  ENDIF

  IF(NPTIND+2 <= NLENIND) THEN
    NPTIND=NPTIND+1
    NINDPG(NPTIND)=NPTRE
    NPTIND=NPTIND+1
    NINDPG(NPTIND)=KLEN
  ELSE
    WRITE(NULOUT,*)' INDEX SECTION FOR PSEUDO-GRIB TO SHORT'
    CALL ABOR1('ADDPGR')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.  ADD DATA INTO REAL PART, UPDATE POINTER.
!            ----------------------------------------

IF(NPTRE+KLEN-1 <= NLENREA) THEN
  DO J=1,KLEN
    REAPG(NPTRE+J-1)=PDATA(J)
  ENDDO
  NPTRE=NPTRE+KLEN
ELSE
  WRITE(NULOUT,*) ' REAL SECTION FOR PSEUDO-GRIB TO SHORT'
  CALL ABOR1('ADDPGR')
ENDIF

!     ------------------------------------------------------------------

!*       3.  ADD UP LENTGTHS FOR BLOCK 2
!            (THIS BLOCK IS ASSUMED TO BE THE SAME LENGTH AS ALL
!             THE FOLLOWING ONES AND THE FIRST ONE EXCEPT FOR
!             HEADER RECORDS)
!            ---------------------------------------------------

IF(KBLOCK == 2) THEN
  NLENRBL2=NLENRBL2+KLEN
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ADDPGR',1,ZHOOK_HANDLE)
END SUBROUTINE ADDPGR

