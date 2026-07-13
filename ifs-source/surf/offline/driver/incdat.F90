SUBROUTINE INCDAT(KIN,KD,KOU)
! (C) Copyright 1996- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *INCDAT*  - Increments a date by a positive number of days

!     Purpose.
!     --------
!     GIVEN A DATE KIN, INCREMENTS IT BY THE NUMBER OF DAYS, KD,
!     TO OBTAIN KOU

!**   Interface.
!     ----------
!        *CALL* *DATTIM

!        Explicit arguments :
!        --------------------

!     INPUT
!     -----

!     KIN   INITIAL DATE (YYYYMMDD)
!     KD    NUMBER OF DAYS TO INCREMENT

!     OUTPUT
!     ------

!     KOU   FINAL DATE (YYYYMMDD)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------

!     Author.
!     -------
!        Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        Original : 96-03-03

!     ------------------------------------------------------------------


USE YOMLUN1S , ONLY : NULOUT
USE YOMCST   , ONLY : RDAY
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

!* Arguments:
!
INTEGER(KIND=JPIM), INTENT(IN) :: KIN
INTEGER(KIND=JPIM), INTENT(IN) :: KD
INTEGER(KIND=JPIM), INTENT(OUT) :: KOU

!* Local variables:
!
INTEGER(KIND=JPIM) :: J,IYN,IMN,IDN

#include "fcttim.h"


!     ------------------------------------------------------------------

!*       1.   TRIVIAL CASE.
!             -------------

IF (KD.EQ.0) THEN
  KOU=KIN 

!*       2.   NEGATIVE CASE.
!             -------------

ELSEIF (KD.LT.0) THEN
  WRITE(NULOUT,*) " Negative increments not yet allowed in INCDAT"
  WRITE(NULOUT,*) " Increment ",KD
  CALL ABORT

!*       3.   REAL WORK.
!             ----------

ELSE
  IYN=NCCAA(KIN)
  IMN=NMM(KIN)
  IDN=NDD(KIN)
  DO J=1,KD
    IF (IDN.EQ.31.AND.IMN.EQ.12) THEN
      IYN=IYN+1
      IMN=1
      IDN=1
    ELSEIF (IDN.EQ.31.AND.&
     &            (IMN.EQ.1.OR.&
     &             IMN.EQ.3.OR.&
     &             IMN.EQ.5.OR.&
     &             IMN.EQ.7.OR.&
     &             IMN.EQ.8.OR.&
     &             IMN.EQ.10)) THEN
      IMN=IMN+1
      IDN=1
    ELSEIF (IDN.EQ.30.AND.&
     &            (IMN.EQ.2.OR.&
     &             IMN.EQ.4.OR.&
     &             IMN.EQ.6.OR.&
     &             IMN.EQ.9.OR.&
     &             IMN.EQ.11)) THEN
      IMN=IMN+1
      IDN=1
    ELSEIF (IDN.EQ.28.AND.IMN.EQ.2.AND..NOT.&
     &            (MOD(IYN,4).EQ.0.AND.&
     &            (MOD(IYN,400).NE.100.OR.&
     &             MOD(IYN,400).NE.200.OR.&
     &             MOD(IYN,400).NE.300))) THEN
      IMN=IMN+1
      IDN=1
    ELSEIF (IDN.EQ.29.AND.IMN.EQ.2.AND.&
     &            (MOD(IYN,4).EQ.0.AND.&
     &            (MOD(IYN,400).NE.100.OR.&
     &             MOD(IYN,400).NE.200.OR.&
     &             MOD(IYN,400).NE.300))) THEN
      IMN=IMN+1
      IDN=1
    ELSE
      IDN=IDN+1
    ENDIF
  ENDDO
  KOU=10000*IYN+100*IMN+IDN
ENDIF

RETURN
END SUBROUTINE INCDAT
