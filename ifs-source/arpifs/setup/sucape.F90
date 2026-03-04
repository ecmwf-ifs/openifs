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

SUBROUTINE SUCAPE(YDPHY,YDVAB,YDDIMV,KULOUT)

!**** *SUCAPE  * - ROUTINE TO INITIALIZE THE VARIABLES FOR CAPE 
!                  AND CIN COMPUTATION

!     PURPOSE.
!     --------
!        SET DEFAULT VALUES, THEN READS NAMELIST NAMCAPE

!**   INTERFACE.
!     ----------
!        *CALL* *SUCAPE(KULOUT)*

!         EXPLICIT ARGUMENTS :  KULOUT
!         --------------------

!         IMPLICIT ARGUMENTS :
!         --------------------
!            COMMON  YOMCAPE
!            COMMON  YOMPHY
!            COMMON  YOMLUN

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        N.Pristov 03/2001

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03/2001
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!   2010-03-11, J.M. Piriou: change NETAPES default value.
!   2018-09, R. Brozkova: compute NCAPEPSD.
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMVERT  , ONLY : TVAB     ,TVETA    ,TVFE
USE YOMDIMV  , ONLY : TDIMV
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCAPE  , ONLY :  NCAPEITER   ,NETAPES  ,NCAPEPSD ,GCAPERET   ,GCAPEPSD
USE YOMPHY   , ONLY : TPHY
USE YOMLUN   , ONLY :  NULNAM

IMPLICIT NONE

TYPE(TPHY)        ,INTENT(IN) :: YDPHY
TYPE(TVAB)        ,INTENT(IN) :: YDVAB
TYPE(TDIMV)       ,INTENT(IN) :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN) :: KULOUT

INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZVETAF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namcape.nam.h"

!      ----------------------------------------------------------------
!*       1.    SET DEFAULT VALUES.
!              -------------------
IF (LHOOK) CALL DR_HOOK('SUCAPE',0,ZHOOK_HANDLE)
ASSOCIATE(NBITER=>YDPHY%NBITER,NFLEVG=>YDDIMV%NFLEVG)
NCAPEITER=NBITER
NETAPES=1
GCAPERET=0._JPRB
GCAPEPSD=0.7_JPRB

!      ----------------------------------------------------------------
!*       2.    MODIFIES DEFAULT VALUES. READ NAMELIST.
!              ------------------------

CALL POSNAM(NULNAM,'NAMCAPE')
READ(NULNAM,NAMCAPE)

IF (NCAPEITER <= 0 ) CALL ABOR1('SUCAPE:  INVALID VALUE FOR NCAPEITER')
IF (NETAPES <= 0 ) CALL ABOR1('SUCAPE:  INVALID VALUE FOR NETAPES')
IF ((GCAPERET < 0.0_JPRB ) .OR. (GCAPERET > 1.0_JPRB ))&
 & CALL ABOR1('SUCAPE:  INVALID VALUE FOR GCAPERET')  
IF (GCAPEPSD <= 0.0_JPRB ) CALL ABOR1('SUCAPE:  INVALID VALUE FOR GCAPEPSD')

DO JLEV=1,NFLEVG
  ZVETAF=(YDVAB%VALH(JLEV)+YDVAB%VBH(JLEV)+YDVAB%VALH(JLEV-1)+YDVAB%VBH(JLEV-1))*0.5_JPRB
  IF (ZVETAF <= GCAPEPSD) THEN
    NCAPEPSD=JLEV
  ENDIF
ENDDO

!      -----------------------------------------------------------
!*       3.    PRINT FINAL VALUES.
!              -------------------

WRITE(KULOUT,'('' NCAPEITER = '',I2,'' NETAPES = '',I2,'' NCAPEPSD = '',I2, &
 & '' GCAPERET = '',E13.6,'' GCAPEPSD = '',E13.6)')                         &
 & NCAPEITER, NETAPES, NCAPEPSD, GCAPERET, GCAPEPSD
 
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCAPE',1,ZHOOK_HANDLE)
END SUBROUTINE SUCAPE
