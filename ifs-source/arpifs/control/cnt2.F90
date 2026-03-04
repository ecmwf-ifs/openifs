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

SUBROUTINE CNT2(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDGOM5,YDODB,YDFPOS)

!**** *CNT2*  -  Controls integration job at level 2.

!     Purpose.
!     --------
!           Controls integration job at level 2.

!**   Interface.
!     ----------
!        *CALL* *CNT2

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    SU2YOM - initialize level 2 commons
!     ----------    CNT3   - controls direct integration on level 3

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      Modified : 03-09-30 C. Fischer remove call to suellt
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      A. Geer      27 Jul 2015  More OOPS cleaning: VARBC by argument
!      P. Lean      17 Aug 2016  OOPS cleaning: ODB by argument
!      P. Lean      22 Mar 2017  OOPS cleaning: Jo-table by argument
!      M. Chrust     3 Jan 2020  OOPS cleaning: Model error config by argument
!     ------------------------------------------------------------------

USE TYPE_MODEL    , ONLY : MODEL
USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE FIELDS_MOD    , ONLY : FIELDS
USE MTRAJ_MOD     , ONLY : MTRAJ
USE PARKIND1      , ONLY : JPRB, JPIM
USE YOMHOOK       , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCT0        , ONLY : LECMWF, NFPOS
USE SUPPDATE_MOD  , ONLY : IDATEF_CFNISH
USE YOMRIP0       , ONLY : NINDAT, NSSSSS
USE YOMOPH0       , ONLY : CFNISH, CFNIGG
USE YOMLUN        , ONLY : NINISH
USE JO_TABLE_MOD  , ONLY : JO_TABLE
USE YOMVAR        , ONLY : LMODERR
USE VARBC_CLASS   , ONLY : CLASS_VARBC
USE TOVSCV_MOD    , ONLY : TOVSCV
USE SUPERGOM_CLASS, ONLY : CLASS_SUPERGOM
USE DBASE_MOD     , ONLY : DBASE
USE FULLPOS       , ONLY : TFPOS
USE YOMMODERR     , ONLY : YGMODERRCONF

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(INOUT) :: YDGEOMETRY  !! INOUT needed for call to CNT3
TYPE(FIELDS)        ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)         ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
TYPE(JO_TABLE)      ,INTENT(INOUT) :: YDJOT
TYPE(CLASS_VARBC)   ,INTENT(INOUT), OPTIONAL :: YDVARBC
TYPE(TOVSCV)        ,INTENT(IN),    OPTIONAL :: YDTCV
TYPE(CLASS_SUPERGOM),INTENT(INOUT), OPTIONAL :: YDGOM5
CLASS(DBASE)        ,INTENT(INOUT), OPTIONAL :: YDODB
TYPE (TFPOS)        ,INTENT(IN)   , OPTIONAL :: YDFPOS

INTEGER (KIND=JPIM) :: ICNT3

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cnt3_wait.intfb.h"
#include "cnt3.intfb.h"
#include "su2yom.intfb.h"
#include "suinimoderr.intfb.h"
#include "filedate.intfb.h"
!     ------------------------------------------------------------------

!*       1.    Initialization.
!              ---------------

!*       1.1   Initialize YOMCT2.


IF (LHOOK) CALL DR_HOOK('CNT2',0,ZHOOK_HANDLE)
CALL SU2YOM(YDMODEL%YRML_GCONF%YRRIP)

IF (LMODERR) CALL SUINIMODERR(YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,YGMODERRCONF)

!     ------------------------------------------------------------------

!*       2.    Call level 3 control routine.
!              -----------------------------

IF (LECMWF.OR.NFPOS==0) THEN

  CALL CNT3(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDGOM5,YDODB,YDFPOS=YDFPOS)

ELSE

  ICNT3 = 0

  DO

! Wait for next file; ICNT3 is updated :
! - ICNT3 =  0 : regular case : single iteration
! - ICNT3 =  N : N-th iteration
! - ICNT3 = -1 : stop now

    CALL CNT3_WAIT (ICNT3)


    IF (ICNT3 > 0) THEN

    ! Reset model date

      IDATEF_CFNISH = 0

      CALL FILEDATE(CDFILESH=CFNISH,CDFILEGG=CFNIGG,KUNIT=NINISH,KUDATE=NINDAT,KUSSSS=NSSSSS)

    ENDIF

! Stop

    IF (ICNT3 < 0) EXIT

    CALL CNT3(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDJOT,YDVARBC,YDTCV,YDGOM5,YDODB,YDFPOS=YDFPOS)

! Single iteration (usual case)

    IF (ICNT3 == 0) EXIT

  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CNT2',1,ZHOOK_HANDLE)

END SUBROUTINE CNT2
