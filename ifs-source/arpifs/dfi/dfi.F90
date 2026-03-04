! (C) Copyright 1989- Meteo-France.

SUBROUTINE DFI(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDVARBC,YDTCV,YDJOT,YDGOM5,YDODB)

!**** *DFI*  - Controls integration for digital filter initialization
!                on highest level

!     Purpose.
!     --------
!     CONTROLS  DIGITAL FILTER INITIALIZATION ON HIGHEST LEVEL

!**   Interface.
!     ----------
!        *CALL* *DFI

!        Explicit arguments : None
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!     ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!      Gabor Radnoti GMAP/MICECO
!      Original : 92-12-24

!     Modifications.
!     --------------
!      M.Hamrud   01-Oct-2003 CY28 Cleaning
!      K.Yessad : Sep 2010 : simplify organigramme.
!      A. Geer      27-Jul-2015   VarBC is now an object passed by argument, for OOPS
!      F.Suzat    apr2018  YDGOM5 by argument
!     ------------------------------------------------------------------

USE TYPE_MODEL   ,ONLY : MODEL
USE GEOMETRY_MOD ,ONLY : GEOMETRY
USE FIELDS_MOD   ,ONLY : FIELDS
USE MTRAJ_MOD    ,ONLY : MTRAJ
USE PARKIND1     ,ONLY : JPIM     ,JPRB
USE YOMHOOK      ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE JO_TABLE_MOD ,ONLY : JO_TABLE
USE YOMDFI       ,ONLY : NEDFI    ,NSTDFI   ,NSTDFIA  ,RTDFI    ,RTDFIA   ,LADIFH
USE YOMINI       ,ONLY : LINITER
USE VARBC_CLASS  ,ONLY : CLASS_VARBC
USE SUPERGOM_CLASS     , ONLY : CLASS_SUPERGOM
USE TOVSCV_MOD   ,ONLY : TOVSCV
USE DBASE_MOD, ONLY : DBASE
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)   ,INTENT(INOUT) :: YDGEOMETRY  !! INOUT needed for DFI2->SUSPEC
TYPE(FIELDS)     ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)      ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)      ,INTENT(INOUT) :: YDMODEL
TYPE(CLASS_VARBC),INTENT(INOUT) :: YDVARBC
TYPE(CLASS_SUPERGOM), INTENT(INOUT), OPTIONAL :: YDGOM5
CLASS(DBASE), INTENT(INOUT), OPTIONAL         :: YDODB
TYPE(TOVSCV)     ,INTENT(IN),    OPTIONAL :: YDTCV
TYPE(JO_TABLE)   ,INTENT(INOUT), OPTIONAL :: YDJOT
CALL ABOR1('OIFS - DFI should never be called - EXIT')

END SUBROUTINE DFI
