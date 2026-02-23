! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ACCUM1C(&
  &  YDMODEL, YDDIMV,  PRUISP  ,PRUISS&
  & ,PVDIS   ,PVDISG  ,PUSTRG  ,PVSTRG&
  & ,PFCLL   ,PFCLN&
  & ,PCPR    ,PCPS    ,PSPR    ,PSPS    ,PIRUNOF&
  & ,PISLHF  ,PSWRF0  ,PSWRFN  ,PLWRF0  ,PLWRFN&
  & ,PDIFTQ  ,PDIFTS&
  & ,PFPLCL  ,PFPLCN  ,PFPLSL  ,PFPLSN&
  & ,PFRSO   ,PFRTH&
  &,PSTRTU  ,PSTRTV  )

!**** *ACCUM1C*  - prepares and accumulates diagnostics

!     Purpose.
!     --------
!        prepares and accumulates diagnostic variables

!**   Interface.
!     ----------
!        *CALL* *ACCUM1C

!        Explicit arguments :
!        --------------------


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
!        ECMWF Research Department documentation of the one column model

!     Author.
!     -------
!        Joao Teixeira   *ECMWF*

!     Modifications.
!     --------------
!        Original : 94-01-26
!        J.Teixeira  Jan-95   new output files.
!        M. Ko"hler  6-6-2006 Single Column Model integration within IFS 

!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT3   , ONLY : NSTEP
USE YOMGPD1C , ONLY : VDLSP    ,VDCP     ,VDBLD    ,VDSSHF   ,&
           &VDSLHF   ,VDSSR    ,VDSTR    ,VDTSR    ,VDTTR    ,&
           &VDEWSS   ,VDNSSS   ,VDE      ,VDLGWS   ,VDMGWS   ,&
           &VDGWD    ,VDRO     ,VDIEWSS  ,VDINSSS  ,VDISSHF  ,&
           &VDIE     ,VDCSF    ,VDLSSF

IMPLICIT NONE

!     DUMMY REAL SCALARS
TYPE(MODEL), INTENT(INOUT) :: YDMODEL
TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
REAL(KIND=JPRB) :: PFCLL ,PFCLN
REAL(KIND=JPRB) :: PRUISP,PRUISS
REAL(KIND=JPRB) :: PUSTRG,PVSTRG
REAL(KIND=JPRB) :: PVDIS ,PVDISG
REAL(KIND=JPRB) :: PCPR  ,PCPS  ,PSPR  ,PSPS  ,PIRUNOF
REAL(KIND=JPRB) :: PISLHF,PSWRF0,PSWRFN,PLWRF0,PLWRFN

REAL(KIND=JPRB) :: PFPLCL(0:YDDIMV%NFLEVG), PFPLCN(0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PFPLSL(0:YDDIMV%NFLEVG), PFPLSN(0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PSTRTU(0:YDDIMV%NFLEVG), PSTRTV(0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PFRSO (0:YDDIMV%NFLEVG), PFRTH (0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PDIFTS(0:YDDIMV%NFLEVG), PDIFTQ(0:YDDIMV%NFLEVG)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACCUM1C',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & TSTEP=>YDMODEL%YRML_GCONF%YRRIP%TSTEP)
 
!*       1.   PREPARING FOR WRITE OUT.
!             ------------------------

PCPR    = PFPLCL(NFLEVG)
PCPS    = PFPLCN(NFLEVG)
PSPR    = PFPLSL(NFLEVG)
PSPS    = PFPLSN(NFLEVG)

PSWRF0  = PFRSO(0)
PSWRFN  = PFRSO(NFLEVG)
PLWRF0  = PFRTH(0)
PLWRFN  = PFRTH(NFLEVG)
PIRUNOF = PRUISP+PRUISS

VDIE    = PDIFTQ(NFLEVG)
VDISSHF = PDIFTS(NFLEVG)
PISLHF  = PFCLL+PFCLN
VDIEWSS = PSTRTU(NFLEVG)
VDINSSS = PSTRTV(NFLEVG)

IF(NSTEP == 0) THEN

  VDE    = 0.0_JPRB
  VDSSHF = 0.0_JPRB
  VDSLHF = 0.0_JPRB
  VDEWSS = 0.0_JPRB
  VDNSSS = 0.0_JPRB
  VDSSR  = 0.0_JPRB
  VDSTR  = 0.0_JPRB
  VDTSR  = 0.0_JPRB
  VDTTR  = 0.0_JPRB
  VDLSP  = 0.0_JPRB
  VDCP   = 0.0_JPRB
  VDCSF  = 0.0_JPRB
  VDLSSF = 0.0_JPRB
  VDRO   = 0.0_JPRB
  VDBLD  = 0.0_JPRB
  VDGWD  = 0.0_JPRB
  VDLGWS = 0.0_JPRB
  VDMGWS = 0.0_JPRB

ENDIF

VDE    = VDE    + TSTEP*VDIE
VDSSHF = VDSSHF + TSTEP*VDISSHF
VDSLHF = VDSLHF + TSTEP*PISLHF
VDEWSS = VDEWSS + TSTEP*VDIEWSS
VDNSSS = VDNSSS + TSTEP*VDINSSS
VDSSR  = VDSSR  + TSTEP*PSWRFN
VDSTR  = VDSTR  + TSTEP*PLWRFN
VDTSR  = VDTSR  + TSTEP*PSWRF0
VDTTR  = VDTTR  + TSTEP*PLWRF0
VDLSP  = VDLSP  + TSTEP*PSPR
VDCP   = VDCP   + TSTEP*PCPR
VDCSF  = VDCSF  + TSTEP*PCPS
VDLSSF = VDLSSF + TSTEP*PSPS
VDRO   = VDRO   + TSTEP*PIRUNOF
VDBLD  = VDBLD  + TSTEP*PVDIS
VDGWD  = VDGWD  + TSTEP*PVDISG
VDLGWS = VDLGWS + TSTEP*PUSTRG
VDMGWS = VDMGWS + TSTEP*PVSTRG

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACCUM1C',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------

END SUBROUTINE ACCUM1C
