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

SUBROUTINE GPGW(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,LDNHDYN,KFLEV,KPROMA,KSTART,KEND,LDGWF,LDGDWI,&
 & POROGL,POROGM,PLNPR,PALPH,PUS,PVS,&
 & PRT,PDVER,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PGWH,PGWF,&
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & LDVFE,PRNHPPI,&
 ! --- OPTIONAL OUTPUT -------------------------------------------------------
 & PGDW&
 & )

! GPGW - Diagnoses "Gw" from the vertical divergence "dver" or from "-G dw".

! Purpose
! -------
!   Diagnoses "Gw" from the vertical divergence "dver" if LDGDWI=F, from "-G dw" if LDGDWI=T.
!   For the finite element vertical discretization (lvertfe=T + lvfe_gw=T), "Gw" is computed only at full levels.
!   For the finite difference vertical discretization, "Gw" is computed at both half and full levels.
!   Calculation is done by vertical integration of the formula:
!    dver = - G/(RT) (pre/prehyd) (d w / d log(prehyd))
!   with the following bottom condition:
!    w_surf = V_surf grad[Phi_s].

!   This routine can be used in a NHEE model: PRNHPPI should be present in this case.

!   This routine can be used in a NHQE or a hydrostatic model: in this case
!   the ratio (prehyd/pre) is equal to 1 and PRNHPPI should not be used.

! Interface
! ---------
!   * INPUT:
!   YDGEOMETRY   : structure containing all geometry.
!   KFLEV        : number of levels.
!   KPROMA       : horizontal dimension.
!   KSTART       : start of work.
!   KEND         : end of work.
!   LDGWF        : calculation of "Gw" at full layers asked for
!                  if finite difference vertical discretization.
!   LDGDWI       : T vs F: input content of PDVER is "-G dw" vs "dver".
!   POROGL       : zonal component of grad[Phi_s].
!   POROGM       : meridian component of grad[Phi_s].
!   PLNPR        : "delta" at full layers (computed in GPXYB).
!   PALPH        : "alpha" at full layers (computed in GPXYB).
!   PUS          : surface U wind.
!   PVS          : surface V wind.
!   PRT          : (RT) at full levels, with the version of R used to define vertical divergence "dver".
!                  R may be Rdry or Rmoist according to definition of vertical divergence "dver".
!   PDVER        : vertical divergence "dver" at full layers.

!   * OUTPUT:
!   PGWH         : G times vertical velocity w at half layers.
!                  (computed only if lvertfe=F)
!   PGWF         : G times vertical velocity w at full layers.

!   * OPTIONAL INPUT:
!   LDVFE        : T if VFE discretisation is used in this routine.
!   PRNHPPI      : (prehyd/pre) at full layers, required for NHEE.

!   * OPTIONAL OUTPUT:
!   PGDW         : contains 'G (dw)'.

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. Yessad, Dec 2004 (after GNHSVD2GW)

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)    :: YDGEOMETRY
LOGICAL           ,INTENT(IN)    :: LDNHDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDGWF
LOGICAL           ,INTENT(IN)    :: LDGDWI
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVER(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWF(KPROMA,KFLEV) 
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LDVFE
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PRNHPPI(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGDW(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
LOGICAL            :: LLVFE
REAL(KIND=JPRB) :: ZGDW(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZGWH(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZIN(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB) :: ZRNHPPI(KPROMA,KFLEV) 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGW',0,ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

!*      0. PRELIMINARY CALCULATIONS
!       ---------------------------

IF( PRESENT(PRNHPPI) ) THEN
  ! * NHEE
  ZRNHPPI(KSTART:KEND,1:KFLEV)=PRNHPPI(KSTART:KEND,1:KFLEV)
ELSE
  ! * NHQE, HYD
  ZRNHPPI(KSTART:KEND,1:KFLEV)=1.0_JPRB
ENDIF

IF( PRESENT(LDVFE) ) THEN
  LLVFE=LDVFE
ELSE
  IF (LDNHDYN) THEN
    ! * NHEE and NHQE
    LLVFE=YDGEOMETRY%YRCVER%LVERTFE.AND.YDGEOMETRY%YRCVER%LVFE_GW
  ELSE
    LLVFE=YDGEOMETRY%YRCVER%LVERTFE
  ENDIF
ENDIF

! -----------------------------------------------------------------------------

!*      1. Case LDGDWI=F: computes "Gw" from "dver"
!       -------------------------------------------

IF (.NOT.LDGDWI) THEN

  ! * Transform "dver" into "-G.dw".
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      ZGDW(JROF,JLEV)=PDVER(JROF,JLEV)* &
       & PRT(JROF,JLEV)*PLNPR(JROF,JLEV)*ZRNHPPI(JROF,JLEV)
    ENDDO
  ENDDO

  ! * Compute "Gw" at the surface (free slip boundary condition)
  DO JROF=KSTART,KEND
      PGWH(JROF,KFLEV)=PUS(JROF)*POROGL(JROF)+PVS(JROF)*POROGM(JROF)
  ENDDO

  IF (LLVFE) THEN

    ! * Store "G dw" at full levels.
    IF (PRESENT(PGDW)) PGDW(KSTART:KEND,1:KFLEV)=-ZGDW(KSTART:KEND,1:KFLEV)

    ! * Compute "Gw" at full levels.

    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        ZGDW(JROF,JLEV)=-ZGDW(JROF,JLEV)*YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
      ENDDO
    ENDDO

    IF (YDGEOMETRY%YRCVER%NVFE_INTBC==0.OR.YDGEOMETRY%YRCVER%NVFE_INTBC==1) THEN
      DO JROF=KSTART,KEND
        ZIN(JROF,0)        = 0.0_JPRB
        ZIN(JROF,1:KFLEV)  = ZGDW(JROF,1:KFLEV)
        ZIN(JROF,KFLEV+1)  = 0.0_JPRB
      ENDDO
      CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'IBOT','11',KPROMA,KSTART,KEND,KFLEV,ZIN,ZGWH)
    ELSE
      DO JROF=KSTART,KEND
        ZIN(JROF,0)        = ZGDW(JROF,1)
        ZIN(JROF,1:KFLEV)  = ZGDW(JROF,1:KFLEV)
        ZIN(JROF,KFLEV+1)  = ZGDW(JROF,KFLEV) ! not applied with INGW
      ENDDO
      CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'INGW','00',KPROMA,KSTART,KEND,KFLEV,ZIN,ZGWH)
    ENDIF

    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        PGWF(JROF,JLEV)=ZGWH(JROF,JLEV-1)+PGWH(JROF,KFLEV)
      ENDDO
    ENDDO

  ELSE

    ! * Compute "Gw" at half levels.

    ! transform -G.dw into Gw
    DO JLEV=KFLEV,1,-1
      DO JROF=KSTART,KEND
        PGWH(JROF,JLEV-1)=PGWH(JROF,JLEV)+ZGDW(JROF,JLEV)
      ENDDO
    ENDDO

    ! * Compute "Gw" at full levels.

    IF (LDGWF) THEN
      ! k.y.: formula pgwf(jlev)=pgwh(jlev)(1-palph(jlev)/plnpr(jlev))
      ! +pgwh(jlev-1)(palph(jlev)/plnpr(jlev)) must be equivalent.
      DO JLEV=1,KFLEV
        DO JROF=KSTART,KEND
          PGWF(JROF,JLEV)=PGWH(JROF,JLEV)+PDVER(JROF,JLEV)* &
           & PRT(JROF,JLEV)*PALPH(JROF,JLEV)*ZRNHPPI(JROF,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ENDIF ! test on LLVFE

ENDIF

! -----------------------------------------------------------------------------

!*      2. Case LDGDWI=T: computes "Gw" from -G.dw 
!       ------------------------------------------

IF (LDGDWI) THEN

  ! * Copy "-G.dw".
  ZGDW(KSTART:KEND,1:KFLEV)=PDVER(KSTART:KEND,1:KFLEV)

  ! * Compute "Gw" at the surface (free slip boundary condition)
  DO JROF=KSTART,KEND
    PGWH(JROF,KFLEV)=PUS(JROF)*POROGL(JROF)+PVS(JROF)*POROGM(JROF)
  ENDDO

  IF (LLVFE) THEN

    ! * Store "G dw" at full levels.
    IF (PRESENT(PGDW)) PGDW(KSTART:KEND,1:KFLEV)=-ZGDW(KSTART:KEND,1:KFLEV)

    ! * Compute "Gw" at full levels.

    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        ZGDW(JROF,JLEV)=-ZGDW(JROF,JLEV)*YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
      ENDDO
    ENDDO

    IF (YDGEOMETRY%YRCVER%NVFE_INTBC==0.OR.YDGEOMETRY%YRCVER%NVFE_INTBC==1) THEN
      DO JROF=KSTART,KEND
        ZIN(JROF,0)        = 0.0_JPRB
        ZIN(JROF,1:KFLEV)  = ZGDW(JROF,1:KFLEV)
        ZIN(JROF,KFLEV+1)  = 0.0_JPRB
      ENDDO
      CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'IBOT','11',KPROMA,KSTART,KEND,KFLEV,ZIN,ZGWH)
    ELSE
      DO JROF=KSTART,KEND
        ZIN(JROF,0)        = ZGDW(JROF,1)
        ZIN(JROF,1:KFLEV)  = ZGDW(JROF,1:KFLEV)
        ZIN(JROF,KFLEV+1)  = ZGDW(JROF,KFLEV)
      ENDDO
      ! Apply RINTBF00, constructed from bottom
      ! Implicit BC not allowed here (set KBC=3)
      CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'ITOP','00',KPROMA,KSTART,KEND,KFLEV,ZIN,ZGWH,KBC=3)
    ENDIF

    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        PGWF(JROF,JLEV)=ZGWH(JROF,JLEV-1)+PGWH(JROF,KFLEV)
      ENDDO
    ENDDO

  ELSE

    ! * Compute "Gw" at half levels.

    ! transform -G.dw into Gw
    DO JLEV=KFLEV,1,-1
      DO JROF=KSTART,KEND
        PGWH(JROF,JLEV-1)=PGWH(JROF,JLEV)+ZGDW(JROF,JLEV)
      ENDDO
    ENDDO

    ! * Compute "Gw" at full levels.

    IF (LDGWF) THEN
      ! k.y.: formula pgwf(jlev)=pgwh(jlev)(1-palph(jlev)/plnpr(jlev))
      ! +pgwh(jlev-1)(palph(jlev)/plnpr(jlev)) must be valid.
      CALL ABOR1(' GPGW: compute "Gw" at full levels: case not coded')
    ENDIF

  ENDIF ! test on LLVFE

ENDIF

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGW',1,ZHOOK_HANDLE)

END SUBROUTINE GPGW

