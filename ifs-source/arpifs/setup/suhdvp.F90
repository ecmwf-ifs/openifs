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

#ifdef RS6K
@PROCESS NOOPTIMIZE
#endif
!OPTION! -O nochg
!OCL  NOUNROLL,NOPREEX,NOEVAL
SUBROUTINE SUHDVP(YDGEOMETRY,YDDYN,YDDYNA,PDILEV,PDILEV_SLD,PDILEVS,LDMESO,KSTDV2)

!**** *SUHDVP*   - Initialize vertical profiles for horizontal diffusion (3D models).

!     Purpose.
!     --------

!         COMPUTES HORIZONTAL DIFFUSION COEFFICIENTS BY USE OF SHAPE FUNCTION
!          (copied from *SUDYN* without any change)

!         Some remarks:
!          - Several options are available for the vertical profile:
!            * NPROFILEHD=1 is designed for ECMWF, but in theory it can
!              be used if LECMWF=F too.
!            * NPROFILEHD=2 and 3 are designed for LECMWF=F, but in theory they can
!              be used if LECMWF=T too.
!            * NPROFILEHD=4 uses and old code (not available in the LAM model),
!              the results of which can be nearly be reproduced by NPROFILEHD=3.

!**   Interface.
!     ----------
!        *CALL* *SUHDVP

!        Explicit arguments :
!        --------------------
!         OUTPUT:
!          PDILEV     : vertical profile applied when no SLHD ('RDI' diffusion)
!          PDILEV_SLD : vertical profile applied when SLHD ('RDI' diffusion)
!          PDILEVS    : vertical profile applied when SLHD ('RDS' diffusion)
!          LDMESO     : T if level nr 1 is in the mesosphere.
!          KSTDV2     : threshold level.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!      K. YESSAD (CNRM/GMAP) after SUHDF.
!      Original : Jan 2012
!      Code was partially in SUHDF until CY38.

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMDYNA  , ONLY : TDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMVERT  , ONLY : VP00

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)  :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT)  :: YDDYN
TYPE(TDYNA)    ,INTENT(INOUT)  :: YDDYNA
REAL(KIND=JPRB),INTENT(OUT)    :: PDILEV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT)    :: PDILEV_SLD(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT)    :: PDILEVS(YDGEOMETRY%YRDIMV%NFLEVG)
LOGICAL,INTENT(OUT)            :: LDMESO
INTEGER(KIND=JPIM),INTENT(OUT) :: KSTDV2

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZPRESH(0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEV_SLD(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEVS(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: ISTDV2, JLEV, JLEVG

LOGICAL :: LLMESO

REAL(KIND=JPRB) :: ZA, ZB, ZAS, ZBS, ZRPROFHDBT, ZRPROFHDTP, ZRAPTP, ZRAP
REAL(KIND=JPRB) :: ZEPS, ZP1DH, ZP1DHS, ZP3DH, ZDILEV, ZENSTDF, ZPRES

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHDVP',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & NPROFILEHD=>YDDYN%NPROFILEHD, &
 & RPROFHDBT=>YDDYN%RPROFHDBT, RPROFHDEX=>YDDYN%RPROFHDEX, &
 & RPROFHDMX=>YDDYN%RPROFHDMX, RPROFHDTP=>YDDYN%RPROFHDTP, &
 & SDRED=>YDDYN%SDRED, SLEVDH=>YDDYN%SLEVDH, SLEVDH3=>YDDYN%SLEVDH3, &
 & SLEVDHS=>YDDYN%SLEVDHS, STPRE=>YDSTA%STPRE, &
 & YDDIMV=>YDGEOMETRY%YRDIMV)
!     ------------------------------------------------------------------

!*       1.    VERTICAL PROFILES.
!              ------------------

LLMESO=0.5_JPRB*YDVAB%VAH(1) <= 10.5_JPRB
IF (NPROFILEHD == 1) THEN
  ! * vertical profile used by default at ECMWF (varies between 1 and 16).
  IF(NFLEVG == 19) THEN
    ZENSTDF=2._JPRB
    ISTDV2=6
  ELSE
    ZENSTDF=1.3_JPRB
    ISTDV2=MIN(9,NFLEVG)
  ENDIF
  IF(LLMESO) THEN
    ! * Remark: why not use STPRE instead of ZPRES?
    ! FOR MODELS WITH TOP ABOVE 0.105hPa PDILEV INCREASES LINEARLY
    ! IN LOG10(PRESSURE) FROM 1. AT 10hPa TO 16. AT 0.1hPa.
    ZPRESH(0)=YDVAB%VAH(0)+YDVAB%VBH(0)*VP00
    DO JLEV=1,NFLEVG
      ZPRESH(JLEV)=YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*VP00
      ZPRES=0.5_JPRB*(ZPRESH(JLEV)+ZPRESH(JLEV-1))
      ZDILEV=1.0_JPRB+7.5_JPRB*(3._JPRB-LOG10(ZPRES))
      Z_PDILEV(JLEV)=MAX(1.0_JPRB,ZDILEV)
      Z_PDILEV_SLD(JLEV)=MAX(MAX(0.0_JPRB,1.0_JPRB-SDRED),ZDILEV)
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
      Z_PDILEV(JLEV)=MIN(16._JPRB,ZENSTDF**MAX(0,ISTDV2-JLEV))
      Z_PDILEV_SLD(JLEV)=MIN(16._JPRB,MAX(0.0_JPRB,ZENSTDF**MAX(0,ISTDV2-JLEV) -SDRED))
    ENDDO
  ENDIF
  Z_PDILEVS(:)=Z_PDILEV(:)
ELSEIF (NPROFILEHD == 2) THEN
  ! * flexible profile using SLEVDH, SLEVDHS, and also RPROFHD...
  ! - first guess "1/prehyd" profile like for NPROFILEHD=3:
  ISTDV2=NFLEVG
  ZEPS=1.E-6_JPRB
  ZP1DH=MAX(2.0_JPRB*ZEPS,SLEVDH*VP00)
  ZP1DHS=MAX(2.0_JPRB*ZEPS,SLEVDHS*VP00)
  DO JLEVG=1,NFLEVG
    IF(STPRE(JLEVG) >= ZP1DH) THEN
      Z_PDILEV(JLEVG)=1.0_JPRB
    ELSE
      Z_PDILEV(JLEVG)=ZP1DH/STPRE(JLEVG)
    ENDIF
  ENDDO
  Z_PDILEV_SLD(:)=Z_PDILEV(:)
  DO JLEVG=1,NFLEVG
    IF(STPRE(JLEVG) >= ZP1DHS) THEN
      Z_PDILEVS(JLEVG)=1.0_JPRB
    ELSE
      Z_PDILEVS(JLEVG)=ZP1DHS/STPRE(JLEVG)
    ENDIF
  ENDDO
  ! - modify these profiles for "prehyd < RPROFHDBT":
  ZRPROFHDTP=MAX(2.0_JPRB*ZEPS,RPROFHDTP)
  ZRPROFHDBT=MAX(2.0_JPRB*ZEPS,RPROFHDBT)
  ZRAPTP=(ZRPROFHDTP/ZRPROFHDBT)**RPROFHDEX
  ZA=((ZRAPTP+1.0_JPRB)*RPROFHDMX-2.0_JPRB*MAX(ZP1DH/ZRPROFHDBT,1.0_JPRB))/(ZRAPTP-1.0_JPRB)
  ZB=2.0_JPRB*MAX(ZP1DH/ZRPROFHDBT,1.0_JPRB)-ZA
  DO JLEVG=1,NFLEVG
    IF (STPRE(JLEVG) < RPROFHDBT .AND. STPRE(JLEVG) > RPROFHDTP) THEN
      ZRAP=(STPRE(JLEVG)/RPROFHDBT)**RPROFHDEX
      Z_PDILEV(JLEVG)=(ZA*ZRAP+ZB)/(ZRAP+1.0_JPRB)
    ELSEIF (STPRE(JLEVG) <= RPROFHDTP) THEN
      Z_PDILEV(JLEVG)=RPROFHDMX
    ENDIF
  ENDDO
  DO JLEVG=1,NFLEVG
    Z_PDILEV_SLD(JLEVG)=MAX(0.0_JPRB,Z_PDILEV(JLEVG)-SDRED)
  ENDDO
  ZAS=((ZRAPTP+1.0_JPRB)*RPROFHDMX-2.0_JPRB*MAX(ZP1DHS/ZRPROFHDBT,1.0_JPRB))/(ZRAPTP-1.0_JPRB)
  ZBS=2.0_JPRB*MAX(ZP1DHS/ZRPROFHDBT,1.0_JPRB)-ZAS
  DO JLEVG=1,NFLEVG
    IF (STPRE(JLEVG) < RPROFHDBT .AND. STPRE(JLEVG) > RPROFHDTP) THEN
      ZRAP=(STPRE(JLEVG)/RPROFHDBT)**RPROFHDEX
      Z_PDILEVS(JLEVG)=(ZAS*ZRAP+ZBS)/(ZRAP+1.0_JPRB)
    ELSEIF (STPRE(JLEVG) <= RPROFHDTP) THEN
      Z_PDILEVS(JLEVG)=RPROFHDMX
    ENDIF
  ENDDO
ELSEIF (NPROFILEHD == 3 .OR. NPROFILEHD == 4) THEN
  ISTDV2=NFLEVG
  ZEPS=1.E-6_JPRB
  ZP1DH=MAX(2.0_JPRB*ZEPS,SLEVDH*VP00)
  ZP1DHS=MAX(2.0_JPRB*ZEPS,SLEVDHS*VP00)
  ZP3DH=MAX(2.0_JPRB*ZEPS,SLEVDH3*VP00)
  DO JLEVG=1,NFLEVG
    IF(STPRE(JLEVG) >= ZP1DH) THEN
      Z_PDILEV(JLEVG)=1.0_JPRB
    ELSEIF (STPRE(JLEVG) < ZP1DH .AND. STPRE(JLEVG) >= ZP3DH) THEN
      Z_PDILEV(JLEVG)=ZP1DH/STPRE(JLEVG)
    ELSE
      Z_PDILEV(JLEVG)=ZP1DH/ZP3DH
    ENDIF
  ENDDO
  DO JLEVG=1,NFLEVG
    Z_PDILEV_SLD(JLEVG)=MAX(0.0_JPRB,Z_PDILEV(JLEVG)-SDRED)
  ENDDO
  DO JLEVG=1,NFLEVG
    IF(STPRE(JLEVG) >= ZP1DHS) THEN
      Z_PDILEVS(JLEVG)=1.0_JPRB
    ELSEIF (STPRE(JLEVG) < ZP1DHS .AND. STPRE(JLEVG) >= ZP3DH) THEN
      Z_PDILEVS(JLEVG)=ZP1DHS/STPRE(JLEVG)
    ELSE
      Z_PDILEVS(JLEVG)=ZP1DHS/ZP3DH
    ENDIF
  ENDDO
ENDIF

PDILEV(:)=Z_PDILEV(:)
PDILEV_SLD(:)=Z_PDILEV_SLD(:)
PDILEVS(:)=Z_PDILEVS(:)
LDMESO=LLMESO
KSTDV2=ISTDV2

!     ------------------------------------------------------------------

!*       2.    PRINTINGS.
!              ----------

WRITE(NULOUT,'(A)')
WRITE(NULOUT,'(A)') ' --- SUHDVP:'
WRITE(NULOUT,'(A)') ' VERTICAL PROFILES FOR HORIZONTAL DIFFUSION'
WRITE(NULOUT,'(A)') ' * PDILEV:'
WRITE(NULOUT,'(5(1X,E23.14))')  (PDILEV(JLEV),JLEV=1,NFLEVG)
IF ((.NOT.YDDYNA%LSLHD_STATIC) .AND. YDDYNA%LSLHD) THEN
  WRITE(NULOUT,'(A)') ' * PDILEV_SLD:'
  WRITE(NULOUT,'(5(1X,E23.14))')  (PDILEV_SLD(JLEV),JLEV=1,NFLEVG)
  WRITE(NULOUT,'(A)') ' * PDILEVS:'
  WRITE(NULOUT,'(5(1X,E23.14))')  (PDILEVS(JLEV),JLEV=1,NFLEVG)
ENDIF
WRITE(NULOUT,'(A)')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHDVP',1,ZHOOK_HANDLE)
END SUBROUTINE SUHDVP
