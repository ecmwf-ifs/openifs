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
SUBROUTINE SUHDVPN(YDGEOMETRY,YDDYN,PREXPDH,PSLEVDH1,PSLEVDH2,KNSHD)

!**** *SUHDVPN*   - Initialize vertical profiles N(l) for horizontal diffusion (3D models).

!     Purpose.
!     --------

!      Compute N(l) between NSMAX and NSREFDH.

!**   Interface.
!     ----------
!        *CALL* *SUHDVPN

!     Explicit arguments :
!     --------------------

!      INPUT:
!       PREXPDH    : order of the diffusion.
!       PSLEVDH1   : first threshold for the pressure dependency scaled by VP00.
!       PSLEVDH2   : second threshold for the pressure dependency scaled by VP00 (PSLEVDH2<=PSLEVDH1).

!      OUTPUT:
!       KNSHD     : vertical profile N(l).

!     Implicit arguments :
!     --------------------

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

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMDYN   , ONLY : TDYN
USE YOMVERT  , ONLY : VP00

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)  :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT)  :: YDDYN
REAL(KIND=JPRB),INTENT(IN)     :: PREXPDH
REAL(KIND=JPRB),INTENT(IN)     :: PSLEVDH1
REAL(KIND=JPRB),INTENT(IN)     :: PSLEVDH2
INTEGER(KIND=JPIM),INTENT(OUT) :: KNSHD(YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEVG
INTEGER(KIND=JPIM) :: INSHDI(YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZEPS, ZP1DH, ZP2DH
REAL(KIND=JPRB) :: ZNSMAX,ZNSREFDH,ZNSREFDHR,ZNSMAXR,ZNSR

LOGICAL :: LLVCST ! if T, final profile for KNSHD is vertically constant.

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHDVPN',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
  & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NSREFDH=>YDDYN%NSREFDH, &
 & STPRE=>YDSTA%STPRE)
!     ------------------------------------------------------------------

!*       1.    VERTICAL PROFILES.
!              ------------------

ZEPS=1.E-6_JPRB

ZNSREFDH=REAL(NSREFDH,JPRB)
ZNSMAX=REAL(NSMAX,JPRB)
ZNSREFDHR=ZNSREFDH**PREXPDH
ZNSMAXR=ZNSMAX**PREXPDH

IF (PSLEVDH2 < PSLEVDH1 .AND. NSREFDH < NSMAX) THEN
  ZP1DH=MAX(2.0_JPRB*ZEPS,PSLEVDH1*VP00)
  ZP2DH=MIN(MAX(ZEPS,PSLEVDH2*VP00),ZP1DH-ZEPS)
  DO JLEVG=1,NFLEVG
    IF(STPRE(JLEVG) >= ZP1DH) THEN
      INSHDI(JLEVG)=NSMAX
    ELSEIF(STPRE(JLEVG) >= ZP2DH) THEN
      ZNSR=((ZP1DH*(STPRE(JLEVG)-ZP2DH))/(STPRE(JLEVG)*(ZP1DH-ZP2DH)))*ZNSMAXR&
       & +((ZP2DH*(STPRE(JLEVG)-ZP1DH))/(STPRE(JLEVG)*(ZP2DH-ZP1DH)))*ZNSREFDHR
      INSHDI(JLEVG)=NINT(ZNSR**(1._JPRB/PREXPDH))
    ELSE
      INSHDI(JLEVG)=NSREFDH
    ENDIF
  ENDDO
  LLVCST=.FALSE.
ELSE
  INSHDI(1:NFLEVG)=NSMAX
  LLVCST=.TRUE.
ENDIF

KNSHD(:)=INSHDI(:)

!     ------------------------------------------------------------------

!*       2.    PRINTINGS.
!              ----------

IF (LLVCST) THEN
  WRITE(NULOUT,'(A)') ' KNSHD=NSMAX at all levels '
ELSE
  WRITE(NULOUT,'(A)') ' * KNSHD:'
  WRITE(NULOUT,'(25(1X,I4))') (KNSHD(JLEVG),JLEVG=1,NFLEVG)
ENDIF
WRITE(NULOUT,'(A)')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHDVPN',1,ZHOOK_HANDLE)
END SUBROUTINE SUHDVPN
