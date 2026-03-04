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

SUBROUTINE SPC2(YDGEOMETRY,YDRIP,YDDYN,CDCONF,KM,KMLOC,YDSPEC)

!**** *SPC2* - SPECTRAL SPACE COMPUTATIONS MULTITASKED ON ZONAL WAVE NUMBER.

!                Caution! Only works for a triangular troncature, i.e.
!                         NASRN(n)=n.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPC2(...)

!        Explicit arguments :
!        --------------------
!          INPUT:
!            CDCONF   - configuration of work.
!            KM       - wavenumber.
!            KMLOC    - wavenumber (DM-local numbering)

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-11-24

!     Modifications.
!     --------------
!      J.PADILLA-BARBOSA 2001-03-12 numerical component of Fourier
!                                   Horizontal Diffusion removed.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad (Dec 2003): cleaning in horizontal diffusion
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden (May 2016): Remove redundant geometry argument
!     ------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0             , ONLY : LRSHW, NCONF, N2DINI
USE YOMCST             , ONLY : ROMEGA
USE YOMRIP             , ONLY : TRIP
USE YOMVODCST          , ONLY : CSTVO, CSTDIV
USE YOMDYN             , ONLY : TDYN
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)          ,INTENT(IN)    :: YDDYN
TYPE(TRIP)          ,INTENT(IN)    :: YDRIP
CHARACTER(LEN=1)    ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KM 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KMLOC
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPDIV(KM:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPVOR(KM:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPDIVG(KM:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPSP(KM:YDGEOMETRY%YRDIM%NSMAX,2)

REAL(KIND=JPRB) :: ZSDIV(KM:YDGEOMETRY%YRDIM%NSMAX,2)

REAL(KIND=JPRB) :: ZSPX(3,KM:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPY(3,KM:YDGEOMETRY%YRDIM%NSMAX,2)

REAL(KIND=JPRB) :: ZALPHA (KM:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZDENIM (KM:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZEPSI  (KM:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZFPLUS (KM:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZFMINUS(KM:YDGEOMETRY%YRDIM%NSMAX+1)

REAL(KIND=JPRB) :: ZY(2,KM:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZX(2,KM:YDGEOMETRY%YRDIM%NSMAX)

INTEGER(KIND=JPIM) :: II, ILL, ILO, IS0, IS02, ISE, JI, JL, JN

REAL(KIND=JPRB) :: ZAL, ZBDT, ZBDT2, ZEM, ZEN, ZF, ZTEMP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "mxptma.h"
#include "mxture.h"
#include "mxturhd.h"
#include "mxturs.h"
#include "simplico.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPC2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, YDLAP=>YDGEOMETRY%YRLAP,YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & LIMPF=>YDDYN%LIMPF, LSIDG=>YDDYN%LSIDG, LSTRHD=>YDDYN%LSTRHD, &
 & RBTS2=>YDDYN%RBTS2, RDHI=>YDDYN%RDHI, RDIDIV=>YDDYN%RDIDIV, &
 & RDISP=>YDDYN%RDISP, RDIVOR=>YDDYN%RDIVOR, SIHEG=>YDDYN%SIHEG, &
 & SIHEG2=>YDDYN%SIHEG2, SIVP=>YDDYN%SIVP, &
 & RSTRET=>YDGEM%RSTRET, &
 & TDT=>YDRIP%TDT, &
 & SCGMAP=>YDSPGEOM%SCGMAP)
!     ------------------------------------------------------------------

!     II=1 if KM=0, II=2 if KM different from 0.
II=1+MIN(1,IABS(KM))
IS0=YDLAP%NSE0L(KMLOC)
IS02=0

IF(.NOT.LRSHW) THEN
  CALL ABOR1(' = SPC2 = ERROR NCONF ')
ENDIF

!     ------------------------------------------------------------------

!*       1.    INITIAL MEMORY TRANSFER.
!              ------------------------

DO JI=1,II
  DO JN=KM,NSMAX
    ISE=YDLAP%NASM0(KM)+(JN-KM)*2+JI-1
    ZSPDIV(JN,JI)=YDSPEC%DIV(1,ISE)
    ZSPVOR(JN,JI)=YDSPEC%VOR(1,ISE)
    ZSPSP(JN,JI)=YDSPEC%SP(ISE)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS FOR SHALLOW WATER MODEL.
!              ------------------------------------------------------------

IF(CDCONF == 'A'.OR.CDCONF == 'I')THEN
  IF (LRSHW) THEN

!*      2.1  Preliminary initialisations.

    ZBDT=RBTS2*TDT
    ZBDT2=(ZBDT*RSTRET)**2

!            Set up helper arrays for implicit Coriolis case

    IF (LIMPF) THEN
      ZEM=REAL(KM,JPRB)
      ZAL=2.0_JPRB*ZBDT*ROMEGA*ZEM
      ILO=KM
      IF (KM == 0) THEN
        ZALPHA(0)=0.0_JPRB
        ZDENIM(0)=0.0_JPRB
        ZEPSI(0)=0.0_JPRB
        ILO=1
      ENDIF
      DO JL=ILO,NSMAX
        ZEN=REAL(JL,JPRB)
        ZALPHA(JL)=ZAL/(ZEN*(ZEN+1.0_JPRB))
        ZDENIM(JL)=1.0_JPRB/(1.0_JPRB+ZALPHA(JL)**2)
        ZEPSI(JL)=SQRT((ZEN*ZEN-ZEM*ZEM)/(4.0_JPRB*ZEN*ZEN-1.0_JPRB))
      ENDDO
      ZALPHA(NSMAX+1)=0.0_JPRB
      ZDENIM(NSMAX+1)=0.0_JPRB

      IF (KM == 0) THEN
        ZFPLUS(0)=0.0_JPRB
        ZFMINUS(0)=0.0_JPRB
      ENDIF
      ZF=2.0_JPRB*ZBDT*ROMEGA
      DO JL=ILO,NSMAX-1
        ZEN=REAL(JL,JPRB)
        ZFPLUS(JL)=ZF*ZEN*ZEPSI(JL+1)/(ZEN+1.0_JPRB)
        ZFMINUS(JL)=ZF*(ZEN+1.0_JPRB)*ZEPSI(JL)/ZEN
      ENDDO
      ZEN=REAL(NSMAX,JPRB)
      ZFPLUS(NSMAX)=0.0_JPRB
      ZFMINUS(NSMAX)=ZF*(ZEN+1.0_JPRB)*ZEPSI(NSMAX)/ZEN
      ZFPLUS(NSMAX+1)=0.0_JPRB
      ZFMINUS(NSMAX+1)=0.0_JPRB
    ENDIF

!*      2.2  Computes right-hand side of Helmholtz equation.

    IF (LSIDG) THEN
      DO JI=1,II
        IF (KM > 0) THEN
          DO JN=KM,NSMAX
            ZSDIV(JN,JI)=YDLAP%RLAPIN(JN)*ZSPDIV(JN,JI)-ZBDT*ZSPSP(JN,JI)
          ENDDO
        ELSEIF (KM == 0) THEN
          DO JN=KM,NSMAX
            ZSDIV(JN,JI)=ZSPDIV(JN,JI)-ZBDT*YDLAP%RLAPDI(JN)*ZSPSP(JN,JI)
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO JI=1,II
        DO JN=KM,NSMAX
          ZSDIV(JN,JI)=ZSPDIV(JN,JI)-ZBDT*YDLAP%RLAPDI(JN)*ZSPSP(JN,JI)
        ENDDO
      ENDDO
    ENDIF

!        For implicit Coriolis case, multiply rhs of vorticity equn by INV[J]
!        Careful : LIMPF assumes all n-values are available !

    IF (LIMPF) THEN
      IF (KM > 0) THEN
        DO JN=KM,NSMAX
          ZTEMP=ZDENIM(JN)*(ZSPVOR(JN,2)+ZALPHA(JN)*ZSPVOR(JN,1))
          ZSPVOR(JN,1)=ZDENIM(JN)*(ZSPVOR(JN,1)-ZALPHA(JN)*ZSPVOR(JN,2))
          ZSPVOR(JN,2)=ZTEMP
        ENDDO
      ENDIF

!        Add [F] * result to rhs of Helmholtz equation

      DO JI=1,II
        DO JN=KM+1,NSMAX
          ZSDIV(JN,JI)=ZSDIV(JN,JI)+ZFMINUS(JN)*ZSPVOR(JN-1,JI)
        ENDDO
        DO JN=KM,NSMAX-1
          ZSDIV(JN,JI)=ZSDIV(JN,JI)+ZFPLUS(JN)*ZSPVOR(JN+1,JI)
        ENDDO
      ENDDO
    ENDIF

!*      2.3  Solve Helmholtz equation

    IF (LSIDG) THEN
!              Inversion of two tridiagonal systems (Helmholtz equation)
!                 --> (DIVprim(t+dt)).
      IF (KM > 0) THEN
        CALL MXTURS(NSMAX+1-KM,1,1,II,&
         & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
         & ZSDIV(KM,1),ZSPDIV(KM,1))  
      ELSEIF (KM == 0) THEN
        CALL MXTURE(NSMAX+1-KM,1,1,II,-2,.TRUE.,&
         & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
         & ZSDIV(KM,1),ZSPDIV(KM,1))  
        CALL MXTURE(NSMAX+1-KM,1,1,II,3,.FALSE.,&
         & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
         & SIHEG2(1,IS02+1,3),ZSDIV(KM,1),ZSPDIV(KM,1))  
      ENDIF
    ELSEIF (LIMPF) THEN

!              Solve complex pentadiagonal system

!           (first have to transpose right-hand side)
      DO JI=1,II
        DO JN=KM,NSMAX
          ZY(JI,JN)=ZSDIV(JN,JI)
        ENDDO
      ENDDO
      CALL SIMPLICO(KM,NSMAX,1,1,ZALPHA(KM),ZDENIM(KM),&
       & ZFPLUS(KM),ZFMINUS(KM),SIVP(1),YDLAP%RLAPDI(0:),&
       & ZBDT2,ZY,ZX)  
!           (now transpose result)
      DO JI=1,II
        DO JN=KM,NSMAX
          ZSPDIV(JN,JI)=ZX(JI,JN)
        ENDDO
      ENDDO
    ELSE
!              Inversion of a diagonal system (Helmholtz equation)
!                 --> (DIVprim(t+dt)).
      DO JI=1,II
        DO JN=KM,NSMAX
          ZSPDIV(JN,JI)=ZSDIV(JN,JI)/(1.0_JPRB-ZBDT2*SIVP(1)*YDLAP%RLAPDI(JN))
        ENDDO
      ENDDO
    ENDIF

    IF (LSIDG) THEN
!           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .
      CALL MXPTMA(NSMAX+1-KM,1,1,II,SCGMAP(IS0+1,1),&
       & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
       & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
       & ZSPDIV(KM,1),ZSPDIVG(KM,1))  
    ELSE
!           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GMBAR**2 * DIVprim(t+dt)) .
      DO JI=1,II
        DO JN=KM,NSMAX
          ZSPDIVG(JN,JI)=ZSPDIV(JN,JI)*RSTRET*RSTRET
        ENDDO
      ENDDO
    ENDIF

!*      2.4  Increment equivalent height.

    DO JI=1,II
      DO JN=KM,NSMAX
        ZSPSP(JN,JI)=ZSPSP(JN,JI)-ZBDT*SIVP(1)*ZSPDIVG(JN,JI)
      ENDDO
    ENDDO

!*       2.5  Increment vorticity

    IF (LIMPF) THEN
      IF (KM == 0) THEN
        DO JN=2,NSMAX
          ZSPVOR(JN,1)=ZSPVOR(JN,1)-ZDENIM(JN)*ZFMINUS(JN)*ZSPDIVG(JN-1,1)
        ENDDO
        DO JN=1,NSMAX-1
          ZSPVOR(JN,1)=ZSPVOR(JN,1)-ZDENIM(JN)*ZFPLUS(JN)*ZSPDIVG(JN+1,1)
        ENDDO
      ELSE
        DO JN=KM+1,NSMAX
          ZSPVOR(JN,1)=ZSPVOR(JN,1)&
           & -ZDENIM(JN)*ZFMINUS(JN)*(ZSPDIVG(JN-1,1)&
           & -ZALPHA(JN)*ZSPDIVG(JN-1,2))  
          ZSPVOR(JN,2)=ZSPVOR(JN,2)&
           & -ZDENIM(JN)*ZFMINUS(JN)*(ZSPDIVG(JN-1,2)&
           & +ZALPHA(JN)*ZSPDIVG(JN-1,1))  
        ENDDO
        DO JN=KM,NSMAX-1
          ZSPVOR(JN,1)=ZSPVOR(JN,1)&
           & -ZDENIM(JN)*ZFPLUS(JN)*(ZSPDIVG(JN+1,1)&
           & -ZALPHA(JN)*ZSPDIVG(JN+1,2))  
          ZSPVOR(JN,2)=ZSPVOR(JN,2)&
           & -ZDENIM(JN)*ZFPLUS(JN)*(ZSPDIVG(JN+1,2)&
           & +ZALPHA(JN)*ZSPDIVG(JN+1,1))  
        ENDDO
      ENDIF
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       5.    HORIZONTAL DIFFUSION ON WORK SPHERE.
!              ------------------------------------

IF (CDCONF == 'A') THEN

  IF (LSTRHD) THEN

!             5.1  Diffusion when stretching.

    DO JI=1,II
      DO JN=KM,NSMAX
        ZSPX(1,JN,JI)=ZSPVOR(JN,JI)*YDLAP%RLAPIN(JN)
        ZSPX(2,JN,JI)=ZSPDIV(JN,JI)*YDLAP%RLAPIN(JN)
        ZSPX(3,JN,JI)=ZSPSP (JN,JI)
      ENDDO
    ENDDO

    ILL=3
    CALL MXTURHD(NSMAX+1-KM,ILL,ILL,-2,.TRUE.,&
     & RDHI(1,IS0+1,1),RDHI(1,IS0+1,2),ZSPX(1,KM,1),ZSPY(1,KM,1))
    CALL MXTURHD(NSMAX+1-KM,ILL,ILL,3,.FALSE.,&
     & RDHI(1,IS0+1,1),RDHI(1,IS0+1,3),ZSPX(1,KM,1),ZSPY(1,KM,1))

    IF (KM > 0) THEN
      CALL MXTURHD(NSMAX+1-KM,ILL,ILL,-2,.TRUE.,&
       & RDHI(1,IS0+1,1),RDHI(1,IS0+1,2),ZSPX(1,KM,2),ZSPY(1,KM,2))
      CALL MXTURHD(NSMAX+1-KM,ILL,ILL,3,.FALSE.,&
       & RDHI(1,IS0+1,1),RDHI(1,IS0+1,3),ZSPX(1,KM,2),ZSPY(1,KM,2))
    ENDIF

    DO JI=1,II
      DO JN=KM,NSMAX
        ZSPVOR(JN,JI)=ZSPY(1,JN,JI)*YDLAP%RLAPDI(JN)
        ZSPDIV(JN,JI)=ZSPY(2,JN,JI)*YDLAP%RLAPDI(JN)
        ZSPSP (JN,JI)=ZSPY(3,JN,JI)
      ENDDO
    ENDDO

  ELSE

!             5.2  Diffusion when no stretching.

    DO JI=1,II
      DO JN=KM,NSMAX
        ZSPVOR(JN,JI)=ZSPVOR(JN,JI)/(1.0_JPRB+TDT*RDIVOR(1,JN))
        ZSPDIV(JN,JI)=ZSPDIV(JN,JI)/(1.0_JPRB+TDT*RDIDIV(1,JN))
        ZSPSP(JN,JI)=ZSPSP(JN,JI)/(1.0_JPRB+TDT*RDISP(JN))
      ENDDO
    ENDDO

  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       6.    FINAL MEMORY TRANSFER.
!              ----------------------

DO JI=1,II
!OCL NOVREC
  DO JN=KM,NSMAX
    ISE=YDLAP%NASM0(KM)+(JN-KM)*2+JI-1
    YDSPEC%DIV(1,ISE)=ZSPDIV(JN,JI)
    YDSPEC%VOR(1,ISE)=ZSPVOR(JN,JI)
    YDSPEC%SP(ISE)=ZSPSP(JN,JI)
  ENDDO
ENDDO

IF(NCONF == 201 .AND. N2DINI == 11) THEN
  DO JI=1,II
!OCL NOVREC
    DO JN=KM,NSMAX
      ISE=YDLAP%NASM0(KM)+(JN-KM)*2+JI-1
      YDSPEC%DIV(1,ISE)=CSTDIV(1,ISE)
      YDSPEC%VOR(1,ISE)=CSTVO(1,ISE)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPC2',1,ZHOOK_HANDLE)
END SUBROUTINE SPC2
