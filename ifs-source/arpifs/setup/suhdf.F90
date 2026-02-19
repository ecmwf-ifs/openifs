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
SUBROUTINE SUHDF(YDGEOMETRY,YDML_GCONF,YDDYN,YDDYNA)

!**** *SUHDF*   - Initialize horizontal diffusion for 3D models (NFLEVG>1).

!     Purpose.
!     --------

!         COMPUTES HORIZONTAL DIFFUSION COEFFICIENTS BY USE OF SHAPE FUNCTION
!          (copied from *SUDYN* without any change)

!         Some remarks:
!          - The code of suhdf.F90 is valid for ARPEGE/IFS, the corresponding
!            code for LAM models is computed in suehdf.F90 (project ALD).
!          - Several options are available for the vertical profile:
!            * NPROFILEHD=1 is designed for ECMWF, but in theory it can be used if LECMWF=F too.
!            * NPROFILEHD=2 and 3 are designed for LECMWF=F, but they must work for LECMWF=T too.
!            * NPROFILEHD=4 saves an old code which can be nearly reproduced by NPROFILEHD=3.
!              Users (especially MF climate team) must switch on NPROFILEHD=3.
!              This obsolescent code must disappear in the future.
!              It is not available in LAM model.
!          - Two options are available for the horizontal function:
!            * LRDISPE_EC=T: uses PDISPEL, PDISPEE, PDISPE and PDISPEX and is designed for ECMWF.
!              But in theory it must work for LECMWF=F too.
!            * LRDISPE_EC=F: uses PDISPE and PDISPVOR and is designed for MF.
!              But in theory it must work for LECMWF=T too.
!          - For GFL variables, the old passive scalars are now considered
!            as obsolete features, and there are not any longer diffused;
!            GFL variables other than "q" and "O3" are not diffused.
!          - The equivalent code for 2D models (shallow-water and vorticity models)
!            is now in suhdf2.F90: LAM version is not coded.

!**   Interface.
!     ----------
!        *CALL* *SUHDF

!        Explicit arguments :
!        --------------------

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
!      G. Radnoti CNRM/GMAP MICECO
!      Original : 92-12-24

!     Modifications.
!     --------------
!      M.Hamrud     : 03-08-01 GFL introduction
!      F. Vana and K. YESSAD: 03-08-25 new vertical profile of diffusion
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.YESSAD (Dec 2003): cleaning horizontal diffusion.
!      F. Vana and K. Yessad: modify diffusion when LSLHD=T.
!      K. Yessad: (jul 2004) bound the increase of hor diffusion
!       in the upper stratosphere.
!      F. Vana : (Nov-2004) avoiding double division of RRDXTAU
!      F. Vana : (Feb-2005) split of the global LSLHD key + cleaning
!      K. Yessad 08-Feb-2005: horizontal diffusion for NH variables (case c=1)
!      K. Yessad (Aug 2005): calculation of RCORDI... moved in SURCORDI.
!      K. Yessad (Aug 2005): calculation of RKRF moved in SURAYFRIC.
!      R. El Khatib 28-Jul-2006 Porting to NEC
!      N. Wedi (Nov 2007): adapt HD for NH model at ECMWF
!      F. Vana  13-Jan-2009 : removed argument + recognized LSLHD_STATIC
!      F. Vana  22-Feb-2011 : SLHD extensnion to LECMWF=.T. code
!      K. Yessad (Jan 2012): shallow-water code moved in SUHDF2, more flexibility in vert profiles.
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMDYNA  , ONLY : TDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMVERT  , ONLY : VP00

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA)    ,INTENT(INOUT) :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
REAL(KIND=JPRB) :: Z_PDILEV_SLD(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEVS(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEV_VOR(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: INSHDI(YDGEOMETRY%YRDIMV%NFLEVG),INSHDS(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IL, ILEV, INSDH, ISTDV2, JL, JLEV, JLEVG, JN

LOGICAL :: LLSLHD_W ,LLSLHD_T ,LLSLHD_Q ,LLSLHD_O3 ,LLSLHD_SVD ,LLSLHD_SPD ,LLSLHD
LOGICAL :: LLGRID, LLMESO

REAL(KIND=JPRB) :: ZEPS, ZP1DH
REAL(KIND=JPRB) :: ZBDH, ZBDHVOR, ZP2DH
REAL(KIND=JPRB) :: ZEFOLD, ZEFOLV, ZLSL
REAL(KIND=JPRB) :: ZPKDH, ZPKDH_REF, ZPKDH_SLHD, ZPKDHVOR, ZPKDHVOR_REF, ZPKDHVOR_SLHD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "suhdvp.intfb.h"
#include "suhdvpn.intfb.h"

!     ------------------------------------------------------------------

! * FUNCTIONS AND FUNCTION ARGUMENTS:
REAL(KIND=JPRB) :: PDISPE
REAL(KIND=JPRB) :: PDISPEE
REAL(KIND=JPRB) :: PDISPEL
REAL(KIND=JPRB) :: PDISPEX
REAL(KIND=JPRB) :: PDISPEX_L91
REAL(KIND=JPRB) :: PDISPVOR
INTEGER(KIND=JPIM) :: KL, KMAX, KN
REAL(KIND=JPRB) :: PEXPDH
REAL(KIND=JPRB) :: PLSL

! * HORIZONTAL DIFFUSION SHAPE FUNCTIONS:
PDISPE(KN,KMAX,PEXPDH)=( MAX ( 0.0_JPRB ,&
 & (SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB))-YDDYN%FRANDH)&
 & /(1.0_JPRB-YDDYN%FRANDH) ) )**PEXPDH  
PDISPVOR(KN,KMAX,PEXPDH)=( MAX ( 0.0_JPRB ,&
 & (SQRT(MAX(0.0_JPRB,REAL(KN*(KN+1)-MIN(REAL(KMAX,JPRB),2.0_JPRB),JPRB)&
 & /REAL(KMAX*(KMAX+1)-MIN(REAL(KMAX,JPRB),2.0_JPRB),JPRB)))&
 & -YDDYN%FRANDH)&
 & /(1.0_JPRB-YDDYN%FRANDH) ) )**PEXPDH  
PDISPEL(KN,KMAX,PEXPDH,KL,PLSL)=&
 & 1.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **PEXPDH      *(Z_PDILEV(KL)-1.0_JPRB)+&
 & 0.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **(PEXPDH*3/2)*(Z_PDILEV(KL)-1.0_JPRB)+&
 & 1.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **(PEXPDH    )*MAX(0.0_JPRB,1.0_JPRB-PLSL*YDDYN%SDRED)
PDISPEE(KN,KMAX,PEXPDH,KL)=&
 & 1.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **PEXPDH      *(Z_PDILEV(KL)-1.0_JPRB)+ 0.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **(PEXPDH*3/2)*(Z_PDILEV(KL)-1.0_JPRB)  
PDISPEX(KN,KMAX,PEXPDH,KL)=&
& 5._JPRB*&
& (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))**2&
& *MAX(9._JPRB-KL,0.0_JPRB)  
PDISPEX_L91(KN,KMAX,PEXPDH,KL)=&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & *MAX(11._JPRB-KL,0.0_JPRB)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHDF',0,ZHOOK_HANDLE)
ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE, &
 & YDSTA=>YDGEOMETRY%YRSTA, YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
 & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM, &
 & YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)

ASSOCIATE(YO3=>YGFL%YO3, YO3_NL=>YGFL%YO3_NL, YQ=>YGFL%YQ, YQ_NL=>YGFL%YQ_NL, &
 & NDLON=>YDGEOMETRY%YRDIM%NDLON, NSMAX=>YDGEOMETRY%YRDIM%NSMAX, &
 & LSPT=>YDDIMF%LSPT, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, NFLEVL=>YDGEOMETRY%YRDIMV%NFLEVL, &
 & FRANDH=>YDDYN%FRANDH, HDIRDIV=>YDDYN%HDIRDIV, HDIRO3=>YDDYN%HDIRO3, &
 & HDIRPD=>YDDYN%HDIRPD, HDIRQ=>YDDYN%HDIRQ, HDIRSP=>YDDYN%HDIRSP, &
 & HDIRT=>YDDYN%HDIRT, HDIRVOR=>YDDYN%HDIRVOR, HRDIRDIV=>YDDYN%HRDIRDIV, &
 & HRDIRO3=>YDDYN%HRDIRO3, HRDIRPD=>YDDYN%HRDIRPD, HRDIRQ=>YDDYN%HRDIRQ, &
 & HRDIRSP=>YDDYN%HRDIRSP, HRDIRT=>YDDYN%HRDIRT, HRDIRVD=>YDDYN%HRDIRVD, &
 & HRDIRVOR=>YDDYN%HRDIRVOR, HRDSRDIV=>YDDYN%HRDSRDIV, HRDSRVD=>YDDYN%HRDSRVD, &
 & HRDSRVOR=>YDDYN%HRDSRVOR, LRDISPE_EC=>YDDYN%LRDISPE_EC, &
 & NPROFILEHD=>YDDYN%NPROFILEHD, NSREFDH=>YDDYN%NSREFDH, RDIDIV=>YDDYN%RDIDIV, &
 & RDIGFL=>YDDYN%RDIGFL, RDIPD=>YDDYN%RDIPD, RDISP=>YDDYN%RDISP, &
 & RDITG=>YDDYN%RDITG, RDIVD=>YDDYN%RDIVD, RDIVOR=>YDDYN%RDIVOR, &
 & RDSDIV=>YDDYN%RDSDIV, RDSVD=>YDDYN%RDSVD, RDSVOR=>YDDYN%RDSVOR, &
 & REXPDH=>YDDYN%REXPDH, REXPDHS=>YDDYN%REXPDHS, SDRED=>YDDYN%SDRED, &
 & SLEVDH=>YDDYN%SLEVDH, SLEVDH1=>YDDYN%SLEVDH1, SLEVDH2=>YDDYN%SLEVDH2, &
 & SLEVDH3=>YDDYN%SLEVDH3, SLEVDHS=>YDDYN%SLEVDHS, SLEVDHS1=>YDDYN%SLEVDHS1, &
 & SLEVDHS2=>YDDYN%SLEVDHS2, &
 & LSLHD_W=>YDDYNA%LSLHD_W,LSLHD_T=>YDDYNA%LSLHD_T,LSLHD_SVD=>YDDYNA%LSLHD_SVD, &
 & LSLHD_SPD=>YDDYNA%LSLHD_SPD,LSLHD_STATIC=>YDDYNA%LSLHD_STATIC,LSLHD=>YDDYNA%LSLHD, &
 & MYLEVS=>YDGEOMETRY%YRMP%MYLEVS, &
 & STPRE=>YDSTA%STPRE)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS AND VERTICAL PROFILES.
!              --------------------------------------------------

! * 1.1: Initialization of SLHD keys

LLSLHD    = (.NOT.LSLHD_STATIC) .AND. LSLHD
LLSLHD_W  = (.NOT.LSLHD_STATIC) .AND. LSLHD_W
LLSLHD_T  = (.NOT.LSLHD_STATIC) .AND. LSLHD_T
LLSLHD_Q  = (.NOT.LSLHD_STATIC) .AND. YQ_NL%LSLHD
LLSLHD_O3 = (.NOT.LSLHD_STATIC) .AND. YO3_NL%LSLHD
LLSLHD_SVD= (.NOT.LSLHD_STATIC) .AND. LSLHD_SVD
LLSLHD_SPD= (.NOT.LSLHD_STATIC) .AND. LSLHD_SPD

! * 1.2: Vertical profiles Z_PDILEV..

CALL SUHDVP(YDGEOMETRY,YDDYN,YDDYNA,Z_PDILEV,Z_PDILEV_SLD,Z_PDILEVS,LLMESO,ISTDV2)

! * 1.3: Vertical profiles INSHD..

WRITE(NULOUT,'(A)')
WRITE(NULOUT,'(A)') ' --- SUHDVPN:'
WRITE(NULOUT,'(A)') ' VERTICAL PROFILES N(l) FOR HORIZONTAL DIFFUSION: HDI'
CALL SUHDVPN(YDGEOMETRY,YDDYN,REXPDH,SLEVDH1,SLEVDH2,INSHDI)
IF (LLSLHD) THEN
  WRITE(NULOUT,'(A)') ' VERTICAL PROFILES N(l) FOR HORIZONTAL DIFFUSION: HDS'
  CALL SUHDVPN(YDGEOMETRY,YDDYN,REXPDHS,SLEVDHS1,SLEVDHS2,INSHDS)
ELSE
  ! not used in this case.
  INSHDS(:)=NSMAX
ENDIF

!     ------------------------------------------------------------------

!*       2.    MAIN HORIZONTAL DIFFUSION.
!              --------------------------

LLGRID=NSMAX > (NDLON+3)/3
WRITE(NULOUT,*) ' SUHDF, LLGRID=',LLGRID,' NSMAX=',NSMAX,' NDLON=',NDLON
IF(ALLOCATED(YDDYN%RDIGFL)) RDIGFL(:,:,:) = 0.0_JPRB

!*       2.1   Computation of RDI[X] and RDS[X] for ECMWF type set-up:

IF(LRDISPE_EC) THEN

  WAVE_LOOP_EC : DO JN=0,NSMAX

    !*       2.1.1 3D variables:

    LEVEL_LOOP1_EC : DO JLEV=1,NFLEVL
      ILEV=MYLEVS(JLEV)
      IF (LLSLHD_W) THEN
        ZLSL=1._JPRB
      ELSE
        ZLSL=0._JPRB
      ENDIF
      IF(HDIRVOR > 1.0_JPRB) THEN
        RDIVOR(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,ZLSL)/HDIRVOR
      ELSE
        RDIVOR(JLEV,JN)=0.0_JPRB
      ENDIF
      IF(HDIRDIV > 1.0_JPRB) THEN
        IF(LLMESO) THEN
          IF(NFLEVG >= 91) THEN
            RDIDIV(JLEV,JN)=(PDISPEL(JN,NSMAX,REXPDH,ILEV,ZLSL)+&
             & PDISPEX_L91(JN,NSMAX,REXPDH,ILEV))/HDIRDIV  
          ELSE
            IF(NSMAX > 200) THEN
              RDIDIV(JLEV,JN)=(PDISPEL(JN,NSMAX,REXPDH,ILEV,ZLSL)+&
               & PDISPEX(JN,NSMAX,REXPDH,ILEV))/HDIRDIV
            ELSE
              RDIDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,ZLSL)/HDIRDIV
            ENDIF  
          ENDIF
        ELSE
          RDIDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,ZLSL)/HDIRDIV
        ENDIF
      ELSE
        IF(LLMESO) THEN
          IF(NFLEVG >= 91) THEN
            RDIDIV(JLEV,JN)=(PDISPEE(JN,NSMAX,REXPDH,ILEV)+&
             & PDISPEX_L91(JN,NSMAX,REXPDH,ILEV))/HDIRVOR  
          ELSE
            IF(NSMAX > 200) THEN
              RDIDIV(JLEV,JN)=(PDISPEE(JN,NSMAX,REXPDH,ILEV)+&
               & PDISPEX(JN,NSMAX,REXPDH,ILEV))/HDIRVOR
            ELSE
              RDIDIV(JLEV,JN)=PDISPEE(JN,NSMAX,REXPDH,ILEV)/HDIRVOR              
            ENDIF
          ENDIF
        ELSE
          RDIDIV(JLEV,JN)=PDISPEE(JN,NSMAX,REXPDH,ILEV)/HDIRVOR
        ENDIF
      ENDIF
      IF (YQ%LSP) THEN
        IF(HDIRQ  > 1.0_JPRB) THEN
          IF (LLSLHD_Q) THEN
            RDIGFL(JLEV,JN,YQ%MPSP)=&
             & Z_PDILEV_SLD(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRQ
          ELSE
            RDIGFL(JLEV,JN,YQ%MPSP)=&
             & Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRQ
          ENDIF
        ENDIF
      ENDIF
      IF (YO3%LSP) THEN
        IF(HDIRO3 > 1.0_JPRB) THEN
          IF (LLSLHD_O3) THEN
            RDIGFL(JLEV,JN,YO3%MPSP)=&
             & Z_PDILEV_SLD(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRO3
          ELSE
            RDIGFL(JLEV,JN,YO3%MPSP)=&
             & Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRO3
          ENDIF
        ENDIF
      ENDIF
      IF(YDDYNA%LNHDYN) THEN
        !ori IF(HDIRPD > 1.0_JPRB) THEN
        !ori   RDIPD(JLEV,JN)=Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRPD
        !ori ELSE
        !ori   RDIPD(JLEV,JN)=0.0_JPRB
        !ori ENDIF
        IF (LLSLHD_SPD) THEN
          ZLSL=1._JPRB
        ELSE
          ZLSL=0._JPRB
        ENDIF
        IF(HDIRPD > 1.0_JPRB) THEN
          RDIPD(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,ZLSL)/HDIRPD
        ELSE
          RDIPD(JLEV,JN)=0.0_JPRB
        ENDIF
        !ori IF(HDIRVD > 1.0_JPRB) THEN
        !ori   RDIVD(JLEV,JN)=Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRVD
        !ori ELSE
        !ori   RDIVD(JLEV,JN)=0.0_JPRB
        !ori ENDIF
        ! This is not very nice with respect to SLHD as it implicitly
        !   assumes LSLHD_W=LSLHD_SVD (which is in a way the most logical
        !   option but still bypassing totally the setup.
        RDIVD(JLEV,JN)=RDIDIV(JLEV,JN)
      ENDIF
      ! Supporting diffusion when SLHD
      !  here only complying with the modern setup
      IF (LLSLHD_W) THEN
        RDSVOR(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,0._JPRB)*HRDSRVOR
        RDSDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,0._JPRB)*HRDSRDIV
      ENDIF
      IF (LLSLHD_SVD) THEN
        RDSVD(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,0._JPRB)*HRDSRVD
      ENDIF
    ENDDO LEVEL_LOOP1_EC

    IF (LLSLHD_T) THEN
      ZLSL=1._JPRB
    ELSE
      ZLSL=0._JPRB
    ENDIF
    DO JLEVG=1,NFLEVG
      IF(HDIRT  > 1.0_JPRB) THEN
        RDITG (JLEVG,JN)=PDISPEL(JN,NSMAX,REXPDH,JLEVG,ZLSL)/HDIRT
      ELSE
        RDITG(JLEVG,JN)=PDISPEE(JN,NSMAX,REXPDH,JLEVG)/HDIRVOR
      ENDIF
    ENDDO

    !*       2.1.2 2D variables:

    IF (HDIRSP > 1.0_JPRB) THEN
      RDISP(JN)=PDISPE(JN,NSMAX,REXPDH)/HDIRSP
    ELSE
      RDISP(JN)=0.0_JPRB
    ENDIF

  ENDDO WAVE_LOOP_EC

ENDIF ! LRDISPE_EC=T

!*       2.2   Computation of RDI[X] for METEO-FRANCE type set-up:

IF (.NOT.LRDISPE_EC) THEN

  IF (NPROFILEHD == 1 .OR. NPROFILEHD == 2 .OR. NPROFILEHD == 3) THEN

    DO JN=0,NSMAX

      !*       2.2.1 3D variables:

      !   * GMV variables (other than T for numerical diffusion), and GFL variables.

      DO JLEV=1,NFLEVL

        !     All "RDI..." arrays dimensioned with NFLEVL.

        ILEV=MYLEVS(JLEV)

        ! * vorticity
        IF (LLSLHD_W) THEN
          RDIVOR(JLEV,JN)=Z_PDILEV_SLD(ILEV)*PDISPVOR(JN,INSHDI(ILEV),REXPDH)*HRDIRVOR
        ELSE
          RDIVOR(JLEV,JN)=Z_PDILEV(ILEV)*PDISPVOR(JN,INSHDI(ILEV),REXPDH)*HRDIRVOR
        ENDIF

        ! * divergence
        IF (LLSLHD_W) THEN
          RDIDIV(JLEV,JN)=Z_PDILEV_SLD(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRDIV
        ELSE
          RDIDIV(JLEV,JN)=Z_PDILEV(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRDIV
        ENDIF

        ! * NH pressure departure variable:
        IF(YDDYNA%LNHDYN) THEN
          IF (LLSLHD_SPD) THEN
            RDIPD(JLEV,JN)=Z_PDILEV_SLD(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRPD
          ELSE
            RDIPD(JLEV,JN)=Z_PDILEV(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRPD
          ENDIF
        ENDIF

        ! * NH vertical divergence variable:
        IF(YDDYNA%LNHDYN) THEN
          IF (LLSLHD_SVD) THEN
            RDIVD(JLEV,JN)=Z_PDILEV_SLD(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRVD
          ELSE
            RDIVD(JLEV,JN)=Z_PDILEV(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRVD
          ENDIF
        ENDIF

        ! -- GFL variables:

        ! * specific humidity:
        IF (YQ%LSP) THEN
          IF (LLSLHD_Q) THEN
            RDIGFL(JLEV,JN,YQ%MPSP)=Z_PDILEV_SLD(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRQ
          ELSE
            RDIGFL(JLEV,JN,YQ%MPSP)=Z_PDILEV(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRQ
          ENDIF
        ENDIF

        ! * ozone:
        IF (YO3%LSP) THEN
          IF (LLSLHD_O3) THEN
            RDIGFL(JLEV,JN,YO3%MPSP)=Z_PDILEV_SLD(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRO3
          ELSE
            RDIGFL(JLEV,JN,YO3%MPSP)=Z_PDILEV(ILEV)*PDISPE(JN,INSHDI(ILEV),REXPDH)*HRDIRO3
          ENDIF
        ENDIF

        ! * other GFL variables: they are currently not diffused.

        ! -- end of GFL variables:

      ENDDO ! JLEV

      ! * Temperature for numerical diffusion:

      IF (LSPT) THEN

        DO JLEVG=1,NFLEVG
          !     All "RDI..." arrays dimensioned with NFLEVG.
          IF (LLSLHD_T) THEN
            RDITG(JLEVG,JN)=Z_PDILEV_SLD(JLEVG)*PDISPE(JN,INSHDI(JLEVG),REXPDH)*HRDIRT
          ELSE
            RDITG(JLEVG,JN)=Z_PDILEV(JLEVG)*PDISPE(JN,INSHDI(JLEVG),REXPDH)*HRDIRT
          ENDIF
        ENDDO ! JLEVG

      ENDIF ! LSPT

      !*       2.2.2 2D variables:

      ! Remark: no diffusion is generally done on this variable,
      ! so the SLHD scheme has not been implemented.

      RDISP(JN)=PDISPE(JN,NSMAX,REXPDH)*HRDIRSP

    ENDDO ! JN

  ELSEIF (NPROFILEHD == 4) THEN

    DO JN=0,NSMAX

      !*       2.2.3 3D variables:

      !   * GMV variables (other than T for numerical diffusion), and GFL variables.

      DO JLEV=1,NFLEVL

        !     All "RDI..." arrays dimensioned with NFLEVL.

        ILEV=MYLEVS(JLEV)

        ZEPS=1.E-6_JPRB

        ZP1DH=MAX(2.0_JPRB*ZEPS,SLEVDH*VP00)
        ZP2DH=MIN(MAX(ZEPS,SLEVDH2*VP00),ZP1DH-ZEPS)
        INSDH=MIN(NSREFDH,NSMAX)
        IF(STPRE(ILEV) >= ZP1DH) THEN
          ! * LOW LAYER (BELOW SLEVDH): THE DIFFUSION COEFFICIENT IS A(JN) AND THE TRUNCATION IS NSMAX.
          Z_PDILEV_VOR(ILEV)=1.0_JPRB
          Z_PDILEV(ILEV)=1.0_JPRB
          ZPKDHVOR_REF=PDISPVOR(JN,NSMAX,REXPDH)
          ZPKDH_REF=PDISPE(JN,NSMAX,REXPDH)
          IF (LLSLHD) THEN
            ZPKDHVOR_SLHD=PDISPVOR(JN,NSMAX,REXPDH)*MAX(0.0_JPRB,1.0_JPRB-SDRED)
            ZPKDH_SLHD=PDISPE(JN,NSMAX,REXPDH)*MAX(0.0_JPRB,1.0_JPRB-SDRED)
          ELSE
            ZPKDHVOR_SLHD=ZPKDHVOR_REF
            ZPKDH_SLHD=ZPKDH_REF
          ENDIF
        ELSEIF(STPRE(ILEV) >= ZP2DH) THEN
          ! * MEDIUM LAYER: THE DIFFUSION COEFFICIENT IS A(JN)/(PRESSURE-B(JN))
          ZBDHVOR=ZP2DH*ZP1DH*&
           & (PDISPVOR(JN,INSDH,REXPDH)-PDISPVOR(JN,NSMAX,REXPDH))&
           & /MAX(ZEPS,ZP1DH*PDISPVOR(JN,INSDH,REXPDH)&
           & -ZP2DH*PDISPVOR(JN,NSMAX,REXPDH))  
          ZBDH=ZP2DH*ZP1DH*&
           & (PDISPE(JN,INSDH,REXPDH)-PDISPE(JN,NSMAX,REXPDH))/&
           & MAX(ZEPS,ZP1DH*PDISPE(JN,INSDH,REXPDH)&
           & -ZP2DH*PDISPE(JN,NSMAX,REXPDH))  
          Z_PDILEV_VOR(ILEV)=(ZP1DH-ZBDHVOR)/(STPRE(ILEV)-ZBDHVOR)
          Z_PDILEV_VOR(ILEV)=MIN(Z_PDILEV_VOR(ILEV),1.0_JPRB/SLEVDH3)
          Z_PDILEV(ILEV)=(ZP1DH-ZBDH)/(STPRE(ILEV)-ZBDH)
          Z_PDILEV(ILEV)=MIN(Z_PDILEV(ILEV),1.0_JPRB/SLEVDH3)
          ZPKDHVOR_REF=PDISPVOR(JN,NSMAX,REXPDH)*Z_PDILEV_VOR(ILEV)
          ZPKDH_REF=PDISPE(JN,NSMAX,REXPDH)*Z_PDILEV(ILEV)
          IF (LLSLHD) THEN
            ZPKDHVOR_SLHD=PDISPVOR(JN,NSMAX,REXPDH)*MAX(0.0_JPRB,Z_PDILEV_VOR(ILEV)-SDRED)
            ZPKDH_SLHD=PDISPE(JN,NSMAX,REXPDH)*MAX(0.0_JPRB,Z_PDILEV(ILEV)-SDRED)
          ELSE
            ZPKDHVOR_SLHD=ZPKDHVOR_REF
            ZPKDH_SLHD=ZPKDH_REF
          ENDIF
        ELSE
          ! * HIGH LAYER (ABOVE SLEVDH2): THE DIFFUSION COEFFICIENT IS A(JN)/PRESSURE AND THE TRUNCATION IS MIN(NSMAX,NSREFDH)
          Z_PDILEV_VOR(ILEV)=ZP1DH/STPRE(ILEV)
          Z_PDILEV_VOR(ILEV)=MIN(Z_PDILEV_VOR(ILEV),1.0_JPRB/SLEVDH3)
          Z_PDILEV(ILEV)=ZP1DH/STPRE(ILEV)
          Z_PDILEV(ILEV)=MIN(Z_PDILEV(ILEV),1.0_JPRB/SLEVDH3)
          ZPKDHVOR_REF=PDISPVOR(JN,INSDH,REXPDH)*Z_PDILEV_VOR(ILEV)
          ZPKDH_REF=PDISPE(JN,INSDH,REXPDH)*Z_PDILEV(ILEV)
          IF (LLSLHD) THEN
            ZPKDHVOR_SLHD=PDISPVOR(JN,INSDH,REXPDH)*MAX(0.0_JPRB,Z_PDILEV_VOR(ILEV)-SDRED)
            ZPKDH_SLHD=PDISPE(JN,INSDH,REXPDH)*MAX(0.0_JPRB,Z_PDILEV(ILEV)-SDRED)
          ELSE
            ZPKDHVOR_SLHD=ZPKDHVOR_REF
            ZPKDH_SLHD=ZPKDH_REF
          ENDIF
        ENDIF

        ! * vorticity and divergence
        IF (LLSLHD_W) THEN
          ZPKDHVOR=ZPKDHVOR_SLHD
          ZPKDH=ZPKDH_SLHD
        ELSE
          ZPKDHVOR=ZPKDHVOR_REF
          ZPKDH=ZPKDH_REF
        ENDIF
        RDIVOR(JLEV,JN)=ZPKDHVOR*HRDIRVOR
        RDIDIV(JLEV,JN)=ZPKDH*HRDIRDIV

        ! * NH pressure departure variable:
        IF (LLSLHD_SPD) THEN
          ZPKDH=ZPKDH_SLHD
        ELSE
          ZPKDH=ZPKDH_REF
        ENDIF
        IF(YDDYNA%LNHDYN) THEN
          RDIPD(JLEV,JN)=ZPKDH*HRDIRPD
        ENDIF

        ! * NH vertical divergence variable:
        IF (LLSLHD_SVD) THEN
          ZPKDH=ZPKDH_SLHD
        ELSE
          ZPKDH=ZPKDH_REF
        ENDIF
        IF(YDDYNA%LNHDYN) THEN
          RDIVD(JLEV,JN)=ZPKDH*HRDIRVD
        ENDIF

        ! -- GFL variables:

        ! * specific humidity:
        IF (LLSLHD_Q) THEN
          ZPKDH=ZPKDH_SLHD
        ELSE
          ZPKDH=ZPKDH_REF
        ENDIF
        IF (YQ%LSP) THEN
          RDIGFL(JLEV,JN,YQ%MPSP)=ZPKDH*HRDIRQ
        ENDIF

        ! * ozone:
        IF (LLSLHD_O3) THEN
          ZPKDH=ZPKDH_SLHD
        ELSE
          ZPKDH=ZPKDH_REF
        ENDIF
        IF (YO3%LSP) THEN
          RDIGFL(JLEV,JN,YO3%MPSP)=ZPKDH*HRDIRO3
        ENDIF

        ! * other GFL variables: they are currently not diffused.

        ! -- end of GFL variables:

      ENDDO ! JLEV

      ! * Temperature for numerical diffusion:

      IF (LSPT) THEN

        DO JLEVG=1,NFLEVG

          !     All "RDI..." arrays dimensioned with NFLEVG.

          ZEPS=1.E-6_JPRB

          ZP1DH=MAX(2.0_JPRB*ZEPS,SLEVDH*VP00)
          ZP2DH=MIN(MAX(ZEPS,SLEVDH2*VP00),ZP1DH-ZEPS)
          INSDH=MIN(NSREFDH,NSMAX)
          IF(STPRE(JLEVG) >= ZP1DH) THEN
            ! * LOW LAYER (BELOW SLEVDH): THE DIFFUSION COEFFICIENT IS A(JN) AND THE TRUNCATION IS NSMAX.
            Z_PDILEV(JLEVG)=1.0_JPRB
            IF (LLSLHD_T) THEN
              ZPKDH=PDISPE(JN,NSMAX,REXPDH)*MAX(0.0_JPRB,1.0_JPRB-SDRED)
            ELSE
              ZPKDH=PDISPE(JN,NSMAX,REXPDH)
            ENDIF
          ELSEIF(STPRE(JLEVG) >= ZP2DH) THEN
            ! * MEDIUM LAYER: THE DIFFUSION COEFFICIENT IS A(JN)/(PRESSURE-B(JN))
            ZBDH=ZP2DH*ZP1DH*&
             & (PDISPE(JN,INSDH,REXPDH)-PDISPE(JN,NSMAX,REXPDH))/&
             & MAX(ZEPS,ZP1DH*PDISPE(JN,INSDH,REXPDH)&
             & -ZP2DH*PDISPE(JN,NSMAX,REXPDH))  
            Z_PDILEV(JLEVG)=(ZP1DH-ZBDH)/(STPRE(JLEVG)-ZBDH)
            Z_PDILEV(JLEVG)=MIN(Z_PDILEV(JLEVG),1.0_JPRB/SLEVDH3)
            IF (LLSLHD_T) THEN
              ZPKDH=PDISPE(JN,NSMAX,REXPDH)*MAX(0.0_JPRB,Z_PDILEV(JLEVG)-SDRED)
            ELSE
              ZPKDH=PDISPE(JN,NSMAX,REXPDH)*Z_PDILEV(JLEVG)
            ENDIF
          ELSE
            ! * HIGH LAYER (ABOVE SLEVDH2): THE DIFFUSION COEFFICIENT IS A(JN)/PRESSURE AND THE TRUNCATION IS MIN(NSMAX,NSREFDH)
            Z_PDILEV(JLEVG)=ZP1DH/STPRE(JLEVG)
            Z_PDILEV(JLEVG)=MIN(Z_PDILEV(JLEVG),1.0_JPRB/SLEVDH3)
            IF (LLSLHD_T) THEN
              ZPKDH=PDISPE(JN,INSDH,REXPDH)*MAX(0.0_JPRB,Z_PDILEV(JLEVG)-SDRED)
            ELSE
              ZPKDH=PDISPE(JN,INSDH,REXPDH)*Z_PDILEV(JLEVG)
            ENDIF
          ENDIF

          RDITG  (JLEVG,JN)=ZPKDH*HRDIRT

        ENDDO ! JLEV

      ENDIF

      !*       2.2.4 2D variables:

      ! Remark: no diffusion is generally done on this variable,
      ! so the SLHD scheme has not been implemented.

      RDISP(JN)=PDISPE(JN,NSMAX,REXPDH)*HRDIRSP

    ENDDO ! JN

  ENDIF ! NPROFILEHD

ENDIF ! .NOT.LRDISPE_EC

!*       2.3   Computation of RDS[X] for VOR,DIV,VD for METEO-FRANCE type set-up:

IF ((.NOT.LRDISPE_EC).AND.(LLSLHD_W.OR.LLSLHD_SVD)) THEN

  IF (NPROFILEHD == 1 .OR. NPROFILEHD == 2 .OR. NPROFILEHD == 3) THEN

    DO JN=0,NSMAX

      DO JLEV=1,NFLEVL

        !     All "RDS..." arrays dimensioned with NFLEVL.

        ILEV=MYLEVS(JLEV)

        ! * vorticity and divergence:
        IF (LLSLHD_W) THEN
          RDSVOR(JLEV,JN)=Z_PDILEVS(ILEV)*PDISPVOR(JN,INSHDS(ILEV),REXPDHS)*HRDSRVOR
          RDSDIV(JLEV,JN)=Z_PDILEVS(ILEV)*PDISPE(JN,INSHDS(ILEV),REXPDHS)*HRDSRDIV
        ENDIF

        ! * NH vertical divergence variable:
        IF(YDDYNA%LNHDYN.AND.LLSLHD_SVD) RDSVD(JLEV,JN)=Z_PDILEVS(ILEV)*PDISPE(JN,INSHDS(ILEV),REXPDHS)*HRDSRVD

      ENDDO ! JLEV

    ENDDO ! JN

  ELSEIF (NPROFILEHD == 4) THEN

    DO JN=0,NSMAX

      DO JLEV=1,NFLEVL

        !     All "RDS..." arrays dimensioned with NFLEVL.

        ILEV=MYLEVS(JLEV)

        ZEPS=1.E-6_JPRB

        ZP1DH=MAX(2.0_JPRB*ZEPS,SLEVDHS*VP00)
        ZP2DH=MIN(MAX(ZEPS,SLEVDHS2*VP00),ZP1DH-ZEPS)
        INSDH=MIN(NSREFDH,NSMAX)
        IF(STPRE(ILEV) >= ZP1DH) THEN
          ! * LOW LAYER (BELOW SLEVDH): THE DIFFUSION COEFFICIENT IS A(JN) AND THE TRUNCATION IS NSMAX.
          Z_PDILEV_VOR(ILEV)=1.0_JPRB
          Z_PDILEV(ILEV)=1.0_JPRB
          ZPKDHVOR=PDISPVOR(JN,NSMAX,REXPDHS)
          ZPKDH=PDISPE(JN,NSMAX,REXPDHS)
        ELSEIF(STPRE(ILEV) >= ZP2DH) THEN
          ! * MEDIUM LAYER: THE DIFFUSION COEFFICIENT IS A(JN)/(PRESSURE-B(JN))
          ZBDHVOR=ZP2DH*ZP1DH*&
           & (PDISPVOR(JN,INSDH,REXPDHS)-PDISPVOR(JN,NSMAX,REXPDHS))&
           & /MAX(ZEPS,ZP1DH*PDISPVOR(JN,INSDH,REXPDHS)&
           & -ZP2DH*PDISPVOR(JN,NSMAX,REXPDHS))  
          ZBDH=ZP2DH*ZP1DH*&
           & (PDISPE(JN,INSDH,REXPDHS)-PDISPE(JN,NSMAX,REXPDHS))/&
           & MAX(ZEPS,ZP1DH*PDISPE(JN,INSDH,REXPDHS)&
           & -ZP2DH*PDISPE(JN,NSMAX,REXPDHS))  
          Z_PDILEV_VOR(ILEV)=(ZP1DH-ZBDHVOR)/(STPRE(ILEV)-ZBDHVOR)
          Z_PDILEV_VOR(ILEV)=MIN(Z_PDILEV_VOR(ILEV),1.0_JPRB/SLEVDH3)
          Z_PDILEV(ILEV)=(ZP1DH-ZBDH)/(STPRE(ILEV)-ZBDH)
          Z_PDILEV(ILEV)=MIN(Z_PDILEV(ILEV),1.0_JPRB/SLEVDH3)
          ZPKDHVOR=PDISPVOR(JN,NSMAX,REXPDHS)*Z_PDILEV_VOR(ILEV)
          ZPKDH=PDISPE(JN,NSMAX,REXPDHS)*Z_PDILEV(ILEV)
        ELSE
          ! * HIGH LAYER (ABOVE SLEVDH2): THE DIFFUSION COEFFICIENT IS A(JN)/PRESSURE AND THE TRUNCATION IS MIN(NSMAX,NSREFDH)
          Z_PDILEV_VOR(ILEV)=ZP1DH/STPRE(ILEV)
          Z_PDILEV_VOR(ILEV)=MIN(Z_PDILEV_VOR(ILEV),1.0_JPRB/SLEVDH3)
          Z_PDILEV(ILEV)=ZP1DH/STPRE(ILEV)
          Z_PDILEV(ILEV)=MIN(Z_PDILEV(ILEV),1.0_JPRB/SLEVDH3)
          ZPKDHVOR=PDISPVOR(JN,INSDH,REXPDHS)*Z_PDILEV_VOR(ILEV)
          ZPKDH=PDISPE(JN,INSDH,REXPDHS)*Z_PDILEV(ILEV)
        ENDIF

        ! * vorticity and divergence:
        IF (LLSLHD_W) THEN
          RDSVOR(JLEV,JN)=ZPKDHVOR*HRDSRVOR
          RDSDIV(JLEV,JN)=ZPKDH   *HRDSRDIV
        ENDIF

        ! * NH vertical divergence variable:
        IF (YDDYNA%LNHDYN.AND.LLSLHD_SVD) RDSVD(JLEV,JN)=ZPKDH*HRDSRVD

      ENDDO ! JLEV

    ENDDO ! JN

  ENDIF ! NPROFILEHD

ENDIF ! .NOT.LRDISPE_EC

!*       2.4   Diagnostic of e-folding time:

WRITE(NULOUT,*) ' '
DO JL=NFLEVL,1,-1
  IL=MYLEVS(JL)
  IF (IL <= ISTDV2) THEN
    IF(RDIVOR(JL,NSMAX) > 1E-7_JPRB) THEN
      ZEFOLV=1.0_JPRB/RDIVOR(JL,NSMAX)/3600._JPRB
    ELSE
      ZEFOLV=0._JPRB
    ENDIF
    IF(RDIDIV(JL,NSMAX) > 1E-7_JPRB) THEN
      ZEFOLD=1.0_JPRB/RDIDIV(JL,NSMAX)/3600._JPRB
    ELSE
      ZEFOLD=0._JPRB
    ENDIF
    WRITE(NULOUT,*) '  E-FOLDING TIME FOR VORTICITY=',ZEFOLV,' AND FOR DIVERGENCE=',ZEFOLD,' HOURS AT LEVEL ',IL  
  ENDIF
ENDDO
WRITE(NULOUT,*) ' '

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHDF',1,ZHOOK_HANDLE)
END SUBROUTINE SUHDF
