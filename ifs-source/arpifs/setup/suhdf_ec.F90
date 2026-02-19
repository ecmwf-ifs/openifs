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
SUBROUTINE SUHDF_EC(YDGEOMETRY,YDML_GCONF,YDDYN,YDDYNA)

!**** *SUHDF_EC*   - Initialize horizontal diffusion for ECMWF

!     Purpose.
!     --------

!         COMPUTES HORIZONTAL DIFFUSION COEFFICIENTS BY USE OF SHAPE FUNCTION
!         ECMWF version.

!**   Interface.
!     ----------
!        *CALL* *SUHDF_EC

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
!      ECMWF, after SUHDF, 2012

!     Modifications.
!     --------------
!      N. Semane+P. Bechtold  04-10-2012 replace 3600s by RHOUR and add RPLRADI for small planet
!      F. Vana  15-Mar-2013 : Cleaning, simplification and consistent treatment of top with settings bellow
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      I. Polichtchouk (Nov 2020): Introduce sponge that is a level number & start at any height.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RHOUR !, RA
USE YOMVERT  , ONLY : VP00
USE YOMLUN   , ONLY : NULOUT
USE YOMDYNA  , ONLY : TDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMDYNCORE,ONLY : RPLRADI


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA)    ,INTENT(INOUT) :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
REAL(KIND=JPRB) :: ZPRESH(0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEV(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IL, ILEV, ISTDV2, JL, JLEV, JLEVG, JN

LOGICAL :: LLGRID, LLMESO, LLSLHD_W  ,LLSLHD_T  ,LLSLHD_Q  ,LLSLHD_O3 ,&
 & LLSLHD_SVD ,LLSLHD_SPD ,LLSLHD

REAL(KIND=JPRB) :: ZDILEV, ZEFOLD, ZEFOLV, ZEFOLT, ZENSTDF, ZEPS, ZPRES, ZZDIR
REAL(KIND=JPRB) :: ZP1, ZP2, ZMASK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZFAC, ZQN, ZSPDIFF
INTEGER(KIND=JPIM) :: JNLOC, ISHIFT, IGG, ISMAX

!     ------------------------------------------------------------------

! * FUNCTIONS AND FUNCTION ARGUMENTS:
REAL(KIND=JPRB) :: PDISPE
REAL(KIND=JPRB) :: PDISPEE
REAL(KIND=JPRB) :: PDISPEL
REAL(KIND=JPRB) :: PDISPEX
REAL(KIND=JPRB) :: PDISPEX_L91
REAL(KIND=JPRB) :: PDISPEX_VOR_L137
INTEGER(KIND=JPIM) :: KL, KMAX, KN, KBOT
REAL(KIND=JPRB) :: PEXPDH, BSP, KLSP

! * HORIZONTAL DIFFUSION SHAPE FUNCTIONS:
PDISPE(KN,KMAX,PEXPDH)=( MAX ( 0.0_JPRB ,&
 & (SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB))-YDDYN%FRANDH)&
 & /(1.0_JPRB-YDDYN%FRANDH) ) )**PEXPDH  
! modified
PDISPEL(KN,KMAX,PEXPDH,KL,BSP)=&
! switch on/off linear diffusion, also at model top
 & 1.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **PEXPDH      *(Z_PDILEV(KL)-BSP)+&
 & 0.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **(PEXPDH*3/2)*(Z_PDILEV(KL)-BSP)+&
! switch on/off linear diffusion
 & 1.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **(PEXPDH    )
PDISPEE(KN,KMAX,PEXPDH,KL,BSP)=&
 & 1.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **PEXPDH      *(Z_PDILEV(KL)-BSP)+ 0.0_JPRB*&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & **(PEXPDH*3/2)*(Z_PDILEV(KL)-BSP)  
PDISPEX(KN,KMAX,PEXPDH,KL,KLSP)=&
& 5._JPRB*&
& (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))**2&
& *MAX(KLSP-KL,0.0_JPRB)  
PDISPEX_L91(KN,KMAX,PEXPDH,KL,KLSP)=&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB)))))&
 & *MAX(KLSP-KL,0.0_JPRB)
PDISPEX_VOR_L137(KN,KMAX,PEXPDH,KBOT,KL)=&
 & (MAX(0.0_JPRB,(SQRT(REAL(KN*(KN+1),JPRB)/ &
 &                     REAL(KMAX*(KMAX+1),JPRB)))))**PEXPDH &
 & *MAX(REAL(KBOT-KL),0.0_JPRB)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHDF_EC',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA, &
 & YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
 & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM, &
 & YGFL=>YDML_GCONF%YGFL,YDRIP=>YDML_GCONF%YRRIP)

ASSOCIATE(YO3=>YGFL%YO3, YO3_NL=>YGFL%YO3_NL, YQ=>YGFL%YQ, YQ_NL=>YGFL%YQ_NL, &
 & NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NDLON=>YDGEOMETRY%YRDIM%NDLON, NSMAX=>YDGEOMETRY%YRDIM%NSMAX, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, NFLEVL=>YDGEOMETRY%YRDIMV%NFLEVL, &
 & FRANDH=>YDDYN%FRANDH, HDIRDIV=>YDDYN%HDIRDIV, HDIRO3=>YDDYN%HDIRO3, &
 & HDIRPD=>YDDYN%HDIRPD, HDIRQ=>YDDYN%HDIRQ, HDIRSP=>YDDYN%HDIRSP, &
 & HDIRT=>YDDYN%HDIRT, HDIRVOR=>YDDYN%HDIRVOR, LGPSTRESS=>YDDYN%LGPSTRESS, &
 & LHDIFFM=>YDDYN%LHDIFFM, LSPECVIS=>YDDYN%LSPECVIS, NDIFFACT=>YDDYN%NDIFFACT, &
 & LTOP_VOR=>YDDYN%LTOP_VOR, NTOP_VOR_TRUNC=>YDDYN%NTOP_VOR_TRUNC, &
 & NTOP_VOR_BOT=>YDDYN%NTOP_VOR_BOT, KLSPONGE=>YDDYN%KLSPONGE, RBEGSPONGE=>YDDYN%RBEGSPONGE, &
 & RDAMPDIVS=>YDDYN%RDAMPDIVS, RDAMPVDS=>YDDYN%RDAMPVDS, &
 & RDAMPVORS=>YDDYN%RDAMPVORS, RDIDIV=>YDDYN%RDIDIV, RDIGFL=>YDDYN%RDIGFL, &
 & RDIPD=>YDDYN%RDIPD, RDISP=>YDDYN%RDISP, RDITG=>YDDYN%RDITG, &
 & RDIVD=>YDDYN%RDIVD, RDIVOR=>YDDYN%RDIVOR, RDSDIV=>YDDYN%RDSDIV, &
 & RDSVD=>YDDYN%RDSVD, RDSVOR=>YDDYN%RDSVOR, REXPDH=>YDDYN%REXPDH, &
 & REXPDHS=>YDDYN%REXPDHS, SDRED=>YDDYN%SDRED, &
 & LSLHD_W=>YDDYNA%LSLHD_W,LSLHD_T=>YDDYNA%LSLHD_T,LSLHD_SVD=>YDDYNA%LSLHD_SVD, &
 & LSLHD_SPD=>YDDYNA%LSLHD_SPD,LSLHD_STATIC=>YDDYNA%LSLHD_STATIC,LSLHD=>YDDYNA%LSLHD, &
 & MYLEVS=>YDGEOMETRY%YRMP%MYLEVS, &
 & TSTEP=>YDRIP%TSTEP, &
 & STPRE=>YDSTA%STPRE)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------


! * 1.1: Initialization of SLHD keys

LLSLHD    = (.NOT.LSLHD_STATIC) .AND. LSLHD
LLSLHD_W  = (.NOT.LSLHD_STATIC) .AND. LSLHD_W
LLSLHD_T  = (.NOT.LSLHD_STATIC) .AND. LSLHD_T
LLSLHD_Q  = (.NOT.LSLHD_STATIC) .AND. YQ_NL%LSLHD
LLSLHD_O3 = (.NOT.LSLHD_STATIC) .AND. YO3_NL%LSLHD
LLSLHD_SVD= (.NOT.LSLHD_STATIC) .AND. LSLHD_SVD
LLSLHD_SPD= (.NOT.LSLHD_STATIC) .AND. LSLHD_SPD

! To ensure the supporting diffusion is initialized everywhere
IF (LLSLHD_W) THEN
  RDSVOR(:,:)=0._JPRB
  RDSDIV(:,:)=0._JPRB
ENDIF
IF (LLSLHD_SVD) THEN
  RDSVD(:,:)=0._JPRB
ENDIF

! * 1.2: Vertical profiles for ECMWF

  IF(NFLEVG == 19) THEN
    ZENSTDF=2._JPRB
    ISTDV2=6
  ELSE
    ZENSTDF=1.3_JPRB
    ISTDV2=MIN(9,NFLEVG)
    !     ZENSTDF=SQRT(2.)
    !     ISTDV2=MIN(10,NFLEVG)
  ENDIF
  LLMESO=0.5_JPRB*YDVAB%VAH(1) <= 10.5_JPRB
  IF(LLMESO) THEN
    ! * Remark: this piece of code is fully consistent with
    !   lvertfe=F, ndlnpr=0, laprxpk=T
    !   where prehyd(l)=0.5*prehyd(lbar-1)+0.5*prehyd(lbar)
    !   and not completely consistent with the other cases
    !   for lvertfe,ndlnpr,laprxpk.
    WRITE(NULOUT,*) ' SUHDF_EC: ',' Z_PDILEV COMPUTED AS FUNCTION OF PRESSURE: '
    ! FOR MODELS WITH TOP ABOVE 0.105hPa PDILEV INCREASES LINEARLY 
    ! IN LOG10(PRESSURE) FROM 1. AT 10hPa TO 16. AT 0.1hPa.
    ZPRESH(0)=YDVAB%VAH(0)+YDVAB%VBH(0)*VP00
    DO JLEV=1,NFLEVG
      ZPRESH(JLEV)=YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*VP00
      ZPRES=0.5_JPRB*(ZPRESH(JLEV)+ZPRESH(JLEV-1))
      ZDILEV=1.0_JPRB+7.5_JPRB*(3._JPRB-LOG10(ZPRES))
      Z_PDILEV(JLEV)=MAX(RBEGSPONGE,ZDILEV)
    ENDDO
    WRITE(NULOUT,'(A/(5E24.14))') ' Z_PDILEV: ',(Z_PDILEV(JLEV),JLEV=1,NFLEVG)
  ELSE
    RBEGSPONGE=1._JPRB
    WRITE(NULOUT,*) ' SUHDF_EC: RBEGSPONGE RESET TO 1.0 AS LLMESO=F'

    DO JLEV=1,NFLEVG
      Z_PDILEV(JLEV)=MIN(16._JPRB,ZENSTDF**MAX(0,ISTDV2-JLEV))
    ENDDO
    WRITE(NULOUT,'(A/(5E24.14))') ' Z_PDILEV: ',(Z_PDILEV(JLEV),JLEV=1,NFLEVG)
  ENDIF

!     ------------------------------------------------------------------

!*       2.    MAIN HORIZONTAL DIFFUSION.
!              --------------------------

LLGRID=NSMAX > (NDLON+3)/3
WRITE(NULOUT,*) ' SUHDF_EC, LLGRID=',LLGRID,' NSMAX=',NSMAX,' NDLON=',NDLON
IF(ALLOCATED(YDDYN%RDIGFL)) RDIGFL(:,:,:) = 0.0_JPRB

! sponge special for cubic
IF( (NDGLG > NSMAX+1) .AND. YDDYNA%LSLAG ) THEN
  ISMAX=NDGLG-1
ELSE
  ISMAX=NSMAX
ENDIF

!*       2.1   Computation of RDI[X] at ECMWF:

  ! for the case LGRADSP=T and LSPECVIS=F outside the sponge
  ! we want to use less horizontal diffusion
  ZZDIR = REAL(NDIFFACT,JPRB)*TSTEP
  ! Model top (sponge layer) use the previous linear diffusion

IF (NFLEVG > 1) THEN

    DO JN=0,NSMAX

      DO JLEV=1,NFLEVL
        ILEV=MYLEVS(JLEV)
        
        IF( Z_PDILEV(ILEV) > RBEGSPONGE ) THEN
          
          IF(HDIRVOR > 1.0_JPRB) THEN
            RDIVOR(JLEV,JN)=PDISPEL(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRVOR
          ENDIF
          IF(HDIRDIV > 1.0_JPRB) THEN
            IF(LLMESO) THEN
              IF(NFLEVG >= 91) THEN
                RDIDIV(JLEV,JN)=(PDISPEL(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)+&
                 & PDISPEX_L91(JN,ISMAX,REXPDH,ILEV,KLSPONGE))/HDIRDIV  



                IF (LTOP_VOR) THEN
                  RDIVOR(JLEV,JN)=(PDISPEL(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)+&
                   & PDISPEX_VOR_L137(JN,NTOP_VOR_TRUNC,REXPDH,NTOP_VOR_BOT,ILEV))/HDIRDIV
                ENDIF
              ELSE
                IF(REAL(ISMAX,JPRB)/RPLRADI > 200._JPRB) THEN
                  RDIDIV(JLEV,JN)=(PDISPEL(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)+&
                   & PDISPEX(JN,ISMAX,REXPDH,ILEV,KLSPONGE))/HDIRDIV
                ELSE
                  RDIDIV(JLEV,JN)=PDISPEL(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRDIV
                ENDIF
              ENDIF
            ELSE
              RDIDIV(JLEV,JN)=PDISPEL(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRDIV
            ENDIF
          ELSE
            IF(LLMESO) THEN
              IF(NFLEVG >= 91) THEN
                RDIDIV(JLEV,JN)=(PDISPEE(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)+&
                 & PDISPEX_L91(JN,ISMAX,REXPDH,ILEV,KLSPONGE))/HDIRVOR  



                IF (LTOP_VOR) THEN
                  RDIVOR(JLEV,JN)=(PDISPEL(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)+&
                   & PDISPEX_VOR_L137(JN,NTOP_VOR_TRUNC,REXPDH,NTOP_VOR_BOT,ILEV))/HDIRVOR
                ENDIF
              ELSE
                IF(REAL(ISMAX,JPRB)/RPLRADI > 200._JPRB) THEN
                  RDIDIV(JLEV,JN)=(PDISPEE(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)+&
                   & PDISPEX(JN,ISMAX,REXPDH,ILEV,KLSPONGE))/HDIRVOR
                ELSE
                  RDIDIV(JLEV,JN)=PDISPEE(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRVOR              
                ENDIF
              ENDIF
            ELSE
              RDIDIV(JLEV,JN)=PDISPEE(JN,ISMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRVOR
            ENDIF
          ENDIF
          ! Make sure it is not smaller than the setting below the sponge
          IF( LHDIFFM .AND. .NOT.LSPECVIS ) THEN
            IF( NDGLG > NSMAX+1 ) THEN
              ! cubic grid test, reset diffusion to linear grid
              RDIDIV(JLEV,JN)=MAX(RDIDIV(JLEV,JN),PDISPEL(JN,NDGLG-1,REXPDH,ILEV,RBEGSPONGE)/ZZDIR)
              RDIVOR(JLEV,JN)=MAX(RDIVOR(JLEV,JN),PDISPEL(JN,NDGLG-1,REXPDH,ILEV,RBEGSPONGE)/ZZDIR)
            ELSE
              RDIDIV(JLEV,JN)=MAX(RDIDIV(JLEV,JN),PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/ZZDIR)
              RDIVOR(JLEV,JN)=MAX(RDIVOR(JLEV,JN),PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/ZZDIR)
            ENDIF
          ENDIF
        ELSE ! outside sponge 
          ! Note: Whatever is coded here should be followed for the 
          !       security in the sponge layer
          IF( LHDIFFM .AND. .NOT.LSPECVIS ) THEN
            IF( NDGLG > NSMAX+1 ) THEN
              RDIDIV(JLEV,JN)=PDISPEL(JN,NDGLG-1,REXPDH,ILEV,RBEGSPONGE)/ZZDIR
              RDIVOR(JLEV,JN)=PDISPEL(JN,NDGLG-1,REXPDH,ILEV,RBEGSPONGE)/ZZDIR
            ELSE
              RDIDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/ZZDIR
              RDIVOR(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/ZZDIR
            ENDIF
          ELSE
            RDIDIV(JLEV,JN)=0.0_JPRB
            RDIVOR(JLEV,JN)=0.0_JPRB
          ENDIF
        ENDIF

        ! nonhydrostatic
        IF(YDDYNA%LNHDYN) THEN
          RDIPD(JLEV,JN)=0.0_JPRB
          RDIVD(JLEV,JN)=RDIDIV(JLEV,JN)
        ENDIF
          
        ! not used really
        IF (YQ%LSP) THEN
          RDIGFL(JLEV,JN,YQ%MPSP)= 0.0_JPRB
        ENDIF
        IF (YO3%LSP) THEN
          RDIGFL(JLEV,JN,YO3%MPSP)= 0.0_JPRB
        ENDIF
 
      ENDDO

      ! linear diffusion on T        
      DO JLEVG=1,NFLEVG
        IF( Z_PDILEV(JLEVG) > RBEGSPONGE ) THEN
          IF(HDIRT  > 1.0_JPRB) THEN
            RDITG (JLEVG,JN)=PDISPEL(JN,ISMAX,REXPDH,JLEVG,RBEGSPONGE)/HDIRT
          ELSE
            RDITG(JLEVG,JN)=PDISPEE(JN,ISMAX,REXPDH,JLEVG,RBEGSPONGE)/HDIRVOR
          ENDIF
          ! Make sure it is not smaller than the settting bellow the sponge 
          IF (LHDIFFM .AND. .NOT.LSPECVIS) THEN
            IF( NDGLG <= NSMAX+1 ) THEN
              RDITG (JLEVG,JN)=MAX(RDITG (JLEVG,JN),PDISPEL(JN,NSMAX,REXPDH,JLEVG,RBEGSPONGE)/ZZDIR)
            ENDIF
          ENDIF     
        ELSE ! outside sponge 
          ! Note: Whatever is coded here should be followed for the 
          !       security in the sponge layer
          IF( LHDIFFM .AND. .NOT.LSPECVIS) THEN
            IF( NDGLG > NSMAX+1 ) THEN
              ! special, no diffusion on T with cubic grid
              RDITG (JLEVG,JN) = 0.0_JPRB
            ELSE
              RDITG (JLEVG,JN) = PDISPEL(JN,NSMAX,REXPDH,JLEVG,RBEGSPONGE)/ZZDIR
            ENDIF
          ELSE
            RDITG (JLEVG,JN) = 0.0_JPRB
          ENDIF
        ENDIF
      ENDDO

      ! linear diffusion on ps
      RDISP(JN)=0.0_JPRB

    ENDDO
ELSE
  ! shallow water init
  RDIDIV(:,:)=0.0_JPRB
  RDIVOR(:,:)=0.0_JPRB
  RDITG (:,:)=0.0_JPRB
ENDIF
      
IF( LSPECVIS ) THEN
    ! spectral viscosity parameters after Gelb and Gleeson, 2001
    !ZEPS = 2.0_JPRB*RA**3/REAL(NSMAX**3,JPRB)

    ! Pasquelli
    ZEPS = 1.0_JPRB/REAL(2*NSMAX,JPRB)
    ! remove time-scale
    ZEPS = ZEPS/TSTEP

    ! Pasquelli
    IGG = NSMAX/2
    IF( NSMAX >= 1999 ) THEN
      IGG = NSMAX - 800
    ENDIF
    
    ISHIFT=0

    DO JN=0,NSMAX
      
      JNLOC = JN + ISHIFT

      IF( JNLOC > IGG ) THEN

        ZQN=EXP(-0.5_JPRB*REAL(JNLOC-NSMAX,JPRB)**2/REAL(JNLOC-IGG,JPRB)**2)
!        after Gelb and Gleeson, 2001
!        ZFAC=REAL(JNLOC*(JNLOC+1),JPRB)/RA**2
        ZFAC=REAL(JNLOC*(JNLOC+1),JPRB)/REAL(NSMAX*(NSMAX+1),JPRB)
        ZSPDIFF=ZEPS*(ZFAC**2)*(ZQN**2)
      ELSE
        ZSPDIFF=0.0_JPRB
      ENDIF

      DO JLEV=1,NFLEVL
        ! note that shallow water equations only have spectral viscosity now
        RDIVOR(JLEV,JN)=RDIVOR(JLEV,JN)+ZSPDIFF
        RDIDIV(JLEV,JN)=RDIDIV(JLEV,JN)+ZSPDIFF
        
        IF( NFLEVG > 1 ) THEN

          ! nonhydrostatic
          IF(YDDYNA%LNHDYN) THEN
            RDIPD(JLEV,JN)=0.0_JPRB
            RDIVD(JLEV,JN)=RDIDIV(JLEV,JN)
          ENDIF
          
          ! not used really
          IF (YQ%LSP) THEN
            RDIGFL(JLEV,JN,YQ%MPSP)= 0.0_JPRB
          ENDIF
          IF (YO3%LSP) THEN
            RDIGFL(JLEV,JN,YO3%MPSP)= 0.0_JPRB
          ENDIF
        ENDIF
      ENDDO ! end level loop

      ! linear diffusion on T        
      DO JLEVG=1,NFLEVG
        RDITG (JLEVG,JN)=RDITG (JLEVG,JN) + ZSPDIFF
      ENDDO
      
      ! no diffusion on lnsp
      RDISP(JN)=0.0_JPRB

    ENDDO ! end wave loop

ENDIF

! old default configuration
IF( .NOT. LHDIFFM .AND. .NOT.LSPECVIS) THEN

  WAVE_LOOP_EC : DO JN=0,NSMAX

    !*       2.1.1 3D variables:

    LEVEL_LOOP1_EC : DO JLEV=1,NFLEVL
      ILEV=MYLEVS(JLEV)

      ! outside sponge ...
      IF( Z_PDILEV(ILEV) <= RBEGSPONGE ) THEN

        IF (NFLEVG > 1) THEN
          IF(HDIRVOR > 1.0_JPRB) THEN
            RDIVOR(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRVOR
          ELSE
            RDIVOR(JLEV,JN)=0.0_JPRB
          ENDIF
          IF(HDIRDIV > 1.0_JPRB) THEN
            IF(LLMESO) THEN
              IF(NFLEVG >= 91) THEN
                RDIDIV(JLEV,JN)=(PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)+&
                 & PDISPEX_L91(JN,NSMAX,REXPDH,ILEV,KLSPONGE))/HDIRDIV  



              ELSE
                IF(REAL(NSMAX,JPRB)/RPLRADI > 200._JPRB) THEN
                  RDIDIV(JLEV,JN)=(PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)+&
                   & PDISPEX(JN,NSMAX,REXPDH,ILEV,KLSPONGE))/HDIRDIV
                ELSE
                  RDIDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRDIV
                ENDIF
              ENDIF
            ELSE
              RDIDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRDIV
            ENDIF
          ELSE
            IF(LLMESO) THEN
              IF(NFLEVG >= 91) THEN
                RDIDIV(JLEV,JN)=(PDISPEE(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)+&
                 & PDISPEX_L91(JN,NSMAX,REXPDH,ILEV,KLSPONGE))/HDIRVOR  



              ELSE
                IF(REAL(NSMAX,JPRB)/RPLRADI > 200._JPRB) THEN
                  RDIDIV(JLEV,JN)=(PDISPEE(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)+&
                   & PDISPEX(JN,NSMAX,REXPDH,ILEV,KLSPONGE))/HDIRVOR
                ELSE
                  RDIDIV(JLEV,JN)=PDISPEE(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRVOR              
                ENDIF
              ENDIF
            ELSE
              RDIDIV(JLEV,JN)=PDISPEE(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRVOR
            ENDIF
          ENDIF
          IF (YQ%LSP) THEN
            IF(HDIRQ > 1.0_JPRB) THEN
              RDIGFL(JLEV,JN,YQ%MPSP)=&
               & Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRQ
            ENDIF
          ENDIF
          IF (YO3%LSP) THEN
            IF(HDIRO3 > 1.0_JPRB) THEN
              RDIGFL(JLEV,JN,YO3%MPSP)=&
               & Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRO3
            ENDIF
          ENDIF
          IF(YDDYNA%LNHDYN) THEN
            IF(HDIRPD > 1.0_JPRB) THEN
              RDIPD(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDH,ILEV,RBEGSPONGE)/HDIRPD
            ELSE
              RDIPD(JLEV,JN)=0.0_JPRB
            ENDIF
            RDIVD(JLEV,JN)=RDIDIV(JLEV,JN)
          ENDIF
        ELSE
          ! --- Code for shallow-water equations:
          !     SLHD not available
          IF(HDIRVOR > 1.0_JPRB) THEN
            RDIVOR(JLEV,JN)=Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRVOR
          ELSE
            RDIVOR(JLEV,JN)=0.0_JPRB
          ENDIF
          IF(HDIRDIV > 1.0_JPRB) THEN
            RDIDIV(JLEV,JN)=Z_PDILEV(ILEV)*PDISPE(JN,NSMAX,REXPDH)/HDIRDIV
          ELSE
            RDIDIV(JLEV,JN)=0.0_JPRB
          ENDIF
          ! --- end of code for shallow-water equations
        ENDIF
      ENDIF
    ENDDO LEVEL_LOOP1_EC
      
    ! outside sponge 

    DO JLEVG=1,NFLEVG
      IF( Z_PDILEV(JLEVG) <= RBEGSPONGE ) THEN
        IF(HDIRT  > 1.0_JPRB) THEN
          RDITG (JLEVG,JN)=PDISPEL(JN,NSMAX,REXPDH,JLEVG,RBEGSPONGE)/HDIRT
        ELSE
          RDITG(JLEVG,JN)=PDISPEE(JN,NSMAX,REXPDH,JLEVG,RBEGSPONGE)/HDIRVOR
        ENDIF
      ENDIF
    ENDDO
    
    !*       2.1.2 2D variables:
    
    IF (HDIRSP > 1.0_JPRB) THEN
      RDISP(JN)=PDISPE(JN,NSMAX,REXPDH)/HDIRSP
    ELSE
      RDISP(JN)=0.0_JPRB
    ENDIF

  ENDDO WAVE_LOOP_EC

ENDIF ! .NOT. LHDIFFM

!*       2.3   SLHD effect on spectral diffusion
ZP1=4500._JPRB    ! for  45 hPa the diffusion is in full strength
ZP2=15000._JPRB   ! for 150 hPa the diffusion is fully substituted by the SLHD

IF(LSLHD) THEN
  ! Temperature
  IF (LLSLHD_T) THEN
    DO JLEVG=1,NFLEVG
      ZMASK=MIN(1._JPRB,MAX(0._JPRB,ZP2-STPRE(JLEVG))/(ZP2-ZP1)) 
      ! now bounded to <1-SDRED,1>
      ZMASK=1._JPRB - SDRED + SDRED*ZMASK
      RDITG (JLEVG,0:NSMAX)=ZMASK*RDITG (JLEVG,0:NSMAX)
    ENDDO
  ENDIF
  ! All the other variables
  DO JLEV=1,NFLEVL
    ILEV=MYLEVS(JLEV)
    ZMASK=MIN(1._JPRB,MAX(0._JPRB,ZP2-STPRE(ILEV))/(ZP2-ZP1)) 
    ! now bounded to <1-SDRED,1>
    ZMASK=1._JPRB - SDRED + SDRED*ZMASK

    ! Reducing the standard diffusion
    IF (LLSLHD_W) THEN
      RDIDIV(JLEV,0:NSMAX)=ZMASK*RDIDIV(JLEV,0:NSMAX)
      RDIVOR(JLEV,0:NSMAX)=ZMASK*RDIVOR(JLEV,0:NSMAX)
    ENDIF
    IF (YQ%LSP .AND.LLSLHD_Q ) RDIGFL(JLEV,0:NSMAX,YQ%MPSP) = 0._JPRB
    IF (YO3%LSP.AND.LLSLHD_O3) RDIGFL(JLEV,0:NSMAX,YO3%MPSP)= 0._JPRB
    IF (LLSLHD_SPD) RDIPD(JLEV,0:NSMAX)= ZMASK*RDIPD(JLEV,0:NSMAX)
    IF (LLSLHD_SVD) RDIVD(JLEV,0:NSMAX)= ZMASK*RDIVD(JLEV,0:NSMAX)
  ENDDO

  ! Supporting diffusion
  IF (.NOT.LSPECVIS) THEN ! only applied for clasical spectral diffusion
    DO JLEV=1,NFLEVL
      ILEV=MYLEVS(JLEV)

      ! vorticity and divergence
      IF (LLSLHD_W) THEN
        IF( LHDIFFM ) THEN
          IF (RDAMPVORS /= 0._JPRB) THEN
            DO JN=0,NSMAX
              RDSVOR(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,RBEGSPONGE)/ZZDIR/RDAMPVORS
            ENDDO
          ENDIF
          IF (RDAMPDIVS /= 0._JPRB) THEN
            DO JN=0,NSMAX
              RDSDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,RBEGSPONGE)/ZZDIR/RDAMPDIVS
            ENDDO
          ENDIF
        ELSE
          IF ((HDIRVOR > 1._JPRB).AND.(RDAMPVORS /= 0._JPRB)) THEN
            DO JN=0,NSMAX
              RDSVOR(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,RBEGSPONGE)/HDIRVOR/RDAMPVORS
            ENDDO
          ENDIF
          IF ((HDIRDIV > 1._JPRB).AND.(RDAMPDIVS /= 0._JPRB)) THEN
            DO JN=0,NSMAX
              RDSDIV(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,RBEGSPONGE)/HDIRDIV/RDAMPDIVS
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      ! NH vertical divergence
      IF (LLSLHD_SVD) THEN
        IF( LHDIFFM ) THEN
          IF (RDAMPVDS /= 0._JPRB) THEN
            DO JN=0,NSMAX
              RDSVD(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,RBEGSPONGE)/ZZDIR/RDAMPVDS
            ENDDO
          ENDIF
        ELSE
          IF ((HDIRDIV > 1._JPRB).AND.(RDAMPVDS /= 0._JPRB)) THEN
            DO JN=0,NSMAX
              RDSVD(JLEV,JN)=PDISPEL(JN,NSMAX,REXPDHS,ILEV,RBEGSPONGE)/HDIRDIV/RDAMPVDS
            ENDDO
          ENDIF
        ENDIF
      ENDIF

    ENDDO  
  ENDIF

ENDIF

!*       2.4   Diagnostic of e-folding time:

  DO JL=NFLEVL,1,-1
    IL=MYLEVS(JL)
!    IF (IL <= ISTDV2) THEN
      IF(RDIVOR(JL,NSMAX) > 1E-7_JPRB) THEN
        ZEFOLV=1.0_JPRB/RDIVOR(JL,NSMAX)/RHOUR
      ELSE
        ZEFOLV=0._JPRB
      ENDIF
      IF(RDITG(JL,NSMAX) > 1E-7_JPRB) THEN
        ZEFOLT=1.0_JPRB/RDITG(JL,NSMAX)/RHOUR
      ELSE
        ZEFOLT=0._JPRB
      ENDIF
      IF(RDIDIV(JL,NSMAX) > 1E-7_JPRB) THEN
        ZEFOLD=1.0_JPRB/RDIDIV(JL,NSMAX)/RHOUR
      ELSE
        ZEFOLD=0._JPRB
      ENDIF
      WRITE(NULOUT,*) '  E-FOLDING TIME FOR VORTICITY=',ZEFOLV,&
       & ' HOURS AT LEVEL ',IL  
      WRITE(NULOUT,*) '    FOR DIVERGENCE=',ZEFOLD
      WRITE(NULOUT,*) '    FOR T=',ZEFOLT
!    ENDIF
  ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHDF_EC',1,ZHOOK_HANDLE)
END SUBROUTINE SUHDF_EC

