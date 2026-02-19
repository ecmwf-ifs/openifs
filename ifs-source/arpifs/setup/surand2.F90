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

!OPTION! -O nochg -O nodarg -O nodiv -O noiodo -O nomove -O overlap -O noshapeprop -Nfusion -i errchk
SUBROUTINE SURAND2(YDGEOMETRY,YDECUCONVCA,YDML_PHY_STOCH,KSTEP)

!**** *SURAND2*  - Initialize stochastic physic parameters - part 2

!     Purpose.
!     --------
!           Compute random numbers

!**   Interface.
!     ----------
!        *CALL* *SURAND2

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        COMMON STOPH_MIX

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R Buizza and L Isaksen - ECMWF - 

!     Modifications.
!     --------------
!        Original : 97-11-20
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        J.Berner 06-08-11 added cellular automaton backscatter and 
!                          spectral AR backscatter update
!        G. Shutts 09-03-4 spectral backscatter update
!        M. Steinheimer 09-03-05 added option RVP
!        M. Steinheimer 09-10-27 spectral backscatter update (changes to RVP,
!                          changes related to option LSTOPH_SPBS_FAST)
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!        R. El Khatib 22-Mar-2012 Fix uninitialized variables
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        F. Vana        17-Dec-2015 Support for single precision
!        SJ Lock :       Jan  2016  Removed LSTOPH option
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_STOCHAST_MOD , ONLY : MODEL_PHYSICS_STOCHAST_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE RANDOM_NUMBERS_MIX , ONLY : UNIFORM_DISTRIBUTION, GAUSSIAN_DISTRIBUTION
USE YOMLUN   , ONLY : NULOUT
USE YOE_CUCONVCA, ONLY : TECUCONVCA


IMPLICIT NONE

!    ------- more variables for spectral backscatter

TYPE(GEOMETRY)    ,INTENT(IN)  :: YDGEOMETRY
TYPE(TECUCONVCA)  ,INTENT(INOUT):: YDECUCONVCA
TYPE(MODEL_PHYSICS_STOCHAST_TYPE),INTENT(INOUT):: YDML_PHY_STOCH
INTEGER(KIND=JPIM),INTENT(IN)  :: KSTEP

REAL(KIND=JPRB) :: ZCONV, ZLAT, ZLON, ZRAND_NOS_RE, ZRAND_NOS_IM, ZEPS1, ZEPS2
INTEGER(KIND=JPIM) :: ILATBOX, ILONBOX, JROF, J, IIP, JJP, JGL

!     ----------------variables for spectral backscatter

INTEGER(KIND=JPIM) :: ILEV, IM, IN, JMLOC, JN, IJSE, ILEVG, ILEVGPTR
INTEGER(KIND=JPIM) :: INDSP2G_RE_RAN, INDSP2G_IM_RAN

REAL(KIND=JPRB), ALLOCATABLE :: ZRAND_NOS(:)  ! Either one or YDGEOMETRY%YRDIMV%NFLEVG global spectral field of random numbers, 
!                                               dependent on LSTOPH_UNCORR

REAL(KIND=JPRB), ALLOCATABLE :: ZRAND_RVP(:,:)  ! Random vertical profile, used if YDML_PHY_STOCH%YRSTOPH%LSTOPH_RVP == .true.
REAL(KIND=JPRB), ALLOCATABLE :: ZRAND_UNI(:)

!     ----------------variables for random vertical profile
REAL(KIND=JPRB)    :: ZPHI, ZABS, ZRT
INTEGER(KIND=JPIM) :: JI
REAL(KIND=JPRB) :: ZMUL
REAL(KIND=JPRB) :: ZFACT
REAL(KIND=JPRB) :: ZFACT0
REAL(KIND=JPRB) :: ZSTDEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ----------------variables for T-backscatter
REAL(KIND=JPRB), ALLOCATABLE :: ZRAND_NOS_T(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZRAND_RVP_T(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZRAND_UNI_T(:)
REAL(KIND=JPRB)    :: ZPHI_T, ZABS_T, ZRT_T
REAL(KIND=JPRB) :: ZMUL_T
REAL(KIND=JPRB) :: ZRAND_NOS_RE_T, ZRAND_NOS_IM_T

#include "updcelaut.intfb.h"

IF (LHOOK) CALL DR_HOOK('SURAND2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, &
 & YDLAP=>YDGEOMETRY%YRLAP,  YDVAB=>YDGEOMETRY%YRVAB,  &
 & YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA, &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
  & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  YDSPGEOM=>YDGEOMETRY%YSPGEOM,  &
  & YDSTOPH=>YDML_PHY_STOCH%YRSTOPH)

ASSOCIATE(NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, NSPEC2G=>YDDIM%NSPEC2G, &
 & NUMP=>YDDIM%NUMP, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & LCUCONV_CA=>YDECUCONVCA%LCUCONV_CA, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & MYMS=>YDLAP%MYMS, NASM0=>YDLAP%NASM0, &
 & MYLEVS=>YDMP%MYLEVS, &
 & ADLATSTOPH_CA=>YDSTOPH%ADLATSTOPH_CA, ADLONSTOPH_CA=>YDSTOPH%ADLONSTOPH_CA, &
 & ALPHA_STO=>YDSTOPH%ALPHA_STO, ALPHA_STO_T=>YDSTOPH%ALPHA_STO_T, &
 & LSTOPH_CASBS=>YDSTOPH%LSTOPH_CASBS, &
 & LSTOPH_INI=>YDSTOPH%LSTOPH_INI, LSTOPH_JBCOR=>YDSTOPH%LSTOPH_JBCOR, &
 & LSTOPH_RVP=>YDSTOPH%LSTOPH_RVP, LSTOPH_SPBS=>YDSTOPH%LSTOPH_SPBS, &
 & LSTOPH_SPBS_T=>YDSTOPH%LSTOPH_SPBS_T, LSTOPH_TAPER=>YDSTOPH%LSTOPH_TAPER, &
 & LSTOPH_UNCORR=>YDSTOPH%LSTOPH_UNCORR, LSTOPH_UNIFORM=>YDSTOPH%LSTOPH_UNIFORM, &
 & MCELL=>YDSTOPH%MCELL, &
 & NFRSTOPH_SPBS_PAT=>YDSTOPH%NFRSTOPH_SPBS_PAT, NIMRAN=>YDSTOPH%NIMRAN, &
 & NSMAXSPBS=>YDSTOPH%NSMAXSPBS, &
 & ONEMINALPHA_NFRSPBS=>YDSTOPH%ONEMINALPHA_NFRSPBS, &
 & ONEMINALPHA_NFRSPBS_T=>YDSTOPH%ONEMINALPHA_NFRSPBS_T, &
 & RSIGMA2_EPS=>YDSTOPH%RSIGMA2_EPS, RSTOPHCA=>YDSTOPH%RSTOPHCA, &
 & RVP_MUL=>YDSTOPH%RVP_MUL, &
 & RVP_MULFACT=>YDSTOPH%RVP_MULFACT, RVP_MULFACT_T=>YDSTOPH%RVP_MULFACT_T, &
 & RVP_MUL_T=>YDSTOPH%RVP_MUL_T, RWGHT=>YDSTOPH%RWGHT, SPG_AMP=>YDSTOPH%SPG_AMP, &
 & SPG_AMP_T=>YDSTOPH%SPG_AMP_T, SPSTREAM=>YDSTOPH%SPSTREAM, &
 & SPTEMP=>YDSTOPH%SPTEMP, SQRTCORR=>YDSTOPH%SQRTCORR, &
 & TAPER_FACT=>YDSTOPH%TAPER_FACT)
ZCONV = 45._JPRB/ATAN(1.0_JPRB)
ZEPS1 = 1000._JPRB*TINY(ZEPS1)
ZEPS2 = 100._JPRB*EPSILON(ZEPS2)



IF (LSTOPH_CASBS .AND. (.NOT. LCUCONV_CA)) THEN

!     ------------------------
!*    2.0 CASBS
!     ------------------------

!     2.1 UPDATE CELLULAR AUTOMATON  
!     -----------------------------------------------

  IIP= NINT(360._JPRB/ADLONSTOPH_CA)
  JJP= NINT(180._JPRB/ADLATSTOPH_CA)

  CALL UPDCELAUT(YDML_PHY_STOCH,MCELL,RWGHT,IIP,JJP)

!     2.2 generate smooth pattern from CA (on Gaussian grid)
!     -----------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JROF,ZLON,ZLAT,ILATBOX,ILONBOX)
  DO JROF=1,NGPTOT
     ZLON=YDGSGEOM_NB%GELAM(JROF)*ZCONV
     ZLAT=YDGSGEOM_NB%GELAT(JROF)*ZCONV
     ILATBOX=INT((90._JPRB-ZLAT)/(ADLATSTOPH_CA+ZEPS1))+1
     ILONBOX=INT(ZLON/(ADLONSTOPH_CA+ZEPS1))+1
     RSTOPHCA(JROF)= RWGHT(ILONBOX,ILATBOX)*COS(YDGSGEOM_NB%GELAT(JROF))
  ENDDO
!$OMP END PARALLEL DO

ENDIF

IF (LSTOPH_SPBS .AND. MOD(KSTEP,NFRSTOPH_SPBS_PAT) == 0) THEN

!     ------------------------
!*    3.0 SPECTRAL BACKSCATTER
!     ------------------------

  IF (.NOT.LSTOPH_JBCOR.AND..NOT.LSTOPH_UNCORR) THEN
    ! Same random number is used for all vertical levels (fully correlated)
    ! RVP-style correlations may be introduced by setting LSTOPH_RVP=.T.
    ALLOCATE(ZRAND_NOS(NSPEC2G))
    IF (LSTOPH_SPBS_T) ALLOCATE(ZRAND_NOS_T(NSPEC2G))
  ELSE
    ! Different random numbers used for each vertical level (uncorrelated)
    ! Jb-style correlations may be introduced by setting LSTOPH_JBCOR=.T.
    ALLOCATE(ZRAND_NOS(NFLEVG*NSPEC2G))
    IF (LSTOPH_SPBS_T) ALLOCATE(ZRAND_NOS_T(NFLEVG*NSPEC2G))
  ENDIF


  !     3.1 UPDATE MARKOV CHAIN OF SPECTRAL WAVENUMBERS
  !     -----------------------------------------------
  IF (LSTOPH_UNIFORM) THEN
    CALL UNIFORM_DISTRIBUTION(ZRAND_NOS,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOCHPHYS_SPBS)
    IF (LSTOPH_SPBS_T) CALL UNIFORM_DISTRIBUTION(ZRAND_NOS_T,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOPH_SPBS_T)
  ELSE
    CALL GAUSSIAN_DISTRIBUTION(ZRAND_NOS,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOCHPHYS_SPBS)
    ! Change standard deviation of Gaussian distribution to be similar to
    ! standard deviation for a uniform distribution (st.dev=sqrt(1/12.))
    ZSTDEV=SQRT(0.5_JPRB*RSIGMA2_EPS)
    IF (KSTEP == 0)&
     & WRITE(NULOUT,*) ' Standard deviation for Gaussian random distribution: ',ZSTDEV
    ZRAND_NOS(:)=ZSTDEV*ZRAND_NOS(:)
    IF (LSTOPH_SPBS_T) THEN
      CALL GAUSSIAN_DISTRIBUTION(ZRAND_NOS_T,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOPH_SPBS_T)
      ! Change standard deviation of Gaussian distribution to be similar to
      ! standard deviation for a uniform distribution (st.dev=sqrt(1/12.))
      ZSTDEV=SQRT(0.5_JPRB*RSIGMA2_EPS)
      ZRAND_NOS_T(:)=ZSTDEV*ZRAND_NOS_T(:)
    ENDIF
  ENDIF

  !     3.1.1 Random vertical profiles
  !     -----------------------------------------------
  IF (LSTOPH_RVP) THEN
    ALLOCATE(ZRAND_RVP(NFLEVG,NSPEC2))
    ALLOCATE(ZRAND_UNI(NFLEVG*NSPEC2G/2))

    CALL UNIFORM_DISTRIBUTION(ZRAND_UNI,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOCHPHYS_RVP)

    IF (LSTOPH_SPBS_T) THEN
      ALLOCATE(ZRAND_RVP_T(NFLEVG,NSPEC2))
      ALLOCATE(ZRAND_UNI_T(NFLEVG*NSPEC2G/2))
      CALL UNIFORM_DISTRIBUTION(ZRAND_UNI_T,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOPH_RVP_T)
    ENDIF


!$OMP  PARALLEL DO SCHEDULE(STATIC) & 
!$OMP& PRIVATE(JMLOC,IM,JN,IJSE,INDSP2G_RE_RAN,INDSP2G_IM_RAN,ZMUL,ZMUL_T,ZRAND_NOS_RE,ZRAND_NOS_IM,ZPHI,ZPHI_T)&
!$OMP& PRIVATE(ZABS,ZABS_T,JI,ZRT,ZRT_T,ILEVGPTR,ZRAND_NOS_RE_T,ZRAND_NOS_IM_T)
    DO JMLOC=1,NUMP !Number of spectral waves handled by this processor
      IM=MYMS(JMLOC) !actual wave numbers handled by this processor (confiremd: goes from 1-NSMAX)
      DO JN=IM,NSMAXSPBS
        IJSE=NASM0(IM)+2*(JN-IM)                 ! NASM0(0:NSMAX)  : address in a spectral array of (m, n=m)

        INDSP2G_RE_RAN = NIMRAN(IM)+(JN-IM)*2    ! index of real part in a sorted global spectral array of (m, n=m)
        INDSP2G_IM_RAN = NIMRAN(IM)+(JN-IM)*2+1  ! index of imaginary part in a sorted global spectral array of (m, n=m)

        !wavenumber dependent random number scaling
        ZMUL=RVP_MULFACT(JN)

        ZRAND_NOS_RE=ZRAND_NOS(INDSP2G_RE_RAN)
        ZRAND_NOS_IM=ZRAND_NOS(INDSP2G_IM_RAN)

        IF (LSTOPH_UNIFORM) THEN
          ZRAND_NOS_RE=ZRAND_NOS_RE-0.5_JPRB
          ZRAND_NOS_IM=ZRAND_NOS_IM-0.5_JPRB
        ENDIF

        IF (ZRAND_NOS_IM == 0._JPRB .AND. ZRAND_NOS_RE == 0._JPRB) THEN
           ZPHI=0._JPRB
        ELSE
           ZPHI=ATAN2(ZRAND_NOS_IM,ZRAND_NOS_RE)
        ENDIF
        ZABS=(ZRAND_NOS_RE**2+ZRAND_NOS_IM**2)**0.5

        ZRAND_RVP(1,IJSE)  =ZRAND_NOS_RE
        ZRAND_RVP(1,IJSE+1)=ZRAND_NOS_IM

        IF (LSTOPH_SPBS_T) THEN
          ZMUL_T=RVP_MULFACT_T(JN)

          ZRAND_NOS_RE_T=ZRAND_NOS_T(INDSP2G_RE_RAN)
          ZRAND_NOS_IM_T=ZRAND_NOS_T(INDSP2G_IM_RAN)

          IF (LSTOPH_UNIFORM) THEN
            ZRAND_NOS_RE_T=ZRAND_NOS_RE_T-0.5_JPRB
            ZRAND_NOS_IM_T=ZRAND_NOS_IM_T-0.5_JPRB
          ENDIF

          IF (ZRAND_NOS_IM_T == 0._JPRB .AND. ZRAND_NOS_RE_T == 0._JPRB) THEN
            ZPHI_T=0._JPRB
          ELSE
            ZPHI_T=ATAN2(ZRAND_NOS_IM_T,ZRAND_NOS_RE_T)
          ENDIF
          ZABS_T=(ZRAND_NOS_RE_T**2+ZRAND_NOS_IM_T**2)**0.5

          ZRAND_RVP_T(1,IJSE)  =ZRAND_NOS_RE_T
          ZRAND_RVP_T(1,IJSE+1)=ZRAND_NOS_IM_T
        ENDIF


        DO JI=2,NFLEVG
          ILEVGPTR=(JI-1)*NSPEC2G
          ZRT=ZRAND_UNI((ILEVGPTR+INDSP2G_IM_RAN)/2)

          IF ( ZRT > ZEPS1) THEN   !to avoid infinit LOG(0) -> there are for sure better ways to do that
             ZRT=ZRT-0.5_JPRB
             ZRT=SIGN(1._JPRB,ZRT)*LOG(MAX(1._JPRB-2._JPRB*ABS(ZRT),ZEPS1))
          ELSE
              ZRT=0._JPRB
          ENDIF

         !scale random number and calculate new phase
          ZPHI=ZPHI+ZMUL*RVP_MUL(JI)*ZRT

          ZRAND_RVP(JI,IJSE)=ZABS*COS(ZPHI)
          ZRAND_RVP(JI,IJSE+1)=ZABS*SIN(ZPHI)

          IF (LSTOPH_SPBS_T) THEN
            ZRT_T=ZRAND_UNI_T((ILEVGPTR+INDSP2G_IM_RAN)/2)

            IF ( ZRT_T > ZEPS1 ) THEN   !to avoid infinit LOG(0) -> there are for sure better ways to do that
              ZRT_T=ZRT_T-0.5_JPRB
              ZRT_T=SIGN(1._JPRB,ZRT_T)*LOG(MAX(1._JPRB-2._JPRB*ABS(ZRT_T),ZEPS1))
            ELSE
              ZRT_T=0._JPRB
            ENDIF

           !scale random number and calculate new phase
            ZPHI_T=ZPHI_T+ZMUL_T*RVP_MUL_T(JI)*ZRT_T

            ZRAND_RVP_T(JI,IJSE)=ZABS_T*COS(ZPHI_T)
            ZRAND_RVP_T(JI,IJSE+1)=ZABS_T*SIN(ZPHI_T)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF


!$OMP  PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(ILEV,ILEVG,ILEVGPTR,JMLOC,IM,JN,IJSE) &
!$OMP& PRIVATE(ZRAND_NOS_RE,ZRAND_NOS_IM,JGL,INDSP2G_RE_RAN,INDSP2G_IM_RAN) &
!$OMP& PRIVATE(ZFACT,ZFACT0)&
!$OMP& PRIVATE(ZRAND_NOS_RE_T,ZRAND_NOS_IM_T)
  DO ILEV=1,NFLEVL
    ILEVG=MYLEVS(ILEV)   ! To assure perturbations are reproducible if NFLEVL is changed
    ILEVGPTR=(ILEVG-1)*NSPEC2G

    ZFACT=1._JPRB
    IF (LSTOPH_TAPER) THEN
       ZFACT=TAPER_FACT(ILEVG)
    ENDIF

    DO JMLOC=1,NUMP !Number of spectral waves handled by this processor
      IM=MYMS(JMLOC) !actual wave numbers handled by this processor (confiremd: goes from 1-NSMAX)
      DO JN=IM,NSMAXSPBS
        IJSE=NASM0(IM)+2*(JN-IM) ! NASM0(0:NSMAX)  : address in a spectral array of (m, n=m)

        INDSP2G_RE_RAN = NIMRAN(IM)+(JN-IM)*2    ! index of real part in a sorted global spectral array of (m, n=m)
        INDSP2G_IM_RAN = NIMRAN(IM)+(JN-IM)*2+1  ! index of imaginary part in a sorted global spectral array of (m, n=m)

        ZRAND_NOS_RE=0.0_JPRB
        ZRAND_NOS_IM=0.0_JPRB
        IF (LSTOPH_SPBS_T) THEN
          ZRAND_NOS_RE_T=0.0_JPRB
          ZRAND_NOS_IM_T=0.0_JPRB
        ENDIF

        IF (LSTOPH_UNCORR) THEN
          ZRAND_NOS_RE=ZRAND_NOS(ILEVGPTR+INDSP2G_RE_RAN)
          ZRAND_NOS_IM=ZRAND_NOS(ILEVGPTR+INDSP2G_IM_RAN)
          IF (LSTOPH_SPBS_T) THEN
            ZRAND_NOS_RE_T=ZRAND_NOS_T(ILEVGPTR+INDSP2G_RE_RAN)
            ZRAND_NOS_IM_T=ZRAND_NOS_T(ILEVGPTR+INDSP2G_IM_RAN)
          ENDIF
        ELSEIF (LSTOPH_JBCOR) THEN
          DO JGL=1,ILEVG
            ZRAND_NOS_RE=ZRAND_NOS_RE+ZRAND_NOS((JGL-1)*NSPEC2G+INDSP2G_RE_RAN)*&
               & SQRTCORR(JGL,ILEVG)

            ZRAND_NOS_IM=ZRAND_NOS_IM+ZRAND_NOS((JGL-1)*NSPEC2G+INDSP2G_IM_RAN)*&
               & SQRTCORR(JGL,ILEVG)
          ENDDO
          IF (LSTOPH_SPBS_T) THEN
            DO JGL=1,ILEVG
              ZRAND_NOS_RE_T=ZRAND_NOS_RE_T+ZRAND_NOS_T((JGL-1)*NSPEC2G+INDSP2G_RE_RAN)*&
                & SQRTCORR(JGL,ILEVG)

              ZRAND_NOS_IM_T=ZRAND_NOS_IM_T+ZRAND_NOS_T((JGL-1)*NSPEC2G+INDSP2G_IM_RAN)*&
                & SQRTCORR(JGL,ILEVG)
            ENDDO
          ENDIF
        ELSEIF (LSTOPH_RVP) THEN
            ZRAND_NOS_RE=ZRAND_RVP(ILEVG,IJSE)
            ZRAND_NOS_IM=ZRAND_RVP(ILEVG,IJSE+1)
            IF (LSTOPH_SPBS_T) THEN
              ZRAND_NOS_RE_T=ZRAND_RVP_T(ILEVG,IJSE)
              ZRAND_NOS_IM_T=ZRAND_RVP_T(ILEVG,IJSE+1)
            ENDIF
        ELSE
            ZRAND_NOS_RE=ZRAND_NOS(INDSP2G_RE_RAN)
            ZRAND_NOS_IM=ZRAND_NOS(INDSP2G_IM_RAN)
            IF (LSTOPH_SPBS_T) THEN
              ZRAND_NOS_RE_T=ZRAND_NOS_T(INDSP2G_RE_RAN)
              ZRAND_NOS_IM_T=ZRAND_NOS_T(INDSP2G_IM_RAN)
            ENDIF
        ENDIF
        IF (LSTOPH_UNIFORM.AND.(.NOT.LSTOPH_RVP)) THEN  
          ! it has to be verified if using LSTOPH_JBCOR is "compatible" with using LSTOPH_UNIFORM 
          ! i.e. if the current implementation is doing what is expected for both set to true
          ZRAND_NOS_RE=ZRAND_NOS_RE-0.5_JPRB
          ZRAND_NOS_IM=ZRAND_NOS_IM-0.5_JPRB
          IF (LSTOPH_SPBS_T) THEN
            ZRAND_NOS_RE_T=ZRAND_NOS_RE_T-0.5_JPRB
            ZRAND_NOS_IM_T=ZRAND_NOS_IM_T-0.5_JPRB
          ENDIF
        ENDIF

        IF (KSTEP == 0 .AND. LSTOPH_INI) THEN
         !
         !initialization of Markov Chain
         !
           !WRITE(NULOUT,*) ' Initializing Markov Chain '
           ZFACT0=1._JPRB / MAX(SQRT( 1._JPRB -(1._JPRB-ALPHA_STO(JN))**(2*NFRSTOPH_SPBS_PAT) ),ZEPS1)
           ZRAND_NOS_RE=ZFACT0*ZRAND_NOS_RE
           ZRAND_NOS_IM=ZFACT0*ZRAND_NOS_IM
           IF (LSTOPH_SPBS_T) THEN
             ZFACT0=1._JPRB / MAX(SQRT( 1._JPRB -(1._JPRB-ALPHA_STO_T(JN))**(2*NFRSTOPH_SPBS_PAT) ),ZEPS1)
             ZRAND_NOS_RE_T=ZFACT0*ZRAND_NOS_RE_T
             ZRAND_NOS_IM_T=ZFACT0*ZRAND_NOS_IM_T
           ENDIF
        ENDIF

        IF (IM==0) THEN
          !special treatment for m=0 -> imaginary part=0, variance of real part must be inflated to keep
          !"complex" variance equal to modes with m ne 0
          SPSTREAM(ILEV,IJSE)  = (ONEMINALPHA_NFRSPBS(JN)*SPSTREAM(ILEV,IJSE)  +&
            &                    SQRT(2._JPRB)*SPG_AMP(JN)*ZRAND_NOS_RE)*ZFACT
          SPSTREAM(ILEV,IJSE+1)= 0._JPRB

          IF (LSTOPH_SPBS_T) THEN
            SPTEMP(ILEV,IJSE)  = (ONEMINALPHA_NFRSPBS_T(JN)*SPTEMP(ILEV,IJSE)  +&
            &                    SQRT(2._JPRB)*SPG_AMP_T(JN)*ZRAND_NOS_RE_T)*ZFACT
            SPTEMP(ILEV,IJSE+1)= 0._JPRB
          ENDIF
        ELSE
          SPSTREAM(ILEV,IJSE)  = (ONEMINALPHA_NFRSPBS(JN)*SPSTREAM(ILEV,IJSE)  +&
            &                    SPG_AMP(JN)*ZRAND_NOS_RE)*ZFACT
          SPSTREAM(ILEV,IJSE+1)= (ONEMINALPHA_NFRSPBS(JN)*SPSTREAM(ILEV,IJSE+1) +&
            &                    SPG_AMP(JN)*ZRAND_NOS_IM)*ZFACT

          IF (LSTOPH_SPBS_T) THEN
            SPTEMP(ILEV,IJSE)  = (ONEMINALPHA_NFRSPBS_T(JN)*SPTEMP(ILEV,IJSE)  +&
            &                    SPG_AMP_T(JN)*ZRAND_NOS_RE_T)*ZFACT
            SPTEMP(ILEV,IJSE+1)= (ONEMINALPHA_NFRSPBS_T(JN)*SPTEMP(ILEV,IJSE+1) +&
            &                    SPG_AMP_T(JN)*ZRAND_NOS_IM_T)*ZFACT
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  DEALLOCATE(ZRAND_NOS)
  IF (LSTOPH_RVP) THEN
    DEALLOCATE(ZRAND_RVP)
    DEALLOCATE(ZRAND_UNI)
  ENDIF

  IF (LSTOPH_SPBS_T) THEN
    DEALLOCATE(ZRAND_NOS_T)
    IF (LSTOPH_RVP) THEN
      DEALLOCATE(ZRAND_RVP_T)
      DEALLOCATE(ZRAND_UNI_T)
    ENDIF
  ENDIF

  LSTOPH_INI=.FALSE.    !set to false to avoid multiple initialization
                        !(surand2 is called more than once for KSTEP == 0)

ENDIF


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURAND2',1,ZHOOK_HANDLE)
END SUBROUTINE SURAND2
