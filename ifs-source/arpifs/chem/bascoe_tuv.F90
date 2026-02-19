! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_TUV ( KMXCLY, PALBEDO, KLEVTURBO, PFBEAMR,        & ! input
                      & P_HDENS, PDENS, PTEMPER, PSZA, PZS, P_CP,   & ! input
                      & PDRAT, PWRAT )                                ! output

!**   DESCRIPTION
!     ----------
!
!   Computes photodissociation rates for the stratospheric chemistry
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_TUV* IS CALLED FROM *BASCOE_J_CALC*.
!---------------------------------------------------------------------
!               I N P U T    V A R I A B L E S
!---------------------------------------------------------------------
!
!     KMXCLY       : vertical dimension= nb of absorbing layers= nb_levels-1
!
!     PALBEDO      : bottom-boundary albedo for Lambertian reflecting
!                   surface.
!
!     KLEVTURBO = alt index (between 0=top to KMXCLY=bottom) of turbopause
!
!     PFBEAMR(MXWVN)  constant intensity of incident solar  beam at
!                top boundary. ( no. of photons/( cm**2 s) )
!
!     P_HDENS(0:KMXCLY,IDENS) : scale height of species idens (cm)
!                    idens = 1 = air              air
!                    idens = 2 = O3abs            o3
!                    idens = 3 = O2abs            o2
!                    idens = 4 = NOabs            no
!                    idens = 5 = CO2abs           co2
!                    idens = 6 = NO2abs           no2
!
!    dens(0:KMXCLY,idens): density of species idens (molecules/cm**3)
!
!    PTEMPER(lev) lev = 0 to KMXCLY, temperatures (K) at levels.
!                (note that temperature is specified at levels
!                rather than for layers.)  don't forget to put top
!                temperature in 'PTEMPER(0)', not in 'PTEMPER(1)'.
!                Used to calculate the temperature dependence in
!                cross-sections etc.
!
!    PSZA         Solar zenith angle of incident beam (positive, degrees).
!                ** warning **  if this is close to one of the
!                computational polar angles, serious ill-conditioning
!                and a possible crash of 'TWOSTREAM_*' might
!                result;  hence this is flagged as a fatal error.
!
!    PZS(lev)      lev = 0 to KMXCLY, *geometrical* altitude of level (km) :
!                 PZS(0) = upper boundary altitude , PZS(KMXCLY) = 0km
!
!    P_CP(lev)      lev = 0 to KMXCLY, specific heat of air at cst p (J/g/K)
!
!---------------------------------------------------------------------
!               O U T P U T    V A R I A B L E S
!---------------------------------------------------------------------
!
!    PDRAT(iproc,lev):  photodissociation rate (1/s) of process iproc
!    PWRAT(1,0:KMXCLY) : O2 heating (K/s ?)
!    PWRAT(2,0:KMXCLY) : O3 heating (K/s ?)
!
!
!
!     MODIFICATIONS.
!     -------
!       2017-12-07:    YVES CHRISTOPHE    *BIRA* !YC
!           BASCOE TUV_midatm and co adapted and included in ifs/chem
!       2019-12-18:    YVES CHRISTOPHE    *BIRA* !YC
!           allow date-dependent solar flux
!


USE YOMLUN    , ONLY : NULOUT
USE PARKIND1  ,ONLY  : JPIM, JPRB, JPRD              ! (JPRD needed for fcttim.func)
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_MODULE, only : jprocmax => ndiss, J_NO
USE BASCOE_TUV_MODULE, only : MXWVN, sw_lya, iv_lya, nwvn_srb, iv_srb0, iv_srb1, &
                            & nheatspec, nabspec, air, O2abs, O3abs, raycrs, &
                            & init_done, sza_chap


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM), INTENT(IN)                             :: KMXCLY
  REAL(KIND=JPRB), INTENT(IN)                                :: PALBEDO, PSZA
  REAL(KIND=JPRB), DIMENSION(MXWVN), INTENT(IN)              :: PFBEAMR
  INTEGER(KIND=JPIM), optional, INTENT(IN)                   :: KLEVTURBO
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN)            :: PTEMPER, PZS, P_CP
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY,nabspec), INTENT(IN)    :: P_HDENS, PDENS
  REAL(KIND=JPRB), DIMENSION(jprocmax,0:KMXCLY), INTENT(OUT)  :: PDRAT
  REAL(KIND=JPRB), optional, DIMENSION(nheatspec,0:KMXCLY), INTENT(OUT) :: PWRAT


!-----------------------------------------------------------------------
!*      0.2 Local PARAMETERS
!-----------------------------------------------------------------------
  CHARACTER(LEN=*), PARAMETER    :: CL_MY_NAME   = 'BASCOE_TUV'
  LOGICAL, PARAMETER :: LLSPHERE = .TRUE.                   ! .true. to calc pseudo-spherical
  REAL(KIND=JPRB), PARAMETER :: ZPIR = 3.1415926535897931E0 ! pi
  REAL(KIND=JPRB), PARAMETER :: ZD2R = ZPIR/180.E0          ! degree to radian
  REAL(KIND=JPRB), PARAMETER :: ZSOLID_ANGLEI = 1./(4. * ZPIR)
  REAL(KIND=JPRB), PARAMETER :: ZKBOLTZ = 1.38065E-23       ! Boltzmann constant (J/K/molec)

!-----------------------------------------------------------------------
!       ... Local variables
!---------------------------------------------------------------------
!               I n t e r n a l    v a r i a b l e s
!---------------------------------------------------------------------
!
!    babs(lev,idens)   lev = 1 to KMXCLY, idens=2 to nabspec
!               absorption coeffs of each absorbing species in each
!               alitude layer (unitless). No val for total air (idens=1)
!
!    ZBSCAR(lev)   lev = 1 to KMXCLY,
!                scattering coefficient for rayleigh scattering
!
!    laydens(lev,idens)  lev = 1 to KMXCLY, idens=1 to nabspec,
!               densities of each absorbing species in each alitude layer (cm-2)
!
!    CRS(iproc,iv) iproc = 1, phtmax, iv = ivbegin, ivend
!                reaction cross-sections for the different photo-
!                chemical processes at wavelength interval iv as read
!                from datafiles crs97.dat
!
!    PCRST(iproc,lev,iv) iproc = 1, phtmax ; iv = ivbegin, ivend
!                reaction cross-sections for the different photo-
!                chemical processes at wavelength interval iv ,
!                temperature(lev) and pressure adjusted, ( cm**2 )
!
!    mu2(lev,idens)    =cos(sza) if sza < 65deg, or 1/chapman in which case
!                       mu2 depends on the absorbant in heterosphere
!
!    raycrs(iv)  Rayleigh scattering cross-section (cm2)
!
!    dtauc(lev,iv)   lev = 1 to KMXCLY ; iv =1 to mxwvn ,
!                optical depths of computational layers
!
!    gg(lev,iv)  lev = 1 to KMXCLY, iv = 1 to mxwvn
!                1st coeff. in Legendre polynomial expansions of phase
!                functions for computational layers => single scattering
!                albedo asymmetry factor. Previously named pmom(1,lev) .
!                Set to zero (clouds/aerosols removed) - kept for TWSTR
!
!    iv  :  For wavelength interval indexing
!
!    lev     Index of vertical layer (1 at the top to KMXCLY at the bottom)
!            or of vertical level (0 at the top to KMXCLY at the bottom)
!
!    ssalb(lev,iv)   lev = 1 to KMXCLY ; iv =1 to mxwvn ,
!                single-scatter albedos of computational layers
!
!    uavg(lev,iv)   lev = 0 to KMXCLY ; iv =1 to mxwvn ,
!                mean intensity (including the direct beam)
!
! Used only for pseudo-spherical calculation:
!
!    dsdh = slant path of direct beam through each layer crossed
!            when travelling from the top of the atmosphere to layer i;
!            dsdh(i,j), i = 0..KMXCLY, j = 1..KMXCLY   (for pseudo-spherical
!
!    INID = INTEGER, number of layers crossed by the direct beam when
!            travelling from the top of the atmosphere to layer i;
!            nid(i), i = 0..KMXCLY
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM) :: IDENS,IV, ILEV
  INTEGER(KIND=JPIM), DIMENSION(0:KMXCLY)  :: INID
  REAL(KIND=JPRB), DIMENSION(KMXCLY)   :: ZBSCAR, ZTAUTOTAL, ZDELSC, ZDTAU_O3
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY) :: ZFDR, ZFUP, ZFDN, ZJJNO, ZPMB, ZHEATCON
  REAL(KIND=JPRB), DIMENSION(KMXCLY,MXWVN)     ::  ZDTAUC, ZSSALB, ZGG
  REAL(KIND=JPRB), DIMENSION(KMXCLY,nabspec)   :: ZLAYDENS
  REAL(KIND=JPRB), DIMENSION(KMXCLY,2:nabspec) :: ZBABS
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY,nabspec) :: ZMU2, ZPCOLX
  REAL(KIND=JPRB), DIMENSION(jprocmax,0:KMXCLY,MXWVN) :: ZCRST, ZQYTP
  REAL(KIND=JPRB), DIMENSION(MXWVN,0:KMXCLY)   :: ZUAVG
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY,KMXCLY)  :: ZDSDH
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY)   :: ZRMLYA, ZRO2LYA  ! reduction factors for Lyman-a "line"
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY,nwvn_srb) :: ZRM, ZRO2 ! reduction factors for Schumann-R bands
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY,MXWVN)    :: ZEFSHO2, ZEFSHO3
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  REAL(KIND=JPRB) :: ZCHAPMAN_CORR


!-----------------------------------------------------------------------
#include "abor1.intfb.h"
#include "fcttim.func.h"
!-------------------------------------------------------------------
#include "bascoe_pho_chapman.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV',0,ZHOOK_HANDLE )

      IF( .not. init_done ) THEN
         WRITE(NULOUT,*)CL_MY_NAME//' usage error: BADSCOE_TUV_INI must be called first'
         CALL ABOR1(CL_MY_NAME//' usage error: BASCOE_TUV_INI must be called first')
      ENDIF

      IF( any( PZS(1:KMXCLY) > PZS(0:KMXCLY-1) ) ) THEN
         WRITE(NULOUT,*)CL_MY_NAME//' usage error: PZS is not strictly decreasing'
         CALL ABOR1(CL_MY_NAME//', usage error: PZS is not strictly decreasing')
      ENDIF


!-----------------------------------------------------------------------
!       ...Obtain ZLAYDENS, the densities in each horizontal layer
!-----------------------------------------------------------------------
      call GET_LAYDENS( KMXCLY, PZS, PDENS, ZLAYDENS)

!-----------------------------------------------------------------------
!       ... Correct absorption cross-sections for T dependence
!           and calc T-dependent for quantum yield of J(O3)
!-----------------------------------------------------------------------
      call CRS_TDEP( KMXCLY, PTEMPER, ZCRST )
      
      ZPMB(0:KMXCLY) = 1.E4_JPRB * PDENS(0:KMXCLY,air) * ZKBOLTZ * PTEMPER(0:KMXCLY)
      call CALC_QYTP( KMXCLY, PTEMPER, ZPMB, ZQYTP )
      
      IF( PSZA < sza_chap ) THEN
         ZMU2(0:KMXCLY,1:NABSPEC) = COS( PSZA*ZD2R )
      ELSE
!-----------------------------------------------------------------------
!       ... Chapman correction for solar zenith angles > 65 degrees
!           The result variable mu2 is altitude-dependent
!           It is also species-dependent above the turbopause
!-----------------------------------------------------------------------
         IF( present( KLEVTURBO ) ) THEN
            DO ILEV = 0, KLEVTURBO-1
               DO IDENS = 1, nabspec
                  call BASCOE_PHO_CHAPMAN( PZS(ILEV), PSZA, ZCHAPMAN_CORR, 1.E-5_JPRB*P_HDENS(ILEV,IDENS) )
                  ZMU2(ILEV,IDENS) = 1. / ZCHAPMAN_CORR
               ENDDO
            ENDDO
            DO ILEV = KLEVTURBO, KMXCLY
               call BASCOE_PHO_CHAPMAN( PZS(ILEV), PSZA, ZCHAPMAN_CORR, 1.E-5_JPRB*P_HDENS(ILEV,air) )
               ZMU2(ILEV,1:NABSPEC) = 1. / ZCHAPMAN_CORR
            ENDDO
         ELSE
            DO ILEV = 0, KMXCLY
               call BASCOE_PHO_CHAPMAN( PZS(ILEV), PSZA, ZCHAPMAN_CORR )
               ZMU2(ILEV,1:NABSPEC) = 1. / ZCHAPMAN_CORR
            ENDDO
        ENDIF
      ENDIF

      IF( LLSPHERE ) THEN
!-----------------------------------------------------------------------
!       ... Calculate slant path lengths in spherical geometry
!-----------------------------------------------------------------------
         call SPHERL( KMXCLY, PZS, PSZA, ZDSDH, INID )
      ENDIF

!-----------------------------------------------------------------------
!       ... Calc slant overhead column amounts (cm-2) for SRB, LYA, J_NO_PARAMJGR93
!-----------------------------------------------------------------------
      ZPCOLX(0,1:NABSPEC) = PDENS(0,1:NABSPEC) * P_HDENS(0,1:NABSPEC)
      DO ILEV = 1, KMXCLY
         ZPCOLX(ILEV,1:NABSPEC) = ZPCOLX(ILEV-1,1:NABSPEC) + ZLAYDENS(ILEV,1:NABSPEC)
      ENDDO
      DO IDENS = 1, nabspec
         ZPCOLX(0:KMXCLY,IDENS) = ZPCOLX(0:KMXCLY,IDENS) / ZMU2(0:KMXCLY,IDENS)   ! SLANT column contents
      ENDDO

!-----------------------------------------------------------------------
!       ... Calc reduction factors ro2lya and rmlya for
!           Chabrillat/Kockarts Lyman-alpha line parametrization.
!-----------------------------------------------------------------------
      CALL LYA( KMXCLY, ZPCOLX(0:KMXCLY,O2abs), ZRMLYA, ZRO2LYA )

!-----------------------------------------------------------------------
!       ... Calculate reduction factors ro2 and rm for Kockarts
!           Schumann-Runge bands parametrization.
!-----------------------------------------------------------------------
      CALL SRB( KMXCLY, ZPCOLX(0:KMXCLY,O2abs), ZRM, ZRO2 )

!-----------------------------------------------------------------------
!       ... Loop on wavelength intervals
!-----------------------------------------------------------------------
      DO iv = 1,MXWVN
         !VH change initialization value from zero towards 1E-7 
         !VH in order to prevent spurious division by zero in case running in single-precision.
         ZGG(1:KMXCLY,iv) = 1.E-7_JPRB
!---------------------------------------------------------------------
!       ... Calculate the absorption optical depths due to molec abs.
!---------------------------------------------------------------------
         CALL ABSDEP( KMXCLY, ZCRST(1:JPROCMAX,0:KMXCLY,iv), PDENS, PZS, ZBABS )

!---------------------------------------------------------------------
!       ... Rayleigh scattering
!---------------------------------------------------------------------
         ZBSCAR(1:KMXCLY) = ZLAYDENS(1:KMXCLY,air) * raycrs(iv)

!---------------------------------------------------------------------
!       ... Set total optical depth and single scattering albedo
!---------------------------------------------------------------------
         ZDTAUC(1:KMXCLY,iv) = ZBSCAR(1:KMXCLY) &
     &                     + SUM( ZBABS(1:KMXCLY,2:nabspec), DIM=2 )
         ZSSALB(1:KMXCLY,iv) = ZBSCAR(1:KMXCLY) / ZDTAUC(1:KMXCLY,iv)

!---------------------------------------------------------------------
!       ... Solve the one-dimensional radiative transfer equation
!           to obtain mean intensities for each wavelength.
!---------------------------------------------------------------------
         ZUAVG(iv,0:KMXCLY) = 0.

         IF( iv > iv_srb1 ) THEN
!---------------------------------------------------------------------
!       ... Wavelengths longer than SRB
!           (real criterium is ZTAUTOTAL below >> 1) : call Delta-Eddington
!           two-stream method implemented by Sasha Madronich,
!           pseudo-spherical *or* plane parallel version
!---------------------------------------------------------------------
            IF( LLSPHERE ) THEN
               CALL TWOSTREAM_SPHERE( KMXCLY, PSZA, ZMU2(0:KMXCLY,air), PALBEDO, &
     &                                ZDTAUC(1:KMXCLY,iv), ZSSALB(1:KMXCLY,iv), ZGG(1:KMXCLY,iv), &
     &                                ZDSDH, INID, ZFDR, ZFUP, ZFDN )
             ELSE
               CALL TWOSTREAM_PLANE( KMXCLY, ZMU2(0:KMXCLY,air), PALBEDO, ZDTAUC(1:KMXCLY,iv), &
     &                               ZSSALB(1:KMXCLY,iv), ZGG(1:KMXCLY,iv), &
     &                               ZFDR, ZFUP, ZFDN )
            ENDIF

!---------------------------------------------------------------------
!       ... Total actinic flux= (fdr+dfup+fdn) * fbeamr
!           mean Intensity uavg = actinic flux/4*pi
!---------------------------------------------------------------------
            ZUAVG(iv,:) = (ZFDR + ZFUP + ZFDN)* PFBEAMR(iv) * ZSOLID_ANGLEI

         ELSE
!-----------------------------------------------------------------------
!       ... At SRB and shorter wl: direct absorption and single scattering
!           delta scaling (ZDELSC) for total optical depth ZTAUTOTAL
!-----------------------------------------------------------------------
            ZUAVG(iv,0) = 1.
            ZDELSC(1:KMXCLY) = (1. - ZSSALB(1:KMXCLY,iv) * ZGG(1:KMXCLY,iv) * ZGG(1:KMXCLY,iv))
            IF( iv < iv_srb0 .and. ( .not. sw_lya .or. &
     &                             (sw_lya .and. iv /= iv_lya ) ) ) THEN
               ZTAUTOTAL(1:KMXCLY) = ZDELSC(1:KMXCLY) * &
     &                       ( ZBSCAR(1:KMXCLY) / ZMU2(1:KMXCLY,air) &
     &                        + SUM( ZBABS(1:KMXCLY,2:nabspec) &
     &                              /ZMU2(1:KMXCLY,2:nabspec), DIM=2 ) )
               DO ILEV = 1,KMXCLY
                  ZUAVG(iv,ILEV) = ZUAVG(iv,ILEV-1) * EXP( -ZTAUTOTAL(ILEV) )
               ENDDO
             ELSE
!-----------------------------------------------------------------------
!       ... Schumann-Runge Bands and Lyman-alpha line: use absorption
!           depth with only O3 contribution, as O2 contrib will be
!           PARAMETERized and the other contribs are neglected
!-----------------------------------------------------------------------
               ZDTAU_O3(:) = ZDELSC(:) * ZBABS(:,O3abs) /ZMU2(1:KMXCLY,O3abs)
               DO ILEV = 1,KMXCLY
                  ZUAVG(iv,ILEV) = ZUAVG(iv,ILEV-1) * EXP( -ZDTAU_O3(ILEV) )
               ENDDO
            ENDIF
            ZUAVG(iv,1:KMXCLY) = ZUAVG(iv,1:KMXCLY) * PFBEAMR(iv) * ZSOLID_ANGLEI

                                !  PDRAT will be multiplied by fourpi in JRATES
         ENDIF
      ENDDO                        ! of wavelength loop

      PDRAT = 0.
!-----------------------------------------------------------------------
!       ... Integrate over wavelength to obtain photodissociation
!           and heating rates.
!-----------------------------------------------------------------------
      CALL JRATES( KMXCLY, ZCRST, ZQYTP,&
     &             ZUAVG, ZRMLYA, ZRO2LYA, ZRM, ZRO2, PDRAT )

!-----------------------------------------------------------------------
!       ... J(NO) PARAMETERization  (Minschwaner and Siskind '93)
!-----------------------------------------------------------------------
      CALL J_NO_PARAMJGR93( KMXCLY, ZPCOLX, ZUAVG, PDENS(0:KMXCLY,air), ZJJNO )
      PDRAT(J_NO,0:KMXCLY) = ZJJNO


!-----------------------------------------------------------------------
!       ... Optionally, calculate heating rates using JRATES results
!-----------------------------------------------------------------------
      IF( .not. PRESENT( PWRAT ) ) THEN
        IF (LHOOK) CALL DR_HOOK('BASCOE_TUV',1,ZHOOK_HANDLE )
        RETURN
      ENDIF

      CALL SOLHEAT_EFFIC( KMXCLY, ZPMB, ZEFSHO2, ZEFSHO3 )
      CALL SOLHEAT_RATES( KMXCLY, PZS, ZCRST, ZUAVG, ZEFSHO2, ZEFSHO3, PWRAT )
      ZHEATCON(0:KMXCLY) = 1. / ( 4.8097e-23 * P_CP(0:KMXCLY) )
      PWRAT(1,0:KMXCLY) = ZHEATCON * PWRAT(1,0:KMXCLY) * PDENS(0:KMXCLY,O2abs) / PDENS(0:KMXCLY,air)
      PWRAT(2,0:KMXCLY) = ZHEATCON * PWRAT(2,0:KMXCLY) * PDENS(0:KMXCLY,O3abs) / PDENS(0:KMXCLY,air)

      DO ILEV=0,KMXCLY
                  IF (PWRAT(2,ILEV)<0.0) WRITE(*,*) 'in TUV: ',PWRAT(2,ILEV)
      ENDDO



IF (LHOOK) CALL DR_HOOK('BASCOE_TUV',1,ZHOOK_HANDLE )

END SUBROUTINE BASCOE_TUV




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! from original BASCOE MODULE TUV_MIDATM
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE TWOSTREAM_PLANE( KMXCLY, PMU2, PALBEDO, PDTAUCIV, PSSALBIV, PGIV, &
                          &   PFDR, PFUP, PFDN )
!-----------------------------------------------------------------------
!       ... Two-stream equations for multiple layers based on equations
!           from Toon et al., JGR, Vol 94, #d13  Nov. 20, 1989
!           To be called separately for each wavelength interval iv
!           Plane-parallel + Chapman approximation
!-----------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... PARAMETERS
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), PARAMETER :: ZEPS = 1.E-3, ZPRECIS = 1.E-7
      REAL(KIND=JPRB), PARAMETER :: ZPIFS = 1.
      REAL(KIND=JPRB), PARAMETER :: ZFDN0 = 0.

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
!       mu2 = cosine of solar zenith angle if this angle < 65 degrees
!             ( mu2 independent of altitude ) OR 1/Chapman fct
!             if angle > 65 degrees ( mu2(lev) is then lev-dependent )
!             Using Chapman fct based on scale height H for air
!       albedo = surface albedo
!       PDTAUCIV =  unscaled optical depth of each layer
!       PSSALBIV  =  unscaled single scattering albedo
!       PGIV   =  unscaled asymmetry factor
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN)                   :: KMXCLY
      REAL(KIND=JPRB), INTENT(IN)                      :: PALBEDO
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN)  :: PMU2
      REAL(KIND=JPRB), DIMENSION(KMXCLY),   INTENT(IN)  :: PDTAUCIV, PSSALBIV, PGIV
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(OUT) :: PFUP, PFDN, PFDR

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM) :: ILEV, JL
      INTEGER(KIND=JPIM) :: IROW, IKNROWS
      REAL(KIND=JPRB)  :: ZTAUC
      REAL(KIND=JPRB)  :: ZTEMPG
      REAL(KIND=JPRB)  :: ZSSFC
      REAL(KIND=JPRB)  :: ZTAUG
      REAL(KIND=JPRB)  :: ZEXPON, ZEXPON0, ZEXPON1
      REAL(KIND=JPRB)  :: ZDIVISR, ZUP, ZDN, ZTEMP
      REAL(KIND=JPRB)  :: ZGAM1, ZGAM2, ZGAM3, ZGAM4
      REAL(KIND=JPRB)  :: ZOM, ZTAU, ZF, ZG
      REAL(KIND=JPRB), DIMENSION(KMXCLY) :: ZLAM, ZTAUN, ZBGAM, ZE1, ZE2, ZE3, ZE4, &
     &                          ZCUP, ZCDN, ZCUPTN, ZCDNTN, ZMU1
      REAL(KIND=JPRB), DIMENSION(2*KMXCLY) :: ZA, ZB, ZD, ZY   ! 2*KMXCLY = nb of rows in the matrix
      REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TWOSTREAM_PLANE',0,ZHOOK_HANDLE )
!-----------------------------------------------------------------------
!       ... Initial conditions:  pi*solar flux = 1;
!                                diffuse incidence = 0
!-----------------------------------------------------------------------
!       ... Compute coefficients for each layer:
!           gam1 - gam4 = 2-stream coefficients,
!           different for different approximations
!           expon0 = calculation of e when tau is zero
!           expon1 = calculation of e when tau is taun
!           cup and cdn = calculation when tau is zero
!           cuptn and cdntn = calc. when tau is taun
!           ZDIVISR = prevents division by zero
!-----------------------------------------------------------------------
      IKNROWS = 2*KMXCLY
      ZTAUC = 0.
      DO JL = 1,KMXCLY
         ZG = PGIV(JL)
         ZTAU = PDTAUCIV(JL)
         ZOM = PSSALBIV(JL)
!-----------------------------------------------------------------------
!       ... Stay away from 1 by precision.  for g, also stay away from -1
!-----------------------------------------------------------------------
         ZTEMPG = MIN( ABS(ZG),1. - ZPRECIS )
         ZG = SIGN( ZTEMPG,ZG )
         ZOM = MIN( ZOM,1. - ZPRECIS )
!-----------------------------------------------------------------------
!       ... Delta-scaling
!-----------------------------------------------------------------------
         ZF = ZG*ZG
         ZG = (ZG - ZF) / (1. - ZF)
         ZTAUN(JL) = (1. - ZOM*ZF)*ZTAU
         ZOM = (1. - ZF)*ZOM / (1. - ZOM*ZF)
!-----------------------------------------------------------------------
!       ... The following gamma equations are from pg 16,289, table 1
!-----------------------------------------------------------------------
!       ... Eddington approximation
!-----------------------------------------------------------------------
         ZGAM1 = .25 * (7. - ZOM*(4. + 3.*ZG))
         ZGAM2 = -.25 * (1. - ZOM*(4. - 3.*ZG))
         ZGAM3 = .25 * (2. - 3.*ZG*PMU2(JL))
         ZGAM4 = 1. - ZGAM3
!-----------------------------------------------------------------------
!       ... Hemispheric mean; quadrature
!           save mu1 for use in converting irradiance to actinic flux
!-----------------------------------------------------------------------
         ZMU1(JL) = (1. - ZOM) / (ZGAM1 - ZGAM2)
!-----------------------------------------------------------------------
!       ... lambda = pg 16,290 equation 21
!           big gamma = pg 16,290 equation 22
!-----------------------------------------------------------------------
         ZLAM(JL) = SQRT(ZGAM1*ZGAM1 - ZGAM2*ZGAM2)
         ZBGAM(JL) = (ZGAM1 - ZLAM(JL)) / ZGAM2
         ZEXPON = EXP( -ZLAM(JL)*ZTAUN(JL) )
!-----------------------------------------------------------------------
!       ... e1 - e4 = pg 16,292 equation 44
!-----------------------------------------------------------------------
         ZE1(JL) = 1. + ZBGAM(JL)*ZEXPON
         ZE2(JL) = 1. - ZBGAM(JL)*ZEXPON
         ZE3(JL) = ZBGAM(JL) + ZEXPON
         ZE4(JL) = ZBGAM(JL) - ZEXPON

!-----------------------------------------------------------------------
!       ... The following sets up for the c equations 23, and 24
!           found on page 16,290
!           prevent division by zero
!           (if lambda=1/mu, shift 1/mu^2 by eps = 1.e-3)
!           which is approx equiv to shifting mu by 0.5*eps* (mu)**3
!-----------------------------------------------------------------------
         ZEXPON0 = EXP( -ZTAUC/PMU2(JL) )
         ZEXPON1 = EXP( -(ZTAUC + ZTAUN(JL))/PMU2(JL) )
         ZDIVISR = ZLAM(JL)*ZLAM(JL) - 1./(PMU2(JL)*PMU2(JL))
         ZTEMP = MAX( ZEPS,ABS(ZDIVISR) )
         ZDIVISR = SIGN( ZTEMP,ZDIVISR )

         ZUP = ZOM*ZPIFS*((ZGAM1 - 1./PMU2(JL))*ZGAM3 + ZGAM4*ZGAM2)/ZDIVISR
         ZDN = ZOM*ZPIFS*((ZGAM1 + 1./PMU2(JL))*ZGAM4 + ZGAM2*ZGAM3)/ZDIVISR

!-----------------------------------------------------------------------
!       ... cup and cdn are when tau is equal to zero
!           cuptn and cdntn are when tau is equal to taun
!-----------------------------------------------------------------------
         ZCUP(JL) = ZUP*ZEXPON0
         ZCDN(JL) = ZDN*ZEXPON0
         ZCUPTN(JL) = ZUP*ZEXPON1
         ZCDNTN(JL) = ZDN*ZEXPON1
         ZTAUC = ZTAUC + ZTAUN(JL)
      ENDDO

!-----------------------------------------------------------------------
!       ... Set up matrix
!           ssfc = pg 16,292 equation 37  where ZPIFS is one (unity).
!-----------------------------------------------------------------------
      ZSSFC = PALBEDO*PMU2(KMXCLY)*EXP( -ZTAUC/PMU2(KMXCLY) )*ZPIFS

!-----------------------------------------------------------------------
!       ... The following are from pg 16,292  equations 39 - 43.
!           set up first row of matrix:
!-----------------------------------------------------------------------
      ZA(1) = 0.
      ZB(1) = ZE1(1)
      ZD(1) = -ZE2(1)
      ZY(1) = ZFDN0 - ZCDN(1)

!-----------------------------------------------------------------------
!       ... Set up odd rows 3 thru (IKNROWS - 1):
!-----------------------------------------------------------------------
      JL = 0
      DO IROW = 3,IKNROWS-1,2
         JL = JL + 1
         ZA(IROW) = ZE2(JL)*ZE3(JL) - ZE4(JL)*ZE1(JL)
         ZB(IROW) = ZE1(JL)*ZE1(JL + 1) - ZE3(JL)*ZE3(JL + 1)
         ZD(IROW) = ZE3(JL)*ZE4(JL + 1) - ZE1(JL)*ZE2(JL + 1)
         ZY(IROW) = ZE3(JL)*(ZCUP(JL + 1) - ZCUPTN(JL)) &
     &          + ZE1(JL)*(ZCDNTN(JL) - ZCDN(JL + 1))
      ENDDO

!-----------------------------------------------------------------------
!       ... Set up even rows 2 thru (IKNROWS - 2):
!-----------------------------------------------------------------------
      JL = 0
      DO IROW = 2,IKNROWS-2,2
         JL = JL + 1
         ZA(IROW) = ZE2(JL + 1)*ZE1(JL) - ZE3(JL)*ZE4(JL + 1)
         ZB(IROW) = ZE2(JL)*ZE2(JL + 1) - ZE4(JL)*ZE4(JL + 1)
         ZD(IROW) = ZE1(JL + 1)*ZE4(JL + 1) - ZE2(JL + 1)*ZE3(JL + 1)
         ZY(IROW) = (ZCUP(JL + 1) - ZCUPTN(JL))*ZE2(JL + 1) &
     &            - (ZCDN(JL + 1) - ZCDNTN(JL))*ZE4(JL + 1)
      ENDDO

!-----------------------------------------------------------------------
!       ... Set up last row of matrix at IKNROWS:
!-----------------------------------------------------------------------
      ZA(IKNROWS) = ZE1(KMXCLY) - PALBEDO*ZE3(KMXCLY)
      ZB(IKNROWS) = ZE2(KMXCLY) - PALBEDO*ZE4(KMXCLY)
      ZD(IKNROWS) = 0.
      ZY(IKNROWS) = ZSSFC - ZCUPTN(KMXCLY) + PALBEDO*ZCDNTN(KMXCLY)

!-----------------------------------------------------------------------
!       ... Solve tri-diagonal matrix:
!-----------------------------------------------------------------------
      CALL TRIDLA( IKNROWS, ZA, ZB, ZD, ZY )

!-----------------------------------------------------------------------
!       ... Unfold solution of matrix, compute output fluxes
!-----------------------------------------------------------------------
!       ... The following equations are from pg 16,291  equations 31 & 32
!-----------------------------------------------------------------------
      PFDR(0) = 1.
      PFDN(0) = ZFDN0 / ZMU1(1)
      PFUP(0) =  (ZY(1)*ZE3(1) - ZY(2)*ZE4(1) + ZCUP(1)) / ZMU1(1)

      IROW = 1
      ZTAUG = 0.
      DO ILEV = 1,KMXCLY
         ZTAUG = ZTAUG + ZTAUN(ILEV)
         PFDR(ILEV) = EXP( -ZTAUG/PMU2(ILEV) )
         PFDN(ILEV) =  (ZY(IROW)*ZE3(ILEV) + ZY(IROW+1)*ZE4(ILEV) + ZCDNTN(ILEV)) &
     &               / ZMU1(ILEV)
         PFUP(ILEV) =  (ZY(IROW)*ZE1(ILEV) + ZY(IROW+1)*ZE2(ILEV) + ZCUPTN(ILEV)) &
     &               / ZMU1(ILEV)
         IROW = IROW + 2
      ENDDO

IF (LHOOK) CALL DR_HOOK('TWOSTREAM_PLANE',1,ZHOOK_HANDLE )

END SUBROUTINE TWOSTREAM_PLANE

!=======================================================================

SUBROUTINE TWOSTREAM_SPHERE( KMXCLY, PSZA, PMU2, PALBEDO,               &
                           &   PDTAUCIV, PSSALBIV, PGIV, PDSDH, KNID,    &
                           &   PFDR, PFUP, PFDN )
!-----------------------------------------------------------------------
!       ... Two-stream equations for multiple layers based on equations
!           from Toon et al., JGR, Vol 94, #d13  Nov. 20, 1989
!           To be called separately for each wavelength interval iv
!           A pseudo-spherical correction has been added.
!-----------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... PARAMETERS
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), PARAMETER :: ZEPS = 1.E-3, ZPRECIS = 1.E-7
      REAL(KIND=JPRB), PARAMETER :: ZPIFS = 1.
      REAL(KIND=JPRB), PARAMETER :: ZFDN0 = 0.
      REAL(KIND=JPRB), PARAMETER :: ZLARGEST = HUGE(ZLARGEST)   ! largest number

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
!       mu2 = cosine of solar zenith angle if this angle < 65 degrees
!             ( mu2 independent of altitude ) OR 1/Chapman fct
!             if angle > 65 degrees ( mu2(lev) is then lev-dependent )
!             Using Chapman fct based on scale height H for air
!       albedo = surface albedo
!       PDTAUCIV =  unscaled optical depth of each layer
!       PSSALBIV  =  unscaled single scattering albedo
!       PGIV   =  unscaled asymmetry factor
!     PSZA = solar zenith angle
!     PDSDH = slant path of direct beam through each layer crossed
!            when travelling from the top of the atmosphere to layer i;
!            PDSDH(i,j), i = 0..KMXCLY, j = 1..KMXCLY
!     KNID = INTEGER, number of layers crossed by the direct beam when
!            travelling from the top of the atmosphere to layer i;
!            KNID(i), i = 0..KMXCLY
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), INTENT(IN)  :: PALBEDO, PSZA
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN)  :: PMU2
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,KMXCLY), INTENT(IN)  :: PDSDH
      INTEGER(KIND=JPIM), DIMENSION(0:KMXCLY), INTENT(IN) :: KNID
      REAL(KIND=JPRB), DIMENSION(KMXCLY), INTENT(INOUT) :: PDTAUCIV, PSSALBIV, PGIV
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(OUT) :: PFUP, PFDN, PFDR

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM) :: ILEV, i, j
      INTEGER(KIND=JPIM) :: IROW, IKNROWS
      REAL(KIND=JPRB)  :: ZTEMPG
      REAL(KIND=JPRB)  :: ZSUMTAU
      REAL(KIND=JPRB)  :: ZSSFC
      REAL(KIND=JPRB)  :: ZTAUG
      REAL(KIND=JPRB)  :: ZEXPON, ZEXPON0, ZEXPON1
      REAL(KIND=JPRB)  :: ZDIVISR, ZUP, ZDN, ZTEMP
      REAL(KIND=JPRB)  :: ZGAM1, ZGAM2, ZGAM3, ZGAM4
      REAL(KIND=JPRB)  :: ZOM, ZF, ZG
      REAL(KIND=JPRB), DIMENSION(KMXCLY) :: ZLAM, ZTAUN, ZBGAM, ZE1, ZE2, ZE3, ZE4, &
     &                          ZCUP, ZCDN, ZCUPTN, ZCDNTN, ZMU1
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY) :: ZTAUSLA, ZTAUC, ZMU3
      REAL(KIND=JPRB), DIMENSION(2*KMXCLY) :: ZA, ZB, ZD, ZY  ! 2*KMXCLY = nb of rows in the matrix
      REAL(KIND=JPRB)  :: zero
      REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('TWOSTREAM_SPHERE',0,ZHOOK_HANDLE )


      ZERO = SQRT (TINY(zero))  ! a "safe" tiny real

!-----------------------------------------------------------------------
!       ... Initial conditions:  pi*solar flux = 1;
!                                diffuse incidence = 0
!-----------------------------------------------------------------------
!       ... Compute coefficients for each layer:
!           gam1 - gam4 = 2-stream coefficients,
!           different for different approximations
!           ZEXPON0 = calculation of e when tau is zero
!           ZEXPON1 = calculation of e when tau is taun
!           cup and cdn = calculation when tau is zero
!           cuptn and cdntn = calc. when tau is taun
!           ZDIVISR = prevents division by zero
!-----------------------------------------------------------------------
      IKNROWS = 2*KMXCLY        ! number of rows in the matrix
      ZTAUC(0:KMXCLY) = 0.
      ZTAUSLA(0:KMXCLY) = 0.
      ZMU3(0:KMXCLY) = zero

!-----------------------------------------------------------------------
!       ... Delta-scaling
!-----------------------------------------------------------------------
      DO i = 1,KMXCLY
         ZF = PGIV(i)*PGIV(i)
         PGIV(i) = (PGIV(i) - ZF) / (1._JPRB - ZF)
         ZTAUN(i) = (1._JPRB - PSSALBIV(i)*ZF)*PDTAUCIV(i)
         PSSALBIV(i) = (1._JPRB - ZF)*PSSALBIV(i) / (1._JPRB - PSSALBIV(i)*ZF)
      ENDDO

!-----------------------------------------------------------------------
! Calculate slant optical depth at the top of the atmosphere when zen>90.
! in this case, higher altitude of the top layer is recommended which
! can be easily changed in gridz.f.
!-----------------------------------------------------------------------
      IF( PSZA > 90.E0_JPRB ) THEN
        IF( KNID(0) < 0_JPIM ) THEN
           ZTAUSLA(0) = ZLARGEST
         ELSE
           ZSUMTAU = 0.0_JPRB
           DO j = 1, KNID(0)
              ZSUMTAU = ZSUMTAU + 2.*ZTAUN(j)*PDSDH(0,j)
           ENDDO
           ZTAUSLA(0) = ZSUMTAU
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
!  Loop on layer
!-----------------------------------------------------------------------
      DO i = 1, KMXCLY
         ZG = PGIV(i)
         ZTAUC(i) = ZTAUC(i-1) + ZTAUN(i)
         ZOM = PSSALBIV(i)
!-----------------------------------------------------------------------
!       ... Stay away from 1 by precision.  for g, also stay away from -1
!-----------------------------------------------------------------------
         ZTEMPG = MIN( ABS(ZG), 1.E0_JPRB-ZPRECIS )
         ZG = SIGN( ZTEMPG,ZG )
         ZOM = MIN( ZOM, 1.E0_JPRB-ZPRECIS )
!-----------------------------------------------------------------------
!    ...Calculate slant optical depth
!-----------------------------------------------------------------------
         IF( KNID(i) < 0_JPIM ) THEN
            ZTAUSLA(i) = ZLARGEST
          ELSE
            ZSUMTAU = 0.0
            DO j = 1, MIN(KNID(i),i)
               ZSUMTAU = ZSUMTAU + ZTAUN(j)*PDSDH(i,j)
            ENDDO
            DO j = 1+MIN(KNID(i),i), KNID(i)
               ZSUMTAU = ZSUMTAU + 2.*ZTAUN(j)*PDSDH(i,j)
            ENDDO
            ZTAUSLA(i) = ZSUMTAU
            IF(ZTAUSLA(i) == ZTAUSLA(i-1)) THEN
               ZMU3(i) = SQRT(ZLARGEST)
             ELSE
               ZMU3(i) = (ZTAUC(i)-ZTAUC(i-1))/(ZTAUSLA(i)-ZTAUSLA(i-1))
               ZMU3(i) = SIGN( MAX(ABS(ZMU3(i)),zero), ZMU3(i) )
            ENDIF
         ENDIF

!-----------------------------------------------------------------------
!       ... The following gamma equations are from pg 16,289, table 1
!-----------------------------------------------------------------------
!       ... Eddington approximation
!-----------------------------------------------------------------------
         ZGAM1 =  .25_JPRB * (7._JPRB - ZOM*(4._JPRB + 3._JPRB*ZG))
         ZGAM2 = -.25_JPRB * (1._JPRB - ZOM*(4._JPRB - 3._JPRB*ZG))
         ZGAM3 =  .25_JPRB * (2._JPRB - 3._JPRB*ZG*PMU2(i))
         ZGAM4 = 1._JPRB - ZGAM3
!-----------------------------------------------------------------------
!       ... Hemispheric mean; quadrature
!           save mu1 for use in converting irradiance to actinic flux
!-----------------------------------------------------------------------
         ZMU1(i) = (1._JPRB - ZOM) / (ZGAM1 - ZGAM2)   ! mu1 = 0.5 in JF Muller's version
!-----------------------------------------------------------------------
!       ... lambda = pg 16,290 equation 21
!           big gamma = pg 16,290 equation 22
!-----------------------------------------------------------------------
         ZLAM(i) = SQRT(ZGAM1*ZGAM1 - ZGAM2*ZGAM2)
         ZBGAM(i) = (ZGAM1 - ZLAM(i)) / ZGAM2
         ZEXPON = EXP( -ZLAM(i)*ZTAUN(i) )
!-----------------------------------------------------------------------
!       ... e1 - e4 = pg 16,292 equation 44
!-----------------------------------------------------------------------
         ZE1(i) = 1._JPRB + ZBGAM(i)*ZEXPON
         ZE2(i) = 1._JPRB - ZBGAM(i)*ZEXPON
         ZE3(i) = ZBGAM(i) + ZEXPON
         ZE4(i) = ZBGAM(i) - ZEXPON

!-----------------------------------------------------------------------
!       ... The following sets up for the c equations 23, and 24
!           found on page 16,290
!           prevent division by zero
!           (if lambda=1/mu, shift 1/mu^2 by eps = 1.e-3)
!           which is approx equiv to shifting mu by 0.5*eps* (mu)**3
!-----------------------------------------------------------------------
         ZEXPON0 = EXP( -ZTAUSLA(i-1) )
         ZEXPON1 = EXP( -ZTAUSLA(i) )
         ZDIVISR = ZLAM(i)*ZLAM(i) - 1._JPRB/(ZMU3(i)*ZMU3(i))
         ZTEMP = MAX( ZEPS,ABS(ZDIVISR) )
         ZDIVISR = SIGN( ZTEMP,ZDIVISR )

         ZUP = ZOM*ZPIFS*((ZGAM1 - 1._JPRB/ZMU3(i))*ZGAM3 + ZGAM4*ZGAM2)/ZDIVISR
         ZDN = ZOM*ZPIFS*((ZGAM1 + 1._JPRB/ZMU3(i))*ZGAM4 + ZGAM2*ZGAM3)/ZDIVISR

!-----------------------------------------------------------------------
!       ... cup and cdn are when tau is equal to zero
!           cuptn and cdntn are when tau is equal to taun
!-----------------------------------------------------------------------
         ZCUP(i) = ZUP*ZEXPON0
         ZCDN(i) = ZDN*ZEXPON0
         ZCUPTN(i) = ZUP*ZEXPON1
         ZCDNTN(i) = ZDN*ZEXPON1
      ENDDO

!-----------------------------------------------------------------------
!       ... Set up matrix
!           ssfc = pg 16,292 equation 37  where ZPIFS is one (unity).
!-----------------------------------------------------------------------
      ZSSFC = PALBEDO*PMU2(KMXCLY)*EXP( -ZTAUSLA(KMXCLY) )*ZPIFS

!-----------------------------------------------------------------------
!       ... The following are from pg 16,292  equations 39 - 43.
!           set up first row of matrix:
!-----------------------------------------------------------------------
      ZA(1) = 0.
      ZB(1) = ZE1(1)
      ZD(1) = -ZE2(1)
      ZY(1) = ZFDN0 - ZCDN(1)

!-----------------------------------------------------------------------
!       ... Set up odd rows 3 thru (IKNROWS - 1):
!-----------------------------------------------------------------------
      i = 0
      DO IROW = 3,IKNROWS-1,2
         i = i + 1
         ZA(IROW) = ZE2(i)*ZE3(i) - ZE4(i)*ZE1(i)
         ZB(IROW) = ZE1(i)*ZE1(i + 1) - ZE3(i)*ZE3(i + 1)
         ZD(IROW) = ZE3(i)*ZE4(i + 1) - ZE1(i)*ZE2(i + 1)
         ZY(IROW) = ZE3(i)*(ZCUP(i + 1) - ZCUPTN(i)) &
     &          + ZE1(i)*(ZCDNTN(i) - ZCDN(i + 1))
      ENDDO

!-----------------------------------------------------------------------
!       ... Set up even rows 2 thru (IKNROWS - 2):
!-----------------------------------------------------------------------
      i = 0
      DO IROW = 2,IKNROWS-2,2
         i = i + 1
         ZA(IROW) = ZE2(i + 1)*ZE1(i) - ZE3(i)*ZE4(i + 1)
         ZB(IROW) = ZE2(i)*ZE2(i + 1) - ZE4(i)*ZE4(i + 1)
         ZD(IROW) = ZE1(i + 1)*ZE4(i + 1) - ZE2(i + 1)*ZE3(i + 1)
         ZY(IROW) = (ZCUP(i + 1) - ZCUPTN(i))*ZE2(i + 1) &
     &          - (ZCDN(i + 1) - ZCDNTN(i))*ZE4(i + 1)
      ENDDO

!-----------------------------------------------------------------------
!       ... Set up last row of matrix at IKNROWS:
!-----------------------------------------------------------------------
      ZA(IKNROWS) = ZE1(KMXCLY) - PALBEDO*ZE3(KMXCLY)
      ZB(IKNROWS) = ZE2(KMXCLY) - PALBEDO*ZE4(KMXCLY)
      ZD(IKNROWS) = 0.
      ZY(IKNROWS) = ZSSFC - ZCUPTN(KMXCLY) + PALBEDO*ZCDNTN(KMXCLY)

!-----------------------------------------------------------------------
!       ... Solve tri-diagonal matrix:
!-----------------------------------------------------------------------
      CALL TRIDLA( IKNROWS, ZA, ZB, ZD, ZY )

!-----------------------------------------------------------------------
!       ... Unfold solution of matrix, compute output fluxes
!-----------------------------------------------------------------------
!       ... The following equations are from pg 16,291  equations 31 & 32
!-----------------------------------------------------------------------
      PFDR(0) = EXP(-ZTAUSLA(0))
      PFDN(0) = ZFDN0 / ZMU1(1)
      PFUP(0) =  (ZY(1)*ZE3(1) - ZY(2)*ZE4(1) + ZCUP(1)) / ZMU1(1)

      IROW = 1
      ZTAUG = 0.
      DO ILEV = 1,KMXCLY
         PFDR(ILEV) = EXP( -ZTAUSLA(ILEV) )
         PFDN(ILEV) =  (ZY(IROW)*ZE3(ILEV) + ZY(IROW+1)*ZE4(ILEV) + ZCDNTN(ILEV)) &
     &               / ZMU1(ILEV)
         PFUP(ILEV) =  (ZY(IROW)*ZE1(ILEV) + ZY(IROW+1)*ZE2(ILEV) + ZCUPTN(ILEV)) &
     &               / ZMU1(ILEV)
         IROW = IROW + 2
      ENDDO

IF (LHOOK) CALL DR_HOOK('TWOSTREAM_SPHERE',1,ZHOOK_HANDLE )

END SUBROUTINE TWOSTREAM_SPHERE

!=======================================================================

SUBROUTINE ABSDEP( KMXCLY, PCRSTIV, PDENS, PZS, PBABS )
!-----------------------------------------------------------------------
!     ABSorption optical DEPth - Calculates optical depth due to
!        absorption by all the absorbing species
!        PBABS(NO) negligible - PDENS(NO) used only for J(NO)
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_MODULE, only : jprocmax => ndiss, j_O2_O, j_o3_o, j_no2, j_co2
USE BASCOE_TUV_MODULE, only : nabspec, O2abs, O3abs, NO2abs, CO2abs

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), DIMENSION(jprocmax,0:KMXCLY), INTENT(IN) ::  PCRSTIV
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,nabspec), INTENT(IN)  ::  PDENS
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN):: PZS       ! geometric altitude (km)
      REAL(KIND=JPRB), DIMENSION(KMXCLY,2:nabspec), INTENT(OUT) :: PBABS   ! Total Abs.coeff.

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), DIMENSION(KMXCLY) :: ZDELTAZ
      REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ABSDEP',0,ZHOOK_HANDLE )

      PBABS(1:KMXCLY,2:NABSPEC) = 0._JPRB
      ZDELTAZ(1:KMXCLY) = 5.e4_JPRB * (PZS(0:KMXCLY-1) - PZS(1:KMXCLY))    ! (cm)
      PBABS(1:KMXCLY,O2abs) = ZDELTAZ &
     &               * ( PDENS(0:KMXCLY-1,O2abs)*PCRSTIV(J_o2_o,0:KMXCLY-1) &
     &                 + PDENS(1:KMXCLY,O2abs)*PCRSTIV(J_o2_o,1:KMXCLY) )
      PBABS(1:KMXCLY,O3abs) = ZDELTAZ &
     &               * ( PDENS(0:KMXCLY-1,O3abs)*PCRSTIV(J_o3_o,0:KMXCLY-1) &
     &                 + PDENS(1:KMXCLY,O3abs)*PCRSTIV(J_o3_o,1:KMXCLY) )
      PBABS(1:KMXCLY,NO2abs) = ZDELTAZ &
     &               * ( PDENS(0:KMXCLY-1,NO2abs)*PCRSTIV(J_NO2,0:KMXCLY-1) &
     &                  + PDENS(1:KMXCLY,NO2abs)*PCRSTIV(J_NO2,1:KMXCLY) )
      PBABS(1:KMXCLY,CO2abs) = ZDELTAZ &
     &               * ( PDENS(0:KMXCLY-1,CO2abs)*PCRSTIV(j_co2,0:KMXCLY-1) &
     &                 + PDENS(1:KMXCLY,CO2abs)*PCRSTIV(j_co2,1:KMXCLY) )

IF (LHOOK) CALL DR_HOOK('ABSDEP',1,ZHOOK_HANDLE )

END SUBROUTINE ABSDEP

!=======================================================================

SUBROUTINE GET_LAYDENS( KMXCLY, PZS, PDENS, PLAYDENS  )
!-----------------------------------------------------------------------
!       ...Obtain PLAYDENS, the densities in each horizontal layer
!          using exponential or linear interpolation
!          The factor 1.E+5 converts PZS from km to cm
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE, only : nabspec

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... PARAMETERS
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), PARAMETER :: ZEPSILN = 1.E-5

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN) :: PZS              ! (km)
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,nabspec), INTENT(IN) :: PDENS
      REAL(KIND=JPRB), DIMENSION(KMXCLY,nabspec), INTENT(OUT)  :: PLAYDENS

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM) :: i, ILC
      REAL(KIND=JPRB)  :: ZDELTAZ
      REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GET_LAYDENS',0,ZHOOK_HANDLE )

      DO ILC = 1,KMXCLY
         ZDELTAZ = (PZS(ILC-1) - PZS(ILC)) * 1.E5_JPRB
         DO i = 1, nabspec
            IF( PDENS(ILC-1,i) > 0. .and. PDENS(ILC,i) > 0. .and. &
     &                ABS(1.0_JPRB - PDENS(ILC,i)/PDENS(ILC-1,i)) > ZEPSILN ) THEN
               PLAYDENS(ILC,i) = 1._JPRB / (LOG(PDENS(ILC,i) / PDENS(ILC-1,i))) * &
     &                        (PDENS(ILC,i) - PDENS(ILC-1,i)) * ZDELTAZ
             ELSE
               PLAYDENS(ILC,i) = &
     &              .5 * ABS( (PDENS(ILC-1,i) + PDENS(ILC,i) )*ZDELTAZ )
            ENDIF
         ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('GET_LAYDENS',1,ZHOOK_HANDLE )

END SUBROUTINE GET_LAYDENS

!!=======================================================================

SUBROUTINE LYA( KMXCLY, PCOLO2X, PRMLYA, PRO2LYA )
!------------------------------------------------------------------
!       ... Calculation of reduction factors Rm and Ro2 for
!           Lyman-alpha line: Chabrillat & Kockarts, GRL, Vol 24, p2659, 1997
!------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE, only : sw_lya

      IMPLICIT NONE
!------------------------------------------------------------------
!       ... Dummy args
!------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), INTENT(IN)  :: PCOLO2X(0:KMXCLY) ! slant overhead O2 column (molec/cm2)
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(OUT) :: PRMLYA, PRO2LYA

!------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: ic
      REAL(KIND=JPRB), save :: ZB(3) = (/ 0.392161, 0.380645, 0.226447 /)
      REAL(KIND=JPRB), save :: ZC(3) = (/8.30842e-21, 1.41793e-20, 6.61628e-21/)
      REAL(KIND=JPRB), save :: ZD(3) = (/3.61207e-21, 5.08214e-21, 1.60035e-21/)
      REAL(KIND=JPRB), save :: ZE(3) = (/8.54654e-21, 1.50741e-20, 6.63640e-21/)

IF (LHOOK) CALL DR_HOOK('LYA',0,ZHOOK_HANDLE )

!------------------------------------------------------------------
!       ... Calculation of reduction factors Rm and Ro2 for
!           Lyman-alpha line
!           (Chabrillat and Kockarts, submitted to GRL, 1997)
!------------------------------------------------------------------
      PRMLYA = 0.
      PRO2LYA = 0.
      IF( sw_lya ) THEN
         DO ic = 1, 3
            PRMLYA(0:KMXCLY) = PRMLYA(0:KMXCLY) + ZB(ic) * EXP( -ZC(ic) * PCOLO2X(0:KMXCLY) )
            PRO2LYA(0:KMXCLY) = PRO2LYA(0:KMXCLY) + ZD(ic) * EXP( -ZE(ic) * PCOLO2X(0:KMXCLY) )
         ENDDO
      ENDIF

IF (LHOOK) CALL DR_HOOK('LYA',1,ZHOOK_HANDLE )

END SUBROUTINE LYA

!=======================================================================

SUBROUTINE SOLHEAT_EFFIC( KMXCLY, PPMB, PEFSHO2, PEFSHO3 )
!-----------------------------------------------------------------------
!       ... Heating efficiencies  for O3 Hartley band
!           (Mlynczak and Solomon 1993, p10525) and for O2 SR continuum
!           (p 10527 M&S 1993)
!           Corrected from EFFIC in eff.f in SOCRATES up to v1.8c
!           Vertical grid goes down as in PHODIS (PPMB increases with lev)
!                - validity on pmb intervals CORRECTED, since
!                  efficiencies set to last found value above 1e-4 mb.
!       ... added calculation of denerg1 and denerg2
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE, only : MXWVN, iv_lya

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... PARAMETERS
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), PARAMETER ::  ZC(4) = &
     &            (/ .92621, .133960, -.076863, .006897 /)
      REAL(KIND=JPRB), PARAMETER :: ZB(4) = &
     &            (/ .75349, .0036, .059468, -.022795 /)
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN) ::  PPMB  ! depends on dens and temper !
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,MXWVN), INTENT(OUT) :: PEFSHO2, PEFSHO3

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: IV, ILEV1M4
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY) :: ZX, ZEPS_O2, ZEPS_O3

IF (LHOOK) CALL DR_HOOK('SOLHEAT_EFFIC',0,ZHOOK_HANDLE )

      PEFSHO2(:,:) = 1._JPRB
      PEFSHO3(:,:) = 1._JPRB

      WHERE( PPMB > 1._JPRB )
         ZEPS_O3(:) = 1._JPRB
         ZEPS_O2(:) = 1._JPRB
      END WHERE
      WHERE( PPMB <= 1._JPRB .and. PPMB > 1.e-2_JPRB )
         ZX(:) = LOG10( PPMB(:) ) + 1._JPRB
         ZEPS_O3(:) = ZC(1) + ZX*(ZC(2) + ZX*(ZC(3) + ZX*ZC(4)))
         ZEPS_O2(:) = 1._JPRB
      END WHERE
      WHERE( PPMB <= 1.E-2_JPRB .and. PPMB > 1.e-4_JPRB )
         ZX(:) = LOG10( PPMB(:) ) + 3._JPRB
         ZEPS_O3(:) = ZC(1) + ZX*(ZC(2) + ZX*(ZC(3) + ZX*ZC(4)))
         ZEPS_O2(:) = ZB(1) + ZX*(ZB(2) + ZX*(ZB(3) + ZX*ZB(4)))
      END WHERE
      DO ILEV1M4 = KMXCLY, 0, -1
         IF( PPMB(ILEV1M4) <= 1.E-4_JPRB ) EXIT
      ENDDO
      ZEPS_O3(0:ILEV1M4) = ZEPS_O3(ILEV1M4+1)
      ZEPS_O2(0:ILEV1M4) = ZEPS_O2(ILEV1M4+1)

!-----------------------------------------------------------------
!       ... Heating efficiency for O3 Hartley band
!           (Mlynczak and Solomon 1993, p10525)
!-----------------------------------------------------------------
      DO IV = 60,95
         PEFSHO3(:,IV) = ZEPS_O3(:)
      ENDDO

!-----------------------------------------------------------------
!       ... Heating  efficiency for O2 SR continuum (p 10527 M&S 1993)
!-----------------------------------------------------------------
      DO IV = 28,50
         PEFSHO2(:,IV) = ZEPS_O2(:)
      ENDDO

!-----------------------------------------------------------------
!       ... Heating  efficiency for O2 Lyman alpha band (IV = 8)
!-----------------------------------------------------------------
      PEFSHO2(:,iv_lya) = .95_JPRB

IF (LHOOK) CALL DR_HOOK('SOLHEAT_EFFIC',1,ZHOOK_HANDLE )

END SUBROUTINE SOLHEAT_EFFIC

!=======================================================================

SUBROUTINE SOLHEAT_RATES( KMXCLY, PZS, PCRST, PUAVG, PEFSHO2, PEFSHO3, & ! input
                        &   PWRAT )                                 ! output
!-----------------------------------------------------------------------
!     SUBROUTINE SOLar HEATing RATes
!     extracted and enhanced from PHORAT in SOCRATES up to v1.8c
!
!     Calculates the solar UV heating rates from
!     layer 0 (top layer) to layer KMXCLY (surface layer)
!
!     PWRAT(1,0:KMXCLY) : solar heating rate from O2 absorption (K/s)
!     PWRAT(2,0:KMXCLY) : solar heating rate from O3 absorption (K/s)
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULOUT
USE BASCOE_J_MODULE, only : jprocmax => ndiss, J_o2_o, J_o2_O1D, J_o3_o, J_o3_O1D
USE BASCOE_TUV_MODULE, only : MXWVN, nheatspec, ivbegin, ivend, hv, denerg_O2, denerg_O3

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN) :: PZS
      REAL(KIND=JPRB), DIMENSION(jprocmax,0:KMXCLY,MXWVN), INTENT(IN) :: PCRST
      REAL(KIND=JPRB), DIMENSION(MXWVN,0:KMXCLY), INTENT(IN) :: PUAVG
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,MXWVN), INTENT(IN) :: PEFSHO2, PEFSHO3
      REAL(KIND=JPRB), DIMENSION(nheatspec,0:KMXCLY), INTENT(OUT) ::  PWRAT

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      CHARACTER(LEN=*), PARAMETER    :: CL_MY_NAME   = 'SOLHEAT_RATES'
      REAL(KIND=JPRB), PARAMETER :: ZPIR = 3.1415926535897931E0 ! pi
      REAL(KIND=JPRB), PARAMETER :: ZFOURPI = 4. * ZPIR

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: IK, IWB, IWE, IWAVE
      INTEGER(KIND=JPIM) :: ILEV60
      REAL(KIND=JPRB), DIMENSION(jprocmax,0:KMXCLY,MXWVN) :: ZXSECT   ! NOTE memory waste: will use only 4 among jprocmax

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SOLHEAT_RATES',0,ZHOOK_HANDLE )

      DO ILEV60 = 0, KMXCLY             ! find lev60, level index of 1st
         IF( PZS(ILEV60) < 60._JPRB ) EXIT   ! altitude below 60km
      ENDDO
      IF( ILEV60 == 0_JPIM ) THEN
         WRITE(NULOUT,*) CL_MY_NAME//': ILEV60 not found'
         CALL ABOR1( CL_MY_NAME//': ILEV60 not found' )
      ENDIF

      PWRAT(:,:) = 0._JPRB
      ZXSECT(:,:,:) = -999._JPRB

!---------------------------------------------------------------------
!       ... Prepare cross section for heating calculation
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!       ... o2
!---------------------------------------------------------------------

      DO IWAVE = ivbegin(J_o2_o),ivend(J_o2_o)
         ZXSECT(J_o2_o,0:ILEV60-1,IWAVE)      = PCRST(J_o2_o,0:ILEV60-1,IWAVE) &
     &                                    * denerg_O2(IWAVE) &
     &                                    * PEFSHO2(0:ILEV60-1,IWAVE)
         ZXSECT(J_o2_o,ILEV60:KMXCLY,IWAVE)    = PCRST(J_o2_o,ILEV60:KMXCLY,IWAVE) &
     &                                    * hv(IWAVE)
      ENDDO

      IF( any( ZXSECT(J_o2_o,0:ILEV60-1,ivbegin(J_o2_o):ivend(J_o2_o)) < 0. )) THEN
         CALL ABOR1( CL_MY_NAME//': bug found, o2_o(1a)' )
      ENDIF
      IF( any( ZXSECT(J_o2_o,ILEV60:KMXCLY,ivbegin(J_o2_o):ivend(J_o2_o)) < 0. )) THEN
         CALL ABOR1( CL_MY_NAME//': bug found, o2_o(1b)' )
      ENDIF

      DO IWAVE = ivbegin(J_o2_O1D),ivend(J_o2_O1D)
         ZXSECT(J_o2_O1D,0:ILEV60-1,IWAVE)    = PCRST(J_o2_O1D,0:ILEV60-1,IWAVE) &
     &                                    * denerg_O2(IWAVE) &
     &                                    * PEFSHO2(0:ILEV60-1,IWAVE)
         ZXSECT(J_o2_O1D,ILEV60:KMXCLY,IWAVE)  = PCRST(J_o2_O1D,ILEV60:KMXCLY,IWAVE) &
     &                                    * hv(IWAVE)
      ENDDO

!---------------------------------------------------------------------
!       ... o3
!---------------------------------------------------------------------
      DO IWAVE = ivbegin(J_o3_o),ivend(J_o3_o)
         ZXSECT(J_o3_o,0:ILEV60-1,IWAVE)     = PCRST(J_o3_o,0:ILEV60-1,IWAVE) &
     &                                   * denerg_O3(IWAVE) &
     &                                   * PEFSHO3(0:ILEV60-1,IWAVE)
         ZXSECT(J_o3_o,ILEV60:KMXCLY,IWAVE)   = PCRST(J_o3_o,ILEV60:KMXCLY,IWAVE) &
     &                                   * hv(IWAVE)

      ENDDO

      DO IWAVE = ivbegin(J_o3_O1D),ivend(J_o3_O1D)
         ZXSECT(J_o3_O1D,0:ILEV60-1,IWAVE)   = PCRST(J_o3_O1D,0:ILEV60-1,IWAVE) &
     &                                   * denerg_O3(IWAVE) &
     &                                   * PEFSHO3(0:ILEV60-1,IWAVE)
         ZXSECT(J_o3_O1D,ILEV60:KMXCLY,IWAVE) = PCRST(J_o3_O1D,ILEV60:KMXCLY,IWAVE) &
     &                                   * hv(IWAVE)
      ENDDO

!---------------------------------------------------------------------
!       ... o2 heating
!---------------------------------------------------------------------
      IWB = ivbegin(J_o2_o)
      IWE = ivend(J_o2_o)
      DO IK = 0,KMXCLY
         PWRAT(1,IK) = DOT_PRODUCT( ZXSECT(J_o2_o,IK,IWB:IWE), &
     &                            PUAVG(IWB:IWE,IK) )
      ENDDO
      IWB = ivbegin(J_o2_O1D)
      IWE = ivend(J_o2_O1D)
      DO IK = 0,KMXCLY
         PWRAT(1,IK) = DOT_PRODUCT( ZXSECT(J_o2_O1D,IK,IWB:IWE), &
     &                            PUAVG(IWB:IWE,IK) ) &
     &             + PWRAT(1,IK)

      ENDDO

!---------------------------------------------------------------------
!       ... o3 heating
!---------------------------------------------------------------------

      IWB = ivbegin(J_o3_o)
      IWE = ivend(J_o3_o)
      DO IK = 0,KMXCLY
         PWRAT(2,IK) = DOT_PRODUCT( ZXSECT(J_o3_o,IK,IWB:IWE), &
     &                            PUAVG(IWB:IWE,IK) )
      ENDDO

      IWB = ivbegin(J_o3_O1D)
      IWE = ivend(J_o3_O1D)
      DO IK = 0,KMXCLY
         PWRAT(2,IK) = DOT_PRODUCT( ZXSECT(J_o3_O1D,IK,IWB:IWE), &
     &                            PUAVG(IWB:IWE,IK) ) &
     &             + PWRAT(2,IK)

      ENDDO

!---------------------------------------------------------------------
!           Multiply by 4pi because of mean intensity -> actinic flux
!---------------------------------------------------------------------
      PWRAT = ZFOURPI * PWRAT

IF (LHOOK) CALL DR_HOOK('SOLHEAT_RATES',1,ZHOOK_HANDLE )

END SUBROUTINE SOLHEAT_RATES

!=======================================================================

SUBROUTINE SPHERL( KMXCLY, PZS, PSZA, PDSDH, KNID )
!-----------------------------------------------------------------------------
!  PURPOSE:
!  Calculate slant path over vertical depth ds/dh in spherical geometry.
!  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model
!  for computing the radiation field available for photolysis and heating
!  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)
!-----------------------------------------------------------------------------
!  PARAMETERS:
!  NZ      - INTEGER, number of specified altitude levels in the working (I)
!            grid
!  Z       - real*8, specified altitude working grid (km)                  (I)
!  PSZA     - real*8, solar zenith angle (degrees)                          (I)
!  PDSDH    - real*8, slant path of direct beam through each layer crossed  (O)
!            when travelling from the top of the atmosphere to layer i;
!            PDSDH(i,j), i = 0..KMXCLY, j = 1..KMXCLY
!  KNID     - INTEGER, number of layers crossed by the direct beam when   (O)
!            travelling from the top of the atmosphere to layer i;
!            KNID(i), i = 0..KMXCLY
!-----------------------------------------------------------------------------
!  EDIT HISTORY: adapted to writej application, franch@oma.be ;
!                translated to fortran90, simonc@oma.be          nov 1999
!-----------------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), INTENT(IN) :: PSZA
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN) :: PZS
      INTEGER(KIND=JPIM), DIMENSION(0:KMXCLY), INTENT(OUT) :: KNID
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,KMXCLY), INTENT(OUT) :: PDSDH

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), PARAMETER :: ZPIR = 3.1415926535897931E0
      REAL(KIND=JPRB), PARAMETER :: ZD2R = ZPIR/180.E0        ! degree to radian
      REAL(KIND=JPRB), PARAMETER :: ZR0  = 6371.0E0          ! Effective earth radius (km)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JPRB) :: ZSZARAD, ZRPSINZ, ZRJ, ZRJP1, ZDSJ, ZDHJ, ZGA, ZGB, ZSM
      INTEGER(KIND=JPIM) :: i, j, id


IF (LHOOK) CALL DR_HOOK('SPHERL',0,ZHOOK_HANDLE )
      ZSZARAD = PSZA * ZD2R
      KNID(:) = 0
      PDSDH(:,:) = 0.

!-----------------------------------------------------------------------
!        ... calculate ds/dh of every layer
!-----------------------------------------------------------------------
      DO i = 0, KMXCLY

        ZRPSINZ = (ZR0 + PZS(i)) * SIN(ZSZARAD)

        IF ( (PSZA > 90.0E0) .and. (ZRPSINZ < ZR0) ) THEN
           KNID(i) = -1
         ELSE

!-----------------------------------------------------------------------
!        ... Find index of layer in which the screening height lies
!-----------------------------------------------------------------------
           id = i
           IF( PSZA > 90.E0 ) THEN
              DO j = 1, KMXCLY
                 IF( (ZRPSINZ < ( PZS(j-1) + ZR0 ) ) .and. &
     &               (ZRPSINZ >= ( PZS(j) + ZR0 )) ) id = j
              ENDDO
           ENDIF

           DO j = 1, id

             ZSM = 1.0
             IF( j == id .and. id == i .and. PSZA > 90.0 ) ZSM = -1.
             ZRJ = ZR0 + PZS(j-1)
             ZRJP1 = ZR0 + PZS(j)
             ZDHJ = PZS(j-1) - PZS(j)
             ZGA = ZRJ*ZRJ - ZRPSINZ*ZRPSINZ
             ZGB = ZRJP1*ZRJP1 - ZRPSINZ*ZRPSINZ
             IF (ZGA < 0.0) ZGA = 0.0
             IF (ZGB < 0.0) ZGB = 0.0

             IF(id > i .and. j == id) THEN
                ZDSJ = SQRT( ZGA )
              ELSE
                ZDSJ = SQRT( ZGA ) - ZSM*SQRT( ZGB )
             ENDIF

             PDSDH(i,j) = ZDSJ / ZDHJ

           ENDDO

           KNID(i) = id

        ENDIF

      ENDDO
IF (LHOOK) CALL DR_HOOK('SPHERL',1,ZHOOK_HANDLE )

END SUBROUTINE SPHERL

!=======================================================================

SUBROUTINE SRB( KMXCLY, PCOLO2X, PRM, PRO2 )
!-----------------------------------------------------------------------
!       ... Calculate the of reduction factors Rm and Ro2 for
!           Schumann-Runge bands: Kockarts, Ann Geophys, 12, p1207, 1994
!           The 16 wavelength intervals correspond to iv=46-61 in wlmid(mxwvn)
!       ... The necessary PARAMETERS ako,bko MUST have been read
!           from 'srb_kockarts94.dat' by SRB_READ
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE, only : nwvn_srb, ako, bko

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), INTENT(IN) :: PCOLO2X(0:KMXCLY) ! slant overhead O2 column (molec/cm2)
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,nwvn_srb), INTENT(OUT) :: PRM, PRO2

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: iv, ic

IF (LHOOK) CALL DR_HOOK('SRB',0,ZHOOK_HANDLE )

      PRM = 0.
      PRO2 = 0.
      DO iv = 1,nwvn_srb
         DO ic = 1,11,2
            PRM(:,iv) = PRM(:,iv) &
     &               + ako(ic,iv)*EXP( -ako(ic+1,iv)*PCOLO2X(:) )
            PRO2(:,iv) = PRO2(:,iv) &
     &                + bko(ic,iv)*EXP( -bko(ic+1,iv)*PCOLO2X(:) )
         ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('SRB',1,ZHOOK_HANDLE )

END SUBROUTINE SRB

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! from original BASCOE MODULE TUV_MIDATM_CRS (specific for each bascoe chem release)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE CRS_CFC( KMXCLY, PCRST, PTEMPER, PCRS )
!---------------------------------------------------------------------------
!       ... Calculation of cross-sections for halogens (mainly CFCs)
!     as a function of temperature. Source: P.Simon et al., Gillotay et al.
!           These sources are used and documented fully in JPL2002 where
!           temperature range is reported as 210 - 300 K
!---------------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_MODULE, only : jprocmax => ndiss, &
     &                    J_ccl4, J_ch3ccl3, J_ch3cl, J_cfc113, &
     &                    J_hcfc22, J_ha1211, J_ha1301, J_ch3br, J_cfc114, &
     &                    J_cfc11, J_ch2br2
USE BASCOE_TUV_MODULE, only : MXWVN

      IMPLICIT NONE
!---------------------------------------------------------------------------
!       ... Dummy args
!---------------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) ::  KMXCLY
      REAL(KIND=JPRB), INTENT(IN)    ::  PTEMPER(0:KMXCLY)
      REAL(KIND=JPRB), INTENT(IN)    ::  PCRS(jprocmax,MXWVN)
      REAL(KIND=JPRB), INTENT(OUT)   ::  PCRST(jprocmax,0:KMXCLY,MXWVN)
!-----------------------------------------------------------------------
!       ... PARAMETERS
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), PARAMETER  :: ZTMIN = 210.0, ZTMAX = 300.0                      ! T range given by JPL2002
      INTEGER(KIND=JPIM), PARAMETER :: IKNCFC = 11
      INTEGER(KIND=JPIM), PARAMETER :: ILAMIN(IKNCFC) = (/ 60,50,45,52,45,60,42,62,44,45,65 /)
      INTEGER(KIND=JPIM), PARAMETER :: ILAMAX(IKNCFC) = (/ 79,76,67,72,61,93,88,88,68,73,91 /)
      INTEGER(KIND=JPIM), PARAMETER :: I_CFC_INDEX(IKNCFC) =                & ! NAMES PHOTPROCESSES INVOLVED
     &                 (/ J_ccl4, J_ch3ccl3, J_ch3cl, J_cfc113,             & ! NOTE: J_ccl4 was named J_cfc10 until 20091019
     &                    J_hcfc22, J_ha1211, J_ha1301, J_ch3br, J_cfc114,  &
     &                    J_cfc11, J_ch2br2 /)
      REAL(KIND=JPRB), PARAMETER :: ZACFC(IKNCFC,5) = RESHAPE( &
                             &  SOURCE=(/-37.104170, -5.821802E-1, 9.997399E-3, -4.676527E-5, 6.850102E-8, & ! ccl4
                                    &    341.085191, -7.273362,    5.498387E-2, -1.827552E-4, 2.238640E-7, & ! ch3ccl3
                                    &   -299.796165,  5.104685,   -3.363002E-2,  9.580545E-5,-1.013456E-7, & ! ch3cl
                                    &  -1087.881207, 20.004100,   -1.391989E-1,  4.282793E-4,-4.938351E-7, & ! cfc113
                                    &   -106.029241,  1.503771,   -8.247614E-3,  1.420607E-5, 0.,          & ! hcfc22
                                    &   -134.797197,  1.708389,   -9.153990E-3,  2.164407E-5,-1.986293E-8, & ! ha1211
                                    &     62.563060, -2.006832,    1.659204E-2, -5.646547E-5, 6.745870E-8, & ! ha1301
                                    &     46.52,     -1.457962,    1.146929E-2, -3.762666E-5, 4.326408E-8, & ! ch3br
                                    &   -160.495098,  2.480670,   -1.520180E-2,  3.841242E-5,-3.437259E-8, & ! cfc114
                                    &    -84.611,     7.9551E-1,  -2.055E-3,    -4.4812E-6,   1.5838E-8,   & ! cfc11
                                    &    -70.217,     1.9403E-1,   2.7262E-3,   -1.6955E-5,   2.500E-8 /), & ! ch2br2
                                    &   SHAPE=(/IKNCFC,5/), ORDER=(/2,1/))

      REAL(KIND=JPRB), PARAMETER :: ZBCFC(IKNCFC,5) = RESHAPE( &
                             &  SOURCE=(/ 1.073919,  -1.627543E-2, 8.814085E-5,-1.981057E-7, 1.502234E-10, & ! ccl4
                                    &    -1.660090,   3.079969E-2,-2.106719E-4, 6.264984E-7,-6.781342E-10, & ! ch3ccl3
                                    &    -7.172742,   1.483679E-1,-1.146290E-3, 3.918805E-6,-4.999362E-9,  & ! ch3cl
                                    &    12.493465,  -2.393714E-1, 1.714214E-3,-5.439298E-6, 6.454833E-9,  & ! cfc113
                                    &    -1.339882E-1,2.740485E-3,-1.802848E-5, 3.8504E-8,   0.,           & ! hcfc22
                                    &    3.306975E-1,-5.095714E-3, 2.936073E-5,-7.619773E-8, 7.682522E-11, & ! ha1211
                                    &   -9.175482E-1, 1.857479E-2,-1.385710E-4, 4.506561E-7,-5.380311E-10, & ! ha1301
                                    &    9.340858E-1,-1.688734E-2, 1.148689E-4,-3.488086E-7, 3.994462E-10, & ! ch3br
                                    &   -1.529573,    3.524763E-2,-2.995072E-4, 1.112950E-6,-1.525877E-9,  & ! cfc114
                                    &   -5.7912,      1.1689E-1,  -8.8069E-4,   2.9335E-6,  -3.6421E-9,    & ! cfc11
                                    &    2.8993,     -4.3277E-2,   2.3916E-4,  -5.8075E-7,   5.2449E-10 /),& ! ch2br2
                                    &  SHAPE=(/IKNCFC,5/), ORDER=(/2,1/))
!---------------------------------------------------------------------------
!       ... Local variables
!---------------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM)   :: iv, ICFC, iproc
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY) :: ZT
      REAL(KIND=JPRB) ::       ZWFAC1, ZWFAC2
!-----------------------------------------------------------------------
!       ... other PARAMETERS
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), DIMENSION(MXWVN), PARAMETER :: &
               &      zlambd = (/ ( 1.e3/(8.021 - 5.054e-2*FLOAT(iv)), iv = 1,MXWVN ) /)

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('CRS_CFC',0,ZHOOK_HANDLE )

      ZT = MAX( ZTMIN, MIN( PTEMPER, ZTMAX ) )     ! T range given by JPL2002

      DO ICFC = 1,IKNCFC
         iproc = I_CFC_INDEX(ICFC)
         DO iv = 1,MXWVN
            PCRST(iproc,:,iv) = PCRS(iproc,iv)
         ENDDO

         IF( ILAMIN(ICFC) < 42_JPIM .or. ILAMAX(ICFC) > 102_JPIM ) THEN
            CALL ABOR1('CRS_CFC fatal error: wl range out of bound')
         ENDIF  ! zlambd formula represents wl only for intervals 42 to 102, see 'crs.xls'
         DO iv = ILAMIN(ICFC),ILAMAX(ICFC)
            ZWFAC1 = ZACFC(ICFC,1) &
     &            + zlambd(iv)*(ZACFC(ICFC,2) &
     &            + zlambd(iv)*(ZACFC(ICFC,3) &
     &            + zlambd(iv)*(ZACFC(ICFC,4) &
     &            + zlambd(iv)*ZACFC(ICFC,5))))
            ZWFAC2 = ZBCFC(ICFC,1) &
     &            + zlambd(iv)*(ZBCFC(ICFC,2) &
     &            + zlambd(iv)*(ZBCFC(ICFC,3) &
     &            + zlambd(iv)*(ZBCFC(ICFC,4) &
     &            + zlambd(iv)*ZBCFC(ICFC,5))))
            PCRST(iproc,:,iv)  = 10.E0_JPRB **(ZWFAC1 + (ZT(:) - 273._JPRB)*ZWFAC2)
         ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('CRS_CFC',1,ZHOOK_HANDLE )

END SUBROUTINE CRS_CFC

!=======================================================================

SUBROUTINE CRS_TDEP( KMXCLY, PTEMPER, PCRST )
!-----------------------------------------------------------------------
!       ... Correct absorption cross-sections for T dependence
!----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_MODULE, only : jprocmax => ndiss, &
    &                       J_o3_o, J_o3_O1D, J_co2, J_n2o, J_no2, J_hno3, &
    &                       J_cfc12, J_clono2_clo, J_clono2_cl, J_ch2o_hco, J_ch2o_co, &
    &                       J_n2o5, J_h2o2_oh, J_h2o2_ho2, J_ccl4, J_chbr3, J_cl2, J_clo, J_ocs
USE BASCOE_TUV_MODULE, only : MXWVN, wlmid, tb_o3, tc_o3, crs200co2, crs370co2, &
    &                       ta_no2, tb_hno3, ta1_clono2, ta2_clono2, &
    &                       gamma_ch2o, crs295ocs, crs, ivbegin, ivend

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN)    :: PTEMPER
      REAL(KIND=JPRB), DIMENSION(jprocmax,0:KMXCLY,MXWVN), INTENT(OUT) :: PCRST

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: iproc, IN, iv
      REAL(KIND=JPRB) :: ZASUM, ZBSUM, ZSIGMA
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY) :: ZLOG_Y0, ZLOG_Y1, ZTDIFF, ZWORK

      REAL(KIND=JPRB), PARAMETER :: ZA_N2O(0:4) = (/ 68.21023, -4.071805, 4.301146E-02, &
     &                                              -1.777846E-04, 2.520672E-07 /)
      REAL(KIND=JPRB), PARAMETER :: ZB_N2O(0:3) = (/ 123.4014, -2.116255, 1.111572E-02, &
     &                                              -1.881058E-05 /)
      REAL(KIND=JPRB), PARAMETER :: Z_AH2O2(0:7) = &
     &                      (/ 6.4761E+04, -9.2170972E+02, 4.535649, &
     &                          -4.4589016E-03, -4.035101E-05, 1.6878206E-07, &
     &                          -2.652014E-10, 1.5534675E-13 /)
      REAL(KIND=JPRB), PARAMETER :: Z_BH2O2(0:4) = (/ 6.8123E+03, -5.1351E+01, 1.1522E-01, &
     &                              -3.0493E-05, -1.0924E-07 /)

IF (LHOOK) CALL DR_HOOK('CRS_TDEP',0,ZHOOK_HANDLE )

!-----------------------------------------------------------------------
!       ... Only the cross-sections recomputed below are T-dependent
!-----------------------------------------------------------------------
      DO  iv = 1,MXWVN
        DO iproc = 1,jprocmax
          PCRST(iproc,:,iv) = CRS(iproc,iv)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!       ... correct o3 absorption cross-sections
!-----------------------------------------------------------------------
      ZTDIFF = PTEMPER - 230.
      DO iv = ivbegin(J_o3_o),ivend(J_o3_o)
        PCRST(J_o3_o,:,iv) = CRS(J_o3_o,iv) &
     &                 + 1.e-20 * ZTDIFF * (tb_o3(iv) + tc_o3(iv)*ZTDIFF)
        PCRST(J_o3_O1D,:,iv) = PCRST(J_o3_o,:,iv)
      ENDDO

!-----------------------------------------------------------------------
!       ... Calculate CO2 absorption cross-sections from values at 300K
!  (in var crs), 200K and 370K (Dumb linear interpolation) from values in
!  Lewis & Carver, J Quant Spectros Radiat Transfer, vol30, p297, 1983
!-----------------------------------------------------------------------
      DO iv = ivbegin(J_co2),ivend(J_co2)
        WHERE( PTEMPER(:) <= 200. )
           PCRST(J_co2,:,iv) = crs200co2(iv)
        END WHERE
        WHERE( PTEMPER(:) >= 370. )
           PCRST(J_co2,:,iv) = crs370co2(iv)
        END WHERE
        WHERE( (PTEMPER(:) < 300.) .and. (PTEMPER(:) > 200.) )
           PCRST(J_co2,:,iv) = crs200co2(iv) &
     &          + 1.e-2 * (PTEMPER(:)-200.) *(CRS(J_co2,iv)-crs200co2(iv))
        END WHERE
        WHERE( (PTEMPER(:) < 370.) .and. (PTEMPER(:) >= 300.) )
           PCRST(J_co2,:,iv) = crs370co2(iv) &
     &            - (PTEMPER(:)-370.) *(CRS(J_co2,iv)-crs370co2(iv)) / 70.
        END WHERE
      ENDDO

!-----------------------------------------------------------------------
!       ... Calculate n2o absorption cross-sections including T dependence
!       Ref : JPL94,table 16, p.125 : Selwyn & al., 1977,GRL,vol.4,p.427
!-----------------------------------------------------------------------
      PCRST(J_n2o,:,:) = 0.
      DO iv = ivbegin(J_n2o),ivend(J_n2o)
        ZSIGMA = ZA_N2O(4)
        DO IN = 3,0,-1
           ZSIGMA = ZA_N2O(IN) + wlmid(iv)*ZSIGMA
        ENDDO
        ZBSUM = ZB_N2O(3)
        DO IN = 2,0,-1
           ZBSUM = ZB_N2O(IN) + wlmid(iv)*ZBSUM
        ENDDO
        ZBSUM = EXP( ZBSUM )
        PCRST(J_n2o,:,iv) = EXP( ZSIGMA + (PTEMPER - 300.) * ZBSUM )
      ENDDO

!-----------------------------------------------------------------------
!       ... Correct following absorption cross-sections for T dependence
!           See crs97.dat for references
!-----------------------------------------------------------------------
      DO iv = ivbegin(J_no2),ivend(J_no2)
         PCRST(J_no2,:,iv) = CRS(J_no2,iv) + ta_no2(iv)*(PTEMPER - 273.15)
      ENDDO

      DO iv = ivbegin(J_hno3),ivend(J_hno3)
        PCRST(J_hno3,:,iv) = CRS(J_hno3,iv) * EXP( 1.e-3 * tb_hno3(iv) * (PTEMPER - 298.) )
      ENDDO

      DO iv = ivbegin(J_cfc12),ivend(J_cfc12)
        PCRST(J_cfc12,:,iv) = CRS(J_cfc12,iv) &
     &                      * EXP( 4.1e-4*(wlmid(iv) - 184.9)*(PTEMPER - 298.) )
      ENDDO

      ZWORK = MAX( 220., PTEMPER )
      ZWORK = MIN( ZWORK, 298._JPRB )
      ZTDIFF = ZWORK - 296._JPRB
      DO iv = ivbegin(J_clono2_clo),ivend(J_clono2_clo)
        PCRST(J_clono2_clo,:,iv) = CRS(J_clono2_clo,iv) &
     &                       * (1. + ZTDIFF*(ta1_clono2(iv) &
     &                             + ta2_clono2(iv)*ZTDIFF))
         PCRST(J_clono2_cl,:,iv) = PCRST(J_clono2_clo,:,iv)
      ENDDO

!  JPL-06, New formula depending on GAMMA and no more on TA and TB
!          see JPL-06 (eval 15), pp 4-43
      ZWORK = MAX( 223.15, PTEMPER )
      ZWORK = MIN( ZWORK, 293.15_JPRB )
      ZWORK = ZWORK - 298._JPRB
      DO iv = ivbegin(J_ch2o_hco),ivend(J_ch2o_hco)
         PCRST(J_ch2o_hco,:,iv) = CRS(J_ch2o_hco,iv) + gamma_ch2o(iv)*ZWORK*1e-24
         PCRST(J_ch2o_co,:,iv) = PCRST(J_ch2o_hco,:,iv)
      ENDDO

      ZWORK = MAX( 250., PTEMPER )
      ZWORK = MIN( ZWORK, 298._JPRB )
!      DO iv = ivbegin(J_pan),ivend(J_pan)
!         PCRST(J_pan,:,iv) = CRS(J_pan,iv) &
!     &                     * EXP( tb_pan(iv)*(ZWORK - 298.) )
!      ENDDO

!-----------------------------------------------------------------------
!   Calculate n2o5 absorption cross-sections including T dependence in 281-380 nm
!   Source: JPL-2002; WARNING: this param is obsolete, was updated in JPL-2006(=JPL-2011)
!-----------------------------------------------------------------------
      DO iv = ivbegin(J_n2o5),ivend(J_n2o5)
         PCRST(J_n2o5,:,iv) = CRS(J_n2o5,iv)
         IF( wlmid(iv) >= 281. .and. wlmid(iv) <= 380._JPRB )  THEN
           PCRST(J_n2o5,:,iv) = &
     &        1.e-20 * EXP( 2.735_JPRB + ((4728.5_JPRB-17.127_JPRB*wlmid(iv))/PTEMPER) )
         ENDIF
      ENDDO

!-----------------------------------------------------------------------
!       ... Calculate h2o2 absorption cross-sections
!           including T dependence between 259.7 and 352.5 nm
!           See JPL-2011, Table 4B-6, p.259
!-----------------------------------------------------------------------
      ZWORK = 1. / (1. + EXP( -1265./PTEMPER ))
      DO iv = ivbegin(J_h2o2_oh),ivend(J_h2o2_oh)
        PCRST(J_h2o2_oh,:,iv) = CRS(J_h2o2_oh,iv)
        IF( wlmid(iv) >= 259.7_JPRB .and. wlmid(iv) <= 352.5_JPRB ) THEN
            ZASUM = Z_AH2O2(7)
            DO IN = 6,0,-1
               ZASUM = Z_AH2O2(IN) + wlmid(iv)*ZASUM
            ENDDO
            ZBSUM = Z_BH2O2(4)
            DO IN = 3,0,-1
               ZBSUM = Z_BH2O2(IN) + wlmid(iv)*ZBSUM
            ENDDO
            PCRST(J_h2o2_oh,:,iv) = (ZWORK*ZASUM + (1.-ZWORK)*ZBSUM) * 1.e-21_JPRB
         ENDIF
         PCRST(J_h2o2_ho2,:,iv) = PCRST(J_h2o2_oh,:,iv)
      ENDDO

!-----------------------------------------------------------------------
!       ... Correct temperature dependence for the 10 following halogen species :
!                 J_ccl4, J_ch3ccl3, J_ch3cl, J_cfc113,
!                 J_hcfc22, J_ha1211, J_ha1301, J_ch3br, J_cfc114,
!                 J_CFC11 and J_ch2br2
! PCRST(j_CFC11) calculated in CRS_CFC : Feb 2009, sebv@aeronomie.be, JPL-06
! PCRST(j_CH2Br2) calculated in CRS_CFC, Sep 2010, quentin@aeronomie.be
!-----------------------------------------------------------------------
      CALL CRS_CFC( KMXCLY, PCRST, PTEMPER, crs )

!------------------------------------------------------------------
!       ...Special linear-logarithmic interpolation for CCl4 (named CFC10 before 20091019)
!          at wl>250nm (iv>lamax(1)=79) using values found at
!          iv=76 & iv=79                - S.Chabrillat
!------------------------------------------------------------------
      ZLOG_Y0(:) = LOG( PCRST(J_ccl4,:,76) )
      ZLOG_Y1(:) = LOG( PCRST(J_ccl4,:,79) )

      DO iv = 80,ivend(J_ccl4)
        PCRST(J_ccl4,:,iv) = EXP( (ZLOG_Y1 - ZLOG_Y0)*(wlmid(iv) - wlmid(76)) &
     &                           / (wlmid(79) - wlmid(76)) + ZLOG_Y0 )
      ENDDO
!-----------------------------------------------------------------------
!       ... CHBr3 crs are T-dep, range important for atm is wl > 290 nm,
!           JPL2002 gives formula below (p. 4-70) :
!-----------------------------------------------------------------------
      ZWORK = MAX( 210., PTEMPER )
      ZWORK = MIN( ZWORK, 300. )
      DO iv = ivbegin(J_chbr3),ivend(J_chbr3)
         PCRST(J_chbr3,:,iv) = CRS(J_chbr3,iv)
         IF( wlmid(iv) >= 290. .and. wlmid(iv) <= 362. )  THEN
            PCRST(J_chbr3,:,iv) = &
     &       EXP(  (0.06183_JPRB-0.000241_JPRB*wlmid(iv))*(273._JPRB-ZWORK) &
     &                -(2.376_JPRB+0.14757_JPRB*wlmid(iv))   )
         ENDIF
      ENDDO

!-----------------------------------------------------------------------
!       ... Cl2 crs are calculate with 'horrible formula'.
!           JPL2006, pp 4-93
!-----------------------------------------------------------------------
      ZWORK = TANH(402.7_JPRB/PTEMPER)
      DO iv = ivbegin(J_cl2),ivend(J_cl2)
         PCRST(J_cl2,:,iv) = 1e-20 * sqrt(ZWORK) *  &
     &    ( 27.3_JPRB* exp(-99.0_JPRB*ZWORK*(log(329.5_JPRB/wlmid(iv)))**2) + &
     &      0.932_JPRB* exp(-91.5_JPRB*ZWORK*(log(406.5_JPRB/wlmid(iv)))**2) )
      ENDDO

!-----------------------------------------------------------------------
!       ... ClO crs are T-Dep
!           JPL2006, pp 4-93
!-----------------------------------------------------------------------
      ZWORK = MAX( 240., PTEMPER )
      ZWORK = MIN( ZWORK, 298._JPRB )
      DO iv = ivbegin(J_clo),ivend(J_clo)
!         PCRST(J_clo,:,iv) = PCRST(J_clo,:,iv) / (1+0.0036*(ZWORK-298))
!     ! JPL2011 (and 2006) Tref=294 and not 298!
         PCRST(J_clo,:,iv) = PCRST(J_clo,:,iv) / (1._JPRB+0.0036_JPRB*(ZWORK-294))
      ENDDO


!-----------------------------------------------------------------------
!       ... Calculate OCS absorption cross-sections from values 
!  at 225K (in var crs) and 295K
!-----------------------------------------------------------------------
      ZWORK = MAX( 225., PTEMPER )
      ZWORK = MIN( ZWORK, 295. )
      DO iv = ivbegin(J_ocs),ivend(J_ocs)
           PCRST(J_ocs,:,iv) = crs295ocs(iv) &
     &          - (ZWORK(:)-295.) *(CRS(J_ocs,iv)-crs295ocs(iv)) / 70.
      ENDDO


IF (LHOOK) CALL DR_HOOK('CRS_TDEP',1,ZHOOK_HANDLE )

END SUBROUTINE CRS_TDEP


!=======================================================================

SUBROUTINE J_NO_PARAMJGR93( KMXCLY, PCOLX, PUAVG, PTOTDENS, PJJNO )
!-----------------------------------------------------------------------
!       ... Compute NO photolysis following parameterization by
!           Minschwaner and Sisking, JGR, 98, p20401, 1993
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE, only : MXWVN, nabspec

IMPLICIT NONE

!-----------------------------------------------------------------------
!       ... Dummy args - colx are the SLANT overhead columns for NO and O2
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(in)                            :: KMXCLY
      REAL(KIND=JPRB), dimension(0:KMXCLY,nabspec), INTENT(in) :: PCOLX
      REAL(KIND=JPRB), dimension(MXWVN,0:KMXCLY), INTENT(in)   :: PUAVG
      REAL(KIND=JPRB), dimension(0:KMXCLY), INTENT(in)         :: PTOTDENS
      REAL(KIND=JPRB), dimension(0:KMXCLY), INTENT(out)        :: PJJNO

!-----------------------------------------------------------------------
!       ... Parameters
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), parameter :: INS_O2 = 6           ! number of sub-intervals for O2
      INTEGER(KIND=JPIM), parameter :: INSS_NO = 2          ! number of sub-sub-intervals for NO
      REAL(KIND=JPRB), parameter, dimension(INS_O2) ::           &
     &      ZCS2_5_0 = (/  1.117E-23, 2.447E-23, 7.188E-23,      &
     &                     3.042E-22, 1.748E-21, 1.112E-20 /),   &
     &      ZCS2_9_0 = (/  1.350E-22, 2.991E-22, 7.334E-22,      &
     &                     3.074E-21, 1.689E-20, 1.658E-19 /),   &
     &      ZCS2_10_0 = (/ 2.968E-22, 5.831E-22, 2.053E-21,      &
     &                     8.192E-21, 4.802E-20, 2.655E-19 /)

      !-----------------------------------------------------------------------
      ! Uses:
      !     sub-intervals for O2 5-0 at 265K,
      !         sub-sub-intervals for NO 0-0 at 250K
      !     sub-intervals for O2 9-0 band,
      !         sub-sub-intervals for NO 1-0 at 250 K
      !     sub-intervals for O2 10-0 band,
      !         sub-sub-intervals for NO 1-0 at 250 K
      !       ... Note: first 2 sub-intervals in 10-0 band
      !           slightly modified for accuracy
      !-----------------------------------------------------------------------
      REAL(KIND=JPRB), parameter, dimension(INS_O2,INSS_NO)  ::                               &
     &      ZWTNO_5_0 = RESHAPE(                                                             &
     &          SOURCE=(/   0.00E+00, 5.12E-02, 1.36E-01, 1.65E-01, 1.41E-01, 4.50E-02,     &
     &                      0.00E+00, 5.68E-03, 1.83E-02, 1.52E-02, 1.57E-02, 5.00E-03 /),  &
     &          SHAPE=(/INS_O2,INSS_NO/)),                                                    &
     &      ZCSNO_5_0 = RESHAPE(                                                             &
     &          SOURCE=(/   0.00E+00, 1.32E-18, 6.35E-19, 7.09E-19, 2.18E-19, 4.67E-19,     &
     &                      0.00E+00, 4.41E-17, 4.45E-17, 4.50E-17, 2.94E-17, 4.35E-17 /),  &
     &          SHAPE=(/INS_O2,INSS_NO/)),                                                    &
     &      ZWTNO_9_0 = RESHAPE(                                                             &
     &          SOURCE=(/   0.00E+00, 0.00E+00, 1.93E-03, 9.73E-02, 9.75E-02, 3.48E-02,     &
     &                      0.00E+00, 0.00E+00, 2.14E-04, 1.08E-02, 1.08E-02, 3.86E-03 /),  &
     &          SHAPE=(/INS_O2,INSS_NO/)),                                                    &
     &      ZCSNO_9_0 = RESHAPE(                                                             &
     &          SOURCE=(/   0.00E+00, 0.00E+00, 3.05E-21, 5.76E-19, 2.29E-18, 2.21E-18,     &
     &                      0.00E+00, 0.00E+00, 3.20E-21, 5.71E-17, 9.09E-17, 6.00E-17 /),  &
     &          SHAPE=(/INS_O2,INSS_NO/)),                                                    &
     &      ZWTNO_10_0 = RESHAPE(                                                            &
     &          SOURCE=(/   4.50E-02, 1.80E-01, 2.25E-01, 2.25E-01, 1.80E-01, 4.50E-02,     &
     &                      5.00E-03, 2.00E-02, 2.50E-02, 2.50E-02, 2.00E-02, 5.00E-03 /),  &
     &          SHAPE=(/INS_O2,INSS_NO/)),                                                    &
     &      ZCSNO_10_0 = RESHAPE(                                                            &
     &          SOURCE=(/   1.80E-18, 1.50E-18, 5.01E-19, 7.20E-20, 6.72E-20, 1.49E-21,     &
     &                      1.40E-16, 1.52E-16, 7.00E-17, 2.83E-17, 2.73E-17, 6.57E-18 /),  &
     &          SHAPE=(/INS_O2,INSS_NO/))

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: ILEV
      REAL(KIND=JPRB)  ::  ZJNO_5_0, ZJNO_9_0, ZJNO_10_0

IF (LHOOK) CALL DR_HOOK('J_NO_PARAMJGR93',0,ZHOOK_HANDLE )

      DO ILEV = 0,KMXCLY
!         jno_5_0    = PJNO( 55, cs2_5_0, wtno_5_0, csno_5_0 )
!         jno_9_0    = PJNO( 51, cs2_9_0, wtno_9_0, csno_9_0 )
!         jno_10_0   = PJNO( 50, cs2_10_0, wtno_10_0, csno_10_0 )
         CALL JNORATE(ZJNO_5_0 , 55_JPIM, ZCS2_5_0 , ZWTNO_5_0 , ZCSNO_5_0 ,PCOLX(ILEV,:), PUAVG(55,ILEV),PTOTDENS(ILEV) )
         CALL JNORATE(ZJNO_9_0 , 51_JPIM, ZCS2_9_0 , ZWTNO_9_0 , ZCSNO_9_0 ,PCOLX(ILEV,:), PUAVG(51,ILEV),PTOTDENS(ILEV) )
         CALL JNORATE(ZJNO_10_0, 50_JPIM, ZCS2_10_0, ZWTNO_10_0, ZCSNO_10_0,PCOLX(ILEV,:), PUAVG(50,ILEV),PTOTDENS(ILEV) )
         PJJNO(ILEV) = ZJNO_5_0 + ZJNO_9_0 + ZJNO_10_0
      ENDDO
IF (LHOOK) CALL DR_HOOK('J_NO_PARAMJGR93',1,ZHOOK_HANDLE )
END SUBROUTINE J_NO_PARAMJGR93

!      CONTAINS

    !========================================================================

    SUBROUTINE JNORATE ( PJNO, KW, PCS2, PWTNO, PCSNO, PCOLX,PUAVG,PTOTDENS )
    !-----------------------------------------------------------------------
    !           ... Uses xsec at center of g subinterval for O2
    !           uses mean values for NO
    !-----------------------------------------------------------------------

    USE PARKIND1  , ONLY : JPIM     ,JPRB
    USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
    USE BASCOE_TUV_MODULE, only : NABSPEC, O2abs, NOabs
    implicit none
    !-----------------------------------------------------------------------
    !   ... Parameters
    !-----------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JPRB), PARAMETER :: ZPIR = 3.1415926535897931E0
      REAL(KIND=JPRB), PARAMETER :: ZSOLID_ANGLE = 4._JPRB * ZPIR
      INTEGER(KIND=JPIM), PARAMETER :: INS_O2 = 6           ! number of sub-intervals for O2
      INTEGER(KIND=JPIM), PARAMETER :: INSS_NO = 2          ! number of sub-sub-intervals for NO

    !-----------------------------------------------------------------------
    !   ... Dummy args
    !-----------------------------------------------------------------------
      REAL(KIND=JPRB), INTENT(OUT) :: PJNO
      INTEGER(KIND=JPIM), INTENT(IN) :: KW
      REAL(KIND=JPRB), INTENT(IN) :: PCS2(INS_O2)
      REAL(KIND=JPRB), dimension(INS_O2,INSS_NO), INTENT(IN) :: PCSNO, PWTNO
      REAL(KIND=JPRB), dimension(NABSPEC), INTENT(IN)   :: PCOLX
      REAL(KIND=JPRB), INTENT(IN)   :: Puavg
      REAL(KIND=JPRB), INTENT(IN)   :: Ptotdens

    !-----------------------------------------------------------------------
    !   ... Local variables
    !-----------------------------------------------------------------------
      INTEGER(KIND=JPIM) ::  JJ, JK
      REAL(KIND=JPRB) :: ZTAUNO, ZTRANO, ZTRANS, ZTAU
      REAL(KIND=JPRB) :: ZJNO, ZJNO1

IF (LHOOK) CALL DR_HOOK('JNORATE',0,ZHOOK_HANDLE )

      ZJNO = 0._JPRB
      do JK = 1,INS_O2
         ZTAU = PCOLX(O2abs) * PCS2(JK)
         if( ZTAU < 50. ) then
            ZTRANS = EXP( -ZTAU )
         else
            ZTRANS = 0.
         endif
         ZJNO1 = 0.
         do JJ = 1,INSS_NO
            ZTAUNO = PCOLX(NOabs) * PCSNO(JK,JJ)
            if( ZTAUNO < 50. ) then
               ZTRANO = EXP( -ZTAUNO )
            else
               ZTRANO = 0.
            endif
            ZJNO1 = ZJNO1 + PCSNO(JK,JJ) * PWTNO(JK,JJ) * ZTRANO
         enddo
         ZJNO = ZJNO + ZJNO1*ZTRANS
      enddo

      PJNO = ZJNO * Puavg * ZSOLID_ANGLE

      if( KW == 55_JPIM ) then
         PJNO = PJNO * 1.65e9_JPRB &
     &        / (5.1e7_JPRB + 1.65e9_JPRB + 1.5e-9_JPRB*0.79_JPRB*PTOTDENS)
      endif

IF (LHOOK) CALL DR_HOOK('JNORATE',1,ZHOOK_HANDLE )

    END SUBROUTINE JNORATE
    !-----------------------------------------------------------------------



!=======================================================================

SUBROUTINE JRATES( KMXCLY, PXSECT, PQYTP, &
     &                   PUAVG, PRMLYA, PRO2LYA, PRM, PRO2, PDRAT )
!-----------------------------------------------------------------------
!     Calculates the photodissociation rates from
!     layer 0 (top layer) to layer KMXCLY (surface layer)
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_MODULE, only : jprocmax => ndiss, &
    &                       J_o2_O1D, J_o2_o, J_no
USE BASCOE_TUV_MODULE, only : MXWVN, iv_lya, iv_srb0, iv_srb1, nwvn_srb, &
    &                       sw_lya, ivbegin, ivend

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... PARAMETERs
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), PARAMETER :: ZPIR = 3.1415926535897931E0
      REAL(KIND=JPRB), PARAMETER :: ZFOURPI = 4. * ZPIR

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), INTENT(INOUT)  ::  PXSECT(jprocmax,0:KMXCLY,MXWVN)
      REAL(KIND=JPRB), INTENT(in)     ::  PQYTP(jprocmax,0:KMXCLY,mxwvn)
      REAL(KIND=JPRB), DIMENSION(mxwvn,0:KMXCLY), INTENT(IN) :: PUAVG
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN)       :: PRMLYA, PRO2LYA  ! reduction factors for Lyman-a "line"
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY,nwvn_srb), INTENT(IN) :: PRM, PRO2     ! reduction factors for Schumann-R bands
      REAL(KIND=JPRB), INTENT(OUT)    ::  PDRAT(jprocmax,0:KMXCLY)              ! photorates (1/s)

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: IK, IRATE, IWB, IWE, IWAVE

!-----------------------------------------------------------------------
!       ... Calculate photodissociation rates, eq. (18), ISAK
!           used in the 2d code. Skip JNO (43) calculation since
!           J(no) is calculated in subroutine JNO.
!           First modify the cross section array PXSECT to account for :
!           (1) - Lyman-alpha alteration  (iv = 8)
!           (2) - Schumann Runge Bands (iv_srb0=46 <= iv <= iv_srb1=61)
!           (3) - General wavelength quantum yields
!           (4) - o3 z dependent quantum yields
!---------------------------------------------------------------------
!       ... Lyman-alpha line alterations
!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('JRATES',0,ZHOOK_HANDLE )

      IF( sw_lya ) THEN
         DO IRATE = 1,jprocmax
            IF( IRATE == J_no ) THEN
               CYCLE
            ENDIF
            IF( ivbegin(IRATE) <= iv_lya .and. &
     &                               ivend(IRATE) >= iv_lya ) THEN
               IF( IRATE == J_o2_O1D ) THEN
                     PXSECT(J_o2_O1D,:,iv_lya) = PRO2LYA(:)
               ELSE
                     PXSECT(IRATE,:,iv_lya) = PXSECT(IRATE,:,iv_lya) &
     &                                                 * PRMLYA(:)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!---------------------------------------------------------------------
!       ... Schumann Runge band alterations
!---------------------------------------------------------------------
      DO IRATE = 1,jprocmax
         IF( IRATE == J_no ) THEN
            CYCLE
         ENDIF
         IWB = MAX( ivbegin(IRATE), iv_srb0 )
         IWE = MIN( ivend(IRATE), iv_srb1 )
         IF( IWB <= iv_srb1 .and. IWE >= iv_srb0 ) THEN
            IF( IRATE == J_o2_o ) THEN
               DO IWAVE = IWB,IWE
                  PXSECT(J_o2_o,:,IWAVE) = PRO2(:,IWAVE-iv_srb0+1)
               ENDDO
            ELSE
               DO IWAVE = IWB,IWE
                  PXSECT(IRATE,:,IWAVE) = PXSECT(IRATE,:,IWAVE) &
     &                               * PRM(:,IWAVE-iv_srb0+1)
               ENDDO
            ENDIF
         ENDIF

!---------------------------------------------------------------------
!       ... quantum yields
!---------------------------------------------------------------------
         IWB = ivbegin(IRATE)
         IWE = ivend(IRATE)
         DO IWAVE = IWB,IWE
            PXSECT(IRATE,:,IWAVE) = PXSECT(IRATE,:,IWAVE) * PQYTP(IRATE,:,IWAVE)
         ENDDO

!---------------------------------------------------------------------
!       ... Now compute the photorates
!---------------------------------------------------------------------
         DO IK = 0,KMXCLY
            PDRAT(IRATE,IK) =  DOT_PRODUCT( PXSECT(IRATE,IK,IWB:IWE), &
     &                                           PUAVG(IWB:IWE,IK) )
         ENDDO
      ENDDO
      PDRAT = ZFOURPI * PDRAT

IF (LHOOK) CALL DR_HOOK('JRATES',1,ZHOOK_HANDLE )

END SUBROUTINE JRATES

!=======================================================================

SUBROUTINE CALC_QYTP( KMXCLY, PTEMPER, PPMB, PQYTP )
!-----------------------------------------------------------------------
!       ... Calc. T-dependent (optionally p-dependent) quantum yields 
!    NOTE: until BASCOP 07.01.00 this was split into CALC_QY_O3 and CALC_QY_NO3
!----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_MODULE, only : J_O3_O, J_O3_O1D, J_no3_o, J_no3_o2, J_h2so4
USE BASCOE_TUV_MODULE, only : mxwvn, jprocmax, ivbegin, ivend, wlmid, wlbnd, qy, &
    &                   coeff1, coeff2, coeff3, qy_jno3_o, qy_jno3_o2

      IMPLICIT NONE
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN)    :: PTEMPER, PPMB
      REAL(KIND=JPRB), INTENT(OUT)     ::  PQYTP(jprocmax,0:KMXCLY,mxwvn)

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      INTEGER(KIND=JPIM) :: ILC, iv, iproc

      ! for Q.Y. in J(O3):
      REAL(KIND=JPRB), PARAMETER :: ZR = 0.695  ! cm-1 K-1, see JPL2002 p.4-9
      REAL(KIND=JPRB), PARAMETER :: ZNU1 = 0., ZNU2 = 825.518, ZCC = 0.0765
      REAL(KIND=JPRB), DIMENSION(0:KMXCLY) :: ZT_CLIPPED, ZQ1, ZQ2, ZWORK
      
      ! for Q.Y. in J(NO3):
      REAL(KIND=JPRB) :: ZQ1NO3, ZANO3, ZBNO3

      ! for Q.Y. in J(H2SO4):
      REAL(KIND=JPRB), parameter :: ZKRRKM4=1.0/1.7E-7, ZKRRKM5=1.0/8.9E-9, ZKRRKM6=1.0/1.1E-9
      REAL(KIND=JPRB) :: ZKCOLL_LIN
      
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('CALC_QYTP',0,ZHOOK_HANDLE )

!-----------------------------------------------------------------------
!       ... Only the quantum yields recomputed below depend on T (and/or p)
!-----------------------------------------------------------------------
      DO  iv = 1,mxwvn
        DO iproc = 1,jprocmax
          PQYTP(iproc,:,iv) = qy(iproc,iv)
        ENDDO
      ENDDO

!--------------------------------------------------------------------
!       ... Calc the only T-dependent quantum yields : qy_o3 & qy_o3_1d
!           New param from JPL2002 - http://jpldataeval.jpl.nasa.gov
!           implemented at writej/tags/01.00.00_crs2009 i.e. BASCOP q02
!           Temperature range is 200 - 320K
! Previous versions :
! - JPL2000 update : implemented into writej_crs2002
!                    ( -> J-tables for e.g. BASCOE v06q05, BASCOP q01.01)
!                    was also implemented into SOCRATES v06s02 - v7s18
! - Shetter et al. (JGR, 1996, pp.14631-14641) : SOCRATES v1.1 (frozen March 1998)
! - undocumented : SOCRATES "Base" (archived September 1997)
!--------------------------------------------------------------------
      ZT_CLIPPED(:) = MAX( 200., PTEMPER(:) )
      ZT_CLIPPED(:) = MIN( 320., ZT_CLIPPED(:) )
      ZQ1(:) = 1. ! EXP( - ZNU1 / ( R*ZT_CLIPPED(:) ) ) ! but ZNU1 = 0 anyway
      ZQ2(:) = EXP( - ZNU2 / ( ZR*ZT_CLIPPED(:) ) )
      ZWORK(:) = ZT_CLIPPED(:) / 300.

      DO ILC = 0, KMXCLY
         WHERE( wlmid(:) < 306. )
            PQYTP(J_O3_O1D,ILC,:) = 0.9
          ELSE WHERE( wlmid(:) <= 328. )
            PQYTP(J_O3_O1D,ILC,:) = coeff1(:) * ZQ1(ILC) / (ZQ1(ILC)+ZQ2(ILC)) &
     &       + coeff2(:) * ZQ2(ILC) / (ZQ1(ILC)+ZQ2(ILC)) * ZWORK(ILC)*ZWORK(ILC) &
     &       + coeff3(:) * ZWORK(ILC)**1.5 + ZCC
          ELSE WHERE( wlmid(:) < 340. )
            PQYTP(J_O3_O1D,ILC,:) = 0.08
          ELSE where
            PQYTP(J_O3_O1D,ILC,:) = 0.
         END WHERE
      ENDDO

      PQYTP(J_O3_O,:,:) = 1. - PQYTP(J_O3_O1D,:,:)    ! O3 + hv -> O2 + O(3P)
      IF( ANY( PQYTP(J_O3_O,:,:) < 0. ) .or. ANY( PQYTP(J_O3_O,:,:) > 1. ) ) THEN
         CALL ABOR1('CALC_QYTP (file BASCOE_TUV): Fatal error calculating qy_o3')
      ENDIF

!-----------------------------------------------------------------------
!       ... NO3 QY are T dependent
! JPL-06: Quantum Yield of NO3 depend on T in JPL-06. Three QY are given
!         by wavelength from three T: 298, 230 and 190K. Here, linear
!         interpolation between T range is done. If T outside range, use
!         boundary value.
!  Feb 09, sebv@aeronomie.be
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!       ... NO3 + hv -> NO2 + O
!---------------------------------------------------------------------
      DO ILC = 0, KMXCLY
         DO iv = ivbegin(J_no3_o),ivend(J_no3_o)
            IF (PTEMPER(ILC) <= 190.0) THEN
              ZQ1NO3 = qy_jno3_o(1,iv)
            ELSEIF (PTEMPER(ILC)>190.0 .and. PTEMPER(ILC)<=230.0 ) THEN
              ZANO3 = ( qy_jno3_o(2,iv)-qy_jno3_o(1,iv) )/40.0 ! 40 = 230-190K
              ZBNO3 = qy_jno3_o(1,iv)-ZANO3*190.0
              ZQ1NO3 = ZANO3 * PTEMPER(ILC) + ZBNO3
            ELSEIF (PTEMPER(ILC)>230.0 .and. PTEMPER(ILC)<=298.0 ) THEN
              ZANO3 = ( qy_jno3_o(3,iv)-qy_jno3_o(2,iv) )/68.0 ! 68 = 298-230K
              ZBNO3 = qy_jno3_o(2,iv)-ZANO3*230.0
              ZQ1NO3 = ZANO3 * PTEMPER(ILC) + ZBNO3
            ELSE
              ZQ1NO3 = qy_jno3_o(3,iv)
            ENDIF
            PQYTP(j_no3_o,ILC,iv) = ZQ1NO3
         ENDDO
      ENDDO

!---------------------------------------------------------------------
!       ... NO3 + hv -> NO + O2
!---------------------------------------------------------------------
      DO ILC = 0, KMXCLY
         DO iv = ivbegin(J_no3_o2),ivend(J_no3_o2)
            IF (PTEMPER(ILC) <= 190.0) THEN
              ZQ1NO3 = qy_jno3_o2(1,iv)
            ELSEIF (PTEMPER(ILC)>190.0 .and. PTEMPER(ILC)<=230.0 ) THEN
              ZANO3 = ( qy_jno3_o2(2,iv)-qy_jno3_o2(1,iv) )/40.0 ! 40 = 230-190K
              ZBNO3 = qy_jno3_o2(1,iv)-ZANO3*190.0
              ZQ1NO3 = ZANO3 * PTEMPER(ILC) + ZBNO3
            ELSEIF (PTEMPER(ILC)>230.0 .and. PTEMPER(ILC)<=298.0 ) THEN
              ZANO3 = ( qy_jno3_o2(3,iv)-qy_jno3_o2(2,iv) )/68.0 ! 68 = 298-230K
              ZBNO3 = qy_jno3_o2(2,iv)-ZANO3*230.0
              ZQ1NO3 = ZANO3 * PTEMPER(ILC) + ZBNO3
            ELSE
              ZQ1NO3 = qy_jno3_o2(3,iv)
            ENDIF
            PQYTP(j_no3_o2,ILC,iv) = ZQ1NO3
         ENDDO
      ENDDO

!       IF( ANY(PQYTP(j_no3_o,:,:) < 0.d0 ) .or. ANY( PQYTP(j_no3_o:,:) > 1.d0 ) ) THEN
!          stop 'CALC_QY_NO3 (module TUV_MIDATM): Fatal error calculating qy_no3_o'
!       ENDIF
!       IF( ANY(PQYTP(j_no3_o2,:,:) < 0.d0 ) .or. ANY(PQYTP(j_no3_o2:,:) > 1.d0 ) ) THEN
!          stop 'CALC_QY_NO3 (module TUV_MIDATM): Fatal error calculating qy_no3_o2'
!       ENDIF

!-----------------------------------------------------------------------
!       ... H2SO4 qy: apply pressure-dependence based on re-eval of J(H2SO4) 
!   by Miller et al. (GRL, doi:10.1029/2007GL030529, 2007) taking into account 
!   dependence on collision rate which depends on p (ans also T but we neglect it).  
!
!   Implemented in BASCOE's TUV_MIDATM module (chem scheme sb15bs) 
!      by simonc@oma.be on 2018/11/06    (see crs_h2so4.xls for details)
!                                
!-----------------------------------------------------------------------      
      DO ILC = 0, KMXCLY
        ZKCOLL_LIN = 5.56e6 * PPMB(ILC)   ! a linear approx of Table 1 in Miller et al., probably same as neglecting T-dep
        DO iv = ivbegin(J_h2so4),ivend(J_h2so4)
          if(     wlbnd(iv) <= 734.999 .and. wlbnd(iv+1) > 734.999 ) THEN   !   band v=4  
	    ! NOTE: actually for wl=741.3nm but that is LONGWAVE of last bin!
            PQYTP(J_h2so4,ILC,iv) = ZKRRKM4 / ( ZKRRKM4 + ZKCOLL_LIN )
          ELSEIF( wlbnd(iv) <= 606.28 .and. wlbnd(iv+1) > 606.28  ) THEN   !   band v=5
            PQYTP(J_h2so4,ILC,iv) = ZKRRKM5 / ( ZKRRKM5 + ZKCOLL_LIN )
          ELSEIF( wlbnd(iv) <= 518.9  .and. wlbnd(iv+1) > 518.9   ) THEN   !   band v=6
            PQYTP(J_h2so4,ILC,iv) = ZKRRKM6 / ( ZKRRKM6 + ZKCOLL_LIN )
          ENDIF
        ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('CALC_QYTP',1,ZHOOK_HANDLE )

END SUBROUTINE CALC_QYTP


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! from MODULE NUMERICAL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE TRIDLA( KN, PSUB, PDIAG, PSUP, PX )
!-------------------------------------------------------------------
!
! purpose                trdi computes the solution of the tridiagonal
!                        linear system,
!                            b(1)*x(1)+c(1)*x(2)               = y(1)
!                            a(i)*x(i-1)+b(i)*x(i)+c(i)*x(i+1) = y(i)
!                                i=2,3,...,n-1
!                            a(n)*x(n-1)+b(n)*x(n)             = y(n)
!
! usage                  call trdi (n, a, b, c, x )
!
! arguments
!
! on input               n
!                          the number of unknowns.  this subroutine
!                          requires that  n  be greater than  2.
!
!                        sub
!                          the subdiagonal of the matrix is stored in
!                          locations a(2) through a(n).
!
!                        diag
!                          the diagonal of the matrix is stored in
!                          locations b(1) through b(n).
!
!                        sup
!                          the super-diagonal of the matrix is stored in
!                          locations c(1) through c(n-1).
!
!                        x
!                          the right-hand side of the equations is
!                          stored in y(1) through y(n).
!
! on output              x
!                          an array which contains the solution to the
!                          system of equations.
!
! special conditions     this subroutine executes satisfactorily
!                        if the input matrix is diagonally dominant
!                        and non-singular.  the diagonal elements
!                        are used to pivot, and no tests are made to
!                        determine singularity.  if a singular, or
!                        nearly singular, matrix is used as input,
!                        a divide by zero or floating point overflow
!                        may result.
!
! precision              (as defined)
!
! language               fortran
!
! history                originally written by nancy werner at ncar
!                        in october, 1973. modified by s. walters at
!                        ncar in june 1989 to functionally replace
!                        the cal routine tridla.
!
! portability            fortran 90
!
! algorithm              an lu-decomposition is obtained using the
!                        algorithm described in the reference below.
!
!                        to avoid unnecessary divisions, the alpha
!                        values used in the routine are the
!                        reciprocals of the alpha values described
!                        in the reference below.
!
! accuracy               every component of the residual of the linear
!                        system (i.e. the difference between  y  and
!                        the matrix applied to x) should be less in
!                        magnitude than ten times the machine precision
!                        times the matrix order times the maximum
!                        absolute component of the solution vector
!                        times the largest absolute row sum of the
!                        input matrix.
!
! timing                 the timing is roughly proportional to the
!                        order n of the linear system.
!
! references             analysis of numerical methods by
!                        e. isaacson and h. keller
!                        (john wiley and sons, 1966) pp. 55-58.
!-------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

      IMPLICIT NONE

!-------------------------------------------------------------------
!       ... Dummy args
!-------------------------------------------------------------------
      INTEGER(KIND=JPIM), INTENT(IN) ::  KN
      REAL(KIND=JPRB), INTENT(IN)  ::  PSUB(KN)
      REAL(KIND=JPRB), DIMENSION(KN), INTENT(INOUT) ::  PDIAG, PSUP, PX

!-------------------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------------------
      INTEGER(KIND=JPIM) :: I, INM1
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TRIDLA',0,ZHOOK_HANDLE )

      INM1 = KN - 1
!----------------------------------------------------------------------
!       ... Perform the lu-decomposition
!----------------------------------------------------------------------
      PDIAG(1) = 1. / PDIAG(1)
      PSUP(1) = PSUP(1)*PDIAG(1)
      DO I = 2,INM1
         PDIAG(I) = 1. / (PDIAG(I) - PSUB(I)*PSUP(I-1))
         PSUP(I) = PSUP(I)*PDIAG(I)
      ENDDO
!----------------------------------------------------------------------
!       ... Solve the system
!----------------------------------------------------------------------
      PX(1) = PX(1)*PDIAG(1)
      DO I = 2,INM1
         PX(I) = (PX(I) - PSUB(I)*PX(I-1))*PDIAG(I)
      ENDDO

      PX(KN) = (PX(KN) - PSUB(KN)*PX(INM1)) / (PDIAG(KN) - PSUB(KN)*PSUP(INM1))
      DO I = INM1,1,-1
         PX(I) = PX(I) - PSUP(I)*PX(I+1)
      ENDDO

IF (LHOOK) CALL DR_HOOK('TRIDLA',1,ZHOOK_HANDLE )

END SUBROUTINE TRIDLA
