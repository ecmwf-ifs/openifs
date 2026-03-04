! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_TUV_INI

!**   DESCRIPTION
!     ----------
!
!   Initialize parameters for BASCOE_TUV
!       -read tables of cross-sections, quantum yields, and possible
!           additional variables from data files
!       -precompute wavelength-dependent variables
!
!   input cross-section (ASCII) files are  named
!           'crs_[j_shortname].dat'
!       where j_shortname are HARDCODED in module TUV_NAMES
!
!   auxiliary ASCII files are needed:
!         'solflux_socrates.dat', 'solflux_socrates-lean97.dat' and 'srb_kockarts94.dat'
!   NOTE these filenames are HARDCODED in the subroutines which use them !
!
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_TUV_INI* IS CALLED FROM *J_INI*.
!
!
!     MODIFICATIONS.
!     -------
!       2018-03-22:    YVES CHRISTOPHE    *BIRA* !YC
!           photolysis rates computed by TUV_midatm (jon)
!       2017-12-07:    YVES CHRISTOPHE    *BIRA* !YC
!           photolysis tables computed by TUV_midatm (joff)
!---------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE,    ONLY : MXWVN, WLMID, WLBND, &
                            &    init_done, hv, denerg_O2, denerg_O3, &
                            &    coeff1, coeff2, coeff3

IMPLICIT NONE

!-----------------------------------------------------------------------

    REAL(KIND=JPRB), PARAMETER :: ZCLIGHT = 2.99792458d8       ! speed of light  (m/s)
    REAL(KIND=JPRB), PARAMETER :: ZHPL    = 6.626068d-34       ! Planck constant (J*s)

    !    ... Parameters for Heating rates
    REAL(KIND=JPRB), PARAMETER :: ZEN_O2 = 8.18674e-19   !5.11/6.2418e18  O2
    REAL(KIND=JPRB), PARAMETER :: ZEN_O3 = 1.6822e-19    !1.05/6.2418e18  O3

    !   ... Parameters for CALC_QY_NO3
    REAL(KIND=JPRB), PARAMETER :: ZX1 = 304.225,   ZX2 = 314.957,  ZX3=310.737
    REAL(KIND=JPRB), PARAMETER :: ZOM1 = 5.576,    ZOM2 = 6.601,   ZOM3=2.187
    REAL(KIND=JPRB), PARAMETER :: ZA1 = 0.8036,    ZA2 = 8.9061,   ZA3=0.1192

    ! * LOCAL
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI',0,ZHOOK_HANDLE )

      CALL SOLFLUX_READ( )

      CALL CRS_READ( )
      CALL SRB_READ( )
      CALL RAYLEIGH( )

!-----------------------------------------------------------------------
!   Initializations for Heating rates (wavelength dependent only)
!-----------------------------------------------------------------------
      hv(1:MXWVN) = 1.e9 * ZHPL * ZCLIGHT / WLMID(1:MXWVN)    ! 1e9 for wlmid from nm to m
      denerg_O2(1:MXWVN) = MAX( hv(1:MXWVN) - ZEN_O2, 0. )
      denerg_O3(1:MXWVN) = MAX( hv(1:MXWVN) - ZEN_O3, 0. )

!-----------------------------------------------------------------------
!   Initializations of vars for qy_o3 (wavelength dependent only)
!-----------------------------------------------------------------------
      coeff1(1:MXWVN) = ZA1 * EXP( - ((ZX1-WLMID(1:MXWVN))/ZOM1)**4 )
      coeff2(1:MXWVN) = ZA2 * EXP( - ((ZX2-WLMID(1:MXWVN))/ZOM2)**2 )
      coeff3(1:MXWVN) = ZA3 * EXP( - ((ZX3-WLMID(1:MXWVN))/ZOM3)**2 )

      init_done = .true.

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI',1,ZHOOK_HANDLE )

CONTAINS 


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! from MODULE TUV_MIDATM
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE RAYLEIGH( )
!-----------------------------------------------------------------------
!       ... Calculates Rayleigh scattering cross-sections (cm2) for
!           Earth following Nicolet(1984): Planet. Space Sci., 32, p1467
!           Although the formula is valid only above 200nm, it is used
!           below as well. Assumes a constant air molecular weight - ie
!           not exact in heterosphere.
!-----------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE,    ONLY : MXWVN, WLMID, RAYCRS

      implicit none

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      REAL(KIND=JPRB), DIMENSION(MXWVN) :: ZLUM, ZX
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:RAYLEIGH',0,ZHOOK_HANDLE )
      ZLUM = 1.e-3 * WLMID                        ! wavelengths in micrometers
      ZX(1:MXWVN) = 0.04                ! BUGFIX 20180605: was 0.4 until BASCOP 6.5, resulting in raycrs too large by ~10% ! 
      WHERE( ZLUM(1:MXWVN) < 0.55 )
         ZX(1:MXWVN) = 0.389 * ZLUM(1:MXWVN) + 0.09426 / ZLUM(1:MXWVN) - 0.3228
      END WHERE
      RAYCRS(1:MXWVN) = 4.02e-28 / ZLUM(1:MXWVN)**(4.+ZX(1:MXWVN))
IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:RAYLEIGH',1,ZHOOK_HANDLE )

END SUBROUTINE RAYLEIGH

!=======================================================================

SUBROUTINE SOLFLUX_READ( )
!-----------------------------------------------------------------------
! ... Read in wavelength grid and solar flux from files
!
!  reads 'solflux_socrates.dat' and
!      *partly* overwrites it with 'solflux_socrates-lean97.dat'
!      at SOLAR MIN  (see fbeamr = fluxavg/( 1+0.5*solvar )
!
!  if BASCOE_SOLFLUX.dat exists, allocate fbeamr2d and read the
!       time dependent solar flux (check that the wavelength grid is 
!       the same as in solflux_socrates.dat)
!       
!                Note HARDCODED filenames !!
!-----------------------------------------------------------------------
USE YOMLUN    , ONLY : NULOUT
USE PARKIND1  , ONLY : JPIM,    JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE,    ONLY : MXWVN, WLMID, WLBND, FBEAMR,   &
                                & DAILY_SOLFLUX, FBEAM_DATES, FBEAMR2D

      implicit none

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: CL_MY_NAME = 'SOLFLUX_READ'
      INTEGER(KIND=JPIM), PARAMETER  :: IUTMP=77

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: IMXWVN_FILE, IV, IOS, IDUMMY, INL, il
      REAL(KIND=JPRB)  :: ZSOLVAR, ZAL, ZBL, ZFMAX
      REAL(KIND=JPRB), DIMENSION(MXWVN) :: ZLAMB0, ZLAMB1, ZFLUXAVG
      REAL(KIND=JPRB), DIMENSION(MXWVN) :: ZWL
      CHARACTER(LEN=128) :: CL_FILENM, CLLINE
      LOGICAL :: LLFILEOK

#include "abor1.intfb.h"
#include "bascoe_cskip.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:SOLFLUX_READ',0,ZHOOK_HANDLE )

      !-----------------------------------------------------------------------
      !     Open file with default solmin flux + variability...
      !-----------------------------------------------------------------------
      CL_FILENM = 'solflux_socrates.dat'
      OPEN( unit = IUTMP, &
     &      file = CL_FILENM, &
     &      status = 'OLD', form = 'formatted', iostat = ios )
      IF( ios /= 0 ) THEN
        CALL ABOR1(CL_MY_NAME//' Error opening '//TRIM(CL_FILENM))
      ENDIF

      !-----------------------------------------------------------------------
      !     ... Read header, wavelength grid and solflux for TUV...
      !-----------------------------------------------------------------------
      CALL BASCOE_CSKIP( '!', IUTMP)
      READ( IUTMP, * ) IMXWVN_FILE
      CALL BASCOE_CSKIP( '!', IUTMP)
      IF( IMXWVN_FILE /= MXWVN ) THEN
        CALL ABOR1(CL_MY_NAME//' Error Number of Wavelength mismatch in '//TRIM(CL_FILENM))
      ENDIF
      DO IV = 1,MXWVN
        READ(IUTMP,*) IDUMMY, ZLAMB0(IV), ZLAMB1(IV), ZFLUXAVG(IV), ZSOLVAR
        FBEAMR(IV) = ZFLUXAVG(IV) / ( 1. + 0.5*ZSOLVAR )
      ENDDO
      CLOSE( unit = IUTMP )
      WLMID = .5 * ( ZLAMB0 + ZLAMB1 )
      WLBND(1:mxwvn) = ZLAMB0(1:MXWVN)
      WLBND(mxwvn+1) = ZLAMB1(mxwvn)

      !-----------------------------------------------------------------------
      !   ...overwrite with solar flux at solmin & solmax (122-417.5 nm)  by Lean/et-al-1997
      !-----------------------------------------------------------------------
      CL_FILENM = 'solflux_socrates-lean97.dat'
      OPEN( unit = IUTMP, &
     &      file = CL_FILENM, &
     &      status = 'OLD', form = 'formatted', iostat = ios )
      IF( ios /= 0 ) THEN
        CALL ABOR1(CL_MY_NAME//' Error opening '//TRIM(CL_FILENM))
      ENDIF
      CALL BASCOE_CSKIP( '!', IUTMP)
      READ(IUTMP,*) INL
      CALL BASCOE_CSKIP( '!', IUTMP)
      DO il = 1, INL
        READ(IUTMP,*) IV, ZAL, ZBL, FBEAMR(IV), ZFMAX
        IF( ZAL /= ZLAMB0(IV) .or. ZBL /= ZLAMB1(IV) ) THEN
          CALL ABOR1(CL_MY_NAME//' Error Wavelength grid mismatch '//TRIM(CL_FILENM))
        ENDIF
      ENDDO
      CLOSE( unit = IUTMP )

      !-----------------------------------------------------------------------
      !   ...if a daily solar flux file BASCOE_SOLFLUX.dat exists, open it
      !-----------------------------------------------------------------------
      CL_FILENM = 'BASCOE_SOLFLUX.dat'
      INQUIRE(FILE=CL_FILENM, EXIST=LLFILEOK)
      IF(.NOT. LLFILEOK)THEN
        WRITE(NULOUT,*) CL_MY_NAME // ': No daily solflux file found: ' // CL_FILENM
        DAILY_SOLFLUX = .FALSE.
      ELSE
        WRITE(NULOUT,*) CL_MY_NAME // ': Reading daily solflux file: ' // CL_FILENM
        OPEN( unit = IUTMP, &
       &      file = CL_FILENM, &
       &      status = 'OLD', form = 'formatted', iostat = ios )
        IF( ios /= 0 ) THEN
          CALL ABOR1(CL_MY_NAME//' Error opening '//TRIM(CL_FILENM))
        ENDIF
        DAILY_SOLFLUX = .TRUE.
        
        !-----------------------------------------------------------------------
        !    expected format:
        !       - comment lines begin with '#' are ignored
        !       - WL_BOUNDS     2   <number of wavelengths>
        !       - table of wavelength bounds (nm): 
        !               row 1: lower bounds
        !               row 2: upper bounds
        !       - SOLFLUX       <number of dates>  <number of wavelengths>
        !       - table of solar fluxes (photons/cm2/s)
        !               column 1: date (YYYYmmdd)
        !               
        !-----------------------------------------------------------------------
        CALL BASCOE_CSKIP( '#', IUTMP)
        READ( IUTMP, '(a)' , IOSTAT=IOS) clline
        clline = ADJUSTL(clline)
        IF( IOS /= 0 .OR. (clline(1:9)) /= 'WL_BOUNDS' ) THEN
          CALL ABOR1(CL_MY_NAME//' Error reading WL_BOUNDS in '//TRIM(CL_FILENM))
        ENDIF
        READ( clline(10:), *, IOSTAT=IOS) INL, IMXWVN_FILE
        IF( IOS /= 0 .OR. INL /= 2 .OR. IMXWVN_FILE /= MXWVN ) THEN
          CALL ABOR1(CL_MY_NAME//' Error in WL_BOUNDS size in '//TRIM(CL_FILENM))
        ENDIF
        CALL BASCOE_CSKIP( '#', IUTMP)
        READ(IUTMP,* , IOSTAT=IOS) (ZWL(IV),  IV=1,MXWVN)
        IF (  IOS /= 0 .OR. ANY(ZWL /= ZLAMB0) ) THEN
          CALL ABOR1(CL_MY_NAME//' Error Wavelength grid mismatch '//TRIM(CL_FILENM))
        ENDIF
        READ(IUTMP,* , IOSTAT=IOS) (ZWL(IV),  IV=1,MXWVN)
        IF (  IOS /= 0 .OR. ANY(ZWL /= ZLAMB1) ) THEN
          CALL ABOR1(CL_MY_NAME//' Error Wavelength grid mismatch '//TRIM(CL_FILENM))
        ENDIF

        CALL BASCOE_CSKIP( '#', IUTMP)
        READ( IUTMP, '(a)' , IOSTAT=IOS) clline
        clline = ADJUSTL(clline)
        IF( IOS /= 0 .OR. (clline(1:7)) /= 'SOLFLUX' ) THEN
          CALL ABOR1(CL_MY_NAME//' Error reading SOLFLUX in '//TRIM(CL_FILENM))
        ENDIF
        READ( clline(8:), *, IOSTAT=IOS) INL, IMXWVN_FILE
        IF( IOS /= 0 .OR. IMXWVN_FILE /= MXWVN ) THEN
          CALL ABOR1(CL_MY_NAME//' Error in WL_BOUNDS size in '//TRIM(CL_FILENM))
        ENDIF
        WRITE(NULOUT,*)  CL_MY_NAME // ': Reading ', INL, 'daily solflux in file: ' // CL_FILENM
        ALLOCATE( FBEAM_DATES(INL), FBEAMR2D(INL,MXWVN) )
        DO IL = 1, INL
          READ(IUTMP,*, IOSTAT=IOS) FBEAM_DATES(IL), (FBEAMR2D(IL, IV),  IV=1,MXWVN)
          IF( IOS /= 0 ) THEN
            CALL ABOR1(CL_MY_NAME//' Error Wavelength grid mismatch '//TRIM(CL_FILENM))
          ENDIF
        ENDDO

      ENDIF
      
    CLOSE(IUTMP)


IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:SOLFLUX_READ',1,ZHOOK_HANDLE )

END SUBROUTINE SOLFLUX_READ

!=======================================================================

SUBROUTINE SRB_READ( )

!-----------------------------------------------------------------------
!       ... Read from 'srb_kockarts94.dat' the parameters ako, bko
!           that will be used by SRB to calc the reduction factors
!           Rm and Ro2 for Schumann-Runge bands:
!                  Kockarts, Ann Geophys, 12, p1207, 1994
!           The 16 wavelength intervals correspond to iv=46-61 in lamb(MXWVN)
!
!       ... note HARDCODED filename 'srb_kockarts94.dat'
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_TUV_MODULE,    ONLY : nwvn_srb, ako, bko

      implicit none

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: CL_MY_NAME = 'SRB_READ'
      INTEGER(KIND=JPIM), PARAMETER  :: IUTMP=77

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JPIM) :: ios, IV, ic
      CHARACTER(LEN=128) :: CL_FILENM

#include "abor1.intfb.h"
#include "bascoe_cskip.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:SRB_READ',0,ZHOOK_HANDLE )

      CL_FILENM = 'srb_kockarts94.dat'
      OPEN( unit = IUTMP, &
     &      file = CL_FILENM, &
     &      form = 'formatted', &
     &      status = 'OLD', &
     &      iostat = ios )
      IF( ios /= 0 ) THEN
        CALL ABOR1(CL_MY_NAME//' Error opening '//TRIM(CL_FILENM))
      ENDIF
      CALL BASCOE_CSKIP( '!', IUTMP)
      DO IV = nwvn_srb,1,-1
         READ(IUTMP,*)
         READ(IUTMP,*) (ako(ic,IV),ic=1,12)
      ENDDO
      CALL BASCOE_CSKIP( '!', IUTMP)
      DO IV = nwvn_srb,1,-1
         READ(IUTMP,*)
         READ(IUTMP,*) (bko(ic,IV),ic=1,12)
      ENDDO
      CLOSE(IUTMP)

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:SRB_READ',1,ZHOOK_HANDLE )

END SUBROUTINE SRB_READ

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! from MODULE TUV_MIDATM_CRS (sb14b/sb15b specific)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE CRS_READ( )
!-----------------------------------------------------------------------
! For *all* photodissociation process: read cross-sections (crs)
!   and quantum yields from crs_*.dat (this includes reading
!   the *absorption* cross-section for O2 and O3)
!
! Read, for some specified processes, additional parameters from the file:
!   T dependence of O3, CO2, NO2, HNO3, CLONO2 cross sections
!   T dependence of NO3 quantum yields
!     JPL-06: Quantum Yield of NO3 depend on T in JPL-06. Three QY are given
!             by wavelength from three T: 298, 230 and 190K.
!       Feb 09, sebv@aeronomie.be
!   T dependence of CH2O cross sections
!     JPL-06, New formula depending on GAMMA and no more on TA and TB
!          see JPL-06 (eval 15), pp 4-43
!-----------------------------------------------------------------------
USE YOMLUN    , ONLY : NULOUT
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_J_MODULE, only : jprocmax => ndiss, jnames, J_O3_O, J_NO2, J_HNO3, &
     &                          J_CLONO2_CLO, J_CH2O_HCO, J_NO3_O, J_NO3_O2, J_CO2, J_OCS
USE BASCOE_TUV_MODULE, only : crs, qy, ivbegin, ivend,      &
     &                      tb_o3, tc_o3, ta_no2, tb_hno3,  &
     &                      ta1_clono2, ta2_clono2,         &
     &                      crs200co2, crs370co2,           &
     &                      qy_jno3_o, qy_jno3_o2,          &
     &                      gamma_ch2o, crs295ocs !, ta_n2o5, tb_n2o5

      implicit none

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: CL_MY_NAME = 'CRS_READ'
      INTEGER(KIND=JPIM), PARAMETER  :: IUTMP=77

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      INTEGER(KIND=JPIM)    ::  iproc, ios, IV, IDUMMY
      CHARACTER(LEN=128)    :: CL_FILENM
      REAL(KIND=JPRB)       ::  ZDUMMY

#include "abor1.intfb.h"
#include "bascoe_cskip.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:CRS_READ',0,ZHOOK_HANDLE )

!-----------------------------------------------------------------------
!       ... Initializations
!-----------------------------------------------------------------------
      crs(:,:) = 0.
      qy(:,:) = 1.
      
      tb_o3(:) = 0.
      tc_o3(:) = 0.
      crs200co2(:) = 0.
      crs370co2(:) = 0.
      
      ta_no2(:) = 0.
      qy_jno3_o(:,:) = 0.
      qy_jno3_o2(:,:) = 0.
      
      tb_hno3(:) = 0.
      ta1_clono2(:) = 0.
      ta2_clono2(:) = 0.
      
      gamma_ch2o(:) = 0.
      crs295ocs(:)= 0.
      
      write(NULOUT,*) CL_MY_NAME//': starts'
      DO iproc = 1, jprocmax
        !-----------------------------------------------------------------------
        !       ... open cross-section file
        !-----------------------------------------------------------------------
        CL_FILENM = 'crs_' // TRIM(jnames(iproc))// '.dat'
        open( unit = IUTMP, &
       &      file = CL_FILENM, &
       &      status = 'OLD', &
       &      form = 'formatted', &
       &      iostat = ios )
        IF( ios /= 0 ) THEN
          CALL ABOR1(CL_MY_NAME//' Error opening '//TRIM(CL_FILENM))
        ENDIF        
        CALL BASCOE_CSKIP( '!',IUTMP )
        !-----------------------------------------------------------------------
        !       ... read range for index of wavelengths
        !-----------------------------------------------------------------------

        READ(IUTMP,'(31x,i3,31x,i3)',iostat=ios) ivbegin(iproc), ivend(iproc)
        IF (ios /= 0 .or. ivbegin(iproc) > ivend(iproc)) THEN
          CALL ABOR1( CL_MY_NAME//' ERROR reading wavelength range in '//TRIM(CL_FILENM) )
        ENDIF
        CALL BASCOE_CSKIP( '!',IUTMP )
        !-----------------------------------------------------------------------
        !       ... read cross-sections, quantum yields; and depending on
        !              photolysis reaction, additional parameters
        !-----------------------------------------------------------------------
        DO IV = ivbegin(iproc),ivend(iproc)
          IF (iproc == J_o3_o) THEN             ! T dependence of o3 cross sections
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy(iproc,IV), &
            &               tb_o3(IV), tc_o3(IV)
          ELSEIF (iproc == J_co2) THEN          ! CO2 cross sections at 200K and 370K
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy(iproc,IV), &
            &               crs200co2(IV), crs370co2(IV)
          ELSEIF (iproc == J_no2) THEN          ! T dependence of no2 cross sections
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy(iproc,IV), &
            &               ta_no2(IV)
          ELSEIF (iproc == J_no3_o) THEN        ! T dependence of no3 yield: no3-o
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy_jno3_o(1,IV), &
            &               qy_jno3_o(2,IV), qy_jno3_o(3,IV)
          ELSEIF (iproc == J_no3_o2) THEN       ! T dependence of no3 yield: no3-o2
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy_jno3_o2(1,IV), &
            &               qy_jno3_o2(2,IV), qy_jno3_o2(3,IV)
          ELSEIF (iproc == J_hno3) THEN         ! T dependence of hno3 cross sections
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy(iproc,IV), &
            &               tb_hno3(IV)
          ELSEIF (iproc == J_clono2_clo) THEN   ! T dependence of clono2 cross sections
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy(iproc,IV), &
            &               ta1_clono2(IV), ta2_clono2(IV)
          ELSEIF (iproc == J_ch2o_hco) THEN     ! T dependence of ch2o cross sections
            READ(IUTMP,*,iostat=ios) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy(iproc,IV), &
            &               gamma_ch2o(IV)
          ELSEIF (iproc == J_ocs) THEN     ! OCS cross-sections at 225K and 295K
            READ(IUTMP,*) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,iv), qy(iproc,iv), &
            &               crs295ocs(iv)
          ELSE                                  ! general case
            READ(IUTMP,*) IDUMMY,ZDUMMY,ZDUMMY, crs(iproc,IV), qy(iproc,IV)
          ENDIF
          IF (ios /= 0) THEN
            CALL ABOR1( CL_MY_NAME//' ERROR reading data in '//TRIM(CL_FILENM) )
          ENDIF
        ENDDO
        CLOSE( IUTMP )
        
      ENDDO

IF (LHOOK) CALL DR_HOOK('BASCOE_TUV_INI:CRS_READ',1,ZHOOK_HANDLE )

END SUBROUTINE CRS_READ

END SUBROUTINE BASCOE_TUV_INI
