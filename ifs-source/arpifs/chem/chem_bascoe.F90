! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_BASCOE &
 &    (YDML_GCONF, YDCHEM, KSTEP, KIDIA, KFDIA, KLON, KLEV,             &
 &     PTSTEP ,PDELP, PRS1, PRSF1, PGEOH, PQP, PTP,  PALB, PCSZA,       &
 &     PGELAT, PGELAM,  PCEN , PTENC1, POUT)

!**   DESCRIPTION
!     ----------
!
!   routine for C-IFS-BASCOE stratospheric chemistry
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_BASCOE* IS CALLED FROM *CHEM_MAIN*.

! INPUTS:
! -------
! KSTEP                      : Time step
! KIDIA                      : Start of Array
! KFDIA                      : End  of Array
! KLON                       : Length of Arrays
! KLEV                       : Number of Levels
! PTSTEP                     : Time step in seconds
! PDELP   (KLON,KLEV)        : PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PRS1    (KLON,0:KLEV)      : HALF-LEVEL PRESSURE           (Pa)
! PRSF1   (KLON,KLEV)        : FULL-LEVEL PRESSURE           (Pa)
! PGEOH   (KLON,KLEV)        : Geopotential                  (m*m/s*s)
! PQP     (KLON,KLEV)        : SPECIFIC HUMIDITY             (kg/kg)
! PTP     (KLON,KLEV)        : TEMPERATURE                   (K)
! PALB    (KLON)             : Surface albedo
! PCSZA   (KLON)             : COS of Solar Zenit Angle
! PGELAT  (KLON)             : LATITUDE (RADIANS)
! PGELAM  (KLON)             : LONGITUDE (RADIANS)
! PCEN    (KLON,KLEV,NCHEM)  : CONCENTRATION OF TRACERS           (kg/kg)
!
! OUTPUTS:
! -------
! PTENC1  (KLON,KLEV,NCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS BECAUSE OF CHEMISTRY (kg/kg s-1), no update
! POUT    (KLON,KLEV,5)         : additional output, e.g. UBC contribution , Photolysis rates O3 , NO2, tau for output
!
! INOUTPUTS:
! -------
! YDML_GCONF                    : Model general configuration
! YDCHEM                        : Model chem type
!
! LOCAL:
! -------
!
! ZCVM0(KLON,NCHEM)       : initial volume ratios OF TRACERS           (molec/cm3)
! ZCVM (KLON,NCHEM+3)     : final   volume ratios OF TRACERS           (molec/cm3)
!
!
!     AUTHOR.
!     -------
!        JOHANNES FLEMMING  *ECMWF*
!        VINCENT HUIJNEN    *KNMI*
!        Quentin Errera     *BIRA*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2014-02-01
!        BASCOE-implementation of BASCOE scheme            : 2014-03-03
!        BASCOE sb15b J online                             : 2018-03-25
!
!     NOTES
!     -----
!   POUT is used for extra fields to save:
!       consistency of its usage should be carefully checked !
!
!-----------------------------------------------------------------------



USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! NCHEM : number of chemical species
! YCHEM : Data structure with Chemistry meta data
USE YOMLUN   , ONLY : NULERR, NULOUT
USE YOMCHEM  , ONLY : TCHEM
USE YOMCST   , ONLY : RMD, RG , RPI , RNAVO
USE YOMRIP0  , ONLY : NINDAT
USE BASCOE_MODULE, ONLY : IO3, INO, ICO2, INO2, IH2O, ISTRATAER, NHET, NBINS, NAER, &
  &                       IN2O5, IHCL, IHOCL, ICLONO2,IHOBR, IHBR,IBRONO2, IHNO3,   &
  &                       NBC, BASCOE_BC

! General KPP settings
USE CIFS_KPP_INTPARAM  , ONLY : HMIN,HSTART,RTOLS_G,IAUTONOM,IROSMETH, VMR_BAD_LARGE

! BASCOE chemistry...
USE BASCOE_LBC_MODULE   , ONLY : NLATBOUND_LBC, XLATBOUND_LBC,  &
  &                              MONTH_LBC, VALUES_LBC
USE BASCOE_J_MODULE , ONLY : NDISS

! Optionally declare here J rates which may be written in extra fields POUT 
! USE BASCOE_J_MODULE , ONLY : J_O2_O, J_O3_O, J_O3_O1D, J_NO2


USE BASCOE_TUV_MODULE    , ONLY : mxwvn, NABSPEC, fbeamr, fbeamr2d, fbeam_dates, daily_solflux, &
  & AIR, O2ABS, O3ABS, NOABS, CO2ABS, NO2ABS
USE BASCOE_KPP_PARAMETERS, ONLY : NREACT,NVAR,NFIX !,NHEXIT
USE BASCOE_KPP_GLOBAL    , ONLY : RTOL,ATOL,ROUNDOFF_STORE

!-----------------------------------------------------------------------

IMPLICIT NONE

!*       0.1  ARGUMENTS
!             ---------

TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TCHEM)                  ,INTENT(INOUT):: YDCHEM
INTEGER(KIND=JPIM),INTENT(IN) :: KSTEP, KIDIA , KFDIA , KLON , KLEV
REAL(KIND=JPRB)   ,INTENT(IN) :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN) :: PDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PTP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PALB(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PCSZA(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGELAT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PCEN(KLON,KLEV,YDML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB)   ,INTENT(OUT):: PTENC1(KLON,KLEV,YDML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB)   ,INTENT(OUT):: POUT(KLON,KLEV,5)


!*       0.2 Local PARAMETERS
!           ----------------

  CHARACTER(LEN=*), PARAMETER    :: CL_MY_NAME   = 'CHEM_BASCOE'


!*       0.5   LOCAL VARIABLES
!              ---------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * Lat /Lon time
REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLAT
REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLON
INTEGER(KIND=JPIM)                  :: IMONTH0, IDAY0, IYEAR0, ILMONTHS(12)
INTEGER(KIND=JPIM)                  :: IMONTH, IDAY, IYEAR, IYYYYMM, IYYYYMMDD, IJUL, IJULMAX
REAL(KIND=JPRB)                     :: ZJUL1, ZUT

! * counters
INTEGER(KIND=JPIM) :: JK, JL, JT, JLEV, JB
!INTEGER(KIND=JPIM) :: IDRYDEP

! * chemical data 
REAL(KIND=JPRB) , DIMENSION(KLON,YDML_GCONF%YGFL%NCHEM+3)   :: ZCVM
REAL(KIND=JPRB) , DIMENSION(KLON,YDML_GCONF%YGFL%NCHEM)     :: ZCVM0 
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)                      :: ZDENS
REAL(KIND=JPRB)                                             :: ZAIRDM1
REAL(KIND=JPRD) :: ZDENS_DP


! * Photolysis data:
!   ...correction factor for photo rates in function of sun-earth distance
REAL(KIND=JPRB)                          :: ZDISTFAC
!   ...absorbing species, lowest level to compute, solar flux
REAL(KIND=JPRB), DIMENSION(KLEV,NABSPEC) :: ZDENSA
INTEGER(KIND=JPIM)                       :: JBOTJ
REAL(KIND=JPRB), DIMENSION(MXWVN)        :: ZFBEAMR
INTEGER(KIND=JPIM)                       :: IDXDATE(1)

! * reaction rates strato

! BASCOE Conversions  : 1.255e+21 *.767e-17   from kg/m2 --> mol/cm2 --> DU
!REAL(KIND=JPRB), PARAMETER   :: ZTODU = 47296. ! O3 from kg/m2 --> DU
REAL(KIND=JPRB)                           :: ZCST  ! molec/cm2 -> DU
REAL(KIND=JPRB)                           :: ZVMRO3,ZSZA, ZHGT, ZH2O_ECMWF
REAL(KIND=JPRB)                           :: ZCOLO3_DU(KLON,0:KLEV)
INTEGER(KIND=JPIM), DIMENSION(KLON)       :: JBOT_STRATO     ! bottom level strato chemistry
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,NDISS):: ZAJVAL
REAL(KIND=JPRB), DIMENSION(NDISS)         :: ZAJVALH
REAL(KIND=JPRB), DIMENSION(NHET)          :: ZRHET
REAL(KIND=JPRB),PARAMETER                 :: ZSMALL_VMR=1.0E-30                        


!KPP related
INTEGER(KIND=JPIM),DIMENSION(20)         :: ICNTRL, ISTATUS
REAL(KIND=JPRD),   DIMENSION(20)         :: ZCNTRL, ZCNTRL_P, ZSTATE
REAL(KIND=JPRD),   DIMENSION(NREACT)     :: ZRCONST
REAL(KIND=JPRD),   DIMENSION(NVAR)       :: ZVAR
REAL(KIND=JPRD),   DIMENSION(NFIX)       :: ZFIX
INTEGER(KIND=JPIM)                       :: IERR

! Strato. PSC / aerosol related 
LOGICAL,            DIMENSION(KLON)        :: LL_PSC_POSSIBLE
INTEGER(KIND=JPIM), DIMENSION(KLON)        :: JBOT_PSC,JTOP_PSC
REAL(KIND=JPRB),    DIMENSION(KLON,KLEV,NBINS)   :: ZSA_SIZEDIST
REAL(KIND=JPRB),    DIMENSION(KLON,KLEV,NAER)    :: ZAER
REAL(KIND=JPRB),    DIMENSION(KLON,KLEV)  :: ZAER_INFO
INTEGER(KIND=JPIM), DIMENSION(KLON)       :: JTROPOP        ! top level tropo
INTEGER(KIND=JPIM), DIMENSION(9)          :: JHET_TRACER
INTEGER(KIND=JPIM), DIMENSION(2)          :: JPSC_TRACER

! Surface boundary conditions (BASCOE)
REAL(KIND=JPRB) , DIMENSION(KLON)           :: ZTENBC
INTEGER(KIND=JPIM)                          :: JBC
REAL(KIND=JPRB), DIMENSION(NLATBOUND_LBC-1) :: ZBCVAL      ! for every latitude band
INTEGER(KIND=JPIM)                          :: IMONTH_LBC(1), JLAT_LBC

! Variables used to compute altitude
REAL(KIND=JPRB)       :: ZPSURF_STD,ZSURF_H, ZTHKNESS,ZHGT_BASCOE(KLON,KLEV)

INTEGER(KIND=JPIM)                        :: IRANGE_TROPOP


! ------------------------------------------------------------------
#include "fcttim.func.h"
!-------------------------------------------------------------------
#include "bascoe_gs_liq.intfb.h"
! #include "bascoe_zenith_fct.intfb.h"
#include "bascoe_hetconst.intfb.h"
#include "bascoe_j_interp.intfb.h"
#include "bascoe_j_calc.intfb.h"
#include "sundistcorr.intfb.h"
#include "bascoe_psc_param.intfb.h"
#include "bascoe_psc_possible.intfb.h"
#include "bascoe_tropopause.intfb.h"
!#include "bascoe_wetdep.intfb.h"
! KPP code
#include "bascoe_kpp_rates.intfb.h"
#include "bascoe_v0_kpp_initialize.intfb.h"
#include "bascoe_kpp_integrator.intfb.h"
#include "bascoe_v0_kpp_update_cifs_conc.intfb.h"
#include "cifs_kpp_wlamch.intfb.h"
! CH4 boundary condition...
! #include "tm5_boundary_ch4.intfb.h"
!-----------------------------------------------------------------------
! chemistry scheme name - this will later also come from external input
IF (LHOOK) CALL DR_HOOK('CHEM_BASCOE',0,ZHOOK_HANDLE )
!-----------------------------------------------------------------------
! chemistry scheme name - this will later also come from external input
ASSOCIATE(LCHEM_JOUT=>YDCHEM%LCHEM_JOUT, NCHEM=>YDML_GCONF%YGFL%NCHEM, YCHEM=>YDML_GCONF%YGFL%YCHEM, &
 & YDRIP=>YDML_GCONF%YRRIP, LCHEM_BASCOE_JON=>YDCHEM%LCHEM_BASCOE_JON)
!-----------------------------------------------------------------------


! Preparation for kpp - solver

! Set kpp parameters to default, taken from cifs_kpp_IntParam module 
! See comments in Integrator module for
! a list of the defaults.

!VH RTOLS = RTOLS_G

RTOL(1:NVAR) = 0.05_JPRB * RTOLS_G

ATOL(1:NVAR) = 10._JPRB  !  was 1.e-16*cfactor before v3s04 or v3d06

ICNTRL(1:20) = 0_JPIM
ICNTRL(1) = IAUTONOM
! Change some parameters from the default to new values
! Select Integrator
!    ICNTRL(3)  -> selection of a particular method.
!   For Rosenbrock, options are:
!       = 0 :  default method is Rodas3
!       = 1 :  method is  Ros2
!       = 2 :  method is  Ros3 
!       = 3 :  method is  Ros4 
!       = 4 :  method is  Rodas3
!       = 5 :  method is  Rodas4

! ----------------------------------------------------------------------
!  Set Integrator input parameters. Values set in chem_IntParam.f90
! ----------------------------------------------------------------------
ICNTRL(3) = IROSMETH
!VH ICNTRL(4) = 200 ! set max. no of steps ?
ICNTRL(7) = 1 ! Currently no adjoint

ZCNTRL(1:20) = 0._JPRB
ZCNTRL(1) = HMIN
ZCNTRL(2) = PTSTEP
ZCNTRL(3) = HSTART

! Set range for tropopause 4 for 60 and 91 level version, 8 for 137 level version
IRANGE_TROPOP=4
IF (KLEV > 130) THEN 
  IRANGE_TROPOP=8
ENDIF

! Preferred time step - should be made flexible (see HSAVE_KPP)
!ZSTATE(Nhexit) = 0._JPRB

! Lat / Lon
DO JL=KIDIA,KFDIA
  ZLAT(JL)=(180.0_JPRB/RPI)*PGELAT(JL)
  ZLON(JL)=(180.0_JPRB/RPI)*PGELAM(JL)
ENDDO

! Initialize appropriate tracers for strat. heterogeneous and psc chemistry
JHET_TRACER(1)= IH2O
JHET_TRACER(2)= IN2O5
JHET_TRACER(3)= IHCL
JHET_TRACER(4)= IHOCL
JHET_TRACER(5)= ICLONO2
JHET_TRACER(6)= IHOBR
JHET_TRACER(7)= IHBR
JHET_TRACER(8)= IBRONO2
JHET_TRACER(9)= IHNO3

JPSC_TRACER(1)= IH2O
JPSC_TRACER(2)= IHNO3


! Time stuff:
! ...start date...
IYEAR0=NCCAA(NINDAT)
IMONTH0=NMM(NINDAT)
IDAY0=NDD(NINDAT)
! ...current date...
CALL UPDCAL (IDAY0,IMONTH0,IYEAR0,YDRIP%NSTADD,IDAY,IMONTH,IYEAR,ILMONTHS,-1)
IYYYYMM =IYEAR*100+IMONTH
IYYYYMMDD = IYYYYMM*100+IDAY
! ...day of year, number of days in year 
ZJUL1 = RJUDAT(IYEAR,1,1)
IJUL = NINT(RJUDAT(IYEAR,IMONTH,IDAY) - ZJUL1) + 1
IJULMAX = NINT(RJUDAT(IYEAR,12,31) - ZJUL1) + 1

! RHGMT: GMT time of model - between 0 and 86400.
! ZUT: Time of day in hours.
ZUT = YDRIP%RHGMT /3600.

! Compute correction factor for photo rates
CALL SUNDISTCORR(IMONTH,IDAY,ZDISTFAC)

! Compute here 'roundoff' number - there is a paralellization
! issue when calling WLMACH as part of kpp-code.
IF ( KSTEP == 0_JPIM ) THEN
    CALL CIFS_KPP_WLAMCH(ROUNDOFF_STORE , 'E')
ENDIF


IF (.NOT. LCHEM_BASCOE_JON) THEN
  !  1.1  Compute full level O3 TC, needed for offline J-rate computation

  ! BASCOE variant to compute overhead column [DU]
  !-----------------------------------------------------------------------
  ! Calculate overhead ozone columns *at* levels:  Zcst * SUM( vmr * delta_p )
  ! where delta_p is between levels and vmr are at mid-levels
  !-----------------------------------------------------------------------
  ZCST  = (1.0E-4/RG)*( RNAVO / (1.0E-3* RMD) ) * 1.0E3 / 2.687E19  ! molec/cm2 -> DU
  DO JL=KIDIA,KFDIA
    ! convert to mixing ratio
    ZVMRO3=MAX(PCEN(JL,1,IO3) / YCHEM(IO3)%RMOLMASS *RMD ,0._JPRB) 
    ! convert to DU
    ZCOLO3_DU(JL,1)  = ZCST * ZVMRO3 * ( PRSF1(JL,1) - 0._JPRB )
  ENDDO

  DO JLEV=2,KLEV
    DO JL=KIDIA,KFDIA
      ! convert to mixing ratio
      ZVMRO3=MAX(0.5_JPRB*(PCEN(JL,JLEV,IO3)+PCEN(JL,JLEV-1,IO3)) / YCHEM(IO3)%RMOLMASS *RMD ,0._JPRB) 
      ! convert to DU
      ZCOLO3_DU(JL,JLEV)=ZCOLO3_DU(JL,JLEV-1)+ZCST * ZVMRO3* ( PRSF1(JL,JLEV) - PRSF1(JL,JLEV-1))
    ENDDO
  ENDDO
ENDIF



! BASCOE way to compute model altitude, first surface model level:
ZPSURF_STD=101325._JPRB ! std p at surf (Pa)
DO JL=KIDIA,KFDIA

  IF( PRS1(JL,KLEV) < ZPSURF_STD ) THEN
    ZSURF_H = 7._JPRB*LOG( ZPSURF_STD / PRS1(JL,KLEV) )
  ELSE
    ZSURF_H=0.0_JPRB
  ENDIF
  ZTHKNESS = PTP(JL,KLEV)*287./9.806*LOG(PRS1(JL,KLEV)/PRSF1(JL,KLEV))
  ZHGT_BASCOE(JL,KLEV) = ZSURF_H + 1.E-3*ZTHKNESS
ENDDO

DO JLEV=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
     ZTHKNESS=0.5*(PTP(JL,JLEV+1)+PTP(JL,JLEV))*287./9.806*&
     &                 LOG(PRSF1(JL,JLEV+1)/PRSF1(JL,JLEV))
     ZHGT_BASCOE(JL,JLEV)=ZHGT_BASCOE(JL,JLEV+1)+1E-3*ZTHKNESS
  ENDDO
ENDDO

!--------
! Find tropopause level
CALL BASCOE_TROPOPAUSE(KIDIA,KFDIA,KLON,KLEV,PTP,PRSF1,PGEOH,ZLAT,JTROPOP)

! ------ GET SAD climatology - 

IF (KSTEP /= 0_JPIM ) THEN
  ! ZAER_INFO will be used in BASCOE_GS_LIQ
  DO JLEV=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAER_INFO(JL,JLEV)=PCEN(JL,JLEV,ISTRATAER)
    ENDDO
  ENDDO

ELSE
  ! ZAER_INFO will be updated in BASCOE_GS_LIQ
  DO JLEV=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAER_INFO(JL,JLEV)=0._JPRB
    ENDDO
  ENDDO
ENDIF

CALL BASCOE_GS_LIQ(KSTEP, IMONTH, KIDIA, KFDIA, KLON, KLEV,  JTROPOP, PRSF1, ZLAT, PTP, ZAER, ZSA_SIZEDIST, ZAER_INFO)

IF (KSTEP == 0_JPIM) THEN
  ! Create 'tendency' to arrive at aerosol field
  PTENC1(KIDIA:KFDIA,1:KLEV,ISTRATAER)=(ZAER_INFO(KIDIA:KFDIA,1:KLEV) - PCEN(KIDIA:KFDIA,1:KLEV,ISTRATAER))/ PTSTEP
ENDIF


! --------------------
!  Find range where PSC can be present
! --------------------

CALL BASCOE_PSC_POSSIBLE(KIDIA, KFDIA, KLON, KLEV,IJUL,PTP,PRSF1,ZLAT,&
                        & JTROPOP,LL_PSC_POSSIBLE,JTOP_PSC,JBOT_PSC) 

! ---------------------
!  Define range of levels to apply strato chemistry
!       (above 400hPa)

DO JL=KIDIA,KFDIA
  DO JK=KLEV,1,-1
    IF ( PRSF1(JL,JK) < 40000_JPRB ) THEN
      JBOT_STRATO(JL) = JK
      EXIT
    ENDIF
  ENDDO
ENDDO

! ---------------------
! initialize air densities (molec/cm3)

DO JLEV=1,KLEV
  DO JL=KIDIA,KFDIA
    ZDENS(JL,JLEV) = 7.24291e16_JPRB*PRSF1(JL,JLEV)/PTP(JL,JLEV)
  ENDDO
ENDDO

! ---------------------
!  Compute strato photolysis rates
!   but ONLY where needed, to avoid unnecessary processing:
!   - at same levels where strato chemistry is applied
!   - for solar zenith angle < 96 degrees)

ZAJVAL(KIDIA:KFDIA,1:KLEV,1:NDISS) = 0._JPRB

IF (LCHEM_BASCOE_JON) THEN
  ! Use online calculations...
  !
  !   ...obtain solflux for current date if daily solflux is requested...
  !
  IF (daily_solflux) THEN
    ! check if date in range
    !
    IF ( IYYYYMMDD < fbeam_dates(1) .OR. IYYYYMMDD > fbeam_dates(SIZE(fbeam_dates)) ) THEN
      CALL ABOR1(CL_MY_NAME//' error: date outside range of solflux dates read from file')
    ENDIF
    ! get closest date
    !
    IDXDATE = MINLOC( ABS(fbeam_dates - IYYYYMMDD ))
    IF ( fbeam_dates(IDXDATE(1)) /= IYYYYMMDD ) THEN
        WRITE(NULOUT,*) CL_MY_NAME//' warning solflux date ',fbeam_dates(IDXDATE(1)), &
        &                           ' differs from current date ',IYYYYMMDD
    ENDIF
    ZFBEAMR(:) = FBEAMR2D(IDXDATE(1),:)


  ELSE
    ZFBEAMR = FBEAMR
  ENDIF

  DO JL=KIDIA,KFDIA
    ZSZA = ACOS(PCSZA(JL))*180_JPRB/RPI

    !VH propose to take over BASCOE version (double check that it works for all (FC-)times !
    ! CALL BASCOE_ZENITH_FCT(ZLAT(JL),ZLON(JL),IJUL,IJULMAX,ZUT,ZSZA)

    IF( ZSZA < 96._JPRB ) THEN
      JBOTJ = JBOT_STRATO(JL)         ! lowest level to compute photolysis rates

      ! compute number densities profiles for the absorbing species
  
      DO JK = 1, JBOTJ
        ZDENSA(JK,AIR)   = ZDENS(JL,JK)
        ZDENSA(JK,O2ABS) = ZDENS(JL,JK) * 0.209_JPRB
        ZAIRDM1 = ZDENS(JL,JK) * RMD
        ZDENSA(JK,O3ABS) = PCEN(JL,JK,IO3)  / YCHEM(IO3)%RMOLMASS  *ZAIRDM1
        ZDENSA(JK,NOABS) = PCEN(JL,JK,INO)  / YCHEM(INO)%RMOLMASS  *ZAIRDM1
        ZDENSA(JK,CO2ABS)= PCEN(JL,JK,ICO2) / YCHEM(ICO2)%RMOLMASS *ZAIRDM1
        ZDENSA(JK,NO2ABS)= PCEN(JL,JK,INO2) / YCHEM(INO2)%RMOLMASS *ZAIRDM1
      ENDDO

      CALL BASCOE_J_CALC( JBOTJ, ZHGT_BASCOE(JL,1:JBOTJ), ZSZA, PTP(JL,1:JBOTJ), &
                      &   ZDENSA(1:JBOTJ,1:NABSPEC), PALB(JL), ZFBEAMR, ZAJVAL(JL,1:JBOTJ,1:NDISS) )
    ENDIF
  ENDDO
ELSE
  ! use offline (faster lookup-table) approach
  DO JL=KIDIA,KFDIA
    ZSZA = ACOS(PCSZA(JL))*180_JPRB/RPI

    !VH propose to take over BASCOE version (double check that it works for all (FC-)times !
    ! CALL BASCOE_ZENITH_FCT(ZLAT(JL),ZLON(JL),IJUL,IJULMAX,ZUT,ZSZA)

    IF( ZSZA < 96._JPRB ) THEN
      JBOTJ = JBOT_STRATO(JL)         ! lowest level to compute photolysis rates
  
      ! compute number densities profiles for the absorbing species
  
      DO JK = 1, JBOTJ
        !VH ZHGT = PGEOH(JL,JK-1)  * ZRGI *1e-3 ! height in km
        ! ZHGT = 0.5*(PGEOH(JL,JK-1) + PGEOH(JL,JK))  * ZRGI *1e-3 ! height in km
        ! 
        ZHGT = ZHGT_BASCOE(JL,JK)

        ! Set maximimum height to ~110 km altitude 
        ZHGT = MIN(ZHGT,109.9_JPRB)
        ! ZCOLO3_DU: o3 overhead column in DU, converted from kg/m2
        ! compute photolysis rates (ZAJVAL) 
        ZAJVALH(1:NDISS)=0._JPRB
        CALL BASCOE_J_INTERP( ZSZA, ZHGT  , ZCOLO3_DU(JL,JK), ZAJVALH(1:NDISS) )
        ZAJVAL(JL,JK,1:NDISS)=ZAJVALH(1:NDISS)
      ENDDO
    ENDIF
  ENDDO
ENDIF

! adjust photolysis rate to sun earth distance variation
ZAJVAL(KIDIA:KFDIA,1:KLEV,1:NDISS) = ZAJVAL(KIDIA:KFDIA,1:KLEV,1:NDISS) * ZDISTFAC

!-----------------------------------------------------------------------
!2.0 loop over levels for solving chemistry
DO JK=1,KLEV
 
!2.1 Loop over lon/lat
  DO JL=KIDIA,KFDIA
     
    !*  Air density mutiplied with RMD (dry air molar mass) for efficiency
    ZAIRDM1 = ZDENS(JL,JK) * RMD


    !*  convert tracer concentrations from kg/kg to molec/cm3
    DO JT=1,NCHEM
      !*     assure positivity for initial concentrations
      ZCVM0(JL,JT) = MAX(PCEN(JL,JK,JT) / YCHEM(JT)%RMOLMASS *ZAIRDM1 ,ZSMALL_VMR) 
      !*     initialize final (ZCVM) and intermediate (ZCVM1) concentrations
      !*     not strictly needed for kpp, but necessary to initialize tropospheric conc.
      !*     in case trop. chemistry is switched off
      ZCVM(JL,JT)  = ZCVM0(JL,JT)
    ENDDO

!*   Overwrite ECMWF H2O in troposphere as well as tropopause region
    IF (JK > JTROPOP(JL)-IRANGE_TROPOP) THEN
!* Overwrite only in tropopause region, where Strat chemistry is solved:
      ZCVM0(JL,IH2O)=PQP(JL,JK)/ YCHEM(IH2O)%RMOLMASS * ZAIRDM1
      ZCVM(JL,IH2O)=ZCVM0(JL,IH2O)
    ENDIF
  

   ! Try Only solve stratospheric chemistry from tropopause onwards?!

   !VH IF (JK < JTROPOP(JL)+4) THEN

    ! * solve chemistry in strato domain
    IF ( JK <= JBOT_STRATO(JL) ) THEN

      !* Fixed concentrations
      ZFIX(1) = 0.209*ZDENS(JL,JK)    ! O2 number density
      ZFIX(2) = 0.781*ZDENS(JL,JK)    ! N2 number density


      ! ----------------------------------------------------------------------
      !  Compute heterogeneous reaction rates (ZRHET)
      ! ----------------------------------------------------------------------

      IF (YDCHEM%LCHEM_BASCOE_HETCHEM) THEN

        CALL BASCOE_HETCONST(YDML_GCONF%YGFL,JHET_TRACER,PTP(JL,JK),PRSF1(JL,JK),ZDENS(JL,JK), &
          & LL_PSC_POSSIBLE(JL),JTOP_PSC(JL), JBOT_PSC(JL),JK,ZCVM0(JL,1:NCHEM), &
          & ZSA_SIZEDIST(JL,JK,1:NBINS),ZAER(JL,JK,1:NAER),PTSTEP,ZRHET)

      ELSE

        ! in case of emergency switch off: 
        ZRHET(:) = 0._JPRB


      ENDIF


      ! Output EXTRA fields (handle with care)
      ! IF (LCHEM_JOUT) THEN
      !  POUT(JL,JK,2) =  ZAJVAL(JL,JK,J_O3_O1D)
      !  POUT(JL,JK,3) =  ZAJVAL(JL,JK,J_NO2)
      ! ENDIF
   

       ! Initialize concentrations to KPP (fill VAR) 
       CALL BASCOE_V0_KPP_INITIALIZE(YDML_GCONF%YGFL,ZCVM0(JL,1:NCHEM),ZVAR(1:NVAR))

       ! Initialize all kpp reaction rates... (fill ZRCONST) 
       CALL BASCOE_KPP_RATES(ZAJVAL(JL,JK,1:NDISS),ZRHET, PTP(JL,JK), ZDENS(JL,JK), ZRCONST, ZVAR(1:NVAR))  

       ! starting value for integration time step
       ZCNTRL_P=ZCNTRL
       !VH Maybe to be switched on...    ZCNTRL_P(3) = ZHSAVE_KPP(JL,JLEV)

       !VH - Overwrite initial timestep  - not needed here
       !ZCNTRL_P(3) = MIN(PTSTEP, ZHSTART)


       ! ----------------------------------------------------------------------
       !  Now call the chem box solver
       ! ----------------------------------------------------------------------
       ! Call kpp integrator... (provide 'ZVAR' and 'ZRCONST' !)
       ZDENS_DP=ZDENS(JL,JK)
       CALL BASCOE_KPP_INTEGRATOR(0._JPRD, PTSTEP, ICNTRL,ZCNTRL_P, ISTATUS,ZSTATE,IERR,&
          & ZVAR,ZFIX,ZRCONST, ZDENS_DP)
       !- Filter error due to bad concentrations
       IF (IERR>0) THEN

         !No errors: update concentrations...
         CALL BASCOE_V0_KPP_UPDATE_CIFS_CONC(YDML_GCONF%YGFL,ZVAR(1:NVAR),ZCVM(JL,1:NCHEM))

       ELSEIF( IERR==-9) THEN
         WRITE(NULERR,'(a)') '     ZVAR below are out of range:'
         DO JT= 1, NVAR
           IF( ZVAR(JT)/ZDENS(JL,JK) > VMR_BAD_LARGE) THEN
             WRITE(NULERR,'(a,2(i5),a,es12.5)') '  vmr-idx ',JT,JK, ' ; reached value: ', ZVAR(JT)/ZDENS(JL,JK)
           ENDIF
         ENDDO
         WRITE(NULERR,*) '  -> chem integrator skipped'
        ENDIF

    ENDIF ! Within stratosphere

    !VH make sure that within troposphere and around tropopause
    !VH we get the tendencies to arrive at ECMWF-H2O
    !VH so impose initial and final concentration fields.
    !    IF (JK > JTROPOP(JL)-8) THEN
    !Test for lower altitude wrt H2O fix
    IF (JK > JTROPOP(JL)-IRANGE_TROPOP) THEN
      ZCVM0(JL,IH2O)=PCEN(JL,JK,IH2O) / YCHEM(IH2O)%RMOLMASS *ZAIRDM1

!* Nudge towards ECMWF-H2O, but don't enforce it...
      ZH2O_ECMWF=PQP(JL,JK)/ YCHEM(IH2O)%RMOLMASS * ZAIRDM1
      IF (JK < JTROPOP(JL)+IRANGE_TROPOP) THEN
        ! 3 hour decay time near tropopause
        ZCVM(JL,IH2O) = ZCVM0(JL,IH2O)+(1._JPRB-EXP(-PTSTEP/(10800._JPRB)))*(ZH2O_ECMWF-ZCVM0(JL,IH2O))
      ELSE
        ! 1 day decay time in toposphere
        ZCVM(JL,IH2O) = ZCVM0(JL,IH2O)+(1._JPRB-EXP(-PTSTEP/(86400._JPRB)))*(ZH2O_ECMWF-ZCVM0(JL,IH2O))
      ENDIF
    ENDIF

  ENDDO ! loop over JL

!4.0 Compute sedimentation in stratosphere... (Here? Or before chemistry?)
  CALL BASCOE_PSC_PARAM(YDML_GCONF%YGFL,KIDIA,KFDIA,KLON,JPSC_TRACER,PTSTEP,JTROPOP,JK,PTP(1:KLON,JK),&
  & PRSF1(1:KLON,JK),ZCVM(1:KLON,1:NCHEM))

!4.5 Apply extremely simple wet dep, in case one wants to run C-IFS (BASCOE)
!    as similar as possible as BASCOE
!  CALL BASCOE_WETDEP( KIDIA, KFDIA, KLON,PTSTEP, PRSF1(:,JK), ZCVM(:,1:NCHEM))

!5.0 convert concentration tendencies to mass mixing ratio
  DO JL=KIDIA,KFDIA
!*  Air density mutiplied with RMD (dry air molar mass) for efficiency  
    ZAIRDM1 = ZDENS(JL,JK) * RMD
    DO JT=1,NCHEM 
      IF (JT /= ISTRATAER) THEN 
        PTENC1(JL,JK,JT) = (ZCVM(JL,JT) - ZCVM0(JL,JT)) * YCHEM(JT)%RMOLMASS / (ZAIRDM1 * PTSTEP)
      ENDIF  
    ENDDO
  ENDDO

ENDDO ! loop over levels


! ----------------------------------------------------------------------
! Simple boundary conditions at surface for a list of BASCOE species - 
! Check this when introducing emissions!
!
! YC 20180830: allow monthly varying, latitude band dependent LBC
!
! get the index of current month
!   assume current month in LBC coordinates, otherwise nearest neighbour
!   *this could be refined, as trend and seasonality is lost*
IMONTH_LBC = MINLOC( ABS(MONTH_LBC - IYYYYMM) )

DO JB=1,NBC
  JBC    = BASCOE_BC(JB)
  ! Select appropriate BC and convert to mass ratio [kg/kg]
  !
  ZBCVAL = VALUES_LBC(JB, IMONTH_LBC(1), :) * YCHEM(JBC)%RMOLMASS / RMD

  DO JL=KIDIA,KFDIA
    DO JLAT_LBC = 1, NLATBOUND_LBC-1
      IF( ZLAT(JL) >= XLATBOUND_LBC(JLAT_LBC) .AND. ZLAT(JL) <= XLATBOUND_LBC(JLAT_LBC+1) ) EXIT
    ENDDO
    ZTENBC(JL) = ZBCVAL(JLAT_LBC) - PCEN(JL,KLEV,JBC)
    PTENC1(JL,KLEV,JBC) =  ZTENBC(JL)  / PTSTEP
    ! store this special LBC budget contribution in POUT(:,x,x)
    POUT(JL,JB+1,1)=PTSTEP *PTENC1(JL,KLEV,JBC)*PDELP(JL,KLEV) / RG
  ENDDO
ENDDO


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_BASCOE',1,ZHOOK_HANDLE )
END SUBROUTINE CHEM_BASCOE
