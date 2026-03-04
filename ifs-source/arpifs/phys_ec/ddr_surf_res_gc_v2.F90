! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DDR_SURF_RES_GC_V2 (KCHEM_DRYDEP, PTS, PITM, PFRSO, &
                        &  PRSTO, KTILE, KTILE_NOWET, KVEG_GC, PLAI, &
                        &  CDNMS, PXCHEN, PXCHENXP, PXDIMO, PXCF0, &
                        &  PCSZA, PCFRAC,PXMW, PRESS,KDEBUG, &
                        &  PWRC)
!!    ---------
!!         *COMP_SURF_RES_GC_V2* IS CALLED FROM *DEPVEL_GC*
!!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to cumpute the surface resistances for 
!!    the chemical species
!!    Computation is performed for a given date, forecast start time, 
!!    forecast valid time, and cells of a given domain. 
!!   
!!    METHOD
!!    ------
!!    See [Seinfeld et Pandis, 1998, "Atmospheric Chemistry and
!!    Physics", chap. 19, pp.958-996] based primarily on [Wesely, 1989]. 
!!    Note that in [Wesely,1989] only
!!    14 gaseous species are listed, against 17 in [Seinfeld et Pandis, 1998].
!!
!!    REFERENCE
!!    ---------
!!    See [Seinfeld et Pandis, 1998, "Atmospheric Chemistry and
!!    Physics", chap. 19, pp.958-996] based primarily on [Wesely, 1989]. 
!!
!!
!! INPUTS:
!! -------
!! PTS        : Surface temperature
!! PITM       : Land/sea mask
!! PFRSO      : Solar radiation
!! PRSTO      : Stomatal resistances for water vapor(s.m-1)
!! KTILE      : IFS TILE (1... 9) 
!! KTILE_NOWET: IFS TILE (1... 9) no wet skin - reassigned to vegetation or bare ground
!! KVEG_GC    : Vegetation type
!! CDNMS       : Species name
!! PXCHEN      : Henry's specific constant
!! PXCHENXP    : Henry's specific constant exponent
!! PXDIMO      : Molecular diffusivities
!! PXCF0       : Reactivity
!!
!! OUTPUTS:
!! -------
!! ZWRC       : Surface resistances (s.m-1)
!!
!! LOCAL:
!! -------
!!
!!    AUTHOR
!!    ------
!!    D. Finch 
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    original    November 2020


USE PARKIND1 ,ONLY : JPIM, JPRB
USE YOMCST, ONLY : R
USE DRYDEP_PAR_GC  ,ONLY :    & ! GEOS-Chem specific parameters
&     IRI, IRLU, IRAC, IRGSS, IRGSO, IRCLS, IRCLO
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! USE YOMLUN, ONLY: NULERR

! !PRIVATE MEMBER FUNCTIONS:
!PRIVATE :: BIO_FIT

IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN)    :: KCHEM_DRYDEP ! Extra switch
  REAL (KIND=JPRB), INTENT(IN)      :: PTS
  REAL (KIND=JPRB), INTENT(IN)      :: PITM
  REAL (KIND=JPRB), INTENT(IN)      :: PFRSO
  REAL (KIND=JPRB), INTENT(IN)      :: PRSTO
  INTEGER (KIND=JPIM), INTENT(IN)   :: KVEG_GC,KTILE,KTILE_NOWET
  REAL (KIND=JPRB), INTENT(IN)      :: PLAI
  CHARACTER (LEN=*), INTENT(IN)     :: CDNMS
  REAL (KIND=JPRB), INTENT(IN)      :: PXCHEN
  REAL (KIND=JPRB), INTENT(IN)      :: PXCHENXP
  REAL (KIND=JPRB), INTENT(IN)      :: PXDIMO ! Ratio between D_H2O / D_I MOLECULAR diffusivity
  REAL (KIND=JPRB), INTENT(IN)      :: PXCF0
  REAL (KIND=JPRB), INTENT(IN)      :: PCSZA
  REAL (KIND=JPRB), INTENT(IN)      :: PCFRAC
  REAL (KIND=JPRB), INTENT(IN)      :: PXMW ! Molecular weight (in kg/mol)
  REAL (KIND=JPRB), INTENT(IN)      :: PRESS ! Surface pressure ! check - currently mid level pres
  INTEGER(KIND=JPIM), INTENT(IN)    :: KDEBUG ! debug output info if needed
  REAL (KIND=JPRB), INTENT(OUT)     :: PWRC

  !! LOCAL VARIABLES
  REAL (KIND=JPRB)                  :: ZWRCVEG, ZWRCNVEG
  INTEGER (KIND=JPIM)               :: IKVEG_WE ! Mapping to Wesely types (?)
  REAL (KIND=JPRB)                  :: ZI
  REAL (KIND=JPRB)                  :: ZCLO, ZCLS, ZAC, ZGSO, ZGSS
  REAL (KIND=JPRB)                  :: ZDC, ZLU, ZTCOR,ZTS
  REAL (KIND=JPRB)                  :: ZHENRY
  REAL(KIND=JPRB)                   :: ZTR,ZKA1,ZKW
  REAL (KIND=JPRB)                  :: ZGFACT,ZGFACI ! Insolation & temperature factors
  REAL(KIND=JPRB)                   :: ZD_TEMP1, ZD_TEMP2, ZD_TEMP3, ZD_TEMP4 ! Temp calculation steps
  REAL(KIND=JPRB)                   :: ZIXX, ZLUXX, ZGSX, ZCLX ! Intermediate calc steps
  REAL(KIND=JPRB),PARAMETER         :: ZHPLUS =3.16227E-6_JPRB ! [H+]. For now assume soil pH=5.5
  REAL(KIND=JPRB),PARAMETER         :: ZXMWH2O = 18E-3_JPRB ! Molecula weight of H20 (kg/mol)

  REAL(KIND=JPHOOK) ::  ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC_V2',0,ZHOOK_HANDLE)


  ! In a first step, all cells are considered as water cells to optimise the
  ! execution time for the global domain

ZTCOR = 0.0_JPRB
ZWRCNVEG = 0.0_JPRB
ZWRCVEG = 0.0_JPRB

PWRC=-99999.9_JPRB
! PWRC_GC=-99999.9_JPRB

! surface independent treatment of HNO3, H2O2 
 IF  (  ANY( (/ "HNO3","H2O2" /)   ==  CDNMS ) ) THEN 

  ! no temperature correction over sea / water 
  !  ZTCOR = 0.0_JPRB
  !  IF (PTS < 271.0_JPRB .AND. ANY( KTILE == (/ 2 , 4, 5, 6, 7, 8/) ))  THEN
  !    ZTCOR = MIN(PPRMAX, 1000.0_JPRB * EXP(-PTS + 269.0_JPRB) )
  !  ENDIF 
  !  PWRC = MAX(10.0_JPRB, ZTCOR)
    PWRC = 10.0

ELSE 
! all other species 
    !Wesely [1989] gives ZTCOR = 1000.0*EXP(-TEMPC-4.0) (here PTS is in K)
    ZTS=MAX(PTS,260._JPRB)
    ZTCOR = 1000.0_JPRB * EXP(-ZTS + 269.15_JPRB)

    !      1.2  Adjusts Henry's constant according to species and surface
    !             temperature. Later, could also vary with pH.
    ZHENRY = PXCHEN * EXP(PXCHENXP * ( 1.0_JPRB / PTS  - 1.0_JPRB / 298.0_JPRB   )   )
    IF  ( CDNMS == "NH3" ) THEN
      ZTR=(1.0_JPRB/PTS-1.0_JPRB/298.0_JPRB)
      ZKA1=1.7E-5_JPRB *EXP( -450._JPRB*ZTR )
      ZKW =1.0E-14_JPRB*EXP( -6718._JPRB*ZTR )
      ZHENRY = ZHENRY*(1._JPRB + ZKA1*ZHPLUS/ZKW)
    ENDIF



    IKVEG_WE=KVEG_GC
    ! If the IFS tile is snow or ice
    IF (  ANY( (/ 2, 5, 7/)   ==  KTILE ) ) THEN
        IKVEG_WE = 1 ! Re-assign dry dep type to snow parameter
    ENDIF

    !** If ZRLU is '9999' it means there are no cuticular surfaces
    !** on which to deposit so we impose a very large value for ZLU.
    IF ( IRLU(IKVEG_WE) >= 9999 .OR. PLAI <= 0.0 ) THEN
        ZLU = 1.e+6_JPRB
    ELSE
    ! CORRECT FOR SEASON BY USING LAI - cuticular resistances are per unit area of leaf
        ZLU = IRLU(IKVEG_WE) / PLAI
        ! Additional resistance at low temperatures.
        ! Limit increase to a factor of 2.
        ! V. Shah 23 Oct 18
        ! Ref Jaegle et al. 2018
        ZLU = MIN( ZLU + ZTCOR, 2._JPRB * ZLU)
    ENDIF

    ZAC = MAX(1.0_JPRB, IRAC(IKVEG_WE))
    IF (ZAC >= 9999.0_JPRB) ZAC  = 1.E+12_JPRB

    ZGSS = MAX(1.0_JPRB, IRGSS(IKVEG_WE))
    ! Additional resistance at low temperatures.
    ! Limit increase of ZGSS, ZGSO, ZCLS, & ZCLO to a factor of 2.
    ! Ref Jaegle et al. 2018
    ZGSS = MIN(ZGSS + ZTCOR, 2._JPRB * ZGSS)
    IF (ZGSS >= 9999.0_JPRB) ZGSS = 1.E+12_JPRB

    ZGSO = MAX(1.0_JPRB, IRGSO(IKVEG_WE))
    ZGSO = MIN(ZGSO + ZTCOR, 2._JPRB * ZGSO)
    IF (ZGSO >= 9999.0_JPRB) ZGSO = 1.E+12_JPRB

    ZCLS = IRCLS(IKVEG_WE)
    ZCLS = MIN(ZCLS + ZTCOR, 2._JPRB * ZCLS)
    IF (ZCLS >= 9999.0_JPRB) ZCLS = 1.E+12_JPRB

    ZCLO = IRCLO(IKVEG_WE)
    ZCLO = MIN(ZCLO + ZTCOR, 2._JPRB * ZCLO)
    IF (ZCLO >= 9999.0_JPRB) ZCLO = 1.E+12_JPRB


    !compute stomatal resistances
    IF (KCHEM_DRYDEP == 3) THEN
      ! GEOS-Chem type implementation
      ZI = IRI(IKVEG_WE)
      IF (ZI >= 9999.0_JPRB) THEN
          ZI = 1.E+12_JPRB
      ELSE
      ! Adjust stomatal resistances for insolation (ZGFACI) and temperature (ZGFACT):
      ! Temperature adjustment is from Wesely [1989], equation (3).
  
          ZGFACT = 100._JPRB
          IF (PTS > 273.16_JPRB .AND. PTS < 313.14_JPRB) THEN
              ZGFACT = 400_JPRB/(PTS-273.15)/(313.15 - PTS)
          ENDIF
          ZGFACI = 100._JPRB
          ! If sun radiationa and LAI above zero then apply stomatal correction
          IF ((PFRSO > 0._JPRB) .AND. (PLAI > 0._JPRB)) THEN
              ! COEFF = Polynomial coefficients for dry deposition
              ! LAI = Leaf area index
              ! PCSZA / SUNCOS = Cosine of solar zenith angle
              ! PCFRAC = cloud cover fraction
              ZGFACI = 1._JPRB / BIO_FIT(PLAI, PCSZA, PCFRAC)
          ENDIF
          ZI = ZI * ZGFACT * ZGFACI
      ENDIF
      ZIXX = ZI * (DIFFG(PTS,PRESS,ZXMWH2O))/(DIFFG(PTS,PRESS,PXMW)) &
      & + 1._JPRB / (ZHENRY/3000._JPRB + 100._JPRB * PXCF0)

    ELSEIF(KCHEM_DRYDEP ==4) THEN
    
    IF (  ANY( (/ 4, 6, 7 /)   ==  KTILE_NOWET ) ) THEN
      ! Bulk canopy stomatal and mesophyll resistances eq. 19.25 from Seinfeld and Pandis
      ! rely on CTESSEL stomatal resistance for H2O (PRSTO)
      ZIXX = PRSTO * PXDIMO + 1.0_JPRB / (100.0_JPRB * PXCF0 + 3.3E-4_JPRB * ZHENRY)
    ELSE
      ! Very high resistance..
      ZIXX = 1E5_JPRB
    ENDIF
    
    ELSE 
      CALL ABOR1("KCHEM_DRYDEP option not supported") 
    ENDIF


    !Mixing forced by buoyant convection (ZDC)
    ZDC = 100.0_JPRB * (1.0_JPRB + 1000._JPRB / (PFRSO + 10.0_JPRB))


    ZLUXX = 1E+12_JPRB
    IF (ZLU < 9999._JPRB) THEN
        ZLUXX = ZLU/(ZHENRY/1.E+5_JPRB + PXCF0)
    ENDIF

    ! total ground conductance, combining sulfur and ozone-type 
    ZGSX = 1.0e+0_JPRB/(ZHENRY/1.E+5_JPRB/ZGSS + PXCF0/ZGSO)

    ! total lower canopy conductance
    ZCLX = 1.0e+0_JPRB/(ZHENRY/1.E+5_JPRB/ZCLS + PXCF0/ZCLO)

    ZD_TEMP1 = 1._JPRB / ZIXX
    ZD_TEMP2 = 1._JPRB / ZLUXX
    ZD_TEMP3 = 1._JPRB / (ZAC + ZGSX)
    ZD_TEMP4 = 1._JPRB / (ZDC + ZCLX)

    ! IF (KDEBUG==1) WRITE(NULERR,'(a10,6es12.5)')'WRC contrib=',ZIXX, ZLUXX, ZAC, ZGSX,ZDC,ZCLX
    PWRC = 1.0e+0_JPRB/(ZD_TEMP1 + ZD_TEMP2 + ZD_TEMP3 + ZD_TEMP4)

 ENDIF
 ! If ozone and water then limit between 1 and 999999
IF ((CDNMS == "O3") .AND. (IKVEG_WE == 11)) THEN
    PWRC = MAX(1._JPRB, MIN(PWRC,999999._JPRB))
ELSE
    PWRC = MAX(1._JPRB, MIN(PWRC,9999._JPRB))
ENDIF

! Set Rc for strong acids (HNO3,HCl,HBr) to 1 s/m
! V. Shah (23 Oct 18)
! Ref. Jaegle et al. 2018, cf. Erisman,van Pul,Ayers 1994
!IF ( ZHENRY > 1.E+10_JPRB ) THEN
!    PWRC = 1._JPRB
!ENDIF

! Resistances limited to a maximum value (1e5)
!IF (PWRC > PPRMAX) PWRC = PPRMAX

IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC_V2',1,ZHOOK_HANDLE)
CONTAINS

FUNCTION DIFFG( PTS, PRESS, PXM) 
    ! Calculate the molecular diffustivity based on temperature and pressure


    REAL (KIND=JPRB), INTENT(IN)       :: PTS     ! Temp
    REAL (KIND=JPRB), INTENT(IN)       :: PRESS   ! Pres
    REAL (KIND=JPRB), INTENT(IN)       :: PXM     ! molecular mass
    REAL (KIND=JPRB)                   :: DIFFG
    
    REAL(KIND=JPRB), PARAMETER  :: ZXMAIR  = 28.8E-3_JPRB ! Moist air molec wt
    REAL(KIND=JPRB), PARAMETER  :: ZRADAIR = 1.2E-10_JPRB
    REAL(KIND=JPRB), PARAMETER  :: ZRADX   = 1.5E-10_JPRB
    REAL(KIND=JPRB), PARAMETER  :: ZAVOGAD = 6.023E+23_JPRB ! Avogadros constant
    REAL(KIND=JPRB), PARAMETER  :: ZPI     = 3.1415926535897932_JPRB ! Pi
    !! Local variables
    REAL(KIND=JPRB)             :: ZAIRDEN, ZZ, ZDIAM, ZFRPATH, ZSPEED 

    REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC_V2:DIFFG',0,ZHOOK_HANDLE)

    ! Calculate molecular diffusivity of H20 for calculation
    ! Air density [molec/m3]
    ZAIRDEN = ( PRESS * ZAVOGAD ) / ( R * PTS )
    ! DIAM is the collision diameter for gas X with air.
    ZDIAM   = ZRADX + ZRADAIR
    !Calculate the mean free path for gas X in air:
    ! eq. 8.5 of Seinfeld [1986];
    ZZ      = PXM  / ZXMAIR
    ZFRPATH = 1._JPRB /( ZPI * SQRT( 1._JPRB + ZZ ) * ZAIRDEN * ( ZDIAM**2_JPIM ) )

    ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
    ZSPEED  = SQRT( 8_JPRB * R * PTS / ( ZPI * PXM ) )

    ! Calculate diffusion coefficient of gas X in air;
    ! eq. 8.9 of Seinfeld [1986]
    DIFFG = ( 3_JPRB * ZPI / 32_JPRB ) * ( 1_JPRB + ZZ ) * ZFRPATH * ZSPEED

    IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC_V2:DIFFG',1,ZHOOK_HANDLE)

END FUNCTION DIFFG

FUNCTION BIO_FIT( PPXLAI, PPSUNCOS, PPCFRAC) 
    ! Function to compute the light correction for surface deposition

    USE DRYDEP_PAR_GC, ONLY: COEFF, NPOLY      ! Baldocchi drydep coefficients

    REAL (KIND=JPRB), INTENT(IN)       :: PPXLAI     ! Leaf area index [cm2/cm2]
    REAL (KIND=JPRB), INTENT(IN)       :: PPSUNCOS   ! Cosine( Solar Zenith Angle )
    REAL (KIND=JPRB), INTENT(IN)       :: PPCFRAC    ! Cloud fraction [unitless]

    REAL (KIND=JPRB)                   :: BIO_FIT ! Resultant light correction

    INTEGER (KIND=JPIM), PARAMETER     :: IKK = 4_JPIM

    REAL (KIND=JPRB), DIMENSION (IKK)  :: ZTERM
    REAL (KIND=JPRB), DIMENSION (NPOLY):: ZREALTERM
    REAL (KIND=JPRB)                   :: ZLAI,ZCFRAC,ZSUNCOS
    INTEGER (KIND=JPIM)                :: IK, IK1, IK2, IK3

    REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

    !    Wang, Y., D.J. Jacob, and J.A. Logan, "Global simulation of tropospheric
    !     O3-NOx-hydrocarbon chemistry, 1. Model formulation", J. Geophys. Res.,
    !     103/D9, 10,713-10,726, 1998.

    IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC_V2:BIO_FIT',0,ZHOOK_HANDLE)

    ! NPOLY = SIZE(COEFF)

    ! Make sure the variables are within normal range
    ZLAI = MIN(PPXLAI,11._JPRB) ! Do not have a reference as to why this is set to 11
    ZSUNCOS = MIN(PPSUNCOS, 1._JPRB)
    ZCFRAC = MIN(PPCFRAC,1._JPRB)

    ! Make sure each variable is above a minimum
    ZLAI = MAX(ZLAI, 0.2_JPRB)
    ZSUNCOS = MAX(ZSUNCOS,0.05_JPRB)
    ZCFRAC = MAX(ZCFRAC,0._JPRB)

    ! Scaling factor LAI
    ZLAI = ZLAI / 11._JPRB

    ZTERM(1) = 1._JPRB
    ZTERM(2) = ZLAI
    ZTERM(3) = ZSUNCOS
    ZTERM(4) = ZCFRAC

    IK = 0
    DO IK3 = 1, IKK
        DO IK2 = IK3, IKK
            DO IK1 = IK2, IKK
                IK = IK + 1
                ZREALTERM(IK) = ZTERM(IK1) * ZTERM(IK2) * ZTERM(IK3)
            ENDDO
        ENDDO
    ENDDO

    BIO_FIT = 0._JPRB
    DO IK = 1, NPOLY
        BIO_FIT = BIO_FIT + COEFF(IK) * ZREALTERM(IK)
    ENDDO

    IF (BIO_FIT < 0.1_JPRB) BIO_FIT = 0.1_JPRB

    IF (LHOOK) CALL DR_HOOK('DDR_SURF_RES_GC_V2:BIO_FIT',1,ZHOOK_HANDLE)

END FUNCTION BIO_FIT

END SUBROUTINE DDR_SURF_RES_GC_V2
