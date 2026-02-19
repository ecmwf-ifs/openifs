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

MODULE YOMFA

USE PARKIND1  ,ONLY : JPIM
USE TYPE_FADS, ONLY : FAD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    GRIB packing options / Options du compactage GRIB

!     NVGRIB  :  Level of GRIB packing (0:no packing, 1:v0 GRIB, 2:v0mod GRIB)
!                Niveau de codage GRIB (0:pas de cod., 1:v0 GRIB, 2:v0mod GRIB)
!     NBITPG  :  Number of bits to code each grid point
!                Nombre de bits pour coder les points de grille
!     NBITCS  :  Number of bits to code each spectral coefficient
!                Nombre de bits pour coder les coefficients spectraux
!     NSTRON  :  Non-packed sub-truncation
!                Niveau de sous-troncature non compactee
!     NPULAP  :  Laplacian power
!                Puissance de Laplacien
!     NVERBOSE:  FA software verbosiry (0,1 or 2)
!     LSUPPDATE: Read input file date each time suppdate is called
!     JPXTRO   : Maximum truncation
!     JPXLAT   : Maximim number of latitudes
!     JPXNIV   : Maximum number of levels
!     JPNXCA   : Maximum number of "cadres"
!     LEXTERN  : Write GRIB1 fields in an external file

INTEGER(KIND=JPIM) :: NVGRIB
INTEGER(KIND=JPIM) :: NBITPG
INTEGER(KIND=JPIM) :: NBITCS
INTEGER(KIND=JPIM) :: NSTRON
INTEGER(KIND=JPIM) :: NPULAP
INTEGER(KIND=JPIM) :: NVERBOSE
LOGICAL            :: LSUPPDATE
CHARACTER(LEN=64)  :: CMODEL   = ''
INTEGER(KIND=JPIM) :: NIDCEN   = 85
LOGICAL            :: LEXTERN  = .TRUE.


!     Fields descriptors

TYPE(FAD) :: YFAOROG ! Surface geopotential
TYPE(FAD) :: YFASP   ! Surface pressure
TYPE(FAD) :: YFAPSI  ! Velocity potential
TYPE(FAD) :: YFAKHI  ! Stream function
TYPE(FAD) :: YFAUGEO ! U-Geographical wind
TYPE(FAD) :: YFAVGEO ! V-Geographical wind
TYPE(FAD) :: YFAT    ! Temperature
TYPE(FAD) :: YFAPD   ! Pressure departure (NH)
TYPE(FAD) :: YFAVD   ! Vertical divergence (NH)
TYPE(FAD) :: YFAQ    ! Specific humidity
TYPE(FAD) :: YFAL    ! Liquid water
TYPE(FAD) :: YFALCONV! Liquid water (CONV.PART)
TYPE(FAD) :: YFAICONV! Ice    water (CONV.PART)
TYPE(FAD) :: YFARCONV! Rain         (CONV.PART)
TYPE(FAD) :: YFASCONV! Snow         (CONV.PART)
TYPE(FAD) :: YFAI    ! Ice
TYPE(FAD) :: YFALRAD ! Radiative Liquid water
TYPE(FAD) :: YFAIRAD ! Radiative Ice
TYPE(FAD) :: YFAS    ! Snow
TYPE(FAD) :: YFAR    ! Rain
TYPE(FAD) :: YFAG    ! Graupels
TYPE(FAD) :: YFAH    ! Hail
TYPE(FAD) :: YFATKE  ! Turbulent Kinetic Energy
TYPE(FAD) :: YFAO3   ! Ozone mixing ratio
TYPE(FAD) :: YFACLF  ! Cloud fraction
TYPE(FAD) :: YFACPF  ! Convective precipitation flux
TYPE(FAD) :: YFASPF  ! Stratiform precipitation flux
TYPE(FAD) :: YFACVGQ ! CVGQ for French physics
TYPE(FAD) :: YFASDSAT ! Standard deviation of the Saturation Depression
TYPE(FAD) :: YFACVV  ! Convective Vertical Velocity
TYPE(FAD) :: YFARKTH ! Rasch-Kristjansson - tendency of enthalpy
TYPE(FAD) :: YFARKTQV ! Rasch-Kristjansson - tendency of water vapour
TYPE(FAD) :: YFARKTQC ! Rasch-Kristjansson - tendency of condensed water&ice

! Pronostic convection variables
TYPE(FAD) :: YFAUOM  ! Updraught vertic velocity
TYPE(FAD) :: YFAUAL  ! Updraught mesh fraction
TYPE(FAD) :: YFADOM  ! Downdraught vertic velocity
TYPE(FAD) :: YFADAL  ! Downdraught mesh fraction
TYPE(FAD) :: YFAUEN  ! Updraught entrainment
TYPE(FAD) :: YFAUNEBH! Pseudo Hist Conv cloud fraction

! EFB variables
TYPE(FAD) :: YFAEFB1 ! First variable for EFB turbulent paramterization
TYPE(FAD) :: YFAEFB2 ! Second variable for EFB turbulent paramterization
TYPE(FAD) :: YFAEFB3 ! Third variable for EFB turbulent paramterization

! TOUCANS variables
TYPE(FAD) :: YFATTE   ! Total turbulent energy   
TYPE(FAD) :: YFAMXL   ! Prognostic mixing length 
TYPE(FAD) :: YFASHTUR ! Shear source term for turbulence
TYPE(FAD) :: YFAFQTUR ! Flux form source term for turbulence -moisture
TYPE(FAD) :: YFAFSTUR ! Flux form source term for turbulence -enthalpy

! Filtered surface ln Ps for monitoring coupling updates:
TYPE(FAD) :: YFAFSP1  ,YFAFSP2  ,YFAFSP3  ,YFAFSP4  ,YFAFSP5
TYPE(FAD) :: YFASRC  ! Second order flux for Arome Turb. Scheme
TYPE(FAD) :: YFAQVA  ! Total humidity amplitude variation of Q+L+I

! Satellite images
TYPE(FAD) :: YFA_MSAT7C1  ! MSAT7   MVIRI  channel 1
TYPE(FAD) :: YFA_MSAT7C2  ! MSAT7   MVIRI  channel 2
TYPE(FAD) :: YFA_MSAT8C1  ! MSAT8   SEVIRI channel 1
TYPE(FAD) :: YFA_MSAT8C2  ! MSAT8   SEVIRI channel 2
TYPE(FAD) :: YFA_MSAT8C3  ! MSAT8   SEVIRI channel 3
TYPE(FAD) :: YFA_MSAT8C4  ! MSAT8   SEVIRI channel 4
TYPE(FAD) :: YFA_MSAT8C5  ! MSAT8   SEVIRI channel 5
TYPE(FAD) :: YFA_MSAT8C6  ! MSAT8   SEVIRI channel 6
TYPE(FAD) :: YFA_MSAT8C7  ! MSAT8   SEVIRI channel 7
TYPE(FAD) :: YFA_MSAT8C8  ! MSAT8   SEVIRI channel 8
TYPE(FAD) :: YFA_MSAT9C1  ! MSAT9   SEVIRI channel 1
TYPE(FAD) :: YFA_MSAT9C2  ! MSAT9   SEVIRI channel 2
TYPE(FAD) :: YFA_MSAT9C3  ! MSAT9   SEVIRI channel 3
TYPE(FAD) :: YFA_MSAT9C4  ! MSAT9   SEVIRI channel 4
TYPE(FAD) :: YFA_MSAT9C5  ! MSAT9   SEVIRI channel 5
TYPE(FAD) :: YFA_MSAT9C6  ! MSAT9   SEVIRI channel 6
TYPE(FAD) :: YFA_MSAT9C7  ! MSAT9   SEVIRI channel 7
TYPE(FAD) :: YFA_MSAT9C8  ! MSAT9   SEVIRI channel 8
TYPE(FAD) :: YFA_GOES11C1 ! GOES11  IMAGER channel 1
TYPE(FAD) :: YFA_GOES11C2 ! GOES11  IMAGER channel 2
TYPE(FAD) :: YFA_GOES11C3 ! GOES11  IMAGER channel 3
TYPE(FAD) :: YFA_GOES11C4 ! GOES11  IMAGER channel 4
TYPE(FAD) :: YFA_GOES12C1 ! GOES12  IMAGER channel 1
TYPE(FAD) :: YFA_GOES12C2 ! GOES12  IMAGER channel 2
TYPE(FAD) :: YFA_GOES12C3 ! GOES12  IMAGER channel 3
TYPE(FAD) :: YFA_GOES12C4 ! GOES12  IMAGER channel 4
TYPE(FAD) :: YFA_MTSAT1C1 ! MTSAT 1 IMAGER channel 1
TYPE(FAD) :: YFA_MTSAT1C2 ! MTSAT 1 IMAGER channel 2
TYPE(FAD) :: YFA_MTSAT1C3 ! MTSAT 1 IMAGER channel 3
TYPE(FAD) :: YFA_MTSAT1C4 ! MTSAT 1 IMAGER channel 4

! Instantaneous fluxes
TYPE(FAD) :: YFA_XTCC    ! SURFNEBUL.TOTALE 
TYPE(FAD) :: YFA_X10U    ! CLSVENT.ZONAL 
TYPE(FAD) :: YFA_X10V    ! CLSVENT.MERIDIEN 
TYPE(FAD) :: YFA_X10NU   ! CLSVENTNEUTRE.U
TYPE(FAD) :: YFA_X10NV   ! CLSVENTNEUTRE.V
TYPE(FAD) :: YFA_X2T     ! CLSTEMPERATURE 
TYPE(FAD) :: YFA_X2TPW   ! CLSTPRIM_W
TYPE(FAD) :: YFA_XMRT    ! CLSMEAN.RAD.TEMP
TYPE(FAD) :: YFA_X2SH    ! CLSHUMI.SPECIFIQ 
TYPE(FAD) :: YFA_X2RH    ! CLSHUMI.RELATIVE 
TYPE(FAD) :: YFA_XCCC    ! SURFNEBUL.CONVEC 
TYPE(FAD) :: YFA_XHCC    ! SURFNEBUL.HAUTE 
TYPE(FAD) :: YFA_XMCC    ! SURFNEBUL.MOYENN 
TYPE(FAD) :: YFA_XLCC    ! SURFNEBUL.BASSE 
TYPE(FAD) :: YFA_XX2T    ! CLSMAXI.TEMPERAT 
TYPE(FAD) :: YFA_XN2T    ! CLSMINI.TEMPERAT 
TYPE(FAD) :: YFA_XCVU    ! S000FL.U CONVEC
TYPE(FAD) :: YFA_XCVV    ! S000FL.V CONVEC
TYPE(FAD) :: YFA_XCVQ    ! S000FL.Q CONVEC
TYPE(FAD) :: YFA_XCVS    ! S000FL.CT CONVEC
TYPE(FAD) :: YFA_XTUU    ! S000FL.U TURBUL
TYPE(FAD) :: YFA_XTUV    ! S000FL.V TURBUL
TYPE(FAD) :: YFA_XTUQ    ! S000FL.Q TURBUL
TYPE(FAD) :: YFA_XTUS    ! S000FL.CT TURBUL
TYPE(FAD) :: YFA_XGDU    ! S000FL.U ONDG.OR
TYPE(FAD) :: YFA_XGDV    ! S000FL.V ONDG.OR
TYPE(FAD) :: YFA_XLSP    ! S000PLUIE STRATI
TYPE(FAD) :: YFA_XCP     ! SURFPREC.EAU.CON 
TYPE(FAD) :: YFA_XLSS    ! S000NEIGE STRATI
TYPE(FAD) :: YFA_XCSF    ! S000NEIGE CONVEC
TYPE(FAD) :: YFA_XLSG    ! S000GRAUP STRATI
TYPE(FAD) :: YFA_XCSG    ! S000GRAUP CONVEC
TYPE(FAD) :: YFA_XLSH    ! S000GRELE STRATI
TYPE(FAD) :: YFA_XCSH    ! S000GRELE CONVEC
TYPE(FAD) :: YFA_XSSR    ! S000RAYT.SOLAIRE 
TYPE(FAD) :: YFA_XSTR    ! S000RAYT.TERREST 
TYPE(FAD) :: YFA_XTSR    ! SOMMRAYT.SOLAIRE 
TYPE(FAD) :: YFA_XTTR    ! SOMMRAYT.TERREST 
TYPE(FAD) :: YFA_X10FF   ! CLSWIND.VELOCITY 
TYPE(FAD) :: YFA_XCAPE   ! SURFCAPE.MOD.XFU 
TYPE(FAD) :: YFA_XCTOP   ! CLPTOPCONV.MOD.XFU
TYPE(FAD) :: YFA_XMOCO   ! CLPMOCON.MOD.XFU 
TYPE(FAD) :: YFA_XCLPH   ! CLPMHAUT.MOD.XFU 
TYPE(FAD) :: YFA_XVEIN   ! CLPVEIND.MOD.XFU 
TYPE(FAD) :: YFA_XGUST   ! CLSRAFAL.MOD.XFU 
TYPE(FAD) :: YFA_XUGST   ! CLSU.RAF.MOD.XFU 
TYPE(FAD) :: YFA_XVGST   ! CLSV.RAF.MOD.XFU 
TYPE(FAD) :: YFA_XUGST2  ! CLSU.RAF2.MOD.XFU 
TYPE(FAD) :: YFA_XVGST2  ! CLSV.RAF2.MOD.XFU 
TYPE(FAD) :: YFA_XX2HU   ! CLSMAXI.HUMI.REL 
TYPE(FAD) :: YFA_XN2HU   ! CLSMINI.HUMI.REL 
TYPE(FAD) :: YFA_XTHW    ! THETAPWP_FLUX 
TYPE(FAD) :: YFA_XXDIAGH ! SURFDIAGHAIL
TYPE(FAD) :: YFA_VISICLD ! VISICLD
TYPE(FAD) :: YFA_VISIHYD ! VISIHYDRO
TYPE(FAD) :: YFA_MXCLWC  ! MAX CLWC
TYPE(FAD) :: YFA_VISICLD2! VISICLD
TYPE(FAD) :: YFA_VISIHYD2! VISIHYDRO
TYPE(FAD) :: YFA_MXCLWC2 ! MAX CLWC
TYPE(FAD) :: YFA_XPTYPE  ! PTYPE FREQUENT
TYPE(FAD) :: YFA_XPTYPESEV! PTYPE SEVERE
TYPE(FAD) :: YFA_XPTYPE2 ! PTYPE FREQUENT
TYPE(FAD) :: YFA_XPTYPESEV2! PTYPE SEVERE


! Cumulated fluxes
TYPE(FAD) :: YFA_CLSP  ! SURFPREC.EAU.GEC
TYPE(FAD) :: YFA_CCP   ! SURFPREC.EAU.CON
TYPE(FAD) :: YFA_CLSS  ! SURFPREC.NEI.GEC
TYPE(FAD) :: YFA_CCSF  ! SURFPREC.NEI.CON
TYPE(FAD) :: YFA_CLSG  ! SURFPREC.GRA.GEC
TYPE(FAD) :: YFA_CCSG  ! SURFPREC.GRA.CON
TYPE(FAD) :: YFA_CLSH  ! SURFPREC.GRE.GEC
TYPE(FAD) :: YFA_CCSH  ! SURFPREC.GRE.CON
TYPE(FAD) :: YFA_CUSS  ! SURFTENS.TURB.ZO
TYPE(FAD) :: YFA_CVSS  ! SURFTENS.TURB.ME
TYPE(FAD) :: YFA_CSSH  ! SURFFLU.CHA.SENS
TYPE(FAD) :: YFA_CSLH  ! SURFCHAL LATENTE
TYPE(FAD) :: YFA_CTSP  ! SURFPRESSION SOL
TYPE(FAD) :: YFA_CTCC  ! ATMONEBUL.TOTALE
TYPE(FAD) :: YFA_CBLD  ! SURFDISSIP SURF 
TYPE(FAD) :: YFA_CSSR  ! SURFFLU.RAY.SOLA
TYPE(FAD) :: YFA_CSTR  ! SURFFLU.RAY.THER
TYPE(FAD) :: YFA_CTSR  ! SOMMFLU.RAY.SOLA
TYPE(FAD) :: YFA_CTTR  ! SOMMFLU.RAY.THER
TYPE(FAD) :: YFA_CCCC  ! ATMONEBUL.CONVEC
TYPE(FAD) :: YFA_CHCC  ! ATMONEBUL.HAUTE 
TYPE(FAD) :: YFA_CMCC  ! ATMONEBUL.MOYENN
TYPE(FAD) :: YFA_CLCC  ! ATMONEBUL.BASSE 
TYPE(FAD) :: YFA_CUGW  ! SURFTENS.DMOG.ZO
TYPE(FAD) :: YFA_CVGW  ! SURFTENS.DMOG.ME
TYPE(FAD) :: YFA_CUTO  ! SURFTENS.TOTA.ZO
TYPE(FAD) :: YFA_CVTO  ! SURFTENS.TOTA.ME
TYPE(FAD) :: YFA_CE    ! SURFFLU.MEVAP.EA
TYPE(FAD) :: YFA_CS    ! SURFFLU.MSUBL.NE
TYPE(FAD) :: YFA_CT    ! SURFFLU.MTOTA.NE
TYPE(FAD) :: YFA_CLHE  ! SURFFLU.LAT.MEVA
TYPE(FAD) :: YFA_CLHS  ! SURFFLU.LAT.MSUB
TYPE(FAD) :: YFA_CLHT  ! SURFFLU.LAT.MTOT
TYPE(FAD) :: YFA_CWS   ! SURFCONTENU EAU 
TYPE(FAD) :: YFA_CSNS  ! SURFRESERV NEIGE
TYPE(FAD) :: YFA_CQTO  ! ATMOHUMI TOTALE 
TYPE(FAD) :: YFA_CTO3  ! ATMOOZONE TOTALE
TYPE(FAD) :: YFA_CTME  ! TOPMESO ENTH    
TYPE(FAD) :: YFA_CICE  ! ATMOHUMI SOLIDE 
TYPE(FAD) :: YFA_CLI   ! ATMOHUMI LIQUIDE
TYPE(FAD) :: YFA_CCVU  ! SURFTENS.CONV.ZO
TYPE(FAD) :: YFA_CCVV  ! SURFTENS.CONV.ME
TYPE(FAD) :: YFA_CCVQ  ! SURFCFU.Q.CONVEC
TYPE(FAD) :: YFA_CCVS  ! SURFCFU.CT.CONVE
TYPE(FAD) :: YFA_CTUQ  ! SURFCFU.Q.TURBUL
TYPE(FAD) :: YFA_CTUS  ! SURFCFU.CT.TURBU
TYPE(FAD) :: YFA_CSOC  ! SURFRAYT SOL CL 
TYPE(FAD) :: YFA_CTHC  ! SURFRAYT THER CL
TYPE(FAD) :: YFA_CSOP  ! SURFRAYT DIR SUR
TYPE(FAD) :: YFA_CDNI  ! SURFDIR NORM IRR
TYPE(FAD) :: YFA_CGNI  ! SURFGLB NORM IRR
TYPE(FAD) :: YFA_CTOP  ! TOPRAYT DIR SOM 
TYPE(FAD) :: YFA_CSOD  ! SURFRAYT SOLA DE
TYPE(FAD) :: YFA_CTHD  ! SURFRAYT THER DE
TYPE(FAD) :: YFA_CFON  ! SURFFONTE NEIGE 
TYPE(FAD) :: YFA_CCHS  ! SURFCHAL. DS SOL
TYPE(FAD) :: YFA_CEAS  ! SURFEAU DANS SOL
TYPE(FAD) :: YFA_CSRU  ! SURFRUISSELLEMEN
TYPE(FAD) :: YFA_CDRU  ! PROFRUISSELLEMEN
TYPE(FAD) :: YFA_CIRU  ! SURFRUISS. INTER
TYPE(FAD) :: YFA_CETP  ! SURFEVAPOTRANSPI
TYPE(FAD) :: YFA_CTP   ! SURFTRANSPIRATIO
TYPE(FAD) :: YFA_CSMR  ! SURFRAYT.LUNE.DE
TYPE(FAD) :: YFA_CTSOC ! SOMMRAYT SOL CL 
TYPE(FAD) :: YFA_CTTHC ! SOMMRAYT THER CL
TYPE(FAD) :: YFA_CDUTP ! SURFTIME.PREC.TO
TYPE(FAD) :: YFA_CFLASH! Lightning density

TYPE(FAD) :: YFA_CAPE    ! Convective Available Potential Energy
TYPE(FAD) :: YFA_CIEN    ! Convective Inhibition
TYPE(FAD) :: YFA_TCVS    ! Temperature of Convection
TYPE(FAD) :: YFA_ASMC    ! Analysed Surface Moisture Convergence
TYPE(FAD) :: YFA_VSMC    ! Surface Moisture Convergence for VARPACK
TYPE(FAD) :: YFA_ACCGREL ! Dummy hail field

TYPE(FAD) :: YFA_SRH     ! Storm relative helicity
TYPE(FAD) :: YFA_STRMMU  ! U component of storm motion
TYPE(FAD) :: YFA_STRMMV  ! V component of storm motion

!     ------------------------------------------------------------------

INTEGER (KIND=JPIM) :: JPXTRO, JPXLAT, JPXNIV, JPNXFA, JPNXCA

END MODULE YOMFA
