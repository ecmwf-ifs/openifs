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

MODULE YOEAERATM

USE PARKIND1,     ONLY : JPIM, JPRB
USE YOE_AERODIAG, ONLY : NPAERODIAG
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERATM* - CONTROL PARAMETERS FOR AEROSOLS IN THE ATMOSPHERE
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NMAXTAER = 10

TYPE :: TYPE_AERO_DESC
SEQUENCE
! Automatically initialised keys to match GFL
CHARACTER(LEN=16)  :: CNAME                 ! ARPEGE field name
INTEGER(KIND=JPIM) :: IGRBCODE              ! GRIB code
! Initialised from namelist
INTEGER(KIND=JPIM) :: NTYP                  ! Aerosol type number
INTEGER(KIND=JPIM) :: NBIN                  ! Aerosol bin number within type
INTEGER(KIND=JPIM) :: IGRIBDIAG(NPAERODIAG) ! GRIB codes of diagnostics
REAL(KIND=JPRB)    :: RDDEPVSEA             ! Dry deposition velocity over sea
REAL(KIND=JPRB)    :: RDDEPVLIC             ! Dry deposition velocity over land or ice
REAL(KIND=JPRB)    :: RSEDIMV               ! Sedimentation velocity
REAL(KIND=JPRB)    :: RSCAVIN               ! In-cloud scavenging fraction
REAL(KIND=JPRB)    :: RSCAVBCR              ! Below-cloud scavenging coefficient for rain
REAL(KIND=JPRB)    :: RSCAVBCS              ! Below-cloud scavenging coefficitne for snow
CHARACTER(LEN=16)  :: COPTCLASS             ! Optical property class
CHARACTER(LEN=16)  :: CHYGCLASS             ! Hygroscopicity class
INTEGER(KIND=JPIM) :: IAEROCV               ! Aerosol control variable (0=none, 1=fine, 2=coarse)
END TYPE TYPE_AERO_DESC

TYPE :: TEAERATM
INTEGER(KIND=JPIM) :: NAERCONF
INTEGER(KIND=JPIM) :: NINIDAY   
INTEGER(KIND=JPIM) :: NXT3DAER
INTEGER(KIND=JPIM) :: NDD1, NSS1
INTEGER(KIND=JPIM) :: NBCOPTP, NDDOPTP, NOMOPTP, NSSOPTP, NSUOPTP
INTEGER(KIND=JPIM) :: NVISWL
INTEGER(KIND=JPIM) :: NTAER
INTEGER(KIND=JPIM) :: NTYPAER(NMAXTAER)
INTEGER(KIND=JPIM) :: NAER_BLNUCL
INTEGER(KIND=JPIM) :: NAERSCAV

REAL(KIND=JPRB) :: RGRATE

REAL(KIND=JPRB) :: REPSCAER

LOGICAL :: LAERCLIMG, LAERCLIMZ, LAERCLIST, LAERDRYDP, LAERHYGRO, LAERLISI 
LOGICAL :: LAERNGAT , LAERSEDIM, LAERSURF , LAERELVS , LAER6SDIA,LAERSEDIMSS
LOGICAL :: LAERGTOP , LAERRAD  , LAERCCN  , LAEROPT(9),LAERINIT , LAERVOL
LOGICAL :: LAERCSTR , LAERDIAG1, LAERDIAG2, LAERRRTM , LAERUVP
LOGICAL :: LAEREXTR , LAERGBUD , LAERPRNT
LOGICAL :: LAERNITRATE
LOGICAL :: LSEASALT_RH80
LOGICAL :: LAERDUSTSOURCE
LOGICAL :: LAERDUSTSIZEVAR
LOGICAL :: LAERDUST_NEWBIN

! Pre-computed conversion factors according to LSEASALT_RH80.
! Computed in SU_AERW after reading namelist and calling SU_AERP
! to initialise growth tables.
REAL(KIND=JPRB) :: RSS_DRY_DIAFAC, RSS_DRY_DENSFAC, RSS_DRY_MASSFAC
REAL(KIND=JPRB) :: RSS_RH80_DIAFAC, RSS_RH80_DENSFAC, RSS_RH80_MASSFAC
REAL(KIND=JPRB) :: RAERDUST_REBOUND

TYPE(TYPE_AERO_DESC), POINTER :: YAERO_DESC(:) => NULL()
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEAERATM
!============================================================================

TYPE(TEAERATM), POINTER :: YREAERATM => NULL()

!     ------------------------------------------------------------------
! NDD1       : location of first bin for desert dust
! NSS1       : location of first bin for sea-salt
! NBCOPTP    : index for choosing the black carbon LW and SW optic.prop. (1: Boucher, 2: Bond, Bergstrom, 3: Stier et al.)
! NDDOPTP    : index for choosing the dust LW and SW optic.prop. (1: Dubovik(SW)+Fouquart(LW), 2: Fouquart, 3: Woodward)
! NOMOPTP    : index for choosing the organic carbon optic.prop.
! NSSOPTP    : index for choosing the sea salt optic.prop.
! NSUOPTP    : index for choosing the sulphate optic.prop.
! NVISWL     : index of wavelength for visibility computations
! RMFMIN     : minimum mass flux for convective aerosol transport
! RGRATE     : transformation rate from hygrophopic to hygrophilic for BC and OM aerosols

! REPSCAER   : security on aerosol concentration: always >= 1.E-18

! LAERCLIMG  : .T. to start prognostic aerosols with geographical monthly 
!                  mean climatology
! LAERCLIMZ  : .T. to start prognostic aerosols with zonal annual mean 
!                  climatology
! LAERCLIST  : .T. to start prognostic aerosols with geographical monthly 
!                  mean climatology for background stratospheric only
! LAERDRYDP  : .T. dry deposition is active
! LAERHYDRO  : .T. hygroscopic effects on BC and OM aerosols
! LAERNGAT   : .T. prevents negative aerosol concentrations
! NAERSCAV   : aerosol scanvenging scheme: 1=historical, 2=from CB05, 3=Luo et al. 2019
! LAERSEDIM  : .T. sedimentation is active
! LAERSEDIMSS  : .T. special sedimentation for sea-salt is active
! LAERSURF   : .T. if surface emissions
! LAERELVS   : .T. if "elevated" source
! LAER6SDIA  : .T. if radiance diagnostics with 6S
! LAERGTOP   : .T. if gas-to-particle conversion for SO2/SO4
! LAERRAD    : .T. if there is any prognostic aerosols used for RT
! LAEROPT(.) : .T. if a given aerosol type is radiatively interactive
! LAERCCN    : .T. if prognostic aerosols are used to define the Re of liq.wat.clds
! LAERUVP    : .T. if prognostic aerosols are used in UV-processor
! LAERCSTR   : .T. if climatological stratospheric aerosols are used in radiation 
!                  schemes with the prognostic tropospheric aerosols.
! LAERRRTM   : .T. if RRTM schemes get the information from the prognostic aerosols
! LAERINIT   : .T. if analysed prognostic aerosols are ONLY used as "climatological" aerosols
! LAERVOL    : .T. if volcanic aerosol is considered

! NMAXTAER   : MAXIMUM TOTAL NUMBER OF AEROSOLS
! NTAER      : TOTAL NUMBER OF AEROSOLS
! NTYPAER( ) : NBINAER
!        (1) :         3/9 FOR SEA-SALT 
!        (2) :         3/9 FOR DESERT DUST
!        (3) :         2 FOR ORGANIC MATTERS
!        (4) :         2 FOR BLACK CARBON
!        (5) :         2 FOR SO2/SO4 (ESSENTIALLY TROPOSPHERIC)
!        (6) :         2 FOR NITRATES
!        (7) :         1 FOR AMMONIUM
!        (8) :         1 FOR FLY ASH
!        (9) :         2 FOR SO2/SO4 OF VOLCANIC ORIGIN

! LSEASALT_RH80: transport sea salt at fixed 80% RH (as done historically) rather than dry (as all other aerosols)
! YAERO_DESC: aerosol descriptors (metadata moved out of GFL and hard-coded source)
!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TEAERATM), INTENT(IN) :: SELF
INTEGER        , INTENT(IN) :: KDEPTH
INTEGER        , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOEAERATM:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yreaeratm : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NAERCONF = ', SELF%NAERCONF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NINIDAY = ', SELF%NINIDAY
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NXT3DAER = ', SELF%NXT3DAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDD1 = ', SELF%NDD1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSS1 = ', SELF%NSS1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NBCOPTP = ', SELF%NBCOPTP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDDOPTP = ', SELF%NDDOPTP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NOMOPTP = ', SELF%NOMOPTP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSSOPTP = ', SELF%NSSOPTP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSUOPTP = ', SELF%NSUOPTP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVISWL = ', SELF%NVISWL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMAXTAER = ', NMAXTAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTAER = ', SELF%NTAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTYPAER SUM = ', SUM(SELF%NTYPAER)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RGRATE = ', SELF%RGRATE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REPSCAER = ', SELF%REPSCAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERCLIMG = ', SELF%LAERCLIMG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERCLIMZ = ', SELF%LAERCLIMZ
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERCLIST = ', SELF%LAERCLIST
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERDRYDP = ', SELF%LAERDRYDP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERHYGRO = ', SELF%LAERHYGRO
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERLISI = ', SELF%LAERLISI
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERNGAT = ', SELF%LAERNGAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NAERSCAV = ', SELF%NAERSCAV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERSEDIM = ', SELF%LAERSEDIM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERSEDIMSS = ', SELF%LAERSEDIMSS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERSURF = ', SELF%LAERSURF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERELVS = ', SELF%LAERELVS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAER6SDIA = ', SELF%LAER6SDIA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERGTOP = ', SELF%LAERGTOP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERRAD = ', SELF%LAERRAD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERCCN = ', SELF%LAERCCN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAEROPT values ', SELF%LAEROPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERINIT = ', SELF%LAERINIT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERVOL = ', SELF%LAERVOL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERCSTR = ', SELF%LAERCSTR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERDIAG1 = ', SELF%LAERDIAG1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERDIAG2 = ', SELF%LAERDIAG2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERRRTM = ', SELF%LAERRRTM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERUVP = ', SELF%LAERUVP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAEREXTR = ', SELF%LAEREXTR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERGBUD = ', SELF%LAERGBUD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERPRNT = ', SELF%LAERPRNT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERNITRATE = ', SELF%LAERNITRATE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSEASALT_RH80 = ', SELF%LSEASALT_RH80
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERDUSTSOURCE = ', SELF%LAERDUSTSOURCE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERDUSTSIZEVAR = ', SELF%LAERDUSTSIZEVAR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERDUST_NEWBIN = ', SELF%LAERDUST_NEWBIN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAERDUST_REBOUND = ', SELF%RAERDUST_REBOUND
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSS_DRY_DIAFAC = ', SELF%RSS_DRY_DIAFAC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSS_DRY_DENSFAC = ', SELF%RSS_DRY_DENSFAC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSS_DRY_MASSFAC = ', SELF%RSS_DRY_MASSFAC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSS_RH80_DIAFAC = ', SELF%RSS_RH80_DIAFAC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSS_RH80_DENSFAC = ', SELF%RSS_RH80_DENSFAC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSS_RH80_MASSFAC = ', SELF%RSS_RH80_MASSFAC
IF (LHOOK) CALL DR_HOOK('YOEAERATM:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEAERATM

