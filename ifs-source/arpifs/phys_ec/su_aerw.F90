! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU_AERW(YDMODEL)

!**** *SU_AERW*   - DEFINES INDICES AND PARAMETERS FOR VARIOUS AEROSOL VARIABLES

!     PURPOSE.
!     --------
!           INITIALIZE YOEAERATM, YOEAERSRC, YOEAERSNK, THE MODULES THAT CONTAINS INDICES
!           ALLOWING TO GET THE AEROSOL PARAMETERS RELEVANT FOR THE PROGNOSTIC AEROSOL
!           CONFIGURATION.

!**   INTERFACE.
!     ----------
!        *CALL* *SU_AERW

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        YOEAERATM, YOEAERSRC

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE *ECMWF*
!        ORIGINAL : 2005-07-08

!     MODIFICATIONS.
!     --------------
!     G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!     H. Hersbach  : 01-04-2011 Initialize extinction coeffs for historical aerosols
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     K. Yessad (July 2014): Move some variables.
!     S.Remy (Jan 2015): Change the dust emissions, RDDUSRC, more larger particles
!     S.Remy (Jan 2016): Reduce dust emissions over Taklamakan and India
!     S.Remy (Jan 2016): LAERSOA_CHEM added
!     S.Remy (Sep 2016): nitrates added
!     S.Remy (Jul 2017): new dust scheme: LAERDUSTSOURCE and LAERDUST_NEWBIN
!     added
!     ------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE PARKIND1   , ONLY : JPIM, JPRB
USE YOMHOOK    , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN     , ONLY : NULNAM, NULOUT
USE YOMRIP0    , ONLY : NINDAT, NSSSSS
USE YOEAERATM  , ONLY : TYPE_AERO_DESC, NMAXTAER





!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)  ,INTENT(INOUT), TARGET :: YDMODEL
INTEGER(KIND=JPIM) :: IAER, ICAER, ITAER
INTEGER(KIND=JPIM) :: J, JAER, JVOLE

!-- map
REAL(KIND=JPRB) :: ZDDUAER(50)
CHARACTER(LEN=45) :: CLAERWND(0:3)
CHARACTER(LEN=100) :: CFMT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ----------------------------------------------------------------

#include "posnam.intfb.h"

#include "su_aerp.intfb.h"
#include "su_aerop.intfb.h"
#include "abor1.intfb.h"





!      ----------------------------------------------------------------

LOGICAL, POINTER :: LAERCLIMG, LAERCLIMZ, LAERCLIST, LAERDRYDP, LAERELVS,&
 & LAEREXTR, LAERGBUD, LAERHYGRO, LAERNGAT, LAERPRNT, &
 & LAERSEDIM, LAERSURF, LAERGTOP, LAER6SDIA, LAERCCN, LAERSEDIMSS, &
 & LAERINIT, LAERCSTR, LAERRRTM, LAERDIAG1, LAERDIAG2, LAERUVP,&
 & LAERVOL, LEPAERO, LAEROMIN, LOCNDMS,&
 & LSEASALT_RH80, LAERDUSTSOURCE, LAERDUST_NEWBIN, LAERDUSTSIZEVAR

INTEGER(KIND=JPIM), POINTER :: NAERCONF, NXT3DAER, NINIDAY, NBCOPTP,&
 & NDDOPTP, NOMOPTP, NSSOPTP, NSUOPTP, NVISWL, NDRYDEP,  NDRYDEPVEL_DYN, &
 & NTYPAER(:), NDDUST, NSO4SCHEME, NSSALT, NDMSO, NPIST, KSTPDBG(:), NSTPDBG,&
 & NAERWND, NAERLISI, NDUSRCP(:), NAERVOLC, NAERVOLE, NVOLDATS(:), NVOLERUP,&
 & NVOLERUZ(:), NVOLOPTP, NVOLHOMO, NINTERPT, NAER_BLNUCL, NVOLDATE(:), NAERSCAV
REAL(KIND=JPRB), POINTER :: RAERDUB, RDDUAER(:), RFCTDUR, RFCTSSR, RLATVOL,&
 & RLONVOL, RSUCV1, RSUCV2, RAERDUST_REBOUND, &
 & RAERVOLC(:,:), RAERVOLE(:,:), RVOLERUZ(:)

TYPE(TYPE_AERO_DESC), POINTER :: YAERO_DESC(:)

#include "naeaer.nam.h"
#include "naevol.nam.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU_AERW',0,ZHOOK_HANDLE)
ASSOCIATE(YDEAERMAP=>YDMODEL%YRML_PHY_AER%YREAERMAP,YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
 & YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, YDEAERVOL=>YDMODEL%YRML_PHY_AER%YREAERVOL)

ASSOCIATE(REPSCAER=>YDEAERATM%REPSCAER, AERO_SCHEME=>YDCOMPO%AERO_SCHEME, &
 & RSS_DRY_DIAFAC=>YDEAERATM%RSS_DRY_DIAFAC, &
 & RSS_DRY_DENSFAC=>YDEAERATM%RSS_DRY_DENSFAC, &
 & RSS_DRY_MASSFAC=>YDEAERATM%RSS_DRY_MASSFAC, &
 & RSS_RH80_DIAFAC=>YDEAERATM%RSS_RH80_DIAFAC, &
 & RSS_RH80_DENSFAC=>YDEAERATM%RSS_RH80_DENSFAC, &
 & RSS_RH80_MASSFAC=>YDEAERATM%RSS_RH80_MASSFAC, &
 & RSSGROWTH_RHTAB=>YDEAERSNK%RSSGROWTH_RHTAB, &
 & RSSDENS_RHTAB=>YDEAERSNK%RSSDENS_RHTAB, &
 & NACTAERO=>YGFL%NACTAERO, NAERO=>YGFL%NAERO, &
 & YAERO_NL=>YGFL%YAERO_NL, &
 & LAERRAD=>YDEAERATM%LAERRAD, &
 & RDDUSRC=>YDEAERSRC%RDDUSRC, RFCTSS=>YDEAERSRC%RFCTSS, &
 & RSSFLX=>YDEAERSRC%RSSFLX, RDMSMIN=>YDEAERSRC%RDMSMIN, &
 & RDUSRCP=>YDEAERMAP%RDUSRCP, RFCTDU=>YDEAERSRC%RFCTDU, &
 & LESO4HIS=>YDERAD%LESO4HIS, NCHEM=>YGFL%NCHEM, &
 & LAERCHEM=>YGFL%LAERCHEM, &
 & LAERNITRATE=>YDCOMPO%LAERNITRATE,&
 & LAERSOA=>YDCOMPO%LAERSOA, &
 & LAERSOA_COUPLED=>YDCOMPO%LAERSOA_COUPLED)

ALLOCATE(YDEAERATM%YAERO_DESC(YGFL%NACTAERO))

! Associate pointers for variables in namelist
LAERCLIMG => YDEAERATM%LAERCLIMG
NDRYDEPVEL_DYN=> YDEAERSNK%NDRYDEPVEL_DYN
LAERCLIMZ => YDEAERATM%LAERCLIMZ
LAERCLIST => YDEAERATM%LAERCLIST
LAERDRYDP => YDEAERATM%LAERDRYDP
LAERELVS  => YDEAERATM%LAERELVS
LAEREXTR  => YDEAERATM%LAEREXTR
LAERGBUD  => YDEAERATM%LAERGBUD
LAERDUSTSOURCE => YDEAERATM%LAERDUSTSOURCE
RAERDUST_REBOUND => YDEAERATM%RAERDUST_REBOUND
LAERDUST_NEWBIN=>YDEAERATM%LAERDUST_NEWBIN
LAERDUSTSIZEVAR=>YDEAERATM%LAERDUSTSIZEVAR
LAERHYGRO => YDEAERATM%LAERHYGRO
LAERNGAT  => YDEAERATM%LAERNGAT
LAERPRNT  => YDEAERATM%LAERPRNT
LAERSEDIM => YDEAERATM%LAERSEDIM
LAERSEDIMSS => YDEAERATM%LAERSEDIMSS
LAERSURF  => YDEAERATM%LAERSURF
LAERGTOP  => YDEAERATM%LAERGTOP
LAER6SDIA => YDEAERATM%LAER6SDIA
LAERCCN   => YDEAERATM%LAERCCN
LAERINIT  => YDEAERATM%LAERINIT
LAERCSTR  => YDEAERATM%LAERCSTR
LAERRRTM  => YDEAERATM%LAERRRTM
LAERDIAG1 => YDEAERATM%LAERDIAG1
LAERDIAG2 => YDEAERATM%LAERDIAG2
LAERUVP   => YDEAERATM%LAERUVP
LAERVOL   => YDEAERATM%LAERVOL
NAERCONF  => YDEAERATM%NAERCONF
NXT3DAER  => YDEAERATM%NXT3DAER
NINIDAY   => YDEAERATM%NINIDAY
NBCOPTP   => YDEAERATM%NBCOPTP
NDDOPTP   => YDEAERATM%NDDOPTP
NOMOPTP   => YDEAERATM%NOMOPTP
NSSOPTP   => YDEAERATM%NSSOPTP
NSUOPTP   => YDEAERATM%NSUOPTP
NVISWL    => YDEAERATM%NVISWL
NAERSCAV => YDEAERATM%NAERSCAV
LSEASALT_RH80=> YDEAERATM%LSEASALT_RH80
NTYPAER   => YDEAERATM%NTYPAER
YAERO_DESC=> YDEAERATM%YAERO_DESC
LEPAERO   => YDEAERSRC%LEPAERO
LAEROMIN  => YDEAERSRC%LAEROMIN
LOCNDMS   => YDEAERSRC%LOCNDMS
NDDUST    => YDEAERSRC%NDDUST
NSO4SCHEME    => YDEAERSRC%NSO4SCHEME
NAER_BLNUCL    => YDEAERATM%NAER_BLNUCL
NSSALT    => YDEAERSRC%NSSALT
NDMSO     => YDEAERSRC%NDMSO
NPIST     => YDEAERSRC%NPIST
KSTPDBG   => YDMODEL%YRML_PHY_AER%YREDBUG%KSTPDBG
NSTPDBG   => YDMODEL%YRML_PHY_AER%YREDBUG%NSTPDBG
NAERWND   => YDEAERSRC%NAERWND
NAERLISI  => YDMODEL%YRML_PHY_AER%YREAERLID%NAERLISI
NDUSRCP   => YDEAERMAP%NDUSRCP
RAERDUB   => YDEAERSRC%RAERDUB
RDDUAER   => YDEAERMAP%RDDUAER
RFCTDUR   => YDEAERSRC%RFCTDUR
RFCTSSR   => YDEAERSRC%RFCTSSR
RLATVOL   => YDEAERSRC%RLATVOL
RLONVOL   => YDEAERSRC%RLONVOL
NDRYDEP   => YDEAERSNK%NDRYDEP
RSUCV1    => YDEAERSNK%RSUCV1
RSUCV2    => YDEAERSNK%RSUCV2
NAERVOLC  => YDEAERVOL%NAERVOLC
NAERVOLE  => YDEAERVOL%NAERVOLE
NVOLDATS  => YDEAERVOL%NVOLDATS
NVOLDATE  => YDEAERVOL%NVOLDATE
NVOLERUP  => YDEAERVOL%NVOLERUP
NVOLERUZ  => YDEAERVOL%NVOLERUZ
NVOLOPTP  => YDEAERVOL%NVOLOPTP
NVOLHOMO  => YDEAERVOL%NVOLHOMO
NINTERPT  => YDEAERVOL%NINTERPT
RAERVOLC  => YDEAERVOL%RAERVOLC
RAERVOLE  => YDEAERVOL%RAERVOLE
RVOLERUZ  => YDEAERVOL%RVOLERUZ

!      ----------------------------------------------------------------

!*       1.       DEFAULT VALUES OF PARAMETERS
!                 ----------------------------

DO JAER=1,NACTAERO
  ! write (*,*) "SUAERW, JAER=",JAER,YAERO_NL(JAER)%CNAME
  ! Copy these from YAERO_NL, not YAERO itself which won't have
  ! been initialised yet if we're being called early from SUGFL1
  YAERO_DESC(JAER)%CNAME = YAERO_NL(JAER)%CNAME
  YAERO_DESC(JAER)%IGRBCODE = YAERO_NL(JAER)%IGRBCODE

  ! These will be read in via the NAEAER namelist
  YAERO_DESC(JAER)%IGRIBDIAG(:) = -9999
  YAERO_DESC(JAER)%RDDEPVSEA = 0._JPRB
  YAERO_DESC(JAER)%RDDEPVLIC = 0._JPRB
  YAERO_DESC(JAER)%RSEDIMV = 0._JPRB
  YAERO_DESC(JAER)%RSCAVIN = 0._JPRB
  YAERO_DESC(JAER)%RSCAVBCR = 0._JPRB
  YAERO_DESC(JAER)%RSCAVBCS = 0._JPRB
  YAERO_DESC(JAER)%COPTCLASS = ''
  YAERO_DESC(JAER)%CHYGCLASS = ''
  YAERO_DESC(JAER)%IAEROCV = 0

  ! These will be initialised based on NTYPAER if not given
  ! in the namelist
  YAERO_DESC(JAER)%NTYP = -9999
  YAERO_DESC(JAER)%NBIN = -9999
ENDDO

! the 10 types and assumed number of bins are:
!  NTYPAER    bins  type
!     1       1- 3  sea-salt  0.03 - 0.5 -  5  - 20 microns
!     2       4- 6  dust      0.03 - 0.5 - 0.9 - 20 microns
!     3       7- 8  POM    hydrophilic, hydrophobic
!     4       9-10  BC     hydrophilic, hydrophobic
!     5      11-12  SO4/SO2 including sulfate prognostic stratospheric aerosols (SO4 is 11)
!     6      13-14  Nitrate fine and coarse
!     7      15     Ammonium
!     8      16-18  SOA Biogenic and Anthropogenic 1/2
!     9      19     fly ash
!     10     20-21   volcanic SO2/SO4

DO JAER=1,NMAXTAER
  NTYPAER(JAER)=0
ENDDO


NSTPDBG=10
DO J=1,NSTPDBG
  KSTPDBG(J)=-999
ENDDO

NBCOPTP = 1
NDDOPTP = 3
NOMOPTP = 1
NSSOPTP = 1
NSUOPTP = 1
NVISWL  = 7

RSUCV1  = 0._JPRB
RSUCV2  = 0._JPRB
RDMSMIN = 0._JPRB

! Default to backwards compatibility for now
LSEASALT_RH80 = .TRUE.

IF (NAERO == 0) THEN
  LEPAERO  =.FALSE.
  LAERCLIMG=.FALSE.
  LAERCLIMZ=.FALSE.
  LAERCLIST=.FALSE.
  LAERRAD  =.FALSE.

  LAERDRYDP=.FALSE.
  LAERGTOP =.FALSE.
  LAERHYGRO=.FALSE.
  LAERNGAT =.FALSE.
  NAERSCAV = 0
  LAERDUSTSOURCE =.FALSE.
  LAERDUST_NEWBIN =.FALSE.
  LAERDUSTSIZEVAR =.FALSE.
  RAERDUST_REBOUND = 0._JPRB
  LAERSEDIM=.FALSE.
  LAERSEDIMSS=.FALSE.
  LAERSURF =.FALSE.
  LAERELVS =.FALSE.
  LAER6SDIA=.FALSE.
  LAERCCN  =.FALSE.
  LAERCSTR =.FALSE.
  LAERRRTM =.FALSE.
  LAERUVP  =.FALSE.
  LOCNDMS  =.FALSE.
  NDRYDEPVEL_DYN= 0
  NDMSO = 0
  NPIST = 0

  NAERLISI =0

  LAERDIAG1=.FALSE.
  LAERDIAG2=.FALSE.

  LAERVOL = .FALSE.
  NAERVOLC= 0
  NAERVOLE= 0
  NINTERPT= 0
  NVOLERUP= 0
  NVOLOPTP= 0
  NVOLHOMO= 0
  NVOLERUZ(:)= 0
  RAERVOLC(:,:)=0._JPRB
  RAERVOLE(:,:)=0._JPRB
  RVOLERUZ(:)= 1._JPRB

  LAERINIT = .FALSE.
  NAERWND = 0
  NSSALT = 0
  NDDUST = 0
  NSO4SCHEME = 0

  RFCTDU = 0._JPRB
  RFCTSS = 0._JPRB
  RFCTDUR = 0._JPRB
  RFCTSSR = 0._JPRB

  NDRYDEP = 0

  RDMSMIN = 0._JPRB

  IF (YDERAD%NAERMACC == 1) THEN
    LAERRRTM=.TRUE.
    CALL SU_AERP(YDEAERATM,YDMODEL%YRML_PHY_AER,YDCOMPO)
    CALL SU_AEROP(YDEAERATM,YDMODEL%YRML_PHY_AER%YREAERLID,YDEAERSRC)
  ENDIF

ELSE

 LEPAERO  =.TRUE.
 LAERRAD  =.TRUE.

 SELECT CASE (TRIM(AERO_SCHEME))
!-- define a default configuration that can be modified through the "naeaer" namelist
  CASE ("glomap")
    CALL ABOR1("OIFS - glomap should never be called from OIFS, EXIT")


  CASE ("aer")

! the 8 types and assumed number of bins are:
!  NTYPAER    bins  type
!     1       1- 3  sea-salt  0.03 - 0.5 -  5  - 20 microns
!     2       4- 6  dust      0.03 - 0.5 - 0.9 - 20 microns
!     3       7- 8  POM    hydrophilic, hydrophobic
!     4       9-10  BC     hydrophilic, hydrophobic
!     5      11-12  SO4/SO2 including sulfate prognostic stratospheric aerosols
!     (SO4 is 11)
!     6      13-14  Nitrate fine and coarse
!     7      15     Ammonium
!     8      16     fly ash
!     9      17-18   volcanic SO2/SO4

  NTYPAER(1)=3
  NTYPAER(2)=3
  NTYPAER(3)=2
  NTYPAER(4)=2

  !LAERCHEM, LAERNITRATE, LAERSOA and LAERSOA_COUPLED are in NAMCOMPO
  !and have already been read.
  IF (LAERCHEM) THEN
    !-- if LAERCHEM is set, we use SO2 oxidation from the chemistry scheme,
    !   and the aerosol scheme does not need its own tracer.
    NTYPAER(5)=1
  ELSE
    NTYPAER(5)=2
  ENDIF
  IF (LAERNITRATE) THEN
    NTYPAER(6)=2            ! Nitrate (fine and coarse)
    NTYPAER(7)=1            ! Ammonium
  ENDIF
  IF (LAERSOA) THEN
    IF (LAERSOA_COUPLED) THEN
      NTYPAER(8)=2          ! Secondary Organics (Bio, Anthro)
    ELSE
      NTYPAER(8)=4          ! Secondary Organics (Bio, Anthro) + two precursors
    ENDIF
  ENDIF


!-- if some volcanic aerosols are present (LAERVOL=true in NAERAD), default
!   configuration NTYPAER((6) and defaults values for volcanic optical
!   properties and volcanic eruption handling are set up (if nothing
!   else has been read in from namelist NAEVOL)
  IF (LAERVOL) THEN
    NTYPAER(9)=1            ! fly ash
    NTYPAER(10)=2           ! SO4 and SO2 of volcanic origin
    NVOLERUP=1              ! a priori no knowledge of the height distribution of the plume
                            ! if such knowledge exists NVOLERUP=2 and provide a proper file ...
  ELSE
    NVOLERUP=0
  ENDIF

  LAERCLIMG=.FALSE.
  LAERCLIMZ=.FALSE.
  LAERCLIST=.FALSE.

  LAERDRYDP=.TRUE.
  LAERGTOP =.TRUE.
  LAERHYGRO=.TRUE.
  LAERNGAT =.TRUE.
  NAERSCAV = 2
  LAERDUSTSOURCE =.FALSE.
  LAERDUST_NEWBIN =.FALSE.
  LAERDUSTSIZEVAR =.FALSE.
  RAERDUST_REBOUND = 0._JPRB
  LAERSEDIM=.TRUE.
  LAERSEDIMSS=.FALSE.
  LAERSURF =.TRUE.
  LAERELVS =.FALSE.
  LAER6SDIA=.FALSE.
  LAERCCN  =.FALSE.
  LAERCSTR =.FALSE.
  LAERRRTM =.FALSE.
  LAERUVP  =.TRUE.
  LOCNDMS  =.TRUE.
  NDRYDEPVEL_DYN= 0
  NDMSO = 2
  NPIST = 1

  NAERLISI  =1

  LAERDIAG1=.FALSE.
  LAERDIAG2=.FALSE.

!-- NB: if volcanic activity is to be accounted for, LAERVOL has to appear in naerad
  LAERVOL = .FALSE.
!-- and the other parameters have to be potentially defined in naervol
  NAERVOLC= 0
  NAERVOLE= 0
  NINTERPT= 0
  NVOLOPTP= 1              ! SO4 as taken default value for optp.prop. of volcanic aerosols
  NVOLHOMO= 0
  NVOLERUZ(:)= 0
  RAERVOLC(:,:)=0._JPRB
  RAERVOLE(:,:)=0._JPRB
  RVOLERUZ(:)= 1._JPRB

  LAERINIT =.FALSE.

!-- default value are for use of "plain" or "gusty" 10-m wind as predictor for SS and DU
  NAERWND  = 2
!--  other values would be: (see *aer_src*)
!- NAERWND = 0 for "plain" 10-m wind as predictor for sea salt and desert dust emissions
!- NAERWND = 1 for wind+gust in sea salt emission
!- NAERWND = 2 for wind+gust in dust emission
!- NAERWND = 3 for wind+gust in sea salt and dust emissions

!-- SSalt : 1 is Monahan et al. 1986
!           2 is a la LSCE
!           3 is Grythe et al (2014) adapted
  NSSALT =3

!-- DDust : 2 is Nabat et al 2012
!           3 is "pre-operational" with filter on StDevOrog, MODIS Albedo
!           4 is Nabat et al 2012 with roughness length of smooth erodible
!           surfaces
  NDDUST =3
!-- Sulphate sceheme : 1 is operational
!           2 is "mocage like" from MF J.Bock
  NSO4SCHEME = 1
  RAERDUB=1.E-11_JPRB   !  dust emission potential modulated by the areas' RDDUAER
!                       !  (default value when gustiness is accounted for in dust emission)
! NB: if no gustiness is accounted for (NAERWND=0, explicit need to enter RAERDUB ~5.E-11

!-- default values are for use of 10-m wind as predictor for SS and DU
  RFCTDU     = 1.0_JPRB
  RFCTSS     = 1.0_JPRB
  RFCTDUR    = 1.0_JPRB   ! in preliminary tests rfctdur=0.40
  RFCTSSR    = 1.0_JPRB   !                      rfctssr=0.52
    LAERDUSTSOURCE =.FALSE.

!-- dry deposition as in GEMS
  NDRYDEP = 1

!-- minimum oceanic production of DMS
  RDMSMIN = 5.E-11_JPRB

!-- from BEN

    LAEROMIN=.FALSE.

!- note that NAERCONF is now irrelevant
    NAERCONF=-99
    NXT3DAER=0
    NINIDAY=19000101

    ! This is no longer used by the negative fixer (as we use CHEM_NEGAT instead)
    ! - only as a minimum value in the radiation code
    REPSCAER=1.E-20_JPRB

  CASE DEFAULT

    CALL ABOR1(" NO AEROSOL SCHEME "//TRIM(AERO_SCHEME) )

 END SELECT
ENDIF

RDDUAER(:) = -9999.0_JPRB
RDDUSRC(:)= 0.0_JPRB
NDUSRCP(:) = 1
RDUSRCP(:,:) = 0.0_JPRB

!* default reference values for threshold speed and reference particle radius
! -- 1 N & S America, Europe
RDUSRCP(1,1) = 6.0_JPRB
RDUSRCP(1,2) = 5.0_JPRB
! -- 2 Russia, Urals
RDUSRCP(2,1) = 6.0_JPRB
RDUSRCP(2,2) = 5.0_JPRB
! -- 3  Africa, Sahara, S. Africa
RDUSRCP(3,1) = 6.0_JPRB
RDUSRCP(3,2) = 5.0_JPRB
! -- 4 Australasia
RDUSRCP(4,1) = 4.0_JPRB
RDUSRCP(4,2) = 5.0_JPRB
! -- 5 Asian deserts
RDUSRCP(5,1) = 3.5_JPRB
RDUSRCP(5,2) = 5.0_JPRB
! -- 6 dry lands of S.America
RDUSRCP(6,1) = 4.0_JPRB
RDUSRCP(6,2) = 5.0_JPRB
! -- 7 the rest (Japan, Greenland, Antarctica)
RDUSRCP(7,1) = 4.0_JPRB
RDUSRCP(7,2) = 5.0_JPRB

!-- default values are for use of 10-m wind as predictor for SS and DU
RFCTDU     = 1.0_JPRB
RFCTSS     = 1.0_JPRB
!* New defaults taken from previous namelist; stj - 27-10-2010
RFCTDUR    = 1.0_JPRB
RFCTSSR    = 1.0_JPRB
CLAERWND(0) = '10-M WIND AS PREDICTOR FOR SS AND DU         '
CLAERWND(1) = 'PREDICTORS: WIND GUST FOR SS, 10M-WIND FOR DU'
CLAERWND(2) = 'PREDICTORS: WIND GUST FOR DU, 10M-WIND FOR SS'
CLAERWND(3) = 'WIND GUST AS PREDICTORS FOR SS AND DU        '


!     ------------------------------------------------------------------

!*       2.       READ VALUES OF PROGNOSTIC AEROSOL CONFIGURATION
!                 -----------------------------------------------


IF(NAERO > 0) THEN
 CALL POSNAM(NULNAM,'NAEAER')
 READ (NULNAM,NAEAER)
 IF (LAERVOL) THEN
   CALL POSNAM(NULNAM,'NAEVOL')
   READ (NULNAM,NAEVOL)
 ENDIF
ENDIF



IF (NINDAT == NINIDAY .AND. NSSSSS == 00000) THEN
  LAERCLIST=.TRUE.
ENDIF

! supplementary definitions for the AER SCHEME

IF (NACTAERO > 0) THEN
  IF (TRIM(AERO_SCHEME) == "aer") THEN
    IF (ALL(YAERO_DESC(1:NACTAERO)%NTYP == -9999) .AND. ALL(YAERO_DESC(1:NACTAERO)%NBIN == -9999)) THEN
    ! Auto-initialise NTYP and NBIN if not set in the namelist
      ICAER=0
      DO JAER=1,NMAXTAER
        IF (NTYPAER(JAER) /= 0) THEN
          ITAER=NTYPAER(JAER)
          DO IAER=1,ITAER
            ICAER=ICAER+1
            IF (ICAER > NACTAERO) CALL ABOR1('SU_AERW: more bins in NTYPAER than NACTAERO')
            YAERO_DESC(ICAER)%NTYP=JAER
            YAERO_DESC(ICAER)%NBIN=IAER
          ENDDO
        ENDIF
      ENDDO
      IF (ICAER > NACTAERO) CALL ABOR1('SU_AERW: fewer bins in NTYPAER than NACTAERO')
    ENDIF
  ENDIF

  !     ------------------------------------------------------------------

  !*       3.       DISTRIBUTE DUST AEROSOL SOURCE FUNCTIONS
  !                 ----------------------------------------

  ! This is why we use ZDDUAER rather than initialising RDDUAER directly:
  !   we need the NAEAER to have been read so that NDDUST, LAERDUST* are
  !   available, but we also need to be able to override individual RDDUAER
  !   elements via the same namelist.
  ZDDUAER(:) = 0._JPRB
  SELECT CASE (NDDUST)
    CASE(3,9)
      !* Default dust emission factors for Ginoux scheme
      ZDDUAER(1)=1.0_JPRB  ! Canada
      ZDDUAER(2)=1.0_JPRB  ! Alaska
      ZDDUAER(3)=0.5_JPRB  ! USA (or 1.0)
      ZDDUAER(4)=0.6_JPRB  ! Central America
      ZDDUAER(5)=0.6_JPRB  ! South America
      ZDDUAER(6)=1.0_JPRB  ! Brazil
      ZDDUAER(7)=1.0_JPRB  ! Iceland
      ZDDUAER(8)=1.0_JPRB  ! Ireland
      ZDDUAER(9)=1.0_JPRB  ! Britian
      ZDDUAER(10)=1.0_JPRB ! Continental Europe
      ZDDUAER(11)=1.0_JPRB ! Russia (Europe)
      ZDDUAER(12)=1.0_JPRB ! Russia (Georgia)
      ZDDUAER(13)=0.0_JPRB ! Northern Sahara (now 34-37)
      ZDDUAER(14)=0.2_JPRB ! Central Africa (West)
      ZDDUAER(15)=0.5_JPRB ! Southern Africa
      ZDDUAER(16)=1.0_JPRB ! Siberia
      ZDDUAER(17)=1.0_JPRB ! Asian deserts
      ZDDUAER(18)=0.5_JPRB ! Saudi Arabia
      ZDDUAER(19)=0.5_JPRB ! Iraq, Iran, Pakistan
      ZDDUAER(20)=0.8_JPRB ! Central Asia Taklamakan (or 0.4)
      ZDDUAER(21)=1.2_JPRB ! India (or 0.6)
      ZDDUAER(22)=1.5_JPRB ! Mongolia and Gobi
      ZDDUAER(23)=1.5_JPRB ! Central China
      ZDDUAER(24)=0.5_JPRB ! South China
      ZDDUAER(25)=1.0_JPRB ! Japan, South Korea
      ZDDUAER(26)=0.5_JPRB ! Padding Asia
      ZDDUAER(27)=1.0_JPRB ! Tropical Pacific Islands
      ZDDUAER(28)=0.5_JPRB ! Australia, New Zealand
      ZDDUAER(29)=1.0_JPRB ! Greenland
      ZDDUAER(30)=1.0_JPRB ! Antarctica
      ZDDUAER(31)=0.1_JPRB ! Atacama and Uyuni
      ZDDUAER(32)=0.1_JPRB ! Pipanco and others
      ZDDUAER(33)=0.3_JPRB ! Argentinian pampas
      ZDDUAER(34)=0.8_JPRB ! Southern Sahara (West 1)
      ZDDUAER(35)=0.5_JPRB ! Southern Sahara (East)
      ZDDUAER(36)=0.4_JPRB ! Northern Sahara (West)
      ZDDUAER(37)=0.6_JPRB ! Northern Sahara (East)
      ZDDUAER(38)=0.2_JPRB ! Central Africa (East)
      ZDDUAER(39)=0.8_JPRB ! Southern Sahara (West 2)

    CASE(2,4)
      ! Regional tuning parameters for P Nabat's scheme
      ! Depend whether surface smoothness is taken into account (NDDUST 4) or not
      ! Depend on the dust bin definition
      ! Depend on whether a dust source function is used (LAERDUSTSOURCE)
      !
      ZDDUAER(:) = 1._JPRB
      IF (.NOT.LAERDUST_NEWBIN) THEN
        ZDDUAER(28) = 0.03_JPRB
        ZDDUAER(31) = 0.1_JPRB
        ZDDUAER(32) = 0.1_JPRB
        ZDDUAER(33) = 0.1_JPRB
        ZDDUAER(14) = 0.5_JPRB
        ZDDUAER(15) = 0.2_JPRB
        ZDDUAER(18) = 0.25_JPRB
        ZDDUAER(20) = 2.0_JPRB
        !ZDDUAER(22) = 1.5_JPRB
        !ZDDUAER(23) = 1.5_JPRB
        ZDDUAER(36) = 0.3_JPRB
        ZDDUAER(37) = 0.3_JPRB
        ZDDUAER(3) = 0.3_JPRB
        ZDDUAER(4) = 0.6_JPRB
        ZDDUAER(19) = 0.8_JPRB
        ZDDUAER(21) = 0.4_JPRB
        ZDDUAER(4) = 0.6_JPRB
        IF (LAERDUSTSOURCE) THEN
          ZDDUAER(28) = 1.0_JPRB
          ZDDUAER(38) = 0.4_JPRB
          ZDDUAER(37) = 0.45_JPRB
          ZDDUAER(18) = 0.4_JPRB
          ZDDUAER(39) = 1.3_JPRB
        ENDIF
      ELSE
        IF (NDDUST == 4) THEN
          IF (LAERDUSTSOURCE) THEN
            ZDDUAER(38) = 0.2_JPRB
            ZDDUAER(31) = 0.2_JPRB
            ZDDUAER(32) = 0.2_JPRB
            ZDDUAER(33) = 0.2_JPRB
            ZDDUAER(20) = 2.0_JPRB
            ZDDUAER(22) = 2.5_JPRB
            ZDDUAER(39) = 0.5_JPRB
            ZDDUAER(18) = 0.4_JPRB
            ZDDUAER(36) = 0.4_JPRB
            ZDDUAER(37) = 0.4_JPRB
            ZDDUAER(3) = 0.4_JPRB
            ZDDUAER(4) = 0.4_JPRB
          ELSE
            ZDDUAER(28) = 0.025_JPRB
            ZDDUAER(31) = 0.05_JPRB
            ZDDUAER(32) = 0.05_JPRB
            ZDDUAER(33) = 0.05_JPRB
            ZDDUAER(14) = 0.3_JPRB
            ZDDUAER(15) = 0.2_JPRB
            ZDDUAER(18) = 0.3_JPRB
            ZDDUAER(20) = 2.0_JPRB
            ZDDUAER(22) = 0.8_JPRB
            !ZDDUAER(23) = 1.5_JPRB
            ZDDUAER(36) = 0.5_JPRB
            ZDDUAER(37) = 0.5_JPRB
            ZDDUAER(39) = 0.3_JPRB
            ZDDUAER(3) = 0.3_JPRB
            ZDDUAER(4) = 0.3_JPRB
            ZDDUAER(19) = 0.8_JPRB
            ZDDUAER(21) = 0.5_JPRB
          ENDIF
        ELSEIF (NDDUST == 2) THEN
          IF (LAERDUSTSOURCE) THEN
            ZDDUAER(38) = 0.2_JPRB
            ZDDUAER(31) = 0.2_JPRB
            ZDDUAER(32) = 0.2_JPRB
            ZDDUAER(33) = 0.2_JPRB
            ZDDUAER(20) = 2.5_JPRB
            ZDDUAER(22) = 3.5_JPRB
            ZDDUAER(39) = 0.6_JPRB
            !ZDDUAER(18) = 0.4_JPRB
            ZDDUAER(18) = 0.3_JPRB
            !ZDDUAER(36) = 0.4_JPRB
            ZDDUAER(36) = 0.1_JPRB
            ZDDUAER(37) = 0.4_JPRB
            !ZDDUAER(3) = 0.4_JPRB
            ZDDUAER(3) = 0.07_JPRB
            !ZDDUAER(4) = 0.4_JPRB
            ZDDUAER(4) = 0.07_JPRB
          ELSE
            ZDDUAER(28) = 0.03_JPRB
            ZDDUAER(31) = 0.1_JPRB
            ZDDUAER(32) = 0.1_JPRB
            ZDDUAER(33) = 0.1_JPRB
            ZDDUAER(14) = 0.5_JPRB
            ZDDUAER(15) = 0.2_JPRB
            ZDDUAER(18) = 0.25_JPRB
            ZDDUAER(20) = 2.0_JPRB
            !ZDDUAER(22) = 1.5_JPRB
            !ZDDUAER(23) = 1.5_JPRB
            ZDDUAER(36) = 0.3_JPRB
            ZDDUAER(37) = 0.3_JPRB
            ZDDUAER(39) = 0.6_JPRB
            ZDDUAER(3) = 0.3_JPRB
            ZDDUAER(4) = 0.6_JPRB
            ZDDUAER(19) = 0.8_JPRB
            ZDDUAER(21) = 0.3_JPRB
            ZDDUAER(4) = 0.6_JPRB
          ENDIF
        ENDIF
      ENDIF
  END SELECT

  WHERE(RDDUAER(:) < 0.0_JPRB) RDDUAER(:) = ZDDUAER(:)

  IF (TRIM(AERO_SCHEME) == "GLOMAP" .OR. NTYPAER(2) /= 9) THEN
    RDDUSRC(1)=0.3_JPRB
    ! online dry deposition increases a lot dry dep for dust bin 2 and 3 => need
    ! to increase sources to compensate
    IF (NDRYDEPVEL_DYN > 0) THEN
      RDDUSRC(2)=1._JPRB
      RDDUSRC(3)=7._JPRB
    ELSE
      RDDUSRC(2)=0.8_JPRB
      RDDUSRC(3)=5.5_JPRB
    ENDIF
  ELSE
    RDDUSRC(:)=1.0_JPRB
  ENDIF

ENDIF ! NACTAERO > 0

!     ------------------------------------------------------------------

!*       4.       DEFINE VALUES OF PROGNOSTIC AEROSOL CONFIGURATION
!                 -------------------------------------------------

WRITE(UNIT=NULOUT,FMT='(''NAERO='',I2,'' NACTAERO='',I2)') NAERO,NACTAERO

IF (NACTAERO > 0) THEN

  WRITE(UNIT=NULOUT,FMT='(''NAERO='',I2,'' NACTAERO='',I2,3X,10I3)') NAERO,NACTAERO,(NTYPAER(JAER),JAER=1,10)

!-- if MACC model includes volcanic aerosol, no startospheric adjustment is included
  IF (LAERVOL) THEN
    LAERCSTR=.FALSE.
  ENDIF
!     ------------------------------------------------------------------

!*       5.    INITIALIZE PROGNOSTIC AEROSOL PHYSICAL AND OPTICAL PARAMETERS
!              -------------------------------------------------------------

  CALL SU_AERP(YDEAERATM,YDMODEL%YRML_PHY_AER,YDCOMPO)
  CALL SU_AEROP(YDEAERATM,YDMODEL%YRML_PHY_AER%YREAERLID,YDEAERSRC)

  IF (TRIM(AERO_SCHEME) == "aer") THEN
! Pre-compute conversion factors for dry vs. 80%-RH sea salt
    IF (LSEASALT_RH80) THEN
  ! Factors to convert transported SS (@ 80% RH) to dry sea salt
  RSS_DRY_DIAFAC = 1._JPRB/RSSGROWTH_RHTAB(9)
  RSS_DRY_DENSFAC = RSSDENS_RHTAB(1)/RSSDENS_RHTAB(9)
  RSS_DRY_MASSFAC = RSS_DRY_DENSFAC * (RSS_DRY_DIAFAC**3._JPRB)
  RSS_RH80_DIAFAC = 1.0_JPRB
  RSS_RH80_DENSFAC = 1.0_JPRB
  RSS_RH80_MASSFAC = 1.0_JPRB
    ELSE
  RSS_DRY_DIAFAC = 1.0_JPRB
  RSS_DRY_DENSFAC = 1.0_JPRB
  RSS_DRY_MASSFAC = 1.0_JPRB
  ! Factors to convert transported SS (dry) to sea salt @ 80% RH
  RSS_RH80_DIAFAC = RSSGROWTH_RHTAB(9)
  RSS_RH80_DENSFAC = RSSDENS_RHTAB(9)/RSSDENS_RHTAB(1)
  RSS_RH80_MASSFAC = RSS_RH80_DENSFAC * (RSS_RH80_DIAFAC**3._JPRB)
    ENDIF
  ENDIF

!      ----------------------------------------------------------------

!*       6.    PRINT FINAL VALUES.
!              -------------------

  WRITE(UNIT=NULOUT,FMT='('' AERO_SCHEME   = '',A10 &
   & ,'' LEPAERO = '',L5 &
   & ,'' LAERNITRATE = '',L5 &
   & ,'' NMAXTAER = '',I2 ,'' NDDUST = '',I1 &
   & ,'' NSSALT = '',I1,'' NSO4SCHEME = '',I1 &
   & ,'' NINIDAY = '',I8)') &
   & AERO_SCHEME,LEPAERO,LAERNITRATE,NMAXTAER,NDDUST,NSSALT,NSO4SCHEME,NINIDAY

  WRITE(UNIT=NULOUT,FMT='('' NAERO = '',I2,'' NACTAERO = '',I2,&
   &'' NXT3DAER = '',I2,'' NAERCONF = '',I3)') &
   & NAERO,NACTAERO,NXT3DAER,NAERCONF

  WRITE(UNIT=NULOUT,FMT='('' NTYPAER = '',*(I3))') (NTYPAER(J), J=1,NMAXTAER)
  WRITE(UNIT=NULOUT,FMT='('' NTYP    = '',*(I3))') (YAERO_DESC(J)%NTYP, J=1,NACTAERO)
  WRITE(UNIT=NULOUT,FMT='('' NBIN    = '',*(I3))') (YAERO_DESC(J)%NBIN, J=1,NACTAERO)

  WRITE(UNIT=NULOUT,FMT='('' LAERINIT = '',L1)') LAERINIT
  WRITE(UNIT=NULOUT,FMT='(" LAERSURF = ",L1 &
   & ,'' LAERELVS = '',L1 &
   & ,'' LOCNDMS  = '',L1 &
   & ,'' RDMSMIN  = '',E10.3 &
   & ,'' NDMSO = '',I1 &
   & ,'' NPIST = '',I1 &
   & )')&
   & LAERSURF,LAERELVS, LOCNDMS, RDMSMIN, NDMSO, NPIST
  WRITE(UNIT=NULOUT,FMT='('' LAERGBUD = '',L1 &
   & ,'' LAERNGAT = '',L1 &
   & ,'' LAERDRYDP= '',L1 &
   & ,'' LAERSEDIM= '',L1 &
   & ,'' LAERSEDIMSS= '',L1 &
   & ,'' LAERGTOP = '',L1 &
   & ,'' LAERHYGRO= '',L1 &
   & ,'' NAERSCAV = '',I1&
   & ,'' LAERCHEM = '',L1&
   & ,'' LAERNITRATE = '',L1&
   & ,'' LAERSOA = '',L1&
   & ,'' LAERSOA_COUPLED = '',L1&
   & ,'' NDRYDEPVEL_DYN = '',I1&
   & ,'' LAERDUSTSOURCE = '',L1&
   & ,'' LAERDUST_NEWBIN = '',L1&
   & ,'' LAERDUSTSIZEVAR = '',L1&
   & ,'' RAERDUST_REBOUND = '',E3.2&
   & ,'' LAER6SDIA= '',L1&
   & ,'' LAERCLIMZ= '',L1&
   & ,'' LAERCLIMG= '',L1&
   & ,'' LAERCLIST= '',L1&
   & ,'' LAERRAD= '',L1&
   & )')&
   & LAERGBUD,LAERNGAT,LAERDRYDP,LAERSEDIM,LAERSEDIMSS,LAERGTOP,LAERHYGRO,NAERSCAV,LAERCHEM,LAERNITRATE,LAERSOA,LAERSOA_COUPLED, &
   & NDRYDEPVEL_DYN,LAERDUSTSOURCE,LAERDUST_NEWBIN,LAERDUSTSIZEVAR,RAERDUST_REBOUND,LAER6SDIA,&
   & LAERCLIMZ,LAERCLIMG,LAERCLIST,LAERRAD
  WRITE(UNIT=NULOUT,FMT='('' NAERLISI = '',I1)') NAERLISI
  WRITE(UNIT=NULOUT,FMT='('' LAERDIAG1 = '',L1&
   & ,'' LAERDIAG2 = '',L1&
   & )')&
   & LAERDIAG1,LAERDIAG2

  WRITE(UNIT=NULOUT,FMT='('' RSSFLX= '',9E10.3)') RSSFLX
  WRITE(UNIT=NULOUT,FMT='('' NAERWND= '',I1,'' RFCTSS= '',F4.1,'' RFCTDU= '',F4.1,2X,A45)') &
   & NAERWND, RFCTSS, RFCTDU, CLAERWND(NAERWND)
  WRITE(UNIT=NULOUT,FMT='('' RAERDUB= '',E10.3)') RAERDUB
  WRITE(UNIT=NULOUT,FMT='('' RDDUAER= '',25F6.3)') (RDDUAER(J),J= 1,25)
  WRITE(UNIT=NULOUT,FMT='('' NDUSRCP= '',25I6  )') (NDUSRCP(J),J= 1,25)
  WRITE(UNIT=NULOUT,FMT='('' RDDUAER= '',25F6.3)') (RDDUAER(J),J=26,50)
  WRITE(UNIT=NULOUT,FMT='('' NDUSRCP= '',25I6  )') (NDUSRCP(J),J=26,50)
  WRITE(UNIT=NULOUT,FMT='('' LAERUVP= '',L3)') LAERUVP

  IF (TRIM(AERO_SCHEME) == "aer") THEN

    WRITE(UNIT=NULOUT,FMT='('' LAERGBUD = '',L1 &
     & ,'' LAERNGAT = '',L1 &
     & ,'' LAERDRYDP= '',L1 &
     & ,'' LAERSEDIM= '',L1 &
     & ,'' LAERSEDIMSS= '',L1 &
     & ,'' LAERGTOP = '',L1 &
     & ,'' LAERHYGRO= '',L1 &
     & ,'' NAERSCAV = '',I1 &
     & ,'' LAER6SDIA= '',L1 &
     & ,'' LAERCLIMZ= '',L1 &
     & ,'' LAERCLIMG= '',L1 &
     & ,'' LAERCLIST= '',L1 &
     & )')&
     & LAERGBUD,LAERNGAT,LAERDRYDP,LAERSEDIM,LAERSEDIMSS,LAERGTOP,LAERHYGRO, &
     & NAERSCAV,LAER6SDIA, &
     & LAERCLIMZ,LAERCLIMG,LAERCLIST
    WRITE(UNIT=NULOUT,FMT='('' NAERLISI = '',I1 &
     & )')&
     & NAERLISI
    WRITE(UNIT=NULOUT,FMT='('' LAERDIAG1 = '',L1 &
     & ,'' LAERDIAG2 = '',L1 &
     & )')&
     & LAERDIAG1,LAERDIAG2

    WRITE(UNIT=NULOUT,FMT='('' RSSFLX= '',9E10.3)') RSSFLX
    WRITE(UNIT=NULOUT,FMT='('' NAERWND= '',I1,'' RFCTSS= '',F4.1,'' RFCTDU= '',F4.1,2X,A45)') &
     & NAERWND, RFCTSS, RFCTDU, CLAERWND(NAERWND)
    WRITE(UNIT=NULOUT,FMT='('' LAERUVP= '',L3)') LAERUVP

    WRITE(UNIT=NULOUT,FMT='('' Interaction prognostic aerosols and radiation: '',&
    &'' LAERCSTR= '',L1,'' LAERRRTM= '',L1)') LAERCSTR,LAERRRTM
    WRITE(UNIT=NULOUT,FMT='('' Progn. aerosol optical properties: NBCOPTP = '',I1 &
     & ,'' NDDOPTP = '',I1 &
     & ,'' NOMOPTP = '',I1 &
     & ,'' NSSOPTP = '',I1 &
     & ,'' NSUOPTP = '',I1 &
   & )')&
   & NBCOPTP,NDDOPTP,NOMOPTP,NSSOPTP,NSUOPTP

    WRITE(UNIT=NULOUT,FMT='('' Interaction prognostic aerosols and eff.radius of liq.wat.clouds: LAERCCN= '',&
    &L1)') LAERCCN
    WRITE(UNIT=NULOUT,FMT='(" Visibility calculations for wavelength NVISWL= ",I2)') NVISWL
    WRITE(UNIT=NULOUT,FMT='(" Aerosol being dry deposited: NDRYDEP=",I2)') NDRYDEP
  ENDIF !AERO_SCHEME=aer

  IF (LAERVOL) THEN
    WRITE(UNIT=NULOUT,FMT='(" LAERVOL= ",L2,"  NVOLOPTP=",I3)') LAERVOL,NVOLOPTP
      WRITE(UNIT=NULOUT,FMT='(" NAERVOLC=",I3,"  NAERVOLE=",I3,"  NVOLERUP =",I3)') &
     & NAERVOLC,NAERVOLE,NVOLERUP
    IF (NAERVOLC /=  0) THEN
      WRITE(UNIT=NULOUT,FMT='(" VOLC lat ",10F12.3)') (RAERVOLC(J,1),J= 1,NAERVOLC)
      WRITE(UNIT=NULOUT,FMT='(" VOLC lon ",10F12.3)') (RAERVOLC(J,2),J= 1,NAERVOLC)
      WRITE(UNIT=NULOUT,FMT='(" VOLC mssa",10F12.1)') (RAERVOLC(J,3),J= 1,NAERVOLC)
      WRITE(UNIT=NULOUT,FMT='(" VOLC mssu",10E12.5)') (RAERVOLC(J,4),J= 1,NAERVOLC)
      WRITE(UNIT=NULOUT,FMT='(" VOLC bsp ",10E12.5)') (RAERVOLC(J,5),J= 1,NAERVOLC)
      WRITE(UNIT=NULOUT,FMT='(" VOLC tpp ",10E12.5)') (RAERVOLC(J,6),J= 1,NAERVOLC)
      WRITE(UNIT=NULOUT,FMT='(" VOLC fAsh",10E12.5)') (RAERVOLC(J,7),J= 1,NAERVOLC)
      WRITE(UNIT=NULOUT,FMT='(" VOLC fSO2",10E12.5)') (RAERVOLC(J,8),J= 1,NAERVOLC)
    ENDIF
    IF (NAERVOLE /=  0) THEN
      WRITE(UNIT=NULOUT,FMT='(" VOLE lat ",10F12.3)') (RAERVOLE(J,1),J= 1,NAERVOLE)
      WRITE(UNIT=NULOUT,FMT='(" VOLE lon ",10F12.3)') (RAERVOLE(J,2),J= 1,NAERVOLE)
      WRITE(UNIT=NULOUT,FMT='(" VOLE mssa",10F12.1)') (RAERVOLE(J,3),J= 1,NAERVOLE)
      WRITE(UNIT=NULOUT,FMT='(" VOLE mssu",10E12.5)') (RAERVOLE(J,4),J= 1,NAERVOLE)
      WRITE(UNIT=NULOUT,FMT='(" VOLE bsp ",10E12.5)') (RAERVOLE(J,5),J= 1,NAERVOLE)
      WRITE(UNIT=NULOUT,FMT='(" VOLE tpp ",10E12.5)') (RAERVOLE(J,6),J= 1,NAERVOLE)
      WRITE(UNIT=NULOUT,FMT='(" VOLE fAsh",10E12.5)') (RAERVOLE(J,7),J= 1,NAERVOLE)
      WRITE(UNIT=NULOUT,FMT='(" VOLE fSO2",10E12.5)') (RAERVOLE(J,8),J= 1,NAERVOLE)
      IF (NVOLERUP /= 0) THEN
        DO JVOLE=1,NAERVOLE
          IF (NVOLERUZ(JVOLE) == 1) THEN
             WRITE(UNIT=NULOUT,FMT='("PROFILE ERUPTION DATA AVAILABLE FOR VOLCANO No ",I3,&
             &" WITH SCALING FACTOR=",F5.2," FROM TIME: ",I10," TO TIME: ",I10)') &
             & JVOLE,RVOLERUZ(JVOLE),NVOLDATS(JVOLE),NVOLDATE(JVOLE)
            WRITE(UNIT=NULOUT,FMT='(" NVOLHOMO=",I3," NINTERPT=",I3)') NVOLHOMO,NINTERPT
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDIF

  WRITE(UNIT=NULOUT,FMT='('' LSEASALT_RH80 = '',L1)') LSEASALT_RH80
  WRITE(UNIT=NULOUT,FMT='('' RSS_DRY_DIAFAC = '',F6.3 &
   & ,'' RSS_DRY_DENSFAC = '',F6.3 &
   & ,'' RSS_DRY_MASSFAC = '',F6.3)') &
   & RSS_DRY_DIAFAC, RSS_DRY_DENSFAC, RSS_DRY_MASSFAC
  WRITE(UNIT=NULOUT,FMT='('' RSS_RH80_DIAFAC= '',F6.3 &
   & ,'' RSS_RH80_DENSFAC= '',F6.3 &
   & ,'' RSS_RH80_MASSFAC= '',F6.3)') &
   & RSS_RH80_DIAFAC, RSS_RH80_DENSFAC, RSS_RH80_MASSFAC

ELSEIF (LESO4HIS) THEN
  WRITE(NULOUT,*)"SU_AERW: INITIALIZE AEROSOL EXTINCTION COEFFICIENTS FOR HISTORICAL AEROSOLS"
  CALL SU_AERP(YDEAERATM,YDMODEL%YRML_PHY_AER,YDCOMPO)
  CALL SU_AEROP(YDEAERATM,YDMODEL%YRML_PHY_AER%YREAERLID,YDEAERSRC)
ENDIF

!     ----------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SU_AERW',1,ZHOOK_HANDLE)
END SUBROUTINE SU_AERW
