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

SUBROUTINE SUMCCLAG(YDGEM,YDML_GCONF,YDML_AOC,YDCOMPO,YDCHEM,YDEAERSRC,YDEPHY,KULOUT,YDSURF)

!**** *SUMCCLAG* * - ROUTINE TO INITIALIZE SWITCHES FOR THE CLIMATE VERSION

!     PURPOSE.
!     --------
!        SET DEFAULT VALUES, THEN READS NAMELIST NAMMCC

!**   INTERFACE.
!     ----------
!        *CALL* *SUMCCLAG(...)*

!     EXPLICIT ARGUMENTS :  KULOUT
!     --------------------

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMMCC

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!      Michel Deque *CNRM*
!      ORIGINAL : 92-09-15

!     MODIFICATIONS.
!     --------------
!      Modified by M.Hamrud     : 04-10-01 CY28 Cleaning
!      Modified by D. Giard     : 04-09-15 new 923 configurations allowed
!      Modified by A.Alias      : GMEGEC/EAC list of modifications
!      JPh. Piedelievre : coupled model in distributed memory mode.
!                                 coupled model in UPDCPL mode (LMCC05,NFRCPL)
!      JJMorcrette, 20060605 MODIS albedo
!      JJMorcrette, 20061026 DU, BC, OM, SU, VOL, SOA climatological fields
!      S. Boussetta/G.Balsamo      May 2009   (Add switch for variable LAI: LLELAIV)
!      A.Alias    ,20091014 new keys LCURR and LGELAT added
!                           MPI1 + ocean currents (J-F. Gueremy)
!                           Print of NOACOMM added (E.Maisonnave)
!      K. Yessad (Jan 2012): remove old 923
!      K. Yessad (July 2014): NSTOP dependent calculations moved from SUMCC to SUMCCLAG.
!      J. Flemming (Jan 2014) Climate run for MACC / composition (input of emissions and dry deposition velocities)
!      F. Vana  05-Mar-2015  Support for single precision
!      R. Hogan 14-Jan-2019  Changed LE4ALB to NALBEDOSCHEME
!      R. Hogan 22-Feb-2019  Added 6 new albedo components
!      ----------------------------------------------------------------

USE MODEL_ATMOS_OCEAN_COUPLING_MOD , ONLY : MODEL_ATMOS_OCEAN_COUPLING_TYPE
USE MODEL_GENERAL_CONF_MOD         , ONLY : MODEL_GENERAL_CONF_TYPE
USE SURFACE_FIELDS_MIX             , ONLY : TSURF
USE YOMGEM   , ONLY : TGEM
USE PARKIND1 , ONLY : JPRD, JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY   , ONLY : TEPHY
USE YOMCT0   , ONLY : LALLOPR
USE YOMMP0   , ONLY : MYPROC, NPRINTLEV
USE YOMLUN   , ONLY : NULOUT
USE YOMCOMPO , ONLY : TCOMPO
USE YOMCHEM  , ONLY : TCHEM
USE YOEAERSRC, ONLY : TEAERSRC
USE YOM_GRIB_CODES, ONLY : NGRBAERLS13, NGRBAERSO2DD
USE CPLNG    , ONLY : CPLNG_INIT

!      ----------------------------------------------------------------

IMPLICIT NONE

TYPE(TGEM)                           , INTENT(IN)    :: YDGEM
TYPE(MODEL_GENERAL_CONF_TYPE)        , INTENT(INOUT) :: YDML_GCONF
TYPE(MODEL_ATMOS_OCEAN_COUPLING_TYPE), INTENT(INOUT) :: YDML_AOC
TYPE(TCOMPO)                         , INTENT(INOUT) :: YDCOMPO
TYPE(TCHEM)                          , INTENT(IN)    :: YDCHEM
TYPE(TEAERSRC)                       , INTENT(IN)    :: YDEAERSRC
TYPE(TEPHY)                          , INTENT(INOUT) :: YDEPHY
INTEGER(KIND=JPIM)                   , INTENT(IN)    :: KULOUT 
TYPE(TSURF)                          , OPTIONAL, INTENT(IN) :: YDSURF
INTEGER(KIND=JPIM)               :: IOPROC, ISTEP, ITIM, IT, IC, ICHEMFLX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL    :: LLP, LLELAIV, LLMODISALBEDO

!      ----------------------------------------------------------------

#include "abor1.intfb.h"
#include "inicou.intfb.h"
#include "su_coup_com.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUMCCLAG',0,ZHOOK_HANDLE)
ASSOCIATE(NACTAERO=>YDML_GCONF%YGFL%NACTAERO, NCHEM=>YDML_GCONF%YGFL%NCHEM, YCHEM=>YDML_GCONF%YGFL%YCHEM, & 
 & NGHG=>YDML_GCONF%YGFL%NGHG, LAERCHEM=>YDML_GCONF%YGFL%LAERCHEM, &
 & NEMIS2D_DESC=>YDCOMPO%NEMIS2D_DESC, YEMIS2D_DESC=>YDCOMPO%YEMIS2D_DESC, &
 & NEMIS2D=>YDML_GCONF%YGFL%NEMIS2D, &
 & NEMIS2DAUX=>YDML_GCONF%YGFL%NEMIS2DAUX, YEMIS2DAUX_DESC=>YDCOMPO%YEMIS2DAUX_DESC, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME, &
 & NDMSO=>YDEAERSRC%NDMSO, LOCNDMS=>YDEAERSRC%LOCNDMS, &
 & NALBEDOSCHEME=>YDEPHY%NALBEDOSCHEME, LECURR=>YDEPHY%LECURR, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & LMCC03=>YDML_AOC%YRMCC%LMCC03, &
 & LMCC04=>YDML_AOC%YRMCC%LMCC04, LMCC05=>YDML_AOC%YRMCC%LMCC05, LMCCEC=>YDML_AOC%YRMCC%LMCCEC, &
 & LMCCIEC=>YDML_AOC%YRMCC%LMCCIEC, LMCCIEC_COMPO=>YDML_AOC%YRMCC%LMCCIEC_COMPO, &
 & LMCC_COMPO=>YDML_AOC%YRMCC%LMCC_COMPO, LNEMO1WAY=>YDML_AOC%YRMCC%LNEMO1WAY, &
 & LNEMOCOUP=>YDML_AOC%YRMCC%LNEMOCOUP, LNEMOLIMGET=>YDML_AOC%YRMCC%LNEMOLIMGET, &
 & LNEMOLIMPUT=>YDML_AOC%YRMCC%LNEMOLIMPUT, NCLIGC=>YDML_AOC%YRMCC%NCLIGC, &
 & NCLIGC_COMPO=>YDML_AOC%YRMCC%NCLIGC_COMPO,NYSDMP_COMPO=>YDML_AOC%YRMCC%NYSDMP_COMPO, NCLIMR=>YDML_AOC%YRMCC%NCLIMR, &
 & NCLIMR_COMPO=>YDML_AOC%YRMCC%NCLIMR_COMPO, NCLIMTOT=>YDML_AOC%YRMCC%NCLIMTOT, &
 & NDIFC=>YDML_AOC%YRMCC%NDIFC, NDIFC_COMPO=>YDML_AOC%YRMCC%NDIFC_COMPO, NFRCPL=>YDML_AOC%YRMCC%NFRCPL, &
 & NJDCR=>YDML_AOC%YRMCC%NJDCR, NJDCR1=>YDML_AOC%YRMCC%NJDCR1, NJDCR1_COMPO=>YDML_AOC%YRMCC%NJDCR1_COMPO, &
 & NJDCR_COMPO=>YDML_AOC%YRMCC%NJDCR_COMPO, NOACOMM=>YDML_AOC%YRMCC%NOACOMM, NP1=>YDML_AOC%YRMCC%NP1, &
 & NP2=>YDML_AOC%YRMCC%NP2, NPCOMPO_1=>YDML_AOC%YRMCC%NPCOMPO_1, NPCOMPO_2=>YDML_AOC%YRMCC%NPCOMPO_2, &
 & NSTOP=>YDML_GCONF%YRRIP%NSTOP, TSTEP=>YDML_GCONF%YRRIP%TSTEP,&
 & NPCO2_1=>YDML_AOC%YRMCC%NPCO2_1, NPCO2_2=>YDML_AOC%YRMCC%NPCO2_2, LEOCLAKE=>YDEPHY%LEOCLAKE, &
 & LEOBC=>YDEPHY%LEOBC, LEOBCICE=>YDEPHY%LEOBCICE, &
 & KCHEM_DRYDEP=>YDCHEM%KCHEM_DRYDEP, &
 & LNEEONLINE=>YDEPHY%LNEEONLINE, YDCOM=>YDML_AOC%YRCOM, YDMCC=>YDML_AOC%YRMCC) 
!      ----------------------------------------------------------------

!*       1.    CALCULATIONS.
!              -------------

LLMODISALBEDO = (NALBEDOSCHEME > 0)

IF (NFRCPL <= 0) THEN
  NFRCPL=NSTOP+1
  WRITE(KULOUT,*) ' IN SUMCCLAG, NFRCPL REPLACED BY : ',NFRCPL
ENDIF

LLP = NPRINTLEV >= 1.OR. LALLOPR
IF(LMCC05.AND.NFRCPL <= NSTOP) THEN
  ALLOCATE(YDCOM%OMLDTH(NGPTOT))
  IF(LLP)WRITE(KULOUT,9) 'OMLDTH    ',SIZE(YDCOM%OMLDTH),SHAPE(YDCOM%OMLDTH)
  ALLOCATE(YDCOM%GTTLIN(NGPTOT))
  IF(LLP)WRITE(KULOUT,9) 'GTTLIN    ',SIZE(YDCOM%GTTLIN),SHAPE(YDCOM%GTTLIN)
  ALLOCATE(YDCOM%SSTPRE(NGPTOT))
  IF(LLP)WRITE(KULOUT,9) 'SSTPRE    ',SIZE(YDCOM%SSTPRE),SHAPE(YDCOM%SSTPRE)
  ALLOCATE(YDCOM%SSTMSK(NGPTOT))
  IF(LLP)WRITE(KULOUT,9) 'SSTMSK    ',SIZE(YDCOM%SSTMSK),SHAPE(YDCOM%SSTMSK)
ENDIF
! To initialize the YOMCOM common used by the COM coupling.

NCLIGC(:)=-999
LLELAIV = .TRUE.

! Some criminal copy-and-paste going on here - ought to be
! rationalized...

IF (LMCCEC) THEN  !BOUNDARY CONDITIONS UPDATED BY ECMWF ROUTINES
  IF (LMCC04.AND.(.NOT.LNEMO1WAY)) THEN !OASIS COUPLER ACTIVE IN ECMWF ROUTINES
    IF((.NOT.LLMODISALBEDO) .AND. (.NOT.LLELAIV) .AND. (.NOT.LECURR) ) THEN !NO( MODIS albedo + variable LAI + ocean current)
        NCLIMTOT=3
        NCLIMR=1
        NCLIGC(1)=174 ! AL - Broadband albedo
        NCLIGC(2)=139 ! STL1 - Soil temperature level 1
        NCLIGC(3)=31  ! CI - Sea ice area fraction
    ENDIF
    IF((.NOT.LLMODISALBEDO) .AND. (.NOT.LLELAIV) .AND. (LECURR) ) THEN !NO( MODIS albedo + variable LAI) + ocean current
        NCLIMTOT=5
        NCLIMR=1
        NCLIGC(1)=174
        NCLIGC(2)=139
        NCLIGC(3)=31
        NCLIGC(4)=131 ! U - U component of wind
        NCLIGC(5)=132 ! V - V component of wind
    ENDIF
    IF((LLMODISALBEDO) .AND. (.NOT.LLELAIV) .AND. (.NOT.LECURR) ) THEN ! MODIS albedo + No(variable LAI + ocean current)
        NCLIMTOT=13
        NCLIMR=11
        NCLIGC(1)=174
        NCLIGC(2)=15 ! ALUVP - UV/Vis albedo for direct radiation
        NCLIGC(3)=16 ! ALUVD - UV/Vis albedo for diffuse radiation
        NCLIGC(4)=17 ! ALNIP - Near-IR albedo for direct radiation
        NCLIGC(5)=18 ! ALNID - Near-IR albedo for diffuse radiation
        NCLIGC(6)=210186  ! ALUVPI (or ALUVI) - Isotropic component of UV/Vis albedo
        NCLIGC(7)=210187  ! ALUVPV (or ALUVV) - Volumetric component of UV/Vis albedo
        NCLIGC(8)=210188  ! ALUVPG (or ALUVG) - Geometric component of UV/Vis albedo
        NCLIGC(9)=210189  ! ALNIPI (or ALNII) - Isotropic component of Near-IR albedo
        NCLIGC(10)=210190 ! ALNIPV (or ALNIV) - Volumetric component of Near-IR albedo
        NCLIGC(11)=210191 ! ALNIPG (or ALNIG) - Geometric component of Near-IR albedo
        NCLIGC(12)=139 
        NCLIGC(13)=31
    ENDIF
    IF((.NOT.LLMODISALBEDO) .AND. (LLELAIV) .AND. (.NOT.LECURR) ) THEN ! variable LAI  + No(MODIS albedo + ocean current)
        NCLIMTOT=5
        NCLIMR=3
        NCLIGC(1)=174
        NCLIGC(2)=139
        NCLIGC(3)=31
        NCLIGC(4)=66 ! LAI_LV - Leaf area index, low vegetation
        NCLIGC(5)=67 ! LAI_HV - Leaf area index, high vegetation
    ENDIF
    IF((LLMODISALBEDO) .AND. (.NOT.LLELAIV) .AND. (LECURR) ) THEN !No variable LAI  + (MODIS albedo + ocean current)
        NCLIMTOT=15
        NCLIMR=11
        NCLIGC(1)=174
        NCLIGC(2)=15
        NCLIGC(3)=16
        NCLIGC(4)=17
        NCLIGC(5)=18
        NCLIGC(6)=210186
        NCLIGC(7)=210187
        NCLIGC(8)=210188
        NCLIGC(9)=210189
        NCLIGC(10)=210190
        NCLIGC(11)=210191
        NCLIGC(12)=139
        NCLIGC(13)=31
        NCLIGC(14)=131
        NCLIGC(15)=132
    ENDIF
    IF((LLMODISALBEDO) .AND. (LLELAIV) .AND. (.NOT.LECURR) ) THEN !( variable LAI  + MODIS albedo) + No ocean current
        NCLIMTOT=15
        NCLIMR=13
        NCLIGC(1)=174
        NCLIGC(2)=15
        NCLIGC(3)=16
        NCLIGC(4)=17
        NCLIGC(5)=18
        NCLIGC(6)=210186
        NCLIGC(7)=210187
        NCLIGC(8)=210188
        NCLIGC(9)=210189
        NCLIGC(10)=210190
        NCLIGC(11)=210191
        NCLIGC(12)=66
        NCLIGC(13)=67
        NCLIGC(14)=139
        NCLIGC(15)=31
    ENDIF
    IF((.NOT.LLMODISALBEDO) .AND. (LLELAIV) .AND. (LECURR) ) THEN !(variable LAI  + ocean current) + No MODIS albedo 
        NCLIMTOT=7
        NCLIMR=3
        NCLIGC(1)=174
        NCLIGC(2)=66
        NCLIGC(3)=67
        NCLIGC(4)=139
        NCLIGC(5)=31
        NCLIGC(6)=131
        NCLIGC(7)=132
    ENDIF
    IF((LLMODISALBEDO) .AND. (LLELAIV) .AND. (LECURR) ) THEN ! variable LAI  + MODIS albedo +  ocean current
        NCLIMTOT=17
        NCLIMR=13
        NCLIGC(1)=174
        NCLIGC(2)=15
        NCLIGC(3)=16
        NCLIGC(4)=17
        NCLIGC(5)=18
        NCLIGC(6)=210186 
        NCLIGC(7)=210187 
        NCLIGC(8)=210188 
        NCLIGC(9)=210189 
        NCLIGC(10)=210190
        NCLIGC(11)=210191
        NCLIGC(12)=66
        NCLIGC(13)=67
        NCLIGC(14)=139
        NCLIGC(15)=31
        NCLIGC(16)=131
        NCLIGC(17)=132
    ENDIF
!
  ELSE  !OASIS COUPLER NOT ACTIVE IN ECMWF ROUTINES
    IF((.NOT.LLMODISALBEDO) .AND. (.NOT.LLELAIV)) THEN !NO( MODIS albedo + variable LAI)
        NCLIMR=3
        NCLIMTOT=NCLIMR
        NCLIGC(1)=31
        NCLIGC(2)=139
        NCLIGC(3)=174
    ENDIF
    IF((LLMODISALBEDO) .AND. (.NOT.LLELAIV)) THEN ! MODIS albedo + No(variable LAI)
        NCLIMR=13
        NCLIMTOT=NCLIMR
        NCLIGC(1)=31
        NCLIGC(2)=139
        NCLIGC(3)=174
        NCLIGC(4)=15
        NCLIGC(5)=16
        NCLIGC(6)=17
        NCLIGC(7)=18
        NCLIGC(8)=210186
        NCLIGC(9)=210187
        NCLIGC(10)=210188
        NCLIGC(11)=210189
        NCLIGC(12)=210190
        NCLIGC(13)=210191
    ENDIF
    IF((.NOT.LLMODISALBEDO) .AND. (LLELAIV)) THEN ! variable LAI  + No(MODIS albedo)
        NCLIMR=5
        NCLIMTOT=NCLIMR
        NCLIGC(1)=31
        NCLIGC(2)=139
        NCLIGC(3)=174
        NCLIGC(4)=66
        NCLIGC(5)=67
    ENDIF
    IF((LLMODISALBEDO) .AND. (LLELAIV)) THEN ! variable LAI  + MODIS albedo
        NCLIMR=15
        NCLIMTOT=NCLIMR
        NCLIGC(1)=31
        NCLIGC(2)=139
        NCLIGC(3)=174
        NCLIGC(4)=15
        NCLIGC(5)=16
        NCLIGC(6)=17
        NCLIGC(7)=18
        NCLIGC(8)=210186
        NCLIGC(9)=210187
        NCLIGC(10)=210188
        NCLIGC(11)=210189
        NCLIGC(12)=210190
        NCLIGC(13)=210191
        NCLIGC(14)=66
        NCLIGC(15)=67
    ENDIF

  ENDIF  !END OASIS COUPLER ACTIVE IN ECMWF ROUTINES

  IF (LMCCIEC) THEN !THE BOUNDARY CONDITIONS (SST) ARE INTERPOLATED IN TIME
    NP1=1
    NP2=2
    ITIM=2
  ELSE !THE BOUNDARY CONDITIONS (SST) ARE NOT INTERPOLATED IN TIME
    NP1=1
    NP2=NP1
    ITIM=1
  ENDIF !END THE BOUNDARY CONDITIONS (SST) ARE INTERPOLATED IN TIME

  IF(NCLIMR == 0) THEN
    NDIFC=-999
    NJDCR=-999
  ENDIF

  NJDCR1=-999

! Allocates the necessary space

  ALLOCATE(YDMCC%CLIMR(NGPTOT,ITIM,NCLIMTOT))
  YDMCC%CLIMR(:,:,:)=-999._JPRD
  IF(NPRINTLEV >= 1.OR. LALLOPR)&
   & WRITE(KULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)")&
   & 'CLIMR     ',SIZE(YDMCC%CLIMR   ),SHAPE(YDMCC%CLIMR    )  

  IF (LEOCLAKE) THEN
    ALLOCATE(YDMCC%ZLAKE(NGPTOT,ITIM,2))
  ENDIF

  IF (LEOBC .OR. LEOBCICE) THEN
    ALLOCATE(YDMCC%OCEANBC(NGPTOT,2))
    WRITE(NULOUT,*) 'OCEANBC ALLOCATED'
  ENDIF

! naj set up mechanism to read in emissions and dry deposition fields during model run for COMPO (AER AND CHEM) 
  IF (LMCC_COMPO) THEN
    IF (LMCCIEC_COMPO) THEN ! Emissions etc. interpolated in time
      NPCOMPO_1=1
      NPCOMPO_2=2
      ITIM=2
    ELSE ! Emissions etc. not interpolated in time
      NPCOMPO_1=1
      NPCOMPO_2=NP1
      ITIM=1
    ENDIF

    NCLIMR_COMPO=0
    IC=0      

    DO IT=1,NEMIS2D
      IF ( ANY( YEMIS2D_DESC(1:NEMIS2D_DESC)%PARAM_INDEX==IT .AND. &
              & YEMIS2D_DESC(1:NEMIS2D_DESC)%TEMPORALITY /= "Constant" ) ) THEN
        IF ( ANY( YEMIS2D_DESC(1:NEMIS2D_DESC)%PARAM_INDEX==IT .AND. &
                  & YEMIS2D_DESC(1:NEMIS2D_DESC)%TEMPORALITY /= "MCC" .AND. &
                  & YEMIS2D_DESC(1:NEMIS2D_DESC)%TEMPORALITY /= "Default" ) ) THEN
          CALL ABOR1('SUMCCLAG: EMIS2D only supports "Constant", "MCC" or "Default" TEMPORALITY, which must be consistent')
        ENDIF
        IC=IC+1
        NCLIGC_COMPO(IC) = YDSURF%YSD_VF%YEMIS2D(IT)%IGRBCODE
        NYSDMP_COMPO(IC) = YDSURF%YSD_VF%YEMIS2D(IT)%MP
      ENDIF
    ENDDO

    DO IT=1,NEMIS2DAUX
      SELECT CASE (YEMIS2DAUX_DESC(IT)%TEMPORALITY)
        CASE ('Constant')
        CASE ('MCC', 'Default')
          IC=IC+1
          NCLIGC_COMPO(IC) = YDSURF%YSD_VF%YEMIS2DAUX(IT)%IGRBCODE
          NYSDMP_COMPO(IC) = YDSURF%YSD_VF%YEMIS2DAUX(IT)%MP
        CASE DEFAULT
          CALL ABOR1('SUMCCLAG: EMIS2DAUX only supports "Constant", "MCC" or "Default" TEMPORALITY')
      END SELECT
    ENDDO

    ! What about 3D emissions??

    ! Chemistry-scheme prescribed deposition velocities only if LCHEM_DVOL=false
    IF ( NCHEM > 0 .AND. KCHEM_DRYDEP == 0 ) THEN 
      ICHEMFLX=0
      DO IT = 1, NCHEM 
        IF (YCHEM(IT)%IGRIBDV > 0 ) THEN
          IC=IC+1
          ICHEMFLX=ICHEMFLX+1
          NCLIGC_COMPO(IC) = YCHEM(IT)%IGRIBDV 
          NYSDMP_COMPO(IC) = YDSURF%YSD_VF%YCHEMDV(ICHEMFLX)%MP  
        ENDIF 
      ENDDO   
    ENDIF ! NCHEM && LCHEM_DVOL

    ! Aerosol-scheme DMS input only if LAERCHEM=false and we have DMS emissions
    !  - Should test LAEROSFC once emission refactoring is complete
    IF (NACTAERO > 0 .AND. .NOT. LAERCHEM .AND. NDMSO > 0 .AND. LOCNDMS) THEN
      IC=IC+1
      NCLIGC_COMPO(IC) = NGRBAERLS13
      NYSDMP_COMPO(IC) = YDSURF%YSD_VF%YDMSO%MP
    ENDIF

    ! Aerosol-scheme SO2 deposition only if LAERCHEM=false (and ideally only if we have sulphates)
    IF (NACTAERO > 0 .AND. .NOT. LAERCHEM) THEN ! Should test LAEROSFC once emission refactoring is complete
      IC=IC+1
      NCLIGC_COMPO(IC) = NGRBAERSO2DD
      NYSDMP_COMPO(IC) = YDSURF%YSD_VF%YSO2DD%MP
    ENDIF
 
    NCLIMR_COMPO = IC 

    IF(NCLIMR_COMPO == 0) THEN
      NDIFC_COMPO=-999
      NJDCR_COMPO=-999
    ENDIF

    NJDCR1_COMPO=-999

! allocate communication array 
    ALLOCATE(YDMCC%CLIMRCOMPO(NGPTOT,ITIM,NCLIMR_COMPO))
    
    WRITE(KULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)")&
    & 'CLIMRCOMPO     ',SIZE(YDMCC%CLIMRCOMPO),SHAPE(YDMCC%CLIMRCOMPO)

    WRITE(KULOUT,"(1X,A20)")  ' EXPECTED ORDER '
    DO IT=1,NCLIMR_COMPO 
      WRITE(KULOUT,"(1X,3I10)") IT,  NCLIGC_COMPO(IT),  NYSDMP_COMPO(IT)
    ENDDO

  ENDIF ! Composition part

ENDIF ! END BOUNDARY CONDITIONS UPDATED BY ECMWF ROUTINES

CALL SU_COUP_COM(YDML_AOC%YRCOM,KULOUT)

! Initialize NEMO coupling if active

IF (LNEMOCOUP) THEN
   CALL CPLNG_INIT(LDACTIVE=.TRUE.)
   YDMCC%CPLNG_ACTIVE=.TRUE.
ENDIF

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!      ----------------------------------------------------------------

!*       2.    PREPARES SIGNAL FOR OCEAN COUPLER.
!              ----------------------------------

! With OASIS, IOPROC must be 1.
IOPROC=1

IF((.NOT.LNEMOCOUP).AND.((LMCC03.AND.NOACOMM /= 5).OR.(LMCC04)))THEN
  ISTEP=INT(TSTEP)
  IF (NOACOMM == 5) THEN
     WRITE(KULOUT,'(A)')'OASIS3 run in SUMCCLAG'
  ELSE
  IF (LMCC03) THEN
      IF (MYPROC == IOPROC) CALL INICOU(YDML_AOC,YDEPHY,NSTOP,NFRCPL,ISTEP)
  ELSE
      CALL INICOU(YDML_AOC,YDEPHY,NSTOP,NFRCPL,ISTEP)
  ENDIF
  ENDIF
ENDIF

!      ----------------------------------------------------------------

!*       3.    PRINTINGS.
!              ----------

WRITE(UNIT=KULOUT,FMT='('' Printings done in SUMCCLAG '')')

WRITE(UNIT=KULOUT,FMT='('' NFRCPL = '',I6,'' NCLIMR = '',I6,'' NCLIMR_COMPO = '',I6)') NFRCPL,NCLIMR,NCLIMR_COMPO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUMCCLAG',1,ZHOOK_HANDLE)
END SUBROUTINE SUMCCLAG
