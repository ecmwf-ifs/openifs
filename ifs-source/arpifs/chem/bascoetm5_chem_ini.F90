! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOETM5_CHEM_INI(YGFL,YDCHEM,YDCOMPO)

!**   DESCRIPTION 
!     ----------
!
!   TM5 and BASCOE routines for IFS chemistry: Checking tracer indices,
!      and initialization of 
!      - chem-rates lookup-table
!      - tropo photolysis preparations
!      - stratospheric boundary conditions for CH4
!      - stratospheric HNO3 BC  
!      - preparation of computation of heterogeneous rates for n2o5 (outdated).
!
!
!**   INTERFACE.
!     ----------
!          *BASCOETM5_CHEM_INI* IS CALLED FROM *CHEM_INIT*.

! INPUTS: none
! -------
!
! OUTPUTS: none
! -------
!
!
!     AUTHOR
!     -------
!     2014-02-01: VINCENT HUIJNEN    *KNMI*

!     MODIFICATIONS.
!     --------------
!     2017-12-07: YVES CHRISTOPHE (YC) *BIRA*
!              strato photolysis lookup tables initialized from a TUV based radiative
!                   transfer code instead of read from exernal file
!     2018-03-27: YVES CHRISTOPHE (YC) *BIRA*
!              strato photolysis initialization is done from BASCOE_J_INI, which only
!                    prepares environment to run TUV (compute J rates) on line 
!     2018-09-14: YVES CHRISTOPHE (YC) *BIRA*
!               time and latband dependent stratospheric species surface boundary conditions
!
!     2020-02-07: JASON WILLIAMS    *KNMI*  
!                update to NOx, Isoprene, XYL, TOL, CH3CN 
!
!     2021-05-20 JONAS DEBOSSCHER (JD) *BIRA*
!                   Aerosol SAD climatology initialization



USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCHEM  , ONLY : TCHEM
USE YOMCOMPO , ONLY : TCOMPO
USE PARKIND1  ,ONLY : JPIM, JPRM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOETM5_MODULE, ONLY : NBC, BASCOE_BCVAL, BASCOE_BCNAME



IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

!INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KLEV


! * LOCAL 
TYPE(TYPE_GFLD),INTENT(INOUT)   :: YGFL
TYPE(TCHEM),    INTENT(IN)      :: YDCHEM
TYPE(TCOMPO),   INTENT(IN)      :: YDCOMPO
REAL(KIND=JPHOOK)                 :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "abor1.intfb.h"
#include "bascoe_j_ini.intfb.h"
#include "bascoe_lbc_ini.intfb.h"
#include "bascoe_setbin.intfb.h"
#include "bascoe_sage_init.intfb.h"
#include "bascoe_tropopause_init.intfb.h"
#include "bascoe_climSAD_ini.intfb.h"

! chemistry scheme name - this will later also come from external input

IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI',0,ZHOOK_HANDLE)

! Some checking
 CALL TRACER_IDX_CHECK_BASCOETM5(YGFL,YDCOMPO)
 
  ! ----------------------------------------------------
  ! Tropospheric chemistry
  ! ----------------------------------------------------

 !Check on settings for chemistry solver:
 ! Combination of KPP-chemistry together with 
 ! 'revised' chemistry version 'tc02b' is so far not supported
 IF (YDCHEM%LCHEM_REVCHEM) THEN
   CALL ABOR1('bascoetm5_chem_ini: LCHEM_REVCHEM=.true. not supported in combination with BASCOE')
 ENDIF
 
  ! Prepare / read in tropo photolysis data
 CALL PHOTOLYSIS_INI

 ! initialize table needed for budget evaluation
 CALL BASCOETM5_INI_BUDGET(YGFL,YDCHEM)

 ! Calculate look up tables rate constants/henry coefficients:
 CALL BASCOETM5_RATES

 ! climatological stratospheric boundary conditions for HNO3 and CH4
 ! No longer needed ! 
 ! CALL BOUNDARY_HNO3
 ! CALL BOUNDARY_CH4STRAT

  ! ----------------------------------------------------
  ! Stratospheric chemistry
  ! ----------------------------------------------------
 
 ! Prepare / read in photolysis table
 IF (YDCHEM%LCHEM_BASCOE_JON) THEN
   ! Online J-rate computation
   ! Prepare data to compute photolysis rates
   CALL BASCOE_J_INI
 ELSE
   ! Initialize J-lookup tables
   CALL BASCOE_J_TABLES_INIT
 ENDIF

  ! Prepare data with surface boundary conditions
 CALL BASCOE_LBC_INI(NBC, BASCOE_BCVAL, BASCOE_BCNAME)

  ! Prepare strat. aerosol particle size distribution 
 CALL BASCOE_SETBIN
 
 ! SAGE initialization
 CALL BASCOE_SAGE_INIT

 ! Tropopause pressure level initialization
 CALL BASCOE_TROPOPAUSE_INIT

 ! JD: Aerosol SAD climatology initialization
 CALL BASCOE_CLIMSAD_INI


IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI',1,ZHOOK_HANDLE)

CONTAINS

!       CHeck correct tracer indices
!       -----------------------------------------------------------
SUBROUTINE TRACER_IDX_CHECK_BASCOETM5(YGFL,YDCOMPO)


USE BASCOETM5_MODULE    , ONLY : & 
 &  IO3,     IH2O2,   ICH4,      ICO,       INO3_A,   IPSC,     IPB210,   IHNO3,     &
 &  ICH3O2H, ICH2O,   IPAR,      IETH,      IOLE,     IALD2,    IPAN,     IROOH,     &
 &  IORGNTR, IISOP,   ISO2,      IDMS,      INH3,     ISO4,     INH4,     IMSA,      &
 &  IMSA,    IMGLY,   IO3S,      IRN222,    INO,      IHO2,     ICH3O2,   ICH3O2NO2, & 
 &  IHONO,   ICH3,    ICH3O,     IHCO,      IOH,      INO2,     INO3,     IN2O5,     &
 &  IHO2NO2, IC2O3,   IROR,      IRXPAR,    IXO2,     IXO2N,    INH2,     ICH3OH,    &
 &  IHCOOH,  IMCOOH,  IC2H6,     IETHOH,    IC3H8,    IC3H6,    ITERP,    IISPD,     &
 &  IACET,   IACO2,   IHYPROPO2, IIC3H7O2,  IHCN,     ICH3CN,   ITOL,     IXYL,      &
 &  IAROO2,  IHPALD1, IHPALD2,   IISOPOOH,  IGLY,     IGLYALD,  IHYAC,    IISOPBO2,  &
 &  IISOPDO2,IN2O,    IH2O,      IOCLO,     IHCL,     ICLONO2,  IHOCL,    ICL2,      &
 &  IHBR,    IBRONO2, ICL2O2,    IHOBR,     IBRCL,    ICFC11,   ICFC12,   ICFC113,   &
 &  ICFC114, ICFC115, ICCL4,     ICLNO2,    ICH3CCL3, ICH3CL,   IHCFC22,  ICH3BR,    &
 &  IHF,     IHA1301, IHA1211,   ICHBR3,    ICLOO,    IO,       IO1D,     IN,        &
 &  ICLO,    ICL,     IBR,       IBRO,      IH,       IH2,      ICO2,     IBR2,      &
 &  ICH2BR2, ISTRATAER, ISO3,    IOCS,      IH2SO4,   ISOG1,    ISOG2A,   ISOG2B

  
 
USE YOMLUN             , ONLY : NULOUT
USE PARKIND1           , ONLY : JPIM  , JPRB
USE YOM_YGFL           , ONLY : TYPE_GFLD
USE YOMCOMPO           , ONLY : TCOMPO
USE YOMHOOK            ,ONLY  : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE

! Local parameters

! * counters
TYPE(TYPE_GFLD),INTENT(INOUT):: YGFL
TYPE(TCOMPO),   INTENT(IN)      :: YDCOMPO
INTEGER(KIND=JPIM) :: JL
LOGICAL            :: LLFOUND, LLFOUND_B, LLFOUND_C
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:TRACER_IDX_CHECK_BASCOETM5',0,ZHOOK_HANDLE)

ASSOCIATE(NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM)
   DO JL = 1,NCHEM  
      LLFOUND = .FALSE.
     SELECT CASE ( TRIM( YCHEM(JL)%CNAME ) ) 
! Tracers mainly of interest in troposphere:
       CASE ('O3')    ; LLFOUND = ( IO3 == JL )
       CASE ('H2O2')  ; LLFOUND = (IH2O2 == JL ) 
       CASE ('CH4')   ; LLFOUND = (ICH4 == JL )
       CASE ('CO')    ; LLFOUND = (ICO == JL )
       CASE ('HNO3')  ; LLFOUND = (IHNO3 == JL )
       CASE ('CH3OOH'); LLFOUND = (ICH3O2H == JL )
       CASE ('CH2O')  ; LLFOUND = (ICH2O == JL )
       CASE ('PAR')   ; LLFOUND = (IPAR == JL )
       CASE ('ETH','C2H4')   ; LLFOUND = (IETH == JL )
       CASE ('OLE')   ; LLFOUND = (IOLE == JL )
       CASE ('ALD2')  ; LLFOUND = (IALD2 == JL)
       CASE ('PAN')   ; LLFOUND = (IPAN == JL )
       CASE ('ROOH')  ; LLFOUND = (IROOH == JL)
       CASE ('ORGNTR','ONIT'); LLFOUND = (IORGNTR == JL )
       CASE ('CH3O2NO2')  ; LLFOUND = (ICH3O2NO2 == JL)
       CASE ('ISOP','C5H8')  ; LLFOUND = (IISOP == JL )
       CASE ('SO2')   ; LLFOUND = (ISO2 == JL )
       CASE ('SO3')   ; LLFOUND = (ISO3 == JL )
       CASE ('OCS')   ; LLFOUND = (IOCS == JL )
       CASE ('DMS')   ; LLFOUND = (IDMS == JL )
       CASE ('NH3')   ; LLFOUND = (INH3 == JL )
       CASE ('SO4')   ; LLFOUND = (ISO4 == JL )
       CASE ('NH4')   ; LLFOUND = (INH4 == JL )
       CASE ('MSA')   ; LLFOUND = (IMSA == JL )
       CASE ('MGLY','CH3COCHO')  ; LLFOUND = (IMGLY == JL )
       CASE ('O3S')   ; LLFOUND = (IO3S == JL )
       CASE ('Rn')    ; LLFOUND = (IRN222 == JL )
       CASE ('Pb')    ; LLFOUND = (IPB210 == JL )
       CASE ('NO')    ; LLFOUND = (INO == JL )
       CASE ('HO2')   ; LLFOUND = (IHO2 == JL )
       CASE ('CH3O2') ; LLFOUND = (ICH3O2 == JL )
       CASE ('OH')    ; LLFOUND = (IOH == JL )
       CASE ('NO2')   ; LLFOUND = (INO2 == JL )
       CASE ('NO3')   ; LLFOUND = (INO3 == JL )
       CASE ('HONO')   ; LLFOUND = (IHONO == JL )
       CASE ('N2O5')  ; LLFOUND = (IN2O5 == JL )
       CASE ('HNO4','HO2NO2')  ; LLFOUND = (IHO2NO2 == JL )
       CASE ('C2O3')  ; LLFOUND = (IC2O3 == JL )
       CASE ('ROR')   ; LLFOUND = (IROR == JL )
       CASE ('RXPAR') ; LLFOUND = (IRXPAR == JL )
       CASE ('XO2')   ; LLFOUND = (IXO2 == JL )
       CASE ('XO2N')  ; LLFOUND = (IXO2N == JL )
       CASE ('NH2')   ; LLFOUND = (INH2 == JL )
       CASE ('PSC')   ; LLFOUND = (IPSC == JL )
       CASE ('CH3OH')   ; LLFOUND = (ICH3OH == JL )
       CASE ('HCOOH')   ; LLFOUND = (IHCOOH == JL )
       CASE ('MCOOH')   ; LLFOUND = (IMCOOH == JL )
       CASE ('C2H6')   ; LLFOUND = (IC2H6 == JL )
       CASE ('ETHOH','C2H5OH')   ; LLFOUND = (IETHOH == JL )
       CASE ('C3H8')   ; LLFOUND = (IC3H8 == JL )
       CASE ('C3H6')   ; LLFOUND = (IC3H6 == JL )
       CASE ('TERP','C10H16')   ; LLFOUND = (ITERP == JL )
       CASE ('ISPD')   ; LLFOUND = (IISPD == JL )
       CASE ('NO3_A')  ; LLFOUND = (INO3_A == JL )
       CASE ('ACET','CH3COCH3')    ; LLFOUND = (IACET == JL )
       CASE ('ACO2')     ; LLFOUND = (IACO2 == JL )
       CASE ('HYPROPO2') ; LLFOUND = (IHYPROPO2 == JL )
       CASE ('IC3H7O2')  ; LLFOUND = (IIC3H7O2 == JL )
       CASE ('ISOPOOH')  ; LLFOUND = (IISOPOOH == JL )
       CASE ('HCN')      ; LLFOUND = (IHCN == JL )
       CASE ('CH3CN')    ; LLFOUND = (ICH3CN == JL )
       CASE ('XYL')      ; LLFOUND = (IXYL == JL )
       CASE ('TOL')      ; LLFOUND = (ITOL == JL )
       CASE ('AROO2')    ; LLFOUND = (IAROO2 == JL ) 
       CASE ('HPALD1')   ; LLFOUND = (IHPALD1 == JL )
       CASE ('HPALD2')   ; LLFOUND = (IHPALD2 == JL )
       CASE ('GLYALD')   ; LLFOUND = (IGLYALD == JL )
       CASE ('GLY')      ; LLFOUND = (IGLY == JL )
       CASE ('HYAC')     ; LLFOUND = (IHYAC == JL )              
       CASE ('ISOPBO2')  ; LLFOUND = (IISOPBO2 == JL )              
       CASE ('ISOPDO2')  ; LLFOUND = (IISOPDO2 == JL )              
       CASE ('SOG1')    ; LLFOUND = (iSOG1 == JL )
       CASE ('SOG2A')   ; LLFOUND = (iSOG2A == JL )
       CASE ('SOG2B')   ; LLFOUND = (iSOG2B == JL )
! Tracers mainly of interest in stratosphere:
       CASE ('H2SO4') ; LLFOUND = (IH2SO4 == JL )
       CASE ('CH3')   ; LLFOUND = (ICH3 == JL )
       CASE ('CH3O')  ; LLFOUND = (ICH3O == JL )
       CASE ('HCO')   ; LLFOUND = (IHCO == JL )
       CASE ('N2O')   ; LLFOUND = (IN2O == JL )
       CASE ('H2O')   ; LLFOUND = (IH2O == JL )
       CASE ('OCLO')  ; LLFOUND = (IOCLO == JL )
       CASE ('HCL')   ; LLFOUND = (IHCL == JL )
       CASE ('CLONO2'); LLFOUND = (ICLONO2 == JL )
       CASE ('HOCL')  ; LLFOUND = (IHOCL == JL )
       CASE ('CL2')   ; LLFOUND = (ICL2 == JL )
       CASE ('HBR')   ; LLFOUND = (IHBR == JL )
       CASE ('BRONO2'); LLFOUND = (IBRONO2 == JL )
       CASE ('CL2O2') ; LLFOUND = (ICL2O2 == JL )
       CASE ('HOBR')  ; LLFOUND = (IHOBR == JL )
       CASE ('BRCL')  ; LLFOUND = (IBRCL == JL )
       CASE ('CFC11') ; LLFOUND = (ICFC11 == JL )
       CASE ('CFC12') ; LLFOUND = (ICFC12 == JL )
       CASE ('CFC113'); LLFOUND = (ICFC113 == JL )
       CASE ('CFC114'); LLFOUND = (ICFC114 == JL )
       CASE ('CFC115'); LLFOUND = (ICFC115 == JL )
       CASE ('CCL4')  ; LLFOUND = (ICCL4 == JL )
       CASE ('CLNO2') ; LLFOUND = (ICLNO2 == JL )
       CASE ('CH3CCL3'); LLFOUND = (ICH3CCL3 == JL )
       CASE ('CH3CL') ; LLFOUND = (ICH3CL == JL )
       CASE ('HCFC22'); LLFOUND = (IHCFC22 == JL )
       CASE ('CH3BR') ; LLFOUND = (ICH3BR == JL )
       CASE ('HF')    ; LLFOUND = (IHF == JL )
       CASE ('HA1301'); LLFOUND = (IHA1301 == JL )
       CASE ('HA1211'); LLFOUND = (IHA1211 == JL )
       CASE ('CHBR3') ; LLFOUND = (ICHBR3 == JL )
       CASE ('CLOO')  ; LLFOUND = (ICLOO == JL )
       CASE ('O')     ; LLFOUND = (IO == JL )
       CASE ('O1D')   ; LLFOUND = (IO1D == JL )
       CASE ('N')     ; LLFOUND = (IN == JL )
       CASE ('CLO')   ; LLFOUND = (ICLO == JL )
       CASE ('CL')    ; LLFOUND = (ICL == JL )
       CASE ('BR')    ; LLFOUND = (IBR == JL )
       CASE ('BRO')   ; LLFOUND = (IBRO == JL )
       CASE ('H')     ; LLFOUND = (IH == JL )
       CASE ('H2')    ; LLFOUND = (IH2 == JL )
       CASE ('CO2')   ; LLFOUND = (ICO2 == JL )
       CASE ('BR2')   ; LLFOUND = (IBR2 == JL )
       CASE ('CH2BR2'); LLFOUND = (ICH2BR2 == JL )
       CASE ('STRATAER')   ; LLFOUND = (ISTRATAER == JL )
! Additional tracers
       CASE ('NOXA')     ; LLFOUND = .TRUE. 
       CASE ('CLXA')     ; LLFOUND = .TRUE. 
       CASE ('BRXA')     ; LLFOUND = .TRUE.
       CASE ('CO_A_50')  ; LLFOUND = .TRUE. 
       CASE ('CO_A_25')  ; LLFOUND = .TRUE. 
       CASE ('PM10')     ; LLFOUND = .TRUE. 
       CASE ('PM25')     ; LLFOUND = .TRUE. 
       CASE ('VSO2')     ; LLFOUND = .TRUE.
       CASE DEFAULT
         WRITE(NULOUT,*) 'ERROR bascoetm5_chem_ini: no matching tracer name for '//TRIM(YCHEM(JL)%CNAME)
         CALL ABOR1('bascoetm5_chem_ini: No matching tracer name available')
     END SELECT
 
     IF (.NOT. LLFOUND ) THEN 
       WRITE(NULOUT,*) 'ERROR bascoetm5_chem_ini: Wrong tracer index or status for  '//TRIM(YCHEM(JL)%CNAME)
       CALL ABOR1('bascoetm5_chem_ini: wrong tracer index tracer name')
     ENDIF   
  ENDDO

! Also various checks in the inverse: tracers present in input table file  
! Particularly, this version also requires H2SO4, SO3 and OCS available in the list of trace gases
  LLFOUND   = .FALSE.
  DO JL = 1, NCHEM
    IF (TRIM(YCHEM(JL)%CNAME) == 'OCS' )    LLFOUND  = .TRUE.
  ENDDO
  IF (.NOT. LLFOUND) THEN
     CALL ABOR1('ERROR bascoetm5_chem_ini: OCS not defined in table-file')
  ENDIF   
  LLFOUND   = .FALSE.
  DO JL = 1, NCHEM
    IF (TRIM(YCHEM(JL)%CNAME) == 'SO3' )    LLFOUND  = .TRUE.
  ENDDO
  IF (.NOT. LLFOUND) THEN
     CALL ABOR1('ERROR bascoetm5_chem_ini: SO3 not defined in table-file')
  ENDIF   
  LLFOUND   = .FALSE.
  DO JL = 1, NCHEM
    IF (TRIM(YCHEM(JL)%CNAME) == 'H2SO4' )    LLFOUND  = .TRUE.
  ENDDO
  IF (.NOT. LLFOUND) THEN
     CALL ABOR1('ERROR bascoetm5_chem_ini: H2SO4 not defined in table-file')
  ENDIF   


! SOA chemistry also requires new XYL and TOL tracers
  LLFOUND   = .FALSE.
  LLFOUND_B = .FALSE.
  DO JL = 1, NCHEM
    IF (TRIM(YCHEM(JL)%CNAME) == 'XYL' )    LLFOUND  = .TRUE.
    IF (TRIM(YCHEM(JL)%CNAME) == 'TOL'  )    LLFOUND_B= .TRUE.
  ENDDO
  
  IF (.NOT. LLFOUND .OR. .NOT. LLFOUND_B) THEN
     CALL ABOR1('ERROR tm5_chem_ini: XYL and/or TOL not defined in table-file')
  ENDIF   
IF (YDCOMPO%LAERSOA .AND. YDCOMPO%LAERSOA_COUPLED) THEN
  LLFOUND   = .FALSE.
  LLFOUND_B = .FALSE.
  LLFOUND_C = .FALSE.
  DO JL = 1, NCHEM
    IF (TRIM(YCHEM(JL)%CNAME) == 'SOG1' )    LLFOUND  = .TRUE.
    IF (TRIM(YCHEM(JL)%CNAME) == 'SOG2A' )   LLFOUND_B= .TRUE.
    IF (TRIM(YCHEM(JL)%CNAME) == 'SOG2B'  )  LLFOUND_C= .TRUE.
  ENDDO
  
  IF (.NOT. LLFOUND .OR. .NOT. LLFOUND_B) THEN
     CALL ABOR1('ERROR tm5_chem_ini: SOG1 and/or SOG2A not defined in table-file, while LAERSOA=true')
  ENDIF   
  IF (.NOT. LLFOUND_C ) THEN
     CALL ABOR1('ERROR tm5_chem_ini: SOG2B not defined in table-file, while LAERSOA=true')
  ENDIF   
ENDIF


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:TRACER_IDX_CHECK_BASCOETM5',1,ZHOOK_HANDLE)
END SUBROUTINE TRACER_IDX_CHECK_BASCOETM5


!       Initialize new photolysis routine
!       -----------------------------------------------------------
SUBROUTINE PHOTOLYSIS_INI
USE TM5_PHOTOLYSIS , ONLY : WAVE_FULL, DWAVE_FULL, WAVE, DWAVE, MAXWAV,  &
     &   NWAV_NO2,NWAV_CO,NWAV_HNO3,NWAV_O3,NWAV_PAN,NWAV_H2O2,NWAV_N2O5,NWAV_CH2O,&
     &   NWAV_NO3,NWAV_ACET, &
     &   CS_RAY, PHI0, RD, FLUX, NTEMP_MAX, &
     &   XS_O3_LOOK,XS_HNO3_LOOK,XS_PAN_LOOK,QY_O3_LOOK,XS_H2O2_LOOK, &
     &   XS_N2O5_LOOK,XS_NO2_LOOK,QY_NO2_LOOK,QY_CO_LOOK,XS_NO3_LOOK, XS_CH2O_LOOK, &
     &   A1_ACET,A2_ACET,A3_ACET,A4_ACET, &
     &   SCA,ABS_EFF, GFAC, EXT , PN_REF 
  
USE YOMLUN             , ONLY : NULOUT
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

! Local parameters
INTEGER(KIND=JPIM), PARAMETER  :: IUTMP=77
INTEGER(KIND=JPIM)             :: IK, IWAV, JL
REAL(KIND=JPRB)                :: ZWL,ZX
LOGICAL                        :: LLFILEOK
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE

!TM5 photolysis info file
CHARACTER(LEN=27) :: CL_TM5_PHOTO_LOOKUP='tropo_look_up_cbmhybrid.dat'
!TM5 extraterrestrial flux:
CHARACTER(LEN=21) :: CL_TM5_EXTRATERRESTRIAL='OMI.data.extraterrest'
!TM5 aerosol:
CHARACTER(LEN=18) :: CL_TM5_AEROSOL='aerosol_reduce.dat'
#include "abor1.intfb.h"


IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:PHOTOLYSIS_INI',0,ZHOOK_HANDLE)

DO IK=1,13 
  WAVE_FULL(IK)=1._JPRB/(56250._JPRB-500._JPRB*IK) 
ENDDO 
          
DO IK=14,45 
 WAVE_FULL(IK)=1._JPRB/(49750._JPRB-(IK-13)*500._JPRB) 
ENDDO 

DO IK=46,68 
  WAVE_FULL(IK)=(266._JPRB+(IK-13))*1.E-7_JPRB
ENDDO 

DO IK=69,71 
  WAVE_FULL(IK)=(320.5_JPRB+2.*(IK-68))*1.E-7_JPRB 
ENDDO 
DO IK=72,176 
  WAVE_FULL(IK)=(325._JPRB+5._JPRB*(IK-71))*1.E-7_JPRB 
ENDDO 
   
DO IK=2,MAXWAV-1 
   DWAVE_FULL(IK)=0.5_JPRB*(WAVE_FULL(IK+1)-WAVE_FULL(IK-1))
ENDDO                 
DWAVE_FULL(1)    = DWAVE_FULL(2) 
DWAVE_FULL(MAXWAV) = DWAVE_FULL(MAXWAV-1) 
! select a subset from the full spctral grid to remove the first band      
WAVE=WAVE_FULL(14:135)
DWAVE=DWAVE_FULL(14:135)


!---------------------------------------------------------------------  
!     Rayleigh scattering cross sections                                
!     emp. formula of Nicolet   plan. space sci.32, 1467f (1984)        
!---------------------------------------------------------------------  
                                                                        
DO  IK=1,MAXWAV           
  ZWL=WAVE(IK)*1.E4_JPRB      
  ZX=0.389*ZWL+0.09426_JPRB/ZWL-0.3228_JPRB      
  CS_RAY(IK)=4.02E-28_JPRB/ZWL**(4._JPRB+ZX)     
ENDDO

!---------------------------------------------------------------------
!        QUANTUM YIELD OF METHYGLYOXAL S. KOCH and G. K. MOORTGAT       
!        J Phys Chem, 102, 9142-9153, 1998.
!
! JEW: the qy_ch3cocho array has been reduced to 52 bins to represent the pressure.
!      dep. spectral bins ONLY !!   
!---------------------------------------------------------------------

! Calculate co-efficients for METHYGLYOXAL qy outside loop
PHI0 = 1. 
RD = 0.      

! for wave < 380nm qy is essentially = 1.0.
DO IK = 1,39
  IF(IK<32) THEN
     PHI0(IK) = 1
  ELSEIF(IK>=32) THEN
     PHI0(IK) = 8.15E-9_JPRB * EXP(7131._JPRB/WAVE(IK+37_JPIM)*1.E-7_JPRB)
  ENDIF
  RD(IK)   = 7.34E-9_JPRB * EXP(8793._JPRB/WAVE(IK+37_JPIM)*1.E-7_JPRB)
ENDDO  
      
DO IK = 40,52
  PHI0(IK) = 3.63E-7_JPRB * EXP(5693._JPRB/WAVE(IK+37_JPIM)*1.E-7_JPRB) 
  RD(IK)   = 1.93E4_JPRB*EXP(-5639._JPRB/WAVE(IK+37_JPIM)*1.E-7_JPRB)
ENDDO  
   
!--------------------------------------------------------------         
!    Read and store temperature dependent cross-section values
!--------------------------------------------------------------         
LLFILEOK=.FALSE.
INQUIRE(FILE=TRIM(CL_TM5_PHOTO_LOOKUP),EXIST=LLFILEOK )
IF(.NOT. LLFILEOK)THEN
  WRITE(NULOUT,*) 'ERROR bascoetm5_chem_ini: Cannot read TM5 file '//TRIM(CL_TM5_PHOTO_LOOKUP)
  CALL ABOR1('tm5_chem_ini: photo inputfile does not exist.')
ELSE
  WRITE(NULOUT,*) 'INFO bascoetm5_chem_ini: Opening TM5 file '//TRIM(CL_TM5_PHOTO_LOOKUP) 
  OPEN(IUTMP,FORM='FORMATTED',ACTION='READ',STATUS='OLD',   &
   & FILE=TRIM(CL_TM5_PHOTO_LOOKUP))          

READ(IUTMP,*)
READ(IUTMP,*)  
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*) (XS_O3_LOOK(IWAV,IK),IWAV=1,MAXWAV)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*) (XS_NO2_LOOK(IWAV,IK),IWAV=1,NWAV_NO2)
ENDDO 
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*) (XS_HNO3_LOOK(IWAV,IK),IWAV=1,NWAV_HNO3)
ENDDO 
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*) (XS_H2O2_LOOK(IWAV,IK),IWAV=1,NWAV_H2O2)
ENDDO 
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
 READ(IUTMP,*) (XS_N2O5_LOOK(IWAV,IK),IWAV=1,NWAV_N2O5)
ENDDO 
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
 READ(IUTMP,*) (XS_CH2O_LOOK(IWAV,IK),IWAV=1,NWAV_CH2O)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
 READ(IUTMP,*) (XS_PAN_LOOK(IWAV,IK),IWAV=1,NWAV_PAN)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
 READ(IUTMP,*) (XS_NO3_LOOK(IWAV,IK),IWAV=1,NWAV_NO3)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
 READ(IUTMP,*) (QY_O3_LOOK(IWAV,IK),IWAV=1,NWAV_O3)
ENDDO 
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*)(QY_NO2_LOOK(IWAV,IK),IWAV=1,NWAV_NO2)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*)(QY_CO_LOOK(IWAV,IK),IWAV=1,NWAV_CO)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*)(A1_ACET(IWAV,IK),IWAV=1,NWAV_ACET)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*)(A2_ACET(IWAV,IK),IWAV=1,NWAV_ACET)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*)(A3_ACET(IWAV,IK),IWAV=1,NWAV_ACET)
ENDDO
READ(IUTMP,*)
DO IK=1,NTEMP_MAX
  READ(IUTMP,*)(A4_ACET(IWAV,IK),IWAV=1,NWAV_ACET)
ENDDO      
CLOSE(IUTMP)
 
ENDIF
      
! put constants in correct units
A1_ACET=A1_ACET*1E-18_JPRB
A2_ACET=A2_ACET*1E-17_JPRB
A4_ACET=A4_ACET*1E-15_JPRB      

!--------------------------------------------------------------         
!    Read and store temperature dependent cross-section values
!    
!----- extraterrestrial flux from input file ------------------
LLFILEOK=.FALSE.
INQUIRE(FILE=TRIM(CL_TM5_EXTRATERRESTRIAL),EXIST=LLFILEOK )
IF(.NOT. LLFILEOK)THEN
  WRITE(NULOUT,*) 'ERROR bascoectm5_chem_ini: Cannot read TM5 file '//TRIM(CL_TM5_EXTRATERRESTRIAL)
  CALL ABOR1('tm5_chem_ini: extraterrestrial inputfile does not exist.')
ELSE
  WRITE(NULOUT,*) 'INFO bascoetm5_chem_ini: Opening TM5 file '//TRIM(CL_TM5_EXTRATERRESTRIAL) 
  OPEN(IUTMP,FORM='FORMATTED',ACTION='READ',STATUS='OLD',   &
   & FILE=TRIM(CL_TM5_EXTRATERRESTRIAL))          
  READ(IUTMP,*)            
  READ(IUTMP,*) FLUX
  CLOSE(IUTMP)
ENDIF
    
!************************************************************************** 
!     aerosol data from: "Models for Aerosols of the Lower Atmosphere   
!     and the Effects of Humidity Variations on their Optical Properties
!     E. P. Shettle and R. W. Fenn (1979), Environmental Research Paper 
!     No. 676        
!************************************************************************** 
LLFILEOK=.FALSE.
INQUIRE(FILE=TRIM(CL_TM5_AEROSOL),EXIST=LLFILEOK )
IF(.NOT. LLFILEOK)THEN
  WRITE(NULOUT,*) 'ERROR bascoetm5_chem_ini: Cannot read TM5 file '//TRIM(CL_TM5_AEROSOL)
  CALL ABOR1('tm5_chem_ini: aerosol inputfile does not exist.')
ELSE
  WRITE(NULOUT,*) 'INFO bascoetm5_chem_ini: Opening TM5 file '//TRIM(CL_TM5_AEROSOL) 
  OPEN(IUTMP,FORM='FORMATTED',ACTION='READ',STATUS='OLD',   &
   & FILE=TRIM(CL_TM5_AEROSOL))          
  DO JL = 1,4
    READ(IUTMP,*)
    DO IK = 1,8
      READ(IUTMP,*)
      READ(IUTMP,'(6(1X,F7.5))')(EXT(IWAV,IK,JL),IWAV=1,MAXWAV) 
      READ(IUTMP,'(6(1X,F7.5))')(ABS_EFF(IWAV,IK,JL),IWAV=1,MAXWAV)
      READ(IUTMP,'(6(1X,F7.5))')(GFAC(IWAV,IK,JL),IWAV=1,MAXWAV)
      DO IWAV=1,MAXWAV
        SCA(IWAV,IK,JL) = EXT(IWAV,IK,JL) - ABS_EFF(IWAV,IK,JL)
      ENDDO
    ENDDO
  ENDDO
  CLOSE(IUTMP) 

  DO IWAV=1,MAXWAV
    DO IK=1,8
      DO JL=1,4
        SCA(IWAV,IK,JL) = SCA(IWAV,IK,JL)/PN_REF(JL)*1.E-5_JPRB
        ABS_EFF(IWAV,IK,JL) = ABS_EFF(IWAV,IK,JL)/PN_REF(JL)*1.E-5_JPRB
      ENDDO 
    ENDDO 
  ENDDO      

ENDIF  

!VH - file format(6(1X,F7.5))

WRITE(NULOUT,*) 'INFO bascoetm5_chem_ini: Finished reading TM5 photolysis input files ' 


IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:PHOTOLYSIS_INI',1,ZHOOK_HANDLE)


END SUBROUTINE PHOTOLYSIS_INI




SUBROUTINE BASCOETM5_INI_BUDGET(YGFL,YDCHEM)

!**   DESCRIPTION 
!     ----------
!
!    routine for initializing the chemical budget information table for trop. chem.
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOETM5_INI_BUDGET* IS CALLED FROM *BASCOETM5_CHEM_INI*.
!
USE TM5_CHEM_MODULE    , ONLY : NRR, NRJ, NREAC


USE BASCOETM5_MODULE    , ONLY : &
 &  IACID,   IAIR,    IH2O,      IO3,       IH2O2,    ICH4,     ICO,     &
 &  IHNO3,   ICH3O2H, ICH2O,     IPAR,      IETH,     IOLE,     IALD2,   &
 &  IPAN,    IROOH,   ICH3O2NO2, IORGNTR,   IISOP,    ISO2,     IDMS,    &
 &  INH3,    IMGLY,   IRN222,    INO,       IHO2,     IHONO,    ICH3O2,  &
 &  IOH,     INO2,    INO3,      IN2O5,     IHO2NO2,  IC2O3,    IROR,    &
 &  IRXPAR,  IXO2,    IXO2N,     INH2,      ICH3OH,   IHCOOH,   IMCOOH,  &
 &  IC2H6,   IETHOH,  IC3H8,     IC3H6,     ITERP,    IISPD,    IACET,   &
 &  IHCN,    ICH3CN,  IACO2,     IHYPROPO2, IIC3H7O2, IO2,      IISOPOOH,&
 &  IHPALD1, IHPALD2, IGLY,      IGLYALD,   IHYAC,    IISOPBO2, IISOPDO2,&
 &  ISO3,    ITOL,    IXYL,      IAROO2, ISOG2A


USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCHEM  , ONLY : TCHEM

IMPLICIT NONE
! * LOCAL 
TYPE(TYPE_GFLD),INTENT(INOUT)   :: YGFL
TYPE(TCHEM),    INTENT(IN)      :: YDCHEM
REAL(KIND=JPHOOK)                 :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:BASCOETM5_INI_BUDGET',0,ZHOOK_HANDLE)
ASSOCIATE(NCHEM=>YGFL%NCHEM)
IACID = NCHEM+1_JPIM
IAIR  = NCHEM+2_JPIM


  NRR = RESHAPE((/& 
  &         INO     , IHO2     , ICH3O2    , INO2     , IHNO3   , INO2     , INO      , INO2     , IN2O5    , IHONO     , &
  &         INO     , IHO2NO2  , INO2      , IHO2NO2  , ICH3O2  , ICH3O2NO2, IAIR     , IH2O     , IO3      , ICO       , &
  &         IO3     , IH2O2    , ICH2O     , ICH4     , ICH3O2H , IROOH    , ICH3O2   , ICH3O2   , ICH3O2   , IHO2      , &
  &         IHO2    , IN2O5    , IN2O5     , IOH      , ICH2O   , IALD2    , IALD2    , IC2O3    , IC2O3    , IPAN      , &
  &         IC2O3   , IC2O3    , IC2O3     , IPAR     , IROR    , IROR     , IOLE     , IO3      , INO3     , IETH      , &
  &         IO3     , IMGLY    , IISOP     , IO3      , INO3    , IXO2     , IXO2     , IXO2N    , IXO2     , IRXPAR    , &
  &         IORGNTR , IXO2N    , IDMS      , IDMS     , IDMS    , ISO2     , INH3     , INH3     , INH2     , INH2      , &
  &         INH2    , INH2     , INH2      , ICH3OH   , IHCOOH  , INO3     , INO3     , INO3     , INO3     , IMCOOH    , &
  &         IC2H6   , IETHOH   , IC3H8     , IC3H6    , IO3     , INO3     , ITERP    , ITERP    , ITERP    , IISPD     , &
  &         IISPD   , IISPD    , IRN222    , IO3      , IAIR    , IACET    , IACO2    , IACO2    , IACO2    , IACO2     , &
  &         IXO2    , IXO2N    , IN2O5     , IHO2     , INO3    , IHO2     , IIC3H7O2 , IIC3H7O2 , IHYPROPO2, IHYPROPO2 , &
  &         ISO2    , ISO3     , 1         , IOH      , IOH     , IOH      , IO3      , INO3     , IOH      , IO3       , &
  &         INO3    , INO      , IAROO2    , IHO2     , IXO2    , IISOPBO2 , IISOPBO2 , IISOPDO2 , IISOPDO2 , IISOPBO2  , &
  &         IISOPBO2, IISOPDO2 , IISOPDO2  , IISOPOOH , IHPALD1 , IHPALD2  , IGLY     , IGLYALD  , IHYAC    , ISOG2A, & 
         !second reaction partner (if monmolecular = 0)
  &         IO3     ,    INO   ,    INO    ,     IOH  ,     IOH ,      IO3 , INO3     ,  INO3    ,      0  ,  IOH  , &
  &         IOH     ,    IOH   ,   IHO2    ,       0  ,    INO2 ,        0 ,    0     ,     0    ,   IHO2  ,  IOH  , &
  &         IOH     ,    IOH   ,    IOH    ,     IOH  ,     IOH ,      IOH , IHO2     ,  IHO2    , ICH3O2  ,  IOH  , &
  &         IHO2    ,      0   ,      0    ,       0  ,    INO3 ,      IOH , INO3     ,   INO    ,   INO2  ,    0  , &
  &         IC2O3   ,   IHO2   ,    IHO2   ,     IOH  ,       0 ,        0 , IOH      ,  IOLE    ,   IOLE  ,  IOH  , &
  &         IETH    ,    IOH   ,    IOH    ,   IISOP  ,   IISOP ,      INO , IXO2     ,   INO    ,   IHO2  , IPAR  , &
  &         IOH     ,   IHO2   ,    IOH    ,     IOH  ,    INO3 ,      IOH , IACID    ,   IOH    ,   IOH   ,  INO  , &
  &         INO2    ,   IHO2   ,      0    ,     IOH  ,     IOH ,     IHO2 , ICH3O2   , IC2O3    ,   IXO2  ,  IOH  , &
  &         IOH     ,    IOH   ,     IOH   ,     IOH  ,   IC3H6 ,    IC3H6 , IOH      ,   IO3    ,   INO3  ,  IOH  , &
  &         IO3     ,   INO3   ,       0   ,    IAIR  ,       0 ,      IOH , IHO2     , ICH3O2   ,    INO  ,  IXO2 , &
  &         IXO2N   ,  IXO2N   ,       0   ,       0  ,       0 ,        0 , INO      ,  IHO2    ,    INO  ,  IHO2 , &
  &         IO3     ,      0   ,     IOH   ,    IHCN  ,  ICH3CN ,     ITOL , ITOL     ,  ITOL    ,    IXYL ,  IXYL , &
  &         IXYL    , IAROO2   ,  IAROO2   ,  IAROO2  ,  IAROO2 ,        0 ,    0     ,     0    ,      0  ,  IHO2 , &
  &         INO     ,   IHO2   ,     INO   ,     IOH  ,     IOH ,      IOH ,  IOH     ,    IOH   ,     IOH ,  IOH/),(/NREAC,2/))
  
  NRJ=(/ &
  &         IO3     , INO2    , IH2O2   , IHNO3  ,   IHO2NO2,  IN2O5, ICH2O, ICH2O,   ICH3O2H, INO3, &
  &         INO3    , IPAN    , IPAN    , IORGNTR, ICH3O2NO2,  IALD2, IMGLY, IROOH,   IO2,     IISPD, &
  &         IACET   , IACET   , IISOPOOH, IHONO  ,   IGLYALD,   IGLY,  IGLY, IHPALD1, IHPALD2, IHYAC/)  

  
  
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:BASCOETM5_INI_BUDGET',1,ZHOOK_HANDLE)
END SUBROUTINE BASCOETM5_INI_BUDGET




SUBROUTINE BASCOETM5_RATES
!----------------------------------------------------------------------
!     
!**   DESCRIPTION 
!     ----------
!      calculation of look up tables for rate constants and
!      Henry's law constants
!
!
!**   INTERFACE.
!     ----------
!          *RATES* IS CALLED FROM *BASCOETM5_CHEM_INI*.
!
!      method
!      ------
!      use known array of temperatures (integers) to calculate rate constants
!      3 body reactions are explicitly calculated in chemistry
!      
!      external
!      --------
!      none
!      
!      reference
!      ---------
!      Williams et al., ACP 2013
!
!------------------------------------------------------------------
USE TM5_CHEM_MODULE, ONLY : &  
  &  KNOO3,      KHO2NO,    KMO2NO,      KNO2O3,       KNONO3,     KN2O5A,     KN2O5B,     KHNO4OH,    KODM,       KH2OOD,    &
  &  KHNO4A,     KHNO4B,    &
  &  KO3HO2,     KO3OH,     KHPOH,       KFRMOH,       KCH4OH,     KOHMPER,    KOHROOH,    KMO2HO2A,   KMO2HO2B,   KMO2MO2,   &
  &  KHO2OH,     KH2OH,     KC41,        KC46,         KC49,       KC50A,      KC50B,      KC52,       KC53,       KC54,      &
  &  KC57,       KC58,      KC59,        KC62,         KC73,       KC76,       KC77,       KC78,       KC79,       KC80,      &
  &  KC81,       KC82,      KC83,        KC84,         KC85,       KDMSOHA,    KDMSOHB,    KDMSNO3,    KNH3OH,     KNH2OH,    &
  &  KNH2NO,     KNH2NO2,   KNH2HO2,     KNH2O2,       KOHCH3OH,   KOHHCOOH,   KOHMCOOH,   KOHC2H6,    KOHETHOH,   KOHC3H8,   &
  &  KO3C3H6,    KNO3C3H6,  KOHTERP,     KO3TERP,      KNO3TERP,   KNO2OHA,    KNO2OHB,    KOHHNO3A,   KOHHNO3B,   KOHHNO3C,  &
  &  KNO2NO3A,   KNO2NO3B,  KNO2HO2A,    KNO2HO2B,     KHO2HO2A,   KHO2HO2B,   KHO2HO2C,   KC47A,      KC47B,      KC48,      &
  &  KC61A,      KC61B,     KSO2OHA,     KSO2OHB,      KDMSOHC,    KCOOHA,     KCOOHB,     KCOOHC,     KCOOHD,     KOHMCHO,   &
  &  KOHMCH2CHO, KNO3MCHO,  KNO3MCH2CHO, KOHMVK,       KOHOLE,     KOHMACR,    KO3MVK,     KO3OLE,     KO3MACR,    KNO3OLE,   &
  &  KOHC3H6A,   KOHC3H6B,  KOHACETA,    KOHACETB,     KACO2HO2,   KACO2MO2,   KACO2NO,    KACO2XO2,   KXO2XO2N,   KXO2N,     &
  &  KO3PO2,     KO3PO3,    KNOHYPROPO2, KHO2HYPROPO2, KNOIC3H7O2, KHO2IC3H7O2,KSO2O3G,    KSO3H2O,    KOCSOH,     KOHHCN,    &
  &  KOHCH3CN,   KOHTOL,    KO3TOL,      KNO3TOL,      KOHXYL,     KCXYLO3A,   KCXYLO3B,   KNO3XYL,    KAROO2NO,   K2AROO2,   &
  &  KAROO2HO2,  KAROO2XO2, KISOPBO2A,   KISOPBO2B,    KISOPDO2A,  KISOPDO2B,  KISOPBO2HO2,KISOPBO2NO, KISOPDO2HO2,KISOPDO2NO,&
  &  KHPALD1OH, KHPALD2OH,   KGLYOH,       KGLYALDOH,  KHYACOH,    KISOPOOHOH, RATES_LUT,  NTEMP,      NTLOW, &
  &  KNO3HO2,KOHHONO,KHONOA, KHONOB, KMENO2A,KMENO2B, KMENO2M
USE BASCOETM5_MODULE, ONLY :   IH2O2,    IHNO3,  ICH3O2H, ICH2O,   IROOH,  IORGNTR, &
  &  ISO4,   INH4, IMSA,  ISO2, INH3,    IO3,    IMGLY,   IALD2,   IHCOOH, ICH3OH,  &
  &  IMCOOH, IETHOH, IO3S, IACET, IHCN,  ICH3CN, IHPALD1, IHPALD2, IHYAC,  IGLY,  HENRY
  
  
USE PARKIND1  ,ONLY : JPIM,JPRB, JPRM, JPRD
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

! local
INTEGER(KIND=JPIM) ::  ITEMP,IK
REAL(KIND=JPRD)    ::  ZRX1
REAL(KIND=JPRB)    ::  ZTREC,ZT3REC,ZTEMP,ZKH1,ZKH2,ZR
REAL(KIND=JPHOOK)    ::  ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:BASCOETM5_RATES',0,ZHOOK_HANDLE)
! start
! JEW : updated rates JPL(2006), incl. Evaluation Number 15 (March, 2006)
!    
DO IK=1,NTEMP
   ITEMP=IK+NTLOW
   ZTEMP=FLOAT(ITEMP)
   ZTREC=1./FLOAT(ITEMP)
   ZT3REC=300./FLOAT(ITEMP)
   !JEW: changed to JPL2006 
   RATES_LUT(KNOO3,IK)=ZFARR(3.E-12_JPRB,-1500._JPRB,ZTREC)
   !JEW: changed to JPL2006   
   RATES_LUT(KHO2NO,IK)=ZFARR(3.3E-12_JPRB,270._JPRB,ZTREC)
   !JEW: changed to JPL2006  
   ! RATES_LUT(KNO2OHA,IK)=1.8E-30_JPRB*ZT3REC**3.0_JPRB
   ! RATES_LUT(KNO2OHB,IK)=2.8E-11_JPRB
   RATES_LUT(KNO2OHA,IK)=3.2E-30_JPRB*ZT3REC**4.5_JPRB            ! IUPAC 2017 
   RATES_LUT(KNO2OHB,IK)=3.0E-11_JPRB
   RATES_LUT(KOHHNO3A,IK)=ZFARR(2.41E-14_JPRB,460._JPRB,ZTREC)    !!!wp!!! new ravi
   RATES_LUT(KOHHNO3B,IK)=ZFARR(6.51E-34_JPRB,1335._JPRB,ZTREC)   !!!wp!!! new ravi
   RATES_LUT(KOHHNO3C,IK)=ZFARR(2.69E-17_JPRB,2199._JPRB,ZTREC)   !!!wp!!! new ravi
   RATES_LUT(KNO2O3,IK)=ZFARR(1.4E-13_JPRB,-2470._JPRB,ZTREC)
   RATES_LUT(KNONO3,IK)=ZFARR(1.8E-11_JPRB,110._JPRB,ZTREC)
   RATES_LUT(KNO3HO2,IK)=4.0E-12_JPRB
   !JEW: changed to JPL2006  
   !RATES_LUT(KNO2NO3A,IK)=2.0E-30_JPRB*ZT3REC**4.4_JPRB
   !RATES_LUT(KNO2NO3B,IK)=1.4E-12_JPRB*ZT3REC**0.7_JPRB
   !VH 12/01/2018: change to JPL-15  
   RATES_LUT(KNO2NO3A,IK)=3.6E-30_JPRB*ZT3REC**4.1_JPRB
   RATES_LUT(KNO2NO3B,IK)=1.9E-12_JPRB*ZT3REC**(-0.2_JPRB)
   !JEW: changed to JPL2006 
   !RATES_LUT(KN2O5,IK)=ZFARR(2.7E-27_JPRB,11000._JPRB,ZTREC)
   !VH 12/01/2018: change to JPL-15  

   RATES_LUT(KN2O5A,IK)=1.3E-3_JPRB*(ZTEMP/300._JPRB)**(-3.5)
   RATES_LUT(KN2O5B,IK)=9.7E14_JPRB*(ZTEMP/300._JPRB)**(0.1)   

   RATES_LUT(KOHHONO,IK)=ZFARR(2.5E-12_JPRB,260._JPRB,ZTREC)      ! IUPAC 2017
   RATES_LUT(KHONOA,IK)=7.4E-31_JPRB*ZT3REC**2.4_JPRB             ! IUPAC 2017
   RATES_LUT(KHONOB,IK)=3.6E-11_JPRB*ZT3REC**(-0.3_JPRB)


   RATES_LUT(KHNO4OH,IK)=ZFARR(3.2E-13_JPRB,690._JPRB,ZTREC)
   !JEW: changed to JPL2006  
   !RATES_LUT(KNO2HO2A,IK)=2.0E-31_JPRB*ZT3REC**3.4_JPRB
   !RATES_LUT(KNO2HO2B,IK)=2.9E-12_JPRB*ZT3REC**1.1_JPRB
   !VH 12/1/2018 changed to JPL-15
   !RATES_LUT(KNO2HO2A,IK)=1.9E-31_JPRB*ZT3REC**3.4_JPRB
   !RATES_LUT(KNO2HO2B,IK)=4.0E-12_JPRB*ZT3REC**0.3_JPRB

   RATES_LUT(KNO2HO2A,IK)=1.4E-31_JPRB*ZT3REC**3.1_JPRB
   RATES_LUT(KNO2HO2B,IK)=4.0E-12_JPRB ! independent of T

   RATES_LUT(KHNO4A,IK)=ZFARR(4.1E-5_JPRB,-10650._JPRB,ZTREC)
   RATES_LUT(KHNO4B,IK)=ZFARR(6.0E15_JPRB,-11170._JPRB,ZTREC)
   !RATES_LUT(KHNO4M,IK)=ZFARR(2.1E-27_JPRB,10900._JPRB,ZTREC)

   RATES_LUT(KMENO2A,IK)=1.0E-30_JPRB*ZT3REC**4.8_JPRB
   RATES_LUT(KMENO2B,IK)=7.2E-12_JPRB*ZT3REC**2.1_JPRB

   RATES_LUT(KMENO2M,IK)=ZFARR(9.5E-29_JPRB,11234._JPRB,ZTREC)

   !JEW: changed to JPL2006  
   RATES_LUT(KODM,IK)=.2095_JPRB*ZFARR(3.3E-11_JPRB,55._JPRB,ZTREC)+ &
    &  .7808*ZFARR(2.15E-11_JPRB,110._JPRB,ZTREC)
   !JEW: changed to JPL2006
   RATES_LUT(KH2OOD,IK)=ZFARR(1.63E-10_JPRB,60._JPRB,ZTREC)
   RATES_LUT(KO3HO2,IK)=ZFARR(1.0E-14_JPRB,-490._JPRB,ZTREC)
   RATES_LUT(KO3OH,IK)=ZFARR(1.7E-12_JPRB,-940._JPRB,ZTREC)
   !JEW: changed to JPL2006  
   RATES_LUT(KHPOH,IK)=1.8E-12_JPRB
   !JEW: changed to JPL2006  
   RATES_LUT(KFRMOH,IK)=ZFARR(5.5E-12_JPRB,125._JPRB,ZTREC)
   !JEW: changed to JPL2006  
   RATES_LUT(KCH4OH,IK)=ZFARR(2.45E-12_JPRB,-1775._JPRB,ZTREC)
   RATES_LUT(KCOOHA,IK)=5.9E-33_JPRB*ZT3REC**1.4_JPRB
   RATES_LUT(KCOOHB,IK)=1.1E-12_JPRB*ZT3REC**(-1.3_JPRB)
   RATES_LUT(KCOOHC,IK)=1.5E-13_JPRB*ZT3REC**(-0.6_JPRB)
   RATES_LUT(KCOOHD,IK)=2.1E9_JPRB*ZT3REC**(-6.1_JPRB)
   !JEW: changed to JPL2006  
   RATES_LUT(KOHMPER,IK)=ZFARR(3.8E-12_JPRB,200._JPRB,ZTREC)
!   rates_lut(kohrooh,IK)=zfarr(3.01e-12,190.,ztrec) ! CB05
!
! modified according to the findings of Archibald et al., AE, 2011.
! the lifetime of ROOH is too long in CB05 compared to other mechanisms
! therefore reduce
!      
   RATES_LUT(KOHROOH,IK)=2.0E-11_JPRB
   RATES_LUT(KHO2OH,IK)=ZFARR(4.8E-11_JPRB,250._JPRB,ZTREC)
   RATES_LUT(KH2OH,IK)=ZFARR(2.8E-12_JPRB,-1800._JPRB,ZTREC)
   !JEW: changed to JPL2006  
   ZR=(1._JPRB/(1._JPRB+(498._JPRB*EXP(-1160._JPRB/ZTEMP))))
   RATES_LUT(KMO2HO2A,IK)=ZFARR(3.8E-13_JPRB,780._JPRB,ZTREC)*(1._JPRB-ZR) 
   RATES_LUT(KMO2HO2B,IK)=ZFARR(3.8E-13_JPRB,780._JPRB,ZTREC)*ZR
   RATES_LUT(KMO2NO,IK)=ZFARR(2.8E-12_JPRB,300._JPRB,ZTREC)      
   RATES_LUT(KMO2MO2,IK)=ZFARR(9.5E-14_JPRB,390._JPRB,ZTREC)     
   !JEW: changed to JPL2006  
   !RATES_LUT(KHO2HO2A,IK)=ZFARR(3.5E-13_JPRB,430._JPRB,ZTREC)
   !RATES_LUT(KHO2HO2B,IK)=ZFARR(1.7E-33_JPRB,1000._JPRB,ZTREC)
   !RATES_LUT(KHO2HO2C,IK)=ZFARR(1.4E-21_JPRB,2200._JPRB,ZTREC)
   !VH 11/1/2018: changed to JPL 15  
   RATES_LUT(KHO2HO2A,IK)=ZFARR(3.0E-13_JPRB,460._JPRB,ZTREC)
   RATES_LUT(KHO2HO2B,IK)=ZFARR(2.1E-33_JPRB,920._JPRB,ZTREC)
   RATES_LUT(KHO2HO2C,IK)=ZFARR(1.4E-21_JPRB,2200._JPRB,ZTREC)

   RATES_LUT(KC41,IK)=5.8E-16_JPRB
   RATES_LUT(KC46,IK)=ZFARR(8.1E-12_JPRB,270._JPRB,ZTREC)
   ! from IUPAC (Atkinson et al, 2006)       
   !RATES_LUT(KC47A,IK)=2.7E-28_JPRB*ZT3REC**7.1_JPRB
   !RATES_LUT(KC47B,IK)=1.2E-11_JPRB*ZT3REC**0.9_JPRB
   ! VH 12/01/2018: IUPAC (Updated PAN formation, as of 2014)       
   RATES_LUT(KC47A,IK)=9.7E-29_JPRB*ZT3REC**5.6_JPRB
   RATES_LUT(KC47B,IK)=9.3E-12_JPRB*ZT3REC**1.5_JPRB
   !RATES_LUT(KC48A,IK)=ZFARR(4.9E-3_JPRB,-12100._JPRB,ZTREC)
   !RATES_LUT(KC48B,IK)=ZFARR(5.4E16_JPRB,-13830._JPRB,ZTREC)
   !VH 12/01/2018: IUPAC (Updated PAN formation, as of 2014)   
   RATES_LUT(KC48,IK)=ZFARR(9.0E-29_JPRB,14000._JPRB,ZTREC)
   !JEW: changed to JPL2006 
   RATES_LUT(KC49,IK)=ZFARR(2.9E-12_JPRB,500._JPRB,ZTREC)
   !VH 12/01/2018 - introduce two rate equations, according to branching probability,
   ! following Gross et al (and MOM)
   RATES_LUT(KC50A,IK)=ZFARR(5.2E-13_JPRB,980._JPRB,ZTREC)*1.507*0.84
   RATES_LUT(KC50B,IK)=ZFARR(5.2E-13_JPRB,980._JPRB,ZTREC)*1.507*0.16
   !------------------------------------------------------  
   RATES_LUT(KC52,IK)=8.1E-13_JPRB
   RATES_LUT(KC53,IK)=ZFARR(1.E15_JPRB,-8000._JPRB,ZTREC)
   RATES_LUT(KC54,IK)=1.6E3_JPRB
   !JEW: changed to JPL2006
   RATES_LUT(KC57,IK)=ZFARR(5.4E-12_JPRB,-610._JPRB,ZTREC)
   !JEW: changed to JPL2006  
   RATES_LUT(KC58,IK)=ZFARR(8.5E-16_JPRB,1520._JPRB,ZTREC)
   RATES_LUT(KC59,IK)=ZFARR(4.6E-14_JPRB,400._JPRB,ZTREC)
   !JEW: changed to JPL2006  
   !RATES_LUT(KC61A,IK)=1.E-28_JPRB*ZT3REC**4.5_JPRB
   !RATES_LUT(KC61B,IK)=8.8E-12_JPRB*ZT3REC**0.85_JPRB
   !VH 12/01/2018: changed to JPL-15  
   RATES_LUT(KC61A,IK)=1.1E-28_JPRB*ZT3REC**3.5_JPRB
   RATES_LUT(KC61B,IK)=8.4E-12_JPRB*ZT3REC**1.75_JPRB
   !JEW: changed to IUPAC 2019 
   RATES_LUT(KC62,IK)=ZFARR(6.82E-15_JPRB,-2500._JPRB,ZTREC)
   !JEW: changed to IUPAC2004    
   !RATES_LUT(KC73,IK)=1.5E-11_JPRB ! IUPAC
   !VH 18/01/2018 CH3COCHO+OH updated following IUPAC
   RATES_LUT(KC73,IK)=ZFARR(1.9E-12_JPRB,575._JPRB,ZTREC)
   RATES_LUT(KC76,IK)=ZFARR(2.7E-11_JPRB,390._JPRB,ZTREC) ! IUPAC 2009
   RATES_LUT(KC77,IK)=ZFARR(1.04E-14_JPRB,-1995._JPRB,ZTREC) ! IUPAC
   RATES_LUT(KC78,IK)=ZFARR(2.95E-12_JPRB,-450._JPRB,ZTREC) ! IUPAC
   !JEW: changed to JPL2006  
   RATES_LUT(KC79,IK)=ZFARR(2.6E-12_JPRB,365._JPRB,ZTREC)
   RATES_LUT(KC80,IK)=ZFARR(1.6E-12_JPRB,-2200._JPRB,ZTREC)
   RATES_LUT(KC81,IK)=ZFARR(2.6E-12_JPRB,365._JPRB,ZTREC) !CB05   
   RATES_LUT(KC82,IK)=ZFARR(7.5E-13_JPRB,700._JPRB,ZTREC) ! CB05
   RATES_LUT(KC83,IK)=8.E-11_JPRB
   RATES_LUT(KC84,IK)=ZFARR(5.9E-13_JPRB,-360._JPRB,ZTREC) ! CB05 temp dep
   RATES_LUT(KC85,IK)=ZFARR(8.0E-12_JPRB,-2060._JPRB,ZTREC) ! CB05
   RATES_LUT(KO3PO2,IK)=6.0E-34_JPRB*ZT3REC**2.4_JPRB
   RATES_LUT(KO3PO3,IK)=ZFARR(8.0E-12_JPRB,-2060._JPRB,ZTREC)       

   ! sulfur and ammonia gas phase reactions
   ! the hynes et al. parameterisation is replaced by chin et al. 1996

   !JEW: changed to JPL2006
   RATES_LUT(KDMSOHA,IK)= 1.11E-11_JPRB*EXP(-240._JPRB/ZTEMP)
   RATES_LUT(KDMSOHB,IK)=1.0E-9_JPRB*EXP(5820._JPRB/ZTEMP)
   RATES_LUT(KDMSOHC,IK)=5.0E-0_JPRB*EXP(6280._JPRB/ZTEMP)
   RATES_LUT(KDMSNO3,IK)=ZFARR(1.9E-13_JPRB,520._JPRB,ZTREC)!at 298 1.01e-12
   !JEW: changed to JPL2006  
   RATES_LUT(KSO2OHA,IK)=3.3E-31_JPRB*(ZTEMP/300._JPRB)**(-4.3_JPRB)
   RATES_LUT(KSO2OHB,IK)= 1.6E-12_JPRB*(ZTEMP/300._JPRB)
   !VH rates for OCS and SO3
   RATES_LUT(KSO2O3G,IK)=ZFARR(3.0E-12_JPRB,-7000._JPRB,ZTREC)

   ! First compute in DP..
   ZRX1=8.5E-41_JPRD*EXP(6540_JPRD*ZTREC)
   RATES_LUT(KSO3H2O,IK)=ZRX1
   !RATES_LUT(KSO3H2O,IK)=ZFARR(8.5E-41_JPRB,6540._JPRB,ZTREC)
   RATES_LUT(KOCSOH,IK) =ZFARR(1.1E-13_JPRB,-1200._JPRB,ZTREC)
   !
   ! fate of ammonia
   !
   RATES_LUT(KNH3OH,IK)=  ZFARR(1.7E-12_JPRB,-710._JPRB,ZTREC)!1.56e-13 at 298k
   RATES_LUT(KNH2OH,IK)=  0.0_JPRB !Reaction /  Reference not reported elsewhere.
   RATES_LUT(KNH2NO,IK)=  ZFARR(4.0E-12_JPRB,+450._JPRB,ZTREC)!1.72e-11
   RATES_LUT(KNH2NO2,IK)= ZFARR(2.1E-12_JPRB,650._JPRB,ZTREC)!1.86e-11
   RATES_LUT(KNH2HO2,IK)= 3.4E-11_JPRB
   RATES_LUT(KNH2O2,IK)= 6.0E-21_JPRB
   !VH RATES_LUT(KNH2O3,IK)= ZFARR(4.3E-12_JPRB,-930._JPRB,ZTREC)!1.89e-13 at 298k
   !
   ! for higher organics
   RATES_LUT(KOHMCHO,IK) = ZFARR(4.4E-12_JPRB,365._JPRB,ZTREC) ! IUPAC
   RATES_LUT(KOHMCH2CHO,IK) = ZFARR(4.9E-12_JPRB,405._JPRB,ZTREC)
   
   RATES_LUT(KNO3MCHO,IK) = ZFARR(1.4E-12_JPRB,-1860._JPRB,ZTREC) 
   RATES_LUT(KNO3MCH2CHO,IK) = 6.4E-15_JPRB
   
   RATES_LUT(KOHOLE,IK) = ZFARR(8.2E-12_JPRB,610._JPRB,ZTREC)   ! IUPAC
   
   RATES_LUT(KO3OLE,IK) = 1.0E-17_JPRB
   
   RATES_LUT(KNO3OLE,IK) = ZFARR(4.0E-14_JPRB,-400._JPRB,ZTREC)
   !
   ! the rates for additional BVOC's 
   !
   RATES_LUT(KOHCH3OH,IK) = ZFARR(2.85E-12_JPRB,-345._JPRB,ZTREC)
   RATES_LUT(KOHHCOOH,IK) = 4.0E-13_JPRB
   !RATES_LUT(KOHMCOOH,IK) = ZFARR(4.2E-14_JPRB,-855._JPRB,ZTREC)
   !VH 12/01/2018 Update following latest IUPAC + bugfix
   RATES_LUT(KOHMCOOH,IK) = ZFARR(4.0E-14_JPRB,850._JPRB,ZTREC)
   RATES_LUT(KOHC2H6,IK) = ZFARR(6.9E-12_JPRB,-1000._JPRB,ZTREC)
   RATES_LUT(KOHETHOH,IK) = ZFARR(3.0E-12_JPRB,20._JPRB,ZTREC)
   RATES_LUT(KOHC3H8,IK) = ZFARR(7.6E-12_JPRB,-585._JPRB,ZTREC)
   RATES_LUT(KOHC3H6A,IK) = 8.0E-27_JPRB*ZT3REC**3.5_JPRB
   RATES_LUT(KOHC3H6B,IK) = 3.0E-11_JPRB*ZT3REC
   RATES_LUT(KO3C3H6,IK) = ZFARR(5.5E-15_JPRB,-1880._JPRB,ZTREC) ! IUPAC
   RATES_LUT(KNO3C3H6,IK) = ZFARR(4.6E-13_JPRB,-1155._JPRB,ZTREC) ! IUPAC 
   !
   ! Oxidation products from C3H6 and C3H8
   !
   ! Following MOZART (Emmons et al., GMD 2010)
   !
   RATES_LUT(KNOHYPROPO2,IK)=4.2E-12_JPRB*EXP(180._JPRB/ZTEMP)
   RATES_LUT(KHO2HYPROPO2,IK)= 7.5E-13_JPRB*EXP(700._JPRB/ZTEMP)

   RATES_LUT(KNOIC3H7O2,IK)=4.2E-12_JPRB*EXP(180._JPRB/ZTEMP)
   RATES_LUT(KHO2IC3H7O2,IK)= 7.5E-13_JPRB*EXP(700._JPRB/ZTEMP)
   !
   ! Terpenes
   ! 
   RATES_LUT(KOHTERP,IK) =  ZFARR(1.2E-11_JPRB,440._JPRB,ZTREC) ! IUPAC
   !RATES_LUT(KO3TERP,IK) =  ZFARR(6.3E-16_JPRB,-580._JPRB,ZTREC) ! IUPAC
   !VH update 18/01/2018 following latest IUPAC for alpha-pinene
   RATES_LUT(KO3TERP,IK) =  ZFARR(8.05E-16_JPRB,-640._JPRB,ZTREC) ! IUPAC
   RATES_LUT(KNO3TERP,IK) = ZFARR(1.2E-12_JPRB,490._JPRB,ZTREC) ! IUPAC      
   !
   ! Acetone
   ! 
   RATES_LUT(KOHACETA,IK) = ZFARR(8.8E-12_JPRB,-1320._JPRB,ZTREC) ! IUPAC
   RATES_LUT(KOHACETB,IK) = ZFARR(1.7E-14_JPRB,423._JPRB,ZTREC) ! IUPAC
   
   RATES_LUT(KACO2HO2,IK) = 1.0E-11_JPRB
   RATES_LUT(KACO2MO2,IK) = 3.8E-12_JPRB ! IUPAC
   RATES_LUT(KACO2NO,IK) =  8.0E-12_JPRB
   
   RATES_LUT(KACO2XO2,IK) = 6.8E-14_JPRB
   RATES_LUT(KXO2XO2N,IK) = 6.8E-14_JPRB
   RATES_LUT(KXO2N,IK) = 6.8E-14_JPRB
   !
   ! Biomass burning tracers (JEW : 2020)
   !
   RATES_LUT(KOHHCN,IK)=ZFARR(1.2E-13_JPRB,-400._JPRB,ZTREC)   
   RATES_LUT(KOHCH3CN,IK)=ZFARR(8.1E-13_JPRB,-1080._JPRB,ZTREC)
   !
   ! Taken from Karl et al., ACP, 2009 (JEW : 2020)
   !
   RATES_LUT(KOHXYL,IK) = 1.7E-11_JPRB
   RATES_LUT(KCXYLO3A,IK) = ZFARR(2.4E-13_JPRB,-5586._JPRB,ZTREC)
   RATES_LUT(KCXYLO3B,IK) = ZFARR(5.37E-13_JPRB,-6039._JPRB,ZTREC)
   RATES_LUT(KNO3XYL,IK) = 3.54E-16_JPRB    
   RATES_LUT(KOHTOL,IK) = 5.96E-12_JPRB
   RATES_LUT(KO3TOL,IK) = ZFARR(2.34E-12_JPRB,-6694._JPRB,ZTREC)
   RATES_LUT(KNO3TOL,IK) = 7.8E-17_JPRB
   RATES_LUT(KAROO2NO,IK) = ZFARR(4.2E-12_JPRB,180._JPRB,ZTREC)
   RATES_LUT(K2AROO2,IK) = ZFARR(1.7E-14_JPRB,1300._JPRB,ZTREC)
   RATES_LUT(KAROO2HO2,IK) = ZFARR(3.5E-13_JPRB,1000._JPRB,ZTREC)
   RATES_LUT(KAROO2XO2,IK) = ZFARR(1.7E-14_JPRB,1300._JPRB,ZTREC)      
   !
   !  MACR, MVK
   !                     
   RATES_LUT(KOHMVK,IK)=ZFARR(2.6E-12_JPRB,610._JPRB,ZTREC)
   RATES_LUT(KOHMACR,IK)=ZFARR(8.0E-12_JPRB,380._JPRB,ZTREC)     
   RATES_LUT(KO3MVK,IK)=ZFARR(8.5E-16_JPRB,-1520._JPRB,ZTREC)
   RATES_LUT(KO3MACR,IK)=ZFARR(1.4E-15_JPRB,-2100._JPRB,ZTREC)
   !
   !   Stavrakou (2010) isoprene mechanism
   !
   ! RATES_LUT(KISOPBO2A,IK)=ZFARR(4.06E+9_JPRB,-7302._JPRB,ZTREC)
   RATES_LUT(KISOPBO2B,IK)=ZFARR(2.08E+11_JPRB,-8993._JPRB,ZTREC)     
   ! RATES_LUT(KISOPDO2A,IK)=ZFARR(8.5E+9_JPRB,-7342._JPRB,ZTREC)
   RATES_LUT(KISOPDO2B,IK)=ZFARR(2.08E+11_JPRB,-8993._JPRB,ZTREC)   

   ! VH Follow JAMC reaction rate expressions (~ slower ISOPxO2->HO2; faster ISOPxO2+HO2->ISOPOOH)
   RATES_LUT(KISOPBO2A,IK)=ZFARR(4.1E+8_JPRB,-7700._JPRB,ZTREC)
   RATES_LUT(KISOPDO2A,IK)=ZFARR(4.1E+8_JPRB,-7700._JPRB,ZTREC)
   RATES_LUT(KISOPBO2HO2,IK)=ZFARR(2.05E-13_JPRB,1300._JPRB,ZTREC)
   RATES_LUT(KISOPDO2HO2,IK)=ZFARR(2.05E-13_JPRB,1300._JPRB,ZTREC)

   ! Stavrakou et al., 2010
   ! RATES_LUT(KISOPBO2HO2,IK)=ZFARR(8.0E-13_JPRB,700._JPRB,ZTREC)
   RATES_LUT(KISOPBO2NO,IK)=ZFARR(4.4E-12_JPRB,180._JPRB,ZTREC)     
   ! RATES_LUT(KISOPDO2HO2,IK)=ZFARR(8.0E-13_JPRB,700._JPRB,ZTREC)
   RATES_LUT(KISOPDO2NO,IK)=ZFARR(4.4E-12_JPRB,180._JPRB,ZTREC)
   
   RATES_LUT(KISOPOOHOH,IK)=ZFARR(1.52E-11_JPRB,200._JPRB,ZTREC)
   RATES_LUT(KHPALD1OH,IK)=ZFARR(1.86E-11_JPRB,175._JPRB,ZTREC)
   RATES_LUT(KHPALD2OH,IK)=ZFARR(1.86E-11_JPRB,175._JPRB,ZTREC)   
   RATES_LUT(KGLYOH,IK)=ZFARR(3.1E-12_JPRB,340._JPRB,ZTREC)
   RATES_LUT(KGLYALDOH,IK)=8.0E-12
   RATES_LUT(KHYACOH,IK)=ZFARR(2.0E-12_JPRB,320._JPRB,ZTREC)    
   !
   ! **** solubility Henry equilibrium
   !      HNO3/so4/nh4 just a very high number to take H and 
   !      dissociation into account
   !
   HENRY(:,IK)=0._JPRB
   HENRY(IH2O2,IK)=1.0E5_JPRB*EXP(6300_JPRB*ZTREC)*6.656E-10_JPRB
   HENRY(IHNO3,IK)=1E7_JPRB 
   HENRY(ICH3O2H,IK)=310._JPRB*EXP(5200_JPRB*ZTREC)*2.664E-8_JPRB
   ! HENRY(ich2o,IK)=3000.*exp(7200*ztrec)*exp(-7200.*(1./298.15))
   HENRY(ICH2O,IK)=3000._JPRB*EXP(7200_JPRB*ZTREC)*3.253E-11_JPRB
   ! account for enhanced solubility (i.e. KHeff)
   HENRY(ICH2O,IK)=HENRY(ICH2O,IK)*37_JPRB
   HENRY(IROOH,IK)=340._JPRB*EXP(6000._JPRB*ZTREC)*1.821E-9_JPRB
   HENRY(IORGNTR,IK)=ZFARR(1.8E-6_JPRB,6000._JPRB,ZTREC) 
   HENRY(ISO4,IK)=1.E7_JPRB
   HENRY(INH4,IK)=1.E7_JPRB
   HENRY(IMSA,IK)=1.E7_JPRB
   !!JW/VH       HENRY(iso2,IK) =1.2*exp(3120.*ZTREC)*3.41e-5 !correction for the 1/298. part
   !!JW : " The second term (exp (-3120* (1/298.15)) should equal 2.85e-5 "
   HENRY(ISO2,IK) =1.2_JPRB*EXP(3120._JPRB*ZTREC)*2.85E-5_JPRB 
   HENRY(INH3,IK) =75.0_JPRB*EXP(3400._JPRB*ZTREC)*1.10E-5_JPRB
   HENRY(IO3,IK)=1.1E-2_JPRB*EXP(2300._JPRB*ZTREC)*4.45E-4_JPRB

  ! JEW add two new scavenging rates for CH3COCHO and ALD2
  ! need KH(eff) due to hydration steps for both species
   HENRY(IMGLY,IK) = 3.2E4_JPRB*48.6_JPRB
   ZKH1=17._JPRB*EXP(5000._JPRB*ZTREC)*EXP(-5000._JPRB*(1._JPRB/298.15_JPRB))
   ZKH2=13._JPRB*EXP(5700._JPRB*ZTREC)*EXP(-5700._JPRB*(1._JPRB/298.15_JPRB))
   HENRY(IALD2,IK)  = ((ZKH1+ZKH2)/2._JPRB)*1.0246_JPRB
   HENRY(IHCOOH,IK) = 8900._JPRB*EXP(6100._JPRB*ZTREC)*EXP(-6100._JPRB*(1._JPRB/298.15_JPRB))
   HENRY(ICH3OH,IK) = 220.0_JPRB*EXP(5200._JPRB*ZTREC)*2.66E-8_JPRB
   HENRY(IMCOOH,IK) = 4100._JPRB*EXP(6300._JPRB*ZTREC)*EXP(-6300._JPRB*(1._JPRB/298.15_JPRB))
   HENRY(IETHOH,IK) = 190.0_JPRB*EXP(6600._JPRB*ZTREC)*EXP(-6600._JPRB*(1._JPRB/298.15_JPRB))
   HENRY(IACET,IK)  = 35._JPRB*EXP(3800._JPRB*ZTREC)*EXP(-3800._JPRB*(1._JPRB/298.15_JPRB))
   HENRY(IHCN,IK)  = 12._JPRB*EXP(5000._JPRB*ZTREC)*EXP(-5000._JPRB*(1._JPRB/298.15_JPRB))
  ! Snider and Davis (1985) values chosen as for ETHOH   
   HENRY(ICH3CN,IK)  = 49._JPRB*EXP(4000._JPRB*ZTREC)*EXP(-4000._JPRB*(1._JPRB/298.15_JPRB))
ENDDO !IK temperature loop
!
! marked tracers:
!
HENRY(IO3S,:) = HENRY(IO3,:)

IF (LHOOK) CALL DR_HOOK('BASCOETM5_CHEM_INI:BASCOETM5_RATES',1,ZHOOK_HANDLE)

END SUBROUTINE BASCOETM5_RATES


FUNCTION ZFARR(PRX1,PER,PTREC)
  !------------------------------------------------------------------
  !
  !****  ZFARR calculation of Arrhenius expression for rate constants
  !
  !------------------------------------------------------------------
  !
  IMPLICIT NONE
  REAL(KIND=JPRB)            :: ZFARR
  REAL(KIND=JPRB),INTENT(IN) :: PRX1,PER
  REAL(KIND=JPRB),INTENT(IN) :: PTREC
  !
  ZFARR=PRX1*EXP(PER*PTREC)
  !
END FUNCTION ZFARR


END SUBROUTINE BASCOETM5_CHEM_INI 


