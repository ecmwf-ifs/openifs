! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_CHEM_INI(YGFL,YDCHEM)

!**   DESCRIPTION 
!     ----------
!
!   BASCOE routine for IFS chemistry : Initialization of chem-rates lookup-table
!       and photolysis lookup-table
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_CHEM_INI* IS CALLED FROM *CHEM_INIT*.
!

! INPUTS: 
! -------
! YGFL
!
! OUTPUTS: none
! -------
!
!
!     AUTHOR
!     -------
!     2014-02-01: VINCENT HUIJNEN    *KNMI*
!
!     MODIFICATIONS.
!     --------------
!     2017-12-07: YVES CHRISTOPHE (YC) *BIRA*
!                   photolysis lookup tables initialized from a TUV based radiative
!                   transfer code instead of read from exernal file
!     2018-03-27: YVES CHRISTOPHE (YC) *BIRA*
!                   photolysis initialization is done from BASCOE_J_INI, which only
!                   prepares environment to run TUV (compute J rates) on line 
!     2021-05-20 JONAS DEBOSSCHER (JD) *BIRA*
!                   Aerosol SAD climatology initialization


USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCHEM  , ONLY : TCHEM
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE BASCOE_MODULE, ONLY : NBC, BASCOE_BCVAL, BASCOE_BCNAME



IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------


! * LOCAL 
TYPE(TYPE_GFLD),INTENT(INOUT)   :: YGFL
TYPE(TCHEM),    INTENT(IN)      :: YDCHEM
REAL(KIND=JPHOOK)                 :: ZHOOK_HANDLE

#include "bascoe_j_ini.intfb.h"
#include "bascoe_lbc_ini.intfb.h"
#include "bascoe_setbin.intfb.h"
#include "bascoe_sage_init.intfb.h"
#include "bascoe_tropopause_init.intfb.h"
#include "bascoe_climSAD_ini.intfb.h"


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('BASCOE_CHEM_INI',0,ZHOOK_HANDLE)

 ! Some checking
 CALL TRACER_IDX_CHECK_BASCOE(YGFL)
 
 IF (YDCHEM%LCHEM_BASCOE_JON) THEN
   ! Online J-rate computation
   ! Prepare data to compute photolysis rates
   CALL BASCOE_J_INI
 ELSE
   ! Initialize J-lookup tables
   CALL BASCOE_J_TABLES_INIT
 ENDIF

  ! Prepare data with lower boundary conditions
 CALL BASCOE_LBC_INI(NBC, BASCOE_BCVAL, BASCOE_BCNAME)

  ! Prepare strat. aerosol particle size distribution 
 CALL BASCOE_SETBIN
 
 ! SAGE initialization
 CALL BASCOE_SAGE_INIT

 ! Tropopause pressure level initialization
 CALL BASCOE_TROPOPAUSE_INIT
 
 ! JD: Aerosol SAD climatology initialization
 CALL BASCOE_CLIMSAD_INI

IF (LHOOK) CALL DR_HOOK('BASCOE_CHEM_INI',1,ZHOOK_HANDLE)

CONTAINS


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   auxiliary routines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!       Check correct tracer indices
!       -----------------------------------------------------------
SUBROUTINE TRACER_IDX_CHECK_BASCOE(YGFL)


USE BASCOE_MODULE  , ONLY : &
 & IO3,  IH2O2,  ICH4, ICO,  IHNO3,  ICH3OOH,  ICH2O, INO, &
 & IHO2 ,ICH3 ,ICH3O ,IHCO ,ICH3O2 ,IOH ,INO2 ,IN2O5 , &
 & IHO2NO2 ,INO3 ,IN2O ,IH2O ,IOCLO ,IHCL ,ICLONO2 ,IHOCL , &
 & ICL2 ,IHBR ,IBRONO2,ICL2O2,IHOBR,IBRCL ,ICFC11 ,ICFC12 ,&
 & ICFC113 ,ICFC114 ,ICFC115 ,&
 & ICCL4 ,ICLNO2 ,ICH3CCL3 ,ICH3CL ,IHCFC22 ,ICH3BR ,IHF ,IHA1301 ,&
 & IHA1211 ,ICHBR3 ,ICLOO ,IO ,IO1D ,IN ,ICLO ,ICL ,&
 & IBR ,IBRO ,IH ,IH2 ,ICO2 ,IBR2, ICH2BR2, ISTRATAER
 
   
USE YOMLUN             , ONLY : NULOUT
USE PARKIND1           , ONLY : JPIM  , JPRB
USE YOM_YGFL           , ONLY : TYPE_GFLD
USE YOMHOOK            , ONLY  : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE

! Local parameters

! * counters
TYPE(TYPE_GFLD),INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM) :: JL
LOGICAL            :: LLFOUND
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BASCOE_CHEM_INI:TRACER_IDX_CHECK_BASCOE',0,ZHOOK_HANDLE)

ASSOCIATE(NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM)
   DO JL = 1,NCHEM  
      LLFOUND = .FALSE.
     SELECT CASE ( TRIM( YCHEM(JL)%CNAME ) ) 
       CASE ('O3')    ; LLFOUND = ( IO3 == JL )
       CASE ('H2O2')  ; LLFOUND = (IH2O2 == JL ) 
       CASE ('CH4')   ; LLFOUND = (ICH4 == JL )
       CASE ('CO')    ; LLFOUND = (ICO == JL )
       CASE ('HNO3')  ; LLFOUND = (IHNO3 == JL )
       CASE ('CH3OOH'); LLFOUND = (ICH3OOH == JL )
       CASE ('CH2O')  ; LLFOUND = (ICH2O == JL )
       CASE ('NO')    ; LLFOUND = (INO == JL )
       CASE ('HO2')   ; LLFOUND = (IHO2 == JL )
       CASE ('CH3')   ; LLFOUND = (ICH3 == JL )
       CASE ('CH3O')  ; LLFOUND = (ICH3O == JL )
       CASE ('HCO')   ; LLFOUND = (IHCO == JL )
       CASE ('CH3O2') ; LLFOUND = (ICH3O2 == JL )
       CASE ('OH')    ; LLFOUND = (IOH == JL )
       CASE ('NO2')   ; LLFOUND = (INO2 == JL )
       CASE ('N2O5')  ; LLFOUND = (IN2O5 == JL )
       CASE ('HO2NO2'); LLFOUND = (IHO2NO2 == JL )
       CASE ('NO3')   ; LLFOUND = (INO3 == JL )
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
       CASE ('NOXA')  ; LLFOUND = .TRUE. 
       CASE ('CLXA')  ; LLFOUND = .TRUE. 
       CASE ('BRXA')  ; LLFOUND = .TRUE.
       CASE ('STRATAER')   ; LLFOUND = (ISTRATAER == JL )
       CASE DEFAULT
         WRITE(NULOUT,*) 'ERROR BASCOE_chem_ini: no matching tracer name for '//TRIM(YCHEM(JL)%CNAME)
         CALL ABOR1('BASCOE_chem_ini: No matching tracer name available')
     END SELECT
 
     IF (.NOT. LLFOUND ) THEN 
       WRITE(NULOUT,*) 'ERROR BASCOE_chem_ini: Wrong tracer index or status for  '//TRIM(YCHEM(JL)%CNAME)
       CALL ABOR1('BASCOE_chem_ini: wrong tracer index tracer name')
     ENDIF   
  ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('BASCOE_CHEM_INI:TRACER_IDX_CHECK_BASCOE',1,ZHOOK_HANDLE)

 END SUBROUTINE TRACER_IDX_CHECK_BASCOE


END SUBROUTINE BASCOE_CHEM_INI 
