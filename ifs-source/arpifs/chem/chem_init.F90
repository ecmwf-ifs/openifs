! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_INIT(YDGEOMETRY,YDML_GCONF,YDDYNA,YDML_CHEM)

!**   DESCRIPTION
!     ----------
!
!   init routine for IFS chemistry
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_INIT* IS CALLED FROM *CNT4*.
!

!     Externals.
!     ---------
!                 CHEM_INIT_MOCAGE
!                 CHEM_INIT_MOZART
!                 CHEM_INIT_TM5

!
!     AUTHOR.
!     -------
!        JOHANNES FLEMMING  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2009-09-11
!        18-09, M. Michou, call to ARPEGE-Climat 6.3 chemistry implemented (scheme arpclim_repro)

USE YOMLUN,       ONLY : NULOUT
USE MODEL_CHEM_MOD , ONLY : MODEL_CHEM_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYNA  , ONLY : TDYNA
USE YOMCST , ONLY : RMV
USE TM5_CHEM_MODULE, ONLY : IFLRN222, IFLPB210, IFLO3, IFLCO

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.5   LOCAL VARIABLES
!              ---------------

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
TYPE(MODEL_CHEM_TYPE),INTENT(INOUT):: YDML_CHEM
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TDYNA),                  INTENT(IN)   :: YDDYNA
REAL(KIND=JPRB)    :: ZMW
CHARACTER(LEN=10)  :: CLNAME
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
INTEGER(KIND=JPIM)    :: ICC1, ICC2, ICC3, ICHEM

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!#include "chem_init_mocage.intfb.h"
#include "tm5_chem_ini.intfb.h"
#include "bascoe_chem_ini.intfb.h"
#include "bascoetm5_chem_ini.intfb.h"
#include "n2o_chem_ini.intfb.h"
#include "suvolc.intfb.h"
#include "linco_chem_ini.intfb.h"




#include "arpclim_chem_ini.intfb.h"
#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('CHEM_INIT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDCOMPO=>YDML_CHEM%YRCOMPO,YDCHEM=>YDML_CHEM%YRCHEM,YGFL=>YDML_GCONF%YGFL,&
 & YDRYDEP=>YDML_CHEM%YRDRYDEP)
ASSOCIATE(NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, NGHG=>YGFL%NGHG, NACTAERO =>YGFL%NACTAERO, &
 & CHEM_SCHEME=>YDCHEM%CHEM_SCHEME, IEXTR_FREE=>YDCHEM%IEXTR_FREE, &
 & LCHEM_DIA=>YDCOMPO%LCHEM_DIA, KCHEM_DRYDEP=>YDCHEM%KCHEM_DRYDEP)
!*             Init chem scheme
!              ---------------

IF(CHEM_SCHEME(1:4) == 'MIN_') THEN

  WRITE(NULOUT,*) 'CHEM_INIT - NO INITIALIZATION FOR MINIMIZATION FOR NOW'

ELSE

   SELECT CASE (TRIM(CHEM_SCHEME))
      CASE ("decay")

        CALL SUVOLC(YDGEOMETRY)
      CASE ("bascoe")

        CALL BASCOE_CHEM_INI(YGFL,YDCHEM)

      CASE ("n2o")

        CALL N2O_CHEM_INI

      CASE ("tm5")

        CALL TM5_CHEM_INI(YGFL,YDCHEM,YDCOMPO)

      CASE ("bascoetm5")

        CALL BASCOETM5_CHEM_INI(YGFL,YDCHEM,YDCOMPO)

      CASE ("linco")

        CALL LINCO_CHEM_INI(YDGEOMETRY,YDML_GCONF,YDDYNA,YDCHEM)

      CASE ("carbontracers")

        CALL LINCO_CHEM_INI(YDGEOMETRY,YDML_GCONF,YDDYNA,YDCHEM)

      CASE ("nwpo3")

        IFLO3=-999
        DO ICHEM=1,NCHEM
          IF ( TRIM(YCHEM(ICHEM)%CNAME) == 'O3') IFLO3=ICHEM
        ENDDO
        IF (IFLO3 == -999) CALL ABOR1(" Missing  O3 chemistry field for  "//TRIM(CHEM_SCHEME) )

      CASE ("RnPb")

        IFLRN222=-999
        IFLPB210=-999
        DO ICHEM=1,NCHEM
          IF ( TRIM(YCHEM(ICHEM)%CNAME) == 'Rn') IFLRN222=ICHEM
          IF ( TRIM(YCHEM(ICHEM)%CNAME) == 'Pb') IFLPB210=ICHEM
        ENDDO
        ! can be run without Pb
        IF (IFLRN222 == -999) CALL ABOR1(" Missing  Rn chemistry field for  "//TRIM(CHEM_SCHEME) )

      CASE ("linco_RnPb")

        IFLRN222=-999
        IFLPB210=-999
        IFLCO=-999
        DO ICHEM=1,NCHEM
          IF ( TRIM(YCHEM(ICHEM)%CNAME) == 'Rn') IFLRN222=ICHEM
          IF ( TRIM(YCHEM(ICHEM)%CNAME) == 'Pb') IFLPB210=ICHEM
          IF ( TRIM(YCHEM(ICHEM)%CNAME) == 'CO') IFLCO=ICHEM
        ENDDO
        ! can be run without Pb
        IF (IFLRN222 == -999) CALL ABOR1(" Missing  Rn chemistry field for  "//TRIM(CHEM_SCHEME) )
        IF (IFLCO == -999) CALL ABOR1(" Missing  CO chemistry field for  "//TRIM(CHEM_SCHEME) )
        CALL LINCO_CHEM_INI(YDGEOMETRY,YDML_GCONF,YDDYNA,YDCHEM)

      CASE ("arpclim_repro")

        CALL ARPCLIM_CHEM_INI(YDGEOMETRY,YGFL,YDCHEM)

      CASE DEFAULT

       CALL ABOR1(" NO KNOWN CHEMISTRY SCHEME "//TRIM(CHEM_SCHEME) )

    END SELECT

ENDIF

! Find un-used levels in extra arrays
IF (LCHEM_DIA) THEN
  ICC1=0
  ICC2=0
  ICC3=0
  IEXTR_FREE(:,:)  = -999_JPIM
  DO ICHEM=1,NCHEM
  ! Find fields without surface emissions or pseudo-emissions
  ! This is not a simple thing to do with flexible emission configuration.
    !IF (.NOT. ( YCHEM(ICHEM)%IGRIBSFC > 0 .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CH4'&
    !& .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'N2O'    .OR.  TRIM(YCHEM(ICHEM)%CNAME) == 'CFC11'&
    !& .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CFC12'  .OR.  TRIM(YCHEM(ICHEM)%CNAME) == 'CFC113'&
    !& .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CFC114' .OR.  TRIM(YCHEM(ICHEM)%CNAME) == 'CCL4'&
    !& .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CH3CCL3'.OR.  TRIM(YCHEM(ICHEM)%CNAME) == 'HCFC22'&
    !& .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'HA1301' .OR.  TRIM(YCHEM(ICHEM)%CNAME) == 'HA1211'&
    !& .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CH3BR'  .OR.  TRIM(YCHEM(ICHEM)%CNAME) == 'CHBR3'&
    !& .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CH3CL'&
    !& )) THEN
    !   ICC1=ICC1+1
    !   IEXTR_FREE(1,ICC1)  = ICHEM
    !   IF (ICC1 <= UBOUND(IEXTR_FREE,2)) IEXTR_FREE(1,ICC1)  = ICHEM
    ! ENDIF
  ! Find fields without dry dep or nudging at surface
    IF( .NOT. YCHEM(ICHEM)%IGRIBDV > 0 .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CH4') THEN
       ICC2=ICC2+1
       IF (ICC2 <= UBOUND(IEXTR_FREE,2)) IEXTR_FREE(2,ICC2)  = ICHEM
     ENDIF
    ! Find fields without wetdep
    IF( .NOT. YCHEM(ICHEM)%HENRYA > 0  ) THEN
      ICC3=ICC3+1
      IF (ICC3 <+ UBOUND(IEXTR_FREE,2)) IEXTR_FREE(3,ICC3)  = ICHEM
    ENDIF
  ENDDO
  ! These ugly-hack fields really shouldn't be relied upon except for testing, but
  ! unfortunately 6 of them are even in standard/operational configuration. Make sure
  ! they actually get defined even if there aren't sufficient gaps.
  ICHEM=NCHEM+NGHG+NACTAERO
  DO WHILE(ICC1 < 6)
    ICC1=ICC1+1
    ICHEM=ICHEM+1
    IEXTR_FREE(1,ICC1) = ICHEM
  ENDDO
ENDIF


! set speciec specific constants for on-line dry depostion
IF (KCHEM_DRYDEP >0) THEN

  ALLOCATE(YDRYDEP%RCHEN(NCHEM))
  ALLOCATE(YDRYDEP%RCHENXP(NCHEM))
  ALLOCATE(YDRYDEP%RCF0(NCHEM))
  ALLOCATE(YDRYDEP%RDIMO(NCHEM))

  DO  ICHEM = 1, NCHEM
    YDRYDEP%RCHEN(ICHEM)=-99999.9_JPRB
    YDRYDEP%RCHENXP(ICHEM)=-99999.9_JPRB
    YDRYDEP%RCF0(ICHEM)=-99999.9_JPRB
    YDRYDEP%RDIMO(ICHEM)=-99999.9_JPRB
    IF (YCHEM(ICHEM)%IGRIBDV > 0 ) THEN
      ZMW=YCHEM(ICHEM)%RMOLMASS
      CLNAME=YCHEM(ICHEM)%CNAME
      IF ( CLNAME  == 'O3S')  CLNAME='O3'
      YDRYDEP%RDIMO(ICHEM) = 1.0_JPRB/((RMV/ZMW)**0.5_JPRB) ! RATIO BETWEEN D_H2O / D_I MOLECULAR DIFFUSIVITIES
      YDRYDEP%RCHEN(ICHEM)=MAX(YCHEM(ICHEM)%HENRYA,1E-5_JPRB)
      YDRYDEP%RCHENXP(ICHEM)=YCHEM(ICHEM)%HENRYB
! SO2 in table is already "effcient" value of 1e5/3000
! not given in table - should we also calculate wet dep ???
      IF ( CLNAME == 'O3' .OR. CLNAME=='OX' ) THEN
        YDRYDEP%RCHEN(ICHEM) =   0.01_JPRB
        YDRYDEP%RCHENXP(ICHEM) = 2300.0_JPRB
      ENDIF
      IF ( CLNAME == 'ISPD' ) THEN
        YDRYDEP%RCHEN(ICHEM) =   6.5_JPRB
        YDRYDEP%RCHENXP(ICHEM) = 5300.0_JPRB
      ENDIF
!  normalised reactivity ! see table 19.3
      SELECT CASE (TRIM(CLNAME))
        CASE ( 'O3','H202','O3S','OX')
          YDRYDEP%RCF0(ICHEM)=1.0_JPRB
        CASE ('PAN','ROOH','ORGNTR','ONIT','N2O5','NO2','NO3','TPAN','MPAN','ONITR','HONO','OP2','PAA','HO2NO2' )
          YDRYDEP%RCF0(ICHEM)=0.1_JPRB
        CASE ('ALKOOH','BENOOH','C2H5OOH','C3H7OOH','ISOPOOH','TERPOOH','TOLOOH','XOOH','XYLOOH' )
          YDRYDEP%RCF0(ICHEM)=0.1_JPRB
        CASE('CH3OOH')
          YDRYDEP%RCF0(ICHEM)=0.3_JPRB
        CASE DEFAULT
          YDRYDEP%RCF0(ICHEM)=0.0_JPRB
       END SELECT
    ENDIF
  ENDDO
ENDIF


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_INIT
