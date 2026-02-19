! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LINCO_CHEM_INI(YDGEOMETRY,YDML_GCONF,YDDYNA,YDCHEM)

!**   DESCRIPTION 
!     ----------
!
!   LINCO routine for IFS chemistry : Initialization of the chemistry 
!
!
!
!**   INTERFACE.
!     ----------
!          *LINCO_CHEM_INI* IS CALLED FROM *CHEM_INI*.

! INPUTS:
! -------
!
! OUTPUTS:
! -------
!
!
! LOCAL:
! -------
!
!     AUTHOR.
!     -------
!        SEBASTIEN MASSART  *CERFACS/ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2011-03-10 
!        2016-11-25 ANNA AGUSTI-PANAREDA: add carbon tracers (including CO with linco scheme)

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYNA  , ONLY : TDYNA
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM,    JPRB     
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCHEM  , ONLY : TCHEM

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(GEOMETRY),               INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT) :: YDML_GCONF
TYPE(TDYNA),                  INTENT(IN)    :: YDDYNA
TYPE(TCHEM),                  INTENT(INOUT) :: YDCHEM

! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "updcoch.intfb.h"
#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('LINCO_CHEM_INI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV, YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NDGSAG=>YDDIM%NDGSAG, NDGENG=>YDDIM%NDGENG, &
 & NCHEM=>YGFL%NCHEM, &
 & CHEM_SCHEME=>YDCHEM%CHEM_SCHEME, NUCOCH1=>YDCHEM%NUCOCH1, NUCOCH2=>YDCHEM%NUCOCH2, & 
 & NCHEM_LCOCOEF=>YDCHEM%NCHEM_LCOCOEF, CHEM_LCOVERS=>YDCHEM%CHEM_LCOVERS, &
 & LCHEM_LCOMESO=>YDCHEM%LCHEM_LCOMESO, RCHEM_LCOCOEFA1=>YDCHEM%RCHEM_LCOCOEFA1, &
 & LCHEM_LCOCSTCLIM=>YDCHEM%LCHEM_LCOCSTCLIM, RCHEM_LCOTAUTOP=>YDCHEM%RCHEM_LCOTAUTOP, &
 & RCHEM_LCOCLIMTOP=>YDCHEM%RCHEM_LCOCLIMTOP, LCHEM_LCOLIMIT=>YDCHEM%LCHEM_LCOLIMIT)

!
! Preliminary validation checks
!
IF (TRIM(CHEM_SCHEME) == 'linco' .OR. TRIM(CHEM_SCHEME) == "carbontracers" .OR. &
  &  TRIM(CHEM_SCHEME) == 'linco_RnPb' ) THEN 
    WRITE(NULOUT,*) TRIM(CHEM_SCHEME)
    WRITE(NULOUT,*) 'LINCO SCHEME'
    WRITE(NULOUT,*) 'Number of coefficients ', NCHEM_LCOCOEF
    WRITE(NULOUT,*) 'Version for the coefficients ', CHEM_LCOVERS
    WRITE(NULOUT,*) 'LCHEM_LCOMESO    = ',    LCHEM_LCOMESO 
    WRITE(NULOUT,*) 'LCHEM_LCOCSTCLIM = ', LCHEM_LCOCSTCLIM
    WRITE(NULOUT,*) 'RCHEM_LCOTAUTOP  = ',  RCHEM_LCOTAUTOP
    WRITE(NULOUT,*) 'RCHEM_LCOCLIMTOP = ', RCHEM_LCOCLIMTOP
    WRITE(NULOUT,*) 'LCHEM_LCOLIMIT   = ',   LCHEM_LCOLIMIT 
    WRITE(NULOUT,*) 'RCHEM_LCOCOEFA1  = ',  RCHEM_LCOCOEFA1
ELSE
    WRITE(NULOUT,*) 'ABORT : THE CHEMICAL SCHEME ', CHEM_SCHEME
    WRITE(NULOUT,*) 'IS NOT CONSISTENT WITH THE LINCO SCHEM!'
    CALL ABOR1("")
ENDIF 

!
!   SETUP THE LOGICAL UNITS WHERE TO READ THE COEFICIENTS
!
    NUCOCH1    = 116_JPIM  
    NUCOCH2    = 117_JPIM
!
!   ALLOCATE THE COEFICIENTS ARRAY
!
    ALLOCATE(YDCHEM%TCO2DG(NFLEVG*NCHEM_LCOCOEF,NDGSAG:NDGENG))
    IF (.NOT. LCHEM_LCOCSTCLIM) ALLOCATE(YDCHEM%TCOTOP(NDGSAG:NDGENG))

!
!   SETUP CARBON MONOXIDE CHEMISTRY
!
    CALL UPDCOCH(YDGEOMETRY,YDDYNA,YDCHEM,YDML_GCONF%YRRIP)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LINCO_CHEM_INI',1,ZHOOK_HANDLE)
END SUBROUTINE LINCO_CHEM_INI 
