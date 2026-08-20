! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GPGRADIENTS(YDGEOMETRY,YDMODEL,YDGMV,YDGFL,PGRADPHY)
use parkind1 , only:&
 & jprb
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL , ONLY : MODEL
USE YOMGMV , ONLY : TGMV
USE YOMGFL , ONLY : TGFL
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(MODEL) ,INTENT(IN) :: YDMODEL
TYPE(TGMV) ,INTENT(IN) :: YDGMV
TYPE(TGFL) ,INTENT(IN) :: YDGFL
REAL(KIND=JPRB) ,INTENT(OUT) :: PGRADPHY(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_MF%YRARPHY%NGRADIENTS,YDGEOMETRY%YRDIM%NGPBLKS)
call abor1("GPGRADIENTS should never be called with this build configuration - EXITING")
END SUBROUTINE GPGRADIENTS
