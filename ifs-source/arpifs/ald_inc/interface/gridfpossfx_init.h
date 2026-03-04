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

INTERFACE
SUBROUTINE GRIDFPOSSFX_INIT(YDAFN,YDFPGEO_DEP, YDMODEL, YDFPUSERGEO, YDGEOMETRY, CDFILE, PSFXMSK)
USE PARKIND1, ONLY : JPRB
USE TYPE_MODEL         , ONLY : MODEL
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMFPGEO, ONLY : TFPGEO
USE YOMAFN, ONLY : TAFN
TYPE (TAFN),  INTENT(IN) :: YDAFN
TYPE (TFPGEO)     ,INTENT(IN) :: YDFPGEO_DEP
TYPE(MODEL)       ,INTENT(IN) :: YDMODEL
TYPE (TFPUSERGEO),  INTENT(IN) :: YDFPUSERGEO(:)
TYPE(GEOMETRY), INTENT(IN)      :: YDGEOMETRY
CHARACTER (LEN=*), INTENT(IN)   :: CDFILE(:)
REAL(KIND=JPRB)   ,INTENT(INOUT), ALLOCATABLE :: PSFXMSK(:,:,:)
END SUBROUTINE GRIDFPOSSFX_INIT
END INTERFACE
