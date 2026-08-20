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
SUBROUTINE RECV_OASIS3_SFX(YDGEOMETRY,PTIMEC,PTSTEP,KSW,PSW_BANDS,ZZENITH,GOASIS_PUT,&
                         ZEMIS,ZTSURF,ZALBDIR,ZALBSCA,ZTSRAD)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1, ONLY : JPRB, JPIM
IMPLICIT NONE
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
REAL(KIND=JPRB),                           INTENT(IN)   :: PTIMEC    ! time of atmospheric model
REAL(KIND=JPRB),                           INTENT(IN)   :: PTSTEP    ! time-step of atmospheric model
INTEGER(KIND=JPIM),                        INTENT(IN)   :: KSW       ! number of SW bands 
REAL(KIND=JPRB), DIMENSION(KSW),           INTENT(IN)   :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL(KIND=JPRB), DIMENSION(YDGEOMETRY%YRGEM%NGPTOT),        INTENT(IN)   :: ZZENITH
LOGICAL,                                   INTENT(OUT)  :: GOASIS_PUT
REAL(KIND=JPRB), DIMENSION(YDGEOMETRY%YRGEM%NGPTOT),        INTENT(OUT)  :: ZEMIS
REAL(KIND=JPRB), DIMENSION(YDGEOMETRY%YRGEM%NGPTOT),        INTENT(OUT)  :: ZTSURF
REAL(KIND=JPRB), DIMENSION(YDGEOMETRY%YRGEM%NGPTOT,KSW),    INTENT(OUT)  :: ZALBDIR
REAL(KIND=JPRB), DIMENSION(YDGEOMETRY%YRGEM%NGPTOT,KSW),    INTENT(OUT)  :: ZALBSCA
REAL(KIND=JPRB), DIMENSION(YDGEOMETRY%YRGEM%NGPTOT),        INTENT(OUT)  :: ZTSRAD
END SUBROUTINE RECV_OASIS3_SFX
END INTERFACE
