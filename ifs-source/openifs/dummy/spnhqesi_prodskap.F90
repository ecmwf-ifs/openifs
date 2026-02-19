! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SPNHQESI_PRODSKAP(LDONEM,YDGEOMETRY,YDDYNA,YDGMV,KDIMV,KDIMH,KSPEC2V,PSPDIVG,PSPPROD)
USE PARKIND1 , ONLY : JPIM, JPRB
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV , ONLY : TGMV
USE YOMDYNA , ONLY : TDYNA
LOGICAL, INTENT(IN) :: LDONEM
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYNA), INTENT(IN) :: YDDYNA
TYPE(TGMV), INTENT(INOUT) :: YDGMV
INTEGER(KIND=JPIM), INTENT(IN) :: KDIMV
INTEGER(KIND=JPIM), INTENT(IN) :: KDIMH
INTEGER(KIND=JPIM), INTENT(IN) :: KSPEC2V
REAL(KIND=JPRB), INTENT(IN) :: PSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB), INTENT(OUT) :: PSPPROD(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
call abor1("spnhqesi_prodskap.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SPNHQESI_PRODSKAP
