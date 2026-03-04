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
SUBROUTINE AROINI_BUDGET(LDBU_ENABLE, HLUOUT,KULOUT, PTSTEP,KSV,KRR, &
                         HRAD,HDCONV,HSCONV,HTURB,HTURBDIM,HCLOUD,&
                         HMET_ADV_SCHEME, HSV_ADV_SCHEME)
!     ##########################################################################
USE PARKIND1  ,ONLY : JPIM     ,JPRB
LOGICAL :: LDBU_ENABLE
CHARACTER (LEN=*), INTENT(IN) :: HLUOUT
INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT
REAL(KIND=JPRB), INTENT(IN) :: PTSTEP
INTEGER(KIND=JPIM), INTENT(IN) :: KSV
INTEGER(KIND=JPIM), INTENT(IN) :: KRR
CHARACTER (LEN=*), INTENT(IN) :: HRAD
CHARACTER (LEN=*), INTENT(IN) :: HDCONV
CHARACTER (LEN=*), INTENT(IN) :: HSCONV
CHARACTER (LEN=*), INTENT(IN) :: HTURB
CHARACTER (LEN=*), INTENT(IN) :: HTURBDIM
CHARACTER (LEN=*), INTENT(IN) :: HCLOUD
CHARACTER (LEN=*), INTENT(IN) :: HMET_ADV_SCHEME
CHARACTER (LEN=*), INTENT(IN) :: HSV_ADV_SCHEME

END SUBROUTINE AROINI_BUDGET
END INTERFACE
