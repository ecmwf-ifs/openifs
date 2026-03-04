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
SUBROUTINE AROINI_NSV(KSV,KSV_CHEMBEG, KSV_CHEMEND, KSV_AERBEG, KSV_AEREND, &
                      KSV_DSTBEG, KSV_DSTEND, KSV_DSTDEPBEG, KSV_DSTDEPEND,&
                      KSV_CO2)
USE PARKIND1  ,ONLY : JPIM

INTEGER(KIND=JPIM), INTENT(IN)::KSV,KSV_CHEMBEG, KSV_CHEMEND, KSV_AERBEG,&
                              & KSV_AEREND , KSV_DSTBEG , KSV_DSTEND,&
                              & KSV_DSTDEPBEG, KSV_DSTDEPEND, KSV_CO2

END SUBROUTINE AROINI_NSV
END INTERFACE
