! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SMAP_SCREEN (KSTART,KEND,KPROMA,KLEVS, &
 & yd_smos_phys, yd_smap_tb)
 use parkind1, only:&
 & jpim
 use yomsmos , only:&
 & smos_phys_type,&
 & smos_tb_type
 INTEGER (KIND=JPIM) , INTENT (IN) :: KSTART, KEND, KPROMA, KLEVS
 TYPE (SMOS_PHYS_TYPE), INTENT (INOUT) :: YD_SMOS_PHYS (KPROMA)
 TYPE (SMOS_TB_TYPE) , INTENT (INOUT) :: YD_SMAP_TB (KPROMA)
call abor1("SMAP_SCREEN should never be called with this build configuration - EXITING")
END SUBROUTINE SMAP_SCREEN
