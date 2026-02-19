! (C) Copyright 2008- ECMWF.
! (C) Copyright 2008- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOEAEROP

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAEROP* - OPTICAL PROPERTIES FOR PROGNOSTIC AEROSOLS
!     ------------------------------------------------------------------

!*---  OPTICAL PROPERTIES FOR OPTICAL DEPTH/RADIANCE SIMULATIONS

REAL(KIND=JPRB) :: ALF_BC(20)     , ASY_BC(20)     , OMG_BC(20)     , RALI_BC(20)
REAL(KIND=JPRB) :: ALF_DD(9,20)   , ASY_DD(9,20)   , OMG_DD(9,20)   , RALI_DD(9,20)
REAL(KIND=JPRB) :: ALF_OM(12,20)  , ASY_OM(12,20)  , OMG_OM(12,20)  , RALI_OM(12,20)
REAL(KIND=JPRB) :: ALF_SS(12,20,9), ASY_SS(12,20,9), OMG_SS(12,20,9), RALI_SS(12,20,9)
REAL(KIND=JPRB) :: ALF_SU(12,20)  , ASY_SU(12,20)  , OMG_SU(12,20)  , RALI_SU(12,20) 
REAL(KIND=JPRB) :: ALF_NI(12,20,2), ASY_NI(12,20,2), OMG_NI(12,20,2), RALI_NI(12,20,2) 
REAL(KIND=JPRB) :: ALF_SOA(12,20,2), ASY_SOA(12,20,2), OMG_SOA(12,20,2), RALI_SOA(12,20,2) 
REAL(KIND=JPRB) :: ALF_AM(12,20)  , ASY_AM(12,20)  , OMG_AM(12,20)  , RALI_AM(12,20) 

!*---  OPTICAL PROPERTIES IN THE 14 SPECTRAL INTERVALS OF THE RRTM_SW 
!        AND 16 SPECTRAL INTERVALS OF THE RRTM_LW RADIATION SCHEMES

REAL(KIND=JPRB) :: ALFS_BC(14)     , ASYS_BC(14)     , OMGS_BC(14)
REAL(KIND=JPRB) :: ALFS_DD(3,14)   , ASYS_DD(3,14)   , OMGS_DD(3,14)
REAL(KIND=JPRB) :: ALFS_FA(14)     , ASYS_FA(14)     , OMGS_FA(14)   
REAL(KIND=JPRB) :: ALFS_OM(12,14)  , ASYS_OM(12,14)  , OMGS_OM(12,14)
REAL(KIND=JPRB) :: ALFS_SS(12,14,3), ASYS_SS(12,14,3), OMGS_SS(12,14,3)
REAL(KIND=JPRB) :: ALFS_SU(12,14)  , ASYS_SU(12,14)  , OMGS_SU(12,14)

REAL(KIND=JPRB) :: ALFL_BC(16)     , ASYL_BC(16)     , OMGL_BC(16)
REAL(KIND=JPRB) :: ALFL_DD(3,16)   , ASYL_DD(3,16)   , OMGL_DD(3,16)
REAL(KIND=JPRB) :: ALFL_FA(16)     , ASYL_FA(16)     , OMGL_FA(16)   
REAL(KIND=JPRB) :: ALFL_OM(12,16)  , ASYL_OM(12,16)  , OMGL_OM(12,16)
REAL(KIND=JPRB) :: ALFL_SS(12,16,3), ASYL_SS(12,16,3), OMGL_SS(12,16,3)
REAL(KIND=JPRB) :: ALFL_SU(12,16)  , ASYL_SU(12,16)  , OMGL_SU(12,16)


!     ------------------------------------------------------------------
! 9 refers to up to 9 bins (3 operationally for SS and DD) 
! 20 to 19 SW channel radiances + 1 channel at 10 micron
! 12 to 12 reference RH

! BC is for black carbon
! DD is for desert dust
! FA is for flying ash
! OM is for organic matter
! SS is for sea-salt
! SU is for sulfate
! NI is for nitrate
! AM is for ammonium
! SOA is for Secondary Organics

! ALF is alpha , the mass extinction coefficient   in m2/g
! ASY is g     , the assymetry factor             ND
! OMG is pizero, the single scattering albedo     ND
! RALI is the lidar ratio  ND
!     ------------------------------------------------------------------
END MODULE YOEAEROP

