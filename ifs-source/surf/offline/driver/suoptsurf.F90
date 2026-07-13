SUBROUTINE SUOPTSURF(KULOUT)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOMLUN1S , ONLY : NULNAM
USE YOMDPHY  , ONLY : NVTYPES

USE YOEPHY   , ONLY : LEOPTSURF, LEFARQUHAR

USE YOEOPTSURF   , ONLY : RVR0VT, RVCMAX25,RHUMREL,RA1,RB1,RG0,RGM25,RE_VCMAX,RE_JMAX,&
                          & OVR0VT,OVCMAX25,OHUMREL,OA1,OB1,OG0,OGM25,OE_VCMAX,OE_JMAX


#ifdef DOC
! (C) Copyright 2021- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUOPTSURF*   - Initialize common YOEOPTSURF optimised parameters in land surface model

!     Purpose.
!     --------
!           Initialize YOEOPTSURF, the common that includes the
!           optimized parameters in land surface model 

!**   Interface.
!     ----------
!        *CALL* *SUOPTSURF(KULOUT) from SU0YOM1S

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOEOPTSURF

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!        or

!        Documentation ARPEGE (depending on which physics will be used)

!     Author.
!     -------
!        A. Agusti-Panareda                    *ECMWF*   02/06/2021
!

!     Modifications.
!     --------------
!
#endif

IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT
INTEGER(KIND=JPIM) :: JSFC, ICO2TYP, ZNOPTDIM

#include "namoptsurf.h"

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!LEOPTSURF=.TRUE. ! read optimized parameters from namelist  NAMOPTSURF
!LEOPTSURF=.FALSE. ! no need to read optimized parameters from namelist  NAMOPTSURF

!ZNOPTDIM=NVTYPES*2
ZNOPTDIM=21

IF (.NOT.ALLOCATED(RVR0VT)) ALLOCATE (RVR0VT(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RVCMAX25)) ALLOCATE (RVCMAX25(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RHUMREL)) ALLOCATE (RHUMREL(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RA1)) ALLOCATE (RA1(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RB1)) ALLOCATE (RB1(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RG0)) ALLOCATE (RG0(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RGM25)) ALLOCATE (RGM25(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RE_VCMAX)) ALLOCATE (RE_VCMAX(NVTYPES,2)) 
IF (.NOT.ALLOCATED(RE_JMAX)) ALLOCATE (RE_JMAX(NVTYPES,2)) 

IF (.NOT.ALLOCATED(OVR0VT)) ALLOCATE (OVR0VT(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OVCMAX25)) ALLOCATE (OVCMAX25(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OHUMREL)) ALLOCATE (OHUMREL(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OA1)) ALLOCATE (OA1(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OB1)) ALLOCATE (OB1(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OG0)) ALLOCATE (OG0(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OGM25)) ALLOCATE (OGM25(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OE_VCMAX)) ALLOCATE (OE_VCMAX(ZNOPTDIM)) 
IF (.NOT.ALLOCATED(OE_JMAX)) ALLOCATE (OE_JMAX(ZNOPTDIM))


OVR0VT(:)=-9999._JPRB 
OVCMAX25(:)=-9999._JPRB 
OHUMREL(:)=-9999._JPRB 
OA1(:)=-9999._JPRB 
OB1(:)=-9999._JPRB 
OG0(:)=-9999._JPRB 
OGM25(:)=-9999._JPRB 
OE_VCMAX(:)=-9999._JPRB 
OE_JMAX(:)=-9999._JPRB 



IF (LEFARQUHAR) THEN

! Optimized values with FLUXNET2015 and BFAS global scaling (2013-2020)
  RVR0VT(1,1) =1.395E-07_JPRB   ! C3 Crops 
  RVR0VT(2,1) =1.332E-07_JPRB   ! C3 Short grass 
  RVR0VT(3,:) =4.617E-07_JPRB   ! Evergreen needle leaf
  RVR0VT(4,:) =4.334E-07_JPRB   ! Decidous  needleleaf trees (no FLUXNET data opt)
  RVR0VT(5,:) =4.060E-07_JPRB   ! Decidious broadleaf 
  RVR0VT(6,:) =3.208E-07_JPRB   ! Evergreen broadleaf
  RVR0VT(7,:) =0.302E-07_JPRB   ! C4 Long grass 
  RVR0VT(8,:) =1.450E-07_JPRB    ! not defined (desert with no vegetation)
  RVR0VT(9,:) =5.230E-07_JPRB   ! Tundra (boreal C3 grass)
  RVR0VT(10,:)=0.557E-07_JPRB  ! Irrigated Crops (C3 crops) (no FLUXNET data opt)
  RVR0VT(11,:)=0.013E-07_JPRB   ! Semidesert C4 grass/bushes (opt params not available) [Default prior value]
  RVR0VT(12,:)=1.450E-07_JPRB    ! not defined (ice, no vegetation)
  RVR0VT(13,:)=3.692E-07_JPRB   ! Bogs and Marshes 
  RVR0VT(14,:)=1.450E-07_JPRB    ! not defined (water, no vegetation)
  RVR0VT(15,:)=1.450E-07_JPRB    ! not defined (ocean, no vegetation)
  RVR0VT(16,:)=1.328E-07_JPRB   ! Evergreen Shrubs      
  RVR0VT(17,:)=0.416E-07_JPRB   ! Deciduous Shrubs
  RVR0VT(18,:)=4.698E-07_JPRB   ! Mixed Forest/woodland
  RVR0VT(19,:)=2.107E-07_JPRB   ! Interrupted Forest (not optimized because not pure PFT) 
  RVR0VT(20,:)=2.700E-07_JPRB    ! water-land mixtures (no vegetation)
  RVR0VT(1,2) =1.038E-07_JPRB   ! C4 Crops
  RVR0VT(2,2) =0.302E-07_JPRB   ! same as PFT7 (long grass)

ELSE 

  !Reference respiration default values (from Boussetta et al., 2013) used in Ags

  RVR0VT(1,:)=0.103E-06_JPRB
  RVR0VT(2,:)=0.080E-06_JPRB
  RVR0VT(3,:)=0.330E-06_JPRB
  RVR0VT(4,:)=0.360E-06_JPRB
  RVR0VT(5,:)=0.280E-06_JPRB
  RVR0VT(6,:)=0.270E-06_JPRB
  RVR0VT(7,:)=0.150E-06_JPRB
  RVR0VT(8,:)=0.145E-06_JPRB!not defined
  RVR0VT(9,:)=0.3595E-06_JPRB
  RVR0VT(10,:)=0.096E-06_JPRB
  RVR0VT(11,:)=0.019E-06_JPRB
  RVR0VT(12,:)=0.145E-06_JPRB!not defined
  RVR0VT(13,:)=0.270E-06_JPRB
  RVR0VT(14,:)=0.145E-06_JPRB!not defined
  RVR0VT(15,:)=0.145E-06_JPRB!not defined
  RVR0VT(16,:)=0.110E-06_JPRB
  RVR0VT(17,:)=0.080E-06_JPRB
  RVR0VT(18,:)=0.420E-06_JPRB
  RVR0VT(19,:)=0.156E-06_JPRB
  RVR0VT(20,:)=0.27E-06_JPRB !water-land mixtures

ENDIF


RVCMAX25(1,1)=49._JPRB    ! Crops (C3 CROPS)
RVCMAX25(2,1)=48.7_JPRB   ! Short grass (C3 GRASS)   
RVCMAX25(3,:)=62.2_JPRB   ! Evergreen needle leaf
RVCMAX25(4,:)=35._JPRB    ! Decidous  needleleaf trees not optimized 
RVCMAX25(5,:)=41.9_JPRB   ! Decidious broadleaf 
!old RVCMAX25(6,:)=33.7_JPRB   ! Evergreen broadleaf 
RVCMAX25(6,:)=33.9_JPRB   ! Evergreen broadleaf  (only using FLUXNET sites from tropics)
RVCMAX25(7,:)=40.3_JPRB   ! C4 Long grass
RVCMAX25(8,:)=50._JPRB    ! Desert (not defined)
RVCMAX25(9,:)=70._JPRB    ! Tundra (boreal C3 GRASS)
RVCMAX25(10,:)=60._JPRB   ! Irrigated Crops (C3 crops) [Default prior value]
RVCMAX25(11,:)=50._JPRB   ! Semidesert C4 grass (opt params not available) [Default prior value]
RVCMAX25(12,:)=50._JPRB   ! Ice caps and glaciers (not defined)
RVCMAX25(13,:)=48.6_JPRB  ! Bogs and Marshes (opt params not available from website) prior
RVCMAX25(14,:)=50._JPRB   ! Inland water (not defined)
RVCMAX25(15,:)=50._JPRB   ! Ocean (not defined)
RVCMAX25(16,:)=58.4_JPRB  ! Evergreen Shrubs
RVCMAX25(17,:)=44.5_JPRB  ! Deciduous Shrubs
RVCMAX25(18,:)=42.2_JPRB  ! Mixed Forest/woodland
RVCMAX25(19,:)=50._JPRB   ! Interrupted Forest (not optimized because not pure PFT) [Default prior value]
RVCMAX25(20,:)=50._JPRB   ! Water land mixtures (not defined)
RVCMAX25(1,2)=70.5        ! C4 Crops
RVCMAX25(2,2)=40.3_JPRB   ! Short grass (C4 GRASS) same as PFT7 (long grass) 

RHUMREL(1,1)=4._JPRB   ! C3 Crops
RHUMREL(2,1)=3.32_JPRB   ! short grass (C3 GRASS)  
RHUMREL(3,:)=2.93_JPRB ! Evergreen needle leaf
RHUMREL(4,:)=1._JPRB   ! Decidious needle leaf trees C3 [default prior value, no FLUXNET data]
RHUMREL(5,:)=3.25_JPRB  ! Decidious broadleaf C3
!old RHUMREL(6,:)=3.73_JPRB ! Evergreen broadleaf C3 (using all FLUXNET sites, most of them outside tropics)
RHUMREL(6,:)=4._JPRB ! Evergreen broadleaf C3 (only using FLUXNET sites from tropics)
RHUMREL(7,:)=3.52_JPRB ! Tall grass C4
RHUMREL(8,:)=1.0_JPRB !  Desert (not defined)
RHUMREL(9,:)=4._JPRB   ! Tundra (C3 GRASS)
RHUMREL(10,:)=1._JPRB  ! Irrigated Crops (C3 GRASS)
RHUMREL(11,:)=1._JPRB  ! Semidesert ! opt par not available in website, user prior
RHUMREL(12,:)=1.0_JPRB !  Ice caps and glaciers (not defined)
RHUMREL(13,:)=4._JPRB   ! Bogs and Marshes (opt param not available from website) use prior
RHUMREL(14,:)=1.0_JPRB !  Inland water (not defined)
RHUMREL(15,:)=1.0_JPRB !  Ocean (not defined)
RHUMREL(16,:)=2.63_JPRB    ! Evergreen Shrubs
RHUMREL(17,:)=4._JPRB    ! Deciduous Shrubs
RHUMREL(18,:)=1.73_JPRB    ! Mixed Forest/woodland
RHUMREL(19,:)=1._JPRB    ! Interrupted Forest [default prior value]
RHUMREL(20,:)=1.0_JPRB !  Water and land mixtures (not defined)
RHUMREL(1,2)=2.06_JPRB   ! C4 Crops
RHUMREL(2,2)=3.52_JPRB   ! short grass (C4 GRASS)  same as PFT7 (long grass) 

RA1(1,1)=0.946_JPRB  ! C3 Crops
RA1(2,1)=0.888_JPRB   ! short grass (C3 GRASS)  
RA1(3,:)=0.68_JPRB ! Evergreen needle leaf
RA1(4,:)=0.85_JPRB   ! Decidious needle leaf trees C3
RA1(5,:)=0.9775_JPRB  ! Decidious broadleaf C3
!old RA1(6,:)=0.68_JPRB ! Evergreen broadleaf C3 (using all FLUXNET sites, most of them outside tropics)
RA1(6,:)=0.9775_JPRB ! Evergreen broadleaf C3 (only using FLUXNET sites from tropics)
RA1(7,:)=0.668_JPRB   ! Tall grass C4
RA1(8,:)=0.85_JPRB     ! Desert (not defined)
RA1(9,:)=0.9775_JPRB     ! Tundra (C3 GRASS)
RA1(10,:)=0.85_JPRB    ! Irrigated Crops [default prior value]
RA1(11,:)=0.72_JPRB    ! Semidesert ! [default prior value]
RA1(12,:)=0.8_JPRB     ! Ice caps and glaciers (not defined)
RA1(13,:)=0.919_JPRB    ! Bogs and Marshes
RA1(14,:)=0.85_JPRB     ! Inland water (not defined)
RA1(15,:)=0.85_JPRB     ! Ocean (not defined)
RA1(16,:)=0.864_JPRB    ! Evergreen Shrubs
RA1(17,:)=0.9775_JPRB    ! Deciduous Shrubs
RA1(18,:)=0.722_JPRB    ! Mixed Forest/woodland
RA1(19,:)=0.85_JPRB    ! Interrupted Forest [default prior value], mixed type
RA1(20,:)=0.85_JPRB    ! Water and land mixtures (not defined)
RA1(1,2)=0.864_JPRB  ! C4 Crops
RA1(2,2)=0.668_JPRB   ! short grass (C3 GRASS) same as PFT7 (long grass) 

RB1(1,1)=0.13_JPRB   ! C3 Crops
RB1(2,1)=0.117_JPRB   ! short grass (C3 GRASS)  
RB1(3,:)=0.168_JPRB ! Evergreen needle leaf
RB1(4,:)=0.14_JPRB   ! Decidious needle leaf trees C3 (not optimized, prior)
RB1(5,:)=0.112_JPRB  ! Decidious broadleaf C3
!old RB1(6,:)=0.168_JPRB ! Evergreen broadleaf C3 (using all FLUXNET sites, most of them outside tropics)
RB1(6,:)=0.138_JPRB ! Evergreen broadleaf C3 (only using FLUXNET sites from tropics)
RB1(7,:)=0.21_JPRB   ! Tall grass C4
RB1(8,:)=0.14_JPRB     ! Desert (not defined)
RB1(9,:)=0.112_JPRB     ! Tundra (C3 GRASS)
RB1(10,:)=0.14_JPRB    ! Irrigated Crops [default prior value]
RB1(11,:)=0.2_JPRB    ! Semidesert [default prior value]
RB1(12,:)=0.14_JPRB     ! Ice caps and glaciers (not defined)
RB1(13,:)=0.126_JPRB    ! Bogs and Marshes
RB1(14,:)=0.14_JPRB     ! Inland water (not defined)
RB1(15,:)=0.14_JPRB     ! Ocean (not defined)
RB1(16,:)=0.206_JPRB    ! Evergreen Shrubs
RB1(17,:)=0.112_JPRB    ! Deciduous Shrubs
RB1(18,:)=0.153_JPRB    ! Mixed Forest/woodland
RB1(19,:)=0.14_JPRB    ! Interrupted Forest [default prior value], mixed type
RB1(20,:)=0.14_JPRB     ! Water and land mixtures (not defined)
RB1(1,2)=0.168_JPRB   ! C4 Crops
RB1(2,2)=0.21_JPRB   ! short grass (C4 GRASS) same as PFT7 (long grass)

RG0(1,1)=0.00699_JPRB ! Crops (C4 CROPS)  NEEDS TO BE CHANGED TO C3!
RG0(2,1)=0.00664_JPRB   ! short grass (C3 GRASS)  
RG0(3,:)=0.00533_JPRB ! Evergreen needle leaf [default prior value, no FLUXNET data]
RG0(4,:)=0.00625_JPRB   ! Decidious needle leaf trees C3 
RG0(5,:)=0.007696_JPRB  ! Decidious broadleaf C3 
!old RG0(6,:)= 0.004375_JPRB ! Evergreen broadleaf C3 (using all FLUXNET sites, most of them outside tropics)
RG0(6,:)=0.007838189 ! Evergreen broadleaf C3 (only using FLUXNET sites from tropics)
RG0(7,:)=0.017_JPRB   ! Tall grass C4
RG0(8,:)=0.00625_JPRB  ! Desert (undefined)
RG0(9,:)=0.008125_JPRB     ! Tundra (C3 GRASS) 
RG0(10,:)=0.00625_JPRB    ! Irrigated Crops (C3 GRASS) [default prior value, no FLUXNET data]
RG0(11,:)=0.01875_JPRB    ! Semidesert [default prior value, no FLUXNET data]
RG0(12,:)=0.00625_JPRB  !  Ice caps and glaciers (undefined)
RG0(13,:)=0.00739_JPRB    ! Bogs and Marshes
RG0(14,:)=0.00625_JPRB  !  Inland waters (undefined)
RG0(15,:)=0.00625_JPRB  !  Ocean (undefined)
RG0(16,:)=0.0142_JPRB    ! Evergreen Shrubs
RG0(17,:)=0.00773_JPRB    ! Deciduous Shrubs
RG0(18,:)=0.00543_JPRB    ! Mixed Forest/woodland
RG0(19,:)=0.00625_JPRB    ! Interrupted Forest [default prior value], mixed type
RG0(20,:)=0.00625_JPRB  !  Water and land mixtures (undefined)
RG0(1,2)=0.0219_JPRB ! Crops (C4 CROPS)  NEEDS TO BE CHANGED TO C3!
RG0(2,2)=0.017_JPRB    ! short grass (C4 GRASS) same as PFT7 (long grass)

RGM25(1,1)=0.513_JPRB   ! C3 Crops
RGM25(2,1)=0.657_JPRB   ! short grass (C3 GRASS)  
RGM25(3,:)=0.1_JPRB ! Evergreen needle leaf
RGM25(4,:)=0.4_JPRB   ! Decidious needle leaf trees C3 [default prior value, no FLUXNET data]
RGM25(5,:)=0.719_JPRB  ! Decidious broadleaf C3
!old RGM25(6,:)=0.55_JPRB ! Evergreen broadleaf C3 (using all FLUXNET sites, most of them outside tropics)
RGM25(6,:)=0.3868752_JPRB ! Evergreen broadleaf C3 (only using FLUXNET sites from tropics)
RGM25(7,:)=0.4_JPRB   ! Tall grass C4 (not used for C4, default value here)
RGM25(8,:)=0.4_JPRB    ! Desert (undefined)
RGM25(9,:)=0.8_JPRB     ! Tundra (C3 GRASS)
RGM25(10,:)=0.4_JPRB    ! Irrigated Crops (C3 GRASS)  [default prior value, no FLUXNET data]
RGM25(11,:)=0.4_JPRB    ! Semidesert   [default prior value, no FLUXNET data]
RGM25(12,:)=0.4_JPRB    ! Ice caps and glaciers  (undefined)
RGM25(13,:)=0.59_JPRB    ! Bogs and Marshes (opt param not available from website) use prior
RGM25(14,:)=0.4_JPRB    ! Inland waters (undefined)
RGM25(15,:)=0.4_JPRB    ! Ocean (undefined)
RGM25(16,:)=0.4_JPRB    ! Evergreen Shrubs    [default prior value, no FLUXNET data]
RGM25(17,:)=0.724_JPRB    ! Deciduous Shrubs
RGM25(18,:)=0.404_JPRB    ! Mixed Forest/woodland
RGM25(19,:)=0.4_JPRB    ! Interrupted Forest [default prior value], mixed type
RGM25(20,:)=0.4_JPRB    ! Water and land mixtures (undefined)
RGM25(1,2)=0.4_JPRB   ! C4 Crops (not used for C4, default value here)
RGM25(2,2)=0.657_JPRB   ! short grass (C3 GRASS) same as PFT7 (long grass)

RE_VCMAX(1,1)=72801.3_JPRB  ! C3 Crops
RE_VCMAX(2,1)=50059.1_JPRB  ! short grass (C3 GRASS)  
RE_VCMAX(3,:)=59880.2_JPRB  ! Evergreen needle leaf
RE_VCMAX(4,:)=71513.0_JPRB   ! Decidious needle leaf trees C3 [default prior value, no FLUXNET data]
RE_VCMAX(5,:)=92966.9_JPRB  ! Decidious broadleaf C3
!old RE_VCMAX(6,:)=50059.1_JPRB  ! Evergreen broadleaf C3 (using all FLUXNET sites, most of them outside tropics)
RE_VCMAX(6,:)=67706.09_JPRB  ! Evergreen broadleaf C3 (only using FLUXNET sites from tropics)
RE_VCMAX(7,:)=47110.0_JPRB   ! Tall grass C4
RE_VCMAX(8,:)=71513.0_JPRB    ! Desert (undefined)
RE_VCMAX(9,:)=50059.1_JPRB  ! Tundra (C3 GRASS) ! not optimized, user prior
RE_VCMAX(10,:)=71513.0_JPRB   ! Irrigated Crops (C3 GRASS) [default prior value, no FLUXNET data]
RE_VCMAX(11,:)=67300.0_JPRB    ! Semidesert  [default prior value, no FLUXNET data]
RE_VCMAX(12,:)=71513.0_JPRB    !Ice caps and glaciers  (undefined)
RE_VCMAX(13,:)=72906.6_JPRB    ! Bogs and Marshes  (opt param not available from website) use prior
RE_VCMAX(14,:)=71513.0_JPRB    ! Inland waters (undefined)
RE_VCMAX(15,:)=71513.0_JPRB    ! Ocean (undefined)
RE_VCMAX(16,:)=47110.0_JPRB    ! Evergreen Shrubs
RE_VCMAX(17,:)=50059.1_JPRB    ! Deciduous Shrubs
RE_VCMAX(18,:)=67899.3_JPRB    ! Mixed Forest/woodland
RE_VCMAX(19,:)=71513.0_JPRB    ! Interrupted Forest [default prior value], mixed type
RE_VCMAX(20,:)=71513.0_JPRB    ! Water and land mixtures (undefined)
RE_VCMAX(1,2)=80090.2_JPRB  ! C4 Crops
RE_VCMAX(2,2)=47110.0_JPRB ! short grass (C4 GRASS)  same as PFT7 (long grass)

RE_JMAX(1,1)=64849.2_JPRB   ! C3 Crops
RE_JMAX(2,1)=41459.5_JPRB   ! short grass (C3 GRASS)  
RE_JMAX(3,:)=35048.2_JPRB   ! Evergreen needle leaf
RE_JMAX(4,:)=49884.0_JPRB    ! Decidious needle leaf trees C3 [default prior value, no FLUXNET data]
RE_JMAX(5,:)=34918.8_JPRB   ! Decidious broadleaf C3
!old RE_JMAX(6,:)=34918.8_JPRB   ! Evergreen broadleaf C3
RE_JMAX(6,:)=34918.8_JPRB   ! Evergreen broadleaf C3 (2 sites)
RE_JMAX(7,:)= 54530.0_JPRB   ! Tall grass C4
RE_JMAX(8,:)=49884.0_JPRB    ! Desert (undefined)
RE_JMAX(9,:)=34918.8_JPRB   ! Tundra (C3 GRASS) 
RE_JMAX(10,:)=49884.0_JPRB    ! Irrigated Crops (C3 GRASS) [default prior value, no FLUXNET data]
RE_JMAX(11,:)=77900.0_JPRB   ! Semidesert [default prior value, no FLUXNET data]
RE_JMAX(12,:)=49884.0_JPRB    !Ice caps and glaciers  (undefined) 
RE_JMAX(13,:)=34918.8_JPRB  ! Bogs and Marshes (opt param not available from website) use prior
RE_JMAX(14,:)=49884.0_JPRB   ! Inland waters (undefined)
RE_JMAX(15,:)=49884.0_JPRB   !  Ocean (undefined)
RE_JMAX(16,:)=54530.0_JPRB   ! Evergreen Shrubs
RE_JMAX(17,:)=52755.6_JPRB  ! Deciduous Shrubs
RE_JMAX(18,:)=34918.8_JPRB  ! Mixed Forest/woodland
RE_JMAX(19,:)=49884.0_JPRB   ! Interrupted Forest[default prior value], mixed type
RE_JMAX(20,:)=49884.0_JPRB   ! Water and land mixtures (undefined)
RE_JMAX(1,2)=101270.0_JPRB   ! C4 Crops
RE_JMAX(2,2)=54530.0_JPRB   ! short grass (C4 GRASS) same as PFT7 (long grass)

IF (LEOPTSURF) THEN
!     ------------------------------------------------------------------

!*       2.    Print final values.
!              -------------------

REWIND(NULNAM)
!READ(NULNAM,NAMPHY)
READ(NULNAM,NAMOPTSURF)


RVR0VT(1,1)=OVR0VT(1)   ! C3 Crops 
RVR0VT(2,1)=OVR0VT(2)   ! C3 short grass 
RVR0VT(3,:)=OVR0VT(3)   ! Evergreen needle leaf
RVR0VT(4,:)=OVR0VT(4)   ! Decidous  needleleaf trees not optimized (no FLUXNET sites)
RVR0VT(5,:)=OVR0VT(5)   ! Decidious broadleaf trees
RVR0VT(6,:)=OVR0VT(6)   ! Evergreen broadleaf trees (mostly in tropics)
RVR0VT(7,:)=OVR0VT(7)   ! C4 Tall grass
RVR0VT(8,:)=OVR0VT(8)   ! Desert
RVR0VT(9,:)=OVR0VT(9)   ! Tundra (C3 grass)
RVR0VT(10,:)=OVR0VT(10) ! Irrigated Crops (C3 grass)
RVR0VT(11,:)=OVR0VT(11) ! Semidesert  (opt params not available from website), use prior
RVR0VT(12,:)=OVR0VT(12) ! Ice caps and glaciers (not defined)
RVR0VT(13,:)=OVR0VT(13) ! Bogs and Marshes (opt params not available from website) prior
RVR0VT(14,:)=OVR0VT(14) ! Inland water (not defined)
RVR0VT(15,:)=OVR0VT(15) ! Ocean (not defined)
RVR0VT(16,:)=OVR0VT(16) ! C4 Evergreen shrubs
RVR0VT(17,:)=OVR0VT(17) ! C3 Deciduous shrubs
RVR0VT(18,:)=OVR0VT(18) ! Mixed forest (wood)
RVR0VT(19,:)=OVR0VT(19) ! Interrupted forest
RVR0VT(20,:)=OVR0VT(20) !water-land mixtures
RVR0VT(1,2)=OVR0VT(21)   ! Crops (C4 CROPS) 
RVR0VT(2,2)=OVR0VT(7)   ! C4 short grass, same as PFT7 (long grasss)

RVCMAX25(1,1)=OVCMAX25(1)    ! C3 Crops 
RVCMAX25(2,1)=OVCMAX25(2)    ! C3 short grass 
RVCMAX25(3,:)=OVCMAX25(3)    ! Evergreen needle leaf
RVCMAX25(4,:)=OVCMAX25(4)    ! Decidous  needleleaf trees not optimized 
RVCMAX25(5,:)=OVCMAX25(5)    ! Decidious broadleaf C3
RVCMAX25(6,:)=OVCMAX25(6)    ! Evergreen broadleaf trees (mostly in tropics)
RVCMAX25(7,:)=OVCMAX25(7)    ! C4 Tall grass
RVCMAX25(8,:)=OVCMAX25(8)    
RVCMAX25(9,:)=OVCMAX25(9)    ! Tundra (C3 GRASS)
RVCMAX25(10,:)=OVCMAX25(10)  ! Irrigated Crops (C3 GRASS)
RVCMAX25(11,:)=OVCMAX25(11)  ! Semidesert  (opt params not available from website) C4 prior
RVCMAX25(12,:)=OVCMAX25(12)  
RVCMAX25(13,:)=OVCMAX25(13)  ! Bogs and Marshes (opt params not available from website) prior
RVCMAX25(14,:)=OVCMAX25(14)  
RVCMAX25(15,:)=OVCMAX25(15)  
RVCMAX25(16,:)=OVCMAX25(16)  ! Evergreen Shrubs
RVCMAX25(17,:)=OVCMAX25(17)  ! Deciduous Shrubs
RVCMAX25(18,:)=OVCMAX25(18)  ! Mixed Forest/woodland
RVCMAX25(19,:)=OVCMAX25(19)  ! Interrupted Forest (not optimized because not pure PFT) prior
RVCMAX25(20,:)=OVCMAX25(20)  
RVCMAX25(1,2)=OVCMAX25(21)   ! Crops (C4 CROPS) 
RVCMAX25(2,2)=OVCMAX25(7)    ! C4 short grass, same as PFT7 (long grasss)

RHUMREL(1,1)=OHUMREL(1)   ! C3 Crops
RHUMREL(2,1)=OHUMREL(2)   ! short grass (C3 GRASS)  
RHUMREL(3,:)=OHUMREL(3) ! Evergreen needle leaf
RHUMREL(4,:)=OHUMREL(4)   ! Decidious needle leaf trees C3
RHUMREL(5,:)=OHUMREL(5)  ! Decidious broadleaf C3
RHUMREL(6,:)=OHUMREL(6) ! Evergreen broadleaf C3
RHUMREL(7,:)=OHUMREL(7) ! Tall grass C4
RHUMREL(8,:)=OHUMREL(8) ! Tall grass C4
RHUMREL(9,:)=OHUMREL(9)   ! Tundra (C3 GRASS)
RHUMREL(10,:)=OHUMREL(10)  ! Irrigated Crops (C3 GRASS)
RHUMREL(11,:)=OHUMREL(11)  ! Semidesert ! opt par not available in website, user prior
RHUMREL(12,:)=OHUMREL(12)  
RHUMREL(13,:)=OHUMREL(13)   ! Bogs and Marshes (opt param not available from website) use prior
RHUMREL(14,:)=OHUMREL(14)
RHUMREL(15,:)=OHUMREL(15)
RHUMREL(16,:)=OHUMREL(16)    ! Evergreen Shrubs
RHUMREL(17,:)=OHUMREL(17)    ! Deciduous Shrubs
RHUMREL(18,:)=OHUMREL(18)    ! Mixed Forest/woodland
RHUMREL(19,:)= OHUMREL(19)   ! Interrupted Forest
RHUMREL(20,:)= OHUMREL(20)  
RHUMREL(1,2)=OHUMREL(21)  ! C4 Crops
RHUMREL(2,2)=OHUMREL(7)   ! short grass (C4 GRASS), same as PFT7 (long grass)   

!DO I=1,OPTDIM 
!   II=MAX(MOD(I,NPFT+1),1)
!   JJ=1+INT(I/(NPFT+1))
!   RHUMREL(II,JJ)=OHUMREL(I)    ! Mixed Forest/woodland

RA1(1,1)=OA1(1)  ! C3 Crops
RA1(2,1)=OA1(2)   ! short grass (C3 GRASS)  
RA1(3,:)=OA1(3) ! Evergreen needle leaf
RA1(4,:)=OA1(4)   ! Decidious needle leaf trees C3
RA1(5,:)=OA1(5)  ! Decidious broadleaf C3
RA1(6,:)=OA1(6) ! Evergreen broadleaf C3
RA1(7,:)=OA1(7) ! Tall grass C4
RA1(8,:)=OA1(8) 
RA1(9,:)=OA1(9)     ! Tundra (C3 GRASS) ! not optimized, user prior
RA1(10,:)=OA1(10)    ! Irrigated Crops (C3 GRASS)
RA1(11,:)=OA1(11)   ! Semidesert ! opt par not available in website, user prior
RA1(12,:)=OA1(12) 
RA1(13,:)=OA1(13)    ! Bogs and Marshes (opt param not available from website) use prior
RA1(14,:)=OA1(14)
RA1(15,:)=OA1(15) 
RA1(16,:)=OA1(16)    ! Evergreen Shrubs
RA1(17,:)=OA1(17)    ! Deciduous Shrubs
RA1(18,:)=OA1(18)    ! Mixed Forest/woodland
RA1(19,:)=OA1(19)    ! Interrupted Forest (not optimized because not pure PFT) prior
RA1(20,:)=OA1(20)    
RA1(1,2)=OA1(21)  ! C4 Crops
RA1(2,2)=OA1(7)   ! short grass (C4 GRASS), same as PFT7 (long grass)   

RB1(1,1)=OB1(1)   ! C3 Crops
RB1(2,1)=OB1(2)   ! short grass (C3 GRASS)  
RB1(3,:)=OB1(3) ! Evergreen needle leaf
RB1(4,:)=OB1(4)   ! Decidious needle leaf trees C3 (not optimized, prior)
RB1(5,:)=OB1(5)  ! Decidious broadleaf C3
RB1(6,:)=OB1(6) ! Evergreen broadleaf C3
RB1(7,:)=OB1(7)   ! Tall grass C4
RB1(8,:)=OB1(8)   
RB1(9,:)=OB1(9)     ! Tundra (C3 GRASS) ! not optimized, user prior
RB1(10,:)=OB1(10)    ! Irrigated Crops (C3 GRASS)
RB1(11,:)=OB1(11)    ! Semidesert ! opt par not available in website, user prior
RB1(12,:)=OB1(12)   
RB1(13,:)=OB1(13)    ! Bogs and Marshes (opt param not available from website) use prior
RB1(14,:)=OB1(14)
RB1(15,:)=OB1(15)   
RB1(16,:)=OB1(16)    ! Evergreen Shrubs
RB1(17,:)=OB1(17)    ! Deciduous Shrubs
RB1(18,:)=OB1(18)   ! Mixed Forest/woodland
RB1(19,:)=OB1(19)    ! Interrupted Forest(not optimized because not pure PFT) prior
RB1(20,:)=OB1(20) 
RB1(1,2)=OB1(21)   ! C4 Crops
RB1(2,2)=OB1(7)   ! short grass (C4 GRASS), same as PFT7 (long grass) 

RG0(1,1)=OG0(1) ! Crops (C3 CROPS) 
RG0(2,1)=OG0(2)   ! short grass (C3 GRASS)  
RG0(3,:)=OG0(3) ! Evergreen needle leaf
RG0(4,:)=OG0(4)  ! Decidious needle leaf trees C3 ! not optimized, user prior
RG0(5,:)=OG0(5)  ! Decidious broadleaf C3
RG0(6,:)=OG0(6) ! Evergreen broadleaf C3
RG0(7,:)=OG0(7)   ! Tall grass C4
RG0(8,:)=OG0(8) 
RG0(9,:)=OG0(9)     ! Tundra (C3 GRASS) ! not optimized, user prior
RG0(10,:)=OG0(10)    ! Irrigated Crops (C3 GRASS)
RG0(11,:)=OG0(11)    ! Semidesert ! opt par not available in website, user prior
RG0(12,:)=OG0(12) 
RG0(13,:)=OG0(13)    ! Bogs and Marshes (opt param not available from website) use prior
RG0(14,:)=OG0(14)
RG0(15,:)=OG0(15) 
RG0(16,:)=OG0(16)    ! Evergreen Shrubs
RG0(17,:)=OG0(17)    ! Deciduous Shrubs
RG0(18,:)=OG0(18)    ! Mixed Forest/woodland
RG0(19,:)=OG0(19) ! Interrupted Forest (not optimized because not pure PFT) prior
RG0(20,:)=OG0(20)
RG0(1,2)=OG0(21) ! Crops (C4 CROPS) 
RG0(2,2)=OG0(7)   ! short grass (C4 GRASS),  same as PFT7 (long grass)

RGM25(1,1)=OGM25(1)   ! C3 Crops
RGM25(2,1)=OGM25(2)   ! short grass (C3 GRASS)  
RGM25(3,:)=OGM25(3) ! Evergreen needle leaf
RGM25(4,:)=OGM25(4)   ! Decidious needle leaf trees C3 ! not optimized, user prior
RGM25(5,:)=OGM25(5)  ! Decidious broadleaf C3
RGM25(6,:)=OGM25(6) ! Evergreen broadleaf C3
RGM25(7,:)=OGM25(7)   ! Tall grass C4
RGM25(8,:)=OGM25(8)  
RGM25(9,:)=OGM25(9)     ! Tundra (C3 GRASS) ! not optimized, user prior
RGM25(10,:)=OGM25(10)    ! Irrigated Crops (C3 GRASS)
RGM25(11,:)=OGM25(11)    ! Semidesert (opt param not available from website) use prior
RGM25(12,:)=OGM25(12)    
RGM25(13,:)=OGM25(13)    ! Bogs and Marshes (opt param not available from website) use prior
RGM25(14,:)=OGM25(14)    
RGM25(15,:)=OGM25(15)    
RGM25(16,:)=OGM25(16)    ! Evergreen Shrubs (currently C4, but needs to be set to C3)?? 
RGM25(17,:)=OGM25(17)    ! Deciduous Shrubs
RGM25(18,:)=OGM25(18)    ! Mixed Forest/woodland
RGM25(19,:)=OGM25(19)    ! Interrupted Forest(not optimized because not pure PFT) prior
RGM25(20,:)=OGM25(20)   
RGM25(1,2)=OGM25(21)   ! C4 Crops
RGM25(2,2)=OGM25(7)   ! short grass (C3 GRASS), same as PFT7 (long grass)

RE_VCMAX(1,1)=OE_VCMAX(1)   ! C3 Crops
RE_VCMAX(2,1)=OE_VCMAX(2)
RE_VCMAX(3,:)=OE_VCMAX(3)
RE_VCMAX(4,:)=OE_VCMAX(4)
RE_VCMAX(5,:)=OE_VCMAX(5)
RE_VCMAX(6,:)=OE_VCMAX(6)
RE_VCMAX(7,:)=OE_VCMAX(7)
RE_VCMAX(8,:)=OE_VCMAX(8)
RE_VCMAX(9,:)=OE_VCMAX(9)
RE_VCMAX(10,:)=OE_VCMAX(10)
RE_VCMAX(11,:)=OE_VCMAX(11)
RE_VCMAX(12,:)=OE_VCMAX(12)
RE_VCMAX(13,:)=OE_VCMAX(13)
RE_VCMAX(14,:)=OE_VCMAX(14)
RE_VCMAX(15,:)=OE_VCMAX(15)
RE_VCMAX(16,:)=OE_VCMAX(16)
RE_VCMAX(17,:)=OE_VCMAX(17)
RE_VCMAX(18,:)=OE_VCMAX(18)
RE_VCMAX(19,:)=OE_VCMAX(19)
RE_VCMAX(20,:)=OE_VCMAX(20)
RE_VCMAX(1,2)=OE_VCMAX(21)
RE_VCMAX(2,2)=OE_VCMAX(7)

RE_JMAX(1,1)=OE_JMAX(1)
RE_JMAX(2,1)=OE_JMAX(2)
RE_JMAX(3,:)=OE_JMAX(3)
RE_JMAX(4,:)=OE_JMAX(4)
RE_JMAX(5,:)=OE_JMAX(5)
RE_JMAX(6,:)=OE_JMAX(6)
RE_JMAX(7,:)=OE_JMAX(7)
RE_JMAX(8,:)=OE_JMAX(8)
RE_JMAX(9,:)=OE_JMAX(9)
RE_JMAX(10,:)=OE_JMAX(10)
RE_JMAX(11,:)=OE_JMAX(11)
RE_JMAX(12,:)=OE_JMAX(12)
RE_JMAX(13,:)=OE_JMAX(13)
RE_JMAX(14,:)=OE_JMAX(14)
RE_JMAX(15,:)=OE_JMAX(15)
RE_JMAX(16,:)=OE_JMAX(16)
RE_JMAX(17,:)=OE_JMAX(17)
RE_JMAX(18,:)=OE_JMAX(18)
RE_JMAX(19,:)=OE_JMAX(19)
RE_JMAX(20,:)=OE_JMAX(20)
RE_JMAX(1,2)=OE_JMAX(21)
RE_JMAX(2,2)=OE_JMAX(7)

ENDIF !LEOPTSURF

WRITE(UNIT=KULOUT,FMT='('' COMMON YOEOPTSURF '')')

DO JSFC=1,20
 DO ICO2TYP=1,2

  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RVR0VT = ',RVR0VT(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RVCMAX25 = ',RVCMAX25(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RHUMREL = ',RHUMREL(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RA1 = ', RA1(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RB1 = ', RB1(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RG0 = ', RG0(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RGM25 = ', RGM25(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' E_VCMAX = ', RE_VCMAX(JSFC,ICO2TYP)
  WRITE(UNIT=KULOUT,FMT=*) 'PFT=',JSFC,' ICO2TYP = ',ICO2TYP,' RE_JMAX = ', RE_JMAX(JSFC,ICO2TYP)

 ENDDO
ENDDO 


RETURN
END SUBROUTINE SUOPTSURF
