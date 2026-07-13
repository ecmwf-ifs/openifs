SUBROUTINE SUDIM1S

USE PARKIND1  ,ONLY : JPIM     ,JPRB,    JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMDPHY  , ONLY : NCSS, NGPP     ,NVPD     ,NGCC     ,NVSO     ,&
            &NVSG     ,NVRS     ,NVSF     ,NVDIAG   ,&
            &NLALO    ,NLAT     ,NLON     ,NTILES   ,NSOTY    ,&
            &NGPA     ,NCOOR    ,NTRAC,  NLEV,& 
            &NVLA     ,&! ENDTURA
            &NCOM     ,NVELO    ,NSCLRO,&   !KPP 
            &NVEG     ,NVHILO, NCLIMDAT,NCSNEC, NVTYPES
USE YOMDIM1S , ONLY : NPROMA
USE YOMLUN1S , ONLY : NULOUT   ,NULNAM
USE YOMLOG1S , ONLY : CFOUT    ,NDIMCDF
USE YOMFORC1S, ONLY : JPSTPFC
USE YOEPHY   , ONLY : LEOCML   ,LEFLAKE  ,LEURBAN, LECTESSEL, LECLIM10D

#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUDIM1S * - Setting up of the dimensions of the surface
!                  multi-column model

!     Purpose.
!     --------
!            dimensioning the surface one-column model

!**   Interface.
!     ----------
!        *CALL* *SUDIM1S*

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------


!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the one-column model

!     Author.
!     -------
!        Pedro Viterbo and Jean-Francois Mahfouf  *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-06-21
!        Bart vd Hurk (KNMI) interface for multi-column model

!     ------------------------------------------------------------------

#endif

IMPLICIT NONE

INTEGER(KIND=JPIM) :: NVS0,NVS1,IVDIAGNTL,IVDIAGTL,NDFORC &
                    &,IVSFNVT,IVSFVT,IVSFVTZ0,IVDIAGVT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "namdim1s.h"

IF (LHOOK) CALL DR_HOOK('SUDIM1S',0,ZHOOK_HANDLE)

NCSS=4      ! soil layers 
NCSNEC=1    ! snow layers 
!! Update from namelist 
REWIND(NULNAM)  
READ(NULNAM,NAMDIM)

!     ------------------------------------------------------------------
!*       1. SET DEFAULT VALUES.
!           -------------------


NSOTY=7		! soil types 
NLEV=1		! number of atmospheric levels

NVS0=3		! soil water, soil temp, ice temp
NVS1=1            ! liquid soil water
NVSG=4*NCSNEC+1		! (snow depth, temp, water,density)*NCSNEC + alb
NVRS=2		! skin temp, skin depth
NVLA=7          ! lake prognostic variable ! ENDUTRA 
NVEG=3          ! lai, Biomstr, Biomstr2


IF (LECLIM10D) THEN !number of seasonally varying data in a year
  NCLIMDAT=38
ELSE
  NCLIMDAT=12
ENDIF

IF (LEFLAKE) THEN
  NTILES=9        ! ENDUTRA - ADDED LAKE TILE 
ELSE
  NTILES=8
ENDIF

IF (LEURBAN) THEN
  NTILES=10
ENDIF

NVHILO=2

NVTYPES=20      ! Number of vegetation types

NTRAC=1		! nr of tracers

IF(LEOCML) THEN
  NCOM=34       !KPP vertical layer
ELSE
  NCOM=2
ENDIF
NVELO=2       !KPP 
NSCLRO=2      !KPP

!NVSF=17
NVSF=11 !test
NVSF=NVSF + 6 ! Add 6-component MODIS albedo parameters (ALUVI,ALUVV,ALUVG,ALNII,ALNIV,ALNIG)
NVSF=NVSF+2   ! INCLUDE LAKE DEPTH AND LAKE COVER (and LAIH LAIL removed; added later in IVSFVT)
NVSF=NVSF+1   ! include CO2type from climate files (C3/C4)
NVSF=NVSF+8+9*NCOM      !KPP 
NVSF=NVSF+1 ! +FWET
NVSF=NVSF+1 ! add VFPBLOB - global gaussian reduced index 
NVSF=NVSF+1 ! add MVFCLAKEF  lake + flood 
NVSF=NVSF+1 ! Include urban

IVSFNVT=9 ! yomgpd1s: z0f,albf,itm,geo,z0h,sst,ci,soty,aror, climatology 
IVSFVT=6 ! yomgpd1s: for each vegetation type: tile fraction, LAI, coverage,type,(rsmin),R0VT
!NVSF=IVSFNVT+IVSFVT*NVHILO (IVSFNVT is may be 17 now)

NVSF=NVSF+IVSFVT*NVHILO

IVDIAGNTL=25 !25 in the original CTESSEL?? (yomgpd1s: physics without VDMSL)
IVDIAGTL=5 ! yomgpd1s: for each tile: fluxes (LHF SHF,E, skin temp
IVDIAGVT=7 ! yomgpd1s: Anday, Anfm, respstr,respstr2,biomass_last,
           !           bloss, bgain

NVDIAG=IVDIAGNTL+IVDIAGTL*NTILES+IVDIAGVT*NVHILO

!*       2.  DIMENSION OF GENERICS GP*, GPD, GPCC, AND GF.
!            -------------------------------------------

NGPP=NVS0*NCSS+NVSG+NVRS
NGPP=NGPP+NVLA  
NGPP=NGPP+(NCOM+1)*(NVELO+NSCLRO) !KPP
NGPP=NGPP+NVEG*NVHILO !CTESSEL

NGPA=NGPP+NVS1*NCSS
NVPD=NVSF+NVDIAG
NGCC=4 !ALBEDO LAIL LAIL FWET

!*       3.  GRID SIZE
!            ---------
NLAT=1
NLON=1 
JPSTPFC=150000 ! max forcing dimension (1d site)
NDFORC=0    ! namelist forcing dimension
NCOOR=0		! selected grid point

!*       4.  OpenMP
!            ------

NPROMA=40

REWIND(NULNAM)  
READ(NULNAM,NAMDIM) 

IF(NDFORC /= 0)JPSTPFC=NDFORC
IF (NDIMCDF == 2)THEN
  NLALO=NLAT*NLON
ELSE
  NLALO=NLON
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUDIM1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SUDIM1S
