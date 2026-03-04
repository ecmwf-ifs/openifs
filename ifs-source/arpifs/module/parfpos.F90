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

MODULE PARFPOS

USE PARKIND1  ,ONLY : JPIM

USE PARDIM, ONLY : JPMXLE, JPMXGL

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! === basic dimensions for Full POST-PROCESSING ===

!     JPOSDOM : Maximum number of horizontal (sub)domains
!     JPOSLEN : Maximum length of a (sub)domain name
!     JPOSDIR : Maximum length of the path (or prefix) for the output files
!     JPOSLE  : Maximum number of eta levels on the output subdomain
!     JPOSGL  : Maximum number of latitude rows of the output gaussian grid

!     JPOSSCVA: Maximum number of post-processable passive scalars
!     JPOSPHY : Maximum number of physical fields
!     JPOSSGP : Maximum number of surface gridpoint fields
!     JPOSCFU : Maximum number of cumulated fluxes
!     JPOSXFU : Maximum number of instantaneous fluxes
!     JPOS3P  : Maximum number of pp. pressure levels
!     JPOS3H  : Maximum number of pp. height (above orography) levels
!     JPOS3TH : Maximum number of pp. potential temperature levels
!     JPOS3PV : Maximum number of pp. potential vorticity levels
!     JPOS3S  : Maximum number of pp. eta levels
!     JPOS3I  : Maximum number of pp. temperature levels
!     JPOSSCVA: Maximum number of passive scalars
!     JPOSVX2 : Maximum number of free gp/sp upper air fields or extra surface fields
!     JPOSFSU : Maximum number of free gp/sp surface fields
!     JPOS3DF : Maximum number of 3D dynamic fields
!     JPOS2DF : Maximum number of 2D dynamic fields
!     JPOSDYN : Maximum number of dynamic fields
!     JPOSGHGFLX : Maximum number of post-processable greenhouse gases
!     JPOSGHG : Maximum number of post-processable greenhouse gases
!     JPOSAERO: Maximum number of post-processable aerosols
!     JPOSAERODIAG: Maximum number of aerosol diagnostics per aerosol type
!     JPOSAEROCLIM: Maximum number of post-processable 3D input aerosol fields
!     JPOSAERO_WVL_DIAG: Maximum number of aerosol diagnostic wavelengths
!     JPOSAERO_WVL_DIAG_TYPES: Maximum number of aerosol diagnostic types per wavelength
!     JPOSEMIS2D: Maximum number of 2D emission fluxes for composition
!     JPOSEMIS2DAUX: Maximum number of 2D emission auxiliary fields for composition
!     JPOSAOD:  Maximum number of post-processable Aerosol Optical Depths
!     JPOSPM:  Maximum number of post-processable Particulate Matter
!     JPOSERA40 : Maximum number of post-processable ERA40 GFL fields
!     JPOSNOGW  : Maximum number of post-processable noro gwd scheme fields
!     JPOSUVP : Maximum number of post-processable fields from UV processor
!     JPOSAERAOT: Maximum number of post-processable aerosol optical thicknesses
!     JPOSAERLISI_VAR: Maximum number of post-processable lidar simulator variables
!     JPOSAERLISI_WVL: Maximum number of post-processable lidar simulator wavelengths
!     JPOSAEROUT: Maximum number of post-processable output aerosol fields
!     JPOSSFX   : Maximum number of PREP fields
!     JPOSEZDIAG: Maximum number of post-processable diag fields (EZDIAG)
!     JPOSEDRP   : Maximum number of post-processable EDRP turbulence fields

INTEGER(KIND=JPIM), PARAMETER :: JPOSERA40=14
INTEGER(KIND=JPIM), PARAMETER :: JPOSNOGW=2
INTEGER(KIND=JPIM), PARAMETER :: JPOSDOM=15
INTEGER(KIND=JPIM), PARAMETER :: JPOSLEN=20
INTEGER(KIND=JPIM), PARAMETER :: JPOSDIR=180
INTEGER(KIND=JPIM), PARAMETER :: JPOSLE=JPMXLE
INTEGER(KIND=JPIM), PARAMETER :: JPOSGL=JPMXGL

INTEGER(KIND=JPIM), PARAMETER :: JPOSCFU=82
INTEGER(KIND=JPIM), PARAMETER :: JPOSXFU=65
INTEGER(KIND=JPIM), PARAMETER :: JPOS3P=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3H=127
INTEGER(KIND=JPIM), PARAMETER :: JPOS3TH=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3PV=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3I=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3F=75
INTEGER(KIND=JPIM), PARAMETER :: JPOS3S=JPOSLE
INTEGER(KIND=JPIM), PARAMETER :: JPOSSCVA=6
INTEGER(KIND=JPIM), PARAMETER :: JPOSGHG=3
INTEGER(KIND=JPIM), PARAMETER :: JPOSGHGFLX=6
INTEGER(KIND=JPIM), PARAMETER :: JPOSCHEM=150
INTEGER(KIND=JPIM), PARAMETER :: JPOSCHEMFLX=50
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERO=42
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERO2=2*JPOSAERO
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERODIAG=8
INTEGER(KIND=JPIM), PARAMETER :: JPOSAEROCLIM=3
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERO_WVL_DIAG=20
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERO_WVL_DIAG_TYPES=5
INTEGER(KIND=JPIM), PARAMETER :: JPOSEMIS2D=500 ! sync with NPEMIS2D in yomcompo.F90
INTEGER(KIND=JPIM), PARAMETER :: JPOSEMIS2DAUX=10 ! sync with NPEMIS2DAUX in yomcompo.F90
INTEGER(KIND=JPIM), PARAMETER :: JPOSEMIS3D=10 ! sync with NPEMIS3D in yomcompo.F90
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERAOT=3
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERLISI_VAR=6
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERLISI_WVL=3
INTEGER(KIND=JPIM), PARAMETER :: JPOSAERLISI=JPOSAERLISI_VAR*JPOSAERLISI_WVL
INTEGER(KIND=JPIM), PARAMETER :: JPOSAEROUT=38
INTEGER(KIND=JPIM), PARAMETER :: JPOSAOD=9
INTEGER(KIND=JPIM), PARAMETER :: JPOSPM=3
INTEGER(KIND=JPIM), PARAMETER :: JPOSUVP=2
INTEGER(KIND=JPIM), PARAMETER :: JPOSVX2=70
INTEGER(KIND=JPIM), PARAMETER :: JPOSEZDIAG=20
INTEGER(KIND=JPIM), PARAMETER :: JPOSEDRP=2
! 15 + 1 (Soiltype) +2 (lai) +2(runoff)
! + 2 (Non-Hydro) + 4 (Modis-Albedo) + 3 (ocean T)
! 14 = 5 + 2 Anton + 5 Philippe + 2 (direct SW)
! + 7 for Tomas + 6 for land carbon
! + 17 for convective and other height level diagnostics
! + 3 downward clearsky SWand LW and direct beam radiation at sfc
! + 9 lakes (7 prognostics, 2 ancillary, GpB)
! + 17 precip(10),ptype,pfrac,gust,potevap,totcolSLW,fzra,tpr (RForbes)
! + 2 GPP/REC flux adjustment coefficients + 4 lightning (Lopez) + 5 ocean
! + 1 totpreciprate + 20 for snowML (4*5 snowML GA)
INTEGER(KIND=JPIM), PARAMETER :: JPOSFSU=28
INTEGER(KIND=JPIM), PARAMETER :: JPOSSFX=1000
INTEGER(KIND=JPIM), PARAMETER :: &
 & JPOSSGP=180+JPOSVX2+JPOSFSU+JPOSGHGFLX+JPOSCHEM+3*JPOSCHEMFLX+(JPOSAERO*JPOSAERODIAG)+ &
 & (JPOSAERO_WVL_DIAG*JPOSAERO_WVL_DIAG_TYPES)+JPOSEMIS2D+JPOSEMIS2DAUX+ &
 & JPOSAOD+JPOSPM+JPOSUVP+10+5+2+2+3+2+7+6+17+3+1+9+17+2+4+5
INTEGER(KIND=JPIM), PARAMETER :: JPOSPHY=JPOSSGP+JPOSCFU+JPOSXFU+JPOSSFX
INTEGER(KIND=JPIM), PARAMETER :: &
 & JPOS3DF=101+JPOSSCVA+JPOSVX2+JPOSAERO+JPOSGHG+JPOSCHEM+JPOSERA40+JPOSNOGW &
 &+JPOSUVP+JPOSAEROUT+JPOSAERAOT+JPOSAERLISI+JPOSAEROCLIM+JPOSEDRP
INTEGER(KIND=JPIM), PARAMETER :: JPOS2DF=100+JPOSFSU+JPOSAERO2
INTEGER(KIND=JPIM), PARAMETER :: JPOSDYN=JPOS3DF+JPOS2DF+JPOSEZDIAG

! Support for quadratic interpolations : 2 rows needed around each point
INTEGER(KIND=JPIM), PARAMETER :: JPROW=2

! Maximum number of post-processing objects
INTEGER(KIND=JPIM), PARAMETER :: JPOSOBJ=7
!     ------------------------------------------------------------------
END MODULE PARFPOS
