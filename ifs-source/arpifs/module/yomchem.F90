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

!****  YOMCHEM
!
!     PURPOSE.
!     --------
!       MODULE CONTAINING VARIABLES FOR CHEMISTRY 
!
!     PARAMETER        DESCRIPTION                                  
!     ---------         -----------                                  !      
!     CHEM_SCHEME      Name of chemistry scheme
!
!     REFERENCE.
!     ----------
!                   
!     AUTHOR.
!     -------
!     2009-08-20   J. Flemming

!     MODIFICATIONS.
!     --------------
!     2012-12-01   M. Diamantakis  Move LTRCMFIX_PS to YOM_YGFL to allow use by other GFL fields 
!     2014-12-01   M. Diamantakis  Move all mass fixers to YOM_YGFL to allow use by other 
!                                  GFL fields
!     2015-12-22   S. Massart     Add the linear carbon monoxide (CO) scheme
!     2018-07-06   V. Huijnen     Add options for ASIS solver for MOC-chemistry, following implementation 
!                                 by Philippe Moinat (CERFACS)
!     2018-09      M. Michou       Add key to ARPEGE-Climat 6.3 chemistry scheme  
MODULE YOMCHEM

USE PARKIND1  ,ONLY : JPIM, JPRB

IMPLICIT NONE
SAVE

! minimum value for chemistry tracer  
     REAL(KIND=JPRB), PARAMETER :: REPSCHEM = 1.0E-25_JPRB 

!   pointer to mass budget fields in extra array 

    ! general fields always required by LCHEM_DIA
    INTEGER(KIND=JPIM),PARAMETER:: IEXTR_EM=1, IEXTR_DD=2, IEXTR_WD=3, IEXTR_CH=4, IEXTR_CHTR=5, IEXTR_NG=6, IEXTR_FE=7

    ! additional fields required by LCHEM_DIA if GLOMAP is used
    INTEGER(KIND=JPIM),PARAMETER:: IEXTR_SDM=8, IEXTR_COND=9, IEXTR_NUC=10, IEXTR_RMOD=11

    ! additional fields required by LCHEM_DIAC
    INTEGER(KIND=JPIM),PARAMETER:: IEXTR_ROH=12,IEXTR_PH=13, IEXTR_ROH_TROP=14, IEXTR_CHEMX=15

TYPE :: TCHEM
     CHARACTER(LEN=20)    :: CHEM_SCHEME 
!  mass diagnostics switch for specific chemistry budget  (OH PHOTO OH_TROP PHOTO_trop ) 
     LOGICAL    :: LCHEM_DIAC
!  switch to use online J-rate ('J-ONline') computation for BASCOE chemistry, or alternatively LUT
     LOGICAL    :: LCHEM_BASCOE_JON
!  switch to use heterogeneous chemistry in BASCOE-stratospheric chemistry. Default is true.
!  Only set to false in case of an emergency.
     LOGICAL    :: LCHEM_BASCOE_HETCHEM
!  switch to use JNO2 from tropo (CBO5) in strato (BASCOE)
     LOGICAL    :: LCHEM_EXTENDJNO2
!  MOZART: use a dummy RSF table for testing without huge lookup table
     LOGICAL    :: LMOZART_RSF_DUMMY
!  output of photolysis rates in extra ) 
     LOGICAL    :: LCHEM_JOUT
!  Lightning emissions active  
     LOGICAL    :: LCHEM_LIGHT
!  Lightning NOx emissions tuning factor
     REAL(KIND=JPRB) :: RCHEM_LINOX_SCALING
!  Air craft emissions active  
     LOGICAL    :: LCHEM_ACNOX
!  NOX (1) /NOy (2) Family advection  
     INTEGER(KIND=JPIM) :: KCHEM_NOXADV
! Chemistry wet deposition
     INTEGER(KIND=JPIM)   :: KCHEM_WETDEP
! Dry deposition option switch
    INTEGER(KIND=JPIM)   :: KCHEM_DRYDEP
!    apply both O3 CBM4 and Cariolle in Stratopshere
     LOGICAL    :: LCHEM_REVCHEM
!    apply revised (updated) CB05 chemistry
     LOGICAL    :: LCHEM_STMA
!    use analysed ozone of ic in startosphere 
     LOGICAL    :: LCHEM_ANAO3
!   mass change in SL pressure fixer
    ! 1: after MF 2: before MF 3: troposphere after 4: troposphere before MF        
    REAL(KIND=JPRB)        :: SMASSCOR (300,4)  
!  surpress NOx in Stratopshere by setting it to zero  
    LOGICAL     :: LCHEM_0NOX
!  use CH4 analysis as IC for methane   
    LOGICAL     :: LCHEM_ANACH4
!  use weak relaxation for CH4 (set true when emissions present in CB05)
    LOGICAL     :: LCHEM_WEAK_CH4_RELAXATION
!  apply scvangening in convection, but not removal   
    LOGICAL     :: LCHEM_CONVSCAV
!  apply C-SHAPE NOx LI profile   
    LOGICAL     :: LCHEM_CSHAPE 
!  apply grid scaling in wet deposition   
    LOGICAL     :: LCHEM_WDFR 
!  apply htap setting   
    LOGICAL     :: LCHEM_HTAP 
!  year if pre-historic application (pre 2000)    
    INTEGER(KIND=JPIM) :: KCHEM_YEARPI 
! fixed CH4 value 
    REAL(KIND=JPRB) :: RCH4CONST  
! variable to select the solver type (CB05-chemistry)
    INTEGER(KIND=JPIM) :: KCHEM_SOLVE
! variable use chemistry o3 in radiation (need NGHGRAD=20,  LEPO3RA=true )
    LOGICAL :: LCHEM_O3RAD
!  Chemistry TL/AD
   LOGICAL  :: LCHEM_TL 
! variable use prognostic aersol for photolysis and heterogenous uptake
    LOGICAL :: LCHEM_AEROI
!   arrays with of unused fields in extra arrays for IEXTR_EM, IEXTR_DD and IEXTR_WD 
    INTEGER(KIND=JPIM) :: IEXTR_FREE(3,160) 

!     NUCOCH1:   unit number for reading the carbon monoxide coefficient file 1
!     NUCOCH2:   unit number for reading the carbon monoxide coefficient file 2

    INTEGER(KIND=JPIM) :: NUCOCH1
    INTEGER(KIND=JPIM) :: NUCOCH2

!   number of coeficients for the linco scheme
    INTEGER(KIND=JPIM) :: NCHEM_LCOCOEF
!   version of the coeficient to read
    CHARACTER(LEN=5)   :: CHEM_LCOVERS
!   TRUE to have a mesospheric source
    LOGICAL             :: LCHEM_LCOMESO
!   climatological (True) or derived from the scheme coeficients (False) values for the  source 
    LOGICAL             :: LCHEM_LCOCSTCLIM
!   climatological relaxation time (in days)
    REAL(KIND=JPRB)     :: RCHEM_LCOTAUTOP
!   climatological relaxation value (in kg/kg)
    REAL(KIND=JPRB)     :: RCHEM_LCOCLIMTOP
!   TRUE to limit the CO tendency and concentration
    LOGICAL             :: LCHEM_LCOLIMIT
!   Coeficient to apply to the A1 coeficient of the scheme
    REAL(KIND=JPRB)     :: RCHEM_LCOCOEFA1
!   ASIS solution method
    CHARACTER(LEN=5)   :: CSOLMET_ASIS             ! DGESV, GMRES or GSEID  or others to be added (NAG methods)
    INTEGER(KIND=JPIM) :: M_soliter_ASIS           ! Max number of iterations, iterative method
    REAL(KIND=JPRB)    :: solcv_ASIS               ! Convergence criterion, iterative method

    REAL(KIND=JPRB)    :: RTOL_ASIS                ! Relative Tolerance parameter of ASIS
    REAL(KIND=JPRB)    :: ATOL_ASIS                ! Absolute Tolerance parameter of ASIS

    INTEGER(KIND=JPIM) :: N_F_RTOL_ASIS            ! function applied to ASIS_RTOL to modulate its value
                                                   ! according eg to the altitude or latitude
                                                   ! (0 = no function applied, 1 = vertical profile with RTOL=0.1 at 100hPa,
                                                   !  2 = vertical profile with RTOL=0.05 at 100hPa)
    INTEGER(KIND=JPIM) :: N_F_ATOL_ASIS            ! function applied to ASIS_RTOL to modulate its value
                                                   ! (0 = no function applied, 1 = ATOL=1.E-13*density of air molecules)

!     GRID POINT ARRAYS FOR CARBON MONOXYDE: DM-GLOBAL VERSIONS

    !REAL(KIND=JPRB),ALLOCATABLE:: TCO1D(:)
    REAL(KIND=JPRB),ALLOCATABLE:: TCO2DG(:,:)
    !REAL(KIND=JPRB),ALLOCATABLE:: TCO3DBG(:,:,:)
    REAL(KIND=JPRB),ALLOCATABLE:: TCOTOP(:)

!   apply ARPEGE-Climat chemistry specificities    
    LOGICAL             :: LCHEM_ARPCLIM
    
END TYPE TCHEM

!!TYPE(TCHEM), POINTER :: YRCHEM => NULL()

END MODULE YOMCHEM
