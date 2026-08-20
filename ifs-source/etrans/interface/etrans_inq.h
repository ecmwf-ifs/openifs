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
SUBROUTINE ETRANS_INQ(KRESOL,KSPEC,KSPEC2,KSPEC2G,KSPEC2MX,KNUMP,&
 & KGPTOT,KGPTOTG,KGPTOTMX,KGPTOTL,&
 & KMYMS,KESM0,KUMPP,KPOSSP,KPTRMS,KALLMS,KDIM0G,&
 & KFRSTLAT,KLSTLAT,KFRSTLOFF,KPTRLAT,&
 & KPTRFRSTLAT,KPTRLSTLAT,KPTRFLOFF,KSTA,KONL,&
 & KULTPP,KPTRLS,&
 & KPRTRW,KMYSETW,KMYSETV,KMY_REGION_NS,KMY_REGION_EW,&
 & LDSPLITLAT,LDLINEAR_GRID,&
 & KSMAX,KMSMAX,KNVALUE,KMVALUE,PLEPINM,KDEF_RESOL,LDLAM,&
 & PMU,PGW,PRPNM,KLEI3,KSPOLEGL,KPMS,KCPL2M,KCPL4M,KPROCM)  

!**** *ETRANS_INQ* - Extract information from the transform package

!     Purpose.
!     --------
!     Interface routine for extracting information from the T.P.

!**   Interface.
!     ----------
!     CALL ETRANS_INQ(...)
!     Explicit arguments : All arguments are optional.
!     -------------------- 
!     KRESOL   - resolution tag for which info is required ,default is the
!                first defined resulution (input)

!                   MULTI-TRANSFORMS MANAGEMENT
!     KDEF_RESOL - number or resolutions defined
!     LDLAM      - .T. if the corresponding resolution is LAM, .F. if it is global

!                   SPECTRAL SPACE
!     KSPEC    - number of complex spectral coefficients on this PE
!     KSPEC2   - 2*KSPEC 
!     KSPEC2G  - global KSPEC2
!     KSPEC2MX - maximun KSPEC2 among all PEs
!     KNUMP    - Number of spectral waves handled by this PE
!     KGPTOT   - Total number of grid columns on this PE
!     KGPTOTG  - Total number of grid columns on the Globe
!     KGPTOTMX - Maximum number of grid columns on any of the PEs
!     KGPTOTL  - Number of grid columns one each PE (dimension N_REGIONS_NS:N_REGIONS_EW)
!     KMYMS    - This PEs spectral zonal wavenumbers
!     KESM0    - Address in a spectral array of (m, n=m)
!     KUMPP    - No. of wave numbers each wave set is responsible for
!     KPOSSP   - Defines partitioning of global spectral fields among PEs
!     KPTRMS   - Pointer to the first wave number of a given a-set
!     KALLMS   - Wave numbers for all wave-set concatenated together 
!                to give all wave numbers in wave-set order
!     KDIM0G   - Defines partitioning of global spectral fields among PEs
!     KSMAX    - spectral truncation - n direction
!     KMSMAX   - spectral truncation - m direction
!     KNVALUE  - n value for each KSPEC2 spectral coeffient
!     KMVALUE  - m value for each KSPEC2 spectral coeffient
!     LDLINEAR_GRID : .TRUE. if the grid is linear

!                 GRIDPOINT SPACE                  
!     KFRSTLAT    - First latitude of each a-set in grid-point space
!     KLSTTLAT    - Last latitude of each a-set in grid-point space
!     KFRSTLOFF   - Offset for first lat of own a-set in grid-point space
!     KPTRLAT     - Pointer to the start of each latitude
!     KPTRFRSTLAT - Pointer to the first latitude of each a-set in 
!                   NSTA and NONL arrays
!     KPTRLSTLAT  - Pointer to the last latitude of each a-set in
!                   NSTA and NONL arrays
!     KPTRFLOFF   - Offset for pointer to the first latitude of own a-set
!                   NSTA and NONL arrays, i.e. nptrfrstlat(myseta)-1
!     KSTA        - Position of first grid column for the latitudes on a 
!                   processor. The information is available for all processors.
!                   The b-sets are distinguished by the last dimension of 
!                   nsta().The latitude band for each a-set is addressed by 
!                   nptrfrstlat(jaset),nptrlstlat(jaset), and 
!                   nptrfloff=nptrfrstlat(myseta) on this processors a-set.
!                   Each split latitude has two entries in nsta(,:) which 
!                   necessitates the rather complex addressing of nsta(,:)
!                   and the overdimensioning of nsta by N_REGIONS_NS.
!     KONL        - Number of grid columns for the latitudes on a processor.
!                   Similar to nsta() in data structure.
!     LDSPLITLAT  - TRUE if latitude is split in grid point space over 
!                   two a-sets

!                FOURIER SPACE
!     KULTPP   - number of latitudes for which each a-set is calculating 
!                the FFT's.
!     KPTRLS   - pointer to first global latitude of each a-set for which
!                it performs the Fourier calculations

!                 LEGENDRE
!     PMU      - sin(Gaussian latitudes)
!     PGW      - Gaussian weights
!     PRPNM    - Legendre polynomials
!     KLEI3    - First dimension of Legendre polynomials
!     KSPOLEGL - Second dimension of Legendre polynomials
!     KPMS     - Adress for legendre polynomial for given M (NSMAX)
!     PLEPINM  - Eigen-values of the inverse Laplace operator

!     Method.
!     -------

!     Externals.  ESET_RESOL - set resolution
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Nmiri       15-Nov-2007 Phasing with TFL 32R3
!        A.Bogatchev   16-Sep-2010 Phasing with TFL 36R4
!        R. El Khatib 08-Aug-2012 KSMAX,KMSMAX,KNVALUE,KMVALUE,PLEPINM,LDLAM,KDEF_RESOL,LDLINEAR_GRID

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KRESOL 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KSPEC 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KSPEC2 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KSPEC2G 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KSPEC2MX 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KNUMP 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KGPTOT 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KGPTOTG 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KGPTOTMX 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KGPTOTL(:,:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KMYMS(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KESM0(0:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KUMPP(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPOSSP(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPTRMS(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KALLMS(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KDIM0G(0:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KFRSTLAT(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KLSTLAT(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KFRSTLOFF 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPTRLAT(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPTRFRSTLAT(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPTRLSTLAT(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPTRFLOFF 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KSTA(:,:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KONL(:,:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KULTPP(:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPTRLS(:) 
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPRTRW
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYSETW
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYSETV
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMY_REGION_NS
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMY_REGION_EW
LOGICAL           ,OPTIONAL,INTENT(INOUT) :: LDSPLITLAT(:) 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PMU(:) 
REAL(KIND=JPRB)   ,OPTIONAL               :: PGW(:) ! Argument NOT used
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PRPNM(:,:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KLEI3 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KSPOLEGL 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPMS(0:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KCPL2M(0:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KCPL4M(0:) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KPROCM(0:) 
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KNVALUE(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMVALUE(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PLEPINM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KDEF_RESOL
LOGICAL           ,OPTIONAL,INTENT(OUT)   :: LDLAM
LOGICAL           ,OPTIONAL,INTENT(OUT)   :: LDLINEAR_GRID

END SUBROUTINE ETRANS_INQ
END INTERFACE
