! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE TRANS_INQ(KRESOL,KSPEC,KSPEC2,KSPEC2G,KSPEC2MX,KNUMP,&
                    &KGPTOT,KGPTOTG,KGPTOTMX,KGPTOTL,&
                    &KMYMS,KASM0,KUMPP,KPOSSP,KPTRMS,KALLMS,KDIM0G,&
                    &KFRSTLAT,KLSTLAT,KFRSTLOFF,KPTRLAT,&
                    &KPTRFRSTLAT,KPTRLSTLAT,KPTRFLOFF,KSTA,KONL,&
                    &KULTPP,KPTRLS,KNMENG,&
                    &KPRTRW,KMYSETW,KMYSETV,KMY_REGION_NS,KMY_REGION_EW,&
                    &LDSPLITLAT,&
                    &KSMAX,PLAPIN,KNVALUE,KDEF_RESOL,LDLAM,&
                    &PMU,PGW,PRPNM,KLEI3,KSPOLEGL,KPMS,KDGLU)

!**** *TRANS_INQ* - Extract information from the transform package

!     Purpose.
!     --------
!     Interface routine for extracting information from the T.P.

!**   Interface.
!     ----------
!     CALL TRANS_INQ(...)
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
!     KASM0    - Address in a spectral array of (m, n=m)
!     KUMPP    - No. of wave numbers each wave set is responsible for
!     KPOSSP   - Defines partitioning of global spectral fields among PEs
!     KPTRMS   - Pointer to the first wave number of a given a-set
!     KALLMS   - Wave numbers for all wave-set concatenated together 
!                to give all wave numbers in wave-set order
!     KDIM0G   - Defines partitioning of global spectral fields among PEs
!     KSMAX    - spectral truncation
!     KNVALUE  - n value for each KSPEC2 spectral coeffient

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
!     KNMENG   - associated (with NLOENG) cut-off zonal wavenumber

!                 LEGENDRE
!     PMU      - sin(Gaussian latitudes)
!     PGW      - Gaussian weights
!     PRPNM    - Legendre polynomials
!     KLEI3    - First dimension of Legendre polynomials
!     KSPOLEGL - Second dimension of Legendre polynomials
!     KPMS     - Adress for legendre polynomial for given M (NSMAX)
!     PLAPIN   - Eigen-values of the inverse Laplace operator
!     KDGLU    - Number of active points in an hemisphere for a given wavenumber "m"

!     Method.
!     -------

!     Externals.  SET_RESOL - set resolution
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        M. Hortal : 2001-03-05 Dimensions of the Legendre polynomials
!        R. El Khatib 08-Aug-2012 KSMAX,PLAPIN,KNVALUE,LDLAM,KDEF_RESOL

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD


IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC2
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC2G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC2MX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KNUMP
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOT
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOTG
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOTMX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOTL(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KFRSTLOFF
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRFLOFF

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYMS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KASM0(0:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KUMPP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPOSSP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRMS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KALLMS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KDIM0G(0:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KFRSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KLSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRFRSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRLSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSTA(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KONL(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPRTRW
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYSETW
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYSETV
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMY_REGION_NS
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMY_REGION_EW
LOGICAL   ,OPTIONAL, INTENT(OUT) :: LDSPLITLAT(:)

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KULTPP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRLS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KNMENG(:)

REAL(KIND=JPRD)    ,OPTIONAL, INTENT(OUT) :: PMU(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PGW(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PRPNM(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KLEI3
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPOLEGL
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPMS(0:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KDGLU(0:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PLAPIN(-1:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KNVALUE(:)

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KDEF_RESOL
LOGICAL           ,OPTIONAL,INTENT(OUT)   :: LDLAM

END SUBROUTINE TRANS_INQ






END INTERFACE
