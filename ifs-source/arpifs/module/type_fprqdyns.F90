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

MODULE TYPE_FPRQDYNS

! Purpose :
! -------
!    To define the internal type "TYPE_FPRQDYN" for the control of 
!    post-processed dynamic fields for a given vertical post-processing level :

!      - %NFPOSDYN  : number of requested fields for the current scan
!      - %NFPTRDYN  : internal codes of the requested fields for the current scan
!      - %LL    : .TRUE. if the computation of the field is needed
!      - %ISF   : spectral computation :
!                 0 : no spectral fit
!                 1 : spectral fit only
!                 2 : spectral fit + gaussian filter in the model spectral space
!                 3 : spectral fit + gaussian filter in the spectral space of 
!                     homogenous resolution (stretched model only)
!        spectral space (1) or not (0)
!      - %ISKP  : the number of fields to be reserved for the computation in gridpoint space
!      - %ILED  : Information about tensor aspect of the field in spectral space :
!        -j : vector 2nd momentum corresponding to the field nr j
!         0 : "scalar" 
!         1 : vector 1st momentum
!         2 : Vorticity as (U,V)
!         3 : Divergence as (U,V)
!      - %ILEV  : the number of vertical levels
!      - %ILVP  : the vertical levels pointers
!      - %IDOM  : the number of subdomains for each level
!      - %IDMP  : the subdomains pointers for each level
!      - %IGPX  : the data array pointers for gridpoint fields (GAUX)
!      - %IGT1(:)  : the data array pointers for spectrally fitted fields : 
!         (0) : pointer after spectral fit (GT0)
!         (1) : main field pointer before spectral fit (GT1)
!         (2) : secondary field pointer (if any) before spectral fit (GT1)
!      - %ISET : the processor (V-)set number for each level of each field
!      - %ILOC : the location of the field in the local V-set spectral array
!      - %ISPD : the number of arrays the field is spread over (for stretching only)

! NFPGT1 : number of fields out of vertical post-processing
! NFPAUXB: number of vertically post-processed fields to remain unfited (Both 2D and 3D)
! NFPVEC : number of vector fields to be fitted
! NFPVECG: global number of vector fields to be fitted
! NFPSCA : number of scalar fields to be fitted
! NFPSCAG: global number of scalar fields to be fitted
! NFPUVMN: number of pairs of true (ie : U,V not Vor, Div) wind components to be fitted
! NFPSPD : number of derived vertically pp. fields to be fited, for one horizontal subdomain.
! NFPSPB : number of derived vertically pp. fields to be fited, for the max. number of horizontal domains.
! NFPGT0B: number of fields to be horizontally post-processed (Both 2D and 3D)
! NFPDYNB: number of post-processed dynamic fields.
! NFPIVEC: number of spectral vector fields to be transformed in gridpoints
! NFPISCA: number of spectral scalar fields to be transformed in gridpoints
! NFPIVECG: global number of spectral vector fields to be transformed in gridpoints
! NFPISCAG: global number of spectral scalar fields to be transformed in gridpoints

! NB : NFPSPD,NFPSPB,NFPSCA,NFPVEC,NFPISCA,NFPIVEC  are local per V-set.

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------
!    Fullpos technical & users guide.

! Author :
! ------
!    Ryad El Khatib *METEO-FRANCE* thanks to Mike Fisher *ECMWF*

! Modifications :
! -------------
! Original : 2000-08-18
! R. El Khatib : 01-03-28 Redefinition of %ILED
! R. El Khatib : 03-02-05 Split module into fullpos_mix + TYPE_FPRQDYN
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARFPOS   , ONLY : JPOSDYN

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC TYPE_FPRQDYN, ALLOCATE_FPRQDYN, DEALLOCATE_FPRQDYN, INQUIRE_FPRQDYN  

TYPE TYPE_FPRQDYN
INTEGER(KIND=JPIM)              :: NFPOSDYN = 0
INTEGER(KIND=JPIM), ALLOCATABLE :: NFPTRDYN(:)
LOGICAL,            ALLOCATABLE :: LL(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ISF(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ISKP(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ILED(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ILEV(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ILVP(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IDOM(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IDMP(:,:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IGPX(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IGT1(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ISET(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ILOC(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ISPD(:,:)
INTEGER(KIND=JPIM) :: NFPGT1 =0
INTEGER(KIND=JPIM) :: NFPAUXB =0
INTEGER(KIND=JPIM) :: NFPVEC =0
INTEGER(KIND=JPIM) :: NFPSCA =0
INTEGER(KIND=JPIM) :: NFPSPD =0
INTEGER(KIND=JPIM) :: NFPSPB =0
INTEGER(KIND=JPIM) :: NFPGT0B =0
INTEGER(KIND=JPIM) :: NFPDYNB =0
INTEGER(KIND=JPIM) :: NFPSCAG =0
INTEGER(KIND=JPIM) :: NFPVECG =0
INTEGER(KIND=JPIM) :: NFPIVEC =0
INTEGER(KIND=JPIM) :: NFPISCA =0
INTEGER(KIND=JPIM) :: NFPIVECG =0
INTEGER(KIND=JPIM) :: NFPISCAG =0
INTEGER(KIND=JPIM) :: NFPUVMN = 0
INTEGER(KIND=JPIM) :: IFIT = 0
END TYPE TYPE_FPRQDYN

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE ALLOCATE_FPRQDYN(YDTYPE,KDOM,KLEV)

TYPE(TYPE_FPRQDYN), INTENT(INOUT) :: YDTYPE
INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV, KDOM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: ILEV

! KDOM : number of subdomains
! KLEV : number of levels

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQDYNS:ALLOCATE_FPRQDYN',0,ZHOOK_HANDLE)

ILEV=MAX(1,KLEV) ! to enable a request with 2D fields only
YDTYPE%NFPOSDYN = 0
ALLOCATE(YDTYPE%NFPTRDYN(JPOSDYN)) ! allocation to nfposdyn has to be enough
ALLOCATE(YDTYPE%LL(JPOSDYN))
ALLOCATE(YDTYPE%ISF(JPOSDYN))
ALLOCATE(YDTYPE%ISKP(JPOSDYN))
ALLOCATE(YDTYPE%ILED(JPOSDYN))
ALLOCATE(YDTYPE%ILEV(JPOSDYN))
ALLOCATE(YDTYPE%ILVP(ILEV,JPOSDYN))
ALLOCATE(YDTYPE%IDOM(ILEV,JPOSDYN))
ALLOCATE(YDTYPE%IDMP(KDOM,ILEV,JPOSDYN))
ALLOCATE(YDTYPE%IGPX(JPOSDYN))
ALLOCATE(YDTYPE%IGT1(JPOSDYN,0:2))
ALLOCATE(YDTYPE%ISET(ILEV,JPOSDYN))
ALLOCATE(YDTYPE%ILOC(ILEV,JPOSDYN))
ALLOCATE(YDTYPE%ISPD(ILEV,JPOSDYN))

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQDYNS:ALLOCATE_FPRQDYN',1,ZHOOK_HANDLE)

END SUBROUTINE ALLOCATE_FPRQDYN

!-----------------------------------------------------------------------------

SUBROUTINE DEALLOCATE_FPRQDYN(YDTYPE)

TYPE(TYPE_FPRQDYN), INTENT(INOUT) :: YDTYPE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQDYNS:DEALLOCATE_FPRQDYN',0,ZHOOK_HANDLE)

YDTYPE%NFPOSDYN = 0
DEALLOCATE(YDTYPE%NFPTRDYN)
DEALLOCATE(YDTYPE%LL)
DEALLOCATE(YDTYPE%ISF)
DEALLOCATE(YDTYPE%ISKP)
DEALLOCATE(YDTYPE%ILED)
DEALLOCATE(YDTYPE%ILEV)
DEALLOCATE(YDTYPE%ILVP)
DEALLOCATE(YDTYPE%IDOM)
DEALLOCATE(YDTYPE%IDMP)
DEALLOCATE(YDTYPE%IGPX)
DEALLOCATE(YDTYPE%IGT1)
DEALLOCATE(YDTYPE%ISET)
DEALLOCATE(YDTYPE%ILOC)
DEALLOCATE(YDTYPE%ISPD)
YDTYPE%NFPGT1 =0
YDTYPE%NFPAUXB =0
YDTYPE%NFPVEC =0
YDTYPE%NFPSCA =0
YDTYPE%NFPSPD =0
YDTYPE%NFPSPB =0
YDTYPE%NFPGT0B =0
YDTYPE%NFPDYNB =0
YDTYPE%NFPSCAG =0
YDTYPE%NFPVECG =0
YDTYPE%NFPIVEC =0
YDTYPE%NFPISCA =0
YDTYPE%NFPIVECG =0
YDTYPE%NFPISCAG =0
YDTYPE%NFPUVMN = 0
YDTYPE%IFIT = 0

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQDYNS:DEALLOCATE_FPRQDYN',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLOCATE_FPRQDYN

!-----------------------------------------------------------------------------

SUBROUTINE INQUIRE_FPRQDYN(YDTYPE,CDNAME,KULOUT)

! YDTYPE  : structure variable
! CDNAME : structure name
! KULOUT : output logical unit number

TYPE(TYPE_FPRQDYN), INTENT(IN) :: YDTYPE
CHARACTER(LEN=*),   INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),          INTENT(IN) :: KULOUT

! CLFMT  : format for output

CHARACTER(LEN=36) :: CLFMT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TYPE_FPRQDYNS:INQUIRE_FPRQDYN',0,ZHOOK_HANDLE)
CLFMT='(1X,''ARRAY '',A,A7,'' ALLOCATED '',8I8)'

WRITE(KULOUT,CLFMT) CDNAME,'%LL    ',SIZE(YDTYPE%LL  ),SHAPE(YDTYPE%LL  )
WRITE(KULOUT,CLFMT) CDNAME,'%ISF   ',SIZE(YDTYPE%ISF) ,SHAPE(YDTYPE%ISF )
WRITE(KULOUT,CLFMT) CDNAME,'%ISKP  ',SIZE(YDTYPE%ISKP),SHAPE(YDTYPE%ISKP)
WRITE(KULOUT,CLFMT) CDNAME,'%ILED  ',SIZE(YDTYPE%ILED),SHAPE(YDTYPE%ILED)
WRITE(KULOUT,CLFMT) CDNAME,'%ILEV  ',SIZE(YDTYPE%ILEV),SHAPE(YDTYPE%ILEV)
WRITE(KULOUT,CLFMT) CDNAME,'%ILVP  ',SIZE(YDTYPE%ILVP),SHAPE(YDTYPE%ILVP)
WRITE(KULOUT,CLFMT) CDNAME,'%IDOM  ',SIZE(YDTYPE%IDOM),SHAPE(YDTYPE%IDOM)
WRITE(KULOUT,CLFMT) CDNAME,'%IDMP  ',SIZE(YDTYPE%IDMP),SHAPE(YDTYPE%IDMP)
WRITE(KULOUT,CLFMT) CDNAME,'%IGPX  ',SIZE(YDTYPE%IGPX),SHAPE(YDTYPE%IGPX)
WRITE(KULOUT,CLFMT) CDNAME,'%IGT1  ',SIZE(YDTYPE%IGT1),SHAPE(YDTYPE%IGT1)
WRITE(KULOUT,CLFMT) CDNAME,'%ISET  ',SIZE(YDTYPE%ISET),SHAPE(YDTYPE%ISET)
WRITE(KULOUT,CLFMT) CDNAME,'%ILOC  ',SIZE(YDTYPE%ILOC),SHAPE(YDTYPE%ILOC)
WRITE(KULOUT,CLFMT) CDNAME,'%ISPD  ',SIZE(YDTYPE%ISPD),SHAPE(YDTYPE%ISPD)
IF (LHOOK) CALL DR_HOOK('TYPE_FPRQDYNS:INQUIRE_FPRQDYN',1,ZHOOK_HANDLE)

END SUBROUTINE INQUIRE_FPRQDYN

!-----------------------------------------------------------------------------

END MODULE TYPE_FPRQDYNS
