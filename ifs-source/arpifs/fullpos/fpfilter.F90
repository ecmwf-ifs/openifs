! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPFILTER(KFPCMAX,KFPDOM,LDFPFIL,KFMAX,PSEL,PBED,PLTF,LDBED,PFPFIL)

!**** *FPFILTER*  - SET UP POST-PROCESSING SPECTRAL FILTER

!     PURPOSE.
!     --------
!        To initialize the profile of the spectral post-processing filter from a global geometry

!**   INTERFACE.
!     ----------
!       *CALL* *FPFILTER*

!        EXPLICIT ARGUMENTS
!        --------------------
!          KFPCMAX: maximum truncation for the post-processing
!          KFPDOM : number of subdomains
!          LDFPFIL : .TRUE. if filter active (for each domain)
!          KFMAX  : maximum truncation of the output subdomains.
!          PSEL : coefficient of selectivity of the low-pass filter
!          PBED : coefficient of the exponential function in the low-pass filter.
!          PLTF : coefficient of the exponential function in the gaussian filter
!          LDBED : .TRUE. to use the low-pass filter on the homogenous resolution space
!                  space ; otherwise the gaussian filter is used.
!          PFPFIL : filters profile for each domain

!       IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 29-May-2013 from SUFPF

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El khatib 03-03-15 read/write filtering matrixes
!      K. Yessad (Nov 2010): minor cleanings
!    P.Marguinaud : 10-Oct-2013 : Fix INI3WRFP arguments
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFPCMAX
INTEGER(KIND=JPIM), INTENT(IN) :: KFPDOM
LOGICAL,            INTENT(IN) :: LDFPFIL(KFPDOM)
INTEGER(KIND=JPIM), INTENT(IN) :: KFMAX(KFPDOM)
REAL(KIND=JPRB),    INTENT(IN) :: PSEL(KFPDOM)
REAL(KIND=JPRB),    INTENT(IN) :: PBED
REAL(KIND=JPRB),    INTENT(IN) :: PLTF
LOGICAL,            INTENT(IN) :: LDBED
REAL(KIND=JPRB), ALLOCATABLE, INTENT(OUT) :: PFPFIL(:,:)

INTEGER(KIND=JPIM) :: J, JN
REAL(KIND=JPRB) :: ZK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPFILTER',0,ZHOOK_HANDLE)

!*       4.    COMPUTE FILTER
!              --------------

ALLOCATE(PFPFIL(0:KFPCMAX,KFPDOM))

IF (LDBED) THEN
! Low-pass (aka "THX") filter
  DO J=1,KFPDOM
    IF (LDFPFIL(J)) THEN
      ZK=EXP(-PBED*PSEL(J))
      DO JN=0,KFPCMAX
        PFPFIL(JN,J)=0.5_JPRB*(1.0_JPRB-TANH(ZK*REAL(JN-KFMAX(J),JPRB)))
      ENDDO
    ENDIF
  ENDDO
ELSE
! Gaussian filter
  DO J=1,KFPDOM
    IF (LDFPFIL(J)) THEN
      ZK=-0.5_JPRB*PLTF/REAL(MIN(KFPCMAX,KFMAX(J))**2,JPRB)
      DO JN=0,KFPCMAX
        PFPFIL(JN,J)=EXP(ZK*REAL(JN,JPRB)**2)
      ENDDO
    ENDIF
  ENDDO
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPFILTER',1,ZHOOK_HANDLE)
END SUBROUTINE FPFILTER
