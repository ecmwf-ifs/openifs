! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUPROCFP(KDISTRIB,KFPDOM,KGP,YDFPUSERGEO,KNUMPROCFP)

!**** *SUPROCFP*  - SET UP MESSAGE PASSINF FOR FULL-POS

!     PURPOSE.
!     --------
!        To initialize control array for the distribution and 
!        to distribute some global fields (DM distribution linked
!        to arrival geometry).

!        Computes the following variables:
!        KNUMPROCFP

!**   INTERFACE.
!     ----------
!       *CALL* *SUPROCFP*

!        EXPLICIT ARGUMENTS
!        ------------------
!         KDISTRIB       : Kind of distribution
!         KFPDOM         : number of subdomains
!         KGP            : total number of output points
!         KNUMPROCFP     : MPI task number for each point on the output grid

!        IMPLICIT ARGUMENTS
!        ------------------
!        See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!       SUPROCFP

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      Original : 98-10-05 from sufpg

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad: 27-Feb-2007 Move old code in SUPROCFP_DEP, adapt to DM-arrival geometry.
!      R. El Khatib  24-Jul-2012 LFPDISTRIB replaced by IDISTRIB
!      R. El Khatib & Tayfun Dalkilic 13-Sep-2012 IDISTRIB=2
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2 (differentiation
!      between interpolation grid and output grid)
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 27-Jul-2016 bugfix for the case NFPBW* > 0 + recode LWIDER_DOM + optimize for Boyd
!      R. El Khatib 16-May-2019 optimize memory access in NGPSET2PE
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMMP0   , ONLY : NPROC, NGPSET2PE, N_REGIONS_NS, N_REGIONS_EW, N_REGIONS
!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)   :: KDISTRIB
INTEGER(KIND=JPIM), INTENT(IN)   :: KFPDOM
INTEGER(KIND=JPIM), INTENT(IN)   :: KGP
TYPE (TFPUSERGEO) ,INTENT(IN)    :: YDFPUSERGEO(KFPDOM)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNUMPROCFP(KGP)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JI
INTEGER(KIND=JPIM) :: JA, JB, IGLOFF, IGL1, IGL2, IOFF, ILAST, IPROC
INTEGER(KIND=JPIM) :: ILOFF, JGL, JLON, IPOINT  ,IRESOL, ILOC
INTEGER(KIND=JPIM) :: IR, IRA, JROC, IST, IEND
INTEGER(KIND=JPIM) :: JDOM

!     Arrays for data transposition :
INTEGER(KIND=JPIM) :: IPTRFRSTLAT(N_REGIONS_NS)
INTEGER(KIND=JPIM) :: IFRSTLAT(N_REGIONS_NS)
INTEGER(KIND=JPIM) :: ILSTLAT(N_REGIONS_NS)
INTEGER(KIND=JPIM) :: IONL(1:YDFPUSERGEO(1)%NLAT+N_REGIONS_NS-1,1:N_REGIONS_EW)
INTEGER(KIND=JPIM) :: ISTA(1:YDFPUSERGEO(1)%NLAT+N_REGIONS_NS-1,1:N_REGIONS_EW)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "trans_inq.h"
#include "etrans_inq.h"

#include "abor1.intfb.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPROCFP',0,ZHOOK_HANDLE)
ASSOCIATE(NLON=>YDFPUSERGEO%NLON, NLAT=>YDFPUSERGEO%NLAT, NFPLUX=>YDFPUSERGEO%NFPLUX, &
 & NFPRESOL=>YDFPUSERGEO%NFPRESOL, NFPSIZEG=>YDFPUSERGEO%NFPSIZEG, &
 & NFPGUX=>YDFPUSERGEO%NFPGUX, NFPBWX=>YDFPUSERGEO%NFPBWX, NFPBWY=>YDFPUSERGEO%NFPBWY, CFPGRID=>YDFPUSERGEO%CFPGRID)

!-----------------------------------------------------------------------

!*     1. COMPUTE CONTROL ARRAY FOR DISTRIBUTION
!         --------------------------------------

IF (KDISTRIB == 1) THEN

!*     1. PRELIMINARY COMPUTATIONS
!         ------------------------

! k.y.: this piece of code (currently commented)
!       assumes that we don't split the different
!       subdomains before distributing.
!       (ex for IFPDOM=2, NPROC=3):
!       ..................********* *********..................
!       1st proc 1st proc 2nd proc  2nd proc 3rd proc 3rd proc
!       dom1     dom1     dom1      dom2     dom2     dom2

!       Is it the right way for NFPSURFEX>1?

!! IR = KGP - NPROC*(KGP/NPROC)
!! IRA = INT(KGP/NPROC)

!! DO JROC=1,IR
!!   IST=1+(JROC-1)*(IRA+1)
!!   IEND=JROC*(IRA+1)
!!   KNUMPROCFP(IST:IEND)=JROC
!! ENDDO

!! DO JROC=IR+1,NPROC
!!   IST=1+IR*(IRA+1)+(JROC-IR-1)*IRA
!!   IEND=IR*(IRA+1)+(JROC-IR)*IRA
!!   KNUMPROCFP(IST:IEND)=JROC
!! ENDDO

! k.y.: alternative:
!       split the different subdomains before distributing.
!       Processor repartition for the NFPRGRG points:
!       (ex for IFPDOM=2, NPROC=3):
!       .........*********......... .........*********.........
!       1st proc 2nd proc 3rd proc  1st proc 2nd proc 3rd proc
!       dom1     dom1     dom1      dom2     dom2     dom2

  IOFF=0
  DO JDOM=1,KFPDOM

    IR = NFPSIZEG(JDOM) - NPROC*(NFPSIZEG(JDOM)/NPROC)
    IRA = INT(NFPSIZEG(JDOM)/NPROC)

    DO JROC=1,IR
      IST=1+(JROC-1)*(IRA+1)
      IEND=JROC*(IRA+1)
      KNUMPROCFP(IOFF+IST:IOFF+IEND)=JROC
    ENDDO

    DO JROC=IR+1,NPROC
      IST=1+IR*(IRA+1)+(JROC-IR-1)*IRA
      IEND=IR*(IRA+1)+(JROC-IR)*IRA
      KNUMPROCFP(IOFF+IST:IOFF+IEND)=JROC
    ENDDO

    IOFF=IOFF+NFPSIZEG(JDOM)

  ENDDO

ELSEIF (KDISTRIB == 2) THEN

  IF (CFPGRID(1)=='LALON') CALL ABOR1 ('SUPROCFP : NOT READY FOR LATLON DOMAIN IN TRANSPOSITION 2')
  IF (SIZE(NFPRESOL) > 1) CALL ABOR1 ('SUPROCFP : NOT READY FOR MULTI-DOMAIN TRANSPOSITION 2')
  IRESOL=NFPRESOL(1)
  IF (CFPGRID(1) /= 'GAUSS') THEN
    CALL ETRANS_INQ(KRESOL=IRESOL,KPTRFRSTLAT=IPTRFRSTLAT,KLSTLAT=ILSTLAT, &
  &  KONL=IONL,KSTA=ISTA,KFRSTLAT=IFRSTLAT)
  ELSE
    CALL TRANS_INQ(KRESOL=IRESOL,KPTRFRSTLAT=IPTRFRSTLAT,KLSTLAT=ILSTLAT, &
  &  KONL=IONL,KSTA=ISTA,KFRSTLAT=IFRSTLAT)
  ENDIF

! Double-loop on target tasks to determine the target distribution
  DO JA=1,N_REGIONS_NS
    IGLOFF=IPTRFRSTLAT(JA)
    IGL1 = IFRSTLAT(JA)
    IGL2 = ILSTLAT(JA)
    IOFF=0
    IF (JA > 1) THEN
      IF ( ILSTLAT(JA-1) == IFRSTLAT(JA) ) THEN
        ILAST=ILSTLAT(JA-1)-1
      ELSE
        ILAST=ILSTLAT(JA-1)
      ENDIF
      DO JI=IFRSTLAT(1),ILAST
        IOFF=IOFF+YDFPUSERGEO(1)%NFPRGRI(JI)
      ENDDO
    ENDIF
    DO JB=1,N_REGIONS(JA)
      IPROC=NGPSET2PE(JB,JA)
      ILOFF=0
      ILOC=0
      DO JGL=IGL1,IGL2
        DO JLON=1,IONL(IGLOFF+JGL-IGL1,JB)
!         IPOINT is the global adress of the local gridpoint in the task IPROC=(JA,JB)
          IPOINT=IOFF+ILOFF+ISTA(IGLOFF+JGL-IGL1,JB)+JLON-1
!         ILOC is the local address of the gridpoint in the task IPROC=(JA,JB)
          ILOC=ILOC+1
          KNUMPROCFP(IPOINT)=IPROC
        ENDDO
        ILOFF=ILOFF+YDFPUSERGEO(1)%NFPRGRI(JGL)
      ENDDO
    ENDDO
  ENDDO

ELSE

  CALL ABOR1('SUPROCFP : INTERNAL ERROR ON KDISTRIB')

ENDIF

! ---------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPROCFP',1,ZHOOK_HANDLE)
END SUBROUTINE SUPROCFP
