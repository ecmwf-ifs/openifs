! (C) Copyright 1989- Meteo-France.

SUBROUTINE WRSFP(YDRQSP,KFPDOM,YDFPOPH,YDVAB,YDFPUSERGEO,PSPEC,CDFPFN,KFPCHKDAT,YDFLDSC,YDFACTX)

!**** *WRSFP*  - WRITE OUT POST-PROCESSED DYNAMIC FIELDS TO ARPEGE/ALADIN FILE
!                AS SPECTRAL COEFFICIENTS OR GRID-POINT FIELDS FOR THOSE WHICH 
!                ARE NOT SPECTRALLY FITTED.

!     PURPOSE.
!     --------
!        To copy the vertically post-processed fields from post-processing
!        spectral arrays or unfitted upper air array to output files

!**   INTERFACE.
!     ----------
!       *CALL* *WRSFP*

!        EXPLICIT ARGUMENTS
!        ------------------
!        CDFPFN : output filenames

!        IMPLICIT ARGUMENTS
!        ------------------
!        See modules above.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-08-03 R. El Khatib

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-03-16 IDATEF(10)=time of previous event; N3DINI for tests
!      R. El Khatib : 01-04-10 Enable post-processing of filtered spectra
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 04-07-22 SPREORD moved inside PREPACKA
!      R. El Khatib : 05-02-22 Change interface to SUPPDATE
!      R. El Khatib : 05-03-03 Prints out limited with epsilon(1.)
!      G. Kerdraon: 10-09-02 IDTIME: read date from model if L_READ_MODEL_DATE=.TRUE.
!      P.Marguinaud : 28-05-10 Change SUMPIOH interface
!      R. El Khatib: 29-Feb-2012 simplified interface to norms computation
!      R. El Khatib: 20-Aug-2012 trwvtof now uses (e)gath_spec
!      P.Marguinaud: 11-Sep-2012 Cleaning + change PREPACKA interface
!      P.Marguinaud: 08-Feb-2013 Add LDCLOSE argument to WRHFP
!      P.Marguinaud: 10-Oct-2013 Use FACTX & WRGATHFLNM
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El Khatib 20-Nov-2015 move out computation of norms
!      R. El Khatib 10-Dec-2015 KDATEF, CDFPFN
!     -----------------------------------------------------------------------

USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMVERT  , ONLY : TVAB
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TYPE_FAOPH, ONLY : TFAOPH
USE YOMFP4L, ONLY : TRQFP
USE YOMMP0   , ONLY : MYPROC   ,NSTROUT ,NPROC
USE FACTX_MOD, ONLY : FACTX
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE IOCPTDESC_MOD, ONLY : IOCPTDESC
USE YOMTAG, ONLY : MTAG_MFIO_WRSFP

IMPLICIT NONE

TYPE (TRQFP),      INTENT(IN)    :: YDRQSP
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPDOM
TYPE (TFAOPH)     ,INTENT(IN)    :: YDFPOPH(KFPDOM)
TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE (TFPUSERGEO) ,INTENT(IN)    :: YDFPUSERGEO(KFPDOM)
REAL(KIND=JPRB),   INTENT(IN)    :: PSPEC(YDRQSP%NFIELDL,YDFPUSERGEO(1)%NSPEC2)
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDFPFN(KFPDOM)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPCHKDAT
TYPE (IOFLDDESC),  INTENT(IN)    :: YDFLDSC(YDRQSP%NFIELDG)
TYPE (FACTX)      ,INTENT(INOUT) :: YDFACTX(KFPDOM)

!     ILEV    : field level to be written out on file
!     IPTR    : pointer of each subdomain in the global overdimensionned
!               output array
!     IFILE  : kind of file to open
!     IDTIME : time indicator ;
!              1 = date of the starting file + model forecast range
!     INFL   : local number of spectral fields in packing distribution
!              (NSTROUT)
!     INFD   : number of spectral fields for each processor in packing
!              distribution (NSTROUT)
!     IFLDOFF: fields offset after distribution among NSTROUT procs.
!     IOPROC : processor number each field belongs
!     IFLDA  : absolute field pointer to be extracted
!     ZVALCO: field-distributed spectral array of packed data
!     ZSPECG : gathered field-distributed spectral array

INTEGER(KIND=JPIM) :: INFD(NPROC), IFLDOFF(NPROC)
INTEGER (KIND=JPIM) :: IMASK(SIZE(YDFLDSC),KFPDOM)

REAL(KIND=JPRB),     ALLOCATABLE :: ZVALCO(:,:), ZSPECG(:,:)
INTEGER(KIND=JPIM) :: IOPROC(YDRQSP%NFIELDG)

INTEGER(KIND=JPIM) :: IRESOL, INFL, JFLD, JROC, JDOM
INTEGER(KIND=JPIM) :: JFLDG1, JFLDG2
INTEGER(KIND=JPIM) :: IPTR(KFPDOM)

TYPE (IOCPTDESC), ALLOCATABLE :: YLCPDSC (:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "wrgathflnm.intfb.h"
#include "prepacka.intfb.h"
#include "sumpioh.intfb.h"
#include "trwvtof.intfb.h"

!     ------------------------------------------------------------------

!*    1. PREPARATIONS
!        ------------

IF (LHOOK) CALL DR_HOOK('WRSFP',0,ZHOOK_HANDLE)
ASSOCIATE(IVSETG=>YDRQSP%IVSETG, &
 & NFPRESOL=>YDFPUSERGEO%NFPRESOL, NSPEC2G=>YDFPUSERGEO%NSPEC2G, &
 & NSPSIZPK=>YDFPOPH%NSPSIZPK)


IF (KFPDOM > 1) CALL ABOR1('WRSFP : NOT READY FOR MULTI-SPECTRAL GATHERING')

!*      2.1 SETUP CONTIGUS DISTRIBUTION AMONG "I/O" PROCS (<=> INFL>0)

CALL SUMPIOH(KPROC=NPROC,KPRIO=NSTROUT,KFLG=YDRQSP%NFIELDG,KFLD=INFD,KFLDOFF=IFLDOFF)
! We can use ioproc retrived from sumpioh as long as infd(:)/ifldoff are used, assuming
! the distribution is contiguous ; therefore we have to recompute locally ioproc :
IOPROC(:)=0
DO JROC=1,NPROC
  IF (INFD(JROC) > 0) THEN
    IOPROC(IFLDOFF(JROC)+1:IFLDOFF(JROC)+INFD(JROC))=JROC
  ENDIF
ENDDO
DO JFLD=1,YDRQSP%NFIELDG
  IF (IOPROC(JFLD) > NPROC .OR. IOPROC(JFLD) < 1) THEN
    CALL ABOR1('WRSFP : INTERNAL ERROR ON IOPROC')
  ENDIF
ENDDO

INFL=INFD(MYPROC)

!*      2.4 TRANSPOSE THE WAVE-DISTRIBUTED SET OF FIELDS TO THE 
!           FIELD-DISTRIBUTED SET OF SPECTRA

ALLOCATE (ZSPECG (NSPEC2G(1), INFL))

IRESOL=NFPRESOL(1)

CALL TRWVTOF(KRESOL=IRESOL,PSPECG=ZSPECG(:,:),KFGATHG=YDRQSP%NFIELDG,KTO=IOPROC(:),&
   & KVSET=IVSETG,PSPEC=PSPEC(:,:))

!       2.5 IO PROCESSORS PRECONDITION DATA

!*        2.5.3 Setup parameters for packing and field name construction

JFLDG1 = IFLDOFF(MYPROC)+1
JFLDG2 = IFLDOFF(MYPROC)+INFL
IPTR(:)=1

IMASK(:,:)=0
DO JDOM=1,KFPDOM
  DO JFLD=1,SIZE(YDFLDSC)
    IMASK(JFLD,JDOM)=IBITS(YDFLDSC(JFLD)%IFPMASK, JDOM-1, 1)
  ENDDO
ENDDO

!         2.5.4 Pack data & buid field names

ALLOCATE(YLCPDSC (INFL, KFPDOM))
ALLOCATE(ZVALCO(NSPSIZPK(1), INFL))
IF (INFL > 0) THEN
  CALL PREPACKA (YDFACTX, 1, YDFLDSC (JFLDG1:JFLDG2), INFL, IMASK (JFLDG1:JFLDG2,:),&
             & YLCPDSC, IPTR, NSPEC2G, YDFPOPH, PFIELD=ZSPECG, PVALCO=ZVALCO)
ENDIF

DEALLOCATE (ZSPECG)

IF (INFL > 0) THEN
  CALL WRGATHFLNM(YDVAB, 1, INFL, NSPSIZPK(1), INFD, IFLDOFF, ZVALCO, YLCPDSC, YDFPOPH, CDFPFN, KPTR=IPTR,&
               & KMASK=IMASK(JFLDG1:JFLDG2,:), YDFACTXS=YDFACTX, KTAG=MTAG_MFIO_WRSFP, &
               & YDFPUSERGEO=YDFPUSERGEO,KFPCHKDAT=KFPCHKDAT)
ELSE
  CALL WRGATHFLNM(YDVAB, 1, INFL, NSPSIZPK(1), INFD, IFLDOFF, ZVALCO, YLCPDSC, YDFPOPH, CDFPFN, KPTR=IPTR,&
               & YDFACTXS=YDFACTX, KTAG=MTAG_MFIO_WRSFP,YDFPUSERGEO=YDFPUSERGEO,KFPCHKDAT=KFPCHKDAT)
ENDIF

DEALLOCATE (ZVALCO)
DEALLOCATE (YLCPDSC)



!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRSFP',1,ZHOOK_HANDLE)

END SUBROUTINE WRSFP
