! (C) Copyright 1989- Meteo-France.

SUBROUTINE WRHFP(KFPXFLD,YDFPGEO,KFPDOM,YDFPOPH,YDVAB,YDFPUSERGEO,KNFG,CDCONF,CDFPFN,KFPCHKDAT,LDPADDING,YDFLDSC,PFPBUF,YDFACTX,PTSTEP)

!**** *WRHFP*  - WRITE OUT THE HORIZONTALLY POST-PROCESSED FIELDS TO
!                ARPEGE/ALADIN FILES

!     PURPOSE.
!     --------
!        Extract the post-processed fields, sort out the output points for
!        each domain, then write records on files.

!**   INTERFACE.
!     ----------
!       *CALL* *WRHFP(...)*

!        EXPLICIT ARGUMENTS
!        --------------------
!        KFPXFLD : maximum number of fields to extract at a time
!        KFPDOM : number of subdomains
!        CDCONF : configuration of work
!        LDPADDING : padding domains wider than the model one

!        IMPLICIT ARGUMENTS
!        --------------------
!         See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION
!     Core data is extracted from all processors 
!              and distributed to NSTROUT processors 
!     (Core data + Extension zone) is packed by NSTROUT processors
!     (Core data + Extension zone) is written out to file by 1 processor

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-07-07

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-03-16 IDATEF(10)=time of previous event; N3DINI for tests
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvments
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      R. El Khatib : 05-02-22 Change interface to SUPPDATE
!      R. El Khatib : 05-03-03 Prints out limited with epsilon(1.)
!      R. El Khatib : 12-Jan-2009 Bugfix on top of the prints above
!      G. Kerdraon: 10-09-02 IDTIME: read date from model if L_READ_MODEL_DATE=.TRUE.
!      P. Marguinaud: 27-May-2010 One output file per NSTROUT proc (NDISTIO==1)
!      P. Marguinaud: 27-May-2010 One output file per NSTROUT proc (NDISTIO(1)==1)
!      D. Degrauwe (Feb 2012) : Boyd periodization
!      P. Marguinaud: 26-Apr-2012 Call IO server
!      P. Marguinaud: 06-Jun-2012 Set date in FA units used for compression
!      P. Marguinaud: 11-Sep-2012 Refactor using IOFLDDESC_MOD, WRGP2FAFP
!      R. El Khatib 20-Aug-2012 GAUXBUF removed and replaced by HFPBUF
!      R. El Khatib 31-Aug-2012 new E-zone management
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!      P. Marguinaud: 08-Feb-2013 Add LDCLOSE argument
!      R. El Khatib 10-Jul-2013 case N3DINI > 1
!      P. Marguinaud: 10-Oct-2013 Use FACTX
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El Khatib 20-Nov-2015 move out computation of norms unless use of IO-server
!      R. El Khatib 10-Dec-2015 KDATEF, CDFPFN
!      R. El Khatib 12-Aug-2016 cleaning
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMMP0   , ONLY : NPROC
USE TYPE_FAOPH, ONLY : TFAOPH
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMFPGEO, ONLY : TFPGEO
USE YOMVERT  , ONLY : TVAB
USE YOMIO_SERV, ONLY : IO_SERV_C001
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE YOMTAG, ONLY : MTAG_MFIO_WRHFP

USE EXTFPSELECT_MOD, ONLY : EXTFPSELECT
USE YOMIO_SERV_HDR,      ONLY : NIO_SERV_HDR_IDOM_FPA
USE YOMIO_SERV_MAP_PLAN, ONLY : IO_SERV_SEND_PLAN
USE MPL_MODULE, ONLY : MPL_BARRIER
USE FACTX_MOD, ONLY : FACTX
USE FA_MOD, ONLY : JD_TST

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFPXFLD
TYPE (TFPGEO),      INTENT(IN) :: YDFPGEO
INTEGER(KIND=JPIM), INTENT(IN) :: KFPDOM
TYPE (TFAOPH),      INTENT(IN) :: YDFPOPH(KFPDOM)
TYPE(TVAB),         INTENT(IN) :: YDVAB
TYPE (TFPUSERGEO),  INTENT(IN) :: YDFPUSERGEO(KFPDOM)
INTEGER(KIND=JPIM), INTENT(IN) :: KNFG
CHARACTER (LEN=1),  INTENT(IN) :: CDCONF 
CHARACTER (LEN=*),  INTENT(IN) :: CDFPFN(KFPDOM)
INTEGER(KIND=JPIM), INTENT(IN) :: KFPCHKDAT
LOGICAL,            INTENT(IN) :: LDPADDING
TYPE (IOFLDDESC),   INTENT(IN) :: YDFLDSC(KNFG)
REAL(KIND=JPRB),    INTENT(IN) :: PFPBUF(YDFPGEO%NFPROMA,KNFG,YDFPGEO%NFPBLOCS)
TYPE (FACTX),       INTENT(INOUT) :: YDFACTX(KFPDOM)
REAL (KIND=JPRB),   INTENT(IN), OPTIONAL :: PTSTEP

#include "io_serv_log.intfb.h"
#include "wrgp2fafp.intfb.h"
#include "io_serv_map_send_part1.intfb.h"
#include "io_serv_map_send_part2.intfb.h"

!     KNFG   : total number of fields
!     IDTIME : time indicator ;
!              1 = date of the starting file + model forecast range

TYPE (IO_SERV_SEND_PLAN)         :: YLIOSMPP

INTEGER(KIND=JPIM) :: JFLDG
INTEGER(KIND=JPIM) :: INPC   ! Number of fields per chunk
INTEGER(KIND=JPIM) :: JFLDG1, JFLDG2
REAL(KIND=JPRB) :: ZTSTEP
REAL (KIND=JPRB), ALLOCATABLE :: ZGPBUFL (:,:)

LOGICAL            :: LLUSE_IOSERV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! PREPARE FOR EXTRACTING DATA AND WRITING OUT TO FILE

IF (LHOOK) CALL DR_HOOK('WRHFP',0,ZHOOK_HANDLE)
ASSOCIATE(NFPRGPL=>YDFPGEO%NFPRGPL, NFPRGPG=>YDFPGEO%NFPRGPG)

LLUSE_IOSERV = IO_SERV_C001%LIO_SERV_WR

IF (LLUSE_IOSERV) THEN

  IF (IO_SERV_C001%LDBUG) CALL IO_SERV_LOG (IO_SERV_C001, 'wrhfp', 1_JPIM) 

  IO_SERV_C001%IDATEF = YDFACTX (1)%IDATEF
  ZTSTEP=0._JPRB
  IF (PRESENT (PTSTEP)) THEN
    ZTSTEP = PTSTEP
  ENDIF
  CALL IO_SERV_MAP_SEND_PART1(IO_SERV_C001, NFPRGPL, NFPRGPG, YDFLDSC, MTAG_MFIO_WRHFP+ICHAR(CDCONF), YLIOSMPP, & 
   & KDOM_TYPE=NIO_SERV_HDR_IDOM_FPA,PIOPROC1=IO_SERV_C001%PIOPROC1_FLP,PIOPROC2=IO_SERV_C001%PIOPROC2_FLP,PTSTEP=ZTSTEP)
  CALL EXTFPSELECT (YDFPGEO,LDPADDING,YDFLDSC,YLIOSMPP%YLBUFD, PFPBUF)
  CALL IO_SERV_MAP_SEND_PART2 (IO_SERV_C001, YLIOSMPP)

  IF (IO_SERV_C001%LDBUG) CALL IO_SERV_LOG (IO_SERV_C001, 'wrhfp', 2_JPIM)

ELSE

  ! Maximum size of fields chunks 
  IF (KFPXFLD <= 0) THEN
    INPC = KNFG
  ELSE
    INPC = KFPXFLD
  ENDIF

  CALL MPL_BARRIER (CDSTRING='WRHFP:')

  DO JFLDG = 1, KNFG, INPC

    JFLDG1 = JFLDG
    JFLDG2 = MIN (JFLDG+INPC-1, KNFG)

    IF (.NOT.ALLOCATED(ZGPBUFL)) THEN
      ALLOCATE (ZGPBUFL (NFPRGPL, JFLDG2-JFLDG1+1))
    ELSEIF( ALLOCATED(ZGPBUFL) .AND. (JFLDG2-JFLDG1+1) /= SIZE(ZGPBUFL,DIM=2)) THEN
      DEALLOCATE(ZGPBUFL)
      ALLOCATE (ZGPBUFL (NFPRGPL, JFLDG2-JFLDG1+1))
    ENDIF
    
    CALL EXTFPSELECT (YDFPGEO,LDPADDING,YDFLDSC, ZGPBUFL, PFPBUF, KFLDGOFF=JFLDG1-1)

    CALL WRGP2FAFP (YDFPGEO,KFPDOM, YDFPOPH, YDVAB, YDFPUSERGEO, ZGPBUFL, YDFLDSC (JFLDG1:JFLDG2), &
     & YDFACTXS=YDFACTX, KSIZEGS=YDFPUSERGEO%NFPSIZEG, KTAG=MTAG_MFIO_WRHFP, &
     & LDNORM=LLUSE_IOSERV, CDFPFN=CDFPFN, KFPCHKDAT=KFPCHKDAT)
    
    IF (((JFLDG+INPC) < KNFG) .AND. (NPROC > 1)) CALL MPL_BARRIER (CDSTRING='WRHFP:')

  
  ENDDO
  DEALLOCATE (ZGPBUFL)

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRHFP',1,ZHOOK_HANDLE)

END SUBROUTINE WRHFP
