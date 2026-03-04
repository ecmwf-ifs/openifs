! (C) Copyright 1989- Meteo-France.

SUBROUTINE WRPLFP(YDRQSP,YDRQGP,KFPRESOL,KSPEC2,KSPEC2G,KFPRGPL,KFPROMA,YDGEOMETRY,YDOFN,PFPBUF,PSPEC,CDPREF,PTSTEP)

!**** *WRPLFP*  - WRITES OUT THE PRESSURE LEVEL POST-PROCESSED FIELDS
!                  IN GRIB

!     PURPOSE.
!     --------
!        CODE IN GRIB FORMAT AND WRITE OUT THE PRESSURE LEVEL POST-PROCESSED FIELDS.

!**   INTERFACE.
!     ----------
!        *CALL* *WRPLFP(...)

!        EXPLICIT ARGUMENTS :     
!        --------------------

!           KFPRESOL : pp resolution tag
!           KSPEC2   : pp local number of waves
!           KSPEC2G  : pp global number of waves
!           KFPRGPL  : pp local number of gridpoints
!           KFPROMA  : pp blocking factor

!        IMPLICIT ARGUMENTS :      THE SPECTRAL ARRAY.
!        --------------------      the post-processed array(s)

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
!      NILS WEDI  *ECMWF*
!      ORIGINAL : 97-02-20 (From WRPLPPG)

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-04-10 Enable post-processing of filtered spectra
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      D.Dent     :02-05-15 remove references to extra field 099
!      D. Salmond   : 05-22-02 Fix to stop Orography overwriting of Geopotential
!      R. El Khatib : 03-04-17 Fullpos improvments
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib: 29-Feb-2012 simplified interface to norms computation
!      R. El Khatib: 19-Jul-2012 LFIT* => NFIT*
!      R. El Khatib: 10-Aug-2012 prepare for Fullpos-2
!      R. El Khatib 20-Aug-2012 GAUXBUF removed and replaced by HFPBUF
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMMP0   , ONLY : LSPLITOUT,MYSETV   
USE YOMFP4L      , ONLY : TRQFP
USE TYPE_FPOFN, ONLY : TFPOFN
USE IOSPECE_MOD,  ONLY : IOSPECE_PL_SELECTF,  IOSPECE_PL_SELECTD,  IOSPECE_PL_COUNT
USE IOGRIDUE_MOD, ONLY : IOGRIDUE_PL_SELECTF, IOGRIDUE_PL_SELECTD, IOGRIDUE_PL_COUNT
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE YOMIO_SERV, ONLY : IO_SERV_C001

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TRQFP),  INTENT(IN) :: YDRQSP
TYPE (TRQFP),  INTENT(IN) :: YDRQGP
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPRESOL
INTEGER(KIND=JPIM),INTENT(IN)  :: KSPEC2
INTEGER(KIND=JPIM),INTENT(IN)  :: KSPEC2G
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPRGPL
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPROMA
TYPE(GEOMETRY)  ,INTENT(IN)    :: YDGEOMETRY
TYPE (TFPOFN)   ,INTENT(IN)    :: YDOFN
REAL(KIND=JPRB), INTENT(IN)    :: PFPBUF(:,:,:)
REAL(KIND=JPRB), INTENT(IN)    :: PSPEC(:,:)
CHARACTER(LEN=1),INTENT (IN)   :: CDPREF
REAL(KIND=JPRB), INTENT(IN), OPTIONAL    :: PTSTEP

!     ------------------------------------------------------------------

CHARACTER (LEN = 1) ::  CLMODE

REAL(KIND=JPRB) ,ALLOCATABLE    :: ZPTRFLD(:,:)
REAL(KIND=JPRB) ,ALLOCATABLE    :: ZREAL(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE  :: IGRIBIOSH(:,:),IGRIBIOGG(:,:)

TYPE (IOFLDDESC),   ALLOCATABLE :: YLFLDSC (:)

INTEGER(KIND=JPIM) :: IFNUM_SP, IFNUM_GRUE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "wroutgpgb.intfb.h"
#include "wroutspgb.intfb.h"
#include "wrplfp_io_serv.intfb.h"

#include "fcttim.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRPLFP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF (IO_SERV_C001%NPROC_IO > 0) THEN
  CALL WRPLFP_IO_SERV(YDRQSP,YDRQGP,YDGEOMETRY,KFPRGPL,KFPROMA,KSPEC2G,PFPBUF,PSPEC,CDPREF,PTSTEP)
ELSE

!*       1.    PREPARATIONS.
!              --------------

!      -----------------------------------------------------------

!*       2.    WRITE SPECTRAL DATA ON PRESSURE LEVELS.
!              ---------------------------------------

!*       2.1   EXTRACT AND COPY SPECTRAL FIELDS
!              --------------------------------

IF (YDRQSP%NFIELDG > 0) THEN

  IFNUM_SP = 0
  CALL IOSPECE_PL_COUNT (YDRQSP,IFNUM_SP,CDPREF=CDPREF)

  ALLOCATE (YLFLDSC (IFNUM_SP))
  CALL IOSPECE_PL_SELECTD (YDRQSP,KSPEC2G,YLFLDSC,CDPREF=CDPREF)

  ALLOCATE (ZPTRFLD (KSPEC2, YDRQSP%NFIELDL))
  CALL IOSPECE_PL_SELECTF (PSPEC,ZPTRFLD, PACK (YLFLDSC, MASK=YLFLDSC%IVSET == MYSETV),CDPREF=CDPREF)

! Initialise igribiosh data structure
!   1 - grib code
!   2 - level 
!   3 - bset to send this field
!   4 - local level for sending bset

  ALLOCATE (IGRIBIOSH (4, IFNUM_SP))
  IGRIBIOSH (1,1:IFNUM_SP) = YLFLDSC%IGRIB
  IGRIBIOSH (2,1:IFNUM_SP) = YLFLDSC%ILEVG
  IGRIBIOSH (3,1:IFNUM_SP) = YLFLDSC%IVSET
  IGRIBIOSH (4,1:IFNUM_SP) = YLFLDSC%ILEVL
  DEALLOCATE (YLFLDSC)

!  Gather global fields and write out

  CALL WROUTSPGB(YDGEOMETRY%YRDIM,ZPTRFLD,IFNUM_SP,IGRIBIOSH,CDPREF,KFPRESOL,YDOFN%CSH,PTSTEP)

  DEALLOCATE (ZPTRFLD)
  DEALLOCATE (IGRIBIOSH)

ENDIF

!      -----------------------------------------------------------

!*       3.0   WRITE OUT UPPER AIR GRID POINT FIELDS
!              --------------------------------------

IF (YDRQGP%NFIELDG > 0) THEN

  IF( LSPLITOUT )THEN
    CLMODE = 'w'
  ELSE
    CLMODE = 'a'
  ENDIF

  IFNUM_GRUE = 0
  CALL IOGRIDUE_PL_COUNT(YDRQGP,YDGEOMETRY%YRGEM,KFPRGPL,KFPROMA,IFNUM_GRUE,CDPREF=CDPREF)

  IF(IFNUM_GRUE > 0) THEN

    ALLOCATE (YLFLDSC (IFNUM_GRUE))
    CALL IOGRIDUE_PL_SELECTD(YDRQGP,YDGEOMETRY%YRGEM,KFPRGPL,KFPROMA,YLFLDSC,CDPREF=CDPREF)
   
    ALLOCATE (ZREAL (KFPRGPL, IFNUM_GRUE))
    CALL IOGRIDUE_PL_SELECTF(YDRQGP,YDGEOMETRY%YRGEM,KFPRGPL,KFPROMA,ZREAL, YLFLDSC, PFPBUF, CDPREF=CDPREF)
   
    ALLOCATE (IGRIBIOGG (2, IFNUM_GRUE))
    IGRIBIOGG (1,1:IFNUM_GRUE) = YLFLDSC%IGRIB
    IGRIBIOGG (2,1:IFNUM_GRUE) = YLFLDSC%ILEVG
    DEALLOCATE (YLFLDSC)

    CALL WROUTGPGB(ZREAL,KFPRGPL,IFNUM_GRUE,IGRIBIOGG,CDPREF,CLMODE,KFPRESOL,YDOFN%CUA,PTSTEP)
    DEALLOCATE (IGRIBIOGG)
    DEALLOCATE (ZREAL)

  ENDIF

ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('WRPLFP',1,ZHOOK_HANDLE)

END SUBROUTINE WRPLFP
