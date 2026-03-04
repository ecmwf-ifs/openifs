! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPINTDYN(YDQTYPE,YDNAMFPINT,YDFPWSTD,KASLB1,KFIELDS,KGPST,KGPEND,KINTER,PBUF,KFPROMA,KFLDBUF,KBLOCK,KFPNUMD_DEP,PROW)

!**** *FPINTDYN*  - FULL-POS horizontal interpolations of dynamic fields

!     PURPOSE.
!     --------
!        Performs horizontal interpolations or buffer transfers.

!**   INTERFACE.
!     ----------
!       *CALL* *FPINTDYN*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KASLB1   : size of interpolation buffer (core + halo)
!         KFIELDS  : number of fields in row
!         KGPST    : first output point in row.
!         KGPEND   : last output point in row.
!         KINTER   : Kind of interpolation for each fullpos field  : 4, 12 or nearest point (0)
!         PBUF     : buffer which contain the fields to interpolate.
!         KFPROMA  : length of the output row.
!         KFLDBUF  : number of fields in input buffer
!         KBLOCK   : current KFPROMA-sized block
!         KFPNUMD_DEP : subdomain index for each interpolation point
!        OUTPUT:
!         PROW     : interpolated fields.

!        IMPLICIT ARGUMENTS
!        --------------------
!          See #include below. 

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!      See below

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad: 28-02-2007 DM-features optimisations in FULL-POS
!      K. Yessad (oct 2010): multi-linear interpolations.
!     -----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMFPC   , ONLY : TNAMFPINT
USE YOMWFPB  , ONLY : TFPWSTD
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     -----------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE(TNAMFPINT)   ,INTENT(IN)    :: YDNAMFPINT
TYPE(TFPWSTD)     ,INTENT(IN)    :: YDFPWSTD
INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDBUF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINTER(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUF(KASLB1*KFLDBUF) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBLOCK
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPNUMD_DEP(KGPEND-KGPST+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PROW(KFPROMA,KFIELDS) 

!     -----------------------------------------------------------------------
!     LLMONO : .TRUE. for monotonic interpolations
!     LLMNMX : .TRUE. to replace interpolations by the search of the 
!              nearest point
!     LLMASK : logical mask : .TRUE. => no interpolation
!     IDER   : Number of "boxes" for each field
!     IDMP   : subdomain number for each field in input buffer (if relevent)

LOGICAL :: LLMONO(KFIELDS), LLMNMX(KFIELDS), LLMASK(KFIELDS)
LOGICAL :: LLML

INTEGER(KIND=JPIM) :: JFLD, IJ, IPTR, INC, JL, JD, IDOM, IPTRBUF, IFPROW
INTEGER(KIND=JPIM) :: ISTB, IENB
INTEGER(KIND=JPIM) :: IDER(KFIELDS)
INTEGER(KIND=JPIM) :: IDMP(KFLDBUF)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "fpint12.intfb.h"
#include "fpint4.intfb.h"

!     -----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPINTDYN',0,ZHOOK_HANDLE)
ASSOCIATE(NFPINDYN=>YDNAMFPINT%NFPINDYN, LFPML_STD=>YDNAMFPINT%LFPML_STD, &
 & LWSTD12=>YDFPWSTD%LWSTD12, LWSTD04=>YDFPWSTD%LWSTD04, WSTD12=>YDFPWSTD%WSTD12, &
 & WSTD04=>YDFPWSTD%WSTD04, WSTDML=>YDFPWSTD%WSTDML, ML0=>YDFPWSTD%ML0)

!     ------------------------------------------------------------------

!*       1. PREPARATIONS
!           ------------

!     Search for fields with missing points / "derivatives"
IDMP(:)=0
IPTR=1
IPTRBUF=1
DO JFLD=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(JFLD)
  INC=YDQTYPE%ILEV(IJ)
  DO JL=1,INC
    IDER(IPTR+JL-1)=YDQTYPE%ISPD(JL,IJ)
    LLMNMX(IPTR+JL-1)=(NFPINDYN == 0).OR.(KINTER(IJ) == 0)
  ENDDO
  DO JL=1,INC
    IDOM=YDQTYPE%ISPD(JL,IJ)
    IF (IDOM > 1) THEN
      DO JD=1, IDOM
        IDMP(IPTRBUF+JD-1)=YDQTYPE%IDMP(JD,JL,IJ)
      ENDDO
    ENDIF
    IPTRBUF=IPTRBUF+IDOM
  ENDDO
  IPTR=IPTR+INC
ENDDO
IF (IPTR-1 /= KFIELDS) CALL ABOR1('FPINTDYN : INTERNAL ERROR IPTR')
IF (IPTRBUF-1 /= KFLDBUF) CALL ABOR1('FPINTDYN : INTERNAL ERROR IPTRBUF')

LLMONO(:)=.FALSE.

!*       2. INTERPOLATIONS
!           --------------

!     Interpolate everything :
LLMASK(:)=.FALSE.
LLML=LFPML_STD
ISTB=KFPROMA*(KBLOCK-1)+1
IENB=KFPROMA*(KBLOCK-1)+KGPEND-KGPST+1
IFPROW=SIZE(ML0,DIM=2)
IF ((NFPINDYN == 12 .OR. NFPINDYN == 0) .AND. LWSTD12) THEN
  IF (KBLOCK==1) WRITE (NULOUT,*) ' FPINTDYN: FPINT12 STD is called'
  CALL FPINT12(KASLB1,IFPROW,KFIELDS,KGPST,KGPEND,KFPROMA,KFLDBUF,IDER,IDMP,LLML,LLMNMX,LLMONO,&
   & KFPNUMD_DEP,ML0(:,:,KBLOCK),WSTD12(:,:,KBLOCK),WSTDML(:,5:16,KBLOCK),LLMASK,PBUF,PROW)  
ELSEIF ((NFPINDYN == 4 .OR. NFPINDYN == 0) .AND. LWSTD04) THEN
  IF (KBLOCK==1) WRITE (NULOUT,*) ' FPINTDYN: FPINT4 STD is called'
  CALL FPINT4(KASLB1,IFPROW,KFIELDS,KGPST,KGPEND,KFPROMA,KFLDBUF,IDER,IDMP,LLML,LLMNMX,&
   & KFPNUMD_DEP,ML0(:,:,KBLOCK),WSTD04(:,:,KBLOCK),WSTDML(:,5:16,KBLOCK),LLMASK,PBUF,PROW)  
ELSE
  CALL ABOR1('FPINTDYN : NO INTERPOLATOR FOUND !')
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPINTDYN',1,ZHOOK_HANDLE)
END SUBROUTINE FPINTDYN


