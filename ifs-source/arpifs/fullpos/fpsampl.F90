! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPSAMPL(YDNAMFPINT,YDFPWSTD,KASLB1,KFIELDS,KGPST,KGPEND,KFPROMA,KFLDBUF,PBUF,PROW,LDONE,&
 & LDSAMPL,KBLOCK)

!**** *FPSAMPL*  - FULL-POS horizontal sampling

!     PURPOSE.
!     --------
!        To sample input buffer (when the output grid is a subset of the 
!        input grid or the same grid)
!        Numbering of the neighbouring points (I is the target point):

!                      13      5       6      14

!                      7       1       2       8
!                                 (I)
!                      9       3       4      10

!                      15     11      12      16

!**   INTERFACE.
!     ----------
!       *CALL* *FPSAMPL*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KASLB1   : size of interpolation buffer (core + halo)
!         KFIELDS  : number of fields in row
!         KGPST    : first output point in row.
!         KGPEND   : last output point in row.
!         PBUF     : buffer which contain the fields to interpolate.
!         KFPROMA  : length of the output row.
!         KFLDBUF  : number of fields in input buffer
!         LDONE    : Logical mask : signal fields (.TRUE.) already computed
!         KBLOCK   : current KFPROMA-sized block
!        OUTPUT:
!         PROW     : interpolated fields.
!         LDSAMPL  : .TRUE. if sampling has been done

!        IMPLICIT ARGUMENTS
!        --------------------
!          See #include below. 

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!      See below.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-09-09

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib  20-May-2005 NFPWIDE moved to YOMWFPDS
!      K. Yessad: 28-02-2007 DM-features optimisations in FULL-POS
!      K. Yessad (oct 2010): multi-linear interpolations.
!     --------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN    ,ONLY : NULOUT

USE YOMFPC, ONLY : TNAMFPINT
USE YOMWFPB, ONLY : TFPWSTD

!     --------------------------------------------------------------------

IMPLICIT NONE

TYPE (TNAMFPINT),  INTENT(IN) :: YDNAMFPINT
TYPE (TFPWSTD)    ,INTENT(IN)    :: YDFPWSTD
INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDBUF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUF(KASLB1*KFLDBUF) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PROW(KFPROMA,KFIELDS) 
LOGICAL           ,INTENT(IN)    :: LDONE(KFIELDS) 
LOGICAL           ,INTENT(OUT)   :: LDSAMPL(KFIELDS) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBLOCK

!     --------------------------------------------------------------------
!     IDER   : number of "boxes" for each field (for derivatives)
!     LLMONO : .TRUE. for monotonic interpolations
!     LLMNMX : .TRUE. to replace interpolations by the search of the 
!              nearest point

LOGICAL :: LLMONO(KFIELDS), LLMNMX(KFIELDS)
LOGICAL :: LLML

INTEGER(KIND=JPIM) :: IBOX(KGPEND), IDER(KFIELDS), IDMP(KFLDBUF)

INTEGER(KIND=JPIM) :: IADD, JF, JFLD, JI, IFPROW
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     --------------------------------------------------------------------

#include "abor1.intfb.h"
#include "fpint12.intfb.h"
#include "fpint4.intfb.h"

!     --------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPSAMPL',0,ZHOOK_HANDLE)
ASSOCIATE(LFPML_STD=>YDNAMFPINT%LFPML_STD, &
 & LWSTD04=>YDFPWSTD%LWSTD04, LWSTD12=>YDFPWSTD%LWSTD12, WSTD04=>YDFPWSTD%WSTD04, &
 & WSTD12=>YDFPWSTD%WSTD12, WSTDML=>YDFPWSTD%WSTDML, ML0=>YDFPWSTD%ML0)

!     --------------------------------------------------------------------

!*    1. BUFFER TRANSFER OR ACTUAL SAMPLING
!        ----------------------------------

IFPROW=SIZE(ML0,DIM=2)

IF (IFPROW == 0) THEN

  ! Buffer copy of all fields
  DO JF = 1, KFIELDS
    IF (.NOT.LDONE(JF)) THEN
      IADD=KASLB1*(JF-1)
      DO JI=KGPST,KGPEND
        PROW(JI,JF)=PBUF(ML0(JI,1,KBLOCK)+IADD)
      ENDDO
    ENDIF
  ENDDO

ELSE

  ! No interpolations => :  
  IDER(:)=1
  IBOX(:)=1
  IDMP(:)=0
  LLMONO(:)=.FALSE.
  ! Search for the nearest point in all cases : 
  LLMNMX(:)=.TRUE.
  LLML=LFPML_STD
  IF (LWSTD04) THEN
    IF (KBLOCK==1) WRITE (NULOUT,*) ' FPSAMPL: FPINT4 STD is called'
    CALL FPINT4(KASLB1,IFPROW,KFIELDS,KGPST,KGPEND,KFPROMA,KFLDBUF,IDER,IDMP,LLML,LLMNMX,&
     & IBOX,ML0(:,:,KBLOCK),WSTD04(:,1:4,KBLOCK),WSTDML(:,5:16,KBLOCK),LDONE,PBUF,PROW)
  ELSEIF (LWSTD12) THEN
    IF (KBLOCK==1) WRITE (NULOUT,*) ' FPSAMPL: FPINT12 STD is called'
    CALL FPINT12(KASLB1,IFPROW,KFIELDS,KGPST,KGPEND,KFPROMA,KFLDBUF,IDER,IDMP,LLML,LLMNMX,&
     & LLMONO,IBOX,ML0(:,:,KBLOCK),WSTD12(:,:,KBLOCK),WSTDML(:,5:16,KBLOCK),LDONE,PBUF,PROW)
  ELSE
    CALL ABOR1('FPSAMPL: NO INTERPOLATOR FOUND !')
  ENDIF

ENDIF

!*    2. CONTROL
!        -------

!     Signal sampling has been done
DO JFLD=1,KFIELDS
  LDSAMPL(JFLD)=.NOT.LDONE(JFLD)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPSAMPL',1,ZHOOK_HANDLE)
END SUBROUTINE FPSAMPL
