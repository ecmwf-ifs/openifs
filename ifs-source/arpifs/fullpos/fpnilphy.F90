! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPNILPHY(YDRQCLI,YDAFN,KFIELDS,KCOD,KST,KEND,KFPROMA,KFLDAUX,LDINTER,YDRQAUX,PAUX,PROW,LDNIL)

!**** *FPNILPHY*  - FULL-POS writer of the non-interpolatable fields

!     PURPOSE.
!     --------
!        To compute the fields that are not directly computable 
!        (ie : an auxilary computation has been necessary or a 
!        straight value is applied, unless abort is required).

!**   INTERFACE.
!     ----------
!       *CALL* *FPNILPHY*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KFIELDS: number of fields in row.
!         KST    : first output point in row.
!         KEND   : last output point in row.
!         KFPROMA: length of the output row.
!         KFLDAUX: total number of fields in auxilary fields row
!         LDINTER: .TRUE. if horizontal interpolation active
!         YDRQAUX : pointer of fields in auxilary row

!        INPUT/OUTPUT:
!         PAUX   : auxilary fields row
!         PROW   : output fields row

!        OUTPUT:
!         LDNIL : logical mask : .TRUE. if the non-interpolatable 
!                  field has been computed

!        IMPLICIT ARGUMENTS
!        --------------------
!          See module below

!     METHOD.
!     -------
!          Scan the request ; if it is a non-interpolatable field
!          then fill the row as possible or abort ; set the logical mask
!          to .TRUE. if the field row is filled.

!     EXTERNALS.
!     ----------
!      DECVLANISO

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-09-03 (partially from FPINTPHY)

!     MODIFICATIONS.
!     --------------
!      R. El Khatib 01-03-26 Modified test on KIND to prevent abort in Aladin E-zone
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad: 31-01-2007 Cleanings + optional calcul. of PROW,PAUX,LDNIL
!      R. El Khatib 13-Dec-2012 Fix erroneous tests
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMFP4L, ONLY : TRQFP, IFPSEARCH
USE YOMAFN   , ONLY : TAFN
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TRQFP),  INTENT(IN) :: YDRQCLI
TYPE (TAFN),       INTENT(IN)    :: YDAFN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOD(KFIELDS) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDAUX 
LOGICAL           ,INTENT(IN)    :: LDINTER 
TYPE (TRQFP),  INTENT(IN), OPTIONAL :: YDRQAUX
REAL(KIND=JPRB)   ,INTENT(IN)   , OPTIONAL :: PAUX(KFPROMA,KFLDAUX) 
REAL(KIND=JPRB)   ,INTENT(INOUT), OPTIONAL :: PROW(KFPROMA,KFIELDS) 
LOGICAL           ,INTENT(OUT)  , OPTIONAL :: LDNIL(KFIELDS) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDUMM(KFPROMA)
REAL(KIND=JPRB) :: ZU(KFPROMA), ZV(KFPROMA), ZSDO(KFPROMA)

INTEGER(KIND=JPIM) :: IERR, JFLD, ICOD, JI
INTEGER(KIND=JPIM) :: IASDO, IADOU, IADOV

LOGICAL :: LLALFA, LLGAMA, LLN, LLR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "decvlaniso.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPNILPHY',0,ZHOOK_HANDLE)
ASSOCIATE(GFP=>YDAFN%GFP, GFP_PHYDS=>YDAFN%GFP_PHYDS)

!     ------------------------------------------------------------------

!*       1. Initializations
!           ---------------

LLN=PRESENT(LDNIL)
LLR=(PRESENT(PROW).AND.PRESENT(PAUX))


!     ------------------------------------------------------------------

!*       2. LOOP ON PHYSICAL FIELDS
!           -----------------------

IERR=0

DO JFLD=1,KFIELDS

  ICOD=KCOD(JFLD)
  IF (IFPSEARCH(YDRQCLI,ICOD) <= 0) THEN

    IF (ICOD == GFP%ACOT%ICOD) THEN    
      ! Anisotropy of orography
      !  => Deconvolution of vector (sigma**2)*(1-gamma)/(1+gamma)
      IF (PRESENT(YDRQAUX)) THEN
        IASDO=IFPSEARCH(YDRQAUX,GFP%SDOG%ICOD)
        IADOU=IFPSEARCH(YDRQAUX,GFP%PADOU%ICOD)
        IADOV=IFPSEARCH(YDRQAUX,GFP%PADOV%ICOD)
      ELSE
        IASDO=0
        IADOU=0
        IADOV=0
      ENDIF
      IF (IASDO > 0.AND.IADOU > 0.AND.IADOV > 0) THEN
        IF (LLR) THEN
          LLGAMA=.TRUE.
          LLALFA=.FALSE.
          DO JI=KST,KEND
            ZU(JI)=PAUX(JI,IADOU)
            ZV(JI)=PAUX(JI,IADOV)
            ZSDO(JI)=PAUX(JI,IASDO)
          ENDDO
          CALL DECVLANISO(KST,KEND,KFPROMA,LLGAMA,LLALFA,PROW(1,JFLD), &
           & ZDUMM,ZU,ZV,ZSDO)  
        ENDIF
        IF (LLN) LDNIL(JFLD)=.TRUE.
      ELSE
        CALL ABOR1('FPNILPHY : INTERNAL ERROR IASDO/IADOU/IADOV')
      ENDIF
    ELSEIF(ICOD == GFP%DPAT%ICOD) THEN
      ! Main axis of orography
      !  => Deconvolution of vector (sigma**2)*(1-gamma)/(1+gamma)
      IF (PRESENT(YDRQAUX)) THEN
        IASDO=IFPSEARCH(YDRQAUX,GFP%SDOG%ICOD)
        IADOU=IFPSEARCH(YDRQAUX,GFP%PADOU%ICOD)
        IADOV=IFPSEARCH(YDRQAUX,GFP%PADOV%ICOD)
      ELSE
        IASDO=0
        IADOU=0
        IADOV=0
      ENDIF
      IF (IADOU > 0.AND.IADOV > 0) THEN
        IF (LLR) THEN
          LLGAMA=.FALSE.
          LLALFA=.TRUE.
          DO JI=KST,KEND
            ZU(JI)=PAUX(JI,IADOU)
            ZV(JI)=PAUX(JI,IADOV)
          ENDDO
          CALL DECVLANISO(KST,KEND,KFPROMA,LLGAMA,LLALFA,ZDUMM, &
           & PROW(1,JFLD),ZU,ZV,ZDUMM)  
        ENDIF
        IF (LLN) LDNIL(JFLD)=.TRUE.
      ELSE
        CALL ABOR1('FPNILPHY : INTERNAL ERROR IADOU/IADOV')
      ENDIF
    ELSEIF(LDINTER.AND.GFP_PHYDS(ICOD)%INTER > 0) THEN
      IF (ICOD == GFP%IVEG%ICOD  .OR. &
         & ICOD == GFP%RSMIN%ICOD .OR. &
         & ICOD == GFP%ARG%ICOD   .OR. &
         & ICOD == GFP%SAB%ICOD   .OR. &
         & ICOD == GFP%D2%ICOD    .OR. &
         & ICOD == GFP%LAI%ICOD   .OR. &
         & ICOD == GFP%BSR%ICOD   .OR. &
         & ICOD == GFP%BAAL%ICOD  .OR. &
         & ICOD == GFP%ALS%ICOD   .OR. &
         & ICOD == GFP%ALV%ICOD   .OR. &
         & ICOD == GFP%Z0H%ICOD) THEN  
        WRITE(NULOUT,'(''INTERPOLATION NOT ALLOWED ON FIELD '',A16,&
         & ''. TRY AGAIN WITH THE NEW OPTION NFPCLI=2'')') GFP_PHYDS(ICOD)%CLNAME  
        IERR=IERR+1
      ENDIF 
    ENDIF
  ENDIF

ENDDO

IF (IERR > 0) THEN
  CALL ABOR1('FPNILPHY : ABOR1 CALLED (IERR>0)')
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPNILPHY',1,ZHOOK_HANDLE)
END SUBROUTINE FPNILPHY
