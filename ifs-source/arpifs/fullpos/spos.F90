! (C) Copyright 1989- Meteo-France.

SUBROUTINE SPOS(YDQTYPE,YDFPFILTERS,KRESOL,KSPEC2,PSPBFP,LDACTIVE)

!**** *SPOS*  - SPECTRAL POST-PROCESSING - FULL POS

!     PURPOSE.
!     --------
!        To modify the post-processed fields in the spectral space :
!        -Overtruncation of each field if required
!        -Horizontal filter for each derivative

!**   INTERFACE.
!     ----------
!       *CALL* *SPOS*

!        EXPLICIT ARGUMENTS
!        --------------------
!          YDFPFILTERS : filters
!          KRESOL : Resolution indicator
!          PSPBFP : spectral data arrray
!          LDACTIVE : .TRUE. if at least one field may have been filtered

!        IMPLICIT ARGUMENTS
!        --------------------
!          See #include below.

!     METHOD.
!     -------
!        See documentation about FULL-POS.

!     EXTERNALS.
!     ----------
!        Calls MXMAOP.
!        Is called by SPCM and ESPCM.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-04-10 Enable post-processing of filtered spectra + cleaning
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      R. El Khatib : 03-08-26 LFPMAT (bugfix)
!      M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!      F. Taillefer : 06-04-01 called by ALADIN too (instead of ESPOS)
!      K. Yessad    : Nov 2010 new treatment for stretched geometry.
!      R. El Khatib 08-Aug-2012 Move away the computation of psi,khi and the gaussian filters
!                               Replace model dimensionning by transforms variables
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV
USE YOMFPFILTERS, ONLY : TFPFILTERS
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE(TFPFILTERS), INTENT(IN)   :: YDFPFILTERS
INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
INTEGER(KIND=JPIM), INTENT(IN) :: KSPEC2
REAL(KIND=JPRB), INTENT(INOUT) :: PSPBFP(YDQTYPE%NFPSPB,KSPEC2)
LOGICAL, INTENT(OUT) :: LDACTIVE

!     II    : =1 if IM=0 ; =2 if not
!     ZSPD  : ALL input derivatives for the working wavenumber
!     ZSPD  : ALL output derivatives of ONE subdomain for the working wavenumber
!     ZSPB  : output derivatives for ONE field for the working wavenumber
!     IDIM  : first dimension of ZSPD, ZSPB
!     ISM   : start adress of matrix in RFPMAT

REAL(KIND=JPRB),ALLOCATABLE :: ZSPD(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSPB(:,:,:)

INTEGER(KIND=JPIM) :: IN, JDOM, JSE, J, IJ, JL, ISPB, IM, IDOM, INUMP
INTEGER(KIND=JPIM) :: IDIM, II, ISE, ISM, JN, ISPD, JD, JMLOC, ISMAX
INTEGER(KIND=JPIM) :: INVALUE(SIZE(PSPBFP,DIM=2))
INTEGER(KIND=JPIM), ALLOCATABLE :: IMYMS(:), IASM0(:)
INTEGER(KIND=JPIM) :: ISRANGE
LOGICAL :: LLETRANS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "trans_inq.h"
#include "etrans_inq.h"
#include "mxmaop.h"

#include "abor1.intfb.h"
#include "updtrans.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPOS',0,ZHOOK_HANDLE)
ASSOCIATE(LFPFIL=>YDFPFILTERS%LFPFIL ,RFPFIL=>YDFPFILTERS%RFPFIL ,LFPMAT=>YDFPFILTERS%LFPMAT ,RFPMAT=>YDFPFILTERS%RFPMAT, &
 & NFPSPD=>YDQTYPE%NFPSPD,NFPSPB=>YDQTYPE%NFPSPB)
!     ------------------------------------------------------------------

LDACTIVE=.FALSE.
CALL UPDTRANS(KRESOL,LLETRANS)

!*       1.    COMPUTATIONS FOR NON-DILATED OR NON-STRETCHED FIELDS
!              ----------------------------------------------------

IF (LLETRANS) THEN
  CALL ETRANS_INQ(KRESOL=KRESOL,KNVALUE=INVALUE(:))
ELSE
  CALL TRANS_INQ(KRESOL=KRESOL,KNVALUE=INVALUE(:))
ENDIF

ISPB=1
DO J=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  SELECT CASE (YDQTYPE%ISF(IJ))
  CASE (1:)
    DO JL=1,YDQTYPE%ILEV(IJ)
      IF (YDQTYPE%ISET(JL,IJ)==MYSETV) THEN
        DO JDOM=1,YDQTYPE%ISPD(JL,IJ) 
          SELECT CASE (YDQTYPE%ISF(IJ))
          CASE (3)
            IDOM=YDQTYPE%IDMP(JDOM,JL,IJ)
            IF ((.NOT.LFPMAT).AND.LFPFIL(IDOM)) THEN
              LDACTIVE=.TRUE.
              IF (LLETRANS) THEN
                PSPBFP(ISPB,:)=PSPBFP(ISPB,:)*RFPFIL(:,IDOM)
              ELSE
                DO JSE=1,KSPEC2
                  IN=INVALUE(JSE)
                  PSPBFP(ISPB,JSE)=PSPBFP(ISPB,JSE)*RFPFIL(IN,IDOM)
                ENDDO
              ENDIF
            ENDIF
          END SELECT
          ISPB=ISPB+1
        ENDDO
      ENDIF
    ENDDO
  END SELECT
ENDDO
IF (ISPB-1 /= NFPSPB) CALL ABOR1('SPOS:INTERNAL ERROR IN ISPB')

!     ------------------------------------------------------------------

!*       2.    COMPUTATIONS FOR DILATED STRETCHED FIELDS
!              -----------------------------------------

IF (LFPMAT .AND. ANY(LFPFIL)) THEN

  LDACTIVE=.TRUE.

  IF (LLETRANS) THEN
    CALL ABOR1('SPOS: WASNT IT A BUG USING NASM0 INSTEAD OF NESM0 ')
  ELSE
    CALL TRANS_INQ(KRESOL=KRESOL,KNUMP=INUMP,KSMAX=ISMAX)
    ALLOCATE(IMYMS(INUMP),IASM0(0:ISMAX))
    CALL TRANS_INQ(KRESOL=KRESOL,KMYMS=IMYMS(:),KASM0=IASM0(:))
  ENDIF

  ISM=1
  DO JMLOC=1,INUMP

    IM=IMYMS(JMLOC)
    ISRANGE=ISMAX+1-IM

!*            2.1 Allocations

    II=1+MIN(1,IABS(IM))
    IDIM=2*(ISRANGE/2)+1
    ALLOCATE(ZSPD(IDIM,NFPSPD,2),ZSPB(IDIM,NFPSPD,2))

!*            2.2 Input memory transfer

    ISPB=1
    ISPD=1
    DO J=1, YDQTYPE%NFPOSDYN
      IJ=YDQTYPE%NFPTRDYN(J)
      IF (YDQTYPE%ISF(IJ) >= 1) THEN
        DO JL=1,YDQTYPE%ILEV(IJ)
          IF (YDQTYPE%ISET(JL,IJ)==MYSETV) THEN
            IF (YDQTYPE%ISF(IJ)==3) THEN
              DO JN=IM,ISMAX
                ISE=IASM0(IM)+2*(JN-IM)
                ZSPD(JN+1-IM,ISPD,1)=PSPBFP(ISPB,ISE  )
                ZSPD(JN+1-IM,ISPD,2)=PSPBFP(ISPB,ISE+1)
              ENDDO
              ISPD=ISPD+1
            ENDIF
            ISPB=ISPB+YDQTYPE%ISPD(JL,IJ)
          ENDIF
        ENDDO 
      ENDIF
    ENDDO
    IF (ISPD-1 /= NFPSPD) CALL ABOR1('SPOS:INTERNAL ERROR IN TRANSFER - ISPD')
    IF (ISPB-1 /= NFPSPB) CALL ABOR1('SPOS:INTERNAL ERROR IN TRANSFER - ISPB')

!*             2.3 Filter then output memory transfer

!       Prevent unitialized values :
    ZSPB(:,:,2)=0.0_JPRB
    DO JDOM=1,SIZE(LFPFIL)
      IF (LFPFIL(JDOM)) THEN
        CALL MXMAOP(RFPMAT(ISM,JDOM),1,ISRANGE,ZSPD,1,IDIM,ZSPB,1,IDIM,&
         & ISRANGE,ISRANGE,II*NFPSPD)  
        ISPB=1
        ISPD=1
        DO J=1, YDQTYPE%NFPOSDYN
          IJ=YDQTYPE%NFPTRDYN(J)
          IF (YDQTYPE%ISF(IJ) >= 1) THEN
            DO JL=1,YDQTYPE%ILEV(IJ)
              IF (YDQTYPE%ISET(JL,IJ)==MYSETV) THEN
                IF (YDQTYPE%ISF(IJ)==3) THEN
                  DO JD=1,YDQTYPE%IDOM(JL,IJ)
                    IF (YDQTYPE%IDMP(JD,JL,IJ) == JDOM) THEN
                      DO JN=IM,ISMAX
                        ISE=IASM0(IM)+2*(JN-IM)
                        PSPBFP(ISPB+JD-1,ISE  )=ZSPB(JN+1-IM,ISPD,1)
                        PSPBFP(ISPB+JD-1,ISE+1)=ZSPB(JN+1-IM,ISPD,2)
                      ENDDO
                      EXIT
                    ENDIF
                  ENDDO
                  ISPD=ISPD+1
                ENDIF
                ISPB=ISPB+YDQTYPE%ISPD(JL,IJ)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IF (ISPD-1 /= NFPSPD) CALL ABOR1('SPOS:INTERNAL ERROR ON ISPD')
        IF (ISPB-1 /= NFPSPB) CALL ABOR1('SPOS:INTERNAL ERROR IN ISPB')
      ENDIF
    ENDDO

!*            2.4 Deallocations

    DEALLOCATE(ZSPD,ZSPB)

    ISM=ISM+ISRANGE**2

  ENDDO

  DEALLOCATE(IMYMS,IASM0)

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPOS',1,ZHOOK_HANDLE)
END SUBROUTINE SPOS
