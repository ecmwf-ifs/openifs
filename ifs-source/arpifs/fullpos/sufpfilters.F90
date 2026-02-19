! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPFILTERS(YDNAMFPSCI,YDNAMFPF,YDGEOMETRY,KLEVSNOW,KFPDOM,KFPCONF,KSF,CDNAME,YDFPFILTERS)

!**** *SUFPFILTERS*  - SET UP POST-PROCESSING SPECTRAL FILTER

!     PURPOSE.
!     --------
!        To initialize the profile of the spectral post-processing filter

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPFILTERS*

!        EXPLICIT ARGUMENTS
!        --------------------
!          YDGEOMETRY : model geometry
!          KFPDOM : number of subdomains
!          KSF    : kind of spectral filter for all fullpos fields
!          CDNAME : name of all fullpos fields
!          YDFPFILTERS : filters

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
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!     R. El Khatib : 01-04-09 : enable filtering for ECMWF and for spectra + cleanings
!     R. El Khatib : 01-08-07 Pruning options
!     R. El Khatib : 03-04-17 Fullpos improvments
!     R. El Khatib : 03-08-26 LFPMAT (bugfix)
!     M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!     R. El Khatib : 04-03-25 LFPWRFIL ,LFPRDFIL
!     F. Taillefer : 06-04-01 prepare ALADIN filters
!     R. El Khatib : 03-Feb-2009 Bugfix Filters in Arpege applied only if the
!                    target geometry is coarser than the model one.
!     K. Yessad    : Nov 2010 new treatment for stretched geometry.
!    R. El Khatib : 19-Jul-2012 NSPFIL*
!    R. El Khatib : 08-Aug-2012 Move the setup of gaussian filters to sposgf
!                               + replace model variables by transforms variables
!    R. El Khatib : 22-Aug-2012 LFPFIL switched off via the computation of NFMAX
!                   if spectral outputs because then the role of low-pass filter
!                   should be played by the target truncations. Anyway this
!                   filter can be re-activated by setting NFMAX in namelist
!    R. El Khatib : 12-Apr-2012 NFPREADALL
!    R. El Khatib : 27-May-2013 Compute matrixes here
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!    R. El Khatib : 09-Aug-2013 Protection against useless matrix usage
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : LELAM
USE YOMFPC   , ONLY : LTRACEFP, TNAMFPL, TNAMFPSCI
USE YOMFPF   , ONLY : TNAMFPF
USE YOMFPFILTERS, ONLY : TFPFILTERS
USE TYPE_FPFIELDS, ONLY : TFPFIELDS
USE YOMFP_SERV, ONLY : FP_SERV_C001
USE YOMVERT , ONLY : TVAB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TNAMFPSCI)   , INTENT(IN) :: YDNAMFPSCI
TYPE(TNAMFPF),      INTENT(IN) :: YDNAMFPF
TYPE(GEOMETRY),     INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN) :: KFPDOM
INTEGER(KIND=JPIM), INTENT(IN) :: KFPCONF
INTEGER(KIND=JPIM), INTENT(IN) :: KSF(:)
CHARACTER(LEN=*),   INTENT(IN) :: CDNAME(:)
TYPE(TFPFILTERS),   INTENT(OUT) :: YDFPFILTERS
INTEGER(KIND=JPIM), INTENT(IN) :: KLEVSNOW

LOGICAL :: LLFILTER_ACTIVE
INTEGER(KIND=JPIM) :: J, IFPCMAX, IREADALL, JFLD, ISIZE
INTEGER(KIND=JPIM) :: ICMAX(KFPDOM) ! what is actually got from matrixes files
INTEGER(KIND=JPIM) :: IFPMAX(KFPDOM) ! cutoff wavelength
TYPE(TFPFIELDS) :: YLFPFIELDS
TYPE(TNAMFPL) :: YLNAMFPL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "sufpc.intfb.h"
#include "sufpfields.intfb.h"
#include "cpfpfilter.intfb.h"
#include "wrfpfilter.intfb.h"
#include "rdfpfilter.intfb.h"
#include "fpfilter.intfb.h"
#include "efpfilter.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPFILTERS',0,ZHOOK_HANDLE)
ASSOCIATE(LFPBED=>YDNAMFPF%LFPBED, LFPWRFIL=>YDNAMFPF%LFPWRFIL, LFPRDFIL=>YDNAMFPF%LFPRDFIL, &
 & NFPGLFI=>YDNAMFPF%NFPGLFI, NFPCMAX=>YDNAMFPF%NFPCMAX, NFMAX=>YDNAMFPF%NFMAX, RFPMAXDEV=>YDNAMFPF%RFPMAXDEV, &
 & RFPLTF=>YDNAMFPF%RFPLTF, RFPBED=>YDNAMFPF%RFPBED, RFPSEL=>YDNAMFPF%RFPSEL, CFPDILA=>YDNAMFPF%CFPDILA, &
 & CFPCONT=>YDNAMFPF%CFPCONT, CFPMATRD=>YDNAMFPF%CFPMATRD, CFPMATWR=>YDNAMFPF%CFPMATWR, NFPREADALL=>YDNAMFPF%NFPREADALL, &
 & RSTRET=>YDGEOMETRY%YRGEM%RSTRET, NSMAX=>YDGEOMETRY%YRDIM%NSMAX, NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)


!*       1.    FIND OUT IF FILTERS ARE ACTIVE OR NOT
!              -------------------------------------

IF (.NOT.LELAM) THEN
  ! All this mess because their computation in stretched model is expensive :
  CALL SUFPC(YDNAMFPL=YLNAMFPL,LDPRINT=.FALSE.)
  CALL SUFPFIELDS(YLNAMFPL,NFLEVG,KLEVSNOW,YLFPFIELDS,LDPRINT=.FALSE.)
  LLFILTER_ACTIVE = &
   & (YDNAMFPSCI%NSPFILP == 3 .AND. YLFPFIELDS%NFP3P  > 0) .OR. &
   & (YDNAMFPSCI%NSPFILT == 3 .AND. YLFPFIELDS%NFP3TH > 0) .OR. &
   & (YDNAMFPSCI%NSPFILV == 3 .AND. YLFPFIELDS%NFP3PV > 0) .OR. &
   & (YDNAMFPSCI%NSPFILI == 3 .AND. YLFPFIELDS%NFP3I  > 0) .OR. &
   & (YDNAMFPSCI%NSPFILF == 3 .AND. YLFPFIELDS%NFP3F  > 0) .OR. &
   & (YDNAMFPSCI%NSPFILH == 3 .AND. YLFPFIELDS%NFP3H  > 0) .OR. &
   & (YDNAMFPSCI%NSPFILS == 3 .AND. YLFPFIELDS%NFP3S  > 0)
ENDIF

ISIZE=SIZE(KSF)

IF (.NOT.LLFILTER_ACTIVE) THEN
  DO JFLD=1,YLFPFIELDS%NFP3DF
    DO J=1, ISIZE
      IF (YLFPFIELDS%CFP3DF(JFLD) == CDNAME(J)(1:12)) EXIT
    ENDDO
    IF (J <= ISIZE) THEN
      IF (KSF(J) == 3) THEN
        LLFILTER_ACTIVE=.TRUE.
        EXIT
      ENDIF
    ENDIF
  ENDDO
  IF (.NOT.LLFILTER_ACTIVE) THEN
    DO JFLD=1,YLFPFIELDS%NFP2DF
      DO J=1,ISIZE
        IF (YLFPFIELDS%CFP2DF(JFLD) == CDNAME(J)(1:16)) EXIT
      ENDDO
      IF (J <= ISIZE) THEN
        IF (KSF(J) == 3) THEN
          LLFILTER_ACTIVE=.TRUE.
          EXIT
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF


!*       4.    COMPUTE FILTER
!              --------------
IF (KFPCONF==1) THEN
  IFPMAX(1:KFPDOM)=NFMAX(1:KFPDOM)
ELSE
  IFPMAX(:)=0
ENDIF
ALLOCATE(YDFPFILTERS%LFPFIL(KFPDOM))
IF (LELAM) THEN
  DO J=1,KFPDOM
  YDFPFILTERS%LFPFIL(J)=(IFPMAX(J) > 0)
  ENDDO
  YDFPFILTERS%LFPMAT=.FALSE.
ELSE
  IFPCMAX=MAX(NFPCMAX,NINT(REAL(NSMAX,KIND=JPRB)*RSTRET))
  DO J=1,KFPDOM
    YDFPFILTERS%LFPFIL(J)=(IFPMAX(J) > 0 .AND. IFPMAX(J) < IFPCMAX .AND. LLFILTER_ACTIVE)
  ENDDO
  YDFPFILTERS%LFPMAT=(IFPCMAX > NSMAX)
  IREADALL=MIN(1,MAX(0,NFPREADALL))
ENDIF

!*       5.    PRINT OUT FINAL VALUES
!              ----------------------

IF(LTRACEFP) THEN
  WRITE(UNIT=NULOUT,FMT='('' MODULE YOMFPF'')')
  WRITE(UNIT=NULOUT,FMT='('' LFPBED = '',L2,'' RFPBED = '',E14.7,'' RFPLTF = '',E14.7)') &
   & LFPBED, RFPBED, RFPLTF
  IF (.NOT.LELAM) THEN
    WRITE(UNIT=NULOUT,FMT='('' NFPCMAX = '',I5,'' RFPMAXDEV = '',E14.7)') IFPCMAX, RFPMAXDEV
    IF (ANY(YDFPFILTERS%LFPFIL)) THEN
      WRITE(UNIT=NULOUT,FMT='('' LFPWRFIL = '',L2,'' LFPRDFIL = '',L2, &
       & '' NFPGLFI = '',I3, '' NFPREADALL = '',I2)') LFPWRFIL, LFPRDFIL, NFPGLFI, IREADALL
    ENDIF
  ENDIF
  WRITE(UNIT=NULOUT,FMT='('' ( J   NFMAX          RFPSEL) '')')
  DO J=1,KFPDOM
  WRITE(UNIT=NULOUT,FMT='(10('' ('',I2,1X,I5,4X,E20.14'')''))') J,IFPMAX(J),RFPSEL(J)  
  ENDDO
  WRITE(UNIT=NULOUT,FMT='('' MODULE YOMFPFILTERS'')')
  WRITE(UNIT=NULOUT,FMT='('' LFPMAT = '',L2)') YDFPFILTERS%LFPMAT
  WRITE(UNIT=NULOUT,FMT='('' (J LFPFIL) '')')
  DO J=1,KFPDOM
    WRITE(UNIT=NULOUT,FMT='(15(''( '',I2,1X,L2,'' ) ''))') J,YDFPFILTERS%LFPFIL(J)
  ENDDO
ENDIF

IF (ANY(YDFPFILTERS%LFPFIL)) THEN
  IF (LELAM) THEN
    CALL EFPFILTER(YDGEOMETRY,KFPDOM,YDFPFILTERS%LFPFIL,IFPMAX,RFPLTF,YDFPFILTERS%RFPFIL)
  ELSE
    IF (YDFPFILTERS%LFPMAT) THEN
      ! stretched filters
      IF (LFPRDFIL) THEN
        CALL RDFPFILTER(YDGEOMETRY,IREADALL,KFPDOM,IFPCMAX,IFPMAX,NFPGLFI,LFPBED,RFPBED, &
         & RFPLTF,CFPMATRD,ICMAX,YDFPFILTERS%LFPFIL,YDFPFILTERS%RFPMAT)
        IF (LFPWRFIL) THEN
          WRITE(NULOUT,'('' FILTERING MATRIXES WILL NOT BE WRITTEN AFTER BEING READ'')')
        ENDIF
      ELSE
        IF ((.NOT. FP_SERV_C001%LFP_SERVER) .OR. .NOT.FP_SERV_C001%LFP_SERVER_FPMTS) THEN
          CALL CPFPFILTER(YDGEOMETRY,CFPDILA,CFPCONT,KFPDOM,IFPCMAX,RFPSEL,IFPMAX, &
           & RFPBED,RFPLTF,LFPBED,ICMAX,RSTRET,RFPMAXDEV,YDFPFILTERS%LFPFIL,YDFPFILTERS%RFPFIL,YDFPFILTERS%RFPMAT)
        ENDIF
        IF (LFPWRFIL) THEN
          CALL WRFPFILTER(YDGEOMETRY,KFPDOM,SIZE(YDFPFILTERS%RFPMAT,DIM=1),ICMAX,IFPMAX,NFPGLFI,LFPBED,RFPBED,RFPLTF,CFPMATWR,&
           & YDFPFILTERS%LFPFIL,YDFPFILTERS%RFPMAT)
        ENDIF
      ENDIF
    ELSE
      CALL FPFILTER(NSMAX,KFPDOM,YDFPFILTERS%LFPFIL,IFPMAX,RFPSEL,RFPBED,RFPLTF,LFPBED,YDFPFILTERS%RFPFIL)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPFILTERS',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPFILTERS
