! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPF(YDNAMFPF,KFPCONF,YDFPUSERGEO)

!**** *SUFPF*  - SET UP POST-PROCESSING SPECTRAL FILTER

!     PURPOSE.
!     --------
!        To initialize the profile of the spectral post-processing filter

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPF*

!        EXPLICIT ARGUMENTS
!        --------------------
!           KFPCONF : configuration of the post-processing :
!                     0 : vertical interpolation only (<CFPFMT='MODEL'>)
!                     1 : gridpoint post-processing, possibly with spectral filters (<NFPOS=1>)
!                     2 : gridpoint/spectral post-processing (spectral outputs possible) (<NFPOS=2>)

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
!    R. El Khatib : 27-02-2019 Support for simple precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULNAM, NULOUT
USE YOMCT0   , ONLY : LELAM
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMFPF   , ONLY : TNAMFPF
USE YOMCST   , ONLY : RPI, RA 

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TNAMFPF), TARGET, INTENT(OUT) :: YDNAMFPF
INTEGER(KIND=JPIM), INTENT(IN) :: KFPCONF
TYPE (TFPUSERGEO),  INTENT(IN) :: YDFPUSERGEO(:)

LOGICAL, POINTER :: LFPBED
LOGICAL, POINTER :: LFPWRFIL
LOGICAL, POINTER :: LFPRDFIL
INTEGER(KIND=JPIM), POINTER :: NFPREADALL
INTEGER(KIND=JPIM), POINTER :: NFPGLFI
INTEGER(KIND=JPIM), POINTER :: NFPCMAX
INTEGER(KIND=JPIM), POINTER :: NFMAX(:)
REAL(KIND=JPRB), POINTER :: RFPMAXDEV
REAL(KIND=JPRB), POINTER :: RFPLTF
REAL(KIND=JPRB), POINTER :: RFPBED
REAL(KIND=JPRB), POINTER :: RFPSEL(:)
CHARACTER(LEN=256), POINTER :: CFPDILA
CHARACTER(LEN=256), POINTER :: CFPCONT
CHARACTER(LEN=256), POINTER :: CFPMATRD(:)
CHARACTER(LEN=256), POINTER :: CFPMATWR(:)

INTEGER(KIND=JPIM) :: J
REAL(KIND=JPRB) :: Z1, Z2, ZEPS, ZNFMAX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "posnam.intfb.h"

#include "namfpf.nam.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUFPF',0,ZHOOK_HANDLE)
ASSOCIATE(NFPMAX=>YDFPUSERGEO%NFPMAX, FPSTRET=>YDFPUSERGEO%FPSTRET, CFPGRID=>YDFPUSERGEO%CFPGRID, &
 & NLAT=>YDFPUSERGEO%NLAT, NLON=>YDFPUSERGEO%NLON, NFPLUX=>YDFPUSERGEO%NFPLUX, NFPGUX=>YDFPUSERGEO%NFPGUX, &
 & RDELX=>YDFPUSERGEO%RDELX, RDELY=>YDFPUSERGEO%RDELY, CFPDOM=>YDFPUSERGEO%CFPDOM)
!     ------------------------------------------------------------------

WRITE(NULOUT,'('' == Full-Pos : setup spectral filters == '')')

!*       0. POINTERS TO THE NAMELIST STRUCTURE
!           ----------------------------------

LFPBED=>YDNAMFPF%LFPBED
LFPWRFIL=>YDNAMFPF%LFPWRFIL
LFPRDFIL=>YDNAMFPF%LFPRDFIL
NFPREADALL=>YDNAMFPF%NFPREADALL
NFPGLFI=>YDNAMFPF%NFPGLFI
NFPCMAX=>YDNAMFPF%NFPCMAX
NFMAX=>YDNAMFPF%NFMAX(:)
RFPMAXDEV=>YDNAMFPF%RFPMAXDEV
RFPLTF=>YDNAMFPF%RFPLTF
RFPBED=>YDNAMFPF%RFPBED
RFPSEL=>YDNAMFPF%RFPSEL(:)
CFPDILA=>YDNAMFPF%CFPDILA
CFPCONT=>YDNAMFPF%CFPCONT
CFPMATRD=>YDNAMFPF%CFPMATRD(:)
CFPMATWR=>YDNAMFPF%CFPMATWR(:)

!*       1.    SET DEFAULT VALUES
!              ------------------

LFPBED=.TRUE.

NFPCMAX=0
IF (JPRB == JPRD) THEN
  ! Double precision version
  RFPMAXDEV=1.E-10_JPRB
ELSE
  RFPMAXDEV=0.2E-5_JPRB
ENDIF
NFPREADALL=1
LFPWRFIL=.FALSE.
LFPRDFIL=.FALSE.
CFPDILA='MATDILA'
CFPCONT='MATCONT'
NFPGLFI=12

IF (LELAM) THEN
  RFPBED=0.0_JPRB
  RFPLTF=6._JPRB
ELSE
! Default value of RFPBED for continuity with the previous formulation ;
! gives roughly RFPBED = 3.08
  RFPBED=-LOG(LOG(9._JPRB)/48._JPRB)
  RFPLTF=4._JPRB
ENDIF


IF (KFPCONF==1) THEN
  ZEPS=1.E-10_JPRB
  DO J=1, SIZE(YDFPUSERGEO)
    IF (CFPGRID(J) == 'GAUSS') THEN
      NFMAX(J)=NINT(REAL(NFPMAX(J),JPRB)*FPSTRET(J))
    ELSEIF (CFPGRID(J) == 'LALON') THEN
      ! This default value considers the local equivalent of a quadratic gaussian grid
      Z1=(360._JPRB*REAL(NLON(J),JPRB)/ &
       & MAX(RDELX(J)*REAL(NLON(J)-1,JPRB),ZEPS)-1.0_JPRB)/3._JPRB
      Z2=(360._JPRB*REAL(NLAT(J),JPRB)/ &
       & MAX(RDELY(J)*REAL(NLAT(J)-1,JPRB),ZEPS)-1.0_JPRB)/3._JPRB
      ZNFMAX=MIN(Z1,Z2)
      NFMAX(J)=INT(ZNFMAX)
    ELSE
      ! This default value considers the local equivalent of a quadratic gaussian grid
      Z1=(2._JPRB*RPI*RA*REAL(NFPLUX(J),JPRB)/ &
       & MAX(RDELX(J)*REAL(NFPLUX(J)-1,JPRB),ZEPS)-1.0_JPRB)/3._JPRB
      Z2=(2._JPRB*RPI*RA*REAL(NFPGUX(J),JPRB)/ &
       & MAX(RDELY(J)*REAL(NFPGUX(J)-1,JPRB),ZEPS)-1.0_JPRB)/3._JPRB
      ZNFMAX=MIN(Z1,Z2)
      NFMAX(J)=INT(ZNFMAX)
    ENDIF
    RFPSEL(J)=1.0_JPRB
    CFPMATRD(J)='matrix.fil.'//CFPDOM(J)
    CFPMATWR(J)='matrix.fil.'//CFPDOM(J)
  ENDDO
ELSE
  NFMAX(:)=0
  RFPSEL(:)=1.0_JPRB
  CFPMATRD(:)=' '
  CFPMATWR(:)=' '
ENDIF

!*       2.    READ NAMELIST
!              -------------

CALL POSNAM(NULNAM,'NAMFPF')
READ(NULNAM,NAMFPF)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPF',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPF
