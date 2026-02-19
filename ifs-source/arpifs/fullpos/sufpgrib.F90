! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPGRIB(CDFPCLIFNAME,YDFPGEO,KFPDOM,YDFPUSERGEO,YDFPOFN)

!**** *SUFPGRIB* - Routine to intitialize parameters for GRIB coding
!                - used for full-pos output file

!     Purpose.
!     --------
!         Set up parameters for GRIB coding including GRIB codes
!     for all parameters.

!**   Interface.
!     ----------
!        *CALL* *SUFPGRIB*

!        Explicit arguments :
!        --------------------

!          CDFPCLIFNAME : filenames of climatology file on target geometry
!          YDFPOFN : partial output filename (filename without time stamp)

!     Method.  
!     -------  

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-03-01

!     Modifications.
!     --------------
!      Modified : 01-05-31  D. Richardson - multi-analysis stream
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      S. Serrar    : 05-06-23 missing values to -9999.
!                     (in case section 3 defined for GRIB)
!      K. Yessad (Jan 2010): remove useless variables.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMFPGEO, ONLY : TFPGEO
USE TYPE_FPOFN, ONLY : TFPOFN
USE IOSTREAM_MIX , ONLY : INI_IOSTREAM
USE YOMMP0   , ONLY : NPROC
!     ------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*),   INTENT(IN)    :: CDFPCLIFNAME(:)
TYPE (TFPGEO),      INTENT(IN)    :: YDFPGEO
INTEGER(KIND=JPIM), INTENT(IN)    :: KFPDOM
TYPE (TFPUSERGEO) , INTENT(IN)    :: YDFPUSERGEO(KFPDOM)
TYPE (TFPOFN),      INTENT(INOUT) :: YDFPOFN(KFPDOM)

INTEGER(KIND=JPIM) :: ISTEP, ISTOP, IDIGITS, J, JROC, JI
REAL(KIND=JPRB) :: ZTSTEP
LOGICAL :: LLINC
INTEGER(KIND=JPIM) :: INUMPROCFP(YDFPGEO%NFPRGPG)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "suecfname.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUFPGRIB',0,ZHOOK_HANDLE)
ASSOCIATE(NFPRGPL=>YDFPGEO%NFPRGPL, NFPRGPLX=>YDFPGEO%NFPRGPLX, NFPRGPG=>YDFPGEO%NFPRGPG, &
 & NFPRGPNUM=>YDFPGEO%NFPRGPNUM, NFPRGPIND=>YDFPGEO%NFPRGPIND, &
 & NFPMAX=>YDFPUSERGEO%NFPMAX,CFPDOM=>YDFPUSERGEO%CFPDOM)
!     ------------------------------------------------------------------

!*       0.    SETUP OUTPUT FILE NAMES
!              -----------------------
ISTEP=5 ! dummy here
ZTSTEP=47._JPRB ! dummy here
ISTOP=6 ! dummy here
IDIGITS=6 ! is hard-coded in suecfname
LLINC=.FALSE. ! Dummy, anyway

DO J=1,KFPDOM

  CALL SUECFNAME(LLINC,'s',ISTEP,ZTSTEP,ISTOP,YDFPOFN(J)%CGG,YDFPOFN(J)%CSH)
  CALL SUECFNAME(LLINC,'m',ISTEP,ZTSTEP,ISTOP,YDFPOFN(J)%CUA,YDFPOFN(J)%CSH)
  YDFPOFN(J)%CSH=YDFPOFN(J)%CSH(1:LEN_TRIM(YDFPOFN(J)%CSH)-IDIGITS-1)//' '
  YDFPOFN(J)%CGG=YDFPOFN(J)%CGG(1:LEN_TRIM(YDFPOFN(J)%CGG)-IDIGITS-1)//' '
  YDFPOFN(J)%CUA=YDFPOFN(J)%CUA(1:LEN_TRIM(YDFPOFN(J)%CUA)-IDIGITS-1)//' '

  YDFPOFN(J)%CLI=CDFPCLIFNAME(J)

ENDDO

!     ------------------------------------------------------------------

!*       2.    SETUP IO STREAM
!              ---------------

CALL INI_IOSTREAM(KFPRGPL=NFPRGPL,KFPRGPG=NFPRGPG,KFPRGPLX=NFPRGPLX)
CALL INI_IOSTREAM(KFPMAX=NFPMAX)

! perhaps it is better to save numprocfp than nfprgpind ... REK
DO JROC=1,NPROC
  DO JI=1,NFPRGPNUM(JROC)
    INUMPROCFP(NFPRGPIND(JI,JROC))=JROC
  ENDDO
ENDDO
CALL INI_IOSTREAM(KNUMPROCFP=INUMPROCFP)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPGRIB',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPGRIB
