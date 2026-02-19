! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE WROUTGPGB(PGGFLD,KGPTOT,KGGMAX,KGRIBIO,CDLEVTYPE,CDMODE,KRESOL,CDFN,PTSTEP)

!**** *WROUTGPGB* - GRIB codes and writes out grid-point fields

!     Purpose.
!     --------
!     Write out gridpoint fields in GRIB

!**   Interface.
!     ----------
!        *CALL* *WROUTGPGB*(...)

!        Explicit arguments : PGGFLD  - fields to be written out   
!        -------------------- KGPTOT  - number of gridpoints on your proc.
!                             KGGMAX  - number of fields to write out
!                             KGRIBIO - help array for GRIB coding
!                             CDLEVTYPE - type of level for fields
!                             CDMODE - mode to open files 
!                             KRESOL - resolution tag
!                             CDFN - output file name
!                             PTSTEP - model time step

!        Implicit arguments :      
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    - PRGRIBENC - modify GRIB headers, output resolution
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*
!       Original : 01-12-14 Adapted from WR****

!     Modifications.
!     --------------
!       J.Hague(IBM)/M.Hamrud 15-03-02  - Improvements in message passing
!       P.Towers 13-11-02 - Fixes for fewer writers than gatherers
!       John Hague : 13-01-03 WAIT if NPROC=1
!       P.Towers 23-04-03 - Added ISETFIELDCOUNTFDB logic
!       M.Hamrud      01-Oct-2003 CY28 Cleaning
!       M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!       M.Hamrud      01-Dec-2005 Generalized IO scheme
!       R. El Khatib : 01-Mar-2012 LFPOS => LECFPOS
!       T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!       K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      Ph. Lopez 24-Feb-2017 Added option to write out trajectory GP fields 
!                            in TL evolution experiments (conf 501 with LDTLEVOL=T). 
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMGRIB  , ONLY : NSTEPLPP  
USE YOMMP0   , ONLY : NOUTTYPE, NPROC
USE IOSTREAM_MIX , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,&
 & CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST, Y_IOSTREAM_FDB,&
 & CLOSE_IOREQUEST
USE MPL_MODULE , ONLY : MPL_ALLREDUCE

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KGPTOT 
INTEGER(KIND=JPIM), INTENT(IN) :: KGGMAX 
REAL(KIND=JPRB)   , INTENT(IN) :: PGGFLD(KGPTOT,KGGMAX) 
INTEGER(KIND=JPIM), INTENT(IN) :: KGRIBIO(2,KGGMAX) 
CHARACTER(LEN=1)  , INTENT(IN) :: CDLEVTYPE 
CHARACTER(LEN=1)  , INTENT(IN) :: CDMODE
INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
CHARACTER(LEN=*)  , INTENT(IN) :: CDFN
REAL(KIND=JPRB)   , INTENT(IN), OPTIONAL :: PTSTEP

!     ------------------------------------------------------------------

CHARACTER :: CLEVT*3

TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST
REAL(KIND=JPRB) :: ZGGFLD(KGPTOT,KGGMAX,1)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WROUTGPGB',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!  PREPARATIONS

IF(CDLEVTYPE == 'm') THEN
  CLEVT='ML'
ELSEIF(CDLEVTYPE == 'p') THEN
  CLEVT='PL'
ELSEIF(CDLEVTYPE == 'v') THEN
  CLEVT='PV'
ELSEIF(CDLEVTYPE == 't') THEN
  CLEVT='TH'
ELSEIF(CDLEVTYPE == 's') THEN
  CLEVT='SFC'
ELSEIF(CDLEVTYPE == 'o') THEN  !KPP
  CLEVT='OML'                  !KPP
ELSE
  CALL ABOR1('WROUTGPGB:UNKNOWN LEVEL TYPE')
ENDIF
CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,KRESOL=KRESOL,&
 & KGRIB2D=KGRIBIO(1,:),KLEVS2D=KGRIBIO(2,:),CDLEVTYPE=CLEVT,PTSTEP=PTSTEP)
ZGGFLD(:,:,1)=PGGFLD(:,:)
IF(NOUTTYPE == 2) THEN
  CALL IO_PUT(Y_IOSTREAM_FDB,YL_IOREQUEST,PR3=ZGGFLD)
ELSE
  CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CDFN),CDMODE=CDMODE,KIOMASTER=1)
  CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR3=ZGGFLD)
  CALL CLOSE_IOSTREAM(YL_IOSTREAM)
ENDIF
 
CALL CLOSE_IOREQUEST(YL_IOREQUEST)

!Make sure everyone knows when a surface parameter was last postprocessed
IF(NPROC > 1 .AND. CLEVT == 'SFC') THEN
  CALL MPL_ALLREDUCE(NSTEPLPP(:,2),'MAX',CDSTRING='WROUTGPGB:')
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WROUTGPGB',1,ZHOOK_HANDLE)
END SUBROUTINE WROUTGPGB
