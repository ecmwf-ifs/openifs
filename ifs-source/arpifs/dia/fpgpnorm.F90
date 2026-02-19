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

SUBROUTINE FPGPNORM(KFPXFLD,YDFPGEO,KFPDOM,CDFPDOM,KFPRESOL,KFPSIZEG,PGP,LDSURF,CDPREF,PLEV,KDOM,KDOMPTR,CDNAME,KGRIB,CDCONF)

!**** *FPGPNORM*  - Fullpos gridpoint norms

!     Purpose.
!     --------
!        Compute and print norms of gridpoint post-processed fields

!**   Interface.
!     ----------
!        *CALL* *FPGPNORM(CDCONF)

!        Explicit arguments :
!        --------------------
!         KFPXFLD : maximum number of fields to extract at a time
!         YDFPGEO: target mixed geometry
!         KFPDOM : number of subdomains
!         CDFPDOM: names of the subdomains
!         PGP    : distributed data array
!         LDSURF : .TRUE. if surface field (not 3D) for each field
!         CDPREF : prefix / kind of level for each field
!         PLEV   : level for each field
!         KDOM   : number of subdomains for each field
!         KDOMPTR: subdomains pointers for each field
!         CDNAME : FA name for each field
!         KGRIB  : GRIB code name for each field
!         CDCONF : configuration of I/Os

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    see includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Ryad El Khatib  *Meteo-France*
!      Original : 19-Nov-2015 from miscellaneous parts

! Modifications
! -------------
! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULOUT
USE YOMMP0    , ONLY : MYPROC, NPROC
USE YOMFPGEO  , ONLY : TFPGEO
USE MPL_MODULE, ONLY : MPL_BARRIER

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFPXFLD
TYPE (TFPGEO),  INTENT(IN) :: YDFPGEO
INTEGER(KIND=JPIM), INTENT(IN) :: KFPDOM
CHARACTER(LEN=*),   INTENT(IN) :: CDFPDOM(KFPDOM)
INTEGER(KIND=JPIM), INTENT(IN) :: KFPRESOL(KFPDOM)
INTEGER(KIND=JPIM), INTENT(IN) :: KFPSIZEG(KFPDOM)
REAL(KIND=JPRB),    INTENT(IN) :: PGP(:,:,:)
LOGICAL,            INTENT(IN) :: LDSURF(:)
CHARACTER(LEN=*),   INTENT(IN) :: CDPREF(:)
REAL(KIND=JPRB),    OPTIONAL, INTENT(IN) :: PLEV(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KDOM(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KDOMPTR(:,:)
CHARACTER(LEN=*),   INTENT(IN)  :: CDNAME(:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: KGRIB(:)
CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: CDCONF

INTEGER(KIND=JPIM), PARAMETER :: JPLEN=16
CHARACTER(LEN=JPLEN) :: CLEVEL, CLNOMA
INTEGER(KIND=JPIM) :: JFL, INCHIF, ILPRFU, ILSUFU, IOMASTER, ILEV, JDOM, JD, IFIELDS
INTEGER(KIND=JPIM) :: IFLDSCH, ICHUNKS, INFG, IBFA, IEFA, JCH
INTEGER(KIND=JPIM) :: INFD(NPROC), IFLDOFF(NPROC)
! ZAVEG,ZMING,ZMAXG  : global norms of the fields (mean val, min val, max val)
REAL(KIND=JPRB) :: ZAVEG(SIZE(KDOMPTR,DIM=1),SIZE(LDSURF))
REAL(KIND=JPRB) :: ZMING(SIZE(KDOMPTR,DIM=1),SIZE(LDSURF))
REAL(KIND=JPRB) :: ZMAXG(SIZE(KDOMPTR,DIM=1),SIZE(LDSURF))
REAL(KIND=JPRB) :: ZAVE(SIZE(LDSURF)),ZMIN(SIZE(LDSURF)),ZMAX(SIZE(LDSURF))
LOGICAL :: LLETRANS, LLAVE_ONLY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------

#include "gpnorm_trans.h"
#include "egpnorm_trans.h"

#include "abor1.intfb.h"
#include "updtrans.intfb.h"
#include "sumpioh.intfb.h"
#include "extfpnorm.intfb.h"

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPGPNORM',0,ZHOOK_HANDLE)

!*       0. CONTROL

IF (SIZE(PGP,DIM=2) /= SIZE(LDSURF)) CALL ABOR1('FPGPNORM: SIZES OF PGP(DIM=2) AND LDSURF MISMATCH')
IF (SIZE(CDPREF) /= SIZE(LDSURF)) CALL ABOR1('FPGPNORM: SIZES OF CDPREF AND LDSURF MISMATCH')
IF (SIZE(KDOM) /= SIZE(LDSURF)) CALL ABOR1('FPGPNORM: SIZES OF KDOM AND LDSURF MISMATCH')
IF (SIZE(KDOMPTR,DIM=2) /= SIZE(LDSURF)) CALL ABOR1('FPGPNORM: SIZES OF KDOMPTR(DIM=2) AND LDSURF MISMATCH')
IF (SIZE(CDNAME) /= SIZE(LDSURF)) CALL ABOR1('FPGPNORM: SIZES OF CDNAME AND LDSURF MISMATCH')
IF (PRESENT(KGRIB)) THEN
  IF (SIZE(KGRIB) /= SIZE(LDSURF)) CALL ABOR1('FPGPNORM: SIZES OF KGRIB AND LDSURF MISMATCH')
ENDIF
IF (.NOT.ALL(LDSURF(:))) THEN
  IF (PRESENT(PLEV)) THEN
    IF (SIZE(PLEV) /= SIZE(LDSURF)) CALL ABOR1('FPGPNORM: SIZES OF LDSURF AND PLEV MISMATCH')
  ELSE
    CALL ABOR1('FPGPNORM: PLEV IS NEEDED FOR UPPER AIR FIELDS')
  ENDIF
ENDIF

!*       1. SETUP

IOMASTER=1 ! notice : it is hard-coded to 1 inside (e)gpnorm_trans

IFIELDS=SIZE(ZAVEG,DIM=2)

IF (.FALSE.) THEN ! for tests
CALL UPDTRANS(KFPRESOL(1),LLETRANS)
LLAVE_ONLY=.FALSE.
IF (LLETRANS) THEN
  CALL EGPNORM_TRANS(PGP,IFIELDS,SIZE(PGP,DIM=1),ZAVE(:),ZMIN(:),ZMAX(:),LLAVE_ONLY,KRESOL=KFPRESOL(1))
ELSE
  CALL GPNORM_TRANS(PGP,IFIELDS,SIZE(PGP,DIM=1),ZAVE(:),ZMIN(:),ZMAX(:),LLAVE_ONLY,KRESOL=KFPRESOL(1))
ENDIF
DO JFL=1,IFIELDS
  ZAVEG(1,JFL)=ZAVE(JFL)
  ZMING(1,JFL)=ZMIN(JFL)
  ZMAXG(1,JFL)=ZMAX(JFL)
ENDDO
ENDIF

!     Maximum size of fields chunks :
IF (KFPXFLD <= 0) THEN
  IFLDSCH=IFIELDS
ELSE
  IFLDSCH=KFPXFLD
ENDIF
ICHUNKS=(IFIELDS-1)/IFLDSCH+1

DO JCH=1,ICHUNKS

!*       2.1 SETUP DISTRIBUTION AMONG PROCESSOR

  IBFA=(JCH-1)*IFLDSCH+1
  IEFA=MIN(JCH*IFLDSCH,IFIELDS)
  INFG=IEFA-IBFA+1
  CALL SUMPIOH(NPROC,NPROC,INFG,INFD,IFLDOFF)

!*       2.2 EXTRACT A CHUNK OF NORMS, ONLY FOR DOUBLE PREC.
  IF (JPRB == JPRD) THEN
    CALL EXTFPNORM(YDFPGEO,PGP,KDOM,KDOMPTR,KFPSIZEG,INFD,IFLDOFF,INFG,IBFA,IOMASTER, &
     & ZAVEG(:,IBFA:IEFA),ZMING(:,IBFA:IEFA),ZMAXG(:,IBFA:IEFA))
  ENDIF
ENDDO
CALL MPL_BARRIER (CDSTRING='FPGPNORMS:')

!*       3. PRINT OUT
!           only in double prec
IF (MYPROC == IOMASTER.AND.JPRB == JPRD)  THEN
! Print title the first time : 
  WRITE (NULOUT,FMT='(/,'' FULL-POS GPNORMS'')')
  DO JDOM=1,SIZE(KDOMPTR,DIM=1)
    IF (JDOM==1) THEN
      IF (.NOT. PRESENT(KGRIB)) THEN
        WRITE (NULOUT, FMT='(29X,3(6X,A7,9X))') 'AVERAGE','MINIMUM','MAXIMUM'
      ELSE
        WRITE (NULOUT, FMT='(26X,A7,3(9X,A7,6X))') 'PARAMID','AVERAGE','MINIMUM','MAXIMUM'
      ENDIF
    ENDIF
    DO JFL= 1, IFIELDS
      ILPRFU=LEN_TRIM(CDPREF(JFL))
      IF (LDSURF(JFL)) THEN
        ILEV=0
        INCHIF=0
      ELSEIF(PRESENT(CDCONF)) THEN
        CALL FA_FIELD_FORMAT(PLEV(JFL),CDCONF,INCHIF,ILEV)
      ELSE
        IF (NINT(PLEV(JFL)) /= 0) THEN
          ILEV=NINT(PLEV(JFL))
        ELSE
!         Most likely PV level => deci-PVU
          ILEV=NINT(PLEV(JFL)*10000000._JPRB)
        ENDIF
        INCHIF=MAX(INT(LOG10(ABS(PLEV(JFL))))+1,3)
      ENDIF
      DO JD=1,KDOM(JFL)
        IF (KDOMPTR(JD,JFL) == JDOM) THEN
!         A security while comparing listings of two experiments :
          IF (ABS(ZAVEG(JD,JFL)) < EPSILON(1._JPRB)) ZAVEG(JD,JFL) = 0._JPRB
          IF (ABS(ZMING(JD,JFL)) < EPSILON(1._JPRB)) ZMING(JD,JFL) = 0._JPRB
          IF (ABS(ZMAXG(JD,JFL)) < EPSILON(1._JPRB)) ZMAXG(JD,JFL) = 0._JPRB
          IF (INCHIF /= 0) THEN
            WRITE (UNIT=CLEVEL,FMT='(I8.8)') ILEV
            ILSUFU=MIN(JPLEN-ILPRFU-INCHIF,LEN_TRIM(CDNAME(JFL)))
            CLNOMA=CDPREF(JFL)(1:ILPRFU)//CLEVEL(9-INCHIF:8)//CDNAME(JFL)(1:ILSUFU)
            IF (.NOT.PRESENT(KGRIB)) THEN
              WRITE (NULOUT, FMT='(1X,A16,''/'',A7,'' : '',3(E21.15,1X))') &
               & CLNOMA, CDFPDOM(JDOM), ZAVEG(JD,JFL), ZMING(JD,JFL), ZMAXG(JD,JFL)
            ELSE
              CLNOMA=CDPREF(JFL)(1:ILPRFU)//CLEVEL(9-INCHIF:8)
              WRITE(NULOUT,FMT='(1X,A8,1X,A16,1X,I6.0,'' : '',3(E21.15,1X))') &
               & ADJUSTL(CLNOMA), CDNAME(JFL),KGRIB(JFL), &
               & ZAVEG(JD,JFL), ZMING(JD,JFL), ZMAXG(JD,JFL)
            ENDIF
          ELSE
            ILSUFU=MIN(JPLEN-ILPRFU-INCHIF,LEN_TRIM(CDNAME(JFL)))
            CLNOMA=CDPREF(JFL)(1:ILPRFU)//CDNAME(JFL)(1:ILSUFU)
            IF (.NOT.PRESENT(KGRIB)) THEN
              WRITE (NULOUT, FMT='(1X,A16,''/'',A7,'' : '',3(E21.15,1X))') &
             & CLNOMA, CDFPDOM(JDOM), ZAVEG(JD,JFL), ZMING(JD,JFL), ZMAXG(JD,JFL)
            ELSE
              WRITE (NULOUT, FMT='(1X,A8,1X,A16,1X,I6.0,'' : '',3(E21.15,1X))')  &
               & ADJUSTL(CLNOMA),CDNAME(JFL),KGRIB(JFL), &
               & ZAVEG(JD,JFL), ZMING(JD,JFL), ZMAXG(JD,JFL)
            ENDIF
          ENDIF
        ENDIF
      ENDDO ! JD
    ENDDO ! JFL
  ENDDO !JDOM
ENDIF

IF (LHOOK) CALL DR_HOOK('FPGPNORM',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE FA_FIELD_FORMAT(PLEV,CDCONF,KNCHIF,KLEV)

! Helper to define the format for printing the norms. Inspired from "fanfar".

REAL(KIND=JPRB),    INTENT(IN)  :: PLEV
CHARACTER(LEN=*),   INTENT(IN)  :: CDCONF
INTEGER(KIND=JPIM), INTENT(OUT) :: KNCHIF
INTEGER(KIND=JPIM), INTENT(OUT) :: KLEV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FPGPNORM:FA_FIELD_FORMAT',0,ZHOOK_HANDLE)

! PLEV   : raw level
! CDNAME : field name
! CDCONF : configuration of I/Os (kind of level)
 
! KNCHIF : number of digits for formatted level
! KLEV   : formatted level as integer

SELECT CASE (CDCONF)
CASE ('V')
!   Potential vorticity levels are written in deci-PVU :
  KLEV=NINT(PLEV*10000000._JPRB)
CASE ('K')
  KLEV=NINT(PLEV)
  IF (KLEV < 0) KLEV=-KLEV
CASE DEFAULT
  KLEV=NINT(PLEV)
END SELECT
IF (KLEV <= 0 ) THEN
  KNCHIF=0
ELSE
  SELECT CASE (CDCONF)
    CASE ('P','H')
      KNCHIF=5
    CASE ('F')
      KNCHIF=4
    CASE DEFAULT
      KNCHIF=3
  END SELECT
ENDIF

IF (LHOOK) CALL DR_HOOK('FPGPNORM:FA_FIELD_FORMAT',1,ZHOOK_HANDLE)

END SUBROUTINE FA_FIELD_FORMAT

END SUBROUTINE FPGPNORM
