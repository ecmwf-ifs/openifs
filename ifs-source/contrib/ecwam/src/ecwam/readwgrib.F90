! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE READWGRIB(IU06, FILNM, IPARAM, CDATE,        &
     &                     BLK2LOC,                           &
     &                     NXS, NXE, NYS, NYE, FIELDG,        &
     &                     FIELD, KZLEV, LLONLYPOS, IREAD )

!-----------------------------------------------------------------------

!**** *READWGRIB*  READS FROM GRIB WAVE MODEL FIELD

!     J. BIDLOT   ECMWF   OCTOBER 1997 

!*    PURPOSE.
!     --------

!       SPECIAL INPUT FROM GRIB WAVE FIELD (SEE METHOD BELOW !!!) 

!**   INTERFACE.
!     ----------

!       *CALL* *READWGRIB*(IU06, FILNM, IPARAM, CDATE,
!    &                     IFROMIJ, JFROMIJ,
!    &                     NXS, NXE, NYS, NYE, FIELDG,
!    &                     FIELD, KZLEV, LLONLYPOS, IREAD )

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *FILNM*     DATA INPUT FILENAME.
!      *IPARAM*    INTEGER   PARAMETER IDENTIFIER OF FIELD
!      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD 
!      *IFROMIJ*   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!      *JFROMIJ*   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!      *NXS:NXE*   FIRST DIMENSION OF FIELDG
!      *NYS:NYE*   SECOND DIMENSION OF FIELDG
!      *FIELDG*    INPUT FORCING FIELDS ON THE WAVE MODEL GRID

!      *FIELD*     REAL      WAVE FIELD IN BLOCK FORMAT 
!      *KZLEV*     INTEGER   REFERENCE LEVEL IN full METER
!                           (SHOULD BE 0 EXCEPT FOR 233, 245 AND 249 WHERE IT
!                           MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
!                           CODED, SEE *GRIBPAC*
!      *LLONLYPOS* LOGICAL   TAKE ONLY THE POSITIVE VALUES
!      *IREAD*     INTEGER  PROCESSOR WHICH WILL ACCESS THE FILE ON DISK


!     METHOD.
!     -------
!      READS GRIB DATA USING GRIB_API, CHECK WHETHER GRID DEFINITION ARE 
!      COMPATIBLE. IN CASE THE INPUT DATA IS DEFINED ON A REGULAR LAT-
!      LON GRID, BUT THE MODEL USES AN IRREGULAR GRID, THE DATA WILL BE 
!      INTERPOLATED TO THE IRREGULAR GRID BY TAKING THE NEAREST GRID
!      POINTS.
!      THE GRID DATA ARE PUT INTO BLOCK FORMAT PROVIDED THEY ARE
!      DIFFERENT THAN ZMISS TO AVOID THE USE OF MISSING DATA INDICATOR
!      FOR LAND POINT FOLLOWING A POSSIBLE CHANGE IN THE LAND-SEA
!      MASK . FOR THAT REASON, 
!      THE BLOCK VALUE OF FIELD SHOUD BE INITIALISED PRIOR TO THE CALL
!      TO THIS ROUTINE.

!     EXTERNALS.
!     ----------

!      *ABORT1*
!      *INWGRIB*

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS, WVGRIDLOC

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWMAP   , ONLY : NLONRGG  , NGX      ,NGY      ,NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPCONS , ONLY : ZMISS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!-----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "inwgrib.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06
      CHARACTER(LEN=24), INTENT(IN) :: FILNM
      INTEGER(KIND=JWIM), INTENT(IN) :: IPARAM
      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: FIELD 
      INTEGER(KIND=JWIM), INTENT(INOUT) :: KZLEV
      LOGICAL, INTENT(IN) :: LLONLYPOS
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD


      INTEGER(KIND=JWIM) :: KPARAM
      INTEGER(KIND=JWIM) :: ICHNK, IJ, IX, JY
      INTEGER(KIND=JWIM) :: KLONRGG(NGY)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: WORK(NXS:NXE, NYS:NYE)

      CHARACTER(LEN=14) :: CCDDATE
      CHARACTER(LEN=40) :: MSG

!-----------------------------------------------------------------------

!*    1. INPUT OF GRIB DATA.
!     -----------------------

      IF (LHOOK) CALL DR_HOOK('READWGRIB',0,ZHOOK_HANDLE)

      CALL INWGRIB(FILNM, IREAD, CCDDATE, KPARAM, KZLEV,  &
     &             NXS, NXE, NYS, NYE, FIELDG, WORK)

!*    SIMPLE CHECKS ON THE RETRIEVED DATA 
!     -----------------------------------

      IF (KPARAM /= IPARAM) THEN
        WRITE(IU06,*)'********************************'
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'* FATAL ERROR IN SUB READWGRIB *'
        WRITE(IU06,*)'* ===========================  *'
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'* GRIB PARAMETER  ',KPARAM
        WRITE(IU06,*)'* WAS READ INSTEAD OF ',IPARAM
        WRITE(IU06,*)'* IN FILE: ',FILNM 
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'********************************'
        CALL ABORT1
      ENDIF
      IF (CCDDATE /= CDATE) THEN
        WRITE(IU06,*)'**********************************'
        WRITE(IU06,*)'*                                *'
        WRITE(IU06,*)'* FATAL ERROR IN SUB READWGRIB   *'
        WRITE(IU06,*)'* ===========================    *'
        WRITE(IU06,*)'*                                *'
        WRITE(IU06,*)'* REQUESTED DATE IS NOT EQUAL TO *' 
        WRITE(IU06,*)'* RETRIEVED DATE.                *' 
        WRITE(IU06,*)'* IN FILE: ',FILNM 
        WRITE(IU06,*)'* CCDDATE = ',CCDDATE
        WRITE(IU06,*)'* CDATE = ',CDATE
        WRITE(IU06,*)'*                                *'
        WRITE(IU06,*)'**********************************'
        CALL ABORT1
      ENDIF


! TRANSFORM GRID DATA TO BLOCK DATA

      IF (LLONLYPOS) THEN
        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
        DO ICHNK = 1, NCHNK
          DO IJ = 1, NPROMA_WAM
            IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
            JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
            IF (WORK(IX,JY) /= ZMISS .AND. WORK(IX,JY) > 0.0_JWRB) FIELD(IJ, ICHNK) = WORK(IX,JY)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

      ELSE
        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
        DO ICHNK = 1, NCHNK
          DO IJ = 1, NPROMA_WAM
            IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
            JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
            IF (WORK(IX,JY) /= ZMISS) FIELD(IJ, ICHNK) = WORK(IX,JY)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)
      ENDIF

      IF (LHOOK) CALL DR_HOOK('READWGRIB',1,ZHOOK_HANDLE)

      END SUBROUTINE READWGRIB
