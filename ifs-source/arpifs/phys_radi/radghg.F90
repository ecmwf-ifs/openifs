! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADGHG ( YDERAD,YDRIP,KIDIA, KFDIA, KLON, KLEV, &
  & KDLON, &
  & PPRESSURE_HL, PGEMU,&
  & PCO2 , PCH4 , PN2O , PNO2 , PCFC11 , PCFC12, PO3, PHCFC22, PCCL4 )

!***********************************************************************
! CAUTION: THIS ROUTINE WORKS ONLY ON A NON-ROTATED, UNSTRETCHED GRID
!***********************************************************************

!**** *RADGHG* - COMPUTES DISTRIBUTION OF GREENHOUSE GASES FROM CLIMATOLOGY

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *RADGHG* FROM *RADINTG*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!      J.-J. MORCRETTE  E.C.M.W.F.  20080424

!     MODIFICATIONS.
!     --------------
!      M.Hamrud    ECMWF 2003-10-01  CY28 Cleaning
!      JJMorcrette ECMWF 20090217    vertical interpolation from 91 lev climatol.
!      K. Yessad (July 2014): Move some variables.
!      R. Hogan  (Aug  2016): Correct result if pressure >= 1100 hPa
!      R. Hogan  (Jan  2019): Complete rewrite, and use YRADGHG for gas climatology
!-----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM ,   JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMRIP    , ONLY : TRIP
USE YOERADGHG , ONLY : YRADGHG
USE YOERAD    , ONLY : TERAD

!     -----------------------------------------------------------------

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPRESSURE_HL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: &
     & PCO2(KLON,KLEV), PCH4(KLON,KLEV), PN2O(KLON,KLEV), &
     & PNO2(KLON,KLEV), PCFC11(KLON,KLEV), PCFC12(KLON,KLEV), PO3(KLON,KLEV), &
     & PHCFC22(KLON,KLEV), PCCL4(KLON,KLEV)

!     ----------------------------------------------------------------- 

! Weight of ILAT+1 versus ILAT
REAL(KIND=JPRB) :: ZLATWEIGHT(KLON)

INTEGER(KIND=JPIM) :: ILAT(KLON), JI, JL

REAL(KIND=JPRB) :: ZINT, ZSINLAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RADGHG',0,ZHOOK_HANDLE)
ASSOCIATE(NGHGRAD=>YDERAD%NGHGRAD, &
 & NSTART=>YDRIP%NSTART, YDERADGHG=>YDRIP%YRERADGHG)
!     ------------------------------------------------------------------

!*         1.     LATITUDE INDEX WITHIN GAS CLIMATOLOGY
!                 ---------------------------------------

DO JI = KIDIA,KFDIA
  ZSINLAT = PGEMU(JI)
  IF (ZSINLAT > YRADGHG%SINLAT(YRADGHG%NLATITUDE)) THEN
    ! Near the north pole
    ILAT(JI) = YRADGHG%NLATITUDE-1
    ZLATWEIGHT(JI) = 1.0_JPRB ! Weight ILAT+1 entirely, i.e. the final point
  ELSEIF (ZSINLAT < YRADGHG%SINLAT(1)) THEN
    ! Near the south pole
    ILAT(JI) = 1
    ZLATWEIGHT(JI) = 0.0_JPRB ! Weight ILAT entirely, i.e. the first point 
  ELSE
    ILAT(JI) = 0
    DO JL = YRADGHG%NLATITUDE-1,1,-1
      IF (ZSINLAT <= YRADGHG%SINLAT(JL+1) .AND. ZSINLAT >= YRADGHG%SINLAT(JL)) THEN
        ILAT(JI)=JL
        EXIT
      ENDIF
    ENDDO
    IF (ILAT(JI) == 0) THEN
      CALL ABOR1('RADGHG: Error in latitude interpolation')
    ENDIF
    ZLATWEIGHT(JI) = (ZSINLAT-YRADGHG%SINLAT(ILAT(JI))) &
         &         / (YRADGHG%SINLAT(ILAT(JI)+1)-YRADGHG%SINLAT(ILAT(JI)))
  ENDIF
ENDDO

IF (NGHGRAD == 1 .OR. NGHGRAD >= 10) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_CO2,PPRESSURE_HL, PCO2)
ENDIF

IF (NGHGRAD == 2 .OR. NGHGRAD >= 11) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_CH4,PPRESSURE_HL, PCH4)
ENDIF

IF (NGHGRAD == 3 .OR. NGHGRAD >= 12) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_N2O,PPRESSURE_HL, PN2O)
ENDIF

IF (NGHGRAD == 4 .OR. NGHGRAD >= 20) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_NO2,PPRESSURE_HL, PNO2)
ENDIF

IF (NGHGRAD == 5 .OR. NGHGRAD >= 15) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_CFC11,PPRESSURE_HL, PCFC11)
ENDIF

IF (NGHGRAD == 6 .OR. NGHGRAD >= 16) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_CFC12,PPRESSURE_HL, PCFC12)
ENDIF

IF (NGHGRAD == 7 .OR. NGHGRAD == 17 .OR. NGHGRAD >= 21) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_O3,PPRESSURE_HL, PO3)
ENDIF

IF (NGHGRAD == 8 .OR. NGHGRAD >= 18) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_HCFC22,PPRESSURE_HL, PHCFC22)
ENDIF

IF (NGHGRAD == 9 .OR. NGHGRAD >= 19) THEN
  CALL INTERP_GAS(KIDIA,KFDIA,KLON,ILAT,ZLATWEIGHT,YRADGHG%NPRESSURE, KLEV, &
       &          YRADGHG%PRESSURE_HL,YRADGHG%MASS_CCL4,PPRESSURE_HL, PCCL4)
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADGHG',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE INTERP_GAS(KIDIA,KFDIA,KLON,KLAT,PLATWEIGHT,KLEVCLIM,KLEV, &
     &                PPRESS_CLIM_HL,PMASS,PPRESSURE_HL,PMMR)

  ! Index to start and end columns, and number of columns
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA, KFDIA, KLON
  ! Number of levels in climatology
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEVCLIM
  ! Number of levels in model
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
  ! Index to latitude dimension in climatology
  INTEGER(KIND=JPIM), INTENT(IN) :: KLAT(KLON)
  REAL(KIND=JPRB),    INTENT(IN) :: PLATWEIGHT(KLON)
  ! Half-level pressure of the climatology (Pa), indexed from 0 at TOA
  REAL(KIND=JPRB),    INTENT(IN) :: PPRESS_CLIM_HL(0:)
  REAL(KIND=JPRB),    INTENT(IN) :: PMASS(:,:)
  ! Half-level pressure of the model (Pa), indexed from 0 at TOA
  REAL(KIND=JPRB),    INTENT(IN) :: PPRESSURE_HL(KLON,0:KLEV)
  REAL(KIND=JPRB),    INTENT(OUT):: PMMR(KLON,KLEV)

  ! Mass and accumulated mass (both in Pa) in latitude-interpolated climatology
  REAL(KIND=JPRB) :: ZMASSCLIM(KLEVCLIM), ZMASSCLIMACC(0:KLEVCLIM)

  ! Accumulated mass for one model column (Pa)
  REAL(KIND=JPRB) :: ZMASSACC(0:KLEV)

  INTEGER(KIND=JPIM) :: JI,JLEV,ILEVCLIM,JLEVCLIM

  ! Loop over model column
  DO JI = KIDIA,KFDIA

    ! Interpolate in latitude
    ZMASSCLIM(1:KLEVCLIM) = PMASS(KLAT(JI),1:KLEVCLIM) &
      &  + PLATWEIGHT(JI)*(PMASS(KLAT(JI)+1,1:KLEVCLIM)-PMASS(KLAT(JI),1:KLEVCLIM))
 
    ZMASSCLIMACC(0) = 0.0_JPRB
    ! Accumulate climatology from top of atmosphere
    DO JLEV = 1,KLEVCLIM
      ZMASSCLIMACC(JLEV) = ZMASSCLIMACC(JLEV-1) + ZMASSCLIM(JLEV)
    ENDDO

    ! Interpolate accumulated mass in pressure space
    ZMASSACC(0) = 0.0_JPRB
    DO JLEV = 1,KLEV
      ILEVCLIM = KLEVCLIM-1

      DO JLEVCLIM = 0,KLEVCLIM-1
        IF (PPRESSURE_HL(JI,JLEV) >= PPRESS_CLIM_HL(JLEVCLIM) &
             &  .AND. PPRESSURE_HL(JI,JLEV) < PPRESS_CLIM_HL(JLEVCLIM+1)) THEN
          ILEVCLIM = JLEVCLIM
          EXIT
        ENDIF
      ENDDO
      ZINT = (PPRESSURE_HL(JI,JLEV)-PPRESS_CLIM_HL(ILEVCLIM)) &
           & / (PPRESS_CLIM_HL(ILEVCLIM+1)-PPRESS_CLIM_HL(ILEVCLIM))
      ZMASSACC(JLEV) = ZMASSCLIMACC(ILEVCLIM) &
           &         + ZINT * (ZMASSCLIMACC(ILEVCLIM+1) - ZMASSCLIMACC(ILEVCLIM))
    ENDDO

    ! Compute mass mixing ratio as difference between adjacent
    ! accumulated values
    DO JLEV = 1,KLEV
      PMMR(JI,JLEV) = (ZMASSACC(JLEV) - ZMASSACC(JLEV-1)) &
           &        / (PPRESSURE_HL(JI,JLEV) - PPRESSURE_HL(JI,JLEV-1))
    ENDDO

  ENDDO

END SUBROUTINE INTERP_GAS

END SUBROUTINE RADGHG
