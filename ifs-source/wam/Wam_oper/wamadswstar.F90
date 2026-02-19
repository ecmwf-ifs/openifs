! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WAMADSWSTAR(NXS, NXE, NYS, NYE, FIELDG,            &
 &                     BLK2LOC, ADS, WSTAR)
! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!     RE-INITIALISES ADS AND WSTAR TO THE VALUES PROVIDED BY FIELDG

!**   INTERFACE.
!     ----------
!     *CALL* *WAMADSWSTAR*(NXS, NXE, NYS, NYE, FIELDG,
!                          BLK2LOC, ADS, WSTAR)
!     *NXS:NXE*  FIRST DIMENSION OF FIELDG
!     *NYS:NYE*  SECOND DIMENSION OF FIELDG
!     *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!     *BLK2LOC*- POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *ADS*      AIR DENSITY IN KG/M3.
!     *WSTAR*    CONVECTIVE VELOCITy SCALE 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, FORCING_FIELDS

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), DIMENSION(NXS:NXE, NYS:NYE), INTENT(IN) :: FIELDG
      TYPE(WVGRIDLOC), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: BLK2LOC
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(OUT):: ADS, WSTAR


      INTEGER(KIND=JWIM) :: ICHNK, IJ, IX, JY

      REAL(KIND=JPHOOK):: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMADSWSTAR',0,ZHOOK_HANDLE)

ASSOCIATE(IFROMIJ => BLK2LOC%IFROMIJ, &
 &        JFROMIJ => BLK2LOC%JFROMIJ)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
      DO ICHNK = 1, NCHNK
        DO IJ = 1, NPROMA_WAM
          IX = IFROMIJ(IJ, ICHNK)
          JY = JFROMIJ(IJ, ICHNK)
          ADS(IJ, ICHNK) = FIELDG(IX,JY)%AIRD
          WSTAR(IJ, ICHNK)= FIELDG(IX,JY)%WSTAR
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('WAMADSWSTAR',1,ZHOOK_HANDLE)
 
END SUBROUTINE WAMADSWSTAR
