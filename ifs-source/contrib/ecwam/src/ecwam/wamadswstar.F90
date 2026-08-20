! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WAMADSWSTAR(NXS, NXE, NYS, NYE, FIELDG, BLK2LOC, FF_NOW) 
! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!     RE-INITIALISES ADS, WSTAR, USTRA, VSTRA TO THE VALUES PROVIDED BY FIELDG

!**   INTERFACE.
!     ----------
!     *CALL* *WAMADSWSTAR*(NXS, NXE, NYS, NYE, FIELDG,
!                          BLK2LOC, ADS, WSTAR, USTRA, VSTRA)
!     *NXS:NXE*  FIRST DIMENSION OF FIELDG
!     *NYS:NYE*  SECOND DIMENSION OF FIELDG
!     *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!     *BLK2LOC*- POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *ADS*      AIR DENSITY IN KG/M3.
!     *WSTAR*    CONVECTIVE VELOCITy SCALE 
!     *USTRA*  - U-COMPONENT OF THE ATMOSPHERIC SURFACE STRESS
!     *VSTRA*  - V-COMPONENT OF THE ATMOSPHERIC SURFACE STRESS


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, FORCING_FIELDS

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW


      INTEGER(KIND=JWIM) :: ICHNK, IJ, IX, JY

      REAL(KIND=JPHOOK):: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMADSWSTAR',0,ZHOOK_HANDLE)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
      DO ICHNK = 1, NCHNK
        DO IJ = 1, NPROMA_WAM
          IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
          JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
          FF_NOW%AIRD(IJ,ICHNK) = FIELDG%AIRD(IX,JY)
          FF_NOW%WSTAR(IJ,ICHNK)= FIELDG%WSTAR(IX,JY)
          FF_NOW%USTRA(IJ,ICHNK)= FIELDG%USTRA(IX,JY)
          FF_NOW%VSTRA(IJ,ICHNK)= FIELDG%VSTRA(IX,JY)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

IF (LHOOK) CALL DR_HOOK('WAMADSWSTAR',1,ZHOOK_HANDLE)
 
END SUBROUTINE WAMADSWSTAR
