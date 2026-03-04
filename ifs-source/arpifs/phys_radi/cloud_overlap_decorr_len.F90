! (C) Copyright 2016- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE CLOUD_OVERLAP_DECORR_LEN &
     & (YDECLD,KIDIA, KFDIA, KLON, PGEMU, KDECOLAT, &
     &  PDECORR_LEN_EDGES_KM, PDECORR_LEN_WATER_KM, PDECORR_LEN_RATIO)

! CLOUD_OVERLAP_DECORR_LEN
!
! PURPOSE
! -------
!   Calculate the cloud overlap decorrelation length as a function of
!   latitude for use in the radiation scheme
!
! INTERFACE
! ---------
!   CLOUD_OVERLAP_DECORR_LEN is called from RADLSWR and RADIATION_SCHEME
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF (using code extracted from radlswr.F90)
!   Original: 2016-02-16
!
! MODIFICATIONS
! -------------
!
! -------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE YOECLD   , ONLY : TECLD

! -------------------------------------------------------------------

IMPLICIT NONE

! INPUT ARGUMENTS

! *** Array dimensions and ranges
TYPE(TECLD)       ,INTENT(INOUT):: YDECLD
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA    ! Start column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA    ! End column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KLON     ! Number of columns

! *** Configuration variable controlling the overlap scheme
INTEGER(KIND=JPIM),INTENT(IN) :: KDECOLAT

! *** Single-level variables 
REAL(KIND=JPRB),   INTENT(IN) :: PGEMU(KLON) ! Sine of latitude

! OUTPUT ARGUMENTS

! *** Decorrelation lengths for cloud edges and cloud water content,
! *** in km
REAL(KIND=JPRB), INTENT(OUT)           :: PDECORR_LEN_EDGES_KM(KLON)
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PDECORR_LEN_WATER_KM(KLON)
  
! Ratio of water-content to cloud-edge decorrelation lengths
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PDECORR_LEN_RATIO

! LOCAL VARIABLES

REAL(KIND=JPRB) :: ZRADIANS_TO_DEGREES, ZABS_LAT_DEG, ZCOS_LAT

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLOUD_OVERLAP_DECORR_LEN',0,ZHOOK_HANDLE)
  
! -------------------------------------------------------------------

IF (KDECOLAT == 0) THEN

  ! Decorrelation lengths are constant values
  PDECORR_LEN_EDGES_KM(KIDIA:KFDIA) = YDECLD%RDECORR_CF
  IF (PRESENT(PDECORR_LEN_WATER_KM)) THEN
    PDECORR_LEN_WATER_KM(KIDIA:KFDIA) = YDECLD%RDECORR_CW
  ENDIF
  IF (PRESENT(PDECORR_LEN_RATIO)) THEN
    PDECORR_LEN_RATIO = YDECLD%RDECORR_CW / YDECLD%RDECORR_CF
  ENDIF

ELSE

  ZRADIANS_TO_DEGREES = 180.0_JPRB / RPI

  IF (KDECOLAT == 1) THEN
    ! Shonk et al. (2010) Eq. 13 formula
    DO JL = KIDIA,KFDIA
      ZABS_LAT_DEG = ABS(ASIN(PGEMU(JL)) * ZRADIANS_TO_DEGREES)
      PDECORR_LEN_EDGES_KM(JL) = 2.899_JPRB - 0.02759_JPRB * ZABS_LAT_DEG
    ENDDO
  ELSE ! KDECOLAT == 2
    DO JL = KIDIA,KFDIA
      ! Shonk et al. (2010) but smoothed over the equator
      ZCOS_LAT = COS(ASIN(PGEMU(JL)))
      PDECORR_LEN_EDGES_KM(JL) = 0.75_JPRB + 2.149_JPRB * ZCOS_LAT*ZCOS_LAT
    ENDDO
  ENDIF

  ! Both KDECOLAT = 1 and 2 assume that the decorrelation length for
  ! cloud water content is half that for cloud edges
  IF (PRESENT(PDECORR_LEN_WATER_KM)) THEN
    PDECORR_LEN_WATER_KM(KIDIA:KFDIA) = PDECORR_LEN_EDGES_KM(KIDIA:KFDIA) * 0.5_JPRB
  ENDIF

  IF (PRESENT(PDECORR_LEN_RATIO)) THEN
    PDECORR_LEN_RATIO = 0.5_JPRB
  ENDIF

ENDIF

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLOUD_OVERLAP_DECORR_LEN',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUD_OVERLAP_DECORR_LEN
