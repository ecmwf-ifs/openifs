! (C) Copyright 2019- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOE_SPECTRAL_PLANCK

! YOE_SPECTRAL_PLANCK
!
! PURPOSE
! -------
!   Calculate Planck function integrated across user-specified
!   spectral intervals, used in RADHEATN by approximate longwave
!   update scheme to modify longwave fluxes to account for the
!   spectral emissivity on the high-resolution model grid (rather than
!   the lower resolution grid seen by the radiation scheme).
!
! INTERFACE
! ---------
!   Call the INIT member routine to configure the look-up table of the
!   TSPECRALPLANCK type, followed by any number of CALC calls with the
!   temperatures at which the Planck function is required. FREE then
!   deallocates memory.
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF
!   Original: 2019-02-04
!
! MODIFICATIONS
! -------------
!   A Dawson 2019-08-05 avoid single precision overflow in INIT

!-----------------------------------------------------------------------

USE PARKIND1, ONLY : JPRB,JPRD,JPIM
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
! Type for storing Planck function look-up table
TYPE TSPECTRALPLANCK
  ! Number of intervals over which the integrated Planck function is
  ! required. Note that an interval need not be contiguous in
  ! wavelength.
  INTEGER(KIND=JPIM) :: NINTERVALS

  ! Number of temperatures in look-up table
  INTEGER(KIND=JPIM) :: NTEMPS

  ! Start temperature and temperature spacing of look-up table
  REAL(KIND=JPRB) :: TEMP1, DTEMP

  ! Integrated Planck functions in look-up table, dimensioned
  ! (NINTERVALS,NTEMPS)
  REAL(KIND=JPRB),    ALLOCATABLE :: PLANCK_LUT(:,:)

  ! Store interval data
  REAL(KIND=JPRB),    ALLOCATABLE :: WAVLEN_BOUND(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: INTERVAL_MAP(:)

CONTAINS
  PROCEDURE :: INIT
  PROCEDURE :: CALC
  PROCEDURE :: PRINT=> PRINT_SPECTRAL_PLANCK
  PROCEDURE :: FREE => FREE_SPECTRAL_PLANCK

END TYPE TSPECTRALPLANCK

CONTAINS

!-----------------------------------------------------------------------
! Generate a Planck function look-up table consisting of KINTERVALS
! spectral intervals (which need not be contiguous in wavelength),
! whose wavelength bounds are defined by PWAVLEN_BOUND and mapping
! on to KINTERVALS described by KINTERVAL_MAP.
SUBROUTINE INIT(SELF, KINTERVALS, PWAVLEN_BOUND, KINTERVAL_MAP)

  USE YOMCST,   ONLY : RPI, RKBOL, RHPLA, RCLUM
  USE YOMLUN,   ONLY : NULOUT

  CLASS(TSPECTRALPLANCK), INTENT(INOUT) :: SELF
  INTEGER(KIND=JPIM)    , INTENT(IN)    :: KINTERVALS
  REAL(KIND=JPRB)       , INTENT(IN)    :: PWAVLEN_BOUND(:)
  INTEGER(KIND=JPIM)    , INTENT(IN)    :: KINTERVAL_MAP(:)

  ! Current temperature (K)
  REAL(KIND=JPRB) :: ZTEMP

  ! Combinations of constants in the Planck function
  REAL(KIND=JPRB) :: ZCOEFF1, ZCOEFF2

  ! Wavelengths at start and end of range
  REAL(KIND=JPRB) :: ZWAVLEN1, ZWAVLEN2, DWAVLEN

  ! Wavelength, wavelength squared
  REAL(KIND=JPRB) :: ZWAVLEN, ZWAVLEN_SQR

  ! Sum of Planck values, integration weight
  REAL(KIND=JPRB) :: ZSUM, ZWEIGHT

  ! A double-precision temporary to hold the exponential term in Planck's law.
  ! A single precision float can overflow during this calculation.
  REAL(KIND=JPRD) :: ZPLANCKEXP

  ! Number of wavelength ranges represented by PWAVLEN_BOUND and
  ! KINTERVAL_MAP
  INTEGER(KIND=JPIM) :: NRANGES, NWAVLEN

  INTEGER(KIND=JPIM) :: JT, JI, JW, JR

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:INIT',0,ZHOOK_HANDLE)

  IF (KINTERVALS == 1) THEN
    ! We can use Stefan-Boltzmann law without a look-up table
    WRITE(NULOUT,'(a)') 'YOE_SPECTRAL_PLANCK: Single-band look-up table requested: use Stefan-Boltzmann law'
    SELF%NINTERVALS = KINTERVALS
    CALL SELF%FREE
  ELSE
    ! Full look-up table required
    ZCOEFF1 = 2.0_JPRB * RPI * RHPLA * RCLUM * RCLUM
    ZCOEFF2 = RHPLA * RCLUM / RKBOL

    NRANGES = SIZE(KINTERVAL_MAP,1)
    IF (SIZE(PWAVLEN_BOUND,1) /= NRANGES-1) THEN
      CALL ABOR1('YOS_SPECTRAL_PLANCK:INIT: PWAVLEN_BOUND must have one fewer elements than KINTERVAL_MAP')
    ENDIF

    CALL SELF%FREE

    ALLOCATE(SELF%WAVLEN_BOUND(NRANGES-1))
    ALLOCATE(SELF%INTERVAL_MAP(NRANGES))
    SELF%WAVLEN_BOUND(1:NRANGES-1) = PWAVLEN_BOUND(1:NRANGES-1)
    SELF%INTERVAL_MAP(1:NRANGES)   = KINTERVAL_MAP(1:NRANGES)

    SELF%NINTERVALS = KINTERVALS
    ! Temperature in 1-K intervals from 150 K to 350 K
    SELF%TEMP1  = 150.0_JPRB
    SELF%DTEMP  = 1.0_JPRB
    SELF%NTEMPS = 1 + NINT((350.0_JPRB - SELF%TEMP1) / SELF%DTEMP)

    ALLOCATE(SELF%PLANCK_LUT(SELF%NINTERVALS,SELF%NTEMPS))
    SELF%PLANCK_LUT(:,:) = 0.0_JPRB

    ! Print the properties of the look-up table
    WRITE(NULOUT,'(a,i0,a,f5.1,a,f5.1,a)') &
         &  'YOE_SPECTRAL_PLANCK: Generating Planck look-up table with ', &
         &  SELF%NTEMPS, ' temperatures from ', &
         &  SELF%TEMP1, ' to ', SELF%TEMP1+SELF%DTEMP*(SELF%NTEMPS-1), ' K:'
    DO JI = 1,SELF%NINTERVALS
      WRITE(NULOUT,'(a,i0,a)',advance='no') '  Band ', JI, ':'
      DO JR = 1,NRANGES
        IF (KINTERVAL_MAP(JR) == JI) THEN
          IF (JR == 1) THEN
            WRITE(NULOUT,'(a,f0.2)',advance='no') ' 0.00-', &
                 &  PWAVLEN_BOUND(1)*1.0e6_JPRB
          ELSEIF (JR == NRANGES) THEN
            WRITE(NULOUT,'(a,f0.2,a)',advance='no') ' ', &
                 &  PWAVLEN_BOUND(JR-1)*1.0e6_JPRB, '-Inf'
          ELSE
            WRITE(NULOUT,'(a,f0.2,a,f0.2)',advance='no') ' ', &
                 &  PWAVLEN_BOUND(JR-1)*1.0e6_JPRB, '-', &
                 &  PWAVLEN_BOUND(JR)*1.0e6_JPRB
          ENDIF
        ENDIF
      ENDDO
      WRITE(NULOUT,'(a)') ' microns'
    ENDDO

    ! Create the look-up table
    DO JT = 1,SELF%NTEMPS

      ZTEMP = SELF%TEMP1 + (JT-1) * SELF%DTEMP

      DO JI = 1,NRANGES

        IF (JI == 1) THEN
          ZWAVLEN1 = MIN(1.0E-6_JPRB, 0.8_JPRB * PWAVLEN_BOUND(1))
          ZWAVLEN2 = PWAVLEN_BOUND(1)
        ELSEIF (JI == NRANGES) THEN
          ZWAVLEN1 = PWAVLEN_BOUND(NRANGES-1)
          ! Simulate up to at least 200 microns wavelength
          ZWAVLEN2 = MAX(200.0E-6_JPRB, PWAVLEN_BOUND(NRANGES-1)+20.0E-6_JPRB)
        ELSE
          ZWAVLEN1 = PWAVLEN_BOUND(JI-1)
          ZWAVLEN2 = PWAVLEN_BOUND(JI)
        ENDIF

        NWAVLEN = 100
        DWAVLEN = (ZWAVLEN2 - ZWAVLEN1) / NWAVLEN
        ZSUM = 0.0_JPRB
        DO JW = 0,NWAVLEN
          ZWAVLEN = ZWAVLEN1 + DWAVLEN*JW
          ! Weights for trapezoidal rule
          !IF (JW > 0 .AND. JW < NWAVLEN) THEN
          !  ZWEIGHT = 2.0_JPRB
          !ELSE
          !  ZWEIGHT = 1.0_JPRB
          !ENDIF
          ! Weights for Simpson's rule
          IF (JW > 0 .AND. JW < NWAVLEN) THEN
            ZWEIGHT = 2.0_JPRB + 2.0_JPRB * MOD(JW,2)
          ELSE
            ZWEIGHT = 1.0_JPRB
          ENDIF
          ! Planck's law
          !
          ! The exponential term is computed in double precision to avoid
          ! overflow. The final result should still be in the range of a single
          ! precision float.
          ZWAVLEN_SQR = ZWAVLEN*ZWAVLEN
          ZPLANCKEXP = EXP(REAL(ZCOEFF2, JPRD) &
                     &     / (REAL(ZWAVLEN, JPRD) * REAL(ZTEMP, JPRD)))
          ZSUM = ZSUM + ZWEIGHT / (ZWAVLEN_SQR*ZWAVLEN_SQR*ZWAVLEN &
               &  * (ZPLANCKEXP - 1.0_JPRB))
        ENDDO
        SELF%PLANCK_LUT(KINTERVAL_MAP(JI),JT) = SELF%PLANCK_LUT(KINTERVAL_MAP(JI),JT) &
             &  + ZCOEFF1 * ZSUM * DWAVLEN / 3.0_JPRB
      ENDDO

    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:INIT',1,ZHOOK_HANDLE)

END SUBROUTINE INIT


!-----------------------------------------------------------------------
! Calculate Planck function in spectral intervals from temperature
SUBROUTINE CALC(SELF, KIDIA, KFDIA, KLON, PTEMPERATURE, PPLANCK)

  USE YOMCST,   ONLY : RSIGMA

  CLASS(TSPECTRALPLANCK), INTENT(IN)  :: SELF
  ! Process columns KIDIA-KFDIA from total of KLON columns
  INTEGER(KIND=JPIM)    , INTENT(IN)  :: KIDIA, KFDIA, KLON
  ! Temperature in Kelvin
  REAL(KIND=JPRB)       , INTENT(IN)  :: PTEMPERATURE(KLON)
  ! Integrated Planck function as an irradiance, in W m-2
  REAL(KIND=JPRB)       , INTENT(OUT) :: PPLANCK(KLON,SELF%NINTERVALS)

  ! Column loop counter, index to temperature interval
  INTEGER(KIND=JPRB) :: JL, ITEMP

  ! Interpolation weight, highest temperature in look-up table
  REAL(KIND=JPRB) :: ZWEIGHT, ZTEMP2

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:CALC',0,ZHOOK_HANDLE)

  IF (SELF%NINTERVALS == 1) THEN
    ! Stefan-Boltzmann law
    PPLANCK(KIDIA:KFDIA,1) = RSIGMA * PTEMPERATURE(KIDIA:KFDIA)**4

  ELSE
    ! Look-up table
    ZTEMP2 = SELF%TEMP1 + SELF%DTEMP * (SELF%NTEMPS - 1)

    DO JL = KIDIA,KFDIA

      IF (PTEMPERATURE(JL) <= SELF%TEMP1) THEN
        ! Cap the Planck function at the low end
        ITEMP   = 1
        ZWEIGHT = 0.0_JPRB
      ELSEIF (PTEMPERATURE(JL) < ZTEMP2) THEN
        ! Linear interpolation
        ZWEIGHT = 1.0_JPRB + (PTEMPERATURE(JL) - SELF%TEMP1) / SELF%DTEMP
        ITEMP   = NINT(ZWEIGHT)
        ZWEIGHT = ZWEIGHT - ITEMP
      ELSE
        ! Linear extrapolation at high temperatures off the scale
        ITEMP   = SELF%NTEMPS-1
        ZWEIGHT = 1.0_JPRB + (PTEMPERATURE(JL) - SELF%TEMP1) / SELF%DTEMP - ITEMP
      ENDIF

      PPLANCK(JL,:) = SELF%PLANCK_LUT(:,ITEMP) &
           &  + ZWEIGHT * (SELF%PLANCK_LUT(:,ITEMP+1) - SELF%PLANCK_LUT(:,ITEMP))

      ! Force sum to equal Stefan-Boltzmann law
      PPLANCK(JL,:) = PPLANCK(JL,:) * RSIGMA * PTEMPERATURE(JL)**4 / SUM(PPLANCK(JL,:),1)

    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:CALC',1,ZHOOK_HANDLE)

END SUBROUTINE CALC


!-----------------------------------------------------------------------
! Print look-up table to a unit
SUBROUTINE PRINT_SPECTRAL_PLANCK(SELF, IUNIT)

  CLASS(TSPECTRALPLANCK), INTENT(IN) :: SELF
  INTEGER(KIND=JPIM),     INTENT(IN) :: IUNIT

  INTEGER(KIND=JPIM) :: JT

  CHARACTER(len=24)  :: MY_FORMAT
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:PRINT_SPECTRAL_PLANCK',0,ZHOOK_HANDLE)
  IF (SELF%NINTERVALS == 1) THEN

    WRITE(IUNIT,'(A)') 'Spectral Planck in only one interval: using Stefan-Boltzmann law'

  ELSE

    WRITE(IUNIT,'(A,I0,A)') 'Spectral Planck look-up table defined in ', &
         &  SELF%NINTERVALS, ' intervals:'
    WRITE(MY_FORMAT,'(A,I0,A)') '(f7.2,', SELF%NINTERVALS, 'e15.5)'
    DO JT = 1,SELF%NTEMPS
      WRITE(IUNIT,TRIM(MY_FORMAT)) SELF%TEMP1 + (JT-1) * SELF%DTEMP, &
           &  SELF%PLANCK_LUT(:,JT)
    ENDDO

  ENDIF
IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:PRINT_SPECTRAL_PLANCK',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_SPECTRAL_PLANCK


!-----------------------------------------------------------------------
! Free allocated memory
SUBROUTINE FREE_SPECTRAL_PLANCK(SELF)

  CLASS(TSPECTRALPLANCK), INTENT(INOUT) :: SELF
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:FREE_SPECTRAL_PLANCK',0,ZHOOK_HANDLE)
  IF (ALLOCATED(SELF%PLANCK_LUT)) THEN
    DEALLOCATE(SELF%PLANCK_LUT)
  ENDIF
  IF (ALLOCATED(SELF%WAVLEN_BOUND)) THEN
    DEALLOCATE(SELF%WAVLEN_BOUND)
  ENDIF
  IF (ALLOCATED(SELF%INTERVAL_MAP)) THEN
    DEALLOCATE(SELF%INTERVAL_MAP)
  ENDIF
IF (LHOOK) CALL DR_HOOK('YOE_SPECTRAL_PLANCK:FREE_SPECTRAL_PLANCK',1,ZHOOK_HANDLE)

END SUBROUTINE FREE_SPECTRAL_PLANCK


END MODULE YOE_SPECTRAL_PLANCK
