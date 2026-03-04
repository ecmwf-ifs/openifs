! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COS_SZA(YDERAD,YDERIP,YDRIP,KSTART,KEND,KCOL,PGEMU,PGELAM,LDRADIATIONTIMESTEP,PMU0)

!**** *COS_SZA*   

!     Purpose.
!     --------
!        Compute the cosine of the solar zenith angle.  Note that this
!        is needed for three different things: (1) as input to the
!        radiation scheme in which it is used to compute the path
!        length of the direct solar beam through the atmosphere, (2)
!        every timestep to scale the solar fluxes by the incoming
!        solar radiation at top-of-atmosphere, and (3) to compute the
!        albedo of the ocean.  For (1) we ideally want an average
!        value for the duration of a radiation timestep while for (2)
!        we want an average value for the duration of a model
!        timestep.

!**   Interface.
!     ----------
!        *CALL* *COS_SZA(...)

!        Explicit arguments : 
!        ------------------
!            PGEMU - Sine of latitude
!            PGELAM - Geographic longitude in radians
!            LDRadiationTimestep - Is this for a radiation timestep?
!            PMU0 - Output cosine of solar zenith angle

!        Implicit arguments :
!        --------------------
!            YRRIP%RWSOVR, RWSOVRM - Solar time for model/radiation timesteps
!            RCODECM, RSIDECM - Sine/cosine of solar declination
!            YRERAD%LAverageSZA - Average solar zenith angle in time interval?
!            YRRIP%TSTEP - Model timestep in seconds
!            YRERAD%NRADFR - Radiation frequency in timesteps

!     Method.
!     -------
!        Compute cosine of the solar zenith angle, mu0, from lat, lon
!        and solar time using standard formula.  If
!        YRERAD%LAverageSZA=FALSE then this is done at a single time,
!        which is assumed to be the mid-point of either the model or
!        the radiation timestep.  If YRERAD%LAverageSZA=TRUE then we
!        compute the average over the model timestep exactly by first
!        computing sunrise/sunset times. For radiation timesteps, mu0
!        is to be used to compute the path length of the direct solar
!        beam through the atmosphere, and the fluxes are subsequently
!        weighted by mu0.  Therefore night-time values are not used,
!        so we average mu0 only when the sun is above the horizon.

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!
!        See also: Zhou, L., M. Zhang, Q. Bao, and Y. Liu (2015), On
!        the incident solar radiation in CMIP5
!        models. Geophys. Res. Lett., 42, 1930–1935. doi:
!        10.1002/2015GL063239.

!     Author.
!     -------
!      Robin Hogan, ECMWF, May 2015

!     Modifications:
!     --------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI, RDAY
USE YOMRIP   , ONLY : TRIP
USE YOERIP   , ONLY : TERIP
USE YOERAD   , ONLY : TERAD
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(IN) :: YDERAD
TYPE(TERIP)       ,INTENT(IN) :: YDERIP
TYPE(TRIP)        ,INTENT(IN) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN) :: KSTART      ! Start column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KEND        ! Last column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KCOL        ! Number of columns in arrays
REAL(KIND=JPRB),   INTENT(IN) :: PGEMU(KCOL) ! Sine of latitude
REAL(KIND=JPRB),   INTENT(IN) :: PGELAM(KCOL)! Longitude in radians
LOGICAL, INTENT(IN) :: LDRADIATIONTIMESTEP   ! Is this for a radiation timestep?
REAL(KIND=JPRB),  INTENT(OUT) :: PMU0(KCOL)  ! Cosine of solar zenith angle

! Solar time at the start and end of the time interval
REAL(KIND=JPRB) :: ZSOLARTIMESTART, ZSOLARTIMEEND

! The time of half a model/radiation timestep, in radians
REAL(KIND=JPRB) :: ZHALFTIMESTEP

! For efficiency we precompute sin(solar declination)*sin(latitude)
REAL(KIND=JPRB) :: ZSINDECSINLAT(KSTART:KEND)
!...and cos(solar declination)*cos(latitude)
REAL(KIND=JPRB) :: ZCOSDECCOSLAT(KSTART:KEND)
! ...and cosine of latitude
REAL(KIND=JPRB) :: ZCOSLAT(KSTART:KEND)

! Tangent of solar declination
REAL(KIND=JPRB) :: ZTANDEC

! Hour angles (=local solar time in radians plus pi)
REAL(KIND=JPRB) :: ZHOURANGLESTART, ZHOURANGLEEND
REAL(KIND=JPRB) :: ZHOURANGLESUNSET, ZCOSHOURANGLESUNSET

INTEGER(KIND=JPIM) :: JCOL        ! Column index

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COS_SZA',0,ZHOOK_HANDLE)

! An average solar zenith angle can only be computed if the solar time
! is centred on the time interval
IF (YDERAD%LAVERAGESZA .AND. .NOT. YDERAD%LCENTREDTIMESZA) THEN
  WRITE(NULOUT,*) 'ERROR IN COS_SZA: LAverageSZA=TRUE but LCentredTimeSZA=FALSE'
  CALL ABOR1('COS_SZA: ABOR1 CALLED')
ENDIF

DO JCOL = KSTART,KEND
  ZCOSLAT(JCOL) = SQRT(1.0_JPRB - PGEMU(JCOL)**2)
ENDDO

IF (LDRADIATIONTIMESTEP) THEN
  ! Compute the effective cosine of solar zenith angle for a radiation
  ! timestep

  ! Precompute quantities that may be used more than once
  DO JCOL = KSTART,KEND
    ZSINDECSINLAT(JCOL) = YDERIP%RSIDECM * PGEMU(JCOL)
    ZCOSDECCOSLAT(JCOL) = YDERIP%RCODECM * ZCOSLAT(JCOL)
  ENDDO

  IF (.NOT. YDERAD%LAVERAGESZA) THEN
    ! Original method: compute the value at the centre of the
    ! radiation timestep (assuming that LCentredTimeSZA=TRUE - see
    ! updtim.F90)
    DO JCOL = KSTART,KEND
      ! It would be more efficient to do it like this...
      ! PMU0(JCOL)=MAX(0.0_JPRB, ZSinDecSinLat(JCOL) &
      !      & - ZCosDecCosLat(JCOL) * COS(YRERIP%RWSOVRM + PGELAM(JCOL)))
      ! ...but for bit reproducibility with previous cycle we do it
      ! like this:
      PMU0(JCOL) = MAX(0.0_JPRB, ZSINDECSINLAT(JCOL) &
           & - YDERIP%RCODECM*COS(YDERIP%RWSOVRM)*ZCOSLAT(JCOL)*COS(PGELAM(JCOL)) &
           & + YDERIP%RCODECM*SIN(YDERIP%RWSOVRM)*ZCOSLAT(JCOL)*SIN(PGELAM(JCOL)))
    ENDDO

  ELSE
    ! Compute the average MU0 for the period of the radiation
    ! timestep, excluding times when the sun is below the horizon

    ! First compute the sine and cosine of the times of the start and
    ! end of the radiation timestep
    ZHALFTIMESTEP = YDRIP%TSTEP * REAL(YDERAD%NRADFR) * RPI / RDAY
    ZSOLARTIMESTART = YDERIP%RWSOVRM - ZHALFTIMESTEP
    ZSOLARTIMEEND   = YDERIP%RWSOVRM + ZHALFTIMESTEP

    ! Compute tangent of solar declination, with check in case someone
    ! simulates a planet completely tipped over
    ZTANDEC = YDERIP%RSIDECM / MAX(YDERIP%RCODECM, 1.0E-12_JPRB)

    DO JCOL = KSTART,KEND
      ! Sunrise equation: cos(hour angle at sunset) =
      ! -tan(declination)*tan(latitude)
      ZCOSHOURANGLESUNSET = -ZTANDEC * PGEMU(JCOL)&
           &              / MAX(ZCOSLAT(JCOL), 1.0E-12_JPRB)
      IF (ZCOSHOURANGLESUNSET > 1.0) THEN
        ! Perpetual darkness
        PMU0(JCOL) = 0.0_JPRB
      ELSE
        ! Compute hour angle at start and end of time interval,
        ! ensuring that the hour angle of the centre of the time
        ! window is in the range -PI to +PI (equivalent to ensuring
        ! that local solar time = solar time + longitude is in the
        ! range 0 to 2PI)
        IF (YDERIP%RWSOVRM + PGELAM(JCOL) < 2.0_JPRB*RPI) THEN
          ZHOURANGLESTART = ZSOLARTIMESTART + PGELAM(JCOL) - RPI
          ZHOURANGLEEND   = ZSOLARTIMEEND   + PGELAM(JCOL) - RPI 
        ELSE
          ZHOURANGLESTART = ZSOLARTIMESTART + PGELAM(JCOL) - 3.0_JPRB*RPI
          ZHOURANGLEEND   = ZSOLARTIMEEND   + PGELAM(JCOL) - 3.0_JPRB*RPI
        ENDIF

        IF (ZCOSHOURANGLESUNSET >= -1.0) THEN
          ! Not perpetual daylight or perpetual darkness, so we need
          ! to check for sunrise or sunset lying within the time
          ! interval
          ZHOURANGLESUNSET = ACOS(ZCOSHOURANGLESUNSET)
          IF (ZHOURANGLEEND <= -ZHOURANGLESUNSET &
               & .OR. ZHOURANGLESTART >= ZHOURANGLESUNSET) THEN
            ! The time interval is either completely before sunrise or
            ! completely after sunset
            PMU0(JCOL) = 0.0_JPRB
            CYCLE
          ENDIF

          ! Bound the start and end hour angles by sunrise and sunset
          ZHOURANGLESTART = MAX(-ZHOURANGLESUNSET, &
               &                MIN(ZHOURANGLESTART, ZHOURANGLESUNSET))
          ZHOURANGLEEND   = MAX(-ZHOURANGLESUNSET, &
               &                MIN(ZHOURANGLEEND,   ZHOURANGLESUNSET))
        ENDIF

        IF (ZHOURANGLEEND - ZHOURANGLESTART > 1.0E-8) THEN
          ! Compute average MU0 in the interval ZHourAngleStart to
          ! ZHourAngleEnd
          PMU0(JCOL) = ZSINDECSINLAT(JCOL) &
               & + (ZCOSDECCOSLAT(JCOL) &
               & * (SIN(ZHOURANGLEEND) - SIN(ZHOURANGLESTART))) &
               & / (ZHOURANGLEEND - ZHOURANGLESTART)

          ! Just in case...
          IF (PMU0(JCOL) < 0.0_JPRB) THEN
            PMU0(JCOL) = 0.0_JPRB
          ENDIF
        ELSE
          ! Too close to sunrise/sunset for a reliable calculation
          PMU0(JCOL) = 0.0_JPRB
        ENDIF

      ENDIF
    ENDDO
  ENDIF

ELSE
  ! Compute the cosine of solar zenith angle for a model timestep

  ! Precompute quantities that may be used more than once
  DO JCOL = KSTART,KEND
    ZSINDECSINLAT(JCOL) = YDRIP%RSIDEC * PGEMU(JCOL)
    ZCOSDECCOSLAT(JCOL) = YDRIP%RCODEC * ZCOSLAT(JCOL)
  ENDDO

  IF (.NOT. YDERAD%LAVERAGESZA) THEN
    ! Original method: compute the value at the centre of the
    ! model timestep
    DO JCOL = KSTART,KEND
      ! It would be more efficient to do it like this...
      ! PMU0(JCOL) = MAX(0.0_JPRB, ZSinDecSinLat(JCOL)        &
      !      & - ZCosDecCosLat(JCOL)*COS(YRRIP%RWSOVR + PGELAM(JCOL)))
      ! ...but for bit reproducibility with previous cycle we do it
      ! like this:
      PMU0(JCOL) = MAX(0.0_JPRB, ZSINDECSINLAT(JCOL)&
           & - YDRIP%RCODEC*COS(YDRIP%RWSOVR)*ZCOSLAT(JCOL)*COS(PGELAM(JCOL))&
           & + YDRIP%RCODEC*SIN(YDRIP%RWSOVR)*ZCOSLAT(JCOL)*SIN(PGELAM(JCOL)))
    ENDDO

  ELSE
    ! Compute the average MU0 for the period of the model timestep

    ! First compute the sine and cosine of the times of the start and
    ! end of the model timestep
    ZHALFTIMESTEP   = YDRIP%TSTEP * RPI / RDAY
    ZSOLARTIMESTART = YDRIP%RWSOVR - ZHALFTIMESTEP
    ZSOLARTIMEEND   = YDRIP%RWSOVR + ZHALFTIMESTEP

    ! Compute tangent of solar declination, with check in case someone
    ! simulates a planet completely tipped over
    ZTANDEC = YDRIP%RSIDEC / MAX(YDRIP%RCODEC, 1.0E-12_JPRB)

    DO JCOL = KSTART,KEND
      ! Sunrise equation: cos(hour angle at sunset) =
      ! -tan(declination)*tan(latitude)
      ZCOSHOURANGLESUNSET = -ZTANDEC * PGEMU(JCOL) &
           &              / MAX(ZCOSLAT(JCOL), 1.0E-12_JPRB)
      IF (ZCOSHOURANGLESUNSET > 1.0) THEN
        ! Perpetual darkness
        PMU0(JCOL) = 0.0_JPRB
      ELSE
        ! Compute hour angle at start and end of time interval,
        ! ensuring that the hour angle of the centre of the time
        ! window is in the range -PI to +PI (equivalent to ensuring
        ! that local solar time = solar time + longitude is in the
        ! range 0 to 2PI)
        IF (YDRIP%RWSOVR + PGELAM(JCOL) < 2.0_JPRB*RPI) THEN
          ZHOURANGLESTART = ZSOLARTIMESTART + PGELAM(JCOL) - RPI
          ZHOURANGLEEND   = ZSOLARTIMEEND   + PGELAM(JCOL) - RPI 
        ELSE
          ZHOURANGLESTART = ZSOLARTIMESTART + PGELAM(JCOL) - 3.0_JPRB*RPI
          ZHOURANGLEEND   = ZSOLARTIMEEND   + PGELAM(JCOL) - 3.0_JPRB*RPI
        ENDIF

        IF (ZCOSHOURANGLESUNSET >= -1.0) THEN
          ! Not perpetual daylight or perpetual darkness, so we need
          ! to check for sunrise or sunset lying within the time
          ! interval
          ZHOURANGLESUNSET = ACOS(ZCOSHOURANGLESUNSET)
          IF (ZHOURANGLEEND <= -ZHOURANGLESUNSET &
               & .OR. ZHOURANGLESTART >= ZHOURANGLESUNSET) THEN
            ! The time interval is either completely before sunrise or
            ! completely after sunset
            PMU0(JCOL) = 0.0_JPRB
            CYCLE
          ENDIF

          ! Bound the start and end hour angles by sunrise and sunset
          ZHOURANGLESTART = MAX(-ZHOURANGLESUNSET, &
               &                MIN(ZHOURANGLESTART, ZHOURANGLESUNSET))
          ZHOURANGLEEND   = MAX(-ZHOURANGLESUNSET, &
               &                MIN(ZHOURANGLEEND,   ZHOURANGLESUNSET))
        ENDIF

        ! Compute average MU0 in the model timestep, although the
        ! numerator considers only the time from ZHourAngleStart to
        ! ZHourAngleEnd that the sun is above the horizon
        PMU0(JCOL) = (ZSINDECSINLAT(JCOL) * (ZHOURANGLEEND-ZHOURANGLESTART)   &
           & + ZCOSDECCOSLAT(JCOL)*(SIN(ZHOURANGLEEND)-SIN(ZHOURANGLESTART))) &
           & / (2.0_JPRB * ZHALFTIMESTEP)

        ! This shouldn't ever result in negative values, but just in
        ! case
        IF (PMU0(JCOL) < 0.0_JPRB) THEN
          PMU0(JCOL) = 0.0_JPRB
        ENDIF

      ENDIF
    ENDDO
  ENDIF

ENDIF


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COS_SZA',1,ZHOOK_HANDLE)
END SUBROUTINE COS_SZA
