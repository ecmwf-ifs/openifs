! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE EC_DATETIME_MOD

!**** Interface to eclib's date-time routines

!     Purpose.
!     --------
!     Fortran 90 Interface to calling eclib date-time routines

!     Author.
!     -------
!        W.Deconinck     ECMWF

!     Modifications.
!     --------------
!        Original: 2021-05-10

!     ------------------------------------------------------------------

IMPLICIT NONE
PRIVATE

PUBLIC :: cd2date, yd2date, idate2cd, idate2yd, icd2ymd, iymd2cd
PUBLIC :: daydiff, hourdiff, mindiff, secdiff
PUBLIC :: hourincr, minincr, secincr

! ISO-C-BINDING C-interfaces
INTERFACE

    ! void cd2date(const int32_t *const icd, int32_t *const iy, int32_t *const im, int32_t *const id, int32_t *const iret);
    subroutine cd2date( &
      & cd, &
      & year, month, day, &
      & iret) &
      & bind(c, name="cd2date")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in)  :: cd
      integer(c_int32_t), intent(out) :: year, month, day
      integer(c_int32_t), intent(out) :: iret
    end subroutine

    ! void yd2date(const int32_t *const iyd, const int32_t *const iy, int32_t *const im, int32_t *const id, int32_t *const iret);
    subroutine yd2date( &
      & yearday, year, &
      & month, day, &
      & iret) &
      & bind(c, name="yd2date")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in)  :: yearday, year
      integer(c_int32_t), intent(out) :: month, day
      integer(c_int32_t), intent(out) :: iret
    end subroutine

    ! int32_t idate2cd(const int32_t *const iy, const int32_t *const im, const int32_t *const id, int32_t *const iret);
    function idate2cd( &
      & year, month, day, &
      & iret) &
      & bind(c, name="idate2cd")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t) :: idate2cd
      integer(c_int32_t), intent(in)  :: year, month, day
      integer(c_int32_t), intent(out) :: iret
    end function

    ! int32_t idate2yd(const int32_t *const iy, const int32_t *const im, const int32_t *const id, int32_t *const iret);
    function idate2yd( &
      & year, month, day, &
      & iret) &
      & bind(c, name="idate2yd")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t) :: idate2yd
      integer(c_int32_t), intent(in)  :: year, month, day
      integer(c_int32_t), intent(out) :: iret
    end function

    ! int32_t icd2ymd(const int32_t *const icd, int32_t *const iret);
    function icd2ymd( &
      & cd, &
      & iret) &
      & bind(c, name="icd2ymd")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t) :: icd2ymd
      integer(c_int32_t), intent(in)  :: cd
      integer(c_int32_t), intent(out) :: iret
    end function

    ! int32_t  iymd2cd(const int32_t *const iymd, int32_t *const iret);
    function iymd2cd( &
      & ymd, &
      & iret) &
      & bind(c, name="iymd2cd")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t) :: iymd2cd
      integer(c_int32_t), intent(in)  :: ymd
      integer(c_int32_t), intent(out) :: iret
    end function

    ! void daydiff(const int32_t *const year1, const int32_t *const month1, const int32_t *const day1,
    !              const int32_t *const year2, const int32_t *const month2, const int32_t *const day2,
    !              int32_t *const days, int32_t *const iret);
    subroutine daydiff( &
      & year1, month1, day1, &
      & year2, month2, day2, &
      & days, &
      & iret) &
      & bind(c, name="daydiff")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in)  :: year1, month1, day1
      integer(c_int32_t), intent(in)  :: year2, month2, day2
      integer(c_int32_t), intent(out) :: days, iret
    end subroutine

    ! void hourdiff(const int32_t *const year1, const int32_t *const month1, const int32_t *const day1, const int32_t *const hour1,
    !               const int32_t *const year2, const int32_t *const month2, const int32_t *const day2, const int32_t *const hour2,
    !               int32_t *const hours, int32_t *const iret);
    subroutine hourdiff( &
      & year1, month1, day1, hour1, &
      & year2, month2, day2, hour2, &
      & hours, &
      & iret) &
      & bind(c, name="hourdiff")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in)  :: year1, month1, day1, hour1
      integer(c_int32_t), intent(in)  :: year2, month2, day2, hour2
      integer(c_int32_t), intent(out) :: hours, iret
    end subroutine

    ! void mindiff(const int32_t *const year1, const int32_t *const month1, const int32_t *const day1,
    !              const int32_t *const hour1, const int32_t *const min1,
    !              const int32_t *const year2, const int32_t *const month2, const int32_t *const day2,
    !              const int32_t *const hour2, const int32_t *const min2,
    !              int32_t *const minutes, int32_t *const iret);
    subroutine mindiff( &
        & year1, month1, day1, hour1, min1, &
        & year2, month2, day2, hour2, min2, &
        & minutes, &
        & iret) &
        & bind(c, name="mindiff")
        use, intrinsic :: iso_c_binding
        integer(c_int32_t), intent(in)  :: year1, month1, day1, hour1, min1
        integer(c_int32_t), intent(in)  :: year2, month2, day2, hour2, min2
        integer(c_int32_t), intent(out) :: minutes, iret
    end subroutine

    ! void secdiff(const int32_t *const year1, const int32_t *const month1, const int32_t *const day1, const int32_t *const hour1,
    !              const int32_t *const min1, const int32_t *const sec1,
    !              const int32_t *const year2, const int32_t *const month2, const int32_t *const day2, const int32_t *const hour2,
    !              const int32_t *const min2, const int32_t *const sec2,
    !              int32_t *const seconds, int32_t *const iret);
    subroutine secdiff( &
        & year1, month1, day1, hour1, min1, sec1, &
        & year2, month2, day2, hour2, min2, sec2, &
        & seconds, &
        & iret) &
        & bind(c, name="secdiff")
        use, intrinsic :: iso_c_binding
        integer(c_int32_t), intent(in)  :: year1, month1, day1, hour1, min1, sec1
        integer(c_int32_t), intent(in)  :: year2, month2, day2, hour2, min2, sec2
        integer(c_int32_t), intent(out) :: seconds, iret
    end subroutine

    ! void hourincr(const int32_t *const year, const int32_t *const month, const int32_t *const day, const int32_t *const hour,
    !               const int32_t *const hours,
    !               int32_t *const new_year, int32_t *const new_month, int32_t *const new_day, int32_t *const new_hour,
    !               int32_t *const iret);
    subroutine hourincr( &
      & year, month, day, hour, &
      & hours, &
      & new_year, new_month, new_day, new_hour, &
      & iret) &
      & bind(c, name="hourincr")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in)  :: year, month, day, hour, hours
      integer(c_int32_t), intent(out) :: new_year, new_month, new_day, new_hour, iret
    end subroutine

    ! void minincr(const int32_t *const year, const int32_t *const month, const int32_t *const day, const int32_t *const hour,
    !              const int32_t *const min,
    !              const int32_t *const minutes,int32_t *const new_year, int32_t *const new_month, int32_t *const new_day,
    !              int32_t *const new_hour, int32_t *const new_min,  int32_t *const iret);
    subroutine minincr( &
      & year, month, day, hour, min, &
      & minutes, &
      & new_year, new_month, new_day, new_hour, new_min, &
      & iret) &
      & bind(c, name="minincr")
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in)  :: year, month, day, hour, min, minutes
      integer(c_int32_t), intent(out) :: new_year, new_month, new_day, new_hour, new_min, iret
    end subroutine

    ! void secincr(const int32_t *const year, const int32_t *const month, const int32_t *const day, const int32_t *const hour,
    !              const int32_t *const min, const int32_t *const sec, const int32_t *const seconds,
    !              int32_t *const new_year, int32_t *const new_month, int32_t *const new_day, int32_t *const new_hour,
    !              int32_t *const new_min, int32_t *const new_sec,  int32_t *const iret);
    subroutine secincr( &
        & year, month, day, hour, min, sec, &
        & seconds, &
        & new_year, new_month, new_day, new_hour, new_min, new_sec, &
        & iret) &
        & bind(c, name="secincr")
        use, intrinsic :: iso_c_binding
        integer(c_int32_t), intent(in)  :: year, month, day, hour, min, sec, seconds
        integer(c_int32_t), intent(out) :: new_year, new_month, new_day, new_hour, new_min, new_sec, iret
    end subroutine

END INTERFACE

CONTAINS

END MODULE
