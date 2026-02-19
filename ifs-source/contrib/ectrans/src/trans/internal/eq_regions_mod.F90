! (C) Copyright 2006- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE eq_regions_mod
!
!     Purpose.
!     --------
!           eq_regions_mod provides the code to perform a high level
!           partitioning of the surface of a sphere into regions of
!           equal area and small diameter.
!           the type.
!
!     Background.
!     -----------
!     This Fortran version of eq_regions is a much cut down version of the
!     "Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox" of the
!     same name developed by Paul Leopardi at the University of New South Wales.
!     This version has been coded specifically for the case of partitioning the
!     surface of a sphere or S^dim (where dim=2) as denoted in the original code.
!     Only a subset of the original eq_regions package has been coded to determine
!     the high level distribution of regions on a sphere, as the detailed
!     distribution of grid points to each region is left to IFS software.
!     This is required to take into account the spatial distribution of grid
!     points in an IFS gaussian grid and provide an optimal (i.e. exact)
!     distribution of grid points over regions.
!
!     The following copyright notice for the eq_regions package is included from
!     the original MatLab release.
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     + Release 1.10 2005-06-26                                                 +
!     +                                                                         +
!     + Copyright (c) 2004, 2005, University of New South Wales                 +
!     +                                                                         +
!     + Permission is hereby granted, free of charge, to any person obtaining   +
!     + a copy of this software and associated documentation files (the         +
!     + "Software"), to deal in the Software without restriction, including     +
!     + without limitation the rights to use, copy, modify, merge, publish,     +
!     + distribute, sublicense, and/or sell copies of the Software, and to      +
!     + permit persons to whom the Software is furnished to do so, subject to   +
!     + the following conditions:                                               +
!     +                                                                         +
!     + The above copyright notice and this permission notice shall be included +
!     + in all copies or substantial portions of the Software.                  +
!     +                                                                         +
!     + THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         +
!     + EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      +
!     + MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  +
!     + IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    +
!     + CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    +
!     + TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       +
!     + SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  +
!     +                                                                         +
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Author.
!     -------
!        George Mozdzynski *ECMWF*
!
!     Modifications.
!     --------------
!        Original : 2006-04-15
!
!--------------------------------------------------------------------------------
!
USE PARKIND1  ,ONLY : JPIM,   JPRB

IMPLICIT NONE

SAVE

PRIVATE

PUBLIC eq_regions,l_regions_debug,n_regions_ns,n_regions_ew,n_regions,my_region_ns,my_region_ew
PUBLIC eq_regions_t, eq_regions_save, eq_regions_load, eq_regions_free

real(kind=jprb) pi

type eq_regions_t
logical :: l_regions_debug=.false.
integer(kind=jpim) :: n_regions_ns
integer(kind=jpim) :: n_regions_ew
integer(kind=jpim) :: my_region_ns
integer(kind=jpim) :: my_region_ew
integer(kind=jpim),pointer :: n_regions(:) => null ()
end type eq_regions_t

logical :: l_regions_debug=.false.
integer(kind=jpim) :: n_regions_ns
integer(kind=jpim) :: n_regions_ew
integer(kind=jpim) :: my_region_ns
integer(kind=jpim) :: my_region_ew
integer(kind=jpim),pointer :: n_regions(:) => null ()

CONTAINS

subroutine eq_regions_save (yder)
type (eq_regions_t), intent (inout) :: yder

yder%l_regions_debug =  l_regions_debug 
yder%n_regions_ns    =  n_regions_ns    
yder%n_regions_ew    =  n_regions_ew    
yder%my_region_ns    =  my_region_ns    
yder%my_region_ew    =  my_region_ew    
yder%n_regions       => n_regions

nullify (n_regions)

end subroutine

subroutine eq_regions_load (yder)
type (eq_regions_t), intent (inout) :: yder

l_regions_debug =  yder%l_regions_debug 
n_regions_ns    =  yder%n_regions_ns    
n_regions_ew    =  yder%n_regions_ew    
my_region_ns    =  yder%my_region_ns    
my_region_ew    =  yder%my_region_ew    
n_regions       => yder%n_regions

nullify (yder%n_regions)

end subroutine

subroutine eq_regions_free (yder)
type (eq_regions_t), intent (inout) :: yder

if (associated (yder%n_regions)) then
  deallocate (yder%n_regions)
  nullify (yder%n_regions)
endif

end subroutine

subroutine eq_regions(N)
!
! eq_regions uses the zonal equal area sphere partitioning algorithm to partition
! the surface of a sphere into N regions of equal area and small diameter.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
integer(kind=jpim),intent(in) :: N
integer(kind=jpim) :: n_collars,j
real(kind=jprb),allocatable :: r_regions(:)
real(kind=jprb) :: c_polar

pi=2.0_jprb*asin(1.0_jprb)

n_regions(:)=0

if( N == 1 )then

  !
  ! We have only one region, which must be the whole sphere.
  !
  n_regions(1)=1
  n_regions_ns=1

else

  !
  ! Given N, determine c_polar
  ! the colatitude of the North polar spherical cap.
  !
  c_polar = polar_colat(N)
  !
  ! Given N, determine the ideal angle for spherical collars.
  ! Based on N, this ideal angle, and c_polar,
  ! determine n_collars, the number of collars between the polar caps.
  !
  n_collars = num_collars(N,c_polar,ideal_collar_angle(N))
  n_regions_ns=n_collars+2
  !
  ! Given N, c_polar and n_collars, determine r_regions,
  ! a list of the ideal real number of regions in each collar,
  ! plus the polar caps.
  ! The number of elements is n_collars+2.
  ! r_regions[1] is 1.
  ! r_regions[n_collars+2] is 1.
  ! The sum of r_regions is N.
  allocate(r_regions(n_collars+2))
  call ideal_region_list(N,c_polar,n_collars,r_regions)
  !
  ! Given N and r_regions, determine n_regions, a list of the natural number
  ! of regions in each collar and the polar caps.
  ! This list is as close as possible to r_regions.
  ! The number of elements is n_collars+2.
  ! n_regions[1] is 1.
  ! n_regions[n_collars+2] is 1.
  ! The sum of n_regions is N.
  !
  call round_to_naturals(N,n_collars,r_regions)
  deallocate(r_regions)
  if( N /= sum(n_regions(:)) )then
    write(*,'("eq_regions: N=",I10," sum(n_regions(:))=",I10)')N,sum(n_regions(:))
    call abor1('eq_regions: N /= sum(n_regions)')
  endif

endif

if( l_regions_debug )then
  write(*,'("eq_regions: N=",I6," n_regions_ns=",I4)') N,n_regions_ns
  do j=1,n_regions_ns
    write(*,'("eq_regions: n_regions(",I4,")=",I4)') j,n_regions(j)
  enddo
endif
n_regions_ew=maxval(n_regions(:))

return
end subroutine eq_regions

function num_collars(N,c_polar,a_ideal) result(num_c)
!
!NUM_COLLARS The number of collars between the polar caps
!
! Given N, an ideal angle, and c_polar,
! determine n_collars, the number of collars between the polar caps.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
integer(kind=jpim),intent(in) :: N
real(kind=jprb),intent(in) :: a_ideal,c_polar
integer(kind=jpim) :: num_c
logical enough
enough = (N > 2) .and. (a_ideal > 0)
if( enough )then
  num_c = max(1,nint((pi-2.*c_polar)/a_ideal))
else
  num_c = 0
endif
return
end function num_collars

subroutine ideal_region_list(N,c_polar,n_collars,r_regions)
!
!IDEAL_REGION_LIST The ideal real number of regions in each zone
!
! List the ideal real number of regions in each collar, plus the polar caps.
!
! Given N, c_polar and n_collars, determine r_regions, a list of the ideal real
! number of regions in each collar, plus the polar caps.
! The number of elements is n_collars+2.
! r_regions[1] is 1.
! r_regions[n_collars+2] is 1.
! The sum of r_regions is N.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
integer(kind=jpim),intent(in) :: N,n_collars
real(kind=jprb),intent(in) :: c_polar
real(kind=jprb),intent(out) :: r_regions(n_collars+2)
integer(kind=jpim) :: collar_n
real(kind=jprb) :: ideal_region_area,ideal_collar_area
real(kind=jprb) :: a_fitting
r_regions(:)=0.0_jprb
r_regions(1) = 1.0_jprb
if( n_collars > 0 )then
  !
  ! Based on n_collars and c_polar, determine a_fitting,
  ! the collar angle such that n_collars collars fit between the polar caps.
  !
  a_fitting = (pi-2.0_jprb*c_polar)/float(n_collars)
  ideal_region_area = area_of_ideal_region(N)
  do collar_n=1,n_collars
    ideal_collar_area = area_of_collar(c_polar+(collar_n-1)*a_fitting, &
     & c_polar+collar_n*a_fitting)
    r_regions(1+collar_n) = ideal_collar_area / ideal_region_area
  enddo
endif
r_regions(2+n_collars) = 1.
return
end subroutine ideal_region_list

function ideal_collar_angle(N) result(ideal)
!
! IDEAL_COLLAR_ANGLE The ideal angle for spherical collars of an EQ partition
!
! IDEAL_COLLAR_ANGLE(N) sets ANGLE to the ideal angle for the
! spherical collars of an EQ partition of the unit sphere S^2 into N regions.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
integer(kind=jpim),intent(in) :: N
real(kind=jprb) :: ideal
ideal = area_of_ideal_region(N)**(0.5_jprb)
return
end function ideal_collar_angle

subroutine round_to_naturals(N,n_collars,r_regions)
!
! ROUND_TO_NATURALS Round off a given list of numbers of regions
!
! Given N and r_regions, determine n_regions, a list of the natural number
! of regions in each collar and the polar caps.
! This list is as close as possible to r_regions, using rounding.
! The number of elements is n_collars+2.
! n_regions[1] is 1.
! n_regions[n_collars+2] is 1.
! The sum of n_regions is N.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
integer(kind=jpim),intent(in) :: N,n_collars
real(kind=jprb),intent(in) :: r_regions(n_collars+2)
integer(kind=jpim) :: zone_n
real(kind=jprb) :: discrepancy
n_regions(1:n_collars+2) = r_regions(:)
discrepancy = 0.0_jprb
do zone_n = 1,n_collars+2
    n_regions(zone_n) = nint(r_regions(zone_n)+discrepancy);
    discrepancy = discrepancy+r_regions(zone_n)-float(n_regions(zone_n));
enddo
return
end subroutine round_to_naturals

function polar_colat(N) result(polar_c)
!
! Given N, determine the colatitude of the North polar spherical cap.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
integer(kind=jpim),intent(in) :: N
real(kind=jprb) :: area
real(kind=jprb) :: polar_c
if( N == 1 ) polar_c=pi
if( N == 2 ) polar_c=pi/2.0_jprb
if( N > 2 )then
  area=area_of_ideal_region(N)
  polar_c=sradius_of_cap(area)
endif
return
end function polar_colat

function area_of_ideal_region(N) result(area)
!
! AREA_OF_IDEAL_REGION(N) sets AREA to be the area of one of N equal
! area regions on S^2, that is 1/N times AREA_OF_SPHERE.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
integer(kind=jpim),intent(in) :: N
real(kind=jprb) :: area_of_sphere
real(kind=jprb) :: area
area_of_sphere = (2.0_jprb*pi**1.5_jprb/gamma(1.5_jprb))
area = area_of_sphere/float(N)
return
end function area_of_ideal_region

function sradius_of_cap(area) result(sradius)
!
! SRADIUS_OF_CAP(AREA) returns the spherical radius of
! an S^2 spherical cap of area AREA.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
real(kind=jprb),intent(in) :: area
real(kind=jprb) :: sradius
sradius = 2.0_jprb*asin(sqrt(area/pi)/2.0_jprb)
return
end function sradius_of_cap

function area_of_collar(a_top, a_bot) result(area)
!
! AREA_OF_COLLAR Area of spherical collar
!
! AREA_OF_COLLAR(A_TOP, A_BOT) sets AREA to be the area of an S^2 spherical
! collar specified by A_TOP, A_BOT, where A_TOP is top (smaller) spherical radius,
! A_BOT is bottom (larger) spherical radius.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
real(kind=jprb),intent(in) :: a_top,a_bot
real(kind=jprb) area
area = area_of_cap(a_bot) - area_of_cap(a_top)
return
end function area_of_collar

function area_of_cap(s_cap) result(area)
!
! AREA_OF_CAP Area of spherical cap
!
! AREA_OF_CAP(S_CAP) sets AREA to be the area of an S^2 spherical
! cap of spherical radius S_CAP.
!
real(kind=jprb),intent(in) :: s_cap
real(kind=jprb) area
area = 4.0_jprb*pi * sin(s_cap/2.0_jprb)**2
return
end function area_of_cap

function gamma(x) result(gamma_res)
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
real(kind=jprb),intent(in) :: x
real(kind=jprb) :: gamma_res
real(kind=jprb) :: p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
real(kind=jprb) :: w,y
integer(kind=jpim) :: k,n
parameter (&
& p0 =   0.999999999999999990e+00_jprb,&
& p1 =  -0.422784335098466784e+00_jprb,&
& p2 =  -0.233093736421782878e+00_jprb,&
& p3 =   0.191091101387638410e+00_jprb,&
& p4 =  -0.024552490005641278e+00_jprb,&
& p5 =  -0.017645244547851414e+00_jprb,&
& p6 =   0.008023273027855346e+00_jprb)
parameter (&
& p7 =  -0.000804329819255744e+00_jprb,&
& p8 =  -0.000360837876648255e+00_jprb,&
& p9 =   0.000145596568617526e+00_jprb,&
& p10 = -0.000017545539395205e+00_jprb,&
& p11 = -0.000002591225267689e+00_jprb,&
& p12 =  0.000001337767384067e+00_jprb,&
& p13 = -0.000000199542863674e+00_jprb)
n = nint(x - 2)
w = x - (n + 2)
y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) *&
&    w + p9) * w + p8) * w + p7) * w + p6) * w + p5) *&
&    w + p4) * w + p3) * w + p2) * w + p1) * w + p0
if (n .gt. 0) then
  w = x - 1
  do k = 2, n
    w = w * (x - k)
  end do
else
  w = 1
  do k = 0, -n - 1
    y = y * (x + k)
  end do
end if
gamma_res = w / y
return
end function gamma

END MODULE eq_regions_mod
