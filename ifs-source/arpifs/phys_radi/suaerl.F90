! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUAERL(PRAER)

!**** *SUAERL*   - INITIALIZE COMMON YOEAER

!     PURPOSE.
!     --------
!           INITIALIZE YOEAER, THE COMMON THAT CONTAINS THE
!           RADIATIVE CHARACTERISTICS OF THE AEROSOLS

!**   INTERFACE.
!     ----------
!              -----        -----

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOEAER

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "IFS MODEL"

!     AUTHOR.
!     -------
!      JEAN-JACQUES MORCRETTE *ECMWF*
!      ORIGINAL : 88-02-15

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!      ----------------------------------------------------------------

!*       1.    SHORTWAVE COEFFICIENTS
!              ----------------------

!      Initialised in SUAERSn depending on NSW number of SW spectral 
!      intervals

!      ----------------------------------------------------------------

!*       2.    LONGWAVE COEFFICIENTS
!              ---------------------

!=======================================================================
!-- The (old) five aerosol types were respectively:

!  1/ continental average (+desert)       2/ maritime
!  3/ urban                               4/ volcanic active
!  5/ stratospheric background

!      RAER  = RESHAPE((/
!     &          .038520, .037196, .040532, .054934, .038520
!     &        , .126130, .18313 , .10357 , .064106, .126130
!     &        , .012579, .013649, .018652, .025181, .012579
!     &        , .011890, .016142, .021105, .028908, .011890
!     &        , .013792, .026810, .052203, .066338, .013792 /)
!     & ,SHAPE=(/5,5/))

!=======================================================================

!-- The six aerosol types are respectively:

!  1/ continental average                 2/ maritime
!  3/ desert                              4/ urban
!  5/ volcanic active                     6/ stratospheric background

! The quantities given are:
! TAU : ratio of average optical thickness in interval to that at 0.55 
!       micron
! PIZA: average single scattering albedo
! CGA : average asymmetry factor

! computed from Hess and Koepke (con, mar, des, urb)
!          from Bonnel et al.   (vol, str)

!-- data are entered for the 6 spectral intervals of the LW scheme (line)
!   and the different types of aerosols (column)

IMPLICIT NONE
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRAER(6,6) ! LW absorption coefficients
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SUAERL',0,ZHOOK_HANDLE)
PRAER( :, 1)= (/&
 & .036271_JPRB, .030153_JPRB, .017343_JPRB, .015002_JPRB, .008806_JPRB, .006865_JPRB /)  
PRAER( :, 2)= (/&
 & .026561_JPRB, .032657_JPRB, .017977_JPRB, .014210_JPRB, .016775_JPRB, .022123_JPRB /)  
PRAER( :, 3)= (/&
 & .014897_JPRB, .016359_JPRB, .019789_JPRB, .030777_JPRB, .013341_JPRB, .014321_JPRB /)  
PRAER( :, 4)= (/&
 & .001863_JPRB, .002816_JPRB, .002355_JPRB, .002557_JPRB, .001774_JPRB, .001780_JPRB /)  
PRAER( :, 5)= (/&
 & .011890_JPRB, .016142_JPRB, .021105_JPRB, .028908_JPRB, .011890_JPRB, .011890_JPRB /)  
PRAER( :, 6)= (/&
 & .013792_JPRB, .026810_JPRB, .052203_JPRB, .066338_JPRB, .013792_JPRB, .013792_JPRB /)  

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUAERL',1,ZHOOK_HANDLE)
END SUBROUTINE SUAERL
