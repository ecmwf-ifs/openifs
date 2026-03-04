! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE BASCOE_J_TABLES_MODULE

USE PARKIND1, ONLY : JPIM, JPRB
USE BASCOE_J_MODULE, ONLY: NDISS
USE BASCOE_J_EXT_MODULE , ONLY : NABSLAYER

IMPLICIT NONE

SAVE

!**   DESCRIPTION 
!     ----------
!
!   Module defining strato photolysis rates.
!
! **** WARNING: THIS IS FOR THE OFFLINE MODE !
!
!
!     MODIFICATIONS.
!     --------------
!     2017-12-07: YVES CHRISTOPHE (YC) *BIRA*
!       extracted from bascoe_module.F90
!----------------

! Photolysis lookup table
  INTEGER(KIND=JPIM), parameter :: nsza     = 21   &! Nb of zenith anles in the J tables
     &                           , njo3     =  5   &! Nb of o3 profiles used for the tables
     &                           , njlev    = NABSLAYER+1    ! Nb of alt levels in the J tables

  REAL(KIND=JPRB), PARAMETER     :: o3du_target(njo3) = (/ 50., 250., 350., 450., 650. /)
  REAL(KIND=JPRB), PARAMETER     :: sza_tables(nsza) = &
   &                    (/  0. , 30. , 40. , 50. , 60. , 70. , 74. , 78.  , &
   &                       80. , 82. , 84. , 86. , 88. , 89.5, 91.1, 92.5 , &
   &                       94.1, 95.5, 97.1, 98.5, 100.                   /)

  REAL(KIND=JPRB), DIMENSION(njlev)      :: zjlev_tables
  REAL(KIND=JPRB), DIMENSION(njlev,njo3) :: o3col_tables
  REAL(KIND=JPRB), DIMENSION(nsza,njlev,njo3,ndiss) :: alj_tables ! LOG of J tables
    


END MODULE BASCOE_J_TABLES_MODULE
