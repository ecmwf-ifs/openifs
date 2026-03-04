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

MODULE YOMDPRECIPS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE TDPRECIPS
!*
!     ------------------------------------------------------------------
!     CONSTANTES POUR LA DETERMINATION DU TYPE DE PRECIPITATIONS
!     
!       HDPRECIPS   : HAUTEUR DE LA DETERMINATION ( EN GENERAL 2M ).
!                   : HEIGHT ( IN GENERAL 2M ).
!       HDCLWC      : HAUTEUR DU SEUIL DE CLWC.
!                   : HEIGHT OF CLWC THRESHOLD.
!       RDHAIL1     : PREMIER SEUIL POUR LA GRELE / HAIL 
!                   : FIRST THRESHOLD 
!       RDHAIL2     : DEUXIEME SEUIL POUR LA GRELE / HAIL 
!                   : SECOND THRESHOLD 
!       RDSEUIL1    : PREMIER SEUIL POUR LE RATIO PLUIE/PRECIPS
!                   : FIRST THRESHOLD 
!       RDSEUIL2    : DEUXIEME SEUIL POUR LE RATIO PLUIE/PRECIPS
!                   : SECOND THRESHOLD 
!       RDSEUIL3    : TROISIEME SEUIL POUR LE RATIO PLUIE/PRECIPS
!                   : THIRD THRESHOLD 
!       RDSEUIL4    : QUATRIEME SEUIL POUR LE RATIO PLUIE/PRECIPS
!                   : FOURTH THRESHOLD 
!       RDSEUIL5    : CINQUIEME SEUIL POUR LE RATIO PLUIE/PRECIPS
!                   : FIFTH THRESHOLD 
!       RDCLWC      : SEUIL POUR LE CLWC
!                   : THRESHOLD OF CLWC 
!       RPRECSEUIL  : SEUIL DE PRESENCE DE PRECIPS
!                   : THRESHOLD  OF PRECIPITATION
!       RHTOP       : SEUIL DE HAUTEUR AU DESSUS DU SOL
!                   : H THRESHOLD ABOGE GROUND
!       RTPW        : T'W THRESHOLD
!       RAWARM      : AREA WITH T > 0 IN VERTICAL PROFILES USED FOR FREEZING
!                     RAIN 
!       RACOLD      : AREA WITH T < 0 IN VERTICAL PROFILES USED FOR FREEZING
!                     RAIN 
!       NDTPREC     : NB OF TIME STEPS USED FOR DPRECIPS
!       NDTPRECCUR  : CURRENT TIME STEP NUMBER


REAL(KIND=JPRB) :: HDPRECIPS
REAL(KIND=JPRB) :: HDCLWC
REAL(KIND=JPRB) :: RDHAIL1
REAL(KIND=JPRB) :: RDHAIL2
REAL(KIND=JPRB) :: RDSEUIL1
REAL(KIND=JPRB) :: RDSEUIL2
REAL(KIND=JPRB) :: RDSEUIL3
REAL(KIND=JPRB) :: RDSEUIL4
REAL(KIND=JPRB) :: RDSEUIL5
REAL(KIND=JPRB) :: RDCLWC
REAL(KIND=JPRB) :: RPRECSEUIL
REAL(KIND=JPRB) :: RHTOP
REAL(KIND=JPRB) :: RTPW
REAL(KIND=JPRB) :: RAWARM
REAL(KIND=JPRB) :: RACOLD
INTEGER(KIND=JPIM) :: NDTPREC  = 0
INTEGER(KIND=JPIM) :: NDTPREC2 = 0
INTEGER(KIND=JPIM) :: NDTPRECCUR,NDTPRECCUR2

END TYPE TDPRECIPS

!     ------------------------------------------------------------------
END MODULE YOMDPRECIPS
