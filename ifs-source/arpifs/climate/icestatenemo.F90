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

SUBROUTINE ICESTATENEMO(YDMCC,KSTGLO,KIDIA,KFDIA,PTHKICE,PSNTICE,PALBICE,PTEMPSEA,PTEMPICE)
!
!**** *ICESTATENEMO*  - Retrieve the thickness of the sea ice and the snow layer on top.
!
!     Purpose.
!     --------
!       Retrieve the thickness of the sea ice and the snow layer on top.
!       This data is used to compute the conductivity of the skin layer in the surcace energy balance model.
!
!**   Interface.
!     ----------
!       *CALL*  *ICESTATENEMO(KSTCLO,KIDIA,KFDIA,PTHKICE,PSNTICE)
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!     Method:
!     ------
!       
!     Externals:
!     ---------
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!       C. A. Severijns, KNMI, 7-Jul-2008 (for ECEARTH-R1)
!
!     Modifications.
!     --------------
!     Linus Magnusson, Added to IFS
!     F. Vana  05-Mar-2015  Support for single precision
!
!     -----------------------------------------------------------
USE PARKIND1 , ONLY : JPIM,  JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMCC   , ONLY : TMCC

IMPLICIT NONE

TYPE(TMCC)         ,INTENT(IN) :: YDMCC
INTEGER(KIND=JPIM), INTENT(IN) :: KSTGLO, KIDIA, KFDIA
REAL(KIND=JPRB),INTENT(INOUT), DIMENSION(KIDIA:KFDIA), OPTIONAL :: PTHKICE
REAL(KIND=JPRB),INTENT(INOUT), DIMENSION(KIDIA:KFDIA), OPTIONAL :: PSNTICE
REAL(KIND=JPRB),INTENT(INOUT), DIMENSION(KIDIA:KFDIA), OPTIONAL :: PALBICE
REAL(KIND=JPRB),INTENT(INOUT), DIMENSION(KIDIA:KFDIA), OPTIONAL :: PTEMPSEA
REAL(KIND=JPRB),INTENT(INOUT), DIMENSION(KIDIA:KFDIA), OPTIONAL :: PTEMPICE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IST,IEND

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ICESTATENEMO',0,ZHOOK_HANDLE)
ASSOCIATE(LNEMOCOUP=>YDMCC%LNEMOCOUP,CPLNG_FLD=>YDMCC%CPLNG_FLD)

IF (LNEMOCOUP) THEN

   IST  = KSTGLO - 1 + KIDIA
   IEND = KSTGLO - 1 + KFDIA

   IF (PRESENT(PTHKICE)) PTHKICE(:) = &
      & REAL(CPLNG_FLD(YDMCC%IP_A_ICE_THICKNESS)%D(IST:IEND,1,1),JPRB)
   IF (PRESENT(PSNTICE)) PSNTICE(:) = &
      & REAL(CPLNG_FLD(YDMCC%IP_A_SNOW_THICKNESS)%D(IST:IEND,1,1),JPRB)
   IF (PRESENT(PALBICE)) PALBICE(:) = &
      & REAL(CPLNG_FLD(YDMCC%IP_A_ICE_ALBEDO)%D(IST:IEND,1,1),JPRB)
   IF (PRESENT(PTEMPSEA)) PTEMPSEA(:) = &
      & REAL(CPLNG_FLD(YDMCC%IP_A_SST)%D(IST:IEND,1,1),JPRB)
   IF (PRESENT(PTEMPICE)) PTEMPICE(:) = &
      & REAL(CPLNG_FLD(YDMCC%IP_A_ICE_TEMP)%D(IST:IEND,1,1),JPRB)

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ICESTATENEMO',1,ZHOOK_HANDLE)

END SUBROUTINE ICESTATENEMO
