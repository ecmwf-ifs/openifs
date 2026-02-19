! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SSALT_ALBERT &
  &( YDEAERATM,YDEAERSRC,KIDIA, KFDIA, KLON, &
  &  PCI  , PLSM , PCLK , PWIND, PSST, PRH, &
  &  PFLXSS)

!*** * AER_SSALT_ALBERT* - SOURCE TERMS FOR SEA SALT AEROSOLS

!**   INTERFACE.
!     ----------
!          *AER_SSALT_ALBERT* IS CALLED FROM *AER_SRC*.

!     AUTHOR.
!     -------
!        Original version
!        S. Rémy May 2019

!     SOURCE.
!     -------
!     Implementation of the parameterization of Albert et al. (2016)

!     MODIFICATIONS.
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOEAERSRC ,ONLY : TEAERSRC
USE YOEAERATM, ONLY : TEAERATM
IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERATM)    ,INTENT(INOUT) :: YDEAERATM
TYPE(TEAERSRC)    ,INTENT(INOUT) :: YDEAERSRC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLON), PLSM(KLON), PCLK(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWIND(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRH(KLON)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFLXSS(KLON,3)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM)::ISIZE

REAL(KIND=JPRB)::ZEMSALT,ZTSC,ZSALTMAS
REAL(KIND=JPRB)::ZA, ZB, ZR

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZFROC, ZOCEA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SSALT_ALBERT',0,ZHOOK_HANDLE)
ASSOCIATE(RSS_RH80_MASSFAC=>YDEAERATM%RSS_RH80_MASSFAC, &
  & RSSFLX=>YDEAERSRC%RSSFLX )
  PFLXSS(:,:) = 0._JPRB

  DO JL=KIDIA,KFDIA
    ZSALTMAS=0._JPRB

    ! LSM treats lakes as ocean, but (with a few exceptions e.g. the Dead Sea) they don't produce sea salt.
    ! This logic should work for either a fractional LSM or a binary one
    IF (PLSM(JL) == 0._JPRB) THEN
      ! Even with a fractional LSM, this must be all ocean or all lake, because there would be no land to divide them.
      ! Thresholding ensures sensible behaviour if LSM is binary but lake cover fractional.
      IF (PCLK(JL) >= 0.5_JPRB) THEN
        ZOCEA= 0._JPRB
      ELSE
        ZOCEA= 1._JPRB
      ENDIF
    ELSE
      ZOCEA= MAX(1._JPRB-PLSM(JL)-PCLK(JL), 0._JPRB)
    ENDIF
    ! Also exclude sea ice
    ZFROC= ZOCEA*(1._JPRB-PCI(JL))
    !-- flux is considered only over full ocean grids, and for their open ocean
    !  fraction
   IF (ZOCEA == 1._JPRB) THEN
     ! dry diameter = wet radius

     ZTSC = MAX(PSST(JL)-273.15_JPRB,0._JPRB)
     ZA=8.46E-5_JPRB + 1.63E-6_JPRB*ZTSC - 3.35E-8_JPRB*ZTSC**2._JPRB
     ZB=3.354_JPRB - 6.2E-2_JPRB * ZTSC
 
     ISIZE=1
     ZR =  RSSFLX(ISIZE) * 0.26_JPRB
     ZEMSALT =  ZA * ((PWIND(JL)+ZB)**2._JPRB)* ZR
     PFLXSS(JL,1) = ZFROC*ZEMSALT/RSS_RH80_MASSFAC
!     IF ( PFLXSS(JL,1) > 1.E-13_JPRB) THEN
!     write (*,*) "SSALT_ALBERT,",JL,ZTSC,PSST(JL),PWIND(JL),ZA,ZB,ZA*((PWIND(JL)+ZB)**2._JPRB),ZEMSALT,ZR
!     ENDIF

     ISIZE=2
     ZR =  RSSFLX(ISIZE) * 0.26_JPRB
     ZEMSALT =  ZA * ((PWIND(JL) +ZB)**2._JPRB)* ZR
     PFLXSS(JL,2) = ZFROC*ZEMSALT/RSS_RH80_MASSFAC

     ISIZE=3
     ZR =  RSSFLX(ISIZE) * 0.26_JPRB
     ZEMSALT =  ZA * ((PWIND(JL) +ZB)**2._JPRB)* ZR
     PFLXSS(JL,3) = ZFROC*ZEMSALT/RSS_RH80_MASSFAC
    ENDIF
   ENDDO


!-----------------------------------------------------------------------Z
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SSALT_ALBERT',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SSALT_ALBERT

