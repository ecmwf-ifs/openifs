! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SSALT_GRYTHE &
  &( YDEAERATM,YDEAERSNK,KIDIA, KFDIA, KLON, &
  &  PCI  , PLSM , PCLK , PWIND, PSST, PRH, &
  &  PFLXSS)

!*** * AER_SSALT_GRYTHE* - SOURCE TERMS FOR SEA SALT AEROSOLS

!**   INTERFACE.
!     ----------
!          *AER_SSALT_GRYTHE* IS CALLED FROM *AER_SRC*.

!     AUTHOR.
!     -------
!        Original version
!        P. Nabat May 2015 % S. Rémy Nov 2016

!     SOURCE.
!     -------
!     Implementation of the parameterization of Grythe et al. (2014,
!     page 1286 eq. 7), in which sea salt emission fluxes depend both on
!     surface wind speed, and possibly SST.

!     MODIFICATIONS.
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST , ONLY: RPI
USE YOEAERSNK ,ONLY : TEAERSNK
USE YOEAERATM, ONLY : TEAERATM
IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERATM)    ,INTENT(INOUT) :: YDEAERATM
TYPE(TEAERSNK)    ,INTENT(INOUT) :: YDEAERSNK
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

REAL(KIND=JPRB),PARAMETER::ZRHOSS=2160._JPRB
REAL(KIND=JPRB)::ZRAYMD(3)
REAL(KIND=JPRB)::ZRAYM,ZEMSALT,ZTSC,ZSALTMAS
REAL(KIND=JPRB)::ZLIM1,ZLIM2,ZLIM3,ZLIM4       !Bin limits (dry diameter)
REAL(KIND=JPRB)::ZCOEFSST(KLON)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZFROC, ZOCEA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SSALT_GRYTHE',0,ZHOOK_HANDLE)
ASSOCIATE(RRHTAB=>YDEAERSNK%RRHTAB, &
  & RSSGROWTH_RHTAB=>YDEAERSNK%RSSGROWTH_RHTAB, &
  & RSS_DRY_MASSFAC=>YDEAERATM%RSS_DRY_MASSFAC, &
  & RSSDENS_RHTAB=>YDEAERSNK%RSSDENS_RHTAB)
  PFLXSS = 0._JPRB
  ZRAYMD = (/0.3_JPRB,2.0_JPRB,7.0_JPRB/) ! Valeurs médianes classiques

     DO JL=KIDIA,KFDIA
       ZTSC = PSST(JL)-273.15_JPRB
       ZCOEFSST(JL) = 0.3_JPRB + 0.1_JPRB*ZTSC -0.0076_JPRB*ZTSC**2._JPRB + 0.00021_JPRB*ZTSC**3._JPRB
     ENDDO
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
     ZLIM1=0.03_JPRB
     ZLIM2=0.5_JPRB
     ZLIM3=5.0_JPRB
     ZLIM4=20.0_JPRB

     ZSALTMAS=0._JPRB
     ISIZE=1
        ZRAYM=ZRAYMD(ISIZE)
        ZEMSALT = 235._JPRB*PWIND(JL)**3.5_JPRB*EXP(-0.55_JPRB*(LOG(ZRAYM/0.1_JPRB))**2._JPRB)&
           &+0.2_JPRB*PWIND(JL)**3.5_JPRB*EXP(-1.5_JPRB*(LOG(ZRAYM/3._JPRB))**2._JPRB)&
           &+6.8_JPRB*PWIND(JL)**3._JPRB*EXP(-((LOG(ZRAYM/30._JPRB))**2._JPRB))
        ZSALTMAS= ZSALTMAS +&
            &4._JPRB/3._JPRB*RPI*ZEMSALT*ZRHOSS*(ZRAYM/2._JPRB*1.0E-6_JPRB)**3._JPRB*(ZLIM2-ZLIM1)
     PFLXSS(JL,1) = ZFROC*ZSALTMAS*ZCOEFSST(JL)/RSS_DRY_MASSFAC

     ZSALTMAS=0._JPRB
     ISIZE=2
        ZRAYM=ZRAYMD(ISIZE)
        ZEMSALT = 235._JPRB*PWIND(JL)**3.5_JPRB*EXP(-0.55_JPRB*(LOG(ZRAYM/0.1_JPRB))**2._JPRB)&
           &+0.2_JPRB*PWIND(JL)**3.5_JPRB*EXP(-1.5_JPRB*(LOG(ZRAYM/3._JPRB))**2._JPRB)&
           &+6.8_JPRB*PWIND(JL)**3._JPRB*EXP(-((LOG(ZRAYM/30._JPRB))**2._JPRB))
        ZSALTMAS= ZSALTMAS +&
            &4._JPRB/3._JPRB*RPI*ZEMSALT*ZRHOSS*(ZRAYM/2._JPRB*1.0E-6_JPRB)**3._JPRB*(ZLIM3-ZLIM2)
     PFLXSS(JL,2) = ZFROC*ZSALTMAS*ZCOEFSST(JL)/RSS_DRY_MASSFAC

     ZSALTMAS=0._JPRB
     ISIZE=3
        ZRAYM=ZRAYMD(ISIZE)
        ZEMSALT = 235._JPRB*PWIND(JL)**3.5_JPRB*EXP(-0.55_JPRB*(LOG(ZRAYM/0.1_JPRB))**2._JPRB)&
           &+0.2_JPRB*PWIND(JL)**3.5_JPRB*EXP(-1.5_JPRB*(LOG(ZRAYM/3._JPRB))**2._JPRB)&
           &+6.8_JPRB*PWIND(JL)**3._JPRB*EXP(-((LOG(ZRAYM/30._JPRB))**2._JPRB))
        ZSALTMAS= ZSALTMAS +&
            &4._JPRB/3._JPRB*RPI*ZEMSALT*ZRHOSS*(ZRAYM/2._JPRB*1.0E-6_JPRB)**3._JPRB*(ZLIM4-ZLIM3)
     PFLXSS(JL,3) = ZFROC*ZSALTMAS*ZCOEFSST(JL)/RSS_DRY_MASSFAC
    ENDIF
   ENDDO


!-----------------------------------------------------------------------Z
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SSALT_GRYTHE',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SSALT_GRYTHE

