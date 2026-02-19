! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SSALT &
  &( YDEAERATM,YDEAERSRC,KIDIA, KFDIA, KLON, KBINSS, &
  &  PCI  , PLSM , PCLK, PWIND, &
  &  PFLXSS & 
  &)

!*** * AER_SSALT* - SOURCE TERMS FOR SEA SALT AEROSOLS

!**   INTERFACE.
!     ----------
!          *AER_SSALT* IS CALLED FROM *AER_SRC*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        FROM O.BOUCHER's seasalt 

!     SOURCE.
!     -------
!     GONG et al., 1997

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2004-05-11
!        Samuel Rémy, 2016-06-17 : adjust to dry sea-salt emissions (was 80% RH
!        SS emissions)

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEAERATM, ONLY : TEAERATM
USE YOEAERSRC ,ONLY : TEAERSRC

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERATM)    ,INTENT(INOUT) :: YDEAERATM
TYPE(TEAERSRC)    ,INTENT(INOUT) :: YDEAERSRC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBINSS

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLON), PLSM(KLON), PCLK(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWIND(KLON) 

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFLXSS(KLON,KBINSS)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JBIN
REAL(KIND=JPRB) :: ZFROC, ZOCEA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SSALT',0,ZHOOK_HANDLE)
ASSOCIATE(RSSFLX=>YDEAERSRC%RSSFLX, &
 & RSS_RH80_MASSFAC=>YDEAERATM%RSS_RH80_MASSFAC)

! N.B.: RSSFLX in mg m-2 s-1 so PFLXSS is in g cm-2 s-1 after x by 1.E-3 x 1.E-4
! N.B.: RSSFLX in mg m-2 s-1 so PFLXSS is in kg m-2 s-1 after x by         1.E-6

DO JBIN=1,KBINSS
  DO JL=KIDIA,KFDIA
    PFLXSS(JL,JBIN)= 0._JPRB
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
!-- flux is considered only over full ocean grids, and for their open ocean fraction
    IF (ZOCEA == 1._JPRB) THEN
      PFLXSS(JL,JBIN) = ZFROC * RSSFLX(JBIN) * (PWIND(JL)**3.41_JPRB) &
      & * 1.E-06_JPRB / RSS_RH80_MASSFAC
    ENDIF
  ENDDO
ENDDO
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SSALT',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SSALT

