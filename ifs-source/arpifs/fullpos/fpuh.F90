! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPUH(KST,KEND,KPROMA,KLEV,PUH_LOWER_LIMIT,PUH_UPPER_LIMIT,&
                     & PGEOPH,PWVEL0,PVOR,PUH)
                     
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG


!******** FPUH  ************

!      PURPOSE:
!      --------
!      Compute updraft helicity

!      INTERFACE:
!      ----------     
!      *CALL FPUH*
      

!        EXPLICIT ARGUMENTS:
!        -------------------
!          INPUT:
!        KST     : start of work
!        KEND    : end of work
!        KPROMA  : dimension of work
!        KLEV    : number of levels
!        PGEOPH  : helf levels geopotential
!        PWVEL0  : w component of wind on bottom half levels
!        PVOR    : vorticity on full levels


!         OUTPUT:
!        PUH    : updraft helicity        
!

!        IMPLICIT ARGUMENTS:
!        -------------------
!           NONE

!      METHOD:
!      -------
!       

!      AUTHOR:
!      -------
!        Clemens Wastl *ZAMG*

!      MODIFICATIONS:
!      --------------

!        02-2022  C. Wittmann   extract from fpsrh routine  
!


IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KST,KEND,KPROMA,KLEV
REAL(KIND=JPRB),INTENT(IN) :: PUH_LOWER_LIMIT
REAL(KIND=JPRB),INTENT(IN) :: PUH_UPPER_LIMIT
REAL(KIND=JPRB),INTENT(IN) :: PGEOPH(KPROMA,0:KLEV)
REAL(KIND=JPRB),INTENT(IN) :: PWVEL0(KPROMA,KLEV),PVOR(KPROMA,KLEV)
REAL(KIND=JPRB),INTENT(OUT):: PUH(KPROMA)

!local declarations 
INTEGER (KIND=JPIM) :: JLON, JLEV
REAL(KIND=JPHOOK)   :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZONEOVERDZ,ZUHLL,ZUHUL, Z1SRG
REAL(KIND=JPRB) :: ZTEST_LL, ZTEST_UL


IF (LHOOK) CALL DR_HOOK('FPUH',0,ZHOOK_HANDLE)
! ------------------------------------------------

Z1SRG=1._JPRB/RG

PUH(KST:KEND)=0._JPRB

DO JLEV=KLEV,1,-1
  DO JLON=KST, KEND

    ZONEOVERDZ=(PGEOPH(JLON,JLEV-1)-PGEOPH(JLON,JLEV))*Z1SRG
    
    ZUHLL=PGEOPH(JLON,KLEV)*Z1SRG+PUH_LOWER_LIMIT
    ZUHUL=PGEOPH(JLON,KLEV)*Z1SRG+PUH_UPPER_LIMIT

    ZTEST_LL=MIN(MAX(SIGN(1._JPRB,PGEOPH(JLON,JLEV)/RG-ZUHLL),0._JPRB),1.0_JPRB)
    ZTEST_UL=MIN(MAX(SIGN(1._JPRB,ZUHUL-PGEOPH(JLON,JLEV)/RG),0._JPRB),1.0_JPRB)

    PUH(JLON)=PUH(JLON)+MAX(0.0_JPRB,PWVEL0(JLON,JLEV))*PVOR(JLON,JLEV)*ZONEOVERDZ*ZTEST_LL*ZTEST_UL*Z1SRG
    
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('FPUH',1,ZHOOK_HANDLE)

END SUBROUTINE FPUH
