! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPSTRMM(KST,KEND,KPROMA,KLEV,&
                     & PU,PV,PDELP,PPFULL,POROG,PGEOPF,PSTRMMU,PSTRMMV)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE YOMFPC   , ONLY : TNAMFPSCI


!******** FPSTRMM  ************

!      PURPOSE:
!      --------
!      Compute storm motion (wind average in lowest RSTRMMH meters)

!      INTERFACE:
!      ----------     
!      *CALL FPSTRMM*
      

!        EXPLICIT ARGUMENTS:
!        -------------------
!          INPUT:
!        KST     : start of work
!        KEND    : end of work
!        KPROMA  : dimension of work
!        KLEV    : number of levels
!        PU      : u wind
!        PV      : v wind
!        PDELP   : layer thickness (Pa) 
!        PPFULL  : pressure on full levels (Pa)


!          OUTPUT:
!        PSTRMMU : u component of storm motion vector
!        PSTRMMV : v component of storm motion vector

!        IMPLICIT ARGUMENTS:
!        -------------------
!           NONE

!      METHOD:
!      -------
!        Simply computes the average wind of all the layers below RSTRMMH meters

!      AUTHOR:
!      -------
!        Jure Cedilnik *ARSO*

!      MODIFICATIONS:
!      --------------

!        original version: August 2013
!



IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN) :: KST,KEND,KPROMA,KLEV
REAL(KIND=JPRB),INTENT(IN):: PU(KPROMA,KLEV),PV(KPROMA,KLEV),PDELP(KPROMA,KLEV),PPFULL(KPROMA,KLEV) 
REAL(KIND=JPRB),INTENT(IN):: POROG(KPROMA), PGEOPF(KPROMA,KLEV)
REAL(KIND=JPRB),INTENT(OUT):: PSTRMMU(KPROMA), PSTRMMV(KPROMA)

!local declarations
INTEGER (KIND=JPIM) :: JLON, JLEV, ITLEV(KPROMA)
REAL(KIND=JPRB) :: ZOVERDP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB), POINTER :: RSTRMMH

TYPE(TNAMFPSCI), TARGET :: YLNAMFPSCI

IF (LHOOK) CALL DR_HOOK('FPSTRMM',0,ZHOOK_HANDLE)
! ------------------------------------------------

RSTRMMH=>YLNAMFPSCI%RSTRMMH

PSTRMMU(:)=0.0_JPRB
PSTRMMV(:)=0.0_JPRB
ITLEV(:)=KLEV

DO JLEV=KLEV,1,-1
  DO JLON=KST, KEND
    PSTRMMU(JLON)=PSTRMMU(JLON)-PU(JLON,JLEV)*PDELP(JLON,JLEV)*MAX(0.0_JPRB,SIGN(1.0_JPRB,RSTRMMH*RG+POROG(JLON)-PGEOPF(JLON,JLEV)))
    PSTRMMV(JLON)=PSTRMMV(JLON)-PV(JLON,JLEV)*PDELP(JLON,JLEV)*MAX(0.0_JPRB,SIGN(1.0_JPRB,RSTRMMH*RG+POROG(JLON)-PGEOPF(JLON,JLEV)))
    ITLEV(JLON)=ITLEV(JLON)-MAX(0.0_JPRB,SIGN(1.0_JPRB,RSTRMMH*RG+POROG(JLON)-PGEOPF(JLON,JLEV)))
  ENDDO
ENDDO
  
DO JLON=KST,KEND
  ZOVERDP=(PPFULL(JLON,KLEV)-PPFULL(JLON,ITLEV(JLON)))
  IF (ABS(ZOVERDP).GT.0.0_JPRB) THEN 
    PSTRMMU(JLON)=PSTRMMU(JLON)*1.0_JPRB/ZOVERDP
    PSTRMMV(JLON)=PSTRMMV(JLON)*1.0_JPRB/ZOVERDP
  ELSE
    PSTRMMU(JLON)=0.0
    PSTRMMV(JLON)=0.0
  ENDIF

ENDDO


IF (LHOOK) CALL DR_HOOK('FPSTRMM',1,ZHOOK_HANDLE)


END SUBROUTINE FPSTRMM
