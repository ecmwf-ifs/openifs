! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPSRH(KST,KEND,KPROMA,KLEV,&
                     & PU,PV,PDELP,PPFULL,POROG,PGEOPF,PSTRMMU,PSTRMMV,PSRH)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE YOMFPC   , ONLY : TNAMFPSCI


!******** FPSRH  ************

!      PURPOSE:
!      --------
!      Computes storm relative helicity

!      INTERFACE:
!      ----------     
!      *CALL FPSRH*
      

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
!        PSTRMMU : u component of storm motion vector
!        PSTRMMV : v component of storm motion vector

!         OUTPUT:
!        PSRH    : storm relative helicity        
!

!        IMPLICIT ARGUMENTS:
!        -------------------
!           NONE

!      METHOD:
!      -------
!        Computes vertical integral of dot product between storm relative wind and 3d vorticity
!        horizontal derivatives of w are neglected at the moment (in the voriticity part)

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
REAL(KIND=JPRB),INTENT(IN):: PSTRMMU(KPROMA), PSTRMMV(KPROMA)
REAL(KIND=JPRB),INTENT(OUT):: PSRH(KPROMA)

!local declarations
INTEGER (KIND=JPIM) :: JLON, JLEV, ITLEV(KPROMA)
REAL(KIND=JPHOOK)   :: ZHOOK_HANDLE
REAL(KIND=JPRB)     :: ZOVERDP , ZONEOVERDZ, ZDUDZ, ZDVDZ

REAL(KIND=JPRB), POINTER :: RSRHH

TYPE(TNAMFPSCI), TARGET :: YLNAMFPSCI


IF (LHOOK) CALL DR_HOOK('FPSRH',0,ZHOOK_HANDLE)
! ------------------------------------------------

RSRHH=>YLNAMFPSCI%RSRHH

PSRH(:)=0.0_JPRB
ITLEV(:)=KLEV

DO JLEV=KLEV,2,-1
  DO JLON=KST, KEND

    ! first compute vertical derivatives of u and v
    ZONEOVERDZ=(PGEOPF(JLON,JLEV-1)-PGEOPF(JLON,JLEV))/RG
    IF (ABS(ZONEOVERDZ).GT.0.0) THEN
      ZDUDZ=(PU(JLON,JLEV-1)-PU(JLON,JLEV))/ZONEOVERDZ
      ZDVDZ=(PV(JLON,JLEV-1)-PV(JLON,JLEV))/ZONEOVERDZ
    ELSE
      ZDUDZ=0.0_JPRB
      ZDVDZ=0.0_JPRB
    ENDIF
 
    PSRH(JLON)=PSRH(JLON) + ( (PV(JLON,JLEV)-PSTRMMV(JLON)) *ZDUDZ - & 
   &                          (PU(JLON,JLEV)-PSTRMMU(JLON)) *ZDVDZ ) & 
   &        * ZONEOVERDZ *MAX(0.0_JPRB,SIGN(1.0_JPRB,RSRHH*RG+POROG(JLON)-PGEOPF(JLON,JLEV)))
  ENDDO
ENDDO
  

IF (LHOOK) CALL DR_HOOK('FPSRH',1,ZHOOK_HANDLE)


END SUBROUTINE FPSRH
