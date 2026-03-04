FUNCTION ZMODULO(PX,PY)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *zmodulo* fonction modulo périodique sur l'ensemble des reels.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   96-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8) :: ZMODULO
REAL(KIND=8) :: PX
REAL(KIND=8) :: PY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ZMODULO',0,ZHOOK_HANDLE)
ZMODULO=MOD(PX,PY)
IF(ZMODULO < 0) ZMODULO=ZMODULO+PY
IF (LHOOK) CALL DR_HOOK('ZMODULO',1,ZHOOK_HANDLE)
END FUNCTION ZMODULO
FUNCTION ZMODYX(PX,PY,PZ)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *zmodyx* amène px dans l'intervalle [py,pz[ modulo (pz-py).
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   96-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! i.e. zmodyx sera le réel tel que
! py <= zmodyx < pz, et px-zmodyx=k*(pz-py.)
! --------------------------------------------------------------
IMPLICIT NONE
REAL(KIND=8) :: ZMODYX
REAL(KIND=8) :: ZMODULO
REAL(KIND=8) :: PX
REAL(KIND=8) :: PY
REAL(KIND=8) :: PZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ZMODYX',0,ZHOOK_HANDLE)
ZMODYX=PY+ZMODULO(PX-PY,PZ-PY)
IF (LHOOK) CALL DR_HOOK('ZMODYX',1,ZHOOK_HANDLE)
END FUNCTION ZMODYX
FUNCTION IMODULO(KX,KY)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *imodulo* fonction modulo périodique sur l'ensemble des entiers.
! --------------------------------------------------------------
! Sujet:
!       Cette fonction pallie le manque, en f77, 
!       d'une fonction modulo périodique sur l'ensemble des entiers.
!       En f90, ce problème a été résolu par la fonction modulo.
!       Il suffit alors d'appeler modulo et non imodulo.
! Limitations: le deuxième argument ky est supposé > 0.
!       On n'a pas fait ici l'analyse de l'adéquation imodulo <> modulo
!       pour les ky <= 0, pour lesquels la fonction mathématique
!       n'est pas usuellement définie. On ne traite donc l'adéquation 
!       que pour tous les kx relatifs, à ky > 0.
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   2001-01, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
IMPLICIT NONE
INTEGER(KIND=4) :: IMODULO
INTEGER(KIND=4) :: KX
INTEGER(KIND=4) :: KY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('IMODULO',0,ZHOOK_HANDLE)
IMODULO=MOD(KX,KY)
IF(IMODULO < 0) IMODULO=IMODULO+KY
IF (LHOOK) CALL DR_HOOK('IMODULO',1,ZHOOK_HANDLE)
END FUNCTION IMODULO
