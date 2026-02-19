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

SUBROUTINE XUTRII(KS,PS,KND)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  FONCTION REALISEE
!  -----------------
!  CLASSEMENT DES ADRESSES DES ELEMENTS D'UN TABLEAU DE REELS ,
!  CORRESPONDANT AU CLASSEMENT PAR ORDRE CROISSANT DES ELEMENTS DE
!  CE TABLEAU .

!  VARIABLES D'ENTREE (ARGUMENTS)
!  ------------------------------
!  INTEGER    KS    : NOMBRE D'ENTIERS A CLASSER
!  REEEL      PS(KS): CONTIENT LES ENTIERS A CLASSER

!  VARIABLES DE SORTIE (ARGUMENTS)
!  -------------------------------
!  INTEGER    KND(KS) : CONTIENT LES ADRESSES DES ELEMENTS DE IS(LS) ,
!                       DANS L'ORDRE CORRESPONDANT AU CLASSEMENT PAR
!                       ORDRE CROISSANT DES ELEMENTS DE IS(LS)

!  REFERENCES
!  ----------
!  NOM DE L'AUTEUR DU SOUS-PROGRAMME :  X.VILTART
!  DATE                              : SEPTEMBRE 85

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  Purpose :   Ordering of indices of a real array with increasing values
!  --------    of elements

!  Arguments : KS     : number of indices to be ordered
!  ----------  PS(KS) : real array to be ordered
!              KND(KS): integer array containing ordered indices

!  References :
!  ----------
!    From XUTRII written by X.VILTART (SEPTEMBRE 85)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KS 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PS(KS) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KND(KS) 
LOGICAL ::       LLINV

INTEGER(KIND=JPIM) :: I, IS1, J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('XUTRII',0,ZHOOK_HANDLE)
DO J=1,KS
  KND(J) = J
ENDDO
IS1 = KS - 1

LLINV = .TRUE.
DO WHILE (LLINV)
  LLINV = .FALSE.
  DO J=1,IS1
    IF (PS(KND(J)) > PS(KND(J+1))) THEN
      I = KND(J+1)
      KND(J+1) = KND(J)
      KND(J) = I
      LLINV = .TRUE.
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('XUTRII',1,ZHOOK_HANDLE)
END SUBROUTINE XUTRII

