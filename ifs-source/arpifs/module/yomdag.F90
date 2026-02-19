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

MODULE YOMDAG

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! CD YOMDAG : CONTIENT LA DATE, LE RESEAU, L'ECHEANCE ET LE TYPE DU GUESS.
! ---------

! AUTEUR    : PB   LE 25/10/90.      MODIFICATION :  PB   LE 12/11/90
! ------                             ------------

!   NDATGU : DATE DU RUN QUI A FOURNI LE GUESS ( AAAAMMDD )

!   NRESGU : RESEAU DU RUN QUI A FOURNI LE GUESS ( MMNN )

!   NECHGU : ECHEANCE DU GUESS EN HEURES (0 A 30).

!   NTYPGU : TYPE DU GUESS FOURNI AU MODELE:
!            2:CLIM / 3:ARPEGE / 4:EMERAUDE / 5:CEP

!  REMARQUES: DANS LE CAS D'UN GUESS CLIM (NTYPGU = 2):
!  ---------       NECHGU = 0
!                  NRESGU = 0
!                  NDATGU = MM00

!  INITIALISATION: LES VARIABLES SONT INITIALISEES DANS LE S/P SUSPEC
!  --------------  A PARTIR DE LA DATE ET DE L'IDENTIFICATEUR DES
!                  FICHIERS D'EBAUCHE.

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NDATGU
INTEGER(KIND=JPIM) :: NRESGU
INTEGER(KIND=JPIM) :: NECHGU
INTEGER(KIND=JPIM) :: NTYPGU

!     ------------------------------------------------------------------
END MODULE YOMDAG
