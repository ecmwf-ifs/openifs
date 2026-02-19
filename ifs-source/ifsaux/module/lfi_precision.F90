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

MODULE LFI_PRECISION

USE PARKIND1, ONLY : JPRB, JPIB, JPIM, JPRD, JPRM
IMPLICIT NONE

!
!----- DESCRIPTION DES "PARAMETER" DU LOGICIEL DE FICHIERS INDEXES -----
!
!     JPDBLE= PRECISION UTILISE POUR LES ENTIERS:
!     * SI JPDBLE=8 les INTEGER (KIND=JPDBLE) seront a priori en 64 BITS
!     * SI JPDBLE=4 les INTEGER (KIND=JPDBLE) seront a priori en 32 BITS
!
!     JPDBLR= PRECISION UTILISE POUR LES FLOTTANTS (REELS):
!     * SI JPDBLR=8 les REAL (KIND=JPDBLR) seront a priori en 64 BITS
!     * SI JPDBLR=4 les REAL (KIND=JPDBLR) seront a priori en 32 BITS
!
!     (les conventions peuvent dependre de la plate-forme consideree)
!
!     JP_SIMPLE_ENTIER= sous-type entier permettant de representer
!                       l'intervalle +/-(10**9 - 1)
!

INTEGER, PARAMETER :: JPDBLE = JPRB, JPDBLR = JPRB, JPDBLD = JPRD, JPDBLM = JPRM
INTEGER, PARAMETER :: JP_SIMPLE_ENTIER = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: JPLIKB = JPIB, JPLIKM = JPIM


END MODULE LFI_PRECISION
