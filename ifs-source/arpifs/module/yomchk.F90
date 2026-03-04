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

MODULE YOMCHK

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!**************  COMDECK YOMCHK  ***************************

!      CONTROLS FOR GRIDPOINT EVOLUTION DIAGNOSTICS

! These variables are set through namelist NEMCHK :
! LECHKEVO=.TRUE. : diagnostics required
!                   activation du processus
! LECHKTND=.TRUE. : global diagnostics for tendencies required
!                   diagnostics (globaux) sur les tendances
! NFRQCHK         : frequency of diagnostics (in s , default: every timestep)
!                   frequence des diagnostics, en secondes
! NGPCHK          : number of gridpoints required
!                   nombre de points de grille selectionnes
! NXCHK(i), NYCHK(i), i=1,NGPCHK :
!                   longitude and latitude indices of gridpoints
!                   indices de longitude et latitude, respectivement,
!                   des points choisis
!                   1 <= NYCHK(i) <= NDGLG    1 <= NXCHK(i) <= NLOEN(NYCHK(i))
! NFLDCHK         : number of "2d" spectral fields required
!                   nombre de champs spectraux "2d" selectionnes
!                   (champs a un niveau donne ou pression de surface ou ln(ps))

! NNFCHK(i), i=1,NFLDCHK :
!                   fields indices
!                   indices des champs choisis
! 3d :  =(n-1)*NFLEVG+k
!       k= vertical level /niveau vertical  (1 <= k <= NFLEVG)
!       n=1 -> vorticity / tourbillon    n=2 -> divergence
!       n=3 -> u                         n=4 -> v
!       n=4+p -> p-th thermodynamic field / p-ieme variable thermodynamique
!               (n=5 -> T  n=6 -> q)
!       n=4+NFTHER+p -> p-th GFL field / p-ieme variable GFL
!       =(4+NFTHER+ ???? )*NFLEVG+1 -> ps or ln(ps)
!       !!! topic about GFL to be updated, not yet coded currently,
!           this is probably the spectral GFL which must be taken in account.

! LECHKPS =.TRUE. : ps required instead of ln(ps)
!                   diagnostic sur ps et non ln(ps)

! At each time-step, NLENCHK values are stored.
! NCHKTEND is a dimensioning constant used for tendency diagnostics

INTEGER(KIND=JPIM), PARAMETER :: JPGPCHK=100
INTEGER(KIND=JPIM), PARAMETER :: JPFLDCHK=125
LOGICAL :: LECHKEVO
LOGICAL :: LECHKTND
LOGICAL :: LECHKPS
INTEGER(KIND=JPIM) :: NXCHK(JPGPCHK)
INTEGER(KIND=JPIM) :: NYCHK(JPGPCHK)
INTEGER(KIND=JPIM) :: NNFCHK(JPFLDCHK)
REAL(KIND=JPRB),ALLOCATABLE:: TENDCHK(:)

INTEGER(KIND=JPIM) :: NFRQCHK
INTEGER(KIND=JPIM) :: NFLDCHK
INTEGER(KIND=JPIM) :: NGPCHK
INTEGER(KIND=JPIM) :: NLENCHK
INTEGER(KIND=JPIM) :: NCHKTEND
!***********************************************************
END MODULE YOMCHK
