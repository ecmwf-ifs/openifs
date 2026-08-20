! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WVWAMINIT1(LDWCOUIFS, LDWCOU2W, LDWCOURNW, LDWCOUHMF, LDWFLUX, LFDBOPIFS)

! ----------------------------------------------------------------------

!**** *WVWAMINIT1* - WAVE MODEL CONTINUED INITIALISATION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W  ,LWCOURNW, LWCOUHMF, LWFLUX
      USE YOWCOUT  , ONLY : LFDB
      USE YOWMESPAS, ONLY : LFDBIOOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,    NNEXT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL, INTENT(IN) :: LDWCOUIFS, LDWCOU2W, LDWCOURNW, LDWCOUHMF, LDWFLUX, LFDBOPIFS

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',0,ZHOOK_HANDLE)

! RE-INITIALIZE LOGICALS IF COUPLED TO IFS

      IF (LDWCOUIFS) THEN
        LWCOU=LDWCOUIFS
        LWCOU2W=LDWCOU2W
        LWCOURNW=LDWCOURNW
        LWCOUHMF=LDWCOUHMF
        LWFLUX=LDWFLUX
        LFDB=LFDBOPIFS
        LFDBIOOUT=LFDBOPIFS
      ENDIF

      NPREVIOUS=IRANK-1
      IF (IRANK == NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVWAMINIT1
