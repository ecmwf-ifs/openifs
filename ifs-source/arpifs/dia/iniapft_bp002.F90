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

SUBROUTINE INIAPFT_BP002(YDLDDH,YDSDDH,YDPHY,KLUNOUT)
  
!     Purpose.
!     --------
!     Preparing descriptors of fluxes or tendencies, for usage in DDH. 
     

!   Interface.
!   ----------
!      CALL INIAPFT_BP002

!    Explicit arguments
!    ------------------
! KLUNOUT  : logical number of standard uotput unit

!     Author.
!     -------
!      Tomislav Kovacic

!     Modifications.
!     --------------
!      Original : 2007-06-01
!      Sep-2009, M. Hrastinski: TKE and TTE phys. fluxes for ALARO DDH
!-------------------------------------------------------------------------
  
USE PARKIND1, ONLY : JPIM     ,JPRB
USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLDDH,  ONLY : TLDDH
USE YOMSDDH,  ONLY : TSDDH
USE YOMPHY,   ONLY : TPHY
!USE YOMPHY,   ONLY  : LSTRAPRO, L3MT !!!, LRAYFM, &
!                    & LCONDWT  ,LPROCLD  ! Lopez???
USE YOMPHFT,  ONLY  : NDDHFT, YAPFT
!     ------------------------------------------------------------------
IMPLICIT NONE

TYPE(TLDDH)        ,INTENT(INOUT):: YDLDDH
TYPE(TPHY)         ,INTENT(IN)   :: YDPHY
TYPE(TSDDH)        ,INTENT(INOUT):: YDSDDH
INTEGER(KIND=JPIM), INTENT(IN) :: KLUNOUT 

CHARACTER(LEN=1)   :: CLAROFT
INTEGER(KIND=JPIM) :: JI, JJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

!#include "abor1.intfb.h"
#include "addft.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INIAPFT_BP002',0,ZHOOK_HANDLE)
ASSOCIATE(LHDHKS=>YDLDDH%LHDHKS, LONLYVAR=>YDLDDH%LONLYVAR, &
 & LHDQLN=>YDSDDH%LHDQLN,NDPSFI=>YDPHY%NDPSFI, &
 & LPTKE=>YDPHY%LPTKE, LCOEFK_PTTE=>YDPHY%LCOEFK_PTTE)
!     ------------------------------------------------------------------

WRITE(KLUNOUT,*) 'Inside of INIAPFT_BP002' 

NDDHFT = 0

! Hardcoded switches
CLAROFT = 'F'
LONLYVAR = .FALSE.

JJ=0
  
!                 Descriptors of ALARO physical fluxes

IF ( LHDHKS ) THEN         

  NDDHFT = NDDHFT + 26
  IF(LHDQLN) THEN         
    NDDHFT = NDDHFT +10 
!    IF (LCONDWT.AND.LPROCLD) THEN         
      NDDHFT = NDDHFT +10 
!    ENDIF
  ENDIF
  IF (LPTKE) THEN
    NDDHFT = NDDHFT + 4
    IF(LCOEFK_PTTE) THEN
       NDDHFT = NDDHFT + 2
    ENDIF
  ENDIF

  ALLOCATE(YAPFT(NDDHFT))
  CALL ADDFT(KLUNOUT,CLAROFT, 'PP', 'SUMFPL     ' ,JJ)!ZDPSFI*(ZFPLL(JROF,JLEV)
                                                      !       +ZFPLN(JROF,JLEV))
                                                      !          IDHCV+ 1
  CALL ADDFT(KLUNOUT,CLAROFT, 'QV', 'DIFT      ' ,JJ)      ! PDIFTQ   IDHCV+ 2
  CALL ADDFT(KLUNOUT,CLAROFT, 'QV', 'DIFC      ' ,JJ)      ! PDIFCQ   IDHCV+ 3
  CALL ADDFT(KLUNOUT,CLAROFT, 'QR', 'PLS       ' ,JJ)      ! PFPLSL   IDHCV+ 4
  CALL ADDFT(KLUNOUT,CLAROFT, 'QS', 'PLS       ' ,JJ)      ! PFPLSN   IDHCV+ 5
  CALL ADDFT(KLUNOUT,CLAROFT, 'QR', 'PLC       ' ,JJ)      ! PFPLCL   IDHCV+ 6
  CALL ADDFT(KLUNOUT,CLAROFT, 'QS', 'PLC       ' ,JJ)      ! PFPLCN   IDHCV+ 7
  CALL ADDFT(KLUNOUT,CLAROFT, 'QV', 'NG        ' ,JJ)      ! PFCQNG   IDHCV+ 8
  CALL ADDFT(KLUNOUT,CLAROFT, 'UU', 'TUR       ' ,JJ)      ! ZSTTUG   IDHCV+ 9 ?
  CALL ADDFT(KLUNOUT,CLAROFT, 'VV', 'TUR       ' ,JJ)      ! ZSTTVG   IDHCV+10 ?
  CALL ADDFT(KLUNOUT,CLAROFT, 'UU', 'TURCONV   ' ,JJ)      ! ZSTCUG   IDHCV+11 ?
  CALL ADDFT(KLUNOUT,CLAROFT, 'VV', 'TURCONV   ' ,JJ)      ! ZSTCVG   IDHCV+12 ?
  CALL ADDFT(KLUNOUT,CLAROFT, 'UU', 'ONDEGREL  ' ,JJ)      ! ZSTDUG   IDHCV+13 ?
  CALL ADDFT(KLUNOUT,CLAROFT, 'VV', 'ONDEGREL  ' ,JJ)      ! ZSTDVG   IDHCV+14 ?
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'DIFT      ' ,JJ)      ! PDIFTS   IDHCV+15
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'DIFC      ' ,JJ)      ! PDIFCS   IDHCV+16
!  IF (LRAYFM) THEN
    CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'RSO       ' ,JJ)    ! PFRSO    IDHCV+17
    CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'RTH       ' ,JJ)    ! PFRTH    IDHCV+18
!  ENDIF
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HSSL      ' ,JJ)      ! PFHSSL   IDHCV+25
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HSSN      ' ,JJ)      ! PFHSSN   IDHCV+26
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HSCL      ' ,JJ)      ! PFHSCL   IDHCV+23
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HSCN      ' ,JJ)      ! PFHSCN   IDHCV+24
  IF (NDPSFI==1) THEN
    CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'CMP       ' ,JJ)    ! 
    CALL ADDFT(KLUNOUT,CLAROFT, 'QV', 'CMPC      ' ,JJ)    ! 
  ENDIF
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HPSL      ' ,JJ)      ! PFHPSL   IDHCV+19
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HPSN      ' ,JJ)      ! PFHPSN   IDHCV+20
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HPCL      ' ,JJ)      ! PFHPCL   IDHCV+21
  CALL ADDFT(KLUNOUT,CLAROFT, 'CT', 'HPCN      ' ,JJ)      ! PFHPCN   IDHCV+22
  IF(LHDQLN) THEN         
    CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'DIFT      ' ,JJ)    ! PDIFTQL  IDHCV+27
    CALL ADDFT(KLUNOUT,CLAROFT, 'QI', 'DIFT      ' ,JJ)    ! PDIFTQN  IDHCV+28
    CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'DIFC      ' ,JJ)    ! PDIFCQL  IDHCV+29
    CALL ADDFT(KLUNOUT,CLAROFT, 'QI', 'DIFC      ' ,JJ)    ! PDIFCQN  IDHCV+30
    CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'NG        ' ,JJ)    ! PFCQLNG  IDHCV+31
    CALL ADDFT(KLUNOUT,CLAROFT, 'QI', 'NG        ' ,JJ)    ! PFCQNNG  IDHCV+32
    CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'CC        ' ,JJ)    ! PFCQLF   IDHCV+33
    CALL ADDFT(KLUNOUT,CLAROFT, 'QN', 'CC        ' ,JJ)    ! PFCQIF   IDHCV+34
    CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'CS        ' ,JJ)    ! PFSQLF   IDHCV+35
    CALL ADDFT(KLUNOUT,CLAROFT, 'QN', 'CS        ' ,JJ)    ! PFSQIF   IDHCV+36
!    IF (LCONDWT.AND.LPROCLD) THEN         
      CALL ADDFT(KLUNOUT,CLAROFT, 'QR', 'PEVS      ' ,JJ)  ! PFPEVPL  IDHCV+37
      CALL ADDFT(KLUNOUT,CLAROFT, 'QS', 'PEVS      ' ,JJ)  ! PFPEVPN  IDHCV+38
      CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'PFPS      ' ,JJ)  ! PFPFPL   IDHCV+39
      CALL ADDFT(KLUNOUT,CLAROFT, 'QN', 'PFPS      ' ,JJ)  ! PFPFPN   IDHCV+40
      CALL ADDFT(KLUNOUT,CLAROFT, 'QR', 'PEVC      ' ,JJ)  ! PFPEVPCL IDHCV+41
      CALL ADDFT(KLUNOUT,CLAROFT, 'QS', 'PEVC      ' ,JJ)  ! PFPEVPCN IDHCV+42
      CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'PFPC      ' ,JJ)  ! PFPFPCL  IDHCV+43
      CALL ADDFT(KLUNOUT,CLAROFT, 'QN', 'PFPC      ' ,JJ)  ! PFPFPCN  IDHCV+44
      CALL ADDFT(KLUNOUT,CLAROFT, 'QR', 'NG        ' ,JJ)  ! PFCQRNG  IDHCV+45
      CALL ADDFT(KLUNOUT,CLAROFT, 'QS', 'NG        ' ,JJ)  ! PFCQSNG  IDHCV+46
      IF (NDPSFI==1) THEN
        CALL ADDFT(KLUNOUT,CLAROFT, 'QL', 'CMP       ' ,JJ)! 
        CALL ADDFT(KLUNOUT,CLAROFT, 'QN', 'CMP       ' ,JJ)! 
      ENDIF
!    ENDIF         
  ENDIF         
  IF (LPTKE) THEN
!                                                          Index in 
!   FTKESHRPROD; shear production term.                     CPPHDDH
    CALL ADDFT(KLUNOUT,CLAROFT, 'TK', 'ESHRPROD  ' ,JJ)  ! IDHCV+47
!   FTKEBUOYPROD; buoyancy production/destruction term.
    CALL ADDFT(KLUNOUT,CLAROFT, 'TK', 'EBUOYPROD ' ,JJ)  ! IDHCV+48
!   FTKEDISSIP; TKE dissipation term.
    CALL ADDFT(KLUNOUT,CLAROFT, 'TK', 'EDISSIP   ' ,JJ)  ! IDHCV+49
!   FTKETRANSPORT; turbulent transport of TKE
    CALL ADDFT(KLUNOUT,CLAROFT, 'TK', 'ETRANSP   ' ,JJ)  ! IDHCV+50
    IF(LCOEFK_PTTE) THEN
!     FTTEDISSIP; TTE dissipation term.
      CALL ADDFT(KLUNOUT,CLAROFT, 'TT', 'EDISSIP   ' ,JJ)! IDHCV+51
!     FTTETRANSPORT; turbulent transport of TTE
      CALL ADDFT(KLUNOUT,CLAROFT, 'TT', 'ETRANSP   ' ,JJ)! IDHCV+52
    ENDIF
  ENDIF
ENDIF

WRITE(KLUNOUT,*)' DDH prepared in INIAPFT_BP002'
WRITE(KLUNOUT,*)' DDH physical fluxes descriptors' 
DO JI= 1,NDDHFT 
  WRITE(KLUNOUT,*) JI,'  ',YAPFT(JI)%CFT,YAPFT(JI)%CVAR,YAPFT(JI)%CNAME 
ENDDO

WRITE(KLUNOUT,*) 'END of INIAPFT_BP002' 

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INIAPFT_BP002',1,ZHOOK_HANDLE)

END SUBROUTINE INIAPFT_BP002
