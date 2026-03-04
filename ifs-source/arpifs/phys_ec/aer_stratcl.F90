! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_STRATCL&
 &(YDVAB, YDCVER,YDDIMV, YDEAERD, YDERAD, YDEAERSNK, PRTAEVO, KIDIA , KFDIA, KLON  , KLEV , KWAVL, &
 &PAPH  , PGEMU, PRHCL, PTH ,&
 &PAERCLIS&
 &)  

!***********************************************************************
! CAUTION: THIS ROUTINE WORKS ONLY ON A NON-ROTATED, UNSTRETCHED GRID
!***********************************************************************

!**** *AER_STRATCL - STRATOSPHERIC BACKGROUND AND VOLCANIC AEROSOLS
!                    ADAPTED FROM *RADACT*

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *AER_STRATCL* FROM *AER_SRC*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     JJMorcrette  E.C.M.W.F.   20061120

!     MODIFICATIONS.
!     --------------
!     JJMorcrette  E.C.M.W.F.   20061120

!-----------------------------------------------------------------------

USE YOMVERT  , ONLY : TVAB, VP00
USE YOMCVER  , ONLY : TCVER
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RPI
USE YOEAERD  , ONLY : TEAERD
USE YOEAERSNK, ONLY : TEAERSNK
USE YOEAEROP , ONLY : ALF_SU 
USE YOERAD   , ONLY : TERAD

IMPLICIT NONE

!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

TYPE(TVAB)        ,INTENT(IN)   :: YDVAB
TYPE(TCVER)       ,INTENT(IN)   :: YDCVER
TYPE(TDIMV)       ,INTENT(IN)   :: YDDIMV
TYPE(TEAERD)      ,INTENT(INOUT):: YDEAERD
TYPE(TEAERSNK)    ,INTENT(INOUT):: YDEAERSNK
TYPE(TERAD)       ,INTENT(INOUT):: YDERAD
REAL(KIND=JPRB)   ,INTENT(IN)   :: PRTAEVO(46)
INTEGER(KIND=JPIM),INTENT(IN)   :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)   :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)   :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)   :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)   :: KWAVL 

REAL(KIND=JPRB)   ,INTENT(IN)   :: PAPH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PGEMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PRHCL(KLON,KLEV), PTH(KLON,0:KLEV) 

REAL(KIND=JPRB)   ,INTENT(OUT) :: PAERCLIS(KLON,KLEV,2) 

!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------

INTEGER(KIND=JPIM) :: JTAB
INTEGER(KIND=JPIM) :: IRH(KLON,KLEV) 

REAL(KIND=JPRB) :: ZAERVO(KLON)

REAL(KIND=JPRB) :: ZAETRN(KLON),ZAETRO(KLON)
REAL(KIND=JPRB) :: ZAER1(KLON,KLEV), ZAER2(KLON,KLEV)

REAL(KIND=JPRB) :: ZDPN(KLON)  , ZDPO(KLON)
REAL(KIND=JPRB) :: ZGRTH(KLON)

REAL(KIND=JPRB) :: ZSILAT(KLON), ZSINR(46)

INTEGER(KIND=JPIM) :: INLA, JK, JL, JLR, ILATR 

REAL(KIND=JPRB) :: ZAETR, ZDPNMO, ZFACT, ZGRIDR, ZLATR, ZSIN 
REAL(KIND=JPRB) :: ZALF, ZEPSAER

!-- security related quantities
REAL(KIND=JPRB) :: ZEPSBGA, ZCTRBGA, ZCSTBGA
REAL(KIND=JPRB) :: ZETAC(KLEV)    , ZETAHC(0:KLEV)    , ZDPRC(KLEV)
REAL(KIND=JPRB) :: ZPRES(KLEV)    , ZPRESH(0:KLEV)    , ZEPSBG(KLEV)
REAL(KIND=JPRB) :: ZCONV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!#include "legtri.intfb.h"
#include "gphpre.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_STRATCL',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & RCTRPT=>YDEAERD%RCTRPT, RCVOBGA=>YDEAERD%RCVOBGA, &
 & RRHTAB=>YDEAERSNK%RRHTAB, &
 & LHVOLCA=>YDERAD%LHVOLCA)
!*         0.     SECURITY BACKGROUND AEROSOLS

ZEPSBGA=1.E-06_JPRB/101325._JPRB
ZCTRBGA=0.03_JPRB/(101325._JPRB-19330._JPRB)
ZCSTBGA=0.045_JPRB/19330._JPRB
ZPRESH(KLEV)=VP00

!-- convert N.D. optical thickness to concentration in kg m-2, given an
!   mass absorption coefficient in m2 g-1
ZCONV=1.E-03_JPRB

CALL GPHPRE(1,NFLEVG,1,1,YDVAB,YDCVER,ZPRESH,PRESF=ZPRES)

DO JK=0,KLEV
  ZETAHC(JK)=ZPRESH(JK)/ZPRESH(KLEV)
ENDDO
DO JK=1,KLEV
  ZETAC(JK)=ZPRES(JK)/ZPRESH(KLEV)
  ZDPRC(JK)=ZPRES(JK)-ZPRES(JK-1)
  ZEPSBG(JK)=ZDPRC(JK)*ZEPSBGA
ENDDO



!     ------------------------------------------------------------------

!*         1.     "NEW AEROSOL DISTRIBUTION" PARAMETERS COMPUTATIONS
!                 --------------------------------------------------

IF (LHVOLCA) THEN


!*         1.1    VOLCANIC AEROSOL DISTRIBUTION PARAMETERS FROM GISS CLIMATOLOGY
!                 --------------------------------------------------------------

  ILATR=46
  ZGRIDR=180._JPRB/(ILATR-1)
  DO JLR=1,ILATR
    ZLATR=90._JPRB-(JLR-1)*ZGRIDR
    ZSINR(JLR)=SIN(ZLATR*RPI/180._JPRB)
  ENDDO

  DO JL=KIDIA,KFDIA
    INLA=0
    ZSILAT(JL)=-9999._JPRB
    ZSIN=PGEMU(JL)
    DO JLR=1,ILATR-1
      IF (ZSIN <= ZSINR(JLR) .AND. ZSIN > ZSINR(JLR+1)) THEN
        INLA=JLR
        ZSILAT(JL)=(ZSIN-ZSINR(INLA))/(ZSINR(INLA+1)-ZSINR(INLA))
        ZAERVO(JL)=PRTAEVO(INLA)+ZSILAT(JL)*(PRTAEVO(INLA+1)-PRTAEVO(INLA))
      ENDIF
    ENDDO
    IF (ZSIN <= ZSINR(ILATR-1) .AND. ZSIN >= ZSINR(ILATR))THEN
      INLA=ILATR
      ZSILAT(JL)=(ZSIN-ZSINR(INLA-1))/(ZSINR(INLA)-ZSINR(INLA-1))
      ZAERVO(JL)=PRTAEVO(INLA-1)&
       & +ZSILAT(JL)*(PRTAEVO(INLA)-PRTAEVO(INLA-1))  
    ENDIF
    IF (INLA == 0) THEN
      CALL ABOR1(' Problem with lat. interpolation in aer_stratcl!')
    ENDIF
! PRTAEVO is now total optical depth, so need to rescale for this routine
    ZAERVO(JL)=ZAERVO(JL)/19330._JPRB
  ENDDO

ELSE
!*         1.2    TANRE ET AL. CLIMATOLOGY
!                 ------------------------

  DO JL=KIDIA,KFDIA
    ZAERVO(JL)=RCVOBGA
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.      VERTICAL DISTRIBUTION
!*               ---------------------

DO JL=KIDIA,KFDIA
  ZDPO(JL)=PAPH(JL,0)
  ZAETRO(JL)=1.0_JPRB
ENDDO

DO JK=1,KLEV
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA
    PAERCLIS(JL,JK,1) = 0._JPRB
    PAERCLIS(JL,JK,2) = 0._JPRB
    ZGRTH(JL)= PTH(JL,JK-1)/PTH(JL,JK)
    ZDPN(JL)=PAPH(JL,JK)

    IF (0.5_JPRB*(PAPH(JL,JK-1)+PAPH(JL,JK)) < 999._JPRB) THEN
! for models with top above 10hPa
      ZAETRN(JL)=1.0_JPRB
      ZAETRO(JL)=1.0_JPRB
      ZEPSAER=ZEPSBG(JK)
    ELSE
      ZAETRN(JL)=ZAETRO(JL)*(MIN(1.0_JPRB, ZGRTH(JL) ))**RCTRPT
      ZEPSAER=1.E-10_JPRB
      ZEPSAER=ZEPSBG(JK)
    ENDIF

    ZAETR=SQRT  (ZAETRN(JL)*ZAETRO(JL))
    ZDPNMO    =ZDPN(JL)-ZDPO(JL)
!- stratospheric background
    ZAER1(JL,JK)=   ZAETR * ZCSTBGA*ZDPNMO
!- volcanic (either GISS or Tanre et al.)
    ZAER2(JL,JK)=   ZAETR * ZAERVO(JL) * ZDPNMO
   ENDDO
  DO JL=KIDIA,KFDIA
    ZDPO(JL)  =ZDPN(JL)
    ZAETRO(JL)=ZAETRN(JL)
  ENDDO
ENDDO

!-- transform the optical thickness in mass mixing ratio

IRH(:,:)=1
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    DO JTAB=1,12
      IF (PRHCL(JL,JK)*100._JPRB > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO

!-- get the relevant absorption coefficient for sulfate
    ZALF = ALF_SU(IRH(JL,JK),KWAVL)
    ZAER1(JL,JK) = MAX(ZEPSAER, ZAER1(JL,JK))
    ZAER2(JL,JK) = MAX(ZEPSAER, ZAER2(JL,JK))

! - tau = abscoef * mmr * DeltaP / RG

    ZFACT = ZCONV * RG / (ZALF * (PAPH(JL,JK)-PAPH(JL,JK-1)))

    PAERCLIS(JL,JK,1)=ZAER1(JL,JK) * ZFACT
    PAERCLIS(JL,JK,2)=ZAER2(JL,JK) * ZFACT
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_STRATCL',1,ZHOOK_HANDLE)
END SUBROUTINE AER_STRATCL
