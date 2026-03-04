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

SUBROUTINE SUNUD(YDARPHY,KULOUT)

!**** *SUNUD* * - ROUTINE TO INITIALIZE COEFS FOR NUDGING

!     PURPOSE.
!     --------
!        SET DEFAULT VALUES, THEN READS NAMELIST NAMNUD

!**   INTERFACE.
!     ----------
!        *CALL* *SUNUD(...)*

!     EXPLICIT ARGUMENTS :  KULOUT
!     --------------------

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMNUD

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!      Michel Deque *CNRM*
!      ORIGINAL : 96-07-01

!     MODIFICATIONS.
!     --------------
!      M. Deque 01-06-11 : nudging with time-variable coefficients
!      P. Marquet 02-02-14 : LWNUDG in YOMNUD
!      M. Deque 05-06-20 : grid point nudging                
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Jul-2006 Revised surface fields
!      A.Alias       31-Oct-2006 Modified setup of LNUDST(=LNUDG, even if XNUDST=0.0)
!      A.Alias       07-Mar-2007 nudging with wavenumber-variable coefficients(M.Deque)
!      A.Alias       19-Oct-2007 mv setup of NTOTFNUDG to sudim1
!      A.Alias       07-Aug-2008 NFRNUDG (frequency of nudging) : added properly
!      F.Chauvin     08-Aug-2008 bugfix LNUDST=LNUDG even if XNUDST=0.
!      A.Alias       Mar-2011    NFRNUDG setup to 999 instead of 0 and print format modified
!                                No nudging of surface fields when LMSE=.T.
!                                Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMARPHY , ONLY : TARPHY
USE YOMNUD   , ONLY : NFNUDG, NTNUDG, NFRNUDG, LNUDG, LNUDDI, LNUDLP, LNUDRM, &
 & LNUDSD, LNUDSH, LNUDSM, LNUDST, LNUDSV, LNUDTE, LNUDVO, XNUDDI, &
 & XNUDLP, XNUDRM, XNUDSD, XNUDSH, XNUDSM, XNUDST, XNUDSV, XNUDTE, XNUDVO, &
 & LNDRIV, LWNUDG, LNUDTG, LNUDQG, LNUDUG, LNUDVG, LNUDPG, &
 & XNUDTG, XNUDQG, XNUDUG, XNUDVG, XNUDPG, XNUVERT, &
 & NSPNU1, NSPNU2, NTOTFNUDG3, NTOTFNUDG2

!      ----------------------------------------------------------------

IMPLICIT NONE

TYPE(TARPHY)      ,INTENT(INOUT) :: YDARPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!      ----------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "namnud.nam.h"

!      ----------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNUD',0,ZHOOK_HANDLE)
ASSOCIATE(LMSE=>YDARPHY%LMSE)
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

!        1.1 Set implicit default values

LNUDG= .FALSE.
LNUDST =.FALSE.
LWNUDG =.FALSE.
LNDRIV =.FALSE.
NFNUDG=0
NFRNUDG=999
NTNUDG=1
NSPNU1=998
NSPNU2=999
XNUDTE =0._JPRB
XNUDSH =0._JPRB
XNUDDI =0._JPRB
XNUDVO =0._JPRB
XNUDSV =0._JPRB
XNUDLP =0._JPRB
XNUDST =0._JPRB
XNUDSM =0._JPRB
XNUDRM =0._JPRB
XNUDSD =0._JPRB

XNUDTG =0._JPRB
XNUDQG =0._JPRB
XNUDUG =0._JPRB
XNUDVG =0._JPRB
XNUDPG =0._JPRB

XNUVERT =1._JPRB
!      ----------------------------------------------------------------

!*       2.    MODIFIES DEFAULT VALUES.
!              ------------------------

CALL POSNAM(NULNAM,'NAMNUD')
READ       (NULNAM, NAMNUD)

!ZEPS=1.E-12_JPRB

LNUDTE=(XNUDTE /= 0.0_JPRB)
LNUDSH=(XNUDSH /= 0.0_JPRB)
LNUDDI=(XNUDDI /= 0.0_JPRB)
LNUDVO=(XNUDVO /= 0.0_JPRB)
LNUDSV=(XNUDSV /= 0.0_JPRB)
LNUDLP=(XNUDLP /= 0.0_JPRB)
IF(LMSE) THEN
  XNUDST=0.0_JPRB
  XNUDSM=0.0_JPRB
  XNUDRM=0.0_JPRB
  XNUDSD=0.0_JPRB
ENDIF
LNUDST=(XNUDST /= 0.0_JPRB)
LNUDSM=(XNUDSM /= 0.0_JPRB)
LNUDRM=(XNUDRM /= 0.0_JPRB)
LNUDSD=(XNUDSD /= 0.0_JPRB)

LNUDTG=(XNUDTG /= 0.0_JPRB)
LNUDQG=(XNUDQG /= 0.0_JPRB)
LNUDUG=(XNUDUG /= 0.0_JPRB)
LNUDVG=(XNUDVG /= 0.0_JPRB)
LNUDPG=(XNUDPG /= 0.0_JPRB)

IF(.NOT.LMSE) LNUDST=LNUDG

LNUDG=(LNUDG.AND.(LNUDTE.OR.LNUDSH.OR.LNUDDI&
 & .OR.LNUDVO.OR.LNUDSV.OR.LNUDLP&
 & .OR.LNUDST.OR.LNUDSM.OR.LNUDRM.OR.LNUDSD&
 & .OR.LNUDTG.OR.LNUDQG.OR.LNUDUG.OR.LNUDVG.OR.LNUDPG))  

 
NTOTFNUDG3=0
IF(LNUDG) THEN
  IF(LNUDTG)NTOTFNUDG3=NTOTFNUDG3+1
  IF(LNUDQG)NTOTFNUDG3=NTOTFNUDG3+1
  IF(LNUDUG)NTOTFNUDG3=NTOTFNUDG3+1
  IF(LNUDVG)NTOTFNUDG3=NTOTFNUDG3+1
ENDIF

NTOTFNUDG2=0
IF(LNUDG) THEN
  IF(LNUDST)NTOTFNUDG2=NTOTFNUDG2+1
  IF(LNUDSM)NTOTFNUDG2=NTOTFNUDG2+1
  IF(LNUDRM)NTOTFNUDG2=NTOTFNUDG2+1
  IF(LNUDSD)NTOTFNUDG2=NTOTFNUDG2+1
  IF(LNUDPG)NTOTFNUDG2=NTOTFNUDG2+1
ENDIF

!      -----------------------------------------------------------

!*       3.    PRINTS FINAL VALUES.
!              --------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMNUD '')')
WRITE(UNIT=KULOUT,FMT='(&
 & '' LNUDG  ='',L2,9X,'' NFNUDG ='',I2,9X,'' NTNUDG ='',I2,9X,/,&
 & '' NFRNUDG='',I3,9X,'' NTOTFNUDG3='',I4,9X,'' NTOTFNUDG2='',I4,9X,/,&
 & '' XNUDDI ='',E11.4,'' XNUDLP ='',E11.4,'' XNUDRM ='',E11.4,&
 & '' XNUDSD ='',E11.4,'' XNUDSH ='',E11.4,'' XNUDSM ='',E11.4,/,&
 & '' XNUDST ='',E11.4,&
 & '' XNUDSV ='',E11.4,'' XNUDTE ='',E11.4,'' XNUDVO ='',E11.4,/,&
 & '' LNUDDI ='',L2,9X,'' LNUDLP ='',L2,9X,'' LNUDRM ='',L2,9X,&
 & '' LNUDSD ='',L2,9X,'' LNUDSH ='',L2,9X,'' LNUDSM ='',L2,9X,/,&
 & '' LNUDST ='',L2,9X,&
 & '' LNUDSV ='',L2,9X,'' LNUDTE ='',L2,9X,'' LNUDVO ='',L2,9X,/,&
 & '' LWNUDG ='',L2,9X&
 & )')&
 & LNUDG,  NFNUDG, NTNUDG, NFRNUDG, NTOTFNUDG3, NTOTFNUDG2,&
 & XNUDDI, XNUDLP, XNUDRM, XNUDSD, XNUDSH, XNUDSM, XNUDST,&
 & XNUDSV, XNUDTE, XNUDVO,&
 & LNUDDI, LNUDLP, LNUDRM, LNUDSD, LNUDSH, LNUDSM, LNUDST,&
 & LNUDSV, LNUDTE, LNUDVO, LWNUDG  
WRITE(UNIT=KULOUT,FMT='(&
 & '' LNDRIV  ='',L2,9X,/,&
 & '' XNUDTG ='',E11.4,'' XNUDQG ='',E11.4,'' XNUDUG ='',E11.4,&
 & '' XNUDVG ='',E11.4,'' XNUDPG ='',E11.4,/,&
 & '' LNUDTG ='',L2,9X,'' LNUDQG ='',L2,9X,'' LNUDUG ='',L2,9X,&
 & '' LNUDVG ='',L2,9X,'' LNUDPG ='',L2,9X&
 & )')&
 & LNDRIV,&
 & XNUDTG,XNUDQG,XNUDUG,XNUDVG,XNUDPG,&
 & LNUDTG,LNUDQG,LNUDUG,LNUDVG,LNUDPG  
WRITE(UNIT=KULOUT,FMT='('' XNUVERT ='')')
WRITE(UNIT=KULOUT,FMT='(4(1X,E19.12))') XNUVERT
WRITE(UNIT=KULOUT,FMT='(" NSPNU1 =",I5," NSPNU2 =",I5)')NSPNU1,NSPNU2

!      -----------------------------------------------------------

!*       4.    VARIOUS CHECKINGS.
!              ------------------

IF(LNUDG.AND.(NFNUDG == 0)) THEN
  WRITE(KULOUT,FMT='('' WHEN NUDGING, NFNUDG MUST BE >0 '')')
  CALL ABOR1('SUNUD: ABOR1 CALLED')
ENDIF

!      ----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNUD',1,ZHOOK_HANDLE)
END SUBROUTINE SUNUD
