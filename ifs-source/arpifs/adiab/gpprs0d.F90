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

SUBROUTINE GPPRS0D(KPROMA,KSTART,KPROF,KFLEV,PT,PQ,PQCLC,PQCLR,&
 & PQCLI,PQCLS,PQCLG,LDUCONV,LDSREC,LDSREDBC,PSIMRFC,PSIMRFCDB)

USE PARKIND1  , ONLY : JPRD, JPIM,  JPRB
USE YOMCST    , ONLY : RPI,RTT
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCLMICST, ONLY : RXCCR    ,RXAS     ,RXBS     ,RXCCS   ,RXCXS    ,RXBI   &
 &                    ,RXAI     ,RXLBI   ,RXLBEXI  ,RXAG     ,RXBG     ,RXCCG   ,RXCXG    ,RXLBEXR  &
 &                    ,RXLBR    ,RXLBEXS ,RXLBS    ,RXLBEXG  ,RXLBG    ,RXRTMIN  &
 &                    ,RMOMG_RXALPHAR_RXNUR_6      ,RMOMG_RXALPHAI_RXNUI_ZEXP &
 &                    ,RMOMG_RXALPHAS_RXNUS_ZEXP   ,RMOMG_RXALPHAG_RXNUG_ZEXP

!**** *GPPRS0D*  - Compute Simulated reflectivities - 0D -

!     PURPOSE.
!     --------
!        To Compute reflectivities knowing precipitating hydrometeores. 
!        The computation is done in 0Dimension, the geometry of the radar
!        is not considered. 
!**   INTERFACE.
!     ----------
!       *CALL* *GPPRS0D*

!        EXPLICIT ARGUMENTS
!        --------------------
!            INPUT :
!        KPROMA    : Horizontal dimension
!        KSTART    : start of work
!        KPROF     : depth of work
!        KFLEV     : number of vertical levels
!        PT        : temperature (KPROMA,KFLEV)
!        PQ        ; specific humidity (KPROMA,KFLEV)
!        PQCLC      : cloud water (KPROMA,KFLEV)
!        PQCLR      : cloud rain (KPROMA,KFLEV)
!        PQCLI      : cloud ice (KPROMA,KFLEV)
!        PQCLS      : cloud snow (KPROMA,KFLEV)
!        PQCLG      : cloud graupel (KPROMA,KFLEV)
!            OUTPUT:
!        PSIMRFC    : simulated reflectivies in mm/h(KPROMA,KLEV)                            
!        PSIMRFCDB  : simulated reflectivies in dBZ (KPROMA,KLEV)  

!        IMPLICIT ARGUMENTS
!        --------------------
!           NONE

!     METHOD.
!     -------
!        Derived from the radar obs operator (routine reflsim by E. Wattrelot)

!     EXTERNALS.
!     ----------
   
!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        Gwenaelle Hello *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 06-07-31
!        K. Yessad (Sep 2008): make coding norms compliant.
!        R. El Khatib  20-Mar-2009 Optimisation & bugfix
!        R. El Khatib  29-Mar-2010 celebrate last anniversary by more
!                                  optimisations ;-)
!        F. Vana    29-Sep-2014  Support for single precision
!        R. El Khatib  25-May-2016 Optimisation again
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KPROMA,KSTART,KPROF,KFLEV
REAL(KIND=JPRB)   ,INTENT(IN)  :: PT(KPROMA,KFLEV),PQ(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PQCLC(KPROMA,KFLEV),PQCLR(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PQCLI(KPROMA,KFLEV),PQCLS(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PQCLG(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSIMRFC(KPROMA,KFLEV),PSIMRFCDB(KPROMA,KFLEV)
LOGICAL           ,INTENT(IN)  :: LDUCONV
LOGICAL,INTENT(IN) :: LDSREC,LDSREDBC

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZRCLR(KPROMA,KFLEV),ZRCLS(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRCLI(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRCLG(KPROMA,KFLEV),ZQD(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZLBDA,ZREFL(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZREFL_MELT_CONV(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZXR, ZRHOLW, ZDMELT_FACT, ZEXP, ZEQICE
REAL(KIND=JPRB) :: ZFRAC_WATER
REAL(KIND=JPRB) :: ZEQWAT, ZTRESH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JROF,JLEV

!     ------------------------------------------------------------------


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPRS0D',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
!*********************************************************
! 0. Initialization  
!*********************************************************

ZREFL(:,:)=0._JPRB
! Threshold value
IF (JPRB == JPRD) THEN
  ! Maintaining consistency with previous double prec. results
  ZTRESH=1.0E120_JPRD
ELSE
  ! General case
  ZTRESH=10._JPRB**(MAXEXPONENT(ZTRESH)/10)
ENDIF

! From qx to mixing ratio

DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF
    ZQD(JROF,JLEV)=1.0_JPRB-PQ(JROF,JLEV)-PQCLC(JROF,JLEV)-PQCLR(JROF,JLEV) &
     &-PQCLI(JROF,JLEV) - PQCLS(JROF,JLEV) - PQCLG(JROF,JLEV)
  ENDDO
ENDDO

DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF
    ZRCLR(JROF,JLEV)=PQCLR(JROF,JLEV)/ZQD(JROF,JLEV)
    ZRCLS(JROF,JLEV)=PQCLS(JROF,JLEV)/ZQD(JROF,JLEV)
    ZRCLG(JROF,JLEV)=PQCLG(JROF,JLEV)/ZQD(JROF,JLEV)
    ZRCLI(JROF,JLEV)=PQCLI(JROF,JLEV)/ZQD(JROF,JLEV)
  ENDDO
ENDDO 

! Initialisation of clouds constants (Duplication of what is done in

ZXR = -1.0_JPRB !Value of x_r in N_r = C_r lambda_r ** x_r
ZRHOLW = 1000.0_JPRB !density of liquid water in kgm-3
ZEQICE = 0.224_JPRB !?

!**********************************************************************
! 1. Compute reflectivities on each model vertical level (0D operator)*
!**********************************************************************

! 1.1. Rain
! -----------

DO JLEV=1,KFLEV
!DIR$ NOVECTOR
  DO JROF=KSTART,KPROF
    IF (ZRCLR(JROF,JLEV) > RXRTMIN(3)) THEN
      ZLBDA = RXLBR*( ZRCLR(JROF,JLEV) )**RXLBEXR
      IF (ZLBDA < ZTRESH) THEN
        ZREFL(JROF,JLEV) = RXCCR*(ZLBDA**(ZXR-6.0_JPRB))&
         &*RMOMG_RXALPHAR_RXNUR_6
      ENDIF
    ENDIF
  ENDDO
ENDDO

! 1.2. Ice
! --------
!ZDMELT_FACT=(RXBI/3.0_JPRB)*((6.0_JPRB*RXAI)/(RPI*ZRHOLW))**(7.0_JPRB/3.0_JPRB)
!ZEXP = 7.0_JPRB*RXBI/3.0_JPRB-1.0_JPRB
ZDMELT_FACT = ( (6.0_JPRB*RXAI)/(RPI*ZRHOLW) )**2
ZEXP = 2.0_JPRB*RXBI
DO JLEV=1,KFLEV
!DIR$ NOVECTOR
  DO JROF=KSTART,KPROF
    IF (ZRCLI(JROF,JLEV) > RXRTMIN(4)) THEN
      ZLBDA = RXLBI*( ZRCLI(JROF,JLEV) / 800.0_JPRB )**RXLBEXI
      IF (ZLBDA < ZTRESH) THEN
        ZREFL(JROF,JLEV) = ZREFL(JROF,JLEV)+ZEQICE*ZDMELT_FACT  &
         &  *800.0_JPRB*(ZLBDA**(-ZEXP))*RMOMG_RXALPHAI_RXNUI_ZEXP
      ENDIF
    ENDIF
  ENDDO
ENDDO

! To be coded
!!$WHERE(PCIWF5(:,:) > XRTMINCIW)
!!$   ZDMELT_FACT = ( XBI/3.0 )*( (6.0*XAI)/(RPI*XRHOLW) )**(7.0/3.0)
!!$   ZEXP = 7.0*XBI/3.0-1.0
!!$   ZLBDA(:,:) = XLBI*( PCLWF5(:,:)/PCIT(:,:) )**XLBEXI
!!$   ZREFL(:,:) = ZREFL(:,:)+ZEQICE*ZDMELT_FACT                 &
!!$      &  *PCIT(:,:)*(ZLBDA(:,:)**(-ZEXP))*RMOMG_XALPHAI_XNUI_ZEXP
!!$END WHERE

! 1.3. Snow  
! ---------

!ZDMELT_FACT=(RXBS/3.0_JPRB)*((6.0_JPRB*RXAS)/(RPI*ZRHOLW))**(7.0_JPRB/3.0_JPRB)
!ZEXP = 7.0_JPRB*RXBS/3.0_JPRB-1.0_JPRB

ZDMELT_FACT = ( (6.0_JPRB*RXAS)/(RPI*ZRHOLW) )**2
ZEXP = 2.0_JPRB*RXBS

DO JLEV=1,KFLEV
!DIR$ NOVECTOR
  DO JROF=KSTART,KPROF
    IF (ZRCLS(JROF,JLEV) > RXRTMIN(5)) THEN
      ZLBDA = RXLBS*( ZRCLS(JROF,JLEV) )**RXLBEXS
      IF (ZLBDA < ZTRESH) THEN
        ZREFL(JROF,JLEV) = ZREFL(JROF,JLEV)+ZEQICE*ZDMELT_FACT      &
         &  *RXCCS*(ZLBDA**(RXCXS-ZEXP))*RMOMG_RXALPHAS_RXNUS_ZEXP
      ENDIF
    ENDIF
  ENDDO
ENDDO

! 1.4. Graupel
! ------------
ZFRAC_WATER = 0.14_JPRB
!ZDMELT_FACT=(RXBG/3.0_JPRB)*((6.0_JPRB*RXAG)/(RPI*ZRHOLW))**(7.0_JPRB/3.0_JPRB)
!ZEXP = 7.0_JPRB*RXBG/3.0_JPRB-1.0_JPRB
! new computation (Smith, 1984)
ZDMELT_FACT = ( (6.0_JPRB*RXAG)/(RPI*ZRHOLW) )**2
ZEXP = 2.0_JPRB*RXBG
!ZEQWAT=ZFRAC_WATER*(1.0_JPRB-ZEQICE)
ZEQWAT=1._JPRB/((1.0_JPRB-ZFRAC_WATER)*.92_JPRB+ZFRAC_WATER)**2

DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF
!    ZREFL_MELT_CONV(JROF,JLEV)=ZDMELT_FACT* &
!     & (ZEQICE + ZEQWAT*MAX(SIGN(1._JPRB,PT(JROF,JLEV)-RTT),0._JPRB))
     ZREFL_MELT_CONV(JROF,JLEV)=ZDMELT_FACT* &
      & (ZEQICE*(1._JPRB-MAX(SIGN(1._JPRB,PT(JROF,JLEV)-RTT),0._JPRB))&
      & + ZEQWAT*MAX(SIGN(1._JPRB,PT(JROF,JLEV)-RTT),0._JPRB))
  ENDDO
ENDDO

DO JLEV=1,KFLEV
!DIR$ NOVECTOR
  DO JROF=KSTART,KPROF
    IF (ZRCLG(JROF,JLEV) > RXRTMIN(6)) THEN 
      ZLBDA = RXLBG*( ZRCLG(JROF,JLEV) )**RXLBEXG
      IF (ZLBDA < ZTRESH) THEN
        ZREFL(JROF,JLEV) = ZREFL(JROF,JLEV)+ZREFL_MELT_CONV(JROF,JLEV) &
         & *RXCCG*(ZLBDA**(RXCXG-ZEXP))*RMOMG_RXALPHAG_RXNUG_ZEXP
      ENDIF
    ENDIF 
  ENDDO
ENDDO

! 1.5. Unit conversion
! --------------------

IF (LDUCONV.AND.LDSREC) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      PSIMRFC(JROF,JLEV) = (0.005_JPRB*1.E18_JPRB*ZREFL(JROF,JLEV))**0.625_JPRB ! Z_e in mmh-1
    ENDDO
  ENDDO
ELSE
  PSIMRFC(:,:)=ZREFL(:,:)
ENDIF
IF (LDUCONV.AND.LDSREDBC) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      PSIMRFCDB(JROF,JLEV) = 10._JPRB*LOG10( 1.E18_JPRB*MAX(ZREFL(JROF,JLEV),1.E-18_JPRB) ) ! Z_e in dBZ
    ENDDO
  ENDDO
ELSE
  PSIMRFCDB(:,:)=ZREFL(:,:)
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPPRS0D',1,ZHOOK_HANDLE)
END SUBROUTINE  GPPRS0D
