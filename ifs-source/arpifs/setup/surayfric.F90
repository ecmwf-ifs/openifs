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

!OCL NOEVAL
SUBROUTINE SURAYFRIC(YDVAB,YDDIMV,YDDYN)

!     ------------------------------------------------------------------
!**** *SURAYFRIC*   - Initialize Rayleigh friction

!     Purpose.
!     --------

!         Computes the Rayleigh friction coefficient RKRF.

!**   Interface.
!     ----------
!        *CALL* *SURAYFRIC

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!        K. Yessad (CNRM/GMAP), after some code previously in SUHDF
!         (coded at ECMWF by A. Untch).
!        This is the old part 3 of SUHDF.
!        Original : August 2005.

!     Modifications.
!     --------------
!      K. Yessad (Nov 2011): use RRFTAU.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE YOMVERT  , ONLY : TVAB, VP00
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMDYN   , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB) , INTENT(IN)    :: YDVAB
TYPE(TDIMV), INTENT(IN)    :: YDDIMV
TYPE(TDYN) , INTENT(INOUT) :: YDDYN

INTEGER(KIND=JPIM) :: JLEV
REAL(KIND=JPRB) :: ZPRESH(0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPRS, ZRFH, ZRFHH, ZRFP0, ZRFPLM, ZRFTAU, ZRFZ1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURAYFRIC',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NMAXLEVRF=>YDDYN%NMAXLEVRF, RKRF=>YDDYN%RKRF, RRFPLM=>YDDYN%RRFPLM, &
 & RRFTAU=>YDDYN%RRFTAU, RRFZ1=>YDDYN%RRFZ1)
!     ------------------------------------------------------------------

!*       1.    RAYLEIGH FRICTION
!              -----------------

! * Remark: this piece of code is fully consistent with
!   lvertfe=F, ndlnpr=0, laprxpk=T
!   where prehyd(l)=0.5*prehyd(lbar-1)+0.5*prehyd(lbar)
!   and not completely consistent with the other cases
!   for lvertfe,ndlnpr,laprxpk.

NMAXLEVRF=0
IF(YDDYN%LRFRIC) THEN
  ZRFHH=7._JPRB
  ZRFP0=1.E+5_JPRB
! Reference height for profile now set by namelist input
  ZRFZ1=RRFZ1
  ZRFH =7.7_JPRB
  ZRFTAU=RRFTAU

! Rayleigh friction only for p < RRFPLM, set by namelist
  ZRFPLM=RRFPLM

  RKRF(:)=0.0_JPRB
  IF(0.5_JPRB*YDVAB%VAH(1) < ZRFPLM) THEN
    ZPRESH(0)=YDVAB%VAH(0)+YDVAB%VBH(0)*VP00
    DO JLEV=1,NFLEVG
      ZPRESH(JLEV)=YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*VP00
      ZPRS=0.5_JPRB*(ZPRESH(JLEV)+ZPRESH(JLEV-1))
      IF(ZPRS < ZRFPLM) THEN
        NMAXLEVRF=JLEV
        RKRF(JLEV)=(1.0_JPRB-TANH((ZRFZ1-ZRFHH*LOG(ZRFP0/ZPRS))/ZRFH))/ZRFTAU
      ENDIF
    ENDDO
  ELSE
    YDDYN%LRFRIC=.FALSE.
    WRITE(UNIT=NULOUT,FMT='('' NO RAYLEIGH FRICTION APPLIED ''&
     & ,''BECAUSE THE TOP OF THE MODEL IS BELOW '',F4.2&
     & ,''hPa. LRFRIC SET TO FALSE!'')')  ZRFPLM/100._JPRB  
  ENDIF
  WRITE(NULOUT,'(A)')
  WRITE(NULOUT,'(A)') ' --- SURAYFRIC:'
  WRITE(NULOUT,'(A)') ' RAYLEIGH FRICTION: RKRF(JLEV)'
  WRITE(NULOUT,'(5E24.14)')  (RKRF(JLEV),JLEV=1,NFLEVG)
  WRITE(NULOUT,'(A)')
  WRITE(NULOUT,FMT='('' RRFTAU  = '',E20.14)') RRFTAU
  WRITE(NULOUT,FMT='('' RRFZ1   = '',E20.14)') RRFZ1
  WRITE(NULOUT,FMT='('' RRFPLM  = '',E20.14)') RRFPLM
  WRITE(NULOUT,FMT='('' NMAXLEVRF  = '',I4)') NMAXLEVRF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURAYFRIC',1,ZHOOK_HANDLE)
END SUBROUTINE SURAYFRIC
