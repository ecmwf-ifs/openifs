! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUGWWMS(YDSTA,YDDIMV,YDEGWWMS,YDRIP,KSMAX)

! INITIALZE YOEGWWMS, THE MODULE THAT CONTROLS THE WARNER MCINTYRE GW PARAMETRIZATION

!          A.ORR            ECMWF     August 2008

!          INTERFACE
!          ---------
!          CALLED FROM *SUPHEC*

!          MODIFICATIONS
!          -------------
!          P. Bechtold, ECMWF (October 2008) Redefine computation of launch level
!                ""           (October 2011) Switch off LRFRIC if NFLEVG=137
!        N.Semane+P.Bechtold    04-10-2012 replace 3600s by RHOUR for small planet
!          T. Stockdale+P. Bechtold, ECMWF (November 2012) Variable (reduced) launch flux 
!                              for Tropics (NGAUSS=2)
!          T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!          K. Yessad (July 2014): Move some variables.
!--------------------------------------------------------------------------------

USE YOMSTA   , ONLY : TSTA
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB 
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : NCONF
!!USE YOMDYNA  , ONLY : LRFRIC
USE YOMRIP   , ONLY : TRIP
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOEGWWMS , ONLY : TEGWWMS
USE YOMCST   , ONLY : RHOUR

!--------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSTA)          ,INTENT(IN)   :: YDSTA
TYPE(TDIMV)         ,INTENT(IN)   :: YDDIMV
TYPE(TEGWWMS),TARGET,INTENT(INOUT):: YDEGWWMS
TYPE(TRIP)          ,INTENT(IN)   :: YDRIP
INTEGER(KIND=JPIM)  ,INTENT(IN)   :: KSMAX ! horizontal (spectral) resolution of host model

! internal
INTEGER(KIND=JPIM) :: JK,IL
INTEGER(KIND=JPIM),PARAMETER :: ILN=3
REAL(KIND=JPRB)    :: ZLAUNCHP(ILN)   ! launch height of gw spectrum in Pa
REAL(KIND=JPRB)    :: ZFLUXLAUN(ILN)  ! launch flux before adjustment for height
REAL(KIND=JPRB)    :: ZSCAL
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!--------------------------------------------------------------------------------

INTEGER(KIND=JPIM), POINTER :: NLAUNCHLEV, NSLOPE, NGAUSS
REAL(KIND=JPRB), POINTER :: GFLUXLAUN, GFLUXLAUNL(:), GCSTAR,&
 & GPTWO, GGAUSSA, GGAUSSB(:), GCOEFF, GTPHYGWWMS, GMSTAR_L(:)
LOGICAL, POINTER :: LOZPR, LGACALC, LGSATL, LGINDL

#include "namgwwms.nam.h"

#include "posnam.intfb.h"

!--------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGWWMS',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NLAUNCH=>YDEGWWMS%NLAUNCH, NLAUNCHL=>YDEGWWMS%NLAUNCHL, &
 & TSTEP=>YDRIP%TSTEP, &
 & STPRE=>YDSTA%STPRE)
! Associate pointers for variables in namelist
NLAUNCHLEV => YDEGWWMS%NLAUNCHLEV
NSLOPE     => YDEGWWMS%NSLOPE
NGAUSS     => YDEGWWMS%NGAUSS
GFLUXLAUN  => YDEGWWMS%GFLUXLAUN
GFLUXLAUNL => YDEGWWMS%GFLUXLAUNL
GCSTAR     => YDEGWWMS%GCSTAR
GPTWO      => YDEGWWMS%GPTWO
GGAUSSA    => YDEGWWMS%GGAUSSA
GGAUSSB    => YDEGWWMS%GGAUSSB
GCOEFF     => YDEGWWMS%GCOEFF
GTPHYGWWMS => YDEGWWMS%GTPHYGWWMS
GMSTAR_L   => YDEGWWMS%GMSTAR_L
LOZPR      => YDEGWWMS%LOZPR
LGACALC    => YDEGWWMS%LGACALC
LGSATL     => YDEGWWMS%LGSATL
LGINDL     => YDEGWWMS%LGINDL

!--------------------------------------------------------------------------------

!             Set values of parameters/switches
!              --------------------------------

!*  SPECIFICATION OF SOURCE SPECTRUM
!*  -------------------------------
NLAUNCHLEV=1   ! number of launch levels (maximum 3)

GFLUXLAUN=3.75E-3_JPRB  !total launch momentum flux in each azimuth (rho_o x F_o)
! resolution scaling for wavenumber>700 now done for octahedral grid in gwdrag_wms.F90
!ZSCAL=(1.0_JPRB-MIN(1.0_JPRB,ATAN((MAX(KSMAX,700)-700)/REAL(6000-700))))
!write(nulout,*)'GWDSETUP SCAL=',zscal
GFLUXLAUN=GFLUXLAUN !*zscal

GFLUXLAUNL(1)=GFLUXLAUN 
GFLUXLAUNL(2)=0.0_JPRB
GFLUXLAUNL(3)=0.0_JPRB

ZLAUNCHP(1)=45000.0_JPRB  ! launch level in Pa
ZLAUNCHP(2)=70000.0_JPRB
ZLAUNCHP(3)=20000.0_JPRB

GCSTAR=1.0_JPRB         !C* (see McLandress and Scinocca 2005)

GMSTAR_L(1)=2000.0_JPRB    !m* (given as length, in m)
NSLOPE=1                !s (1,0,-1 are valid values) s is the slope at small-m 
                        !end of the launch spectrum 
GPTWO=2.0_JPRB          !2*p (3 or 2 are valid values) p is the exponent of omega 
                        !in the expression of the launch spectrum

!* Extra parameters to introduce spatial variation of the launch spectrum, if LOZPR=TRUE
!* NGAUSS=1  - launch momentum flux is proportional to total precipitation  
!* NGAUSS=2  - launch momentum flux is background plus extra gaussian distribution at equator

LOZPR=.TRUE.            !If .TRUE. then variable launch momemtum flux
NGAUSS=2                

GGAUSSA=20.0_JPRB       !gaussian distribution half-width in degrees (used if NGAUSS=2)
GGAUSSB(1)=-0.25_JPRB   !relative height of gaussian distribution at equator (used if NGAUSS=2)
GGAUSSB(2)=-0.25_JPRB
GGAUSSB(3)=-0.25_JPRB

GCOEFF=2000.0_JPRB      !multiplicative factor (used if NGAUSS=1) 
LGACALC=.FALSE.         !If .TRUE. then recompute equidistant bin spacing
LGSATL=.FALSE.          !If .TRUE. then spectrum saturated up to m* at launch
LGINDL=.FALSE.          !If .TRUE. then saturation and launch spectrum independent of launch location

!*  SETUP TIME FREQUENCY CALL OF SCHEME AS FUNCTION OF MODEL RESOLUTION
!*  -------------------------------------------------------------------

GTPHYGWWMS=RHOUR
IF(TSTEP>=RHOUR) THEN
  GTPHYGWWMS=2.0_JPRB*TSTEP
ENDIF

IF (NCONF == 131.OR.NCONF == 401.OR.NCONF == 501.OR.NCONF == 601 &
 & .OR.NCONF == 801) THEN
  GTPHYGWWMS=2.0_JPRB*TSTEP
ENDIF

CALL POSNAM(NULNAM,'NAMGWWMS')
READ (NULNAM,NAMGWWMS)
GTPHYGWWMS=INT(GTPHYGWWMS/TSTEP)*TSTEP

!*  ADJUST SPECTRUM FOR MULTI_LEVEL LAUNCHES
!*  ----------------------------------------
! Adjust spectrum and amplitude such that each neutrally propogating saturated spectrum
! has same spectral peak at a standard level, and that amplitudes are given as 
! measured at this standard level. Adjustment is valid analytically, but
! discretization and truncation of the spectrum gives errors up to 1 percent.

IF(NLAUNCHLEV>1) THEN
  ZFLUXLAUN(:)=GFLUXLAUNL(:)
  IF(NSLOPE==1) THEN
    GMSTAR_L(2)=GMSTAR_L(1)*(ZLAUNCHP(1)/ZLAUNCHP(2))**0.25
    GMSTAR_L(3)=GMSTAR_L(1)*(ZLAUNCHP(1)/ZLAUNCHP(3))**0.25
    GFLUXLAUNL(2)=GFLUXLAUNL(2)*(ZLAUNCHP(2)/ZLAUNCHP(1))**0.5
    GFLUXLAUNL(3)=GFLUXLAUNL(3)*(ZLAUNCHP(3)/ZLAUNCHP(1))**0.5
  ELSEIF(NSLOPE==2) THEN
    GMSTAR_L(2)=GMSTAR_L(1)*(ZLAUNCHP(1)/ZLAUNCHP(2))**0.2
    GMSTAR_L(3)=GMSTAR_L(1)*(ZLAUNCHP(1)/ZLAUNCHP(3))**0.2
    GFLUXLAUNL(2)=GFLUXLAUNL(2)*(ZLAUNCHP(2)/ZLAUNCHP(1))**0.6
    GFLUXLAUNL(3)=GFLUXLAUNL(3)*(ZLAUNCHP(3)/ZLAUNCHP(1))**0.6
  ENDIF
ELSE
  GFLUXLAUNL(1)=GFLUXLAUN
ENDIF


!*  COMPUTE MODEL LEVEL LAUNCH HEIGHT OF GW SPECTRUM
!*  ------------------------------------------------

NLAUNCHL(:)=NFLEVG-1
DO JK=NFLEVG,2,-1
 DO IL=1,ILN
   IF(STPRE(JK) > ZLAUNCHP(IL)) NLAUNCHL(IL)=JK
 ENDDO
ENDDO
NLAUNCH=NLAUNCHL(1)

!*  SET RAYLEIGH FRICTION TO FALSE IF NFLEVG>100 AND PTOP<1hPa
!*  ----------------------------------------------------------
!!$$IF (LRFRIC.AND.NFLEVG==137) THEN
!!$$   LRFRIC=.FALSE.
!!$$   WRITE(UNIT=NULOUT,FMT='('' SUGWWMS: ATTENTION LRFRIC RESET TO '',L1,'' FOR NFLEVG = '',I4)') LRFRIC,NFLEVG
!!$$ENDIF

WRITE(UNIT=NULOUT,FMT='('' SUGWWMS: ZLAUNCHP = '',F7.1,'' GFLUXLAUN = '',F7.5,&
      & '' GCSTAR = '',F5.1,'' GTPHYGWWMS = '',F6.0)') ZLAUNCHP(1),GFLUXLAUN,GCSTAR,GTPHYGWWMS
IF (LOZPR) THEN
   WRITE(UNIT=NULOUT,FMT='('' SUGWWMS: LOZPR = '',L1,'' NGAUSS = '',I1,&
      & '' GGAUSSA = '',F7.1,'' GGAUSSB = '',F7.2)') LOZPR,NGAUSS,GGAUSSA,GGAUSSB(1)
ENDIF
IF(NLAUNCHLEV > 1) THEN
   WRITE(UNIT=NULOUT,FMT='('' SUGWWMS: Multiple launch levels NLAUNCHLEV = '',I1)') NLAUNCHLEV
   WRITE(UNIT=NULOUT,FMT='('' ZLAUNCHPL = '',3F9.1)')  ZLAUNCHP(:)
   WRITE(UNIT=NULOUT,FMT='('' Specified GFLUXLAUNL = '',3F9.5)') ZFLUXLAUN(:)
   WRITE(UNIT=NULOUT,FMT='('' Specified GMSTAR_L = '',F8.1)') GMSTAR_L(1)
   WRITE(UNIT=NULOUT,FMT='('' Adjusted GFLUXLAUNL = '',3F9.5)') GFLUXLAUNL(:)
   WRITE(UNIT=NULOUT,FMT='('' Adjusted GMSTAR_L = '',3F8.1)') GMSTAR_L(:)
   WRITE(UNIT=NULOUT,FMT='('' GGAUSSB = '',3F7.2)') GGAUSSB(:)
ENDIF
!--------------------------------------------------------------------------- 

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGWWMS',1,ZHOOK_HANDLE)
END SUBROUTINE SUGWWMS

