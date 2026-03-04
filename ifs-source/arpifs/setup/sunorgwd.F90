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

SUBROUTINE SUNORGWD(YDSTA,YDDIMV,YDNORGWD,KULOUT)

!**** *SUNORGWD*

!     Purpose.
!     --------
!           Initialize YOMNORGWD, THE MODULE THAT CONTROLS THE NON-OROGRAPHIC GW PARAMETRIZATION

!**   Interface.
!     ----------
!        CALLED FROM *SUPHMF*

!     Author.
!     -------
!        David Saint-Martin - Jan 2012

!     Modifications.
!     --------------
!        D. St-Martin - Sep 2018 : adapted for CY46
!     ------------------------------------------------------------------

USE PARKIND1,  ONLY : JPIM, JPRB
USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMNORGWD, ONLY : TNORGWD
USE YOMLUN   , ONLY : NULNAM
USE YOMDIMV  , ONLY : TDIMV
USE YOMSTA   , ONLY : TSTA

IMPLICIT NONE

TYPE(TSTA), INTENT(IN)               :: YDSTA
TYPE(TDIMV), INTENT(IN)              :: YDDIMV
TYPE(TNORGWD), TARGET, INTENT(INOUT) :: YDNORGWD
INTEGER(KIND=JPIM), INTENT(IN)       :: KULOUT

! internal
INTEGER(KIND=JPIM) :: JK
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

CHARACTER(LEN=4), POINTER   :: NORGWD_SCHEME 
REAL(KIND=JPRB), POINTER    :: NORGWD_PRMAX
REAL(KIND=JPRB), POINTER    :: NORGWD_DZ
REAL(KIND=JPRB), POINTER    :: NORGWD_PTROPO
REAL(KIND=JPRB), POINTER    :: NORGWD_RUWMAX
REAL(KIND=JPRB), POINTER    :: NORGWD_SAT
REAL(KIND=JPRB), POINTER    :: NORGWD_RDISS
REAL(KIND=JPRB), POINTER    :: NORGWD_DELTAT
REAL(KIND=JPRB), POINTER    :: NORGWD_KMIN
REAL(KIND=JPRB), POINTER    :: NORGWD_KMAX
REAL(KIND=JPRB), POINTER    :: NORGWD_CMIN
REAL(KIND=JPRB), POINTER    :: NORGWD_CMAX
REAL(KIND=JPRB), POINTER    :: NORGWD_PLAUNCH
REAL(KIND=JPRB), POINTER    :: NORGWD_PNOVERDIF
REAL(KIND=JPRB), POINTER    :: NORGWD_DZFRON
REAL(KIND=JPRB), POINTER    :: NORGWD_GFRON
REAL(KIND=JPRB), POINTER    :: NORGWD_GB

#include "posnam.intfb.h"
#include "namnorgwd.nam.h"

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNORGWD',0,ZHOOK_HANDLE)

!Associate for variables not in the include namelists nor allocated in the routine 
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, STPRE=>YDSTA%STPRE, &
 & NORGWD_NLAUNCH => YDNORGWD%NORGWD_NLAUNCH, &
 & NORGWD_NNOVERDIF => YDNORGWD%NORGWD_NNOVERDIF, &
 & NORGWD_NTROPO => YDNORGWD%NORGWD_NTROPO)

! Associate pointers for variables in the include namelists or allocated in the routine
NORGWD_SCHEME => YDNORGWD%NORGWD_SCHEME
NORGWD_RUWMAX => YDNORGWD%NORGWD_RUWMAX
NORGWD_SAT => YDNORGWD%NORGWD_SAT
NORGWD_RDISS => YDNORGWD%NORGWD_RDISS
NORGWD_DELTAT => YDNORGWD%NORGWD_DELTAT
NORGWD_KMIN => YDNORGWD%NORGWD_KMIN
NORGWD_KMAX => YDNORGWD%NORGWD_KMAX
NORGWD_CMIN => YDNORGWD%NORGWD_CMIN
NORGWD_CMAX => YDNORGWD%NORGWD_CMAX
NORGWD_PLAUNCH => YDNORGWD%NORGWD_PLAUNCH
NORGWD_PNOVERDIF => YDNORGWD%NORGWD_PNOVERDIF
NORGWD_PRMAX => YDNORGWD%NORGWD_PRMAX
NORGWD_DZ => YDNORGWD%NORGWD_DZ
NORGWD_PTROPO => YDNORGWD%NORGWD_PTROPO
NORGWD_GB => YDNORGWD%NORGWD_GB
NORGWD_GFRON => YDNORGWD%NORGWD_GFRON
NORGWD_DZFRON => YDNORGWD%NORGWD_DZFRON

!-------------------------------------------------------------------------
!*       1. Set default values
!-------------------------------------------------------------------------

NORGWD_SCHEME = 'RAND'

! Characteristic depth of the source (front and jets source)
NORGWD_DZFRON = 1000.0

 ! Parameter G_0 (~1) that controls the amplitude of the EP flux emitted by fronts and jets
NORGWD_GFRON = 0.60

! Parameter that controls the amplitude of the EP flux emitted by 'background' sources
NORGWD_GB = 0.00

! maximum of rain for which our theory applies (in kg/m^2/s)
NORGWD_PRMAX = 20. / 24. / 3600.

! Characteristic depth of the source
NORGWD_DZ = 1000.

! Tropopasue height for GW spectrum in Pa
NORGWD_PTROPO = 20000.0_JPRB

! Max EP-Flux at Launch altitude
NORGWD_RUWMAX = 0.02_JPRB

! Saturation parameter
NORGWD_SAT = 0.85_JPRB

! Dissipation coefficient
NORGWD_RDISS = 40.0_JPRB

! Time scale of the life cycle of the waves parameterized
NORGWD_DELTAT = 24.*3600.

! Min horizontal wavenumbers
NORGWD_KMIN = 5.E-5

! Max horizontal wavenumbers
NORGWD_KMAX = 3.E-4

! Min absolute ph. vel.
NORGWD_CMIN = 1.

! Max absolute ph. vel.
NORGWD_CMAX = 30.

! Launch height of GW spectrum in Pa
NORGWD_PLAUNCH = 80000.0_JPRB

! Bottom height for no vertical diffusion (Pa)
NORGWD_PNOVERDIF = 10000.0_JPRB

!-------------------------------------------------------------------------
!*       2. Modify default values.
!-------------------------------------------------------------------------
CALL POSNAM(NULNAM,'NAMNORGWD')
READ(NULNAM,NAMNORGWD)

!*  COMPUTE MODEL LEVEL LAUNCH HEIGHT OF GW SPECTRUM
!*  ------------------------------------------------
NORGWD_NLAUNCH = NFLEVG
DO JK = NFLEVG, 1, -1
  IF (STPRE(JK) > NORGWD_PLAUNCH) NORGWD_NLAUNCH = JK
ENDDO
NORGWD_NLAUNCH = NFLEVG - NORGWD_NLAUNCH

NORGWD_NTROPO = NFLEVG
DO JK = NFLEVG, 1, -1
  IF (STPRE(JK) > NORGWD_PTROPO) NORGWD_NTROPO = JK
ENDDO
NORGWD_NTROPO = NFLEVG - NORGWD_NTROPO

NORGWD_NNOVERDIF = NFLEVG
DO JK = NFLEVG, 1, -1
  IF (STPRE(JK) > NORGWD_PNOVERDIF) NORGWD_NNOVERDIF = JK
ENDDO

!-------------------------------------------------------------------------
!*       3. Print final values.
!-------------------------------------------------------------------------
WRITE(UNIT=KULOUT,FMT='('' MODULE YOMNORGWD '')')
WRITE(UNIT=KULOUT,FMT='('' TYPE OF SCHEME                       = '',A)')     NORGWD_SCHEME
WRITE(UNIT=KULOUT,FMT='('' MAX RAIN PRMAX                       = '',E11.4)') NORGWD_PRMAX
WRITE(UNIT=KULOUT,FMT='('' DEPTH OF THE SOURCE                  = '',F11.4)') NORGWD_DZ
WRITE(UNIT=KULOUT,FMT='('' TROPOPAUSE HEIGHT (Pa)               = '',F11.4)') NORGWD_PTROPO
WRITE(UNIT=KULOUT,FMT='('' NORGWD_NTROPO (WITH 1=BOT,KLEV=TOP)  = '',I3)')    NORGWD_NTROPO
WRITE(UNIT=KULOUT,FMT='('' LAUNCH HEIGHT (Pa)                   = '',F11.4)') NORGWD_PLAUNCH
WRITE(UNIT=KULOUT,FMT='('' NORGWD_NLAUNCH (WITH 1=BOT,KLEV=TOP) = '',I3)')    NORGWD_NLAUNCH
WRITE(UNIT=KULOUT,FMT='('' TOTAL LAUNCH MOMENTUM FLUX           = '',E11.4)') NORGWD_RUWMAX
WRITE(UNIT=KULOUT,FMT='('' SATURATION PARAMETER                 = '',F11.4)') NORGWD_SAT
WRITE(UNIT=KULOUT,FMT='('' DISSIPATION COEFFICIENT              = '',F11.4)') NORGWD_RDISS
WRITE(UNIT=KULOUT,FMT='('' TIMESCALE OF THE GW LIFE CYCLE       = '',F11.4)') NORGWD_DELTAT
WRITE(UNIT=KULOUT,FMT='('' MIN HORIZ. WAVENUMBERS               = '',F11.4)') NORGWD_KMIN
WRITE(UNIT=KULOUT,FMT='('' MAX HORIZ. WAVENUMBERS               = '',F11.4)') NORGWD_KMAX
WRITE(UNIT=KULOUT,FMT='('' MIN ABS. PHASE SPEED                 = '',F11.4)') NORGWD_CMIN
WRITE(UNIT=KULOUT,FMT='('' MAX ABS. PHASE SPEED                 = '',F11.4)') NORGWD_CMAX
WRITE(UNIT=KULOUT,FMT='('' BOTTOM HEIGHT FOR NO VER DIFF (Pa)   = '',F11.4)') NORGWD_PNOVERDIF
WRITE(UNIT=KULOUT,FMT='('' BOTTOM HEIGHT FOR NO VER DIFF        = '',I3)')    NORGWD_NNOVERDIF
WRITE(UNIT=KULOUT,FMT='('' DEPTH OF THE SOURCE (FRONTS)         = '',F11.4)') NORGWD_DZFRON
WRITE(UNIT=KULOUT,FMT='('' AMPL. LAUNCHED FRONT FLUX            = '',F11.4)') NORGWD_GFRON
WRITE(UNIT=KULOUT,FMT='('' AMPL. LAUNCHED BACKGROUND FLUX       = '',F11.4)') NORGWD_GB

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNORGWD',1,ZHOOK_HANDLE)

END SUBROUTINE SUNORGWD
