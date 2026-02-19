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

SUBROUTINE SUDYNCORE
!**** *SUDYNCORE*  - Initialize YOMDYNCORE

!     Purpose.
!     --------
!           Initialize ACADEMIC SETUP

!**   Interface.
!     ----------
!        *CALL* *SUDYNCORE

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Nils Wedi *ECMWF*

!     Modifications.
!      N. Semane+P.Bechtold     04-10-2012 Add RPLRG/RPLDARE factors for small planet
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULOUT, NULNAM
USE YOMDYNCORE, ONLY : REFERENCE_LEVEL, RPRESSURE_SCALE, RPHSW_PL, &
 & RSIGMAB, RHS_KS_KA, RPHSW_EQ, RDELTA_THETA, RTHSW_0, RPHSW_EQ_PL_HALF, &
 & RHS_KS, RUSPLRADI, RHS_KA, RPHSW_D_REV, RHS_KF, RHS_KAPPA, RPHSW_D, &
 & LHELDSUAREZ, L_HS_WILLIAMSON,  LAQUA, LAPE, LPPSTEPS, &
 & NTESTCASE, MSSTSCHEME, NOISEVOR, RPLRADI, RCORIOI, RPLRG, RUSCORIOI, &
 & RUSPLRG, RPLDARE, RDELTA_T, RT00_DYN, RP00_DYN, RU00_DYN, RST0_DYN, &
 & RJACSLDIA, RLIPSLDIA, CTEST_FAMILY, RAMP

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IERR
REAL(KIND=JPRB) :: ZDAY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

#include "namdyncore.nam.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUDYNCORE',0,ZHOOK_HANDLE)

! postprocessing in steps (emulating hours for fdb)
LPPSTEPS=.FALSE.

! initial profiles
RT00_DYN=288.15_JPRB  ! T0
RP00_DYN=101315._JPRB ! P0
RU00_DYN=0._JPRB      ! U0
RST0_DYN=0._JPRB      ! DU0/DZ
NOISEVOR=0            ! 1 == initial noise in vorticity
RAMP=0._JPRB

! small planet: Earth radius r = a * RPLRADI 
!               Coriolis acceleration as RDAY = 86400. * RCORIOI
!               gravity rg=rg * RPLRG
! (automatically passed over to surf and wam)
! (needs unfortunately manual updating in trans/module/tpm_constants.F90 )
RUSPLRADI=1._JPRB
RUSCORIOI=1._JPRB
RUSPLRG=1._JPRB
RPLRADI=1.0_JPRB/RUSPLRADI
RCORIOI=1.0_JPRB/RUSCORIOI
RPLRG=1.0_JPRB/RUSPLRG
RPLDARE=1._JPRB

! Test case selection
CTEST_FAMILY='            '
NTESTCASE=-1

! aquaplanet defaults
LAQUA=.FALSE.
LAPE=.FALSE.
! default constant SST everywhere
MSSTSCHEME=9

! heldsuarez defaults
LHELDSUAREZ=.FALSE.
! pole - equator T difference
RDELTA_T=60._JPRB
! tropical heating differential
RDELTA_THETA=10._JPRB

RSIGMAB=0.7_JPRB
REFERENCE_LEVEL=1.0_JPRB-RSIGMAB
RPRESSURE_SCALE=1.E-5_JPRB
RHS_KAPPA=2.0_JPRB/7._JPRB

! additional stratosphere
L_HS_WILLIAMSON=.FALSE.
RTHSW_0=200._JPRB
! RPHSW_PL=200._JPRB
RPHSW_PL=10._JPRB
RPHSW_D = 10000._JPRB

! default thresholds for semi-Lagrangian diagnostics
RJACSLDIA=1.5_JPRB
RLIPSLDIA=5.0_JPRB

CALL POSNAM(NULNAM,'NAMDYNCORE')
READ(NULNAM,NAMDYNCORE)
RPLRADI=1.0_JPRB/RUSPLRADI
IF (RUSCORIOI /= 0.0_JPRB) THEN
   RCORIOI=1.0_JPRB/RUSCORIOI
ELSE
   RCORIOI=0.0_JPRB
ENDIF
RPLRG=1.0_JPRB/RUSPLRG

IF( LHELDSUAREZ ) THEN
  ZDAY=86400._JPRB*RCORIOI
  RHS_KF=1.0_JPRB/ZDAY
  RHS_KA=0.025_JPRB*RHS_KF
  RHS_KS=0.25_JPRB*RHS_KF
  RHS_KS_KA=RHS_KS-RHS_KA  

  RPHSW_EQ = RPHSW_D
  RPHSW_EQ_PL_HALF = 0.5_JPRB*(RPHSW_EQ-RPHSW_PL)
  RPHSW_D_REV = 1.0_JPRB/RPHSW_D
ENDIF

!  some tests

IERR=0

IF( LHELDSUAREZ .AND. (LAQUA .OR. LAPE) ) THEN
  IERR=IERR+1
  WRITE(NULOUT,'(1X,A,I3,A)') ' ! SUDYNCORE: ERROR NR ',IERR,' !!!'
  WRITE(NULOUT,*) 'LHELDSUAREZ .AND. (LAQUA or LAPE) cannot both be true!'
ENDIF

IF( MSSTSCHEME > 9 .OR. MSSTSCHEME < 0 ) THEN
  IERR=IERR+1
  WRITE(NULOUT,'(1X,A,I3,A)') ' ! SUDYNCORE: ERROR NR ',IERR,' !!!'
  WRITE(NULOUT,*) 'INVALID CHOICE FOR MSSTSCHEME: 1..9 !'
ENDIF

IF( NTESTCASE > 12 ) THEN
  IERR=IERR+1
  WRITE(NULOUT,'(1X,A,I3,A)') ' ! SUDYNCORE: ERROR NR ',IERR,' !!!'
  WRITE(NULOUT,*) 'INVALID CHOICE FOR NTESTCASE: -1..12 !'
ENDIF

WRITE(UNIT=NULOUT,FMT='('' RPLRADI     = '',E13.7)') RPLRADI
WRITE(UNIT=NULOUT,FMT='('' RCORIOI     = '',E13.7)') RCORIOI
WRITE(UNIT=NULOUT,FMT='('' LAQUA       = '',L2)') LAQUA
WRITE(UNIT=NULOUT,FMT='('' LAPE        = '',L2)') LAPE
WRITE(UNIT=NULOUT,FMT='('' LHELDSUAREZ = '',L2)') LHELDSUAREZ

WRITE(UNIT=NULOUT,FMT='('' RT00_DYN     = '',E13.7)') RT00_DYN
WRITE(UNIT=NULOUT,FMT='('' RP00_DYN     = '',E13.7)') RP00_DYN
WRITE(UNIT=NULOUT,FMT='('' RU00_DYN     = '',E13.7)') RU00_DYN
WRITE(UNIT=NULOUT,FMT='('' RST0_DYN     = '',E13.7)') RST0_DYN
WRITE(UNIT=NULOUT,FMT='('' NOISEVOR     = '',I1)')    NOISEVOR
WRITE(UNIT=NULOUT,FMT='('' RAMP         = '',E13.7)') RAMP


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUDYNCORE',1,ZHOOK_HANDLE)

END SUBROUTINE SUDYNCORE
