! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUCUMF(YDSTA,YDDIMV,YDECUMF)

!     THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR MASSFLUX SCHEME

!          M.TIEDTKE         E.C.M.W.F.    2/89

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *INIPHY*

!          MODIFICATIONS
!          -------------
!          P. Bechtold 2003-2013       Cleaning and revision of entrainment rates
!                                      options implicit, tracers, perturb, stand atmos
!                                      adding scaling factors for different planet
!                                      (modified gravity)
!                                      add options for diurnal cycle over land
!          P. Lopez, ECMWF (Oct 2007)  Put reading of NAMCUMF back in.
!          R. Forbes, May 2008         Changed factor in RTAUMEL from 
!                                      1.5 to 0.66
!          N. Semane+P.Bechtold     04-10-2012 Add RCORIOI/RPLRG/RPLDARE/RHOUR/RCVRFACTOR for small planet
!          T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!-----------------------------------------------------------------------

USE YOMSTA    , ONLY : TSTA
USE YOMDIMV   , ONLY : TDIMV
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RTT
USE YOMLUN    , ONLY : NULOUT, NULNAM
USE YOECUMF   , ONLY : TECUMF
USE YOETHF    , ONLY : RTWAT, RTICECU, RTWAT_RTICECU_R
USE YOMDYNCORE, ONLY : RPLRG, RPLDARE

IMPLICIT NONE

TYPE(TSTA)        ,INTENT(IN) :: YDSTA
TYPE(TDIMV)       ,INTENT(IN) :: YDDIMV
TYPE(TECUMF)      ,INTENT(INOUT), TARGET :: YDECUMF
INTEGER(KIND=JPIM) :: JLEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB), POINTER :: RMFSOLUV, RMFSOLTQ, RMFSOLCT, RMFSOLRHS, ENTRORG, ENTSHALP,&
 & DETRPEN, RTAUA, RMFDEPS, RPRCON, ENTRDD, RDEPTHS, RMINCAPE, RCAPDCYCL, RHEBC, RCPECONS,&
 & RBASE0, RMINCIN, ENTSTPC1, ENTSTPC2, ENTSTPC3, RMFCFL, RMFLIA, RUVPER, RMFADVW, RMFADVWDD,&
 & RCAPQADV
LOGICAL, POINTER :: LMFPEN, LMFCUCA, LMFSCV, LMFDUDV, LMFDSNOW, LMFWETB, LMFGLAC, LSCVLIQ, LMFENTHCONS
INTEGER(KIND=JPIM), POINTER :: NJKT7

#include "namcumf.nam.h"
#include "posnam.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUCUMF',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & LMFDD=>YDECUMF%LMFDD, LMFMID=>YDECUMF%LMFMID, &
 & LMFPROFP=>YDECUMF%LMFPROFP, LMFSMOOTH=>YDECUMF%LMFSMOOTH, &
 & LMFUVDIS=>YDECUMF%LMFUVDIS, LMFWSTAR=>YDECUMF%LMFWSTAR, NJKT1=>YDECUMF%NJKT1, &
 & NJKT2=>YDECUMF%NJKT2, NJKT3=>YDECUMF%NJKT3, NJKT4=>YDECUMF%NJKT4, &
 & NJKT5=>YDECUMF%NJKT5, NJKT6=>YDECUMF%NJKT6, NJKT8=>YDECUMF%NJKT8,&
 & RCUCOV=>YDECUMF%RCUCOV, RCVRFACTOR=>YDECUMF%RCVRFACTOR, &
 & RMFCMIN=>YDECUMF%RMFCMIN, RTAUMEL=>YDECUMF%RTAUMEL, &
 & STPRE=>YDSTA%STPRE)

! Associate pointers for variables in namelist
NJKT7     => YDECUMF%NJKT7
ENTSTPC1  => YDECUMF%ENTSTPC1
ENTSTPC2  => YDECUMF%ENTSTPC2
ENTSTPC3  => YDECUMF%ENTSTPC3
RMFCFL    => YDECUMF%RMFCFL
RMFSOLUV  => YDECUMF%RMFSOLUV
RMFSOLTQ  => YDECUMF%RMFSOLTQ
RMFSOLCT  => YDECUMF%RMFSOLCT
RMFSOLRHS => YDECUMF%RMFSOLRHS
ENTRORG   => YDECUMF%ENTRORG
ENTSHALP  => YDECUMF%ENTSHALP
DETRPEN   => YDECUMF%DETRPEN
RTAUA     => YDECUMF%RTAUA
RMFDEPS   => YDECUMF%RMFDEPS
RPRCON    => YDECUMF%RPRCON
ENTRDD    => YDECUMF%ENTRDD
RDEPTHS   => YDECUMF%RDEPTHS
LMFPEN    => YDECUMF%LMFPEN
LMFSCV    => YDECUMF%LMFSCV
LMFDUDV   => YDECUMF%LMFDUDV
LMFDSNOW  => YDECUMF%LMFDSNOW
LMFWETB   => YDECUMF%LMFWETB
LSCVLIQ   => YDECUMF%LSCVLIQ
LMFGLAC   => YDECUMF%LMFGLAC
LMFCUCA   => YDECUMF%LMFCUCA
LMFENTHCONS => YDECUMF%LMFENTHCONS
RCAPDCYCL => YDECUMF%RCAPDCYCL
RCAPQADV  => YDECUMF%RCAPQADV
RMINCAPE  => YDECUMF%RMINCAPE
RHEBC     => YDECUMF%RHEBC
RBASE0    => YDECUMF%RBASE0
RMINCIN   => YDECUMF%RMINCIN
RMFLIA    => YDECUMF%RMFLIA
RUVPER    => YDECUMF%RUVPER
RMFADVW   => YDECUMF%RMFADVW
RMFADVWDD => YDECUMF%RMFADVWDD
RCPECONS  => YDECUMF%RCPECONS

!-----------------------------------------------------------------------

!     1.           SPECIFY PARAMETERS FOR MASSFLUX-SCHEME
!                  --------------------------------------
!Nota:     RPLRG is a scaling factor when gravity or scale height of planet
!          is changed (eg for small planet), but for earth =1


!     DETRPEN: AVERAGE DETRAINMENT RATE FOR PENETRATIVE CONVECTION (1/M)
!     -------

DETRPEN=0.75E-4_JPRB*RPLRG

!         NOTA:SHALLOW/DEEP ENTRAINMENT RATES ARE 
!              VERTICALLY SCALED BY FUNCTION  (qs/qsb)**3

!     ENTRORG: ENTRAINMENT FOR POSITIVELY BUOYANT DEEP/SHALLOW CONVECTION 1/(M)
!     -------
ENTRORG =1.75E-3_JPRB*RPLRG

!     ENTSHALP: SHALLOW ENTRAINMENT DEFINED AS ENTSHALP*ENTRORG
ENTSHALP=2.0_JPRB

!     ENTSTPC1,2: SHALLOW ENTRAINMENT CONSTANTS FOR TRIGGER TEST PARCEL ONLY
ENTSTPC1=0.8_JPRB
ENTSTPC2=2.E-4_JPRB
!     ENTSTPC3: ADDITIONAL ENTRAINMENT FACTOR USED ONLY FOR PBL SCHEME INVERSION HEIGHT
ENTSTPC3=3.0_JPRB

!     ENTRDD: AVERAGE ENTRAINMENT RATE FOR DOWNDRAFTS
!     ------

ENTRDD =3.0E-4_JPRB*RPLRG

!     RMFCMIN:   MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     -------

RMFCMIN=1.E-8_JPRB

!     RMFDEPS:   FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     -------

RMFDEPS=0.30_JPRB

!     RDEPTHS:   MAXIMUM ALLOWED SHALLOW CLOUD DEPTH (Pa)
!     -------

RDEPTHS=2.E4_JPRB

!     RPRCON:    COEFFICIENTS FOR DETERMINING CONVERSION FROM CLOUD WATE
!     ------

RPRCON =1.4E-3_JPRB*RPLRG

!                COEFFICIENTS FOR RAIN EVAPORATION BELOW CLOUD
!                AND MELTING
!                ---------------------------------------------
!     RCPECONS:  KESSLER COEFFICIENT
!     RCVRFACTOR:  KESSLER FACTOR FOR EVAPORATION (KESSLER,1969): CONVECTIVE RAIN FLUX (R)--> R/RCVRFACTOR
!     RCUCOV:    ASSUMED CONVECTIVE CLOUD COVER
!     RTAUMEL:   MELTING TIME SCALE
!     RHEBC:     CRITICAL RELATIVE HUMIDITY BELOW CLOUD  FOR EVAPORATION

RCUCOV=0.05_JPRB
RCPECONS=5.44E-4_JPRB/RG
RCVRFACTOR=5.09E-3_JPRB/RPLRG
RTAUMEL=5._JPRB*3.6E3_JPRB/(RPLRG*RPLDARE)*0.66_JPRB
RHEBC=0.92_JPRB !for water, for land it is  modified in cuflxn

!     MODIFY/MULTIPLY ADJUSTMENT TIME SCALE FOR CAPE CLOSURE BY RTAUA
RTAUA=1.0_JPRB

!     LOGICAL SWITCHES
!     ----------------

LMFPEN  =.TRUE.   ! deep convection
LMFSCV  =.TRUE.   ! shallow convection
LMFMID  =.TRUE.   ! mid-level convection
LMFDD   =.TRUE.   ! use downdrafts
LMFDUDV =.TRUE.   ! use convective momentum transport
LMFDSNOW=.TRUE.   ! detrain snow/rain
LMFWETB =.TRUE.   ! use wet bulb T for melting
LMFGLAC =.TRUE.   ! glaciation of precip in updraught
LSCVLIQ =.TRUE.   ! condensation with respect to water for shallow
LMFUVDIS=.TRUE.   ! use kinetic energy dissipation (addit T-tendency)
LMFENTHCONS=.FALSE.! correction for enthalpy conservation in mixed phase
LMFCUCA =.FALSE.  ! modulate cloud base mass flux with CA or other 2D field
LMFPROFP=.FALSE.  ! perturb input T,q profile


!     MASSFLUX SOLVERs FOR MOMEMTUM AND TRACERS
!     0: EXPLICIT 0-1 SEMI-IMPLICIT >=1: IMPLICIT
!     -------------------------------------------

RMFSOLUV=1.0_JPRB  ! mass flux solver for momentum
RMFSOLTQ=1.0_JPRB  ! mass flux solver for T and q 
RMFSOLCT=1.0_JPRB  ! mass flux solver for chemical tracers
RMFSOLRHS=0.0_JPRB ! include (1) or not (0) RHS model tendencies in implicit solver
RMFADVW=0.0_JPRB   ! fraction [0-1] of subsidence from convection to be done by Dynamics
RMFADVWDD=0.0_JPRB ! If RMFADVW>0 keep DDmassflux in convection (0) or include in dynamics (1)
LMFSMOOTH=.FALSE.  ! Smoothing of mass fluxes top/bottom for Tracers

LMFWSTAR=.FALSE.   ! Grant w* closure for shallow convection

RCAPDCYCL=2.0_JPRB ! 0= no CAPE diurnal cycle correction
                   ! 1=    CAPE - surface buoyancy flux
                   ! 2=    CAPE - subcloud CAPE
RCAPQADV=0.8_JPRB  ! 1=    additional total moisture advection term in closure
                   ! 0= no additional total moisture advection term in closure
RMINCAPE=0.05_JPRB ! fraction (0<=RMINCAPE) of CAPE that is always adjusted

RMFLIA=2.0_JPRB*RPLDARE    ! value of absolut mass flux limit


!     RMFCFL:     MASSFLUX MULTIPLE OF CFL STABILITY CRITERIUM
!     -------

IF (RMFSOLTQ > 0.5_JPRB) THEN
    RMFCFL=3.0_JPRB
ELSE
  RMFCFL=1.0_JPRB
ENDIF

!     UPDRAUGHT VELOCITY PERTURBATION FOR IMPLICIT (M/S)
!     --------------------------------------------------

RUVPER=0.3_JPRB


!     LIMITS FOR WHEN TO USE BITMAP IN POST PROCESSING
!     --------------------------------------------------

RBASE0=25.E3_JPRB    ! cloud base height limit
RMINCIN=-1000.0_JPRB ! minimum convective inhibition valule

!IF(LMFDSNOW) RPRCON=1.5E-3_JPRB*RPLRG

CALL POSNAM(NULNAM,'NAMCUMF')
READ(NULNAM,NAMCUMF)

IF(LMFGLAC) THEN
 ! allow for more realistic sublimation range
  RTICECU=RTT-38._JPRB
  RTWAT_RTICECU_R=1.0_JPRB/(RTWAT-RTICECU)
ENDIF

!     TOP INTEGER LEVELS (cheaper computations)
!     -----------------------------------------

NJKT1=2
NJKT2=2
NJKT8=2
NJKT3=NFLEVG-2
DO JLEV=NFLEVG,2,-1
  IF(STPRE(JLEV) > 350.E2_JPRB)NJKT1=JLEV
  IF(STPRE(JLEV) >  60.E2_JPRB)NJKT2=JLEV
  IF(STPRE(JLEV) > 950.E2_JPRB)NJKT3=JLEV
  IF(STPRE(JLEV) > 850.E2_JPRB)NJKT4=JLEV
  IF(STPRE(JLEV) > 500.E2_JPRB)NJKT5=JLEV
  IF(STPRE(JLEV) > 700.E2_JPRB)NJKT6=JLEV
  IF(STPRE(JLEV) > 100.E2_JPRB)NJKT7=JLEV
  IF(STPRE(JLEV) >   1.E2_JPRB)NJKT8=JLEV
ENDDO
NJKT3=MIN(NFLEVG-2,NJKT3)

WRITE(NULOUT,*)'SUCUMF: NJKT1=',NJKT1,' NJKT2=',NJKT2,' NJKT3=',NJKT3
WRITE(UNIT=NULOUT,FMT='('' COMMON YOECUMF '')')
WRITE(UNIT=NULOUT,FMT='('' LMFMID = '',L5 &
 & ,'' LMFDD = '',L5,'' LMFDUDV = '',L5,'' LMFPEN = '',L5 &
 & ,'' RTAUA = '',E12.5,'' s-1'')') &
 & LMFMID,LMFDD,LMFDUDV,LMFPEN,RTAUA
WRITE(UNIT=NULOUT,FMT='('' LMFWETB = '',L5 &
 & ,'' LMFGLAC = '',L5)') &
 & LMFWETB, LMFGLAC 
WRITE(UNIT=NULOUT,FMT='('' RMFSOLUV = '',E12.5 &
 & ,'' RMFSOLTQ = '',E12.5,'' RMFSOLCT = '',E12.5 &
 & ,'' RMFSOLRHS = '',E12.5,'' RMFADVW = '',E12.5)') &
 & RMFSOLUV ,RMFSOLTQ ,RMFSOLCT ,RMFSOLRHS, RMFADVW 
WRITE(UNIT=NULOUT,FMT='('' RCAPDCYCL= '',E12.5 &
 & ,'' RCAPQADV= '',E12.5)') &
 & RCAPDCYCL, RCAPQADV

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCUMF',1,ZHOOK_HANDLE)
END SUBROUTINE SUCUMF
