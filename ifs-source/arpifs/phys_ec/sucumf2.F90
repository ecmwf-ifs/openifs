! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUCUMF2(YDSTA,YDDIMV,YDRIP,YDCUMFS,YDECUMF2,KSMAX)

!     THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR 
!     THE LINEARIZED MASSFLUX SCHEME

!          AUTHOR
!          ------
!     P. Lopez, ECMWF (2006) 

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *INIPHY*

!          MODIFICATIONS
!          -------------
!     P. Lopez, ECMWF (Aug 2006) Separate routine for new linearized convection.
!     R. Forbes, May 2008        Changed factor in RTAUMEL2 from 1.5 to 0.66
!     P. Lopez, ECMWF (Feb 2013) Added options for diurnal cycle over land.
!     P. Lopez, ECMWF Jul 2013   Revision to match non-linear version.
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     P. Lopez, ECMWF Aug 2015   Revised mass flux limiter after Peter Bechtold's mods.
!     P. Lopez, ECMWF Feb 2019   Retuned shallow convection entrainment coefficients
!                                Decrease RMFCFL2 for time steps < 900 seconds.
!     P. Lopez, ECMWF Dec 2020   Added entrainment scaling ENTSTPC32 for PBL height computation.
!
!---------------------------------------------------------------------

USE YOMSTA    , ONLY : TSTA
USE YOMDIMV   , ONLY : TDIMV
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULOUT, NULNAM
USE YOMRIP    , ONLY : TRIP
USE YOECUMF2  , ONLY : TECUMF2
USE YOMCST    , ONLY : RG
USE YOMCUMFS  , ONLY : TCUMFS
USE YOMDYNCORE, ONLY : RPLRG, RPLRADI

IMPLICIT NONE

TYPE(TSTA)          , INTENT(IN)   :: YDSTA
TYPE(TDIMV)         , INTENT(IN)   :: YDDIMV
TYPE(TRIP)          , INTENT(IN)   :: YDRIP
TYPE(TCUMFS) ,TARGET, INTENT(INOUT):: YDCUMFS
TYPE(TECUMF2),TARGET, INTENT(INOUT):: YDECUMF2
INTEGER(KIND=JPIM)  , INTENT(IN)   :: KSMAX

INTEGER(KIND=JPIM) :: JLEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL, POINTER :: LECUMFS, LREGCV, LMFDUDV2, LMFDD2, LMFWETB2, LMFGLAC2, LMFCFL2_SHSTEP
REAL(KIND=JPRB), POINTER :: RMFSOLUV2, RMFSOLTQ2, RMFSOLCT2, RMFSOLRHS2, &
                          & RMFCFL2, RCAPDCYCL2, RTAU02, RMFLIA2

#include "namcumfs.nam.h"
#include "posnam.intfb.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUCUMF2',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & DETRPEN2=>YDECUMF2%DETRPEN2, ENTRDD2=>YDECUMF2%ENTRDD2, &
 & ENTRORG2=>YDECUMF2%ENTRORG2, ENTSHALP2=>YDECUMF2%ENTSHALP2, &
 & ENTSTPC12=>YDECUMF2%ENTSTPC12, ENTSTPC22=>YDECUMF2%ENTSTPC22, &
 & ENTSTPC32=>YDECUMF2%ENTSTPC32, &
 & LMFMID2=>YDECUMF2%LMFMID2, LMFUVDIS2=>YDECUMF2%LMFUVDIS2, &
 & NJKT12=>YDECUMF2%NJKT12, NJKT22=>YDECUMF2%NJKT22, NJKT32=>YDECUMF2%NJKT32, &
 & NJKT52=>YDECUMF2%NJKT52, NJKT62=>YDECUMF2%NJKT62, &
 & RCPECONS2=>YDECUMF2%RCPECONS2, RCUCOV2=>YDECUMF2%RCUCOV2, &
 & RDEPTHS2=>YDECUMF2%RDEPTHS2, RMFCMIN2=>YDECUMF2%RMFCMIN2, &
 & RMFDEPS2=>YDECUMF2%RMFDEPS2, &
 & RPRCON2=>YDECUMF2%RPRCON2, RTAU2=>YDECUMF2%RTAU2, &
 & RTAUMEL2=>YDECUMF2%RTAUMEL2, RUVPER2=>YDECUMF2%RUVPER2, &
 & STPRE=>YDSTA%STPRE, TSTEP=>YDRIP%TSTEP)
! Associate pointers for variables in namelist
! Logical switches.
LECUMFS    => YDCUMFS%LECUMFS
LREGCV     => YDCUMFS%LREGCV
LMFDUDV2   => YDECUMF2%LMFDUDV2
LMFDD2     => YDECUMF2%LMFDD2
LMFWETB2   => YDECUMF2%LMFWETB2
LMFGLAC2   => YDECUMF2%LMFGLAC2
LMFCFL2_SHSTEP => YDCUMFS%LMFCFL2_SHSTEP
! Real parameters.
RMFSOLUV2  => YDECUMF2%RMFSOLUV2
RMFSOLTQ2  => YDECUMF2%RMFSOLTQ2
RMFSOLCT2  => YDECUMF2%RMFSOLCT2
RMFSOLRHS2 => YDECUMF2%RMFSOLRHS2
RMFCFL2    => YDECUMF2%RMFCFL2
RCAPDCYCL2 => YDECUMF2%RCAPDCYCL2
RTAU02     => YDECUMF2%RTAU02
RMFLIA2    => YDECUMF2%RMFLIA2

!     1.           SPECIFY PARAMETERS FOR LINEARIZED MASSFLUX-SCHEME
!                  -------------------------------------------------

!Nota:     RPLRG is a scaling factor when gravity or scale height of planet
!          is changed (eg for small planet), but for earth =1

!     DETRPEN2: AVERAGE DETRAINMENT RATE FOR PENETRATIVE CONVECTION (1/M)

DETRPEN2=0.75E-4_JPRB*RPLRG

!     NOTA:ORGANIZED ENTRAINMENT RATES ARE MULTIPLIED BY rg in cuascn2.F90
!          AND VERTICALLY SCALED BY FUNCTION (qs/qsb)**3.

!     ENTRORG2: ORGANIZED ENTRAINMENT FOR POSITIVELY BUOYANT DEEP CONVECTION 1/(RG M)
!     --------  
ENTRORG2=1.75E-3_JPRB*RPLRG

!     ENTSTPC12,22: SHALLOW ENTRAINMENT CONSTANTS FOR TRIGGER TEST PARCEL ONLY
!     -------------  
ENTSTPC12=0.8_JPRB
ENTSTPC22=2.0E-4_JPRB

!     ENTSTPC32: ADDITIONAL ENTRAINMENT FACTOR USED ONLY FOR PBL SCHEME INVERSION HEIGHT
!     ---------- 
ENTSTPC32=3.0_JPRB

!     ENTSHALP2: SHALLOW ENTRAINMENT DEFINED AS ENTSHALP2*ENTRORG2
!     ---------

ENTSHALP2=2.0_JPRB

!     ENTRDD2: AVERAGE ENTRAINMENT RATE FOR DOWNDRAFTS
!     -------

ENTRDD2=3.0E-4_JPRB*RPLRG

!     RMFCMIN2:   MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     --------

RMFCMIN2=1.E-8_JPRB

!     RMFDEPS2:   FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     --------

RMFDEPS2=0.30_JPRB

!     RDEPTHS2:   MAXIMUM ALLOWED SHALLOW CLOUD DEPTH (Pa)
!     --------

RDEPTHS2=2.E4_JPRB

!     RPRCON2:    COEFFICIENTS FOR DETERMINING CONVERSION FROM CLOUD WATE
!     -------

RPRCON2 =1.4E-3_JPRB*RPLRG

!                COEFFICIENTS FOR RAIN EVAPORATION BELOW CLOUD
!                AND MELTING
!                ---------------------------------------------
!     RCPECONS2:  KESSLER COEFFICIENT
!     RCUCOV2:    ASSUMED CONVECTIVE CLOUD COVER
!     RTAUMEL2:   MELTING TIME SCALE

RCUCOV2=0.05_JPRB
RCPECONS2=5.44E-4_JPRB/RG
RTAUMEL2=5._JPRB*3.6E3_JPRB*0.66_JPRB

!     SET ADJUSTMENT TIME SCALE FOR CAPE CLOSURE AS A FUNCTION
!     OF MODEL RESOLUTION
!     RTAU IS 20 MINUTES FOR RESOLUTIONS HIGHER THAN TL319
!     RTAU IS 10 MINUTES FOR RESOLUTIONS HIGHER THAN TL511
!     RTAU IS 1 HOUR FOR ANY OTHER RESOLUTION
!     --------------------------------------------------------

RTAU2=3600.0_JPRB

IF (KSMAX > 320) RTAU2=1200.0_JPRB
IF (KSMAX > 512) RTAU2=600.0_JPRB

!RTAU2=1.0_JPRB+0.33_JPRB*RPLRADI*REAL(800)/REAL(KSMAX)
RTAU2=0.33_JPRB*RPLRADI
!RTAU2=MIN(3.0_JPRB,RTAU2)

!     LOGICAL SWITCHES
!     ----------------

LMFMID2  =.TRUE.   ! mid-level convection
LMFDD2   =.TRUE.   ! use downdrafts
LMFDUDV2 =.TRUE.   ! use convective momentum transport
LMFUVDIS2=.FALSE.  ! use kinetic energy dissipation (addit T-tendency)
LMFWETB2 =.TRUE.   ! use wet bulb T for melting
LMFGLAC2 =.TRUE.   ! glaciation of precip in updraught

!     UPDRAUGHT VELOCITY PERTURBATION FOR IMPLICIT (M/S)
!     --------------------------------------------------

RUVPER2=0.3_JPRB

!     MASSFLUX SOLVERs FOR MOMEMTUM AND TRACERS
!     0: EXPLICIT 0-1 SEMI-IMPLICIT >=1: IMPLICIT
!     -------------------------------------------

RMFSOLUV2=1.0_JPRB  ! mass flux solver for momentum
RMFSOLTQ2=1.0_JPRB  ! mass flux solver for T and q 
RMFSOLCT2=1.0_JPRB  ! mass flux solver for chemical tracers
RMFSOLRHS2=0.0_JPRB ! include (1) or not (0) RHS model tendencies in implicit solver

RCAPDCYCL2=2.0_JPRB ! 0= no CAPE diurnal cycle correction
                    ! 1=    CAPE - HS
                    ! 2=    CAPE - subcloud CAPE

RMFLIA2=2.0_JPRB   ! value of absolute mass flux limit

!RTAU02=1.0_JPRB/RTAU2

!     TOP INTEGER LEVELS (cheaper computations)
!     -----------------------------------------

NJKT12=2
NJKT22=2
NJKT32=NFLEVG-2
DO JLEV=NFLEVG,2,-1
  IF (STPRE(JLEV) > 350.E2_JPRB) NJKT12=JLEV
  IF (STPRE(JLEV) >  60.E2_JPRB) NJKT22=JLEV
  IF (STPRE(JLEV) > 950.E2_JPRB) NJKT32=JLEV
  IF (STPRE(JLEV) > 500.E2_JPRB) NJKT52=JLEV
  IF (STPRE(JLEV) > 700.E2_JPRB) NJKT62=JLEV
ENDDO
NJKT32=MIN(NFLEVG-2,NJKT32)

!   NEW SIMPLIFIED CONVECTION SCHEME
!   --------------------------------

LECUMFS=.FALSE.
LREGCV=.FALSE.

!     RMFCFL2: MASSFLUX MULTIPLE OF CFL STABILITY CRITERIUM
!     Default value for implicit formulation (resolution dependent)
!     -------------------------------------------------------------

! Mass flux limiter as a fraction of CFL criterion (for stability)
RMFCFL2=1.0_JPRB

!LOP Temporary switch for reduction of CFL for short time steps.
LMFCFL2_SHSTEP = .FALSE.

CALL POSNAM(NULNAM,'NAMCUMFS')
READ(NULNAM,NAMCUMFS)

IF (LMFCFL2_SHSTEP) THEN
  ! Reduce limiter for time steps shorter than 900s (for stability).
  IF (TSTEP < 900._JPRB) RMFCFL2 = RMFCFL2 * TSTEP / 900._JPRB
ENDIF

WRITE(NULOUT,*)'SUCUMF2'
WRITE(UNIT=NULOUT,FMT='('' COMMON YOMCUMFS '')')
WRITE(UNIT=NULOUT,FMT='('' LECUMFS = '',L5&
 & ,'' LREGCV = '',L5)')&
 & LECUMFS,LREGCV  
WRITE(UNIT=NULOUT,FMT='('' COMMON YOECUMF2 '')')
WRITE(UNIT=NULOUT,FMT='('' NJKT12 = '',I4&
 & ,'' NJKT22 = '',I4,'' NJKT32 = '',I4&
 & ,'' NJKT52 = '',I4,'' NJKT62 = '',I4)')&
 & NJKT12,NJKT22,NJKT32,NJKT52,NJKT62
WRITE(UNIT=NULOUT,FMT='('' LMFMID2 = '',L5&
 & ,'' LMFDD2 = '',L5,'' LMFDUDV2 = '',L5&
 & ,'' RTAU2 = '',E12.5,'' s-1'')')&
 & LMFMID2,LMFDD2,LMFDUDV2,RTAU2  
WRITE(UNIT=NULOUT,FMT='('' LMFWETB2 = '',L5&
 & ,'' LMFGLAC2 = '',L5)')&
 & LMFWETB2,LMFGLAC2 
WRITE(UNIT=NULOUT,FMT='('' RMFSOLUV2 = '',E12.5&
 & ,'' RMFSOLTQ2 = '',E12.5,'' RMFSOLCT2 = '',E12.5&
 & ,'' RMFSOLRHS2 = '',E12.5)')&
 & RMFSOLUV2 ,RMFSOLTQ2 ,RMFSOLCT2 ,RMFSOLRHS2
WRITE(UNIT=NULOUT,FMT='('' RMFCFL2 = '',E12.5)') RMFCFL2  
WRITE(UNIT=NULOUT,FMT='('' RCAPDCYCL2 = '',E12.5)') RCAPDCYCL2  
WRITE(UNIT=NULOUT,FMT='('' LMFCFL2_SHSTEP = '',L5)') LMFCFL2_SHSTEP  

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCUMF2',1,ZHOOK_HANDLE)
END SUBROUTINE SUCUMF2
