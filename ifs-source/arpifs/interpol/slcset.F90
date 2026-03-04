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

SUBROUTINE SLCSET(YDGEOMETRY,YDSL,KASET,KBSET,&
 & KDGLG,KDLON,KDGSAG,KDGENG,KDGSAL,KDGENL,KDGSAH,KDGENH,&
 & KDGUXL,KDLUNG,KDLUXG,KDGUNG,KDGUXG,&
 & KDSUR1,KDLSUR,KDGSUR,KGPTOT,&
 & KPTRFLOFF,KFRSTLOFF,KYFRSTACTLAT,KYLSTACTLAT,&
 & KSTA,KONL,KLOENG,KPTRFRSTLAT,KFRSTLAT,KLSTLAT,&
 & PMU,PSQM2,LDCSETONLY,YDSLREP)

!**** *SLCSET * - Routine to initialize parallel environment
!                 necessary for horizontal interpolations. 
!                 Works with a lateral band of NSLWIDE latitudes in the halo. 

!                 Called only if DM code.
!                 For distributed memory computations are DM-local.

!     Purpose.
!     --------
!     initialise communication tables and test for horizontal 
!     interpolations.

!     Additional remarks for ALADIN:
!     - For the semi-Lagrangian or the observations interpolations the
!       calculations concern the sub-domain C+I; the consequence is that
!       KGPTOT_CAP and KDGUXG are used instead of KGPTOT and KDGENG 
!       (KGPTOT_CAP < or = KGPTOT; KDGUXG < or = KDGENG).
!       GM: code was commented out so reverted to KDGENG and KGPTOT
!     - For FULL-POS interpolations the 
!       calculations concern the whole domain C+I+E
!     - This program assumes that in ALADIN: KDGSAG=KDGSAL=1, KDGENG=KDGLG,
!       KLOENG(jgl)=KDLON for all latitudes, so simplications using
!       hard coded KDLON are not made in order to make the code as common
!       as possible between ARPEGE and ALADIN.

!**   Interface.
!     ----------
!        *CALL* *SLCSET *

!        Explicit arguments :
!        --------------------
!         INPUT:
!          YDSL        - SL_STRUCT definition
!          KASET       - north-south set (1..N_REGIONS_NS)
!          KBSET       - west-east set (1..N_REGIONS(KASET))
!          KDGLG       - number of latitude rows
!          KDLON       - length of a row of latitude near equator
!          KDGSAG      - Global version of KDGSA
!          KDGENG      - Global version of KDGEN
!          KDGSAL      - Local version of KDGSA , always 1
!          KDGENL      - No. of latitude rows for which this process has grid points
!          KDGSAH      - 1-YDSL%NSLWIDE
!          KDGENH      - KDGENL+YDSL%NSLWIDE
!          KDGUXL      - local last row in C+I zone in distributed memory Aladin
!          KDLUNG      - first meridian of the area of interest in Aladin
!          KDLUXG      - last  meridian of the area of interest in Aladin
!          KDGUNG      - first row of the area of interest in Aladin
!          KDGUXG      - last  row of the area of interest in Aladin
!          KDSUR1      - over dimensioning of KDLON for technical reasons
!          KDLSUR      - KDLON+KDSUR1
!          KDGSUR      - number of additional rows at each pole for semi-lagrangian
!                        job or calculation of observation equivalents
!          KGPTOT      - Total number of grid columns on a PE
!          KSTA        - Position of first grid column for the latitudes on
!                        a processor
!          KONL        - number of grid columns for the latitudes on a processor
!          KLOENG      - Number of grid points per latitude
!          KPTRFRSTLAT - pointer to the first latitude of each a-set in
!                        KSTA and KONL arrays
!          KFRSTLAT    - first lat of each a-set in grid-point space
!          KLSTLAT     - last lat of each a-set in grid-point space
!          KPTRFLOFF   - offset for pointer to the first latitude of own a-set
!                        KSTA and KONL arrays
!          KFRSTLOFF   - offset for first lat of own a-set in grid-point space,
!                        i.e. kfrstloff=kfrstlat(KASET)-1
!          KYFRSTACTLAT- first actual lat on this PE in grid-point space,
!                        it is kfrstlat(KASET)
!          KYLSTACTLAT - last actual lat on this PE in grid-point space,
!                        it is klstlat(KASET)
!          PMU         - mu              sin(theta)
!          PSQM2       - SQRT(R1MU2)     cos(theta)
!          LDCSETONLY  - optional, if true exit before calling slrset, default false
!                        used in mkglobstab, where some minimal halo info is required

!         OUTPUT:
!          YDSL         - SL_STRUCT definition

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!      Calls SLRSET if NPROC > 1.
!      Is called by SUSC2.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.
!        Documentation about distributed memory.

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original : 95-10-01

!     Modifications.
!     --------------
!      G. Mozdzynski (Dec 2001): support radiation grid
!      A. Bogatchev (April 2002): fixed bug for IALFGLO in LELAM case
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      F. Vana  (August 2007): optimized size of SL buffer for NEC SX
!      Y. Seity (May 2009): bf for B-level parrallelisation for LAMs (ZDIST)
!      R. El Khatib 02-Jul-2009 Bugfix for Aladin
!      K. Yessad (Jun 2009): unified version of SLRSET.
!      G. Mozdzynski (Dec 2010): optimisation to reduce halo volume
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (Aug 2011): support higher order interpolation
!      M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!      G. Mozdzynski (May 2012): further cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      G. Mozdzynski & Nils Wedi (Feb 2015): support for octahedral grids
!      G. Mozdzynski (Feb 2015): improve halo debugging
!      R. El Khatib 27-Jul-2016 halo over C+I+E (cyclic) for Fullpos
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN_IFSAUX, ONLY : NULOUT
USE YOMCST_IFSAUX, ONLY : XRPI, XRA

USE EINT_MOD , ONLY : SL_STRUCT

! arp/ifs dependencies to be solved later.
USE YOMCT0   , ONLY : LALLOPR, LELAM, NCONF
USE YOMMP0   , ONLY : NPRINTLEV, NPROC, LOUTPUT, LMPDIAG, NSLPAD, LSLDEBUG, N_REGIONS_NS, N_REGIONS_EW, MYPROC
USE YOMSLREP , ONLY : TSLREP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)       :: YDGEOMETRY
TYPE(SL_STRUCT),   INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KASET 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBSET 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGLG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSAG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGENG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSAL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGENL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSAH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGENH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUNG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUXG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUNG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPTOT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUXL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUXG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDSUR1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLSUR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSUR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTRFLOFF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFRSTLOFF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KYFRSTACTLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KYLSTACTLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA(KDGSAG:KDGENG+N_REGIONS_NS-1,N_REGIONS_EW) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KONL(KDGSAG:KDGENG+N_REGIONS_NS-1,N_REGIONS_EW) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOENG(KDGSAG:KDGENG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTRFRSTLAT(N_REGIONS_NS) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFRSTLAT(N_REGIONS_NS) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLSTLAT(N_REGIONS_NS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU(KDGSAG:KDGENG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSQM2(KDGSAG:KDGENG) 
LOGICAL,INTENT(IN),OPTIONAL      :: LDCSETONLY
TYPE(TSLREP)   ,INTENT(INOUT), OPTIONAL    :: YDSLREP



!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZLSA(1:KDGENL)
REAL(KIND=JPRB) :: ZLFA(1:KDGENL)

INTEGER(KIND=JPIM) :: IALF, IALFGLO, IALS, IALSGLO,&
 & IDLSUREAST, IDLSURWEST, INSWIDE, IFF, IGLGLO,&
 & IGLLOC, IGPREM, ILF, ILON, ILS, IPERIOD,&
 & IPTF, IPTS, IPTSET, IROFINC, ISLPAD,&
 & ISS, IU, JGLGLO, JGLGLO2, JGLLOC, JGLLOC2,&
 & JLON, JROC, JROF, JLP, IPT, ILPT, IST, IEND,&
 & IGLGLOMIN, IGLGLOMAX, &
 & JSL, JLEV, IGLSLOC, IGLSGLO, ISLWM2W, ISLWM2E, IMAPLEN, IUNUSED, JIND

INTEGER(KIND=JPIM),ALLOCATABLE :: IMAP(:,:)

LOGICAL :: LLCOND, LLP, LLNOR, LLSUD, LLHALO, LLCSETONLY

REAL(KIND=JPRB) :: Z1GP, Z1GPMAX, ZANGLEW,  ZANGLEE, ZDISTW, ZDISTE,&
 & ZDLW, ZDLE, ZDLMAXW, ZDLMAXE, ZDLSW, ZDLSE, ZGPTSW, ZGPTSE,&
 & ZMAXA, ZMINA, ZISLWM2W, ZISLWM2E
REAL(KIND=JPRB) :: ZDLON, ZLAT, ZA, ZC, ZDISTKMW, ZDISTKME
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "slrset.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SLCSET',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDEGEO=>YDGEOMETRY%YREGEO)
ASSOCIATE(NSTENCILWIDE=>YDDIM%NSTENCILWIDE, &
 & NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

! Copy some variables to simplify the SL interface
YDSL%NDGLG=KDGLG
YDSL%NDLON=KDLON
YDSL%NDGSAG=KDGSAG
YDSL%NDGENG=KDGENG
YDSL%NDGSAL=KDGSAL
YDSL%NDGENL=KDGENL
YDSL%NDGSAH=KDGSAH
YDSL%NDGENH=KDGENH
YDSL%NGPTOT=KGPTOT
YDSL%NDGUXL=KDGUXL
YDSL%NDLUNG=KDLUNG
YDSL%NDLUXG=KDLUXG
YDSL%NDGUNG=KDGUNG
YDSL%NDGUXG=KDGUXG
YDSL%NDSUR1=KDSUR1
YDSL%NDLSUR=KDLSUR
YDSL%NPTRFLOFF=KPTRFLOFF
YDSL%NFRSTLOFF=KFRSTLOFF
YDSL%MYFRSTACTLAT=KYFRSTACTLAT
YDSL%MYLSTACTLAT=KYLSTACTLAT

IF(ALLOCATED(YDSL%NLOENG)) DEALLOCATE(YDSL%NLOENG)
ALLOCATE(YDSL%NLOENG(YDSL%NDGSAG:YDSL%NDGENG))
IF (LLP) WRITE(IU,9) 'YDSL%NLOENG ',SIZE(YDSL%NLOENG),SHAPE(YDSL%NLOENG)
YDSL%NLOENG(:)=KLOENG(:)

LLCSETONLY=.FALSE.
IF( PRESENT(LDCSETONLY) ) LLCSETONLY=LDCSETONLY

!*       1. COMPUTATION OF NSLSTA,NSLONL,NSLOFF.
!        ---------------------------------------
IF( YDSL%CVER=='SL' .OR. YDSL%CVER=='AD' .OR. YDSL%CVER=='FP' .OR.&
  & YDSL%CVER=='OB' .OR. YDSL%CVER=='OA' )THEN
  YDSL%NSLGROUP=1
  YDSL%NSLWIDEN=YDSL%NSLWIDE
  YDSL%NSLWIDES=YDSL%NSLWIDE
  YDSL%NSLWIDEW=YDSL%NSLWIDE
  YDSL%NSLWIDEE=YDSL%NSLWIDE
ELSEIF( YDSL%CVER=='RI' .OR. YDSL%CVER=='RO' .OR. YDSL%CVER=='PI' .OR. YDSL%CVER=='PO')THEN
  YDSL%NSLGROUP=2
  YDSL%NSLWIDE=MAX(YDSL%NSLWIDEN,YDSL%NSLWIDES,YDSL%NSLWIDEW,YDSL%NSLWIDEE)
ELSE
  CALL ABOR1(' SLCSET: YDSL%CVER UNKNOWN')
ENDIF

ALLOCATE(YDSL%NSLSTA(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
IF (LLP) WRITE(IU,9) 'YDSL%NSLSTA ',SIZE(YDSL%NSLSTA),SHAPE(YDSL%NSLSTA)

ALLOCATE(YDSL%NSLONL(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
IF (LLP) WRITE(IU,9) 'YDSL%NSLONL ',SIZE(YDSL%NSLONL),SHAPE(YDSL%NSLONL)

ALLOCATE(YDSL%NSLOFF(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
IF (LLP) WRITE(IU,9) 'YDSL%NSLOFF ',SIZE(YDSL%NSLOFF),SHAPE(YDSL%NSLOFF)

IF (LSLDEBUG) THEN
  WRITE(NULOUT,'("SLCSET CALLED: CVER=",A2, " NSLWIDE=",I2)') YDSL%CVER, YDSL%NSLWIDE

  ALLOCATE(YDSL%LCOMPLAT(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
  IF (LLP) WRITE(IU,9) 'YDSL%LCOMPLAT ',SIZE(YDSL%LCOMPLAT),SHAPE(YDSL%LCOMPLAT)
  YDSL%LCOMPLAT(:)=.TRUE.
  
  ALLOCATE(YDSL%NLATGLO(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
  IF (LLP) WRITE(IU,9) 'YDSL%NLATGLO ',SIZE(YDSL%NLATGLO),SHAPE(YDSL%NLATGLO)
  YDSL%NLATGLO(:)=999999

  ALLOCATE(YDSL%DIST1GP(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
  IF (LLP) WRITE(IU,9) 'YDSL%DIST1GP ',SIZE(YDSL%DIST1GP),SHAPE(YDSL%DIST1GP)
  YDSL%DIST1GP(:)=999999.0_JPRB

  ALLOCATE(YDSL%NSLPTSWEST(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
  IF (LLP) WRITE(IU,9) 'YDSL%NSLPTSWEST',SIZE(YDSL%NSLPTSWEST),SHAPE(YDSL%NSLPTSWEST)
  YDSL%NSLPTSWEST(:)=0

  ALLOCATE(YDSL%NSLPTSEAST(YDSL%NDGSAL-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
  IF (LLP) WRITE(IU,9) 'YDSL%NSLPTSEAST',SIZE(YDSL%NSLPTSEAST),SHAPE(YDSL%NSLPTSEAST)
  YDSL%NSLPTSEAST(:)=0
ENDIF

IF( YDSL%CVER=='SL' .OR. YDSL%CVER=='AD' )THEN
  ALLOCATE(IMAP(4,YDSL%NDGLG))
  IMAPLEN=0
ENDIF

IF( YDSL%NSLPAD == 0 .AND. NSLPAD /= 0 )THEN
! only set YDSL%NSLPAD if it is not already set, useful for debugging halos
  YDSL%NSLPAD=NSLPAD
ENDIF

!     IDLSURWEST: number of extra points at start of complete latitude.
!     IDLSUREAST: number of extra points at end of complete latitude.
!     INSWIDE   : maximum north/south stencil width
!     IPERIOD   : extra points at start/end of complete latitude to ensure
!                 interpolation can be performed correctly
IDLSURWEST=NSTENCILWIDE-1
IDLSUREAST=NSTENCILWIDE
INSWIDE = NSTENCILWIDE
IF (YDSL%NDLON == 1) THEN
  IPERIOD=YDSL%NDSUR1
ELSE
  IPERIOD=IDLSURWEST+IDLSUREAST
ENDIF
IF (IPERIOD /= (YDSL%NDLSUR-YDSL%NDLON)) THEN
  WRITE(0,'("SLCSET: ",A," IPERIOD=",I2," YDSL%NDLSUR=",I5," YDSL%NDLON=",I5)')&
    & YDSL%CVER,IPERIOD,YDSL%NDLSUR,YDSL%NDLON
  CALL ABOR1(' SLCSET: IPERIOD AND YDSL%NDLSUR ARE NOT CONSISTENT')
ENDIF

DO JGLLOC=YDSL%NDGSAL-YDSL%NSLWIDEN,YDSL%NDGENL+YDSL%NSLWIDES
  YDSL%NSLSTA(JGLLOC)=-999999999
  YDSL%NSLONL(JGLLOC)=-999999999
  YDSL%NSLOFF(JGLLOC)=-999999999
  IF (LSLDEBUG) YDSL%LCOMPLAT(JGLLOC)=.TRUE.
ENDDO

DO JGLGLO=KFRSTLAT(KASET),KLSTLAT(KASET)
  IF( YDSL%NSLGROUP==1 )THEN
    IGLLOC=MAX(-YDSL%NSLWIDE,JGLGLO-YDSL%NFRSTLOFF)
  ELSE
    IGLLOC=JGLGLO-YDSL%NFRSTLOFF
  ENDIF
  IF (KONL(YDSL%NPTRFLOFF+IGLLOC,KBSET) > 0) THEN
    ISS=KSTA(YDSL%NPTRFLOFF+IGLLOC,KBSET)
    IFF=ISS+KONL(YDSL%NPTRFLOFF+IGLLOC,KBSET)-1
    ZLSA(IGLLOC)=REAL(ISS-1,JPRB)/REAL(YDSL%NLOENG(JGLGLO),JPRB)
    ZLFA(IGLLOC)=REAL(IFF-1,JPRB)/REAL(YDSL%NLOENG(JGLGLO),JPRB)
  ENDIF
ENDDO

IPTSET=0

ISLWM2W=YDSL%NSLWIDEW-KDGSUR+1
ISLWM2E=YDSL%NSLWIDEE-KDGSUR+1
ZISLWM2W=REAL(ISLWM2W,JPRB)
ZISLWM2E=REAL(ISLWM2E,JPRB)

!OCL  NOPREEX
IF (LELAM) THEN
  IF( YDSL%NSLGROUP==1 )THEN
    ZDISTW=ZISLWM2W/REAL(YDSL%NDLON,JPRB)
    ZDISTE=ZISLWM2E/REAL(YDSL%NDLON,JPRB)
  ELSE
    ZDISTW=(MIN(YDEGEO%EDELX,YDEGEO%EDELY)*ZISLWM2W)/XRA
    ZDISTE=(MIN(YDEGEO%EDELX,YDEGEO%EDELY)*ZISLWM2E)/XRA
  ENDIF
ENDIF

IF (LELAM .AND. YDSL%CVER=='FP') THEN
  ILS = YDSL%MYFRSTACTLAT-YDSL%NSLWIDEN-YDSL%NFRSTLOFF
  ILF = YDSL%MYLSTACTLAT +YDSL%NSLWIDES-YDSL%NFRSTLOFF
ELSE
  ILS = MAX(YDSL%NDGSAG,YDSL%MYFRSTACTLAT-YDSL%NSLWIDEN)-YDSL%NFRSTLOFF
  ILF = MIN(YDSL%NDGENG,YDSL%MYLSTACTLAT +YDSL%NSLWIDES)-YDSL%NFRSTLOFF
ENDIF

!     Loop over my processors lats, -/+ NSLWIDE, but bound by YDSL%NDGSAG,YDSL%NDGENG.

DO JGLLOC=ILS,ILF

  IF (LELAM) THEN
    IGLGLO=1
  ELSE
    IGLGLO=JGLLOC+YDSL%NFRSTLOFF
  ENDIF
  IF (LSLDEBUG) YDSL%NLATGLO(JGLLOC)=IGLGLO

  IF( YDSL%NSLGROUP==1 )THEN
    IALS=MAX(YDSL%MYFRSTACTLAT,JGLLOC+YDSL%NFRSTLOFF-YDSL%NSLWIDEN)-YDSL%NFRSTLOFF
    IALF=MIN(YDSL%MYLSTACTLAT ,JGLLOC+YDSL%NFRSTLOFF+YDSL%NSLWIDES)-YDSL%NFRSTLOFF
  ELSE
    IALS=MAX(YDSL%MYFRSTACTLAT,JGLLOC+YDSL%NFRSTLOFF-MAX(YDSL%NSLWIDEN,YDSL%NSLWIDES))-YDSL%NFRSTLOFF
    IALF=MIN(YDSL%MYLSTACTLAT ,JGLLOC+YDSL%NFRSTLOFF+MAX(YDSL%NSLWIDEN,YDSL%NSLWIDES))-YDSL%NFRSTLOFF
  ENDIF
  IF( IALS < 1.OR. IALS > YDSL%NDGENL )THEN
    CALL ABOR1(' SLCSET: INVALID STARTING LOCAL LATITUDE IALS')
  ENDIF
  IF( IALF < 1.OR. IALF > YDSL%NDGENL )THEN
    CALL ABOR1(' SLCSET: INVALID ENDING LOCAL LATITUDE IALF')
  ENDIF

  ZMINA=1.0_JPRB
  ZMAXA=0.0_JPRB

!       For DM-local latitude JGLLOC, loop -/+ YDSL%NSLWIDE latitudes to find 
!       max and min angles for my core region.

  LLHALO=.FALSE.
  DO JGLLOC2=IALS,IALF
    IF (KONL(YDSL%NPTRFLOFF+JGLLOC2,KBSET) > 0) THEN
      ZMINA=MIN(ZMINA,ZLSA(JGLLOC2))
      ZMAXA=MAX(ZMAXA,ZLFA(JGLLOC2))
      LLHALO=.TRUE.
    ENDIF
  ENDDO

!       For DM-local latitude JGLLOC, loop -/+ YDSL%NSLWIDE latitudes to find 
!       max angular separation due to overdimension of KDGSUR relative to
!       the number of northern latitudes actually necessary to do 
!       interpolations, and max angular separation for one grid point 
!       (for interpolation grid requirements).

  IF( LLHALO )THEN

    IF (LELAM) THEN

      Z1GPMAX=1.0_JPRB/REAL(YDSL%NDLON,JPRB)
      ZDLMAXW=ZDISTW
      ZDLMAXE=ZDISTE

    ELSE

      IALSGLO = MAX(YDSL%NDGSAG,JGLLOC+YDSL%NFRSTLOFF-INSWIDE)
      IALFGLO = MIN(YDSL%NDGENG,JGLLOC+YDSL%NFRSTLOFF+INSWIDE)

      Z1GPMAX=0.0_JPRB
      ZDLMAXW=0.0_JPRB
      ZDLMAXE=0.0_JPRB

      DO JGLGLO2=IALSGLO,IALFGLO
        Z1GP=1.0_JPRB/REAL(YDSL%NLOENG(JGLGLO2),JPRB)
        Z1GPMAX=MAX(Z1GPMAX,Z1GP)
        ZDLW=ZISLWM2W/REAL(YDSL%NLOENG(JGLGLO2),JPRB)
        ZDLE=ZISLWM2E/REAL(YDSL%NLOENG(JGLGLO2),JPRB)
        ZDLMAXW=MAX(ZDLMAXW,ZDLW)
        ZDLMAXE=MAX(ZDLMAXE,ZDLE)
      ENDDO

    ENDIF

!       Add 1 grid point for west and east interpolation grid.

    ZANGLEW=ZDLMAXW+Z1GPMAX
    ZANGLEE=ZDLMAXE+Z1GPMAX
    ZGPTSW=ZANGLEW*REAL(YDSL%NLOENG(IGLGLO),JPRB)+0.5_JPRB
    ZGPTSE=ZANGLEE*REAL(YDSL%NLOENG(IGLGLO),JPRB)+0.5_JPRB
    ZANGLEW=ZGPTSW/REAL(YDSL%NLOENG(IGLGLO),JPRB)
    ZANGLEE=ZGPTSE/REAL(YDSL%NLOENG(IGLGLO),JPRB)

    ZMINA=ZMINA-ZANGLEW
    ZMAXA=ZMAXA+ZANGLEE

!       Find latitude fraction and correct to ensure enough points present 
!       for interpolation requirements.

!       note IPTS will be set to the nearest westerly grid point
!       also IPTF will be set to the nearest easterly grid point.

    ZDLSW=REAL(IDLSURWEST,JPRB)
    ZDLSE=REAL(IDLSUREAST,JPRB)
    IPTS=FLOOR  (ZMINA*REAL(YDSL%NLOENG(IGLGLO))-ZDLSW)
    IPTF=CEILING(ZMAXA*REAL(YDSL%NLOENG(IGLGLO))+ZDLSE)

!       Rotate to other side if mimic lats (for ARPEGE only).

    IF( (.NOT.LELAM) .AND. (IGLGLO < 1.OR. IGLGLO > YDSL%NDGLG) )THEN
      IF( IPTS > YDSL%NLOENG(IGLGLO)/2 )THEN
        IPTS=IPTS-YDSL%NLOENG(IGLGLO)/2
        IPTF=IPTF-YDSL%NLOENG(IGLGLO)/2
      ELSE
        IPTS=IPTS+YDSL%NLOENG(IGLGLO)/2
        IPTF=IPTF+YDSL%NLOENG(IGLGLO)/2
      ENDIF
    ENDIF

!       Set NSLSTA and NSLONL to cover halo plus core area, but never more than
!       a whole latitude. Also force whole latitudes for polar latitudes. If
!       halo covers whole latitude then modify ZANGLE and ZGPTS to indicate
!       this on output.

    ISLPAD=YDSL%NSLPAD

!       Leave ISLPAD space either side of start and end points, which
!       will be left initialised to a huge number to trap halo problems
!       should they occur in the future.

    IF (LELAM) THEN
      LLCOND=IPTF-IPTS+1 >= YDSL%NLOENG(IGLGLO)
    ELSE
      LLCOND=IPTF-IPTS+1 >= YDSL%NLOENG(IGLGLO) .OR.&
       & (IGLGLO <= 1+(YDSL%NSLWIDEN-1).OR.&
       & IGLGLO >= YDSL%NDGLG-(YDSL%NSLWIDES-1))  
    ENDIF
  
    IF (LLCOND) THEN
      YDSL%NSLSTA(JGLLOC)=1
      YDSL%NSLONL(JGLLOC)=YDSL%NLOENG(IGLGLO)
      YDSL%NSLOFF(JGLLOC)=IPTSET+ISLPAD+(NSTENCILWIDE-1)
      IF (LSLDEBUG) YDSL%LCOMPLAT(JGLLOC)=.TRUE.
      IPTSET=IPTSET+YDSL%NLOENG(IGLGLO)+IPERIOD+ISLPAD*2
      ZANGLEW=1.0_JPRB
      ZANGLEE=1.0_JPRB
      ZGPTSW=YDSL%NLOENG(IGLGLO)
      ZGPTSE=YDSL%NLOENG(IGLGLO)
    ELSE
      YDSL%NSLSTA(JGLLOC)=IPTS
      YDSL%NSLONL(JGLLOC)=IPTF-IPTS+1
      YDSL%NSLOFF(JGLLOC)=IPTSET+ISLPAD
      IF (LSLDEBUG) YDSL%LCOMPLAT(JGLLOC)=.FALSE.
      IPTSET=IPTSET+YDSL%NSLONL(JGLLOC)+ISLPAD*2
    ENDIF

  ELSE

    YDSL%NSLSTA(JGLLOC)=0
    YDSL%NSLONL(JGLLOC)=0
    YDSL%NSLOFF(JGLLOC)=0
    ZDLMAXW=0.0_JPRB
    ZDLMAXE=0.0_JPRB
    Z1GPMAX=0.0_JPRB
    ZANGLEW=0.0_JPRB
    ZANGLEE=0.0_JPRB
    ZGPTSW=0.0_JPRB
    ZGPTSE=0.0_JPRB

  ENDIF

  IF( YDSL%CVER=='SL' .OR. YDSL%CVER=='AD' )THEN
    IF( IGLGLO>=1.AND.IGLGLO<=YDSL%NDGLG )THEN
      IMAPLEN=IMAPLEN+1
      IF( IMAPLEN > YDSL%NDGLG ) CALL ABOR1('SLCSET: IMAPLEN > YDSL%NDGLG')
      IMAP(1,IMAPLEN)=IGLGLO
      IMAP(2,IMAPLEN)=YDSL%NSLSTA(JGLLOC)
      IMAP(3,IMAPLEN)=YDSL%NSLONL(JGLLOC)
      IMAP(4,IMAPLEN)=YDSL%NSLOFF(JGLLOC)
    ENDIF
  ENDIF
  
!       optional diagnostic output


!       Calculate the distance in KM from ZGPTS and this latitude.
!       The following are only for diagnostic purposes.

  IF (LSLDEBUG) THEN
    Z1GP=1.0_JPRB/REAL(YDSL%NLOENG(IGLGLO),JPRB)
    ZDLON=Z1GP*2.0_JPRB*XRPI
    ZLAT=ASIN(PMU(IGLGLO))
    ZA=COS(ZLAT)*SIN(ZDLON/2.0_JPRB)
    ZC=2.0_JPRB * ASIN( MIN(1.0_JPRB,ZA) )
    YDSL%DIST1GP(JGLLOC)=XRA * ZC /1000.0_JPRB
    ZDISTKMW=ZGPTSW * XRA * ZC /1000.0_JPRB
    ZDISTKME=ZGPTSE * XRA * ZC /1000.0_JPRB
    WRITE(NULOUT,'(" SLCSET INFO: GLOBAL LAT=",I3,&
     & " N",A2,"STA=",I5,&
     & " N",A2,"ONL=",I5,&
     & " N",A2,"OFF=",I8,&
     & " ZDLMAXW=",F10.5," ZDLMAXE=",F10.5,&
     & " Z1GPMAX=",F10.5,&
     & " ZANGLEW=",F10.5," ZANGLEE=",F10.5,&
     & " ZGPTSW=",F10.5," ZGPTSE=",F10.5,&
     & " ZDISTKMW=",F10.1," ZDISTKME=",F10.1)')&
     & IGLGLO,&
     & YDSL%CVER,YDSL%NSLSTA(JGLLOC),YDSL%CVER,YDSL%NSLONL(JGLLOC),YDSL%CVER,YDSL%NSLOFF(JGLLOC),&
     & ZDLMAXW,ZDLMAXE,Z1GPMAX,ZANGLEW,ZANGLEE,ZGPTSW,ZGPTSE,ZDISTKMW,ZDISTKME
  ENDIF

ENDDO

IF( YDSL%CVER=='SL' .OR. YDSL%CVER=='AD' )THEN
  ALLOCATE(YDSL%NSLMAP(4,IMAPLEN))
  YDSL%NSLMAP(:,1:IMAPLEN)=IMAP(:,1:IMAPLEN)
  DEALLOCATE(IMAP)
ENDIF

!     ------------------------------------------------------------------

!*       2. COMPUTATION OF NSLCORE:
!           (position of core points in interpolation buffer).
!        -----------------------------------------------------

ALLOCATE(YDSL%NSLCORE(YDSL%NGPTOT))
IF(LLP)WRITE(IU,9) 'YDSL%NSLCORE  ',SIZE(YDSL%NSLCORE  ),SHAPE(YDSL%NSLCORE)

IROFINC=0
IGLLOC=1
IGPREM=KONL(YDSL%NPTRFLOFF+IGLLOC,KBSET)

DO JROF=1,YDSL%NGPTOT
  IF (IGPREM == 0) THEN
    IGLLOC=IGLLOC+1
    IGPREM=KONL(YDSL%NPTRFLOFF+IGLLOC,KBSET)
    DO WHILE (IGPREM <= 0)
      IGLLOC=IGLLOC+1
      IGPREM=KONL(YDSL%NPTRFLOFF+IGLLOC,KBSET)
    ENDDO
    IROFINC=IROFINC+IPERIOD
  ENDIF
  IGPREM=IGPREM-1
  IROFINC=YDSL%NSLOFF(IGLLOC)+&
   & KONL(YDSL%NPTRFLOFF+IGLLOC,KBSET)-IGPREM-1+&
   & KSTA(YDSL%NPTRFLOFF+IGLLOC,KBSET)-YDSL%NSLSTA(IGLLOC)-JROF  
  YDSL%NSLCORE(JROF)=IROFINC+1+JROF
ENDDO

!     ------------------------------------------------------------------

!*       3. COMPUTATION OF NSLEXT:
!           (data structure used by routines LASCAW, FPSCAW).
!        ----------------------------------------------------

IF (LELAM .AND. YDSL%CVER=='FP') THEN
  ALLOCATE(YDSL%NSLEXT(1-YDSL%NDLON-YDSL%NSLWIDEW:YDSL%NDLON+YDSL%NDLON+YDSL%NSLWIDEE,1-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
ELSE
  ALLOCATE(YDSL%NSLEXT(1-YDSL%NDLON:YDSL%NDLON+YDSL%NDLON,1-YDSL%NSLWIDEN:YDSL%NDGENL+YDSL%NSLWIDES))
ENDIF
IF(LLP)WRITE(IU,9) 'YDSL%NSLEXT   ',SIZE(YDSL%NSLEXT   ),SHAPE(YDSL%NSLEXT )

!Commented out as this initialisation is unnecessary and can be very expensive
!DO JGLLOC= 1-YDSL%NSLWIDEN,YDSL%NDGENL+YDSL%NSLWIDES
!  DO JLON=1-YDSL%NDLON,YDSL%NDLON+YDSL%NDLON
!    YDSL%NSLEXT(JLON,JGLLOC)=-999999999
!  ENDDO
!ENDDO

IF (LELAM .AND. YDSL%CVER=='FP') THEN
  IGLGLOMIN=KFRSTLAT(KASET)-YDSL%NSLWIDEN
  IGLGLOMAX=KLSTLAT(KASET)+YDSL%NSLWIDES
ELSE
  IGLGLOMIN=MAX(YDSL%NDGSAG,KFRSTLAT(KASET)-YDSL%NSLWIDEN)
  IGLGLOMAX=MIN(YDSL%NDGENG,KLSTLAT(KASET)+YDSL%NSLWIDES)
ENDIF
DO JGLGLO=IGLGLOMIN,IGLGLOMAX
  IGLLOC=JGLGLO-KFRSTLAT(KASET)+1
  DO JLON=1-YDSL%NDLON,YDSL%NDLON+YDSL%NDLON
    ILON=JLON
    IF (LELAM) THEN
      ILON=MOD(ILON-YDSL%NSLSTA(IGLLOC)+YDSL%NDLON,YDSL%NDLON)+1
      IF (YDSL%NSLSTA(IGLLOC) == 1.AND.&
         & YDSL%NSLONL(IGLLOC) == YDSL%NDLON.AND.&
         & ILON == YDSL%NDLON) THEN  
        ILON=0
      ENDIF
    ELSE
      IF (JGLGLO < 1.OR. JGLGLO > YDSL%NDGLG) THEN
        ILON=ILON+YDSL%NLOENG(JGLGLO)/2
      ENDIF
      ILON=MOD(ILON-YDSL%NSLSTA(IGLLOC)+YDSL%NLOENG(JGLGLO),YDSL%NLOENG(JGLGLO))+1
      IF (YDSL%NSLSTA(IGLLOC) == 1.AND.&
         & YDSL%NSLONL(IGLLOC) == YDSL%NLOENG(JGLGLO).AND.&
         & ILON == YDSL%NLOENG(JGLGLO)) THEN  
        ILON=0
      ENDIF
    ENDIF
    YDSL%NSLEXT(JLON,IGLLOC)=ILON
!   IF (LOUTPUT.AND.LMPDIAG) THEN
!     WRITE(NULOUT,'(&
!      &" SLCSET: JLON,JGLGLO,IGLLOC,N",A2,"EXT(JLON,IGLLOC)=",&
!      &4(2X,I6))') YDSL%CVER,JLON,JGLGLO,IGLLOC,YDSL%NSLEXT(JLON,IGLLOC)
!   ENDIF
  ENDDO
ENDDO
IF (LOUTPUT.AND.LMPDIAG) CALL FLUSH(NULOUT)

!     ------------------------------------------------------------------

!*       4. COMPUTATION OF NASLB1:
!           THEN CALL SLRSET IF NPROC > 1.
!        ---------------------------------------------

YDSL%NASLB1_TRUE=IPTSET
!     Modify NASLB1 to avoid bad values for bank conflicts
!     choice of '27' is empirical(good for both vpp700 and vpp5000).
!     Similar tests repeated for NEC SX platform peaks for (1024,29).
!     Logically the original (512,27) is replaced by the new couple.
IF(MOD(YDSL%NASLB1_TRUE,1024) == 0) THEN
  YDSL%NASLB1=YDSL%NASLB1_TRUE+29
ELSE
  YDSL%NASLB1=((YDSL%NASLB1_TRUE-1)/1024 +1 ) * 1024 +29
ENDIF

!
! Allocate and compute LSLCORE
!

IF (LSLDEBUG .OR. YDSL%CVER=='AD' ) THEN
  ALLOCATE(YDSL%LSLCORE(YDSL%NASLB1))
  IF(LLP)WRITE(IU,9) 'YDSL%LSLCORE  ',SIZE(YDSL%LSLCORE  ),SHAPE(YDSL%LSLCORE)
  YDSL%LSLCORE(:)=.FALSE.
  DO JROF=1,YDSL%NGPTOT
    YDSL%LSLCORE(YDSL%NSLCORE(JROF))=.TRUE.
  ENDDO
ENDIF

IF( LSLDEBUG .OR. .NOT.LLCSETONLY )THEN
  ALLOCATE(YDSL%MASK_SL1(YDSL%NASLB1+NSTENCILWIDE*2))
  IF(LLP)WRITE(IU,9) 'YDSL%MASK_SL1',SIZE(YDSL%MASK_SL1),SHAPE(YDSL%MASK_SL1)
  YDSL%MASK_SL1(:)=0

  ALLOCATE(YDSL%MASK_SL2(YDSL%NASLB1+NSTENCILWIDE*2))
  IF(LLP)WRITE(IU,9) 'YDSL%MASK_SL2',SIZE(YDSL%MASK_SL2),SHAPE(YDSL%MASK_SL2)
  YDSL%MASK_SL2(:)=0

  ALLOCATE(YDSL%MASK_SLD(YDSL%NASLB1+NSTENCILWIDE*2))
  IF(LLP)WRITE(IU,9) 'YDSL%MASK_SLD',SIZE(YDSL%MASK_SLD),SHAPE(YDSL%MASK_SLD)
  YDSL%MASK_SLD(:)=0
ENDIF

IF (LSLDEBUG )THEN
  DO JGLLOC=YDSL%NDGSAL-YDSL%NSLWIDEN,YDSL%NDGENL+YDSL%NSLWIDES
    IF( YDSL%NSLOFF(JGLLOC) > 0 .AND. .NOT.YDSL%LCOMPLAT(JGLLOC) )THEN
      IUNUSED=0
      DO JIND=YDSL%NSLOFF(JGLLOC)+1,YDSL%NSLOFF(JGLLOC)+YDSL%NSLONL(JGLLOC)
        IF( YDSL%LSLCORE(JIND) )THEN
          EXIT
        ELSE
          IUNUSED=IUNUSED+1
        ENDIF
      ENDDO
      YDSL%NSLPTSWEST(JGLLOC)=IUNUSED
      IUNUSED=0
      DO JIND=YDSL%NSLOFF(JGLLOC)+YDSL%NSLONL(JGLLOC),YDSL%NSLOFF(JGLLOC)+1,-1
        IF( YDSL%LSLCORE(JIND) )THEN
          EXIT
        ELSE
          IUNUSED=IUNUSED+1
        ENDIF
      ENDDO
      YDSL%NSLPTSEAST(JGLLOC)=IUNUSED
      WRITE(NULOUT,'("SLCSET: CVER=",A2," MYPROC=",I6," LAT=",I5," NSLPTSWEST=",I3," NSLPTSEAST=",I3)')&
       & YDSL%CVER,MYPROC,YDSL%NLATGLO(JGLLOC),YDSL%NSLPTSWEST(JGLLOC),YDSL%NSLPTSEAST(JGLLOC)
    ENDIF
  ENDDO
  IF( YDSL%CVER=='OB' .OR. YDSL%CVER=='OA' )THEN
    ! For OB and OA halos only provide mechanism for detecting use to invalid halo points
    !   in routines slint and slintad
    ! Mark MASK_SLD core points with 999999
    ! Mark MASK_SLD halo points and latitude extension points with 888888
    DO JROF=1,YDSL%NGPTOT
      YDSL%MASK_SLD(YDSL%NSLCORE(JROF))=999999
    ENDDO
    DO JGLLOC=YDSL%NDGSAL-YDSL%NSLWIDEN,YDSL%NDGENL+YDSL%NSLWIDES
      IF( YDSL%NSLOFF(JGLLOC) > 0 )THEN
        IF( YDSL%LCOMPLAT(JGLLOC) )THEN
          DO JIND=YDSL%NSLOFF(JGLLOC),YDSL%NSLOFF(JGLLOC)+YDSL%NSLONL(JGLLOC)+IPERIOD-1
            IF( .NOT.YDSL%LSLCORE(JIND) )THEN
              YDSL%MASK_SLD(JIND)=888888
            ENDIF
          ENDDO
        ELSE
          DO JIND=YDSL%NSLOFF(JGLLOC),YDSL%NSLOFF(JGLLOC)+YDSL%NSLONL(JGLLOC)-1
            IF( .NOT.YDSL%LSLCORE(JIND) )THEN
              YDSL%MASK_SLD(JIND)=888888
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF

!
! EXIT IF ONLY BASIC SLCSET HALO DATA STRUCTURES ARE REQUIRED
! i.e. when we don't need to do slcomm* communication with this SL_STRUCT object
!

IF( LLCSETONLY )THEN
  IF (LHOOK) CALL DR_HOOK('SLCSET',1,ZHOOK_HANDLE)
  RETURN
ENDIF

WRITE(NULOUT,'(" YDSL%CVER=",A2," YDSL%NASLB1=",I7," YDSL%NASLB1_TRUE=",I7)') YDSL%CVER,&
  & YDSL%NASLB1,YDSL%NASLB1_TRUE

!     Allocate NADMAP and RSASIGN (required for SL adjoint)

!     NADMAP is initialised to values so that SL buffer points are unique

!     RSASIGN is initialised to handle symmetric and antisymmetric
!     properties of fields at polar latitudes
!       RSASIGN(:,1)= 1.0  (for fields where RPARSL1 = 1.0)
!       RSASIGN(:,2)= 1.0  (for non mirror latitudes)
!                   =-1.0  (for mirror latitudes)

IF ( YDSL%CVER == 'AD'.AND.(NCONF == 131.OR.NCONF == 401.OR.NCONF == 801.OR.NCONF==601 )) THEN
  IF (.NOT.PRESENT(YDSLREP)) CALL ABOR1('SLCSET: YDSLREP NEEDED IN THIS CONFIGURATION')
  ALLOCATE(YDSLREP%NADMAP(YDSL%NASLB1*(NFLEVG+NSTENCILWIDE)))
  IF (LLP) WRITE(NULOUT,9) 'NADMAP ',SIZE(YDSLREP%NADMAP),SHAPE(YDSLREP%NADMAP)
  ALLOCATE(YDSLREP%RSASIGN(YDSL%NASLB1,2))
  IF (LLP) WRITE(NULOUT,9) 'RSASIGN',SIZE(YDSLREP%RSASIGN),SHAPE(YDSLREP%RSASIGN)
  DO JROF=1,YDSL%NASLB1*(NFLEVG+NSTENCILWIDE)
    YDSLREP%NADMAP(JROF)=JROF
  ENDDO
  DO JSL=1,YDSL%NASLB1
    YDSLREP%RSASIGN(JSL,1)=1.0_JPRB
    YDSLREP%RSASIGN(JSL,2)=1.0_JPRB
  ENDDO
!     ** Search through all latitudes within halo region.
  DO JGLLOC=ILS,ILF
    IGLGLO=JGLLOC+YDSL%NFRSTLOFF
    LLNOR=(IGLGLO >= YDSL%NDGSAG).AND.(IGLGLO <= 0)
    LLSUD=(IGLGLO >= YDSL%NDGLG+1).AND.(IGLGLO <= YDSL%NDGENG)
    IF (LLNOR.OR.LLSUD) THEN
      IF (LLNOR) IGLSGLO=1-IGLGLO
      IF (LLSUD) IGLSGLO=2*YDSL%NDGLG+1-IGLGLO
!         * Local view 
      IGLLOC=IGLGLO-YDSL%NFRSTLOFF
      IGLSLOC=IGLSGLO-YDSL%NFRSTLOFF
!         * Determine which mirror points come from this processors core region.
      DO JLP=YDSL%NSLSTA(IGLLOC),YDSL%NSLONL(IGLLOC)+YDSL%NSLSTA(IGLLOC)-1
        ILPT=MOD(JLP-1+YDSL%NLOENG(IGLGLO),YDSL%NLOENG(IGLGLO))+1
        IPT =YDSL%NSLOFF(IGLLOC)+&
         & MOD(ILPT-YDSL%NSLSTA(IGLLOC)+YDSL%NLOENG(IGLGLO),YDSL%NLOENG(IGLGLO))  
        IPTF=YDSL%NSLOFF(IGLSLOC)+&
         & MOD(ILPT-YDSL%NSLSTA(IGLSLOC)+YDSL%NLOENG(IGLSGLO),YDSL%NLOENG(IGLSGLO))  
!           Mark which mirrored latitudes points close to the poles require 
!           correction for symmetric and antisymmetric properties
        YDSLREP%RSASIGN(IPT+1,2)=-1.0_JPRB
!           Remap mirror points that are located on this processor
        IF ( (YDSL%NSLOFF(IGLSLOC) <= IPTF).AND.&
           & (IPTF <= YDSL%NSLOFF(IGLSLOC)+YDSL%NSLONL(IGLSLOC)-1) ) THEN  
          DO JLEV=1,NFLEVG+NSTENCILWIDE
            YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IPT+1)=(JLEV-1)*YDSL%NASLB1+IPTF+1
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  DO JGLLOC=ILS,ILF
    IF (YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF) == YDSL%NSLONL(JGLLOC)) THEN
      IST =YDSL%NSLOFF(JGLLOC)-(NSTENCILWIDE-1)
      IEND=YDSL%NSLOFF(JGLLOC)+YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF)
!         Remap latitudinal extensions
      DO JLEV=1,NFLEVG+NSTENCILWIDE
        IF(     NSTENCILWIDE==2 )THEN
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND+1)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +2)
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND+2)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +3)
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +1)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND  )
        ELSEIF( NSTENCILWIDE==3 )THEN
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND+1)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +3)
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND+2)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +4)
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND+3)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +5)
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +1)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND-1)
          YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IST +2)=YDSLREP%NADMAP((JLEV-1)*YDSL%NASLB1+IEND  )
        ELSE
          CALL ABOR1('SLCSET: NSTENCILWIDE NOT SUPPORTED, NADMAP')
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  DO JGLLOC=ILS,ILF
    IF (YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF) == YDSL%NSLONL(JGLLOC)) THEN
      IST =YDSL%NSLOFF(JGLLOC)-(NSTENCILWIDE-1)
      IEND=YDSL%NSLOFF(JGLLOC)+YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF)
!         Make latitudinal extensions consistent
      IF(     NSTENCILWIDE==2 )THEN
        YDSLREP%RSASIGN(IEND+1,2)=YDSLREP%RSASIGN(IST+2,2)
        YDSLREP%RSASIGN(IEND+2,2)=YDSLREP%RSASIGN(IST+3,2)
        YDSLREP%RSASIGN(IST +1,2)=YDSLREP%RSASIGN(IEND ,2)
      ELSEIF( NSTENCILWIDE==3 )THEN
        YDSLREP%RSASIGN(IEND+1,2)=YDSLREP%RSASIGN(IST+3 ,2)
        YDSLREP%RSASIGN(IEND+2,2)=YDSLREP%RSASIGN(IST+4 ,2)
        YDSLREP%RSASIGN(IEND+3,2)=YDSLREP%RSASIGN(IST+5 ,2)
        YDSLREP%RSASIGN(IST +1,2)=YDSLREP%RSASIGN(IEND-1,2)
        YDSLREP%RSASIGN(IST +2,2)=YDSLREP%RSASIGN(IEND  ,2)
      ELSE
        CALL ABOR1('SLCSET: NSTENCILWIDE NOT SUPPORTED, RSASIGN')
      ENDIF
    ENDIF
  ENDDO
! DO JROF=1,YDSL%NASLB1
!   WRITE(NULOUT,'("SLCSET: SLBUF POINT ",I10," IS REMAPPED TO ",I10)')&
!    &JROF,NADMAP(JROF)
! ENDDO
! DO JROF=1,YDSL%NASLB1
!   WRITE(NULOUT,'("SLCSET: RSASIGN(n,2) POINT ",I10," HAS SIGN ",E10.3)')&
!    &JROF,RSASIGN(JROF,2)
! ENDDO
ENDIF

IF (YDSL%NSLWIDE > 0) THEN

  IF (NPROC > 1) THEN
    CALL SLRSET(YDSL,&
     & KSTA,KONL,KPTRFRSTLAT,KFRSTLAT,KLSTLAT)
!         Check that NASLB1 is big enough to contain communication data
    DO JROC=1,YDSL%NSLPROCS
      IF(YDSL%NSLSENDPTR(JROC+1)-YDSL%NSLSENDPTR(JROC) > YDSL%NASLB1) THEN
        WRITE(NULOUT,*) ' SLCSET: YDSL%NASLB1 too small',&
         & YDSL%NSLSENDPTR(JROC+1)-YDSL%NSLSENDPTR(JROC),&
         & ' required, but the value is ',YDSL%NASLB1  
        CALL ABOR1('SLCSET: NASLB1 too small')
      ENDIF
    ENDDO
  ELSE
!         NSLSENDPOS,NSLRECVPOS,NSENDPTR,NRECVPTR,NSLCOMM are not needed
!         and are therefore zeroed
    YDSL%NSLPROCS=0
    YDSL%NSLRPT=0
    YDSL%NSLSPT=0
    IF (LOUTPUT) THEN
      WRITE(IU,'(" SLCSET: N",A2,"PROCS=",I6)') YDSL%CVER,YDSL%NSLPROCS
      WRITE(IU,'(" SLCSET: N",A2,"RPT=",I9)') YDSL%CVER,YDSL%NSLRPT
      WRITE(IU,'(" SLCSET: N",A2,"SPT=",I9)') YDSL%CVER,YDSL%NSLSPT
    ENDIF
  ENDIF

ENDIF

9 FORMAT(1X,'ARRAY ',A15,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SLCSET',1,ZHOOK_HANDLE)
END SUBROUTINE SLCSET

