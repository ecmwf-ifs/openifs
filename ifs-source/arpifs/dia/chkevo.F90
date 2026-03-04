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

SUBROUTINE CHKEVO(YDGEOMETRY,YDML_GCONF,KSTEP,PSTEP,YDSP)

!**** *CHKEVO*  Writes certain gridpoint values and horizontal RMS tendencies
!                to file for diagnostics. (works from global spectral arrays)

!                !!! Currently GFL fields are not done !!!

!     Purpose.
!     --------
!     Time evolution diagnostics

!**   Interface.
!     ----------
!        *CALL* *CHKEVO(KSTEP,PSTEP)*

!        Explicit arguments :
!        --------------------
!        KSTEP: current step
!        PSTEP: timestep

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------

!     Externals.
!     ----------
!       SPEREE/ESPEREE, SPEUV/ESPVDDUV : inverse spectral transforms

!     Reference.
!     ----------
!       ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!      Gabor Radnoti GMAP/MICECO
!      Original : 92-12-24

!     Modifications.
!     --------------
!      C. Fisher,L. Gaytandjieva: 01-04-02 update for aladin AL15(ntstagp)
!                 - reset nchktend to optimized value for DM
!                 - allocate tendchk here instead of suallo
!      R. El Khatib : 01-08-07 Pruning options
!      O.Spaniel    : 03-04-15 cleaning-empty interface
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.YESSAD (Jan 2004): remove the old passive scalars.
!      S. Ivatek-Sahdan : 29-03-2007 
!                 - bug in comunication if data are on diff CPU
!                 - reordering data for output, eleaning
!      P.Termonia (Jun 2008): fix a bug.
!      K.Yessad (Jan 2011): call ESPEUV instead of ESPUV
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden Sept 2016 : Removed use of SPA3, added explicit spectral argument
!      A. Deckmyn, D. Degrauwe (May 2019): fix u/v bug, allow for E-W distribution 
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN                 , ONLY : NULUSR3
USE YOMCT0                 , ONLY : CNMEXP, LELAM
USE YOMCT2                 , ONLY : NSTAR2
USE YOMOPH0                , ONLY : CNMCA
USE YOMRIP0                , ONLY : NINDAT, NSSSSS
USE YOMCHK                 , ONLY : LECHKTND, LECHKPS, NXCHK, NYCHK, NNFCHK,&
 &                                  TENDCHK, NFRQCHK, NFLDCHK, NGPCHK, NLENCHK, NCHKTEND  
USE YOMMP0                 , ONLY : MYPROC, NPRCIDS, MY_REGION_EW, NPROC, NPRGPNS, MYSETV, NPRGPEW, NPRTRV
USE YOMTAG                 , ONLY : MTAGDISTGP
USE YOMLUN                 , ONLY : NULOUT
USE MPL_MODULE
USE SPECTRAL_FIELDS_MOD    , ONLY : SPECTRAL_FIELD,ASSIGNMENT(=)

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               , INTENT(IN) :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(IN) :: YDML_GCONF
INTEGER(KIND=JPIM)           , INTENT(IN) :: KSTEP 
REAL(KIND=JPRB)              , INTENT(IN) :: PSTEP
TYPE(SPECTRAL_FIELD)         , INTENT(IN) :: YDSP

!     ------------------------------------------------------------------

CHARACTER(LEN= 4) :: CLF2, CLNSTEP
CHARACTER(LEN= 5) :: CLF1
CHARACTER(LEN= 7) :: CLF3
CHARACTER(LEN=16) :: CLFNAM,CLNOMA

INTEGER(KIND=JPIM), ALLOCATABLE :: IDATEF(:), IDOC(:)

INTEGER(KIND=JPIM) :: IXCHK(NGPCHK),IYCHK(NGPCHK),IYCHKG(NGPCHK),&
                      & IXCHKG(NGPCHK), IGPLOC(NGPCHK)
INTEGER(KIND=JPIM) :: ILONG, ILATG, ILON, ILAT, I_KPE

INTEGER(KIND=JPIM) :: IDEC, IERR, IFLD, IGPCHK, IGPCUM, IGPLOC2,&
 & IGPTOT, ILENR, ILEV, ILLCHK, ILLPCHK, ILOSP, INBARI, INBARP,&
 & IND, INDJ, INDJI, INIMES, INR, INZ, IREP, ISNDR,&
 & ITAGR, ITAGS, ITYP,&
 & J, JF, JGL, JI, JLON, JPROC, JSE

REAL(KIND=JPRB), ALLOCATABLE :: ZDIAG(:), ZDIAGL(:), ZDIAGS(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZGEOM(:), ZGEOML(:), ZGEOMS(:)

REAL(KIND=JPRB) :: ZSPEC(YDGEOMETRY%YRDIM%NSPEC2),ZSPVOR(1,YDGEOMETRY%YRDIM%NSPEC2),ZSPDIV(1,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZREEL(YDGEOMETRY%YRGEM%NGPTOT),ZREDU(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZUR(YDGEOMETRY%YRGEM%NGPTOT,1),ZVR(YDGEOMETRY%YRGEM%NGPTOT,1)
REAL(KIND=JPRB) :: ZSPMEANU(1),ZSPMEANV(1)
REAL(KIND=JPRB) :: ZSUMA, ZSUMR

LOGICAL :: LLCOP,  LLSPR
LOGICAL :: LLDONE, LLEXIST, LLWAIT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "esperee.intfb.h"
#include "espeuv.intfb.h"
#include "speree.intfb.h"
#include "speuv.intfb.h"
#include "gl2ll.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHKEVO',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDRIP=>YDML_GCONF%YRRIP,YDDIMF=>YDML_GCONF%YRDIMF)

ASSOCIATE(NDGENL=>YDDIM%NDGENL, NDGSAL=>YDDIM%NDGSAL, NDLON=>YDDIM%NDLON, &
 & NSPEC2=>YDDIM%NSPEC2, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLSUR=>YDDIMV%NFLSUR, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, NSTAGP=>YDGEM%NSTAGP, &
 & NTSTAGP=>YDGEM%NTSTAGP, &
 & MYLATS=>YDMP%MYLATS, NONL=>YDMP%NONL, NPTRFLOFF=>YDMP%NPTRFLOFF, &
 & NSTOP=>YDRIP%NSTOP, NPTRLL=>YDMP%NPTRLL)
!     ------------------------------------------------------------------

!           0. Initial set up of local arrays
!              ------------------------------

! Reset NCHKTEND
IF (LECHKTND) THEN
  NCHKTEND= NFLDCHK*NDLON*(NDGENL - NDGSAL + 1)
ENDIF

! Allocate TENDCHK
IF(.NOT.ALLOCATED(TENDCHK)) THEN
  ALLOCATE(TENDCHK(NCHKTEND))
ENDIF

IGPCHK= 0
IGPLOC2= 0
IGPCUM= 0
IXCHK(:)= 0
IYCHK(:)= 0
IXCHKG(:)=0
IYCHKG(:)=0
IGPLOC(:)=0

! Allocate local arrays
ALLOCATE(ZDIAGL(NLENCHK+1))
ALLOCATE(ZGEOML(7*NGPCHK+1))
ZDIAGL(:)= 0.0_JPRB
ZGEOML(:)= 0.0_JPRB

IDEC= MYLATS(1)-1
DO J=1,NGPCHK
  ILONG=NXCHK(J)
  ILATG=NYCHK(J)
  CALL GL2LL(YDMP,ILATG,ILONG,I_KPE,ILAT,ILON)
    
  IF (MYPROC == I_KPE) THEN
    IGPCHK=IGPCHK+1
    IXCHKG(IGPCHK)= ILONG
    IYCHKG(IGPCHK)= ILATG     
    IF (LELAM) THEN
      IGPLOC(IGPCHK) = ILON+NTSTAGP(ILAT)-1
    ELSE
      IGPLOC(IGPCHK) = ILON+NSTAGP(ILAT)-1
    ENDIF
  ENDIF
ENDDO
ZDIAGL(NLENCHK+1)=REAL(IGPCHK,JPRB)

!     ------------------------------------------------------------------

!           1. Open output file at initial time step
!              -------------------------------------

IF (KSTEP == NSTAR2) THEN

  TENDCHK(:)=0.0_JPRB

  IF (MYPROC == 1) THEN
    ALLOCATE(IDATEF(11))
    ALLOCATE(IDOC(7+NFLDCHK+2*NGPCHK))

    ! Open file
    INIMES = 1
    INBARP = NSTOP+4
    INBARI = 0
    CLF1= 'ICMSH'
    CLF2= CNMEXP(1:4)
    CLF3= 'CHKOUT2'
    WRITE(CLFNAM,'(A5,A4,A7)') CLF1,CLF2,CLF3
    CALL FAITOU(IREP,NULUSR3,.TRUE.,CLFNAM,'UNKNOWN',&
     & .TRUE.,.TRUE.,INIMES,INBARP,INBARI,CNMCA)  

    ! Date
    IDATEF(1) = NINDAT/10000
    IDATEF(2) = (NINDAT-10000*IDATEF(1))/100
    IDATEF(3) = NINDAT-10000*IDATEF(1)-100*IDATEF(2)
    IDATEF(4) = NSSSSS/3600
    IDATEF(5) = 0
    IDATEF(6) = 1
    IDATEF(7) = 0
    IDATEF(8) = 0
    IDATEF(9) = 1
    IDATEF(10)= 0
    IDATEF(11)= 0
    CALL FANDAR(IREP,NULUSR3,IDATEF)

    ! Required diagnostics
    IDOC(1)= NFRQCHK
    IDOC(2)= NFLDCHK
    IDOC(3)= NGPCHK
    IDOC(4)= NFLEVG
    IDOC(5)= NFTHER
    IDOC(6)= 0
    !!! IDOC(6)= will be set to YGFL%NUMSPFLDS             !!!
    !!! when implementing the code for GFL variables       !!!

    ! Check conversion from ln(Ps) to Ps
    ILOSP= 0
    IF (LECHKPS) ILOSP= 1
    IDOC(7)= ILOSP
    LECHKPS= LECHKPS .AND. ILOSP > 0

    ! Fields and gridpoints
    DO JF=1,NFLDCHK
      IDOC(JF+7)= NNFCHK(JF)
    ENDDO
    DO J=1,NGPCHK
      IDOC(2*J-1+7+NFLDCHK)= NXCHK(J)
      IDOC(2*J  +7+NFLDCHK)= NYCHK(J)
    ENDDO
    CLNOMA='REDOCU0000000000'
    CALL FAISAN(IREP,NULUSR3,CLNOMA,IDOC,7+NFLDCHK+2*NGPCHK)

    ! Time-step
    CLNOMA='RESTEP0000000000'
    CALL FAISAN(IREP,NULUSR3,CLNOMA,PSTEP,1)

    DEALLOCATE(IDATEF)
    DEALLOCATE(IDOC)
  ENDIF
ENDIF ! IF (KSTEP == NSTART2) THEN

! Geometry
! Gometry is needed every TSTEP to know where
! point is, but geometry is writen out just once
IGPTOT= 7*NGPCHK+1
DO J=1,IGPCHK
  ZGEOML(J         )= YDGSGEOM_NB%GELAM (INR)
  ZGEOML(J+  NGPCHK)= YDGSGEOM_NB%GELAT (INR)
  ZGEOML(J+2*NGPCHK)= YDGSGEOM_NB%GM    (INR)
  ZGEOML(J+3*NGPCHK)= YDGSGEOM_NB%GNORDL(INR)
  ZGEOML(J+4*NGPCHK)= YDGSGEOM_NB%GNORDM(INR)
	ZGEOML(J+5*NGPCHK)= REAL(IXCHKG(J),JPRB)
  ZGEOML(J+6*NGPCHK)= REAL(IYCHKG(J),JPRB)
ENDDO
ZGEOML(IGPTOT)=REAL(IGPCHK,JPRB)

! Processor communication and file writing
IGPCUM= 0
IF (MYPROC /= 1) THEN
  ! Sending information
  ITAGS= MTAGDISTGP + MYPROC
  IERR= 0
  CALL MPL_SEND(ZGEOML,KDEST=NPRCIDS(1),KTAG=ITAGS,&
   & CDSTRING='CHKEVO:')  
ELSE
  ! Gathering geo information
  ALLOCATE(ZGEOM(7*NGPCHK+1))
  ZGEOM(:)= 0.0_JPRB
  DO J=1,IGPCHK
    ZGEOM(J         )= ZGEOML(J         )
    ZGEOM(J+  NGPCHK)= ZGEOML(J+  NGPCHK)
    ZGEOM(J+2*NGPCHK)= ZGEOML(J+2*NGPCHK)
    ZGEOM(J+3*NGPCHK)= ZGEOML(J+3*NGPCHK)
    ZGEOM(J+4*NGPCHK)= ZGEOML(J+4*NGPCHK)
		ZGEOM(J+5*NGPCHK)= REAL(IXCHKG(J),JPRB)
    ZGEOM(J+6*NGPCHK)= REAL(IYCHKG(J),JPRB)
  ENDDO
  IGPCUM= IGPCHK
  DO JPROC=1,NPROC
    LLDONE= .FALSE.
    DO WHILE (.NOT. LLDONE .AND. JPROC /= 1)
      ITAGS= MTAGDISTGP + JPROC
      IERR= 0
      LLWAIT= .FALSE.
      LLEXIST= .FALSE.
      CALL MPL_PROBE(KSOURCE=JPROC,KTAG=ITAGS,LDWAIT=LLWAIT,&
       & LDFLAG=LLEXIST,CDSTRING='CHKEVO:')  
      IF (LLEXIST) THEN
        LLDONE= .TRUE.
        IERR= 0
        ILENR= 0
        ISNDR= 0
        ITAGR= 0
        ILLCHK=7*NGPCHK+1
        CALL MPL_RECV(ZGEOML(1:ILLCHK),KSOURCE=NPRCIDS(JPROC),KTAG=ITAGS,&
         & KOUNT=ILENR,CDSTRING='CHKEVO:')  
        IF (ILENR /= IGPTOT) CALL ABOR1('CHKEVO/DM : LENGTH ERROR')
        IGPLOC2= NINT(ZGEOML(ILLCHK))
        DO J=1,IGPLOC2
          ZGEOM(J+IGPCUM         )= ZGEOML(J)
          ZGEOM(J+IGPCUM+  NGPCHK)= ZGEOML(J+  NGPCHK)
          ZGEOM(J+IGPCUM+2*NGPCHK)= ZGEOML(J+2*NGPCHK)
          ZGEOM(J+IGPCUM+3*NGPCHK)= ZGEOML(J+3*NGPCHK)
          ZGEOM(J+IGPCUM+4*NGPCHK)= ZGEOML(J+4*NGPCHK)
          ZGEOM(J+IGPCUM+5*NGPCHK)= ZGEOML(J+5*NGPCHK)
          ZGEOM(J+IGPCUM+6*NGPCHK)= ZGEOML(J+6*NGPCHK)
        ENDDO
        IGPCUM= IGPCUM+IGPLOC2
      ELSE
        LLDONE= .FALSE.
      ENDIF
    ENDDO
  ENDDO
  IF (IGPCUM /= NGPCHK) CALL ABOR1('CHKEVO/DM : POINTS LOST')
ENDIF

! Reordering and writing geo information at the beginning
IF (MYPROC == 1 .AND. KSTEP == NSTAR2) THEN
  ALLOCATE(ZGEOMS(7*NGPCHK+1))
  ZGEOMS(:)= 0.0_JPRB
  ZGEOM(7*NGPCHK+1)= 0.0_JPRB
  DO J=1,NGPCHK
    DO JI=1,NGPCHK
      IF( (NXCHK(J) == NINT(ZGEOM(JI+5*NGPCHK))) .AND.&
       &(NYCHK(J) == NINT(ZGEOM(JI+6*NGPCHK))) ) THEN
        ZGEOMS(J         )=ZGEOM(JI)
        ZGEOMS(J+  NGPCHK)=ZGEOM(JI+  NGPCHK)
        ZGEOMS(J+2*NGPCHK)=ZGEOM(JI+2*NGPCHK)
        ZGEOMS(J+3*NGPCHK)=ZGEOM(JI+3*NGPCHK)
        ZGEOMS(J+4*NGPCHK)=ZGEOM(JI+4*NGPCHK)
      ENDIF
    ENDDO
  ENDDO
  ZGEOMS(5*NGPCHK+1)= 0.0_JPRB
  CLNOMA='REGEOM0000000000'
  CALL FAISAN(IREP,NULUSR3,CLNOMA,ZGEOMS,5*NGPCHK+1)
  DEALLOCATE(ZGEOMS)
ENDIF

!     ------------------------------------------------------------------

!           2. Loop over 2d fields
!              -------------------

DO JF=1,NFLDCHK

  ! * 2.1: Preliminar initialisations

  LLCOP= .FALSE.
  LLSPR= .FALSE.

  ZREEL(:)=0.0_JPRB
  ZREDU(:)=0.0_JPRB

  ! * 2.2: Copy required spectral fields

  IFLD= NNFCHK(JF)
  ILEV= MOD(IFLD,NFLEVG)
  ITYP= 1 + (IFLD-1)/NFLEVG
  IF (ILEV == 0) ILEV= NFLEVG

  ! modify ILEV for NPRTRV>1
  IF (ILEV >= NPTRLL(MYSETV) .AND. ILEV<NPTRLL(MYSETV+1)) THEN  ! note that by convention, NPTRLL(NPRTRV+1)=NFLEVG+1
    ILEV=ILEV+1-NPTRLL(MYSETV)  ! ILEV is relative in this B-set!
  ELSE
    ILEV=-1
  ENDIF

  ! vorticity:
  IF (ITYP == 1) THEN
    LLCOP= .TRUE.
    LLSPR= .TRUE.
    IF (ILEV>0) THEN
			DO JSE=1,NSPEC2
				ZSPEC(JSE)= YDSP%VOR(ILEV,JSE)
			ENDDO
    ELSE
      ZSPEC(:)=0._JPRB
    ENDIF

  ! divergence:
  ELSEIF (ITYP == 2) THEN
    LLCOP= .TRUE.
    LLSPR= .TRUE.
    IF (ILEV>0) THEN
			DO JSE=1,NSPEC2
				ZSPEC(JSE)= YDSP%DIV(ILEV,JSE)
			ENDDO
    ELSE
      ZSPEC(:)=0._JPRB
    ENDIF
		
  ! u:
  ELSEIF (ITYP == 3) THEN
    IF (LELAM) THEN
      LLCOP= .TRUE.
      IF (ILEV>0) THEN
				DO JSE=1,NSPEC2
					ZSPVOR(1,JSE)= YDSP%VOR(ILEV,JSE)
					ZSPDIV(1,JSE)= YDSP%DIV(ILEV,JSE)
				ENDDO
				ZSPMEANU(1)= YDSP%UB(ILEV)
				ZSPMEANV(1)= YDSP%VB(ILEV)
      ELSE
        ZSPVOR(1,:)=0._JPRB
        ZSPDIV(1,:)=0._JPRB
        ZSPMEANU(1)=0._JPRB
        ZSPMEANV(1)=0._JPRB
      ENDIF
      CALL ESPEUV(YDGEOMETRY,ZSPVOR,ZSPDIV,ZSPMEANU,ZSPMEANV,ZUR,ZVR,1,1,1) 
      ZREEL(1:NGPTOT)=ZUR(1:NGPTOT,1) 
      ZREDU(1:NGPTOT)=ZVR(1:NGPTOT,1)     
    ELSE
      LLCOP=.TRUE.
      CALL SPEUV(YDGEOMETRY,YDSP%VOR,YDSP%DIV,ZREEL,ZREDU,NFLSUR,1,1)
    ENDIF

  ! v:
  ELSEIF (ITYP == 4) THEN
    IF (LELAM) THEN
      LLCOP= .TRUE.
      IF (ILEV>0) THEN
				DO JSE=1,NSPEC2
					ZSPVOR(1,JSE)= YDSP%VOR(ILEV,JSE)
					ZSPDIV(1,JSE)= YDSP%DIV(ILEV,JSE)
				ENDDO
				ZSPMEANU(1)= YDSP%UB(ILEV)
				ZSPMEANV(1)= YDSP%VB(ILEV)
      ELSE
        ZSPVOR(1,:)=0._JPRB
        ZSPDIV(1,:)=0._JPRB
        ZSPMEANU(1)=0._JPRB
        ZSPMEANV(1)=0._JPRB
      ENDIF
      CALL ESPEUV(YDGEOMETRY,ZSPVOR,ZSPDIV,ZSPMEANU,ZSPMEANV,ZUR,ZVR,1,1,1)     
      ZREEL(1:NGPTOT)=ZVR(1:NGPTOT,1)
      ZREDU(1:NGPTOT)=ZUR(1:NGPTOT,1)
    ELSE
      LLCOP= .TRUE.
      CALL SPEUV(YDGEOMETRY,YDSP%VOR,YDSP%DIV,ZREDU,ZREEL,NFLSUR,1,1)
    ENDIF

  ! non-GFL thermodynamic variables:
  ELSEIF (ITYP <= (4+NFTHER)) THEN
    LLCOP= .TRUE.
    LLSPR= .TRUE.
    IF (ILEV>0) THEN    
			DO JSE=1,NSPEC2
				ZSPEC(JSE)= YDSP%HV(ILEV,JSE,ITYP-4)
			ENDDO
    ELSE
      ZSPEC(:)=0._JPRB
    ENDIF
		
  ! GFL variables:
  ! !!! not yet coded (transfer SPGFL into ZSPEC) !!!

  ! surface pressure:
  ELSEIF (IFLD == 1+(4+NFTHER)*NFLEVG) THEN
    !!! 4+NFTHER to be replaced by 4+NFTHER+YGFL%NUMSPFLDS !!!
    !!! when implementing the code for GFL variables       !!!
    LLCOP= .TRUE.
    LLSPR= .TRUE.
    DO JSE=1,NSPEC2
      ZSPEC(JSE)= YDSP%SP(JSE)
    ENDDO
  ENDIF

  ! * 2.3 Inverse spectral transform

  IF (LLSPR) THEN
    IF (LELAM) THEN
      CALL ESPEREE(YDGEOMETRY,1,1,ZSPEC,ZREEL)
    ELSE
      CALL SPEREE(YDGEOMETRY,1,1,ZSPEC,ZREEL)
    ENDIF
  ENDIF

  ! from ln(Ps) to Ps, if required
  IF (IFLD == 1+(4+NFTHER)*NFLEVG .AND. LECHKPS) THEN
    !!! 4+NFTHER to be replaced by 4+NFTHER+YGFL%NUMSPFLDS !!!
    !!! when implementing the code for GFL variables       !!!
    DO J=1,NGPTOT
      ZREEL(J)= EXP(ZREEL(J))
    ENDDO
  ENDIF

  IF  (LLCOP) THEN

  ! * 2.4 Pick up gridpoint values

    DO J=1,IGPCHK
      IND= 2*NFLDCHK + (JF-1)*NGPCHK + J
      INR=IGPLOC(J)
      ZDIAGL(IND)= ZREEL(INR)
    ENDDO

  ! * 2.5. Global tendency diagnostics
  !   (1/ ROOT-MEAN-SQUARE TENDENCY  2/ MEAN ABSOLUTE TENDENCY)

    IF (LECHKTND) THEN

      IF ( NPRGPEW>1 .OR. NPRTRV>1 ) THEN
        WRITE (NULOUT,*) "LECHKTND is bugged for NPRGPEW>1 or NPRTRV>1"
        CALL ABOR1('LECHKTND bug not yet fixed')
      ENDIF

			ZSUMR= 0.0_JPRB
      ZSUMA= 0.0_JPRB
      DO JGL=NDGSAL,NDGENL
        DO JLON=1,NONL(NPTRFLOFF+JGL,MY_REGION_EW)
          INZ= (JF-1)*(NDGENL-NDGSAL+1)*NDLON + (JGL-1)*NDLON + JLON
          IF(LELAM) THEN
            INR= JLON+NTSTAGP(JGL)-1
          ELSE  
            INR= JLON+NSTAGP(JGL)-1
          ENDIF
          IF (KSTEP > NSTAR2) THEN
            ZSUMR= ZSUMR +    (ZREEL(INR)-TENDCHK(INZ))**2
            ZSUMA= ZSUMA + ABS(ZREEL(INR)-TENDCHK(INZ))
          ENDIF
          TENDCHK(INZ)= ZREEL(INR)
        ENDDO
      ENDDO
      IND= 2*JF-1
      ZDIAGL(IND  )= ZSUMR
      ZDIAGL(IND+1)= ZSUMA
    ENDIF

  ENDIF

ENDDO

!     ------------------------------------------------------------------

!           3. Communication between processors.
!              ---------------------------------

IGPCUM= 0
IF (MYPROC /= 1) THEN
  ! Sending information (on all fields)
  ITAGS= MTAGDISTGP + (KSTEP+1)*NPROC + MYPROC
  IERR= 0
  CALL MPL_SEND(ZDIAGL,KDEST=NPRCIDS(1),KTAG=ITAGS,&
   & CDSTRING='CHKEVO:')  
ELSE
  ! Gathering information (on all fields)
  ALLOCATE(ZDIAG(NLENCHK))
  ZDIAG(:)= 0.0_JPRB
  DO JF=1,NFLDCHK
    DO J=1,IGPCHK
      IND= 2*NFLDCHK + (JF-1)*NGPCHK + J
      ZDIAG(IND)=ZDIAGL(IND)
    ENDDO
    IND=2*JF-1
    ZDIAG(IND  )= ZDIAG(IND  )+ZDIAGL(IND  )
    ZDIAG(IND+1)= ZDIAG(IND+1)+ZDIAGL(IND+1)
  ENDDO
  IGPCUM= IGPCHK
  DO JPROC=1,NPROC
    LLDONE= .FALSE.
    DO WHILE (.NOT. LLDONE .AND. JPROC /= 1)
      ITAGS= MTAGDISTGP+ (KSTEP+1)*NPROC + JPROC
      IERR= 0
      LLWAIT= .FALSE.
      LLEXIST= .FALSE.
      CALL MPL_PROBE(KSOURCE=JPROC,KTAG=ITAGS,LDWAIT=LLWAIT,&
       & LDFLAG=LLEXIST,CDSTRING='CHKEVO:')  
      IF (LLEXIST) THEN
        LLDONE=.TRUE.
        IERR= 0
        ILENR= 0
        ISNDR= 0
        ITAGR= 0
        ILLPCHK=NLENCHK+1
        CALL MPL_RECV(ZDIAGL(1:ILLPCHK),KSOURCE=NPRCIDS(JPROC),KTAG=ITAGS,&
         & KOUNT=ILENR,CDSTRING='CHKEVO:')  
        IF (ILENR /= (NLENCHK+1))CALL ABOR1('CHKEVO/DM : LENGTH ERROR')
        IGPLOC2= NINT(ZDIAGL(NLENCHK+1))
        DO JF=1,NFLDCHK
          DO J=1,IGPLOC2
            IND= 2*NFLDCHK +(JF-1)*NGPCHK + J
            ZDIAG(IND+IGPCUM)= ZDIAGL(IND)
          ENDDO
          IND=2*JF-1
          ZDIAG(IND  )= ZDIAG(IND  ) + ZDIAGL(IND  )
          ZDIAG(IND+1)= ZDIAG(IND+1) + ZDIAGL(IND+1)
        ENDDO
        IGPCUM= IGPCUM+IGPLOC2
      ELSE
        LLDONE= .FALSE.
      ENDIF
    ENDDO
  ENDDO
  IF (IGPCUM /= NGPCHK) CALL ABOR1('CHKEVO/DM : POINTS LOST')
  DO JF=1,2*NFLDCHK,2
    ZDIAG(JF  )= SQRT(ZDIAG(JF)/REAL(NGPTOTG,JPRB))
    ZDIAG(JF+1)= ZDIAG(JF+1)/REAL(NGPTOTG,JPRB)
  ENDDO
  ! Reordering tendencies and data
  ALLOCATE(ZDIAGS(NLENCHK))
  ZDIAGS(:)=0.0_JPRB
  ZDIAGS(1:2*NFLDCHK)=ZDIAG(1:2*NFLDCHK)
  DO J=1,NGPCHK
    DO JI=1,NGPCHK
      IF( (NXCHK(J) == NINT(ZGEOM(JI+5*NGPCHK))) .AND.&
       &(NYCHK(J) == NINT(ZGEOM(JI+6*NGPCHK))) ) THEN
        DO JF=1,NFLDCHK
          INDJ= 2*NFLDCHK +(JF-1)*NGPCHK + J
          INDJI= 2*NFLDCHK +(JF-1)*NGPCHK + JI
          ZDIAGS(INDJ)=ZDIAG(INDJI)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(ZDIAG)
  DEALLOCATE(ZGEOM)
ENDIF

!     ------------------------------------------------------------------

!           4. Write out for diagnostics
!              -------------------------

IF (MYPROC == 1) THEN

  WRITE(CLNSTEP,'(I4)') KSTEP
  CLNOMA='RE'//CLNSTEP//'0000000000'
  CALL FAISAN(IREP,NULUSR3,CLNOMA,ZDIAGS,NLENCHK)
  DEALLOCATE(ZDIAGS)
  IF (KSTEP == NSTOP) CALL FAIRME(IREP,NULUSR3,'UNKNOWN')

ENDIF

! Deallocate TENDCHK
IF (KSTEP == NSTOP .AND. ALLOCATED(TENDCHK)) DEALLOCATE(TENDCHK)
IF(ALLOCATED(ZDIAGL)) DEALLOCATE(ZDIAGL)
IF(ALLOCATED(ZGEOML)) DEALLOCATE(ZGEOML)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHKEVO',1,ZHOOK_HANDLE)
END SUBROUTINE CHKEVO
