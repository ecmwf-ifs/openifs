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

SUBROUTINE SUTRANS0(CDMODE)

!**** *SUTRANS0 * - Initialize the transform package: part which does not depend on geometry

!     Purpose.
!     --------
!       Initialize the transform package: part which does not depend on geometry.
!       In particuliar, some distributed environment present in YOMMP0 and in the
!       transform package too is set-up there.

!**   Interface.  CALL SUTRANS0
!     ---------- 

!     Explicit arguments :
!     --------------------
!        CDMODE    :

!     Externals.
!     ----------
!        SETUP_TRANS0 - basic initialization

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        T.Wilhelmsson 16-Aug-2013
!         (resolution independent setup taken from SUTRANS)

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 02-Jul-2015 Memory management
!       R. El Khatib 01-Sep-2015 Default value of RDISTR_E changed to 1.
!      R. El Khatib 16-May-2019 optimize memory access in NGPSET2PE
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE PARFPOS  , ONLY : JPOSDOM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULERR, NULNAM
USE YOMCT0   , ONLY : NOPT_MEMORY
USE YOMMP0   , ONLY : NPROC, NPRGPNS, NPRGPEW, NPRTRW, LOUTPUT, LMPOFF, &
 & NCOMBFLEN, LEQ_REGIONS, LSYNC_TRANS, NTRANS_SYNC_LEVEL, NPRINTLEV, &
 & N_REGIONS_NS, N_REGIONS_EW, N_REGIONS, MYPROC,MY_REGION_NS,MY_REGION_EW, NGPSET2PE
USE YOMTRANS , ONLY : NPROMATR, NEPROMATR, NMAX_RESOL, RDISTR_E, LMONO_TRANS, LALLOPERM
USE YOMCST   , ONLY : RA

!      -----------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=5)  , OPTIONAL, INTENT(IN)  :: CDMODE

!      -----------------------------------------------------------------

INTEGER(KIND=JPIM) :: I_REGIONS(NPROC+2)
INTEGER(KIND=JPIM) :: JA, JB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
CHARACTER(LEN=5) :: CLMODE

!      -----------------------------------------------------------------

#include "namtrans0.nam.h"

#include "setup_trans0.h"
#include "trans_end.h"

#include "pe2set.intfb.h"
#include "set2pe.intfb.h"
#include "abor1.intfb.h"
#include "posnam.intfb.h"

!      -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUTRANS0',0,ZHOOK_HANDLE)
!      -----------------------------------------------------------------

!*     1.  Default values for namelist variables

CLMODE='SETUP'
IF (PRESENT(CDMODE)) CLMODE=CDMODE

LMONO_TRANS=(CLMODE == 'MONO1')
NPROMATR = 0
NEPROMATR = 0
NMAX_RESOL = 22 
RDISTR_E=1._JPRB
LALLOPERM=(NOPT_MEMORY >= 1)

!      -----------------------------------------------------------------

!*     2.  Read namelist

CALL POSNAM(NULNAM,'NAMTRANS0')
READ(NULNAM,NAMTRANS0)

!      -----------------------------------------------------------------

!*     3.  Checkings

NMAX_RESOL=NMAX_RESOL+JPOSDOM

!      -----------------------------------------------------------------

!*     4.  Setup

! distributed memory environment used in the transform package,
! and N_REGIONS, N_REGIONS_NS, N_REGIONS_EW.
IF (CLMODE == 'MONO1') THEN
  CALL TRANS_END('INTER')
  CALL SETUP_TRANS0(KOUT=NULOUT,KERR=NULERR,KMAX_RESOL=NMAX_RESOL,LDALLOPERM=LALLOPERM) 
ELSEIF (CLMODE == 'RESET' ) THEN
  CALL TRANS_END('INTER')
  CALL SETUP_TRANS0(KOUT=NULOUT,KERR=NULERR,KPRGPNS=NPRGPNS,&
   & KPRGPEW=NPRGPEW,KPRTRW=NPRTRW,KPRINTLEV=NPRINTLEV, &
   & KMAX_RESOL=NMAX_RESOL,&
   & KPROMATR=NPROMATR,KCOMBFLEN=NCOMBFLEN,LDSYNC_TRANS=LSYNC_TRANS,&
   & KTRANS_SYNC_LEVEL=NTRANS_SYNC_LEVEL,LDMPOFF=LMPOFF,LDEQ_REGIONS=LEQ_REGIONS,&
   & K_REGIONS=I_REGIONS(:),K_REGIONS_NS=N_REGIONS_NS,K_REGIONS_EW=N_REGIONS_EW,&
   & PRAD=RA,LDALLOPERM=LALLOPERM)
   ! now that we know the actual size of N_REGIONS we can allocate it and copy
   ! over values from the temporary location
   IF (.NOT. ALLOCATED(N_REGIONS)) ALLOCATE(N_REGIONS(1:N_REGIONS_NS))
   N_REGIONS(1:N_REGIONS_NS)=I_REGIONS(1:N_REGIONS_NS)
ELSEIF (CLMODE == 'SETUP' ) THEN
  CALL SETUP_TRANS0(KOUT=NULOUT,KERR=NULERR,KPRGPNS=NPRGPNS,&
   & KPRGPEW=NPRGPEW,KPRTRW=NPRTRW,KPRINTLEV=NPRINTLEV, &
   & KMAX_RESOL=NMAX_RESOL,&
   & KPROMATR=NPROMATR,KCOMBFLEN=NCOMBFLEN,LDSYNC_TRANS=LSYNC_TRANS,&
   & LDMPOFF=LMPOFF,LDEQ_REGIONS=LEQ_REGIONS,&
   & K_REGIONS=I_REGIONS(:),K_REGIONS_NS=N_REGIONS_NS,K_REGIONS_EW=N_REGIONS_EW,&
   & PRAD=RA,LDALLOPERM=LALLOPERM)
   ! now that we know the actual size of N_REGIONS we can allocate it and copy
   ! over values from the temporary location
  ALLOCATE(N_REGIONS(1:N_REGIONS_NS))
  N_REGIONS(1:N_REGIONS_NS)=I_REGIONS(1:N_REGIONS_NS)
ELSE
  CALL ABOR1('WRONG SUTRANS0 MODE')
ENDIF

! compute MY_REGION_NS and MY_REGION_EW
CALL PE2SET(KPE=MYPROC,KPRGPNS=MY_REGION_NS,KPRGPEW=MY_REGION_EW)

! initialise NGPSET2PE for faster mapping of gp sets to proc
IF (.NOT. ALLOCATED(NGPSET2PE)) ALLOCATE(NGPSET2PE(N_REGIONS_EW,N_REGIONS_NS))
DO JA=1,N_REGIONS_NS
  DO JB=1,N_REGIONS(JA)
    CALL SET2PE(NGPSET2PE(JB,JA),JA,JB,0,0)
  ENDDO
  DO JB=N_REGIONS(JA)+1,N_REGIONS_EW
    NGPSET2PE(JB,JA)=0
  ENDDO
ENDDO

!      -----------------------------------------------------------------

!*     5.  Printings

WRITE(NULOUT,*) 'Printings in SUTRANS0:'
WRITE(NULOUT,*) ' -- NAMTRANS variables:'
WRITE(NULOUT,'('' NMAX_RESOL = '', I3,'' NPROMATR = '',I4, '' NEPROMATR = '',I4, &
 & '' LALLOPERM = '',L2)') NMAX_RESOL,NPROMATR,NEPROMATR,LALLOPERM
WRITE(NULOUT,*) ' RDISTR_E=',RDISTR_E


IF (LOUTPUT) THEN
  WRITE(NULOUT,*) ' -- Some YOMMP0 variables:'
  WRITE(NULOUT,'("N_REGIONS_NS=",I4)')N_REGIONS_NS
  WRITE(NULOUT,'("N_REGIONS_EW=",I4)')N_REGIONS_EW
  WRITE(NULOUT,'("MY_REGION_NS=",I4)')MY_REGION_NS
  WRITE(NULOUT,'("MY_REGION_EW=",I4)')MY_REGION_EW
  IF( LEQ_REGIONS )THEN
    WRITE(NULOUT,'(" ")')
    DO JA=1,N_REGIONS_NS
      WRITE(NULOUT,'("N_REGIONS(",I5,")=",I3)')JA,N_REGIONS(JA)
    ENDDO
  ENDIF
ENDIF

!      -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUTRANS0',1,ZHOOK_HANDLE)
END SUBROUTINE SUTRANS0
