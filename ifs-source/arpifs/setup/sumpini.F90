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

SUBROUTINE SUMPINI

!**** *SUMPINI*   - Preliminary calculations for distributed memory environment.

!     Purpose.
!     --------
!        Preliminary calculations for distributed memory environment.
!        Reads namelist NAMPAR0.
!        Set-up of a subset of YOMMP0 variables.

!**   Interface.
!     ----------
!        *CALL* *SUMPINI

!        explicit arguments : 
!        --------------------

!        implicit arguments :
!        --------------------

!     method.
!     -------
!        see documentation

!     externals.
!     ----------

!     reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     author.
!     -------
!      MPP Group *ECMWF*
!      original : 95-10-01

! Modifications
! -------------
!   09-Jun-2004 J. Masek   NH cleaning (LPC_NOTR, LFULLIMP, LVSLWBC)
!   15-Mar-2005 R. El Khatib setup N_VMASS moved to SUJFH
!   Modified 04-11-25 by Y. Seity Add LAROME
!   Modified 05-05-17 by A. Alias Add NFRCORM
!   Modified 06-06-06 M. Ko"hler   Single Column Model option (LSCMEC) 
!   Modified 01-11-07 G.Mozdzynski Use NSPECRESMIN to control max value of NPRTRW
!   R. El Khatin  11-May-2007 LIMP_NOOLAP=.T. as default
!   R. El Khatib 27-Apr-2007   LIMP_NOOLAP=.F. if no message passing
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   A. Alias (Aug 2009): Change NSHISTS to NSFXHISTS and add NFRSFXHIS, LCALLSFX added
!   S. Saarinen, CSC (Oct 2009) : Added llinfo
!   K. Yessad (Jan 2010): remove useless variables.
!   K. Yessad (Sep 2010): organigramme simplification.
!   K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!   A. Alias (Feb 2011): LSFXLSM added
!   R. El Khatib 10-Aug-2011 LIOLEVG management
!   D. Degrauwe  (Feb 2012): LARPEGEF_WRGP_HIST added
!   R. El Khatib : 01-Mar-2012 (LFPOS,LFPSPEC) => NFPOS
!   R. El Khatib 22-Mar-2012 Fix uninitialized variables
!   P. Marguinaud: 11-Sep-2012 : More namelist parameters
!   P. Marguinaud: 15-May-2013 : Remove LARPEGEF_WRGP_HIST
!   P. Marguinaud: 10-Oct-2013 : Add NHISTSMIN, NSFXHISTSMIN, NPOSTSMIN namelist parameters
!   K. Yessad (July 2013): minor modifications in printings (printing and setup in the same routine).
!   K. Yessad (Oct 2013): cleanings; LOUTPUT set-up there.
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): module+setup refactoring, in order to remove reading NAMCT0 in SUMPINI.
! End Modifications
!------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMLUN   , ONLY : NULNAM
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPROC, NPRGPNS, NPRGPEW, LMPOFF, NSPECRESMIN, NPRTRW, &
 & NPRTRNS, NPRTRV, NPRTRN, NOUTPUT, LOUTPUT, LMPDIAG, LVECADIN, &
 & MYPROC   ,NPRCIDS  ,MP_TYPE  ,MBX_SIZE , &
 & MYSETA   ,MYSETB   ,MYSETV   ,MYSETW   ,MYSETM   ,MYSETN, &
 & LSCMEC, NPRINTLEV, LOPT_SCALAR, LOPT_RS6K
USE YOMMPI   , ONLY : MINTET   ,MREALT
USE YOMTAG   , ONLY : MT_DISTRIBUTED_VECTOR
USE MPL_MODULE, ONLY : MPL_BUFFER_METHOD, MPL_INIT, MPL_MYRANK, MPL_GROUPS_CREATE, &
 & JP_BLOCKING_STANDARD, MPL_CART_COORDS
USE DISTRIBUTED_VECTORS_MIX, ONLY : SETUP_DISTVEC
USE OML_MOD, ONLY : OML_MAX_THREADS

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IDV_CHUNK_SIZE, ISQR, JA, IB, JPRTRV, JJ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IPRTRW,IPRTRV

! - the following ones are normally in different modules.
!   (why not "use" the corresponding module???)
INTEGER(KIND=JPIM) :: NTRACE_STATS
INTEGER(KIND=JPIM) :: NSTATS_MEM
INTEGER(KIND=JPIM) :: NPRNT_STATS
LOGICAL :: LSTATS_MEM,LSTATS_ALLOC
LOGICAL :: LSTATS,LSTATSCPU,LSYNCSTATS,LDETAILED_STATS,LXML_STATS,LLTRACE_STATS
LOGICAL :: LLSTATS_OMP, LLSTATS_COMMS,LBARRIER_STATS,LBARRIER_STATS2
LOGICAL :: LLINFO

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gstats_label_ifs.intfb.h"
#include "posnam.intfb.h"

#include "nampar0.nam.h"
#include "nam_distributed_vectors.nam.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUMPINI',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!     Open NULNAM before going parallel for the benifit of the VPP
!OPEN(NULNAM,ACTION="READ")

!     ------------------------------------------------------------------

!     1. Default values.
!     ------------------

!     1.1: YOMMP0 variables:

#ifdef RS6K
LOPT_SCALAR=.TRUE.
LOPT_RS6K=.TRUE.
#elif _CRAYFTN
LOPT_SCALAR=.TRUE.
LOPT_RS6K=.FALSE.
#elif __INTEL_COMPILER
LOPT_SCALAR=.TRUE.
LOPT_RS6K=.FALSE.
#elif __GFORTRAN__
LOPT_SCALAR=.TRUE.
LOPT_RS6K=.FALSE.
#elif __PGI
LOPT_SCALAR=.TRUE.
LOPT_RS6K=.FALSE.
#elif __NEC__
LOPT_SCALAR=.FALSE.
LOPT_RS6K=.FALSE.
#else
LOPT_SCALAR=.FALSE.
LOPT_RS6K=.FALSE.
#endif

IF(LOPT_SCALAR) THEN
  LVECADIN=.FALSE.
ELSE
  LVECADIN=.TRUE.
ENDIF

LSCMEC=.FALSE.
NPRINTLEV=0

NPROC=0
NPRGPNS=0
NPRGPEW=0
NPRTRW=0
NPRTRV=0
NSPECRESMIN=0
NOUTPUT=1
LMPOFF=.FALSE.
LMPDIAG=.FALSE.
MP_TYPE=JP_BLOCKING_STANDARD
MBX_SIZE=0
IDV_CHUNK_SIZE=1023

!     1.2: "STATS" variables:

LSTATS=.FALSE.
LSTATSCPU=.FALSE.
LSYNCSTATS=.FALSE.
LDETAILED_STATS=.FALSE.
LXML_STATS=.TRUE.
LBARRIER_STATS=.FALSE.
LBARRIER_STATS2=.TRUE.
NSTATS_MEM=0
LSTATS_MEM=.FALSE.
LSTATS_ALLOC=.FALSE.
NTRACE_STATS=0
NPRNT_STATS=2

!     ------------------------------------------------------------------

!     2. Read namelists.
!     ------------------

CALL POSNAM(NULNAM,'NAM_DISTRIBUTED_VECTORS')
READ(NULNAM,NAM_DISTRIBUTED_VECTORS)
CALL POSNAM(NULNAM,'NAMPAR0')
READ(NULNAM,NAMPAR0)

!     ------------------------------------------------------------------

!     3. Checkings and resettings ("STATS" variables).
!     ------------------------------------------------


IF (.NOT.LSTATS) LDETAILED_STATS=.FALSE.
IF (.NOT.LSTATS) LXML_STATS=.FALSE.
IF (.NOT.LSTATS) LBARRIER_STATS=.FALSE.
IF (.NOT.LBARRIER_STATS) LBARRIER_STATS2=.FALSE.

IF (LDETAILED_STATS) THEN
  LLSTATS_OMP=.TRUE.
  LLSTATS_COMMS=.TRUE.
  LLTRACE_STATS=.TRUE.
ELSE
  LLSTATS_OMP=.FALSE.
  LLSTATS_COMMS=.FALSE.
  LLTRACE_STATS=.FALSE.
ENDIF

!     ------------------------------------------------------------------

!     4. Additional actions.
!     ----------------------
  
!     4.1: finalise NPROC, NPRGPNS, NPRGPEW.

IF (NPRGPNS > 0 .AND. NPRGPEW > 0) THEN
  
  ! NPRGPNS and NPRGPEW are known, compute NPROC.
  IF (NPROC > 0 .AND. NPRGPNS*NPRGPEW /= NPROC) &
   & CALL ABOR1('SUMPINI: inconsistency between NPRGPNS, NPRGPEW, and NPROC')  
  NPROC=NPRGPNS*NPRGPEW

ELSE

  ! NPROC is known, compute NPRGPNS and NPRGPEW;
  ! this version selects most square-like distribution.
  IF (NPROC == 0) NPROC = 1
  ISQR=INT(SQRT(REAL(NPROC,JPRB)))
  DO JA=ISQR,NPROC
    IB=NPROC/JA
    IF (JA*IB == NPROC) THEN
      NPRGPNS=MAX(JA,IB)
      NPRGPEW=MIN(JA,IB)
      EXIT
    ENDIF
  ENDDO

ENDIF

!     4.2: finalise LMPOFF.

IF ( NPROC > 1 ) LMPOFF=.FALSE.

!     4.3: close and reopen NULNAM (is it still required?).

! I don't know why we have to close and re-open NULNAM but it fails otherwise on the vpp5000 /MH
!CLOSE(NULNAM)          
IF (.NOT.LMPOFF) THEN
  CALL MPL_INIT(KOUTPUT=NPRINTLEV)
!endif


ENDIF
CALL C_DRHOOK_INIT_SIGNALS(1) ! Get Dr.Hooks signals in case DR_HOOK was OFF
!OPEN(NULNAM,ACTION="READ")

!     4.4: finalise NSPECRESMIN, NPRTRW, NPRTRV, NPRTRNS, NPRTRN.

! Set defaults for NPRTRW and NPRTRV if necessary (where NPROC=NPRTRW*NPRTRV)
! Note NPRTRW will not be allowed to be greater than NSPECRESMIN

IF( NSPECRESMIN == 0 )THEN
  ! NSPECRESMIN was not specified in NAMPAR0 so give it a large default value
  NSPECRESMIN = NPROC
ENDIF

IF(NPRTRW > 0) THEN

  ! NPRTRW must have been set via the namelist read, check it is consistent
  ! i.e. NPROC equals NPRTRW * NPRTRV
  IF (NPROC > 0 .AND. NPRTRW*NPRTRV /= NPROC) &
   & CALL ABOR1(' SUMPINI: inconsistency between NPRTRW, NPRTRV, and NPROC')  
  IF (NPRTRW > NSPECRESMIN ) CALL ABOR1(' SUMPINI: NPRTRW > NSPECRESMIN')

ELSE

  IF(LOPT_SCALAR) THEN
    ! set default values of NPRTRW and NPRTRV for scalar systems
    DO JPRTRV=4,NPROC
      NPRTRV=JPRTRV
      NPRTRW=NPROC/NPRTRV
      IF( NPRTRV*NPRTRW /= NPROC ) CYCLE
      IF( NPRTRV > NPRTRW ) EXIT
      IF( NPRTRW > NSPECRESMIN ) CYCLE
      IF( NPRTRW <= NSPECRESMIN/(2*OML_MAX_THREADS()) ) EXIT
    ENDDO
    ! go for approx square partitioning as backup
    IF( NPRTRV*NPRTRW /= NPROC .OR. NPRTRW > NSPECRESMIN .OR. NPRTRV > NPRTRW ) THEN
      ISQR=INT(SQRT(REAL(NPROC,JPRB)))
      DO JA=ISQR,NPROC
        IB=NPROC/JA
        IF (JA*IB == NPROC) THEN
          NPRTRW=MAX(JA,IB)
          NPRTRV=MIN(JA,IB)
          IF (NPRTRW > NSPECRESMIN ) CALL ABOR1(' SUMPINI: NPRTRW (approx square value) > NSPECRESMIN')
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ELSE
    ! set default values of NPRTRW and NPRTRV for vector systems
    NPRTRW=NPROC
    NPRTRV=1
    DO JPRTRV=1,NPROC
      NPRTRV=JPRTRV
      NPRTRW=NPROC/NPRTRV
      IF( NPRTRV*NPRTRW /= NPROC .OR. NPRTRW > NSPECRESMIN ) THEN
        CYCLE
      ELSE
        EXIT
      ENDIF
    ENDDO
  ENDIF

ENDIF

NPRTRNS=NPRTRW
NPRTRN=NPRTRV

! Create communicators for MPI groups
IF (.NOT.LMPOFF) THEN
  CALL MPL_GROUPS_CREATE(NPRTRW,NPRTRV)
ENDIF

!     4.5: Define processor topology and environment:
!          NPRCIDS, MYPROC, MYSETA, MYSETB, MYSETW, MYSETV, MYSETM, MYSETN.

! compute MYPROC and NPRCIDS:
ALLOCATE(NPRCIDS(NPROC))
IF ( .NOT.LMPOFF ) THEN
  MYPROC=MPL_MYRANK()
  IF (MYPROC > NPROC) CALL ABOR1(' SUMPINI: inconsistency between MYPROC and NPROC')
ELSE
  MYPROC=1
ENDIF
DO JJ=1,NPROC
  NPRCIDS(JJ)=JJ
ENDDO

! compute MYSETA and MYSETB (must be consistent with the content of PE2SET case LEQ_REGIONS=F):
MYSETB=MOD(MYPROC-1,NPRGPEW)+1
MYSETA=(MYPROC-1)/NPRGPEW+1

! compute MYSETW and MYSETV (must be consistent with the content of PE2SET):
IF (LMPOFF) THEN
  MYSETW=(MYPROC-1)/NPRTRV+1
  MYSETV=MOD(MYPROC-1,NPRTRV)+1
ELSE
  CALL MPL_CART_COORDS(MYPROC,MYSETW,MYSETV)
  ! Just checking for now...
  IPRTRV=MOD(MYPROC-1,NPRTRV)+1
  IPRTRW=(MYPROC-1)/NPRTRV+1
  IF (IPRTRV/=MYSETV .OR. IPRTRW/=MYSETW) THEN
    CALL ABOR1('SUMPINI: case LMPOFF=F: inconsistency when computing MYSETW and MYSETV')
  ENDIF
ENDIF

! compute MYSETM and MYSETN:
MYSETM=MYSETW
MYSETN=MYSETV

!     4.6: Compute LOUTPUT:

IF (NOUTPUT <= 0) THEN
  LOUTPUT=.FALSE.
ELSEIF (NOUTPUT == 1) THEN
  IF (MYSETA == 1.AND.MYSETB == 1) THEN
    LOUTPUT=.TRUE.
  ELSE
    LOUTPUT=.FALSE.
  ENDIF
ELSE
  LOUTPUT=.TRUE.
ENDIF

!     4.7: Additional actions:

IF (.NOT.LMPOFF) THEN
  LLINFO=.FALSE.
  IF (MYPROC == 1) LLINFO=.TRUE.
  CALL MPL_BUFFER_METHOD(KMP_TYPE=MP_TYPE,KMBX_SIZE=MBX_SIZE,KPROCIDS=NPRCIDS,LDINFO=LLINFO)
ENDIF

! finalise "STATS" variables.
CALL GSTATS_SETUP(NPROC,MYPROC,NPRCIDS,&
 & LSTATS,LSTATSCPU,LSYNCSTATS,LDETAILED_STATS,LBARRIER_STATS,LBARRIER_STATS2,&
 & LLSTATS_OMP,LLSTATS_COMMS,LSTATS_MEM,NSTATS_MEM,LSTATS_ALLOC,&
 & LLTRACE_STATS,NTRACE_STATS,NPRNT_STATS,LXML_STATS)  
CALL GSTATS_PSUT
CALL GSTATS_LABEL_IFS

CALL SETUP_DISTVEC(NPROC,MYPROC,MINTET,MREALT,&
 & MT_DISTRIBUTED_VECTOR,IDV_CHUNK_SIZE,NPRCIDS)  

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUMPINI',1,ZHOOK_HANDLE)
END SUBROUTINE SUMPINI
