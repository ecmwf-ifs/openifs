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

SUBROUTINE DIWRGRFP(KTAGDISTGP,KNFG,KDIMGP,KGPL,KGPLX,KNFL, &
 & KGPIND,KGPD,KNFD,KFLDOFF,PREALL,PREALG)

!**** *DIWRGRFP* - Routine to send and receive all distributed output 
!                  grid point fields for FULL-POS fields

!     Purpose.
!     --------
!       To transpose a gridpoint-distributed set of fields to a 
!       field-distributed set of global gridpoint fields.

!**   Interface.
!     ----------
!        *CALL* *DIWRGRFP*

!        Explicit arguments :
!        --------------------

!        * INPUT:
!        KTAGDISTGP: tag for processor communication.
!        KNFG   : global number of fields
!        KDIMGP : first dimension of PREALG (over-dimensioning for gridpoints)
!        KGPL   : number of gridpoints on the processor at output distribution
!        KGPLX  : maximum number of gridpoints (taken among all processors)
!        KNFL   : number of local fields on output distribution
!        KGPIND : KGPIND(i,j)=k mean: i-th gridpoint on j-th processor is k-th
!                 global gp
!        KGPD   : number of gridpoints on each processor
!        KNFD   : number of fields on output distribution for each processor
!        KFLDOFF: offset to the first local field for each processor
!        PREALL : gridpoint-distributed fields

!        * OUTPUT:
!        PREALG : field-distributed gridpoints

!        Implicit arguments :
!        --------------------
!           See modules above.

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
!      Valaparaiso Team at Meteo France
!      Original : 98-04    

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad   : 31-Jan-2007: cleaning, tag becomes dummy argument.
!      P.Marguinaud:26-Apr-2012: Add extra argument KHDRS (IO server header size)
!      R. El Khatib 07-Aug-2013 Optimization
!      P.Marguinaud:10-Oct-2013: Remove KHDRS argument
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     -------------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0    , ONLY : MYPROC, NPRCIDS, NPROC
USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_WAIT, JP_NON_BLOCKING_STANDARD

!     -------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KTAGDISTGP
INTEGER(KIND=JPIM),INTENT(IN)    :: KNFG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIMGP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPLX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNFL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPIND(KGPLX,NPROC) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPD(NPROC) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNFD(NPROC) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDOFF(NPROC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREALL(KGPL,KNFG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREALG(:,:)        ! KDIMGP, KNFL

!     -------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ITAG, JROC, IRSND, IRMAX, IOFFPROC
INTEGER(KIND=JPIM) :: II, IFLDLOC, IFLDG, JI
INTEGER(KIND=JPIM) :: ISENDS(NPROC), IRECVS(NPROC), IREQRCV(NPROC)

INTEGER(KIND=JPIM), ALLOCATABLE :: IREQSND(:)

REAL(KIND=JPRB), ALLOCATABLE :: ZBUF(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------------------------

#include "abor1.intfb.h"

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DIWRGRFP',0,ZHOOK_HANDLE)

!     -------------------------------------------------------------------------

!*     0. PREPARATIONS
!         ------------

IF (KNFL /= KNFD(MYPROC)) THEN
  CALL ABOR1('DIWRGRFP : INTERNAL ERROR KNFL/KNFD')
ENDIF
IF (KGPL /= KGPD(MYPROC)) THEN
  CALL ABOR1('DIWRGRFP : INTERNAL ERROR KGPL/KGPD')
ENDIF
IF (KGPLX < MAXVAL(KGPD)) THEN
  CALL ABOR1('DIWRGRFP : INTERNAL ERROR KGPLX/KGPD')
ENDIF


!     -------------------------------------------------------------------------

!*     1. NON-BLOCKING PEER-TO-PEER COMMUNICATIONS
!         ----------------------------------------

!      1.1 Post recv first

DO JROC = 1, NPROC
  IRECVS(JROC)=KGPD(JROC)*KNFL
ENDDO
ALLOCATE(ZBUF(SUM(IRECVS(:))))
IOFFPROC=0
DO JROC = 1, NPROC
  IF (IRECVS(JROC) > 0 .AND. JROC /= MYPROC) THEN
    ITAG = KTAGDISTGP + JROC
    CALL MPL_RECV(ZBUF(IOFFPROC+1:IOFFPROC+IRECVS(JROC)), &
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQRCV(JROC), &
     & KOUNT=IRECVS(JROC),KSOURCE=NPRCIDS(JROC),KTAG=ITAG,CDSTRING='DIWRGRFP:')  
    IOFFPROC=IOFFPROC+IRECVS(JROC)
  ENDIF
ENDDO

!      1.2 Sends

IRMAX=NPROC*KNFG ! Maximum number of send/recv requests
ALLOCATE(IREQSND(IRMAX))
DO JROC = 1, NPROC
  ISENDS(JROC)=KGPL*KNFD(JROC)
ENDDO
IRSND=0
DO JROC = 1, NPROC
  IF (ISENDS(JROC) > 0 .AND. JROC /= MYPROC) THEN
    ITAG = KTAGDISTGP + MYPROC
    IRSND=IRSND+1
    CALL MPL_SEND(PREALL(:,KFLDOFF(JROC)+1:KFLDOFF(JROC)+KNFD(JROC)), &
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQSND(IRSND), &
     & KTAG=ITAG,KDEST=NPRCIDS(JROC),CDSTRING='DIWRGRFP:')  
  ENDIF
ENDDO

!      1.3 Local contribution

IF (KGPL*KNFL > 0) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JI,IFLDLOC,II,IFLDG)
  DO JI=1,KGPL*KNFL
    IFLDLOC=(JI-1)/KGPL+1
    II=JI-((IFLDLOC-1)*KGPL)
    IFLDG=KFLDOFF(MYPROC)+IFLDLOC
    PREALG(KGPIND(II,MYPROC),IFLDLOC)=PREALL(II,IFLDG)
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!      1.4 Reorder data from communication buffers

IOFFPROC=0
DO JROC = 1, NPROC
  IF (IRECVS(JROC) > 0 .AND. JROC /= MYPROC) THEN
    CALL MPL_WAIT(KREQUEST=IREQRCV(JROC),CDSTRING='DIWRGRFP: WAIT FOR RECVS')
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JI,IFLDLOC,II)
    DO JI=1,KGPD(JROC)*KNFL
      IFLDLOC=(JI-1)/KGPD(JROC)+1
      II=JI-((IFLDLOC-1)*KGPD(JROC))
      PREALG(KGPIND(II,JROC),IFLDLOC)=ZBUF(JI+IOFFPROC)
    ENDDO
!$OMP END PARALLEL DO
    IOFFPROC=IOFFPROC+IRECVS(JROC)
  ENDIF
ENDDO
DEALLOCATE(ZBUF)

!      1.5 Wait completion of sends

IF (IRSND > 0) THEN
  CALL MPL_WAIT(KREQUEST=IREQSND(1:IRSND),CDSTRING='DIWRGRFP: WAIT FOR SENDS')
ENDIF
DEALLOCATE(IREQSND)

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DIWRGRFP',1,ZHOOK_HANDLE)
END SUBROUTINE DIWRGRFP
