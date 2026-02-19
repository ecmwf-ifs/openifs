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

SUBROUTINE VERINT(KPROMA,KSTART,KPROF,KFLEV_IN,KFLEV_OUT,PINTE,PIN,POUT,KTYPE,KCHUNK)

!**** *VERINT*   Vertical integral

!     Purpose.
!     --------
!          This subroutine computes the vertical integral (with respect
!          to eta) of a function given at full model
!          levels using a general scheme
!          The integral is given either from the top (KTYPE=0) down
!          (or from the bottom (KTYPE=1) up) to each full model level

!**   Interface.
!     ----------
!        *CALL* *VERINT(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA    - horizontal dimension.
!          KSTART    - first element of work.
!          KPROF     - depth of work.
!          KFLEV_IN  - vertical dimension for array PIN.
!          KFLEV_OUT - vertical dimension for array POUT.
!          PINTE     - matrix operator used to perform vertical integrals.
!          PIN       - Input field
!          KTYPE     - starting point of the integral (0=top, 1=bottom)
!          KCHUNK    - chunking size (to maintain bit reproducibility)

!        OUTPUT:
!                    eta
!                     _
!                    |
!          POUT :    | PIN deta   at each half model level
!                   _|
!                 KTYPE

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none

!     Reference.

!     Author.
!     -------
!         M. Hortal (ECMWF)

!     Modifications.
!     --------------
!        Original : MAY 2000
!        D.SALMOND : APRIL 2002 FORTRAN Matrix multiply replace by BLAS routine
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Y.Seity : MAY 2009 bf for B-Level parallelisation
!        J.Hague : OCT 2012: Parallelise call to DGEMM if not in parallel region
!        P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!        F. Vana  05-Mar-2015  Support for single precision
!        F. Vana  14-Jan-2020  Exclusive usage of double precision
!        P. Gillies & F. Vana  22-Jan-2020  Bit reproducible chunking
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULERR
USE OML_MOD  , ONLY : OML_IN_PARALLEL

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV_IN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV_OUT
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRD)   ,INTENT(IN)    :: PINTE(KFLEV_OUT,KFLEV_IN) 
REAL(KIND=JPRD)   ,INTENT(IN)    :: PIN(KPROMA,KFLEV_IN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(KPROMA,KFLEV_OUT) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCHUNK

!     ------------------------------------------------------------------

REAL(KIND=JPRD) :: ZOUT(KPROMA,KFLEV_OUT)

INTEGER(KIND=JPIM) ::  JLEV, JROF, IDB
!INTEGER(KIND=JPIM) ::  JDB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_XGEMM
INTEGER(KIND=JPIM) :: JLEN,ICHUNK
!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VERINT',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF(KTYPE /= 0.AND.KTYPE /= 1) THEN
  WRITE(NULERR,*) ' INVALID KTYPE IN VERINT =',KTYPE
  CALL ABOR1(' VERINT: ABOR1 CALLED')
ENDIF

! ALTERNATE FORTRAN CODE FOR DGEMM  - Please Don't delete
!DO JLEV=1,KFLEV_OUT
!  DO JROF=KSTART,KPROF
!    ZOUT(JROF,JLEV)=0.0_JPRD
!  ENDDO
!ENDDO
!DO JDB=1,KFLEV_IN
!  DO JLEV=1,KFLEV_OUT
!    DO JROF=KSTART,KPROF
!      ZOUT(JROF,JLEV)=ZOUT(JROF,JLEV)+PINTE(JLEV,JDB)*PIN(JROF,JDB)
!    ENDDO
!  ENDDO
!ENDDO


IDB=SIZE(PINTE,DIM=1)

IF(OML_IN_PARALLEL()) THEN

  IF (KPROF >= KSTART) THEN 
    IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_1',0,ZHOOK_HANDLE_XGEMM)
    CALL DGEMM('N','T',KPROF-KSTART+1,KFLEV_OUT,KFLEV_IN, &
         & 1.0_JPRD,PIN,KPROMA,PINTE,IDB,0.0_JPRD,ZOUT,KPROMA)  
    IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_1',1,ZHOOK_HANDLE_XGEMM)
  ENDIF

ELSEIF (KPROF >= KSTART) THEN ! Not in Parallel Region

  ICHUNK=1  ! Set = 0,1, or 2 to force type of call to DGEMM

  IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_2',0,ZHOOK_HANDLE_XGEMM)

  IF(ICHUNK==0) THEN  ! No chunking
    CALL DGEMM('N','T',KPROF-KSTART+1,KFLEV_OUT,KFLEV_IN, &
         & 1.0_JPRD,PIN,KPROMA,PINTE,IDB,0.0_JPRD,ZOUT,KPROMA)  
  ENDIF

  IF(ICHUNK==1) THEN  ! Chunking across NPROMA
!$OMP PARALLEL DO PRIVATE(JROF,JLEN)
    DO JROF=KSTART,KPROF,KCHUNK
      JLEN=MIN(KCHUNK,KPROF-JROF+1)
      CALL DGEMM('N','T',JLEN,KFLEV_OUT,KFLEV_IN, &
           & 1.0_JPRD,PIN(JROF,1),KPROMA,PINTE,IDB,0.0_JPRD,ZOUT(JROF,1),KPROMA)
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

  IF(ICHUNK==2) THEN  ! Chunkng across KFLEV_OUT
!$OMP PARALLEL DO PRIVATE(JLEV,JLEN)
    DO JLEV=1,KFLEV_OUT,KCHUNK
      JLEN=MIN(KCHUNK,KFLEV_OUT-JLEV+1)
      CALL DGEMM('N','T',KPROF-KSTART+1,JLEN,KFLEV_IN, &
           & 1.0_JPRD,PIN,KPROMA,PINTE(JLEV,1),IDB,0.0_JPRD,ZOUT(1,JLEV),KPROMA)
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_2',1,ZHOOK_HANDLE_XGEMM)

ENDIF

IF(.NOT.OML_IN_PARALLEL()) THEN
  IF(KTYPE == 1) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JROF)
    DO JLEV=1,KFLEV_OUT
      DO JROF=KSTART,KPROF
        POUT(JROF,JLEV)=ZOUT(JROF,JLEV)-ZOUT(JROF,KFLEV_OUT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JROF)
    DO JLEV=1,KFLEV_OUT
      DO JROF=KSTART,KPROF
        POUT(JROF,JLEV)=ZOUT(JROF,JLEV)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
ELSE
  IF(KTYPE == 1) THEN
    DO JLEV=1,KFLEV_OUT
      DO JROF=KSTART,KPROF
        POUT(JROF,JLEV)=ZOUT(JROF,JLEV)-ZOUT(JROF,KFLEV_OUT)
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,KFLEV_OUT
      DO JROF=KSTART,KPROF
        POUT(JROF,JLEV)=ZOUT(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VERINT',1,ZHOOK_HANDLE)
END SUBROUTINE VERINT
