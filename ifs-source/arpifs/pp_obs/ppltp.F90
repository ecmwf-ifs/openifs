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

SUBROUTINE PPLTP(PVCAP,KITERPV,LDISOPV,KPROMA,KSTART,KPROF,KFLEV,&
 & PRPRESF,PTOUP,PXTOUP,LDCP,PCORIO,&
 & PRPRES,PRLNPRES,PCP2D)  

!**** *PPLTP*    Compute the pressures of the TP level

!      PURPOSE.
!      --------
!           COMPUTE THE PRESSURES OF THE TP LEVEL

!**    INTERFACE.
!      ----------
!           *CALL* PPLTP( ... )

!           EXPLICITE ARGUMENTS.
!           --------------------
!           PVCAP       : Minimum pressure of model level to provide an equatorial
!                         cap in the computation of variables on constant PV surfaces
!           LDISOPV     : .T. => new diagnostic for computing an iso-PV level
!           KITERPV     : Nb of vertical iter (1, 2 or 3) used in iso-PV level computing
!           KPROMA      : horizontal dimension.                     (INPUT)
!           KSTART      : start of work.                            (INPUT)
!           KPROF       : depth of work.                            (INPUT)
!           KFLEV       : number of input pressure levels.          (INPUT)
!           PRPRESF(kproma,kflev) : model full level pressures.     (INPUT)
!           PTOUP(kproma,kflev)   : TP.                             (INPUT)
!           PXTOUP      : value of the TP level.                    (INPUT)
!           LDCP        : .true. if post-proc of tropo folding     (INPUT)
!           PCORIO      : Coriolis parameter                       (input)
!           PRPRES(kproma)    : pressures of the TP level.         (OUTPUT)
!           PRLNPRES(kproma)  : ln-pressure of the TP level.       (OUTPUT)
!           PCP2D       : tropopause folding                (OUTPUT)

!      EXTERNALS.
!      ----------
!           NONE

!      AUTHOR.
!      -------
!       C. PERIARD
!       ORIGINAL : 1992

!      MODIFICATIONS.
!      --------------
!       M.Hamrud      01-Oct-2003 CY28 Cleaning
!       R. El Khatib  20-Apr-2007 Optimisation.
!       K.Maynard     17-Apr-2009 Optimisation in case of area folding and 
!                     stratospheric "bulbe" of potential vorticity.
!        P.Marguinaud  15-Jun-2010 Inline ISRCHFGE & ISRCHFLTPV
!      -----------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCAP
INTEGER(KIND=JPIM),INTENT(IN)    :: KITERPV
LOGICAL           ,INTENT(IN)    :: LDISOPV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOUP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXTOUP 
LOGICAL           ,INTENT(IN)    :: LDCP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCORIO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRPRES(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRLNPRES(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCP2D(KPROMA) 
#include "isrchfge.decl.h"
!! INTEGER(KIND=JPIM), EXTERNAL :: ISRCHFGE
INTEGER(KIND=JPIM), EXTERNAL :: ISRCHFLT
#include "isrchfltpv.decl.h"

REAL(KIND=JPRB) :: ZPVPLUS(KPROMA,KFLEV)

INTEGER(KIND=JPIM) :: ILEFLEV, ILEV2, ILEVI, JL, JLEV, JROF
INTEGER(KIND=JPIM) :: ILEVSEC(KPROMA), ILEV(KPROMA), INBLEV(KPROMA), I

REAL(KIND=JPRB) :: ZQNN, ZSECUR, ZPMIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPLTP',0,ZHOOK_HANDLE)
ZSECUR=1.E-10_JPRB
ZPMIN=PVCAP

!     Correct southern hemisphere negative values
DO JLEV = 1, KFLEV
  DO JROF = KSTART, KPROF
    ZPVPLUS(JROF,JLEV) = SIGN(1.0_JPRB,PCORIO(JROF))*PTOUP(JROF,JLEV)
  ENDDO
ENDDO


! X   array name
! LB1 lower bound for dimension 1
! SZ1 size of dimension 1
! I2  indice of dimension 2
#define ARRAY(X,I,LB1,SZ1,I2)  X(1+MOD(I+LB1-2,SZ1),I2+(I+LB1-2)/(SZ1))

IF (LDCP) THEN
  DO JL=KSTART,KPROF
#define ISRCHFGE_RESULT     ILEVSEC(JL)
#define ISRCHFGE_N          KFLEV
#define ISRCHFGE_INC        KPROMA
#define ISRCHFGE_ARRAY(I)   ARRAY(PRPRESF,I,JL,KPROMA,1)
#define ISRCHFGE_TARGET     ZPMIN
#include "isrchfge.body.h"
#undef ISRCHFGE_RESULT
#undef ISRCHFGE_N
#undef ISRCHFGE_INC
#undef ISRCHFGE_ARRAY
#undef ISRCHFGE_TARGET
!   ILEVSEC(JL)=ISRCHFGE(KFLEV,PRPRESF(JL,1),KPROMA,ZPMIN)
    INBLEV(JL) = KFLEV - 1
    IF (LDISOPV) THEN
#define ISRCHFLTPV_N             INBLEV(JL)
#define ISRCHFLTPV_ARRAY(I)      ARRAY(ZPVPLUS,I,JL,KPROMA,2)
#define ISRCHFLTPV_INC           KPROMA
#define ISRCHFLTPV_TARGET        PXTOUP
#define ISRCHFLTPV_NBITER        KITERPV
#define ISRCHFLTPV_RESULT        ILEV(JL)
#include "isrchfltpv.body.h"
#undef ISRCHFLTPV_N    
#undef ISRCHFLTPV_ARRAY
#undef ISRCHFLTPV_INC 
#undef ISRCHFLTPV_TARGET
#undef ISRCHFLTPV_NBITER
#undef ISRCHFLTPV_RESULT
    ELSE
       ILEV(JL)=ISRCHFLT(INBLEV(JL),ZPVPLUS(JL,2),KPROMA,PXTOUP)
    ENDIF
    ILEV(JL) = ILEV(JL) + 1
    ILEVI=MAX(ILEV(JL),ILEVSEC(JL))
    IF (ILEVI == 1.OR. ILEVI >= KFLEV+1) THEN
      ILEVI=MIN(ILEVI,KFLEV)
      PRPRES(JL)=PRPRESF(JL,ILEVI)
      PRLNPRES(JL)=LOG(PRPRES(JL))
    ELSE
      IF (LDISOPV) THEN
         IF ( (ZPVPLUS(JL,ILEVI-1)-ZPVPLUS(JL,ILEVI)) >= 0._JPRB ) THEN
            ZQNN=MAX(ZPVPLUS(JL,ILEVI-1)-ZPVPLUS(JL,ILEVI),ZSECUR)
         ELSE
            ZQNN=MIN(ZPVPLUS(JL,ILEVI-1)-ZPVPLUS(JL,ILEVI),-ZSECUR)
         ENDIF
      ELSE
         ZQNN=MAX(ZPVPLUS(JL,ILEVI-1)-ZPVPLUS(JL,ILEVI),ZSECUR)
      ENDIF
      PRPRES(JL)=(PRPRESF(JL,ILEV(JL)-1)-PRPRESF(JL,ILEV(JL)))*&
       & (PXTOUP-ZPVPLUS(JL,ILEVI))/ZQNN &
       & + PRPRESF(JL,ILEVI)  
      PRPRES(JL)=MAX(PRPRES(JL),PRPRESF(JL,ILEVSEC(JL)))
      PRLNPRES(JL)=LOG(PRPRES(JL))
    ENDIF
    ILEFLEV = MAX(1,KFLEV-ILEV(JL))
#define ISRCHFGE_RESULT    ILEV2
#define ISRCHFGE_N         ILEFLEV
#define ISRCHFGE_INC       KPROMA
#define ISRCHFGE_ARRAY(I)  ARRAY(ZPVPLUS,I,JL,KPROMA,ILEV(JL))
#define ISRCHFGE_TARGET    PXTOUP
#include "isrchfge.body.h"
#undef ISRCHFGE_RESULT
#undef ISRCHFGE_N
#undef ISRCHFGE_INC
#undef ISRCHFGE_ARRAY
#undef ISRCHFGE_TARGET
!   ILEV2=ISRCHFGE(ILEFLEV,ZPVPLUS(JL,ILEV(JL)),KPROMA,PXTOUP)
    IF (ILEV2 >= ILEFLEV+1) THEN
      PCP2D(JL)=0.0_JPRB
    ELSE
      PCP2D(JL)=1.0_JPRB
    ENDIF
  ENDDO
ELSE
  DO JL=KSTART,KPROF
#define ISRCHFGE_RESULT    ILEVSEC(JL)
#define ISRCHFGE_N         KFLEV
#define ISRCHFGE_INC       KPROMA
#define ISRCHFGE_ARRAY(I)  ARRAY(PRPRESF,I,JL,KPROMA,1)
#define ISRCHFGE_TARGET    ZPMIN
#include "isrchfge.body.h"
#undef ISRCHFGE_RESULT
#undef ISRCHFGE_N
#undef ISRCHFGE_INC
#undef ISRCHFGE_ARRAY
#undef ISRCHFGE_TARGET
!   ILEVSEC(JL)=ISRCHFGE(KFLEV,PRPRESF(JL,1),KPROMA,ZPMIN)
  ENDDO
  DO JL=KSTART,KPROF
    INBLEV(JL) = KFLEV-ILEVSEC(JL)+1
  ENDDO
  IF (LDISOPV) THEN
    DO JL=KSTART,KPROF
#define ISRCHFLTPV_N             INBLEV(JL)
#define ISRCHFLTPV_ARRAY(I)      ARRAY(ZPVPLUS,I,JL,KPROMA,ILEVSEC(JL))
#define ISRCHFLTPV_INC           KPROMA
#define ISRCHFLTPV_TARGET        PXTOUP
#define ISRCHFLTPV_NBITER        KITERPV
#define ISRCHFLTPV_RESULT        ILEV(JL)
#include "isrchfltpv.body.h"
#undef ISRCHFLTPV_N    
#undef ISRCHFLTPV_ARRAY
#undef ISRCHFLTPV_INC 
#undef ISRCHFLTPV_TARGET
#undef ISRCHFLTPV_NBITER
#undef ISRCHFLTPV_RESULT
    ENDDO 
  ELSE
    DO JL=KSTART,KPROF
       ILEV(JL)=ISRCHFLT(INBLEV(JL),ZPVPLUS(JL,ILEVSEC(JL)),KPROMA,PXTOUP)
    ENDDO
  ENDIF
  DO JL=KSTART,KPROF
    ILEV(JL)=MAX(ILEV(JL)+ILEVSEC(JL)-1,ILEVSEC(JL))
  ENDDO
  DO JL=KSTART,KPROF
    IF (ILEV(JL) == ILEVSEC(JL) .OR. ILEV(JL) >= KFLEV+1) THEN
      I=MIN(ILEV(JL),KFLEV)
      PRPRES(JL)=PRPRESF(JL,I)
      PRLNPRES(JL)=LOG(PRPRES(JL))
    ELSE
      I=ILEV(JL)
      ZQNN=MAX(ZPVPLUS(JL,I-1)-ZPVPLUS(JL,I),ZSECUR)
      PRPRES(JL)=(PRPRESF(JL,I-1)-PRPRESF(JL,I))*&
       & (PXTOUP-ZPVPLUS(JL,I))/ZQNN + PRPRESF(JL,I)  
      PRPRES(JL)=MAX(PRPRES(JL),PRPRESF(JL,ILEVSEC(JL)))
      PRLNPRES(JL)=LOG(PRPRES(JL))
    ENDIF
  ENDDO
ENDIF

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPLTP',1,ZHOOK_HANDLE)

#undef ARRAY

END SUBROUTINE PPLTP
