! (C) Copyright 1998- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUMPLATB_MOD
CONTAINS
SUBROUTINE SUMPLATB(KDGSA,KDGL,KPROCA,KLOENG,LDSPLIT,LDFOURIER,&
                    &KMEDIAP,KRESTM,KINDIC,KLAST)

!**** *SUMPLATB * - Routine to initialize parallel environment

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUMPLATB *

!     Explicit arguments - input :
!     --------------------
!                          KDGSA      -first latitude (grid-space)
!                                      (may be different from NDGSAG)
!                          KDGL       -last  latitude
!                          KPROCA     -number of processors in A direction
!                          KLOENG     -actual number of longitudes per latitude.
!                          LDSPLIT    -true for latitudes shared between sets
!                          LDFOURIER  -true for fourier space partitioning

!     Explicit arguments - output:
!     --------------------
!                          KMEDIAP    -mean number of grid points per PE
!                          KRESTM     -number of PEs with one extra point
!                          KINDIC     -intermediate quantity for 'sumplat'
!                          KLAST      -intermediate quantity for 'sumplat'

!        Implicit arguments :
!        --------------------


!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD (after old version of sumplat.F).

!     Modifications.
!     --------------
!        Original : 98-12-07
!        G. Mozdzynski (August 2012): rewrite of fourier latitude distribution
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB

USE TPM_DISTR
USE ABORT_TRANS_MOD

IMPLICIT NONE


!     * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KLOENG(KDGSA:KDGL)
LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(IN)  :: LDFOURIER
INTEGER(KIND=JPIM),INTENT(OUT)  :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT)  :: KRESTM
INTEGER(KIND=JPIM),INTENT(OUT)  :: KINDIC(KPROCA)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KLAST(KPROCA)

!     * LOCAL:
INTEGER(KIND=JPIB) :: ICOST(KDGSA:KDGL)
INTEGER(KIND=JPIM) :: ILATS(KPROCA)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ICOMP, IGL, JA, JGL, ILAST, IREST, IA
INTEGER(KIND=JPIM) :: ITOT_TOP, ITOT_BOT, IGL_TOP, IGL_BOT
INTEGER(KIND=JPIB) :: IMEDIA,ITOT
REAL(KIND=JPRB) :: ZLG
LOGICAL   :: LLDONE,LLSIMPLE

!      -----------------------------------------------------------------

!*       1.    COMPUTATION OF KMEDIAP, KRESTM, KINDIC, KLAST.
!              ----------------------------------------------

!     * Computation of KMEDIAP and KRESTM.

IF( LDFOURIER )THEN

! DO JGL=1,KDGL
!   ZLG=LOG(FLOAT(KLOENG(JGL)))
!   ICOST(JGL)=KLOENG(JGL)*ZLG*SQRT(ZLG)
! ENDDO

  DO JGL=1,KDGL
    ICOST(JGL)=KLOENG(JGL)
  ENDDO

ELSE

  DO JGL=1,KDGL
    ICOST(JGL)=KLOENG(JGL)
  ENDDO

ENDIF
  
IMEDIA = SUM(ICOST(KDGSA:KDGL))
KMEDIAP = IMEDIA / KPROCA
KRESTM = IMEDIA - KMEDIAP * KPROCA
IF (KRESTM  >  0) KMEDIAP = KMEDIAP + 1

!     * Computation of intermediate quantities KINDIC and KLAST

KINDIC(:)=0
KLAST(:)=0

IF (LDSPLIT) THEN

  IREST = 0
  ILAST =0
  DO JA=1,KPROCA
    IF (JA  <=  KRESTM .OR. KRESTM  ==  0) THEN
      ICOMP = KMEDIAP
    ELSE
      ICOMP = KMEDIAP - 1
    ENDIF
    ITOT = IREST
    IGL = ILAST+1
    DO JGL=IGL,KDGL
      ILAST = JGL
      IF(ITOT+ICOST(JGL) < ICOMP) THEN
        ITOT = ITOT+ICOST(JGL)
      ELSEIF(ITOT+ICOST(JGL) == ICOMP) THEN
        IREST = 0
        KLAST(JA) = JGL
        KINDIC(JA) = 0
        EXIT
      ELSE
        IREST =  ICOST(JGL) -(ICOMP-ITOT)
        KLAST(JA) = JGL
        KINDIC(JA) = JGL
        EXIT
      ENDIF
    ENDDO
  ENDDO
  
ELSE

  ITOT_TOP=0
  ITOT_BOT=0
  IGL_TOP=1
  IGL_BOT=KDGL
  DO JA=1,(KPROCA-1)/2+1
    IF( JA /= KPROCA/2+1 )THEN
      LLDONE=.TRUE.
      DO WHILE ( LLDONE )
        IF( ITOT_TOP+ICOST(IGL_TOP) < KMEDIAP )THEN
          KLAST(JA)=IGL_TOP
          ITOT_TOP=ITOT_TOP+ICOST(IGL_TOP)
          IGL_TOP=IGL_TOP+1
        ELSE
          ITOT_TOP=ITOT_TOP-KMEDIAP
          LLDONE=.FALSE.
        ENDIF
      ENDDO
      KLAST(KPROCA-JA+1)=IGL_BOT
      LLDONE=.TRUE.
      DO WHILE ( LLDONE )
        IF( ITOT_BOT+ICOST(IGL_BOT) < KMEDIAP )THEN
          ITOT_BOT=ITOT_BOT+ICOST(IGL_BOT)
          IGL_BOT=IGL_BOT-1
        ELSE
          ITOT_BOT=ITOT_BOT-KMEDIAP
          LLDONE=.FALSE.
        ENDIF
      ENDDO
    ELSE
      KLAST(JA)=IGL_BOT
    ENDIF
  ENDDO

  LLSIMPLE=.FALSE.
  DO JA=1,KPROCA
    IF( KLAST(JA)==0 )THEN
      LLSIMPLE=.TRUE.
      EXIT
    ENDIF
  ENDDO
  IF( LLSIMPLE )THEN
!   WRITE(0,'("SUMPLATB_MOD: REVERTING TO SIMPLE LATITUDE DISTRIBUTION")')
    ILATS(:)=0
    IA=0
    DO JGL=1,KDGL
     IA=IA+1
     ILATS(IA)=ILATS(IA)+1
     IF( IA==KPROCA ) IA=0
    ENDDO
    KLAST(1)=ILATS(1)
    DO JA=2,KPROCA
      KLAST(JA)=KLAST(JA-1)+ILATS(JA)
    ENDDO
  ENDIF

ENDIF
  
END SUBROUTINE SUMPLATB
END MODULE SUMPLATB_MOD
