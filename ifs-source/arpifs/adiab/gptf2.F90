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

SUBROUTINE GPTF2(YDGEOMETRY,YDGMV,&
 ! --- INPUT ---------------------------------------------------------
 & YDML_GCONF,YDDYN,YDDYNA,KST,KEN,LDFSTEP,&
 ! --- INPUT-OUTPUT --------------------------------------------------
 & PGMV,PGMVS,PGFL)

!**** *GPTF2* - Timefilter part 2

!     Purpose.
!     --------
!           Performs part 2 of the time-filtering.
!           - leap frog + ldfstep=true : pxt9. = pxt0.
!           - leap frog + ldfstep=false: pxt9. = pxt9. + eps2*pxt0.
!           - sl2tl     + ldfstep=true : put9=put0 and pvt9=pvt0 only.
!           - sl2tl     + ldfstep=false: nothing.

!**   Interface.
!     ----------
!        *CALL* *GPTF2(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!         KST          : start of work
!         KEN          : depth of work
!         LDLSTEP      : check on the last time step.

!        INPUT/OUTPUT:
!         PGMV         : "t" and "t-dt" upper air GMV variables.
!         PGMVS        : "t" and "t-dt" surface GMV variables.
!         PGFL         : "t" and "t-dt" GFL variables.

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 92-02-01

! Modifications.
! --------------
!   Modified 02-09-30 by P. Smolikova (variable d4 in NH)
!   Modified 13-11-02 K. YESSAD : cleanings + improve vectorization.
!   Modified 2003-07-17 C. Fischer - psvdauxt0* come from pgmv
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   08-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
!   Modified Nov 2007 N. Wedi: bug correction
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
! End Modifications
!------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
LOGICAL           ,INTENT(IN)    :: LDFSTEP 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL,JGFL,JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPTF2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & REPS2=>YDDYN%REPS2, REPSM2=>YDDYN%REPSM2, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)
!     ------------------------------------------------------------------

!*       1. PERFORM TIME FILTER (PART 2)
!           ----------------------------

!*        1.1   TWO TIME LEVEL

IF (YDDYNA%LTWOTL) THEN

  IF(LDFSTEP) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        PGMV(JL,JK,YT9%MU) = PGMV(JL,JK,YT0%MU)
        PGMV(JL,JK,YT9%MV) = PGMV(JL,JK,YT0%MV)
      ENDDO
    ENDDO
  ENDIF

ELSE

!*        1.2   THREE TIME LEVEL

  IF(LDFSTEP) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        PGMV(JL,JK,YT9%MU)   = PGMV(JL,JK,YT0%MU)
        PGMV(JL,JK,YT9%MV)   = PGMV(JL,JK,YT0%MV)
        PGMV(JL,JK,YT9%MDIV) = PGMV(JL,JK,YT0%MDIV)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          PGMV(JL,JK,YT9%MT)  = PGMV(JL,JK,YT0%MT)
          PGMV(JL,JK,YT9%MTL) = PGMV(JL,JK,YT0%MTL)
          PGMV(JL,JK,YT9%MTM) = PGMV(JL,JK,YT0%MTM)
        ENDDO
      ENDIF
    ENDDO

    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LT9) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            PGFL(JL,JK,YCOMP(JGFL)%MP9) = PGFL(JL,JK,YCOMP(JGFL)%MP)
          ENDDO
        ENDDO
      ENDIF
    ENDDO

    IF(YDDYNA%LNHDYN) THEN
      DO JK=1,NFLEVG
        DO JL=KST,KEN
          PGMV(JL,JK,YT9%MSPD) = PGMV(JL,JK,YT0%MSPD)
          PGMV(JL,JK,YT9%MSVD) = PGMV(JL,JK,YT0%MSVD)
          PGMV(JL,JK,YT9%MSPDL) = PGMV(JL,JK,YT0%MSPDL)
          PGMV(JL,JK,YT9%MSVDL) = PGMV(JL,JK,YT0%MSVDL)
          PGMV(JL,JK,YT9%MSPDM) = PGMV(JL,JK,YT0%MSPDM)
          PGMV(JL,JK,YT9%MSVDM) = PGMV(JL,JK,YT0%MSVDM)
        ENDDO
      ENDDO

      IF( YDDYNA%NVDVAR==4 ) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            PGMV(JL,JK,YT9%MNHX) = PGMV(JL,JK,YT0%MNHX)
          ENDDO
        ENDDO
      ENDIF

    ENDIF

    DO JL=KST,KEN
      PGMVS(JL,YT9%MSP)  = PGMVS(JL,YT0%MSP)
      PGMVS(JL,YT9%MSPL) = PGMVS(JL,YT0%MSPL)
      PGMVS(JL,YT9%MSPM) = PGMVS(JL,YT0%MSPM)
    ENDDO

  ELSE

    DO JK=1,NFLEVG
      DO JL=KST,KEN
        PGMV(JL,JK,YT9%MU)   = PGMV(JL,JK,YT9%MU)  +REPS2*PGMV(JL,JK,YT0%MU)
        PGMV(JL,JK,YT9%MV)   = PGMV(JL,JK,YT9%MV)  +REPS2*PGMV(JL,JK,YT0%MV)
        PGMV(JL,JK,YT9%MDIV) = PGMV(JL,JK,YT9%MDIV)+REPS2*PGMV(JL,JK,YT0%MDIV)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          PGMV(JL,JK,YT9%MT)  = PGMV(JL,JK,YT9%MT)  +REPS2*PGMV(JL,JK,YT0%MT)
          PGMV(JL,JK,YT9%MTL) = PGMV(JL,JK,YT9%MTL) +REPS2*PGMV(JL,JK,YT0%MTL)
          PGMV(JL,JK,YT9%MTM) = PGMV(JL,JK,YT9%MTM) +REPS2*PGMV(JL,JK,YT0%MTM)
        ENDDO
      ENDIF
    ENDDO

    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LT9) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            PGFL(JL,JK,YCOMP(JGFL)%MP9) = PGFL(JL,JK,YCOMP(JGFL)%MP9)+&
             & REPSM2*PGFL(JL,JK,YCOMP(JGFL)%MP)  
          ENDDO
        ENDDO
      ENDIF
    ENDDO

    IF(YDDYNA%LNHDYN) THEN
      DO JK=1,NFLEVG
        DO JL=KST,KEN
          PGMV(JL,JK,YT9%MSPD)  =&
           & PGMV(JL,JK,YT9%MSPD)+REPS2*PGMV(JL,JK,YT0%MSPD)
          PGMV(JL,JK,YT9%MSVD)  =&
           & PGMV(JL,JK,YT9%MSVD)+REPS2*PGMV(JL,JK,YT0%MSVD)
          PGMV(JL,JK,YT9%MSPDL) =&
           & PGMV(JL,JK,YT9%MSPDL)+REPS2*PGMV(JL,JK,YT0%MSPDL)
          PGMV(JL,JK,YT9%MSVDL) =&
           & PGMV(JL,JK,YT9%MSVDL)+REPS2*PGMV(JL,JK,YT0%MSVDL)
          PGMV(JL,JK,YT9%MSPDM) =&
           & PGMV(JL,JK,YT9%MSPDM)+REPS2*PGMV(JL,JK,YT0%MSPDM)
          PGMV(JL,JK,YT9%MSVDM) =&
           & PGMV(JL,JK,YT9%MSVDM)+REPS2*PGMV(JL,JK,YT0%MSVDM)
        ENDDO
      ENDDO
    ENDIF

    DO JL=KST,KEN
      PGMVS(JL,YT9%MSP)  = PGMVS(JL,YT9%MSP) +REPS2*PGMVS(JL,YT0%MSP)
      PGMVS(JL,YT9%MSPL) = PGMVS(JL,YT9%MSPL)+REPS2*PGMVS(JL,YT0%MSPL)
      PGMVS(JL,YT9%MSPM) = PGMVS(JL,YT9%MSPM)+REPS2*PGMVS(JL,YT0%MSPM)
    ENDDO

  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF2',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF2
