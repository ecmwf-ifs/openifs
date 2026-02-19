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

SUBROUTINE GPTF1(YDGEOMETRY,YDGMV,&
 ! --- INPUT ---------------------------------------------------------
 & YDML_GCONF,YDDYN,YDDYNA,LDTF1,KST,KEN,LDLSTEP,&
 ! --- INPUT/OUTPUT --------------------------------------------------
 & PGMV,PGMVS,PGFL)

!**** *GPTF1* - Timefilter part 1

!     Purpose.
!     --------
!           Performs part 1 of the time-filtering for the t-dt array.
!           Applies to variables and first-order derivatives.

!           - leap-frog and ldlstep=false: px9=eps1*px9+(1-eps1-eps2)*px0
!           - leap-frog and ldlstep=true : px9=0
!           - sl2tl                      : px9=px0

!**   Interface.
!     ----------
!        *CALL* *GPTF1(...)

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

! Modifications
! -------------
!   Modified 15-10-01 D.Salmond FULLIMP mods
!   J.Vivoda 03-2002 PC schemes for NH dynamics (LPC_XXXX keys)
!   Modified 30-09-02 P.Smolikova : variable d4 in NH
!   Modified 13-11-02 K. YESSAD : cleanings + improve vectorization.
!   Modified 01-04-03 M.Hamrud : Revised data flow/GFL
!   Modified 2003-07-17 C. Fischer - psvdauxt0* come from pgmv
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   08-Jun-2004 J. Masek   NH cleaning (LVSLWBC, LFULLIMP)
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (July 2014): Move some variables.
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
LOGICAL           ,INTENT(IN)    :: LDTF1
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
LOGICAL           ,INTENT(IN)    :: LDLSTEP 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZREST, ZRESTM
INTEGER(KIND=JPIM) :: JGFL,JL,JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPTF1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & REPS1=>YDDYN%REPS1, REPS2=>YDDYN%REPS2, REPSM1=>YDDYN%REPSM1, &
 & REPSM2=>YDDYN%REPSM2, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)
!     ------------------------------------------------------------------

!*       1. PERFORM TIME FILTER (PART 1)
!           ----------------------------

IF( LDTF1 ) THEN
  
  IF (YDDYNA%LTWOTL) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        PGMV(JL,JK,YT9%MU)=PGMV(JL,JK,YT0%MU)
        PGMV(JL,JK,YT9%MV)=PGMV(JL,JK,YT0%MV)
      ENDDO

      IF(YDDYNA%LGRADSP)THEN
        DO JL=KST,KEN
          PGMV(JL,JK,YT9%MSGRTL)=PGMV(JL,JK,YT0%MSGRTL)
          PGMV(JL,JK,YT9%MSGRTM)=PGMV(JL,JK,YT0%MSGRTM)
        ENDDO
      ENDIF

      IF(YDDYNA%LPC_FULL)THEN
        DO JL=KST,KEN
          PGMV(JL,JK,YT9%MT)  = PGMV(JL,JK,YT0%MT)
          PGMV(JL,JK,YT9%MTL) = PGMV(JL,JK,YT0%MTL)
          PGMV(JL,JK,YT9%MTM) = PGMV(JL,JK,YT0%MTM)
        ENDDO
      ENDIF
      IF(YDDYNA%LPC_FULL)THEN
        DO JGFL=1,NUMFLDS
          IF(YCOMP(JGFL)%LT9) THEN
            DO JL=KST,KEN
              PGFL(JL,JK,YCOMP(JGFL)%MP9) = PGFL(JL,JK,YCOMP(JGFL)%MP)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    IF(YDDYNA%LPC_FULL)THEN

      IF(YDDYNA%LNHDYN) THEN
        PGMV(KST:KEN,1:NFLEVG,YT9%MSPD) = PGMV(KST:KEN,1:NFLEVG,YT0%MSPD)
        PGMV(KST:KEN,1:NFLEVG,YT9%MSVD) = PGMV(KST:KEN,1:NFLEVG,YT0%MSVD)
        PGMV(KST:KEN,1:NFLEVG,YT9%MSPDL) =PGMV(KST:KEN,1:NFLEVG,YT0%MSPDL)
        PGMV(KST:KEN,1:NFLEVG,YT9%MSVDL) =PGMV(KST:KEN,1:NFLEVG,YT0%MSVDL)
        PGMV(KST:KEN,1:NFLEVG,YT9%MSPDM) =PGMV(KST:KEN,1:NFLEVG,YT0%MSPDM)
        PGMV(KST:KEN,1:NFLEVG,YT9%MSVDM) =PGMV(KST:KEN,1:NFLEVG,YT0%MSVDM)
      ENDIF

      DO JL=KST,KEN
        PGMVS(JL,YT9%MSP)  = PGMVS(JL,YT0%MSP)
        PGMVS(JL,YT9%MSPL) = PGMVS(JL,YT0%MSPL)
        PGMVS(JL,YT9%MSPM) = PGMVS(JL,YT0%MSPM)
      ENDDO

    ENDIF

  ELSE

    ZREST  = 1.0_JPRB-(REPS1+REPS2)
    ZRESTM = 1.0_JPRB-(REPSM1+REPSM2)
    IF(.NOT.LDLSTEP) THEN
      DO JK=1,NFLEVG
        DO JL=KST,KEN
          PGMV(JL,JK,YT9%MU)=REPS1*PGMV(JL,JK,YT9%MU)+ZREST*PGMV(JL,JK,YT0%MU)
          PGMV(JL,JK,YT9%MV)=REPS1*PGMV(JL,JK,YT9%MV)+ZREST*PGMV(JL,JK,YT0%MV)
          PGMV(JL,JK,YT9%MDIV)=&
           & REPS1*PGMV(JL,JK,YT9%MDIV)+ZREST*PGMV(JL,JK,YT0%MDIV)
        ENDDO
        IF(NFTHER >= 1) THEN
          DO JL=KST,KEN
            PGMV(JL,JK,YT9%MT)  =&
             & REPS1*PGMV(JL,JK,YT9%MT) +ZREST*PGMV(JL,JK,YT0%MT)
            PGMV(JL,JK,YT9%MTL) =&
             & REPS1*PGMV(JL,JK,YT9%MTL)+ZREST*PGMV(JL,JK,YT0%MTL)
            PGMV(JL,JK,YT9%MTM) =&
             & REPS1*PGMV(JL,JK,YT9%MTM)+ZREST*PGMV(JL,JK,YT0%MTM)
          ENDDO
        ENDIF
        IF(YDDYNA%LNHDYN) THEN
          DO JL=KST,KEN
            PGMV(JL,JK,YT9%MSPD) =&
             & REPS1*PGMV(JL,JK,YT9%MSPD)+ZREST*PGMV(JL,JK,YT0%MSPD)
            PGMV(JL,JK,YT9%MSVD) =&
             & REPS1*PGMV(JL,JK,YT9%MSVD)+ZREST*PGMV(JL,JK,YT0%MSVD)
            PGMV(JL,JK,YT9%MSPDL) =&
             & REPS1*PGMV(JL,JK,YT9%MSPDL)+ZREST*PGMV(JL,JK,YT0%MSPDL)
            PGMV(JL,JK,YT9%MSVDL) =&
             & REPS1*PGMV(JL,JK,YT9%MSVDL)+ZREST*PGMV(JL,JK,YT0%MSVDL)
            PGMV(JL,JK,YT9%MSPDM) =&
             & REPS1*PGMV(JL,JK,YT9%MSPDM)+ZREST*PGMV(JL,JK,YT0%MSPDM)
            PGMV(JL,JK,YT9%MSVDM) =&
             & REPS1*PGMV(JL,JK,YT9%MSVDM)+ZREST*PGMV(JL,JK,YT0%MSVDM)
          ENDDO
        ENDIF
      ENDDO

      DO JGFL=1,NUMFLDS
        IF(YCOMP(JGFL)%LT9) THEN
          DO JK=1,NFLEVG
            DO JL=KST,KEN
              PGFL(JL,JK,YCOMP(JGFL)%MP9) = REPSM1*PGFL(JL,JK,YCOMP(JGFL)%MP9)+&
               & ZRESTM*PGFL(JL,JK,YCOMP(JGFL)%MP)  
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      DO JL=KST,KEN
        PGMVS(JL,YT9%MSP)  = REPS1*PGMVS(JL,YT9%MSP) +ZREST*PGMVS(JL,YT0%MSP)
        PGMVS(JL,YT9%MSPL) = REPS1*PGMVS(JL,YT9%MSPL)+ZREST*PGMVS(JL,YT0%MSPL)
        PGMVS(JL,YT9%MSPM) = REPS1*PGMVS(JL,YT9%MSPM)+ZREST*PGMVS(JL,YT0%MSPM)
      ENDDO

    ELSE
      PGMV(:,:,YT0%NDIM+1:) = 0.0_JPRB
      PGMVS(:,YT0%NDIMS+1:) = 0.0_JPRB
      DO JGFL=1,NUMFLDS
        IF(YCOMP(JGFL)%LT9) THEN
          DO JK=1,NFLEVG
            DO JL=KST,KEN
              PGFL(JL,JK,YCOMP(JGFL)%MP9) = 0.0_JPRB
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF
ELSE

  IF(YDDYNA%LGRADSP)THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        PGMV(JL,JK,YT9%MSGRTL)=PGMV(JL,JK,YT0%MSGRTL)
        PGMV(JL,JK,YT9%MSGRTM)=PGMV(JL,JK,YT0%MSGRTM)
      ENDDO
    ENDDO
  ENDIF
  
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF1',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF1
