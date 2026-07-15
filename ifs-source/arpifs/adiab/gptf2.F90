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
 & PGMV,PGMVS,PGFL, &
 & P0SP  , P0SPL , P0SPM , P9SP  , P9SPL , P9SPM, &
 & P0DIV , P0NHX , P0SPD , P0SPDL, P0SPDM, P0SVD , P0SVDL, &
 & P0SVDM, P0T   , P0TL  , P0TM  , P0U   , P0V   , P9DIV,  &
 & P9NHX , P9SPD , P9SPDL, P9SPDM, P9SVD , P9SVDL, P9SVDM, &
 & P9T   , P9TL  , P9TM  , P9U   , P9V)

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
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA)   :: P0SP  , P0SPL , P0SPM , P9SP  , P9SPL , P9SPM 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P0DIV , P0NHX , P0SPD , P0SPDL, P0SPDM, P0SVD , P0SVDL
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P0SVDM, P0T   , P0TL  , P0TM  , P0U   , P0V   , P9DIV 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P9NHX , P9SPD , P9SPDL, P9SPDM, P9SVD , P9SVDL, P9SVDM
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P9T   , P9TL  , P9TM  , P9U   , P9V   

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL,JGFL,JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:)   :: Z0SP  , Z0SPL , Z0SPM , Z9SP  , Z9SPL , Z9SPM 
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0DIV , Z0NHX , Z0SPD , Z0SPDL, Z0SPDM, Z0SVD , Z0SVDL
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0SVDM, Z0T   , Z0TL  , Z0TM  , Z0U   , Z0V   , Z9DIV 
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z9NHX , Z9SPD , Z9SPDL, Z9SPDM, Z9SVD , Z9SVDL, Z9SVDM
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z9T   , Z9TL  , Z9TM  , Z9U   , Z9V   

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

Z0SP    => NULL ()
Z0SPL   => NULL ()
Z0SPM   => NULL ()
Z9SP    => NULL ()
Z9SPL   => NULL ()
Z9SPM   => NULL ()

Z0DIV   => NULL ()
Z0NHX   => NULL ()
Z0SPD   => NULL ()
Z0SPDL  => NULL ()
Z0SPDM  => NULL ()
Z0SVD   => NULL ()
Z0SVDL  => NULL ()
Z0SVDM  => NULL ()
Z0T     => NULL ()
Z0TL    => NULL ()
Z0TM    => NULL ()
Z0U     => NULL ()
Z0V     => NULL ()
Z9DIV   => NULL ()
Z9NHX   => NULL ()
Z9SPD   => NULL ()
Z9SPDL  => NULL ()
Z9SPDM  => NULL ()
Z9SVD   => NULL ()
Z9SVDL  => NULL ()
Z9SVDM  => NULL ()
Z9T     => NULL ()
Z9TL    => NULL ()
Z9TM    => NULL ()
Z9U     => NULL ()
Z9V     => NULL ()
           
IF (PRESENT (PGMVS)) THEN           
  Z0SP  => PGMVS(:,YT0%MSP)                       
  Z0SPL => PGMVS(:,YT0%MSPL)                       
  Z0SPM => PGMVS(:,YT0%MSPM)                       
  Z9SP  => PGMVS(:,YT9%MSP)                       
  Z9SPL => PGMVS(:,YT9%MSPL)                       
  Z9SPM => PGMVS(:,YT9%MSPM)                       
ENDIF

IF (PRESENT (P0SP )) Z0SP  => P0SP  
IF (PRESENT (P0SPL)) Z0SPL => P0SPL 
IF (PRESENT (P0SPM)) Z0SPM => P0SPM 
IF (PRESENT (P9SP )) Z9SP  => P9SP  
IF (PRESENT (P9SPL)) Z9SPL => P9SPL 
IF (PRESENT (P9SPM)) Z9SPM => P9SPM 

IF (PRESENT (PGMV)) THEN
  Z0DIV  => PGMV(:,:,YT0%MDIV)                       
  Z0NHX  => PGMV(:,:,YT0%MNHX)                       
  Z0SPD  => PGMV(:,:,YT0%MSPD)                       
  Z0SPDL => PGMV(:,:,YT0%MSPDL)                       
  Z0SPDM => PGMV(:,:,YT0%MSPDM)                       
  Z0SVD  => PGMV(:,:,YT0%MSVD)                       
  Z0SVDL => PGMV(:,:,YT0%MSVDL)                       
  Z0SVDM => PGMV(:,:,YT0%MSVDM)                       
  Z0T    => PGMV(:,:,YT0%MT)                       
  Z0TL   => PGMV(:,:,YT0%MTL)                       
  Z0TM   => PGMV(:,:,YT0%MTM)                       
  Z0U    => PGMV(:,:,YT0%MU)                       
  Z0V    => PGMV(:,:,YT0%MV)                       
  Z9DIV  => PGMV(:,:,YT9%MDIV)                       
  Z9NHX  => PGMV(:,:,YT9%MNHX)                       
  Z9SPD  => PGMV(:,:,YT9%MSPD)                       
  Z9SPDL => PGMV(:,:,YT9%MSPDL)                       
  Z9SPDM => PGMV(:,:,YT9%MSPDM)                       
  Z9SVD  => PGMV(:,:,YT9%MSVD)                       
  Z9SVDL => PGMV(:,:,YT9%MSVDL)                       
  Z9SVDM => PGMV(:,:,YT9%MSVDM)                       
  Z9T    => PGMV(:,:,YT9%MT)                       
  Z9TL   => PGMV(:,:,YT9%MTL)                       
  Z9TM   => PGMV(:,:,YT9%MTM)                       
  Z9U    => PGMV(:,:,YT9%MU)                       
  Z9V    => PGMV(:,:,YT9%MV)                       
ENDIF

IF (PRESENT (P0DIV )) Z0DIV  => P0DIV  
IF (PRESENT (P0NHX )) Z0NHX  => P0NHX  
IF (PRESENT (P0SPD )) Z0SPD  => P0SPD  
IF (PRESENT (P0SPDL)) Z0SPDL => P0SPDL 
IF (PRESENT (P0SPDM)) Z0SPDM => P0SPDM 
IF (PRESENT (P0SVD )) Z0SVD  => P0SVD  
IF (PRESENT (P0SVDL)) Z0SVDL => P0SVDL 
IF (PRESENT (P0SVDM)) Z0SVDM => P0SVDM 
IF (PRESENT (P0T   )) Z0T    => P0T    
IF (PRESENT (P0TL  )) Z0TL   => P0TL   
IF (PRESENT (P0TM  )) Z0TM   => P0TM   
IF (PRESENT (P0U   )) Z0U    => P0U    
IF (PRESENT (P0V   )) Z0V    => P0V    
IF (PRESENT (P9DIV )) Z9DIV  => P9DIV  
IF (PRESENT (P9NHX )) Z9NHX  => P9NHX  
IF (PRESENT (P9SPD )) Z9SPD  => P9SPD  
IF (PRESENT (P9SPDL)) Z9SPDL => P9SPDL 
IF (PRESENT (P9SPDM)) Z9SPDM => P9SPDM 
IF (PRESENT (P9SVD )) Z9SVD  => P9SVD  
IF (PRESENT (P9SVDL)) Z9SVDL => P9SVDL 
IF (PRESENT (P9SVDM)) Z9SVDM => P9SVDM 
IF (PRESENT (P9T   )) Z9T    => P9T    
IF (PRESENT (P9TL  )) Z9TL   => P9TL   
IF (PRESENT (P9TM  )) Z9TM   => P9TM   
IF (PRESENT (P9U   )) Z9U    => P9U    
IF (PRESENT (P9V   )) Z9V    => P9V    
           
!*       1. PERFORM TIME FILTER (PART 2)           
!           ----------------------------           
           
!*        1.1   TWO TIME LEVEL           
           
IF (YDDYNA%LTWOTL) THEN

  IF(LDFSTEP) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        Z9U(JL,JK) = Z0U(JL,JK)
        Z9V(JL,JK) = Z0V(JL,JK)
      ENDDO
    ENDDO
  ENDIF

ELSE

!*        1.2   THREE TIME LEVEL

  IF(LDFSTEP) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        Z9U(JL,JK)   = Z0U(JL,JK)
        Z9V(JL,JK)   = Z0V(JL,JK)
        Z9DIV(JL,JK) = Z0DIV(JL,JK)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          Z9T(JL,JK)  = Z0T(JL,JK)
          Z9TL(JL,JK) = Z0TL(JL,JK)
          Z9TM(JL,JK) = Z0TM(JL,JK)
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
          Z9SPD(JL,JK) = Z0SPD(JL,JK)
          Z9SVD(JL,JK) = Z0SVD(JL,JK)
          Z9SPDL(JL,JK) = Z0SPDL(JL,JK)
          Z9SVDL(JL,JK) = Z0SVDL(JL,JK)
          Z9SPDM(JL,JK) = Z0SPDM(JL,JK)
          Z9SVDM(JL,JK) = Z0SVDM(JL,JK)
        ENDDO
      ENDDO

      IF( YDDYNA%NVDVAR==4 .OR. YDDYNA%NVDVAR==5 ) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            Z9NHX(JL,JK) = Z0NHX(JL,JK)
          ENDDO
        ENDDO
      ENDIF

    ENDIF

    DO JL=KST,KEN
      Z9SP(JL)  = Z0SP(JL)
      Z9SPL(JL) = Z0SPL(JL)
      Z9SPM(JL) = Z0SPM(JL)
    ENDDO

  ELSE

    DO JK=1,NFLEVG
      DO JL=KST,KEN
        Z9U(JL,JK)   = Z9U(JL,JK)  +REPS2*Z0U(JL,JK)
        Z9V(JL,JK)   = Z9V(JL,JK)  +REPS2*Z0V(JL,JK)
        Z9DIV(JL,JK) = Z9DIV(JL,JK)+REPS2*Z0DIV(JL,JK)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          Z9T(JL,JK)  = Z9T(JL,JK)  +REPS2*Z0T(JL,JK)
          Z9TL(JL,JK) = Z9TL(JL,JK) +REPS2*Z0TL(JL,JK)
          Z9TM(JL,JK) = Z9TM(JL,JK) +REPS2*Z0TM(JL,JK)
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
          Z9SPD(JL,JK)  =&
           & Z9SPD(JL,JK)+REPS2*Z0SPD(JL,JK)
          Z9SVD(JL,JK)  =&
           & Z9SVD(JL,JK)+REPS2*Z0SVD(JL,JK)
          Z9SPDL(JL,JK) =&
           & Z9SPDL(JL,JK)+REPS2*Z0SPDL(JL,JK)
          Z9SVDL(JL,JK) =&
           & Z9SVDL(JL,JK)+REPS2*Z0SVDL(JL,JK)
          Z9SPDM(JL,JK) =&
           & Z9SPDM(JL,JK)+REPS2*Z0SPDM(JL,JK)
          Z9SVDM(JL,JK) =&
           & Z9SVDM(JL,JK)+REPS2*Z0SVDM(JL,JK)
        ENDDO
      ENDDO
    ENDIF

    DO JL=KST,KEN
      Z9SP(JL)  = Z9SP(JL) +REPS2*Z0SP(JL)
      Z9SPL(JL) = Z9SPL(JL)+REPS2*Z0SPL(JL)
      Z9SPM(JL) = Z9SPM(JL)+REPS2*Z0SPM(JL)
    ENDDO

  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF2',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF2
