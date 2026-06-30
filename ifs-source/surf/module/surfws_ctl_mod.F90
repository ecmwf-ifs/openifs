MODULE SURFWS_CTL_MOD
CONTAINS
SUBROUTINE SURFWS_CTL( KIDIA, KFDIA, KLON, KLEVSN,  &
                     & PSDOR, &
                     & LSMASK, PFRTI,PMU0,          &
                     & PTSA, PTSKIN, PALBSN,        &
                     & PTSN, PSSN, PRSN, PWSN,      &
                     & YDCST, YDSOIL                )

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_SOIL , ONLY : TSOIL

USE SRFSN_VGRID_MOD
USE SRFSN_REGRID_MOD
USE SURFWS_INIT_ML_MOD
USE SURFWS_INIT_SL_MOD
USE SURFWS_INIT_MLOFF_MOD
USE SURFWS_FGPROF_MOD
USE SURFWS_MASSADJ_MOD
USE SURFWS_TSNADJ_MOD
!USE YOMLUN   , ONLY : NULOUT

USE ABORT_SURF_MOD

! (C) Copyright 2017- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

                   !**** *SURFWS_CTL* - Snow warm start multi-layer 
!     PURPOSE.
!     --------
!          THIS ROUTINE CONTROLS THE WARM START OF MULTI-LAYER SCHEME IN
!          FC MODE. 
!          It initialises/allocate basic fields needed by the parametrizations,
!          and calls routines to setup the profiles (surfws_fgprof) and adjust
!          the mass in the snow layer (surfws_massadj).

!**   INTERFACE.
!     ----------
!          *SURFWS_CTL* IS CALLED FROM *SURFWS* (external).

!     PARAMETER   DESCRIPTION               UNITS
!     ---------   -----------               -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       Length of arrays
!    *KLEVSN*     Snow vertical levels

!     INPUT PARAMETERS (REAL):
!    *PSDOR*      sub grid scale orography   (m)
!    *PFRTI*      tile fractions             (-)
!    *PMU0*       cos solar zenith angle     (-)

!     INPUT PARAMETERS (LOGICAL):
!    *LSMASK*     LAND/SEA MASK         (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1  (REAL):
!    *PTSA*       soil temperature t-1       (K)
!    *PTSKIN*     skin temperature t-1       (K)
!    *PALBSN*     albedo of snow   t-1       (0-1)


!     INPUT/OUTPUT PARAMETERS  (REAL):
!    *PTSN*       snow temperature t-1       (K)
!    *PSSN*       snow mass        t-1       (kg m-2)
!    *PRSN*       snow density     t-1       (kg m-3)
!    *PWSN*       snow liq water   t-1       (kg m-2)

!     METHOD.
!     -------
!     NSNMLWS=1_JPIM ! start from multi-layer top value for temp, avg density
!     NSNMLWS=2_JPIM ! start from single-layer values (e.g. ERA5)
!     NSNMLWS=3_JPIM ! start from multi-layer offline average temperature values

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          

!     Modifications:
!     Original   G. Arduini      ECMWF     28/07/2017


IMPLICIT NONE


! Input variables:
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA, KLON, KLEVSN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDOR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSA(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKIN(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBSN(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: LSMASK(:) 

! Output fields:
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSN(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSSN(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSN(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWSN(:,:) 


TYPE(TCST)     , INTENT(IN) :: YDCST
TYPE(TSOIL)    , INTENT(IN) :: YDSOIL


! Local of warm start Variables
REAL(KIND=JPRB)    :: ZACTDEPTH(KLON)
REAL(KIND=JPRB)    :: ZTHRESWS
REAL(KIND=JPRB)    :: ZSNPERT
REAL(KIND=JPRB)    :: ZDSN(KLON,KLEVSN),ZDSNREAL(KLON,KLEVSN)
REAL(KIND=JPRB)    :: ZSNDEPTH(KLON,KLEVSN)
REAL(KIND=JPRB)    :: ZSNDEPTHREAL(KLON,KLEVSN)

REAL(KIND=JPRB)    :: ZHSN(KLON)
REAL(KIND=JPRB)    :: ZDSNTOT(KLON)
REAL(KIND=JPRB)    :: ZDSNTOTREAL(KLON)
REAL(KIND=JPRB)    :: ZTSN(KLON)
REAL(KIND=JPRB)    :: ZRSN(KLON)
REAL(KIND=JPRB)    :: ZSSN(KLON)
REAL(KIND=JPRB)    :: ZWSN(KLON)

REAL(KIND=JPRB)    :: ZFRSN(KLON) ! snow cover fraction
REAL(KIND=JPRB)    :: ZLEVMAX(KLEVSN),ZLEVMIN(KLEVSN)
REAL(KIND=JPRB)    :: ZLEVMAX_GL(KLEVSN),ZLEVMIN_GL(KLEVSN)
REAL(KIND=JPRB)    :: ZLEVMINA(KLEVSN)
LOGICAL            :: LLNOSNOW(KLON)  ! FALSE to compute snow
INTEGER(KIND=JPRB) :: JL,JK !!,KSNACC

REAL(KIND=JPRB)    :: ZTSNTOP(KLON) 
REAL(KIND=JPRB)    :: ZTSNBOTTOM(KLON)
REAL(KIND=JPRB)    :: ZTSNMIDDLE(KLON)
REAL(KIND=JPRB)    :: ZRSNTOP(KLON)
REAL(KIND=JPRB)    :: ZRSNMAX(KLON)
REAL(KIND=JPRB)    :: ZSADEPTH(KLON)
INTEGER(KIND=JPIM) :: KLEVSNA(KLON)
INTEGER(KIND=JPIM) :: KLEVMID(KLON)
INTEGER(KIND=JPIM) :: ZTMINCL(KLON) 
INTEGER(KIND=JPIM) :: ZRMINCL(KLON)  

LOGICAL            :: LDLAND(KLON) 

! Cluster number
INTEGER(KIND=JPIM)              :: INCL

! Local Warm start arrays
!REAL(KIND=JPRB)    :: ZTSNWS(KLON,0:KLEVSN+1)
REAL(KIND=JPRB)    :: ZTSNWS(KLON,KLEVSN)
REAL(KIND=JPRB)    :: ZRSNWS(KLON,KLEVSN)
REAL(KIND=JPRB)    :: ZSSNWS(KLON,KLEVSN)
REAL(KIND=JPRB)    :: ZWSNWS(KLON,KLEVSN)


! Allocatable
REAL(KIND=JPRB), ALLOCATABLE    :: ZTCONSTAVG(:,:), ZTCONSTSTD(:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZRCONSTAVG(:,:), ZRCONSTSTD(:,:) 


REAL(KIND=JPRB)    :: ZIHCAP

!INTEGER(KIND=JPIM) :: NSNMLWS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_CTL_MOD:SURFWS_CTL',0,ZHOOK_HANDLE)


ASSOCIATE(RTT=>YDCST%RTT,RLMLT=>YDCST%RLMLT, RPI=>YDCST%RPI,    &
        & RLWCSWEA=>YDSOIL%RLWCSWEA, RLWCSWEB=>YDSOIL%RLWCSWEB, &
        & RLWCSWEC=>YDSOIL%RLWCSWEC, RTEMPAMP=>YDSOIL%RTEMPAMP, &
        & RDSNMAX=>YDSOIL%RDSNMAX, RHOMINSND=>YDSOIL%RHOMINSND, &
        & RHOMAXSN_NEW=>YDSOIL%RHOMAXSN_NEW, RDAT=>YDSOIL%RDAT, &
        & NSNMLWS=>YDSOIL%NSNMLWS)


! Define constants 
!****************** 
! 0. Define constants 

! Number of cluster used in k-means clustering
INCL = 27_JPIM 
! snow depth threshold to compute profiles
ZTHRESWS = 0.125_JPRB

! Ice Heat capacity
ZIHCAP     = YDSOIL%RHOCI / YDSOIL%RHOICE

! Glacier mask
ZSNPERT    = 9000.0_JPRB ! permanent snow threshold

! Level of min and max
ZLEVMAX         = YDSOIL%RLEVSNMAX
ZLEVMIN         = YDSOIL%RLEVSNMIN
ZLEVMAX_GL         = YDSOIL%RLEVSNMAX_GL
ZLEVMIN_GL         = YDSOIL%RLEVSNMIN_GL

!*************************
! 1.0 Setup of the routine
!*************************

! Allocate constant for warm start profiles
ALLOCATE( ZTCONSTAVG(KLON,INCL) )
ALLOCATE( ZTCONSTSTD(KLON,INCL) )
ALLOCATE( ZRCONSTAVG(KLON,INCL) )
ALLOCATE( ZRCONSTSTD(KLON,INCL) )

!**********************
! 1.1 Define snow points
DO JL=KIDIA,KFDIA
  ! SNOW COVER FRACTION
  ZFRSN(JL)=MAX(PFRTI(JL,5)+PFRTI(JL,7),YDSOIL%RFRTINY)

! POINTS WITH SNOW
  IF (ZFRSN(JL) < YDSOIL%RFRSMALL) THEN
    LLNOSNOW(JL)=.TRUE.
  ELSE
    LLNOSNOW(JL)=.FALSE.
  ENDIF

  IF (LSMASK(JL) == 1._JPRB) THEN
    LDLAND(JL) = .TRUE.
  ELSE
    LDLAND(JL) = .FALSE.
  ENDIF
ENDDO

!******************************************************************
! 1.2 Define snow vertical grid and mean quantities for computation
CALL SRFSN_VGRID(KIDIA,KFDIA,KLON,KLEVSN, LLNOSNOW, &
               & PSDOR,                             &
               & PSSN,PRSN,                         &
               & ZLEVMIN,ZLEVMAX, &
               & ZLEVMIN_GL,ZLEVMAX_GL, &
               & ZDSN,KLEVSNA)

DO JL=KIDIA,KFDIA
  DO JK=1, KLEVSN
    ZDSNREAL(JL,JK)=ZDSN(JL,JK)
    ZDSN(JL,JK)=MIN( 1._JPRB, ZDSN(JL,JK) )
  ENDDO
ENDDO
! INPUT LAYER HEAT CONTENT 
ZHSN(KIDIA:KFDIA) = ZIHCAP*(PSSN(KIDIA:KFDIA,1))*&
                      (PTSN(KIDIA:KFDIA,1)-RTT)-RLMLT*(PSSN(KIDIA:KFDIA,1)-PWSN(KIDIA:KFDIA,1))

! ZDSNTOT is the total depth, computed with ZDSN...
ZDSNTOT(KIDIA:KFDIA)     = SUM( ZDSN(KIDIA:KFDIA,:), DIM=2 )
ZDSNTOTREAL(KIDIA:KFDIA) = SUM( ZDSNREAL(KIDIA:KFDIA,:), DIM=2 )

! ZTSN is the mean temperature (Single-layer)
ZTSN(KIDIA:KFDIA) = PTSN(KIDIA:KFDIA, 1)

! ZRSN is the mean density (Single-layer)
ZRSN(KIDIA:KFDIA) = PRSN(KIDIA:KFDIA, 1)

! ZSSN is the total water eq (Single-layer, we define it for consistency in the
! routine
ZSSN(KIDIA:KFDIA) = SUM(PSSN(KIDIA:KFDIA, :), DIM=2)

! ZWSN is the Liquid water part (Single-layer, we define it for consistency in the
! routine
ZWSN(KIDIA:KFDIA)    = 0._JPRB 

!**! snow layer to absorbe residuals 
!**!!KSNACC          = MINLOC(YDSOIL%RLEVSNMAX, DIM=1 )
!**!!ZLEVMAX(KSNACC) = 1.E10_JPRB 
!**
!**! accumulated min snow depth 
!**ZLEVMINA(1)    = ZLEVMIN(1)
!**DO JK=2, KLEVSN
!**  ZLEVMINA(JK) = ZLEVMINA(JK-1) + ZLEVMIN(JK)
!**ENDDO
!**
!**! 1.3 Find # of active snow layer
!**DO JL=KIDIA,KFDIA
!**  IF ( LLNOSNOW(JL) ) THEN
!**    KLEVSNA(JL) = 1
!**  ELSE
!**    KLEVSNA(JL) = 1
!**    DO JK=2, KLEVSN
!**    ! IF ( ZDSNTOT(JL) < ZLEVMINA(JK) ) THEN
!**      IF ( ZDSNTOTREAL(JL) < ZLEVMINA(JK) ) THEN
!**        KLEVSNA(JL) = JK-1
!**        EXIT
!**      ELSE
!**        KLEVSNA(JL) = JK 
!**      ENDIF
!**    ENDDO
!**  ENDIF 
!**ENDDO 

! ZSNDEPTH is the depth of each layer with respect to zero.

IF (NSNMLWS == 3_JPIM) THEN
  ZSNDEPTH(KIDIA:KFDIA, :)     = 0.5_JPRB*ZDSN(KIDIA:KFDIA, :)
  ZSNDEPTHREAL(KIDIA:KFDIA, :) = 0.5_JPRB*ZDSNREAL(KIDIA:KFDIA, :)
ELSE
  ZSNDEPTH(KIDIA:KFDIA, 1)     = 0._JPRB
  ZSNDEPTHREAL(KIDIA:KFDIA, 1) = 0._JPRB
ENDIF
DO JK = 2, KLEVSN
  ZSNDEPTH(KIDIA:KFDIA, JK) = ZSNDEPTH(KIDIA:KFDIA, JK-1) + 0.5*(ZDSN(KIDIA:KFDIA,JK-1)+ZDSN(KIDIA:KFDIA,JK))
  ZSNDEPTHREAL(KIDIA:KFDIA, JK) = ZSNDEPTHREAL(KIDIA:KFDIA, JK-1) + &
                                 &0.5*(ZDSNREAL(KIDIA:KFDIA,JK-1)+ZDSNREAL(KIDIA:KFDIA,JK))
ENDDO


!------ SURFWS_INIT -----!
! Basic initialization warm start fields:
DO JL=KIDIA,KFDIA
    ZTSNWS(JL,1:KLEVSN)   = ZTSN(JL)
    ZRSNWS(JL,1:KLEVSN)   = ZRSN(JL)
    ZSSNWS(JL,1)          = ZSSN(JL)
    ZSSNWS(JL,2:KLEVSN)   = 0._JPRB
    ZWSNWS(JL,1:KLEVSN)   = 0._JPRB
ENDDO
ZTCONSTAVG=0._JPRB
ZTCONSTSTD=0._JPRB
ZRCONSTAVG=0._JPRB
ZRCONSTSTD=0._JPRB
ZTMINCL=1_JPIM
ZRMINCL=1_JPIM

! Here is where you want to know which profiles to reconstruct:
! ML exp or SL
!NSNMLWS=1_JPIM ! start from multi-layer top value for temp, avg density
!NSNMLWS=2_JPIM ! start from single-layer values
!NSNMLWS=3_JPIM ! start from multi-layer offline average temperature values
IF (NSNMLWS == 1_JPIM) THEN
  CALL SURFWS_INIT_ML(KIDIA, KFDIA, KLON, KLEVSN,INCL, PMU0,PSDOR,       & ! Input
               & PTSA(:,1), PTSKIN, &
               & ZDSNTOT, ZSNDEPTH,               &
               & ZSNPERT,                & ! Input
               & ZDSNREAL,ZTSN, ZRSN, ZSSN, ZWSN,PALBSN,            & ! Input
               & ZTSNWS,ZSSNWS,ZRSNWS,ZWSNWS,                       & ! Output
               & ZTSNTOP, ZTSNBOTTOM, ZTSNMIDDLE,                   & ! Output
               & ZRSNMAX, ZSADEPTH, KLEVSNA, KLEVMID, ZACTDEPTH,    & ! Output
               & ZTMINCL,ZRMINCL,                                   & ! Output
               & ZTCONSTAVG, ZTCONSTSTD,                            & ! Output
               & ZRCONSTAVG, ZRCONSTSTD, ZRSNTOP,                   & ! Output
               & YDCST, YDSOIL )

ELSE IF (NSNMLWS == 2_JPIM) THEN
  CALL SURFWS_INIT_SL(KIDIA, KFDIA, KLON, KLEVSN,INCL, PMU0,PSDOR,  & ! Input
               & PTSA(:,1), PTSKIN, &
               & ZDSNTOT, ZSNDEPTH,               &
               & ZSNPERT,                & ! Input
               & ZDSNREAL,ZTSN, ZRSN, ZSSN, ZWSN,PALBSN,            & ! Input
               & ZTSNWS,ZSSNWS,ZRSNWS,ZWSNWS,                       & ! Output
               & ZTSNTOP, ZTSNBOTTOM, ZTSNMIDDLE,                   & ! Output
               & ZRSNMAX, ZSADEPTH, KLEVSNA, KLEVMID, ZACTDEPTH,    & ! Output
               & ZTMINCL,ZRMINCL,                                   & ! Output
               & ZTCONSTAVG, ZTCONSTSTD,                            & ! Output
               & ZRCONSTAVG, ZRCONSTSTD, ZRSNTOP,                   & ! Output
               & YDCST, YDSOIL )

ELSE IF (NSNMLWS == 3_JPIM) THEN
  CALL SURFWS_INIT_MLOFF(KIDIA, KFDIA, KLON, KLEVSN,INCL, PMU0,PSDOR,       & ! Input
               & PTSA(:,1), PTSKIN, &
               & ZDSNTOT, ZSNDEPTH,               &
               & ZSNPERT,                & ! Input
               & ZDSNREAL,ZTSN, ZRSN, ZSSN, ZWSN,PALBSN,            & ! Input
               & ZTSNWS,ZSSNWS,ZRSNWS,ZWSNWS,                       & ! Output
               & ZTSNTOP, ZTSNBOTTOM, ZTSNMIDDLE,                   & ! Output
               & ZRSNMAX, ZSADEPTH, KLEVSNA, KLEVMID, ZACTDEPTH,    & ! Output
               & ZTMINCL,ZRMINCL,                                   & ! Output
               & ZTCONSTAVG, ZTCONSTSTD,                            & ! Output
               & ZRCONSTAVG, ZRCONSTSTD, ZRSNTOP,                   & ! Output
               & YDCST, YDSOIL )
ENDIF




CALL SURFWS_FGPROF(KIDIA, KFDIA, KLON, KLEVSN,               &
                 & KLEVSNA, KLEVMID,                         & 
                 & ZTHRESWS, ZSNPERT, LLNOSNOW,              &  ! THESE CAN BE MOVED TO SUSSOIL
                 & ZTSN, ZSSN, ZRSN, PTSA(:,1),              &
                 & ZDSNTOT,ZACTDEPTH, ZSNDEPTH,ZSNDEPTHREAL, ZSADEPTH,    &
                 & ZTSNBOTTOM, ZTSNTOP, ZTSNMIDDLE,          &
                 & ZRSNTOP, ZRSNMAX,                         &
                 & ZTCONSTAVG, ZTMINCL, ZRCONSTAVG, ZRMINCL, &
                 & ZTSNWS, ZRSNWS, YDCST, YDSOIL)


CALL SURFWS_MASSADJ(KIDIA, KFDIA, KLON, KLEVSN,         &
                 &  KLEVSNA,ZTHRESWS,                   &
                 &  ZDSN,ZDSNREAL,ZSNDEPTH,ZSNDEPTHREAL,&
                 &  ZRSN, ZSSN, ZRSNMAX, ZDSNTOT,       &
                 &  ZRCONSTAVG, ZRMINCL,                &
                 &  ZSNPERT,                            &
                 &  ZTSNWS,                             &
                 &  ZRSNWS, ZSSNWS,ZWSNWS, ZRSNTOP,     &
                 &  YDCST, YDSOIL)

IF (NSNMLWS == 3_JPIM) THEN
CALL SURFWS_TSNADJ(KIDIA, KFDIA,KLON,KLEVSN,                  &
                 & KLEVSNA, KLEVMID, ZTHRESWS,                &
                 & ZSNDEPTH, ZTCONSTAVG, ZTCONSTSTD, ZTMINCL, &
                 & ZRSN, ZSSN, ZTSN, ZHSN,                    &
                 & PTSKIN, ZDSNTOT, ZACTDEPTH, ZSADEPTH,      &
                 & ZTSNBOTTOM, ZTSNTOP, ZTSNMIDDLE,           &
                 & ZSNPERT,                                   &
                 & ZTSNWS, ZSSNWS,ZWSNWS,                     &
                 & YDCST, YDSOIL )
ENDIF


!************************************
! 2.4 Final update of snow variables:
DO JL=KIDIA,KFDIA
  IF (.NOT. LLNOSNOW(JL) .AND. KLEVSNA(JL)>1 ) THEN
    ! Safety check to avoid crashes in surfws
    IF ( (.NOT. ANY(ZTSNWS(JL,1:KLEVSN)<100._JPRB)) .AND. (.NOT. ANY(ZTSNWS(JL,1:KLEVSN)>RTT)) ) THEN
      PTSN(JL,1:KLEVSN)=ZTSNWS(JL,1:KLEVSN) 
      PRSN(JL,1:KLEVSN)=ZRSNWS(JL,1:KLEVSN)
      PSSN(JL,1:KLEVSN)=ZSSNWS(JL,1:KLEVSN)
      PWSN(JL,1:KLEVSN)=ZWSNWS(JL,1:KLEVSN)
    ELSE
      WRITE(*,*) 'Very cold snow temperature in surfws_ctl. Reset to initial (not-warm-started) values.'
      WRITE(*,*) 'Tsn warm start:',PTSN(JL,1:KLEVSN),ZTSNWS(JL,1:KLEVSN)
    ENDIF
  ENDIF
ENDDO

DEALLOCATE( ZTCONSTAVG )
DEALLOCATE( ZTCONSTSTD )
DEALLOCATE( ZRCONSTAVG )
DEALLOCATE( ZRCONSTSTD )


END ASSOCIATE

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_CTL_MOD:SURFWS_CTL',1,ZHOOK_HANDLE)

END SUBROUTINE SURFWS_CTL
END MODULE SURFWS_CTL_MOD
