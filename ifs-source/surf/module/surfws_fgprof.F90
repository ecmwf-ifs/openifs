MODULE SURFWS_FGPROF_MOD
CONTAINS

SUBROUTINE SURFWS_FGPROF(KIDIA, KFDIA, KLON, KLEVSN,               &
                       & KLEVSNA, KLEVMID,                         & 
                       & ZTHRESWS, ZSNPERT,LLNOSNOW,               & ! THESE CAN BE MOVED TO SUSSOIL
                       & PTSN, PSSN, PRSN, PTSA,                   &
                       & PDSNTOT,PACTDEPTH, PSNDEPTH,PSNDEPTHR, PSADEPTH,    &
                       & PTSNBOTTOM, PTSNTOP, PTSNMIDDLE,          & 
                       & PRSNTOP, PRSNMAX,    &
                       & PTCONSTAVG, PTMINCL, PRCONSTAVG, PRMINCL, &
                       & PTSNWS, PRSNWS,                           & ! Output
                       & YDCST,YDSOIL )

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_SOIL , ONLY : TSOIL

USE ABORT_SURF_MOD

! (C) Copyright 2017- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURFWS_FGPROF* - Snow warm start multi-layer 
!     PURPOSE.
!     --------
!          THIS ROUTINE GENERATES PARAMETRIZED PROFILES OF 
!          SNOW TEMPERATURE AND DENSITY BASED ON PARAMETERS
!          SETUP IN SURFWS_INIT* 

!**   INTERFACE.
!     ----------
!          *SURFWS_FGPROF* IS CALLED FROM *SURFWS_CTL*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       Length of arrays
!    *KLEVSN*     Snow vertical levels
!    *KLEVSNA*    Snow vertical Active levels
!    *KLEVMID*    Snow middle levels (surfws_init)
!    *PTMINCL*    Cluster index for temperature exp function
!    *PRMINCL*    Cluster index for density exp function


!     INPUT PARAMETERS (REAL):
!    *ZTHRESWS*      snow depth threhsold for warm start
!    *ZSNPERT*       snow depth threshold for glaciers
!    *

!     INPUT PARAMETERS (LOGICAL):
!    *LLNOSNOW*     NOSNOW MASK (TRUE IF NO SNOW)

!     INPUT PARAMETERS  (REAL):
!    *PTSN*       SNOW TEMPERATURE single layer               (K)
!    *PSSN*       SNOW MASS        single layer               (kg m-2)
!    *PRSN*       SNOW DENSITY     single layer               (kg m-3)
!    *PTSA*       soil temperature top layer                  (K)
!    *PDSNTOT*    Total snow depth                            (m)
!    *PACTDEPTH*  Active snow depth                           (m)
!    *PSNDEPTH*   Snow depth per layer wrt 0                  (m)
!    *PSNDEPTHR*  Snow depth per layer (full) wrt 0           (m)
!    *PSADEPTH*   Soil depth top layer                        (m)
!    *PTSNBOTTOM* Snow temperature bottom layer               (K)
!    *PTSNTOP*    Snow temperature top layer                  (K)
!    *PTSNMIDDLE* Snow temperature KLEVMID layer              (K)
!    *PRSNTOP*    Snow density top layer                      (kg m-3)
!    *PRSNMAX*    Snow density MAX allowed                    (kg m-3)
!    *PTCONSTAVG* Constants for temperature exp function
!    *PRCONSTAVG* Constants for density exp function
!    
!     OUTPUT PARAMETERS (REAL):
!    *PTSNWS*        Snow tempeature warm started
!    *PRSNWS*        Snow density    warm started


!     METHOD.
!     -------
!     Compute snow temperature and density profiles using 
!     exponential functions. Parameters of the functions are 
!     pre-computed using a k-clustering algorithm, see Arduini et al. (2019)

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          Arduini et al. (2019)

!     Modifications:
!     Original   G. Arduini      ECMWF     28/07/2017

!     ------------------------------------------------------------------


IMPLICIT NONE

! Input variables:
INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA, KFDIA
INTEGER(KIND=JPIM), INTENT(IN) :: KLON
INTEGER(KIND=JPIM), INTENT(IN) :: KLEVSN
INTEGER(KIND=JPIM),  INTENT(IN):: KLEVSNA(:),KLEVMID(:)
LOGICAL,            INTENT(IN) :: LLNOSNOW(:)

!
REAL(KIND=JPRB),   INTENT(IN)   :: PTSN(:), PRSN(:), PSSN(:), PTSA(:)
REAL(KIND=JPRB),   INTENT(IN)   :: PDSNTOT(:), PACTDEPTH(:), PSADEPTH(:)
REAL(KIND=JPRB),   INTENT(IN)   :: PSNDEPTH(:,:), PSNDEPTHR(:,:)
REAL(KIND=JPRB),   INTENT(IN)   :: PTSNBOTTOM(:), PTSNTOP(:),PTSNMIDDLE(:)
REAL(KIND=JPRB),   INTENT(IN)   :: PRSNTOP(:), PRSNMAX(:)
REAL(KIND=JPRB),   INTENT(IN)   :: PTCONSTAVG(:,:), PRCONSTAVG(:,:)
INTEGER(KIND=JPIM),INTENT(IN)   :: PTMINCL(:), PRMINCL(:)


!REAL(KIND=JPRB),  INTENT(IN)   :: PRMLRA, PRMLRB, PRMLRC, PRMLRD, PRMLRE

REAL(KIND=JPRB),    INTENT(IN) :: ZTHRESWS,ZSNPERT

TYPE(TCST)     , INTENT(IN) :: YDCST
TYPE(TSOIL)    , INTENT(IN) :: YDSOIL

! Output
REAL(KIND=JPRB),    INTENT(INOUT) :: PTSNWS(:,:), PRSNWS(:,:)

! Local variables:
REAL(KIND=JPRB)    :: ZAT, ZBT, ZCT, ZAR, ZBR, ZCR
REAL(KIND=JPRB)    :: ZP1, ZP2, ZP3, ZP4, ZP5 

REAL(KIND=JPRB)    :: ZTMULT
REAL(KIND=JPRB)    :: ZEPSILON

INTEGER(KIND=JPIM) :: JL,JK, ILEVMID


REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE


!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_FGPROF_MOD:SURFWS_FGPROF',0,ZHOOK_HANDLE)

!    -----------------------------------------------------------------

ASSOCIATE(RTT=>YDCST%RTT,RHOMINSND=>YDSOIL%RHOMINSND, NSNMLWS=>YDSOIL%NSNMLWS)
!        & RLMLT=>YDCST%RLMLT, RPI=>YDCST%RPI,    &
!        & RLWCSWEA=>YDSOIL%RLWCSWEA, RLWCSWEB=>YDSOIL%RLWCSWEB, &
!        & RLWCSWEC=>YDSOIL%RLWCSWEC, RTEMPAMP=>YDSOIL%RTEMPAMP, &
!        & RDSNMAX=>YDSOIL%RDSNMAX, RHOMINSND=>YDSOIL%RHOMINSND, &
!        & RHOMAXSN_NEW=>YDSOIL%RHOMAXSN_NEW, RDAT=>YDSOIL%RDAT  )

ZAT = 0._JPRB
ZBT = 0._JPRB
ZCT = 0._JPRB
ZAR = 0._JPRB
ZBR = 0._JPRB
ZCR = 0._JPRB

ZEPSILON = 10E4*EPSILON(ZEPSILON)

DO JL=KIDIA, KFDIA

  IF (.NOT. LLNOSNOW(JL)) THEN 
  IF ( (KLEVSNA(JL) > 1) .AND. (PDSNTOT(JL) >= ZTHRESWS) ) THEN
! Definition of constants:  
    ILEVMID = KLEVMID(JL)
! -- Temperature
    ZCT  = PTCONSTAVG(JL, PTMINCL(JL) )
    IF ( PDSNTOT(JL) < 0.20_JPRB ) THEN  ! IF small snow depth
      ZTMULT=1._JPRB

      ZBT=(PTSNBOTTOM(JL)-PTSNTOP(JL))/SIGN(MAX(ZEPSILON,                   &
                                   &ABS(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCT ))-1) ),&
                                   &(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCT ))-1) )
    !!ZAT=PTSNTOP(JL) - ZBT*MAX(ZEPSILON,EXP( -PSNDEPTH(JL, 1)*ZCT ))
      ZAT=PTSNTOP(JL) - ZBT!!!*MAX(ZEPSILON,EXP( -PSNDEPTH(JL, 1)*ZCT ))

    ELSE  !  Deep snowpack
      ZTMULT=0._JPRB
     !!! CHANGE HERE PTSNOWTOP WITH PTSKIN!! 
     !IF (PTSNTOP(JL)>=PTSN(JL) ) THEN  ! Tskin greater than avg snow temp (typically daytime)
      IF (PTSNTOP(JL)>=PTSNMIDDLE(JL) ) THEN  ! Tskin greater than avg snow temp (typically daytime)
        ZBT=(PTSNBOTTOM(JL)-PTSNTOP(JL))/SIGN(MAX(ZEPSILON,                                        &
                                &ABS((MAX(ZEPSILON,EXP(-PACTDEPTH(JL) * ZCT ))-1.0_JPRB))),&
                                &((MAX(ZEPSILON,EXP(-PACTDEPTH(JL) * ZCT ))-1.0_JPRB) )    )
        IF (PTSNTOP(JL)<=RTT) THEN
          ZAT=PTSNTOP(JL)-ZBT
        ELSE 
          ZAT=RTT-ZBT
        ENDIF
      ELSE ! Tskin smaller than avg snow temp (typically nighttime)
        ZBT=(PTSNBOTTOM(JL)-PTSNTOP(JL))/SIGN(MAX(ZEPSILON,                                 &
                                &ABS(MAX(ZEPSILON,EXP(-(PACTDEPTH(JL)+PSADEPTH(JL)) * ZCT)) - 1._JPRB)),&
                                &(MAX(ZEPSILON,EXP(-(PACTDEPTH(JL)+PSADEPTH(JL)) * ZCT)) - 1._JPRB)     )
      !!ZAT=PTSNTOP(JL) - ZBT*MAX(ZEPSILON,EXP( -PSNDEPTH(JL, 1)*ZCT ))
        ZAT=PTSNTOP(JL) - ZBT!!!!*MAX(ZEPSILON,EXP( -PSNDEPTH(JL, 1)*ZCT ))
      ENDIF ! end if in tskin - tsn sign
    ENDIF   ! end if in snow depth value (deep/thin snowpack)

! -- Density

    ZCR    = PRCONSTAVG(JL, PRMINCL(JL) )

    ZBR=(PRSN(JL)-PRSNTOP(JL))/SIGN(MAX(ZEPSILON,                                                  &
                               &ABS(EXP(-PSNDEPTHR(JL, ILEVMID)*ZCR)-EXP(-PSNDEPTHR(JL,1)*ZCR))),&
                               &(EXP(-PSNDEPTHR(JL, ILEVMID)*ZCR)-EXP(-PSNDEPTHR(JL,1)*ZCR))     )
    ZAR=PRSNTOP(JL) - ZBR*EXP(-PSNDEPTHR(JL,1)*ZCR)
      
! Loop over the vertical levels
    DO JK=1, KLEVSNA(JL)
! -- Temperature profile parametrisation
      PTSNWS(JL, JK) = MIN( RTT,ZAT + ZBT * MAX(ZEPSILON,EXP( -PSNDEPTH(JL, JK) * ZCT )) )
      ! Safety checks to be tested in TSNWS = 2 and 3:
      ! For option 1 better without!!!
      IF (NSNMLWS == 2_JPIM .OR. NSNMLWS == 3_JPIM) THEN
        IF (PTSNWS(JL, JK) > RTT ) PTSNWS(JL, JK) = RTT
        IF (KLEVSNA(JL) >= 2 ) THEN
          IF (PDSNTOT(JL) >= 0.5_JPRB) THEN
            PTSNWS(JL, KLEVSNA(JL)-1)=PTSNMIDDLE(JL)!PTSN(JL)
          ELSEIF (PDSNTOT(JL) >= 0.20_JPRB ) THEN
            PTSNWS(JL, ILEVMID)=PTSNMIDDLE(JL) !PTSN(JL)
          ENDIF
        ENDIF
        IF (PDSNTOT(JL) < 0.20_JPRB ) THEN
          PTSNWS(JL, ILEVMID)=PTSNMIDDLE(JL) !PTSN(JL)
        ENDIF
      ENDIF
      IF (PTSNWS(JL, JK) <100._JPRB) PTSNWS(JL, JK) = PTSNMIDDLE(JL)

! -- RHO profile parametrisation:
      PRSNWS(JL, JK) = ZAR + ZBR *  MAX(ZEPSILON,EXP( - PSNDEPTHR(JL, JK) * ZCR )) 
      PRSNWS(JL, ILEVMID)=PRSN(JL)
      IF (PRSNWS(JL, JK) < RHOMINSND) THEN 
        PRSNWS(JL, JK) = RHOMINSND
      ENDIF
      IF (PRSNWS(JL, JK) > PRSNMAX(JL)) THEN 
        PRSNWS(JL, JK)  = PRSNMAX(JL)
      ENDIF
      IF (PSSN(JL) >= ZSNPERT) THEN
        PRSNWS(JL, KLEVSNA(JL)-1:KLEVSNA(JL))=300._JPRB
      ENDIF
    ENDDO
    IF ( PSSN(JL) >= ZSNPERT ) PTSNWS(JL, KLEVSNA(JL))=MIN(RTT, PTSA(JL))
    IF ( PSSN(JL) >= ZSNPERT ) PRSNWS(JL, 1:KLEVSN)=300._JPRB

    IF (KLEVSNA(JL) < KLEVSN) THEN
      PTSNWS(JL, KLEVSNA(JL)+1:KLEVSN)=PTSNWS(JL, KLEVSNA(JL))
      PRSNWS(JL, KLEVSNA(JL)+1:KLEVSN)=100._JPRB
    ENDIF

  ENDIF ! END IF IN WS THR
  ENDIF
END DO ! END DO IN JL

END ASSOCIATE
!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_FGPROF_MOD:SURFWS_FGPROF',1,ZHOOK_HANDLE)

END SUBROUTINE SURFWS_FGPROF
END MODULE SURFWS_FGPROF_MOD
