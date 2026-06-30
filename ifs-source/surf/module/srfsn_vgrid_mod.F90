MODULE SRFSN_VGRID_MOD
CONTAINS
SUBROUTINE SRFSN_VGRID(KIDIA, KFDIA, KLON, KLEVSN,LLNOSNOW,PSDOR,&
                       PSSN,PRSN,&
                       RLEVSNMIN,RLEVSNMAX, &
                       RLEVSNMIN_GL,RLEVSNMAX_GL, &
                       PDSNOUT,ILEVSNA,LEMIN)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE ABORT_SURF_MOD

! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFSN_VGRID* - Define snow verical grid 
!     PURPOSE.
!     --------
!          THIS ROUTINE DEFINES THE SNOW VERTICAL GRID

!**   INTERFACE.
!     ----------
!          *SRFSN_VGRID* IS CALLED FROM *SRFSN_DRIVER*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*      NUMBER OF HORIZONTAL POINTS
!    *KLEVSN*      NUMBER OF VERTICAL SNOW LAYERS


!     INPUT PARAMETERS (REAL):
!    *PSSN*         SNOW MULTI-LAYER MASS AT T-1, not regridded
!    *PRSN*         SNOW MULTI-LAYER DENSITY AT T-1, not regridded
!    *RLEVSNMIN*    MINIMUM SNOW DEPTH FOR EACH SNOW LAYER
!    *RLEVSNMAX*    MAXIMUM SNOW DEPTH FOR EACH SNOW LAYER
!    *RLEVSNMIN_GL* MINIMUM SNOW DEPTH FOR EACH SNOW LAYER OVER GLACIERS
!    *RLEVSNMAX_GL* MAXIMUM SNOW DEPTH FOR EACH SNOW LAYER OVER GLACIERS

!     INPUT PARAMETERS (LOGICAL):
!    *LLNOSNOW*   NOSNOW MASK (TRUE THEN POINT IS SNOW-FREE)
!    *PSDOR*      OROGRAPHIC PARAMETER                           m

!     OUTPUT PARAMETERS 
!    *PDSNOUT*    NEW VERTICAL THICKNESS OF SNOW LAYERS
!    *ILESVNA*    (optional) NUMBER OF ACTIVE VERTICAL LAYERS
!
!    INPUT PARAMETER, OPTIONAL
!    *LEMIN*    (optional) TRUE IF IT IS INNER LOOP MINIMISATION

!     METHOD.
!     -------
!          

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          

!     Modifications:
!     Original   E. Dutra      ECMWF     04/12/2015
!     G. Arduini: using spatial variable snow thickness depending on
!                 orography and snow depth

!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments 
INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVSN
LOGICAL           , INTENT(IN)   :: LLNOSNOW(:) 
REAL(KIND=JPRB), INTENT(IN) :: PSDOR(:)

REAL(KIND=JPRB), INTENT(IN) :: PSSN(:,:)
REAL(KIND=JPRB), INTENT(IN) :: PRSN(:,:) 

REAL(KIND=JPRB), INTENT(IN) :: RLEVSNMIN(:) 
REAL(KIND=JPRB), INTENT(IN) :: RLEVSNMAX(:) 
REAL(KIND=JPRB), INTENT(IN) :: RLEVSNMIN_GL(:) 
REAL(KIND=JPRB), INTENT(IN) :: RLEVSNMAX_GL(:) 
LOGICAL           , OPTIONAL, INTENT(IN)   :: LEMIN

REAL(KIND=JPRB),    INTENT(OUT) :: PDSNOUT(:,:)
INTEGER(KIND=JPIM), INTENT(OUT):: ILEVSNA(:)

! Local variables 
INTEGER(KIND=JPIM) :: JL,JK
INTEGER(KIND=JPIM) :: KLMAX

REAL(KIND=JPRB)   :: ZLEVMAX(KLEVSN),ZLEVMINA(KLEVSN)
REAL(KIND=JPRB)   :: ZDSNTOT, ZDSN(KLEVSN)
REAL(KIND=JPRB)   :: ZTMP0,ZTMP1,ZEPSILON

REAL(KIND=JPRB)   :: ZRLEVSNMAX(KLEVSN), ZRLEVSNMIN(KLEVSN)
REAL(KIND=JPRB)   :: ZALFA_G
REAL(KIND=JPRB)   :: ZDSNDIFF, ZORODIFF
REAL(KIND=JPRB)   :: ZVAR_LEVSNMIN, ZVAR_LEVSNMAX
REAL(KIND=JPRB)   :: ZDSNTHR, ZOROTHR
REAL(KIND=JPRB)   :: ZOROTHR_MIN,ZOROTHR_MAX
REAL(KIND=JPRB)   :: ZGLACTHR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_VGRID_MOD:SRFSN_VGRID',0,ZHOOK_HANDLE)

!    -----------------------------------------------------------------

! 0. Define constants from grid definition
ZEPSILON=10._JPRB*EPSILON(ZEPSILON)

! Thresholds for snow layer thickness computation:
ZDSNTHR=0.25_JPRB    ! minimum snow depth to use variable thickness
ZOROTHR_MIN=50._JPRB ! std dev oro minimum to use variable thickness

ZALFA_G=0.10_JPRB    ! 

PDSNOUT(KIDIA:KFDIA,1:KLEVSN) = 0._JPRB
!---
IF (PRESENT(LEMIN)) THEN
  ZGLACTHR=2000._JPRB
  ZOROTHR_MIN=80._JPRB ! std dev oro min thr
ELSE
  ZGLACTHR=5000._JPRB
  ZOROTHR_MIN=50._JPRB ! std dev oro min thr
ENDIF

DO JL=KIDIA,KFDIA
   ! layers snow depth 
  ZDSN(1:KLEVSN) = PSSN(JL,1:KLEVSN) / PRSN(JL,1:KLEVSN) 
   ! total snow depth 
  ZDSNTOT = SUM(ZDSN(1:KLEVSN),DIM=1)
  IF (SUM(PSSN(JL,1:KLEVSN),DIM=1)>=ZGLACTHR) THEN ! Glacier check
    ZRLEVSNMIN(1:KLEVSN) = RLEVSNMIN_GL(1:KLEVSN)
    ZRLEVSNMAX(1:KLEVSN) = RLEVSNMAX_GL(1:KLEVSN)
  ELSEIF  (ZDSNTOT>ZDSNTHR ) THEN
    ZDSNDIFF=ZDSNTOT-ZDSNTHR
   ! std dev oro
    IF (PSDOR(JL) < ZOROTHR_MIN) THEN 
    ! we do not change discretization over flat terrain
      ZRLEVSNMIN(1:KLEVSN) = RLEVSNMIN(1:KLEVSN) 
      ZRLEVSNMAX(1:KLEVSN) = RLEVSNMAX(1:KLEVSN) 
    ELSE  
      ZVAR_LEVSNMIN        = 0.10_JPRB+ZALFA_G*ZDSNDIFF
      ZRLEVSNMIN(1:KLEVSN) = MIN(0.25_JPRB,ZVAR_LEVSNMIN)

      ZVAR_LEVSNMAX        = 0.15_JPRB+ZALFA_G*ZDSNDIFF
      ZRLEVSNMAX(1:KLEVSN) = ZVAR_LEVSNMAX
      ZRLEVSNMAX(1)        = ZRLEVSNMIN(1)
      ZRLEVSNMAX(1:KLEVSN) = MIN(0.30_JPRB, ZRLEVSNMAX(1:KLEVSN))
      ZRLEVSNMAX(KLEVSN)   = -1._JPRB ! accumulation layer
    ENDIF
  ELSE 
    ZRLEVSNMAX(1:KLEVSN) = RLEVSNMAX(1:KLEVSN) 
    ZRLEVSNMIN(1:KLEVSN) = RLEVSNMIN(1:KLEVSN) 
  ENDIF

  ! snow layer to absorbe residuals 
  KLMAX=MINLOC(ZRLEVSNMAX(1:KLEVSN),1)  
  ZLEVMAX(:) = ZRLEVSNMAX(:)     
  IF ( ZRLEVSNMAX(KLMAX) > 0._JPRB ) THEN
    CALL ABORT_SURF("SRFSN_VGRID: RLEVSNMAX(KLMAX) should be < 0")
  ENDIF
  ZLEVMAX(KLMAX) = 1.E10_JPRB 
  
! accumulated min snow depth 
  ZLEVMINA(1) = ZRLEVSNMIN(1)
  DO JK=2,KLEVSN
    ZLEVMINA(JK) = ZLEVMINA(JK-1) +  ZRLEVSNMIN(JK)
  ENDDO
  



! 1. Find # of active snow layer
  IF (LLNOSNOW(JL)) THEN
    ILEVSNA(JL) = 1
  ELSE
    ILEVSNA(JL) = 1
    DO JK=2,KLEVSN
      IF ( ZDSNTOT < ZLEVMINA(JK) ) THEN
        ILEVSNA(JL) = JK-1
        EXIT
      ELSE
        ILEVSNA(JL) = JK 
      ENDIF
    ENDDO
  ENDIF

! 2. Apply the layering
  IF (LLNOSNOW(JL)) THEN
    PDSNOUT(JL,1) = ZDSNTOT
  ELSE
    IF (ILEVSNA(JL) > 1 ) THEN
      PDSNOUT(JL,1) = MIN(ZDSNTOT,ZLEVMAX(1))
      ZTMP0= (ZDSNTOT-PDSNOUT(JL,1)) / (ILEVSNA(JL) - 1.0)
      ! 1st screen
      ZTMP1=PDSNOUT(JL,1)
      DO JK=2,ILEVSNA(JL)
        PDSNOUT(JL,JK) = MIN(ZTMP0,ZLEVMAX(JK))
        ZTMP1=ZTMP1+ PDSNOUT(JL,JK)
      ENDDO
      ! Fix total
      IF (ZDSNTOT-ZTMP1 > ZEPSILON ) THEN
        IF (ILEVSNA(JL)<KLMAX) THEN
           PDSNOUT(JL,ILEVSNA(JL)) = PDSNOUT(JL,ILEVSNA(JL)) + ZDSNTOT-ZTMP1
        ELSE
           PDSNOUT(JL,KLMAX) = PDSNOUT(JL,KLMAX) + ZDSNTOT-ZTMP1
        ENDIF
      ENDIF
    ELSE
      PDSNOUT(JL,1) = ZDSNTOT
    ENDIF 
  ENDIF

ENDDO ! END LOOP IN JL
                            
!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_VGRID_MOD:SRFSN_VGRID',1,ZHOOK_HANDLE)

END SUBROUTINE SRFSN_VGRID
END MODULE SRFSN_VGRID_MOD
