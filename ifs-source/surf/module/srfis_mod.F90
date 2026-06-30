MODULE SRFIS_MOD
CONTAINS
SUBROUTINE SRFIS(KIDIA , KFDIA  , KLON  , KLEVS ,&
 & PTMST  , PTIAM1M   , PAHFSTI, PEVAPTI,&
 & PSLRFL , PSSRFLTI  , PTIA   , LDICE , LDNH   ,&
 & YDCST  , YDSOIL)

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST   , ONLY : TCST
USE YOS_SOIL  , ONLY : TSOIL

USE SRFWDIFS_MOD

#ifdef DOC
! (C) Copyright 2011- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFIS* - Computes temperature changes in soil.

!     PURPOSE.
!     --------
!**   Computes temperature evolution of sea ice
!**   INTERFACE.
!     ----------
!          *SRFIS* IS CALLED FROM *SURFS*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*      NUMBER OF SOIL LAYERS
!    *KTILES*     NUMBER OF SURFACE TILES

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDICE*      ICE MASK (TRUE for sea ice)
!    *LDNH*       TRUE FOR NORTHERN HEMISPHERE

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PTIAM1M*    SEA ICE TEMPERATURE                            K
!    *PSLRFL*     NET LONGWAVE  RADIATION AT THE SURFACE         W/M**2
!    *PAHFSTI*    TILE SURFACE SENSIBLE HEAT FLUX                W/M2
!    *PEVAPTI*    TILE SURFACE MOISTURE FLUX                     KG/M2/S
!    *PSSRFLTI*   TILE NET SHORTWAVE RADIATION FLUX AT SURFACE   W/M2

!     UPDATED PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PTIA*       SOIL TEMPERATURE                               K

!     METHOD.
!     -------
!          Parameters are set and the tridiagonal solver is called.

!     EXTERNALS.
!     ----------
!     *SRFWDIFS*

!     REFERENCE.
!     ----------
!          See documentation.

!     Original
!     --------
!            Simplified version based on SRFI
!       M. Janiskova              E.C.M.W.F.     25-07-2011  

!     Modifications
!     -------------

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS

REAL(KIND=JPRB),    INTENT(IN)   :: PTMST
REAL(KIND=JPRB),    INTENT(IN)   :: PTIAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PAHFSTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSLRFL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSSRFLTI(:,:)

LOGICAL,            INTENT(IN)   :: LDICE(:)
LOGICAL,            INTENT(IN)   :: LDNH(:)

TYPE(TCST),         INTENT(IN)   :: YDCST
TYPE(TSOIL),        INTENT(IN)   :: YDSOIL

REAL(KIND=JPRB),    INTENT(OUT)  :: PTIA(:,:)


!      LOCAL VARIABLES

REAL(KIND=JPRB) :: ZSURFL(KLON)
REAL(KIND=JPRB) :: ZRHS(KLON,KLEVS), ZCDZ(KLON,KLEVS),&
 & ZLST(KLON,KLEVS),&
 & ZTIA(KLON,KLEVS)
REAL(KIND=JPRB) :: ZDAI(KLON,KLEVS)
LOGICAL ::LLALLAYS, LLDOICE

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZCONS1, ZCONS2, ZTHFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SRFIS_MOD:SRFIS',0,ZHOOK_HANDLE)
ASSOCIATE(RLSTT=>YDCST%RLSTT, &
 & RCONDSICE=>YDSOIL%RCONDSICE, RDAI=>YDSOIL%RDAI, RDANSICE=>YDSOIL%RDANSICE, &
 & RDARSICE=>YDSOIL%RDARSICE, RRCSICE=>YDSOIL%RRCSICE, RSIMP=>YDSOIL%RSIMP, &
 & RTFREEZSICE=>YDSOIL%RTFREEZSICE, RTMELTSICE=>YDSOIL%RTMELTSICE)

!*    0. INITIALIZATION
!     ------------------

LLDOICE = .FALSE.
DO JL=KIDIA,KFDIA
  IF (LDICE(JL)) THEN
    LLDOICE = .TRUE.  ! if any point is sea ice
  ENDIF

  ZSURFL(JL) = 0.0_JPRB
ENDDO

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    ZLST(JL,JK) = 0.0_JPRB
    ZCDZ(JL,JK) = 0.0_JPRB
    ZRHS(JL,JK) = 0.0_JPRB

    ZTIA(JL,JK) = 0.0_JPRB
    PTIA(JL,JK) = RTFREEZSICE
  ENDDO
ENDDO

!* Computation done for only top or all soil layers

LLALLAYS = .TRUE.    ! done for all layers
!LLALLAYS = .FALSE.   ! done for top layer only

IF (LLDOICE) THEN

!*         1. SET UP SOME CONSTANTS.
!             --- -- ---- ----------

!*    PHYSICAL CONSTANTS.
!     -------- ----------

  DO JK=1,KLEVS-1
    DO JL=KIDIA,KFDIA
      ZDAI(JL,JK)=RDAI(JK)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    IF (LDNH(JL)) THEN
      ZDAI(JL,KLEVS)=RDARSICE-(RDAI(1)+RDAI(2)+RDAI(3))
    ELSE
      ZDAI(JL,KLEVS)=RDANSICE-(RDAI(1)+RDAI(2)+RDAI(3))
    ENDIF
  ENDDO

!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------

  ZCONS1=PTMST*RSIMP*2.0_JPRB
  ZCONS2=1.0_JPRB-1.0_JPRB/RSIMP


!*         2. Compute net heat flux at the surface.
!             -------------------------------------

  DO JL=KIDIA,KFDIA
    IF (LDICE(JL)) THEN
      ZTHFL=PAHFSTI(JL,2)+RLSTT*PEVAPTI(JL,2)
      ZSURFL(JL)=PSSRFLTI(JL,2)+PSLRFL(JL)+ZTHFL
    ENDIF
  ENDDO

!     Layer 1

  JK=1
  DO JL=KIDIA,KFDIA
    IF (LDICE(JL)) THEN
      ZLST(JL,JK)=ZCONS1*RCONDSICE/(ZDAI(JL,JK)+ZDAI(JL,JK+1))
      ZCDZ(JL,JK)=RRCSICE*ZDAI(JL,JK)
      ZRHS(JL,JK)=PTMST*ZSURFL(JL)/ZCDZ(JL,JK)
    ENDIF
  ENDDO

  IF (LLALLAYS) THEN

!     Layers 2 to KLEVS-1
    DO JK=2,KLEVS-1
      DO JL=KIDIA,KFDIA
        IF (LDICE(JL)) THEN
          ZLST(JL,JK)=ZCONS1*RCONDSICE/(ZDAI(JL,JK)+ZDAI(JL,JK+1))
          ZCDZ(JL,JK)=RRCSICE*ZDAI(JL,JK)
        ENDIF
      ENDDO
    ENDDO

!     Layers KLEVS
    JK=KLEVS
    DO JL=KIDIA,KFDIA
      IF (LDICE(JL)) THEN
        ZLST(JL,JK)=ZCONS1*RCONDSICE/(2.*ZDAI(JL,JK))
        ZCDZ(JL,JK)=RRCSICE*ZDAI(JL,JK)
        ZRHS(JL,JK)=(RTFREEZSICE/RSIMP)*ZLST(JL,JK)/ZCDZ(JL,JK)
      ENDIF
     ENDDO
  ENDIF

!*         4. Call tridiagonal solver
!             -----------------------


  CALL SRFWDIFS(KIDIA,KFDIA,KLON,KLEVS,PTIAM1M,ZLST,ZRHS,ZCDZ,ZTIA,&
   & LDICE,LLALLAYS,YDSOIL)

!*         6. New temperatures
!             ----------------

  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      IF (LDICE(JL)) THEN
        PTIA(JL,JK)=PTIAM1M(JL,JK)*ZCONS2+ZTIA(JL,JK)
      ELSE
        PTIA(JL,JK)=RTFREEZSICE
      ENDIF
      PTIA(JL,JK)=MIN(PTIA(JL,JK),RTMELTSICE)
    ENDDO
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFIS_MOD:SRFIS',1,ZHOOK_HANDLE)

END SUBROUTINE SRFIS
END MODULE SRFIS_MOD



