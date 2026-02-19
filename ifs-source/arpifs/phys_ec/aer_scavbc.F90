! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SCAVBC &
 & ( YDEAERSNK,KIDIA  , KFDIA  , KLON , KLEV  ,PSCAVBCR, PSCAVBCS, KSTEP , PDP, &
 &   PFLRAIN, PFLSNOW, PAERO, PTAERI, PTSPHY, &
 &   PTAERO, PFAERO &
 & )

!*** * AER_SCAVBC* - BELOW-CLOUD SCAVENGING OF TRACERS

!**   INTERFACE.
!     ----------
!          *AER_SCAVBC* IS CALLED FROM *CALLPAR*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        FROM O. BOUCHER (LOA, 1999-09) bc_scav

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2004-05-07
!                   2014-10-20  Samuel Remy : add flux output for tracing

!-----------------------------------------------------------------------
! INPUTS:
! -------
! PFLRAIN : liquid precipitation rate (kg m-2 s-1)
! PFLSNOW : solid precipitation rate  (kg m-2 s-1)
! PAERO   : aerosol amount             kg kg-1  
! PTAERI  : initial tendency accumulated from previous processes (kg kg-1)

! OUTPUTS:
!---------
! PTAERO : instantaneous tendency       kg kg-1 s-1

!-----------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST    ,ONLY : RG

USE YOEAERSNK ,ONLY : TEAERSNK

USE YOMLUN   , ONLY : NULOUT
!USE YOEDBUG  , ONLY : KSTPDBG, NSTPDBG

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

TYPE(TEAERSNK)    ,INTENT(INOUT):: YDEAERSNK
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV, KSTEP
REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV)

REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB),INTENT(IN)    :: PSCAVBCR, PSCAVBCS
REAL(KIND=JPRB),INTENT(IN)    :: PFLRAIN(KLON,KLEV+1), PFLSNOW(KLON,KLEV+1)
REAL(KIND=JPRB),INTENT(IN)    :: PAERO(KLON,KLEV), PTAERI(KLON,KLEV)

REAL(KIND=JPRB),INTENT(OUT)   :: PTAERO(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT)   :: PFAERO(KLON)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JK, JL
LOGICAL :: LLPRINT
REAL(KIND=JPRB) :: ZAERI(KLON,KLEV)

REAL(KIND=JPRB) :: ZINCR, ZPR, ZPS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SCAVBC',0,ZHOOK_HANDLE)
ASSOCIATE(RHO_ICE=>YDEAERSNK%RHO_ICE, RHO_WAT=>YDEAERSNK%RHO_WAT, R_R=>YDEAERSNK%R_R, &
 & R_S=>YDEAERSNK%R_S)

LLPRINT=.FALSE.
!DO JL=1,NSTPDBG
!  IF (KSTEP == KSTPDBG(JL)) THEN
!    LLPRINT=.TRUE.
!  ENDIF
!ENDDO
IF (LLPRINT) THEN
  WRITE(NULOUT,'("AER_SCAVBC:",6E12.5)') PSCAVBCR,PSCAVBCS,R_R,RHO_WAT,R_S,RHO_ICE
ENDIF

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZAERI(JL,JK) = PAERO(JL,JK) + PTSPHY * PTAERI(JL,JK)

    ZPR = 0.5_JPRB*(PFLRAIN(JL,JK)+PFLRAIN(JL,JK+1))
    ZPS = 0.5_JPRB*(PFLSNOW(JL,JK)+PFLSNOW(JL,JK+1))
    ZINCR = 0.75_JPRB * (ZPR*PSCAVBCR/R_R/RHO_WAT + ZPS*PSCAVBCS/R_S/RHO_ICE)

    PTAERO(JL,JK) = PTAERI(JL,JK)-ZINCR * ZAERI(JL,JK)
  ENDDO
ENDDO
!---PFAERO in unit of xx m-2 s-1 
DO JL=KIDIA,KFDIA
   PFAERO(JL)=0.0_JPRB
ENDDO

DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      PFAERO(JL) = PFAERO(JL) + (PTAERI(JL,JK)-PTAERO(JL,JK))*(PDP(JL,JK))/RG
   ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SCAVBC',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SCAVBC

