MODULE SRFSN_SSRABS_MOD
CONTAINS
SUBROUTINE SRFSN_SSRABS(KIDIA,KFDIA,KLON,KLEVSN,&
 & LLNOSNOW,PFRTI,PSSRFLTI,&
 & PSSNM1M,PRSNM1M,&
 & YDSOIL,YDCST,&
 & PSNOTRS)


USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOS_SOIL , ONLY : TSOIL 
USE YOS_CST  , ONLY : TCST

USE ABORT_SURF_MOD

! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFSN_SSRABS* - Shortwave radiation absorption by snow
!     PURPOSE.
!     --------
!          THIS ROUTINE COMPUTES SW ABSORBED BY EACH SNOWPACK LAYER 

!**   INTERFACE.
!     ----------
!          *SRFSN_SSRABS* IS CALLED FROM *SRFSN_DRIVER*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF POINTS LON
!    *KLEVSN*     NUMBER OF MAX VERTICAL SNOW LAYERS


!     INPUT PARAMETERS (REAL):
!    *PFRTI*      TILE FRACTION                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LLNOSNOW*     SNOW/NO-SNOW MASK (TRUE IF NO-SNOW)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!     
!    *PSSNM1M*    SNOW WATER EQUIVALENT T-1                     kg/m**2
!    *PRSNM1M*    SNOW DENSITY          T-1                     kg/m**3
!    *PSSRFLTI*   TILED SHORTWAVE RADIATION AT SURFACE           W/m**2

!     OUTPUT FLUX  (UNFILTERED,REAL):
!    *PSNOTRS*        SOLAR RADIATION ABSORBED BY EACH SNOW LAYER W/m**2


!     METHOD.
!     -------
          

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          

!     Modifications:
!     Original   G. Arduini      ECMWF     04/12/2015

!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments 
INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVSN
LOGICAL           , INTENT(IN)   :: LLNOSNOW(:) 

REAL(KIND=JPRB)   , INTENT(IN)   :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSSRFLTI(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PSSNM1M(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PRSNM1M(:,:)

TYPE(TSOIL)       , INTENT(IN)   :: YDSOIL
TYPE(TCST)        , INTENT(IN)   :: YDCST

REAL(KIND=JPRB)   , INTENT(OUT)  :: PSNOTRS(:,:)

! Local variables 
REAL(KIND=JPRB)    :: ZDSN(KLEVSN)    ! actual snow depth
REAL(KIND=JPRB)    :: ZGSNS
REAL(KIND=JPRB)    :: ZSNQRAD
REAL(KIND=JPRB)    :: ZSNEXTCOEFF(KLEVSN)
REAL(KIND=JPRB)    :: ZSNOTRSTMP(KLEVSN+1)
REAL(KIND=JPRB)    :: ZSNSOABS(KLON)
REAL(KIND=JPRB)    :: ZFRSN(KLON)
REAL(KIND=JPRB)    :: ZEPSILON
INTEGER(KIND=JPIM) :: KLACT
INTEGER(KIND=JPIM) :: KSNTILES

INTEGER(KIND=JPIM) :: JL,JK,JT,JTILE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!! INCLUDE FUNCTIONS

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_SSRABS_MOD:SRFSN_SSRABS',0,ZHOOK_HANDLE)

!    -----------------------------------------------------------------
ASSOCIATE(RDSNMAX=>YDSOIL%RDSNMAX, SSAG1=>YDSOIL%SSAG1, SSAG2=>YDSOIL%SSAG2, &
         & SSAG3=>YDSOIL%SSAG3, SSAGSNSMAX=>YDSOIL%SSAGSNSMAX, SSASNEXTMIN=>YDSOIL%SSASNEXTMIN, &
         & SSASNEXTMAX=>YDSOIL%SSASNEXTMAX, SSASNEXTCNST=>YDSOIL%SSASNEXTCNST )

KSNTILES=2
ZEPSILON=10._JPRB*EPSILON(ZEPSILON)

KLACT=1
DO JT=1,KSNTILES
  SELECT CASE(JT)
  CASE(1)
    JTILE=5
    ZFRSN(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,JTILE)
  CASE(2)
    JTILE=7
    ZFRSN(KIDIA:KFDIA)=PFRTI(KIDIA:KFDIA,JTILE)
  END SELECT
  DO JL=KIDIA,KFDIA
    IF (LLNOSNOW(JL)) THEN 
      PSNOTRS(JL,:) = 0.0_JPRB
    ELSE
  
  !! Preparation
      DO JK=1,KLEVSN
        IF (PSSNM1M(JL,JK) > ZEPSILON ) KLACT=JK
      ENDDO
      ZSNOTRSTMP(1:KLEVSN+1)=0._JPRB
      ZSNEXTCOEFF(1:KLEVSN)=25._JPRB
      ZSNQRAD=0._JPRB
  
        DO JK=1, KLACT
          ZDSN(JK)=MAX(0._JPRB, MIN(RDSNMAX,PSSNM1M(JL,JK)/PRSNM1M(JL,JK)))
          ! grain size from Anderson 1976
          ZGSNS=MIN(SSAGSNSMAX,(SSAG1 + SSAG3*PRSNM1M(JL,JK)**(4._JPRB) ))
          ! snow extinction coeff from Jordan 1991
          ZSNEXTCOEFF(JK)=MAX(SSASNEXTMIN, MIN(SSASNEXTMAX, SSASNEXTCNST*PRSNM1M(JL,JK)/ZGSNS**0.5_JPRB ) )
        ENDDO
        ZSNSOABS(JL)=MAX(ZEPSILON, exp(-ZSNEXTCOEFF(1)*ZDSN(1)) )
        ZSNOTRSTMP(1)= PSSRFLTI(JL,JTILE)*ZSNSOABS(JL)
        ZSNQRAD=PSSRFLTI(JL,JTILE)-PSSRFLTI(JL,JTILE)*(1._JPRB-ZSNSOABS(JL))
        IF (KLACT > 1) THEN
          DO JK=2, KLACT
            ZSNSOABS(JL)=exp(-ZSNEXTCOEFF(JK)*ZDSN(JK))
            ZSNOTRSTMP(JK)=ZSNQRAD*(1._JPRB-ZSNSOABS(JL))
            ZSNQRAD=ZSNQRAD-ZSNOTRSTMP(JK)
          ENDDO
        ENDIF
        ZSNOTRSTMP(KLACT+1)=ZSNQRAD
  ! Weighted average with tile fraction:
        DO JK=1, KLACT+1
          PSNOTRS(JL,JK) = PSNOTRS(JL,JK) + ZFRSN(JL)*ZSNOTRSTMP(JK)  
        ENDDO

    ENDIF 
  ENDDO
END DO

END ASSOCIATE
!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_SSRABS_MOD:SRFSN_SSRABS',1,ZHOOK_HANDLE)

END SUBROUTINE SRFSN_SSRABS

END MODULE SRFSN_SSRABS_MOD
