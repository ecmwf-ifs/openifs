! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SEDIMNT &
 !---input
 & ( KIDIA , KFDIA   , KLON, KLEV, PSEDIMV, &
 &   PCI   , PLSM    , & 
 &   PTSPHY, PT      , PAP , PAPH, PVERVEL, &
 !---prognostic fields
 &   PAERO , PTAERI , &
 !---output
 &   PTAERO, PFLUXAER )  

!**** *AER_SEDIMNT* -  ROUTINE FOR PARAMETRIZATION OF AEROSOL SEDIMENTATION

!      Olivier Boucher & Jean-Jacques Morcrette 
!      following the ice sedimentation scheme of Adrian Tompkins

!**   INTERFACE.
!     ----------
!          *AER_SEDIMNT* IS CALLED FROM *CALLPAR*.

! INPUTS:
! -------
! PTSPHY                : TIMESTEP                  (s)
! PAP     (KLON,KLEV)   : LEVEL PRESSURE            (Pa)
! PAPH    (KLON,KLEV+1) : HALF-LEVEL PRESSURE       (Pa)
! PCI     (KLON)        : FRACTION OF SEA ICE
! PLSM    (KLON)        : LAND-SEA MASK   
! PT      (KLON,KLEV)   : LEVEL TEMPERATURE         (K)
! PTAERI   (KLON,KLEV)  : INPUT TENDENCY            (xx kg-1 s-1)
! PVERVEL (KLON,KLEV)   : VERTICAL VELOCITY         (Pa s-1)
! PAERO    (KLON,KLEV)  : CONCENTRATION OF TRACERS  (xx kg-1)

! OUTPUTS:
! --------
! PTAERO (KLON,KLEV)   : TENDENCY                  (xx kg-1 s-1)
! PFLUXAER(KLON)        : SURFACE FLUX              (xx m-2 s-1)

!     EXTERNALS.
!     ----------
!          NONE

!     MODIFICATIONS.
!     -------------

!     SWITCHES.
!     --------

!     MODEL PARAMETERS
!     ----------------

!-----------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    ,ONLY : RG, RD

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

!---input fields
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEDIMV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLON), PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAERI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV)

!---output fields
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUXAER(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAERO(KLON,KLEV) 

!---local fields
REAL(KIND=JPRB) :: ZSEDFLX(KLON), ZAERONWM1(KLON), ZAERI(KLON,KLEV), &
       &           ZSOLAERS, ZSOLAERB, ZGDP, ZDTGDP, ZKK, & 
       &           ZRHO, ZAERONW, ZTAERO(KLON,KLEV)

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SEDIMNT',0,ZHOOK_HANDLE)
!- PAERO     in unit of xx kg-1 (mixing ratio)
!- ZAERONW   in unit of xx kg-1
!- ZSOLAERS, ZSOLAERB in xx kg-1  as m s-2/DeltaP xx m-2 = m s-2/(kg m-1 s-2) xx m-2
!- ZSEDFLX   in unit of xx m-2
!- ZWAIR     in unit of m s-1  = -RD *PVERVEL*PT / (RG*PAP)       vertical speed
!- PVERVEL   in unit of Pa s-1                                    vertical velocity
!- ZRHO      in unit of kg m-3


!--constants
ZKK= -RD/RG    ! in J kg-1 K-1 m-1 s2


!--initialisations of variables carried out from one layer to the next layer
!--actually not needed if (JK>1) test is on
DO JL=KIDIA,KFDIA 
  ZSEDFLX(JL)=0.0_JPRB
  ZAERONWM1(JL)=0.0_JPRB
ENDDO

DO JK=1, KLEV

  DO JL=KIDIA, KFDIA
!--initialisations
    ZSOLAERS=0.0
    ZSOLAERB=0.0
    ZGDP=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
    ZDTGDP=PTSPHY*ZGDP

!- at this stage, with only source fluxes and (simple) dry deposition (by reduction of 
!  source fluxes) having been applied, the initial 3D-field is PAERO=PCEN

!   ZAERI(JL,JK) = PAERO(JL,JK)

!-- In the future (?) with more sophisticated treatment of dry deposition, it will be 
!  necessary to get 
    ZAERI(JL,JK) = PAERO(JL,JK) + PTSPHY * PTAERI(JL,JK) 
!  with PTAERI containing the tendencies from dry deposition

! source from above
    IF (JK>1) THEN 
      ZSEDFLX(JL)=ZSEDFLX(JL)*ZAERONWM1(JL)  
      ZSOLAERS=ZSOLAERS+ZSEDFLX(JL)*ZDTGDP
    ENDIF

! sink to next layer
    ZRHO=PAP(JL,JK)/(RD*PT(JL,JK))
! ZWAIR is W =-RD/RG * PT * PVERVEL / PAP in m s-1 from PVERVEL in Pa s-1
!    ZWAIR=ZKK*PT(JL,JK)*PVERVEL(JL,JK)/PAP(JL,JK)
!    ZSEDFLX(JL)=ZWAIR*ZRHO*PAERO(JL,JK)

    ZSEDFLX(JL)=PSEDIMV(JL,JK)*ZRHO
    ZSOLAERB=ZSOLAERB+ZDTGDP*ZSEDFLX(JL)

!---implicit solver
    ZAERONW=(ZAERI(JL,JK)+ZSOLAERS)/(1.0_JPRB+ZSOLAERB)

!---new time-step AER variable needed for next layer
    ZAERONWM1(JL)=ZAERONW

!---tendency in unit of xx kg-1 s-1
    ZTAERO(JL,JK)=(ZAERONW-ZAERI(JL,JK))/PTSPHY

!-NB: tendencies from sedimentation is added to original tendencies
    PTAERO(JL,JK)=ZTAERO(JL,JK)+PTAERI(JL,JK)
!
  ENDDO
ENDDO

!---sedimentation flux to the surface
!---ZAERONWM1 now contains the surface concentration at the new timestep
!---PFLUXAER in unit of xx m-2 s-1 
DO JL=KIDIA,KFDIA 
  ZRHO=PAP(JL,KLEV)/(RD*PT(JL,KLEV))

!  ZWAIR=ZKK*PT(JL,KLEV)*PVERVEL(JL,KLEV)/PAP(JL,KLEV)
!  PFLUXAER(JL)=ZRHO*ZAERONWM1(JL)*ZWAIR

  PFLUXAER(JL)=ZRHO*ZAERONWM1(JL)*PSEDIMV(JL,KLEV)
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SEDIMNT',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SEDIMNT

