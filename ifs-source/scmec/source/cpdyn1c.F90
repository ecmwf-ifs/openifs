! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CPDYN1C(YDDIMV, LDTWOTL,KFLEV  ,PDT    ,PRESF  ,PVVEL0 ,PRCORI &
         &,PUT0   ,PVT0   ,PTT0   ,PQT0   ,PAT0   ,PLT0   ,PIT0&
         &,PKAP   ,PUG0   ,PVG0   ,PETADOTDPDETA  ,PDELP9&
         &,PUADV  ,PVADV  ,PTADV  ,PQADV&
         &,PUTEND ,PVTEND ,PTTEND ,PQTEND ,PATEND ,PLTEND ,PITEND)

#ifdef DOC
!**** *CPDYN1C* - adiabatic tendencies, Eulerian vertical advection.

!     Purpose.
!     --------
!           Computes the adiabatic tendencies.
!           Uses Eulerian vertical advection.  

!**   Interface.
!     ----------
!        *CALL* *CPDYN1C

!        Explicit arguments :
!        --------------------
!                 KFLEV   - Number of levels   (input)
!                 PUT0    - U time step (t)    (input)
!                 PVT0    - V time step (t)    (input)
!                 PTT0    - T time step (t)    (input)
!                 PQT0    - Q time step (t)    (input)
!                 PAT0    - A time step (t)    (input)
!                 PLT0    - L time step (t)    (input)
!                 PIT0    - I time step (t)    (input)
!                 PRESF   - full level pressure (input)
!                 PVVEL0  - vertical velocity  (input)
!                 petadotdpdeta  - ver. vel. half lev. (input)
!                 PKAP    - K=Rm/Cpm           (input)
!                 PUG0    - U geostrophic      (input)
!                 PVG0    - V geostrophic      (input)
!                 PUADV   - U hor. advection   (input)
!                 PVADV   - V hor. advection   (input)
!                 PTADV   - T hor. advection   (input)
!                 PQADV   - Q hor. advection   (input)
!                 PRCORI  - Coriolis Parameter (input)
!                 PUTEND  - U adiab. tendency (output)
!                 PVTEND  - V adiab. tendency (output)
!                 PTTEND  - T adiab. tendency (output)
!                 PQTEND  - Q adiab. tendency (output)
!                 PATEND  - A adiab. tendency (output)
!                 PLTEND  - L adiab. tendency (output)
!                 PITEND  - I adiab. tendency (output)

!        Implicit arguments :   NONE.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the single column model

!     Author.
!     -------
!        Joao Teixeira   *ECMWF*

!     Modifications.
!     --------------
!        Original : 94-01-12
!        J.Teixeira, May-1995: introduction of cloud variables
!                              and two time-level scheme.
!        M. Ko"hler  Nov-2002  add upwind and Lax schemes
!        M. Ko"hler  6-6-2006  Single Column Model integration within IFS 

!     ------------------------------------------------------------------
#endif

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLOG1C , ONLY : LETADOT, LUPWIND, LWADVCLD

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
LOGICAL, INTENT(IN)        :: LDTWOTL
INTEGER(KIND=JPIM)         :: KFLEV

!     DUMMY REAL SCALARS
REAL(KIND=JPRB) :: PDT
REAL(KIND=JPRB) :: PRCORI

REAL(KIND=JPRB) :: PUT0  (YDDIMV%NFLEVG), PVT0  (YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTT0  (YDDIMV%NFLEVG), PQT0  (YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PAT0  (YDDIMV%NFLEVG), PLT0  (YDDIMV%NFLEVG), PIT0  (YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: PKAP  (YDDIMV%NFLEVG), PRESF (YDDIMV%NFLEVG), PVVEL0(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PETADOTDPDETA(0:YDDIMV%NFLEVG)       , PDELP9(YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: PUG0  (YDDIMV%NFLEVG), PVG0  (YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PUADV (YDDIMV%NFLEVG), PVADV (YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTADV (YDDIMV%NFLEVG), PQADV (YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: PUTEND(YDDIMV%NFLEVG), PVTEND(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PTTEND(YDDIMV%NFLEVG), PQTEND(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PATEND(YDDIMV%NFLEVG), PLTEND(YDDIMV%NFLEVG), PITEND(YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZQVADV(YDDIMV%NFLEVG), ZTVADV(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUVADV(YDDIMV%NFLEVG), ZVVADV(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZAVADV(YDDIMV%NFLEVG), ZLVADV(YDDIMV%NFLEVG), ZIVADV(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUT1  (YDDIMV%NFLEVG), ZVT1  (YDDIMV%NFLEVG), ZTT1  (YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVERVEL(YDDIMV%NFLEVG)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JLEV

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZALFA, ZBETA, ZCONST1, ZCONST2, ZDTCORI, ZA, ZB

LOGICAL :: LLAXMETHOD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPDYN1C',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)
!LLAXMETHOD = .TRUE. 
LLAXMETHOD = .FALSE.


!*       0.   COMPUTATION OF VERTICAL ADVECTION.
!             ----------------------------------

IF (LUPWIND) THEN

 DO JLEV=1,NFLEVG

! interpolate half level petadotdpdeta to full levels

  IF (LETADOT) THEN
    ZVERVEL(JLEV) = ( PETADOTDPDETA(JLEV) + PETADOTDPDETA(JLEV-1) ) / 2.0
  ELSE
    ZVERVEL(JLEV) = PVVEL0(JLEV)

  ENDIF

! upwind differencing

  IF ( ZVERVEL(JLEV) >= 0.0 ) THEN   ! downward motion
   IF ( JLEV .NE. 1 ) THEN
    ZQVADV(JLEV) = ZVERVEL(JLEV) * ( PQT0(JLEV) - PQT0(JLEV-1) ) / PDELP9(JLEV)
    ZTVADV(JLEV) = ZVERVEL(JLEV) * ( PTT0(JLEV) - PTT0(JLEV-1) ) / PDELP9(JLEV)
    ZUVADV(JLEV) = ZVERVEL(JLEV) * ( PUT0(JLEV) - PUT0(JLEV-1) ) / PDELP9(JLEV)
    ZVVADV(JLEV) = ZVERVEL(JLEV) * ( PVT0(JLEV) - PVT0(JLEV-1) ) / PDELP9(JLEV)
    ZAVADV(JLEV) = ZVERVEL(JLEV) * ( PAT0(JLEV) - PAT0(JLEV-1) ) / PDELP9(JLEV)
    ZLVADV(JLEV) = ZVERVEL(JLEV) * ( PLT0(JLEV) - PLT0(JLEV-1) ) / PDELP9(JLEV)
    ZIVADV(JLEV) = ZVERVEL(JLEV) * ( PIT0(JLEV) - PIT0(JLEV-1) ) / PDELP9(JLEV)
   ELSE
    ZQVADV(JLEV) = 0.0
    ZTVADV(JLEV) = 0.0
    ZUVADV(JLEV) = 0.0
    ZVVADV(JLEV) = 0.0
    ZAVADV(JLEV) = 0.0
    ZLVADV(JLEV) = 0.0
    ZIVADV(JLEV) = 0.0
   ENDIF
  ELSE                               ! upward motion
   IF ( JLEV .NE. NFLEVG ) THEN
    ZQVADV(JLEV) = ZVERVEL(JLEV) * ( PQT0(JLEV+1) - PQT0(JLEV) ) / PDELP9(JLEV)
    ZTVADV(JLEV) = ZVERVEL(JLEV) * ( PTT0(JLEV+1) - PTT0(JLEV) ) / PDELP9(JLEV)
    ZUVADV(JLEV) = ZVERVEL(JLEV) * ( PUT0(JLEV+1) - PUT0(JLEV) ) / PDELP9(JLEV)
    ZVVADV(JLEV) = ZVERVEL(JLEV) * ( PVT0(JLEV+1) - PVT0(JLEV) ) / PDELP9(JLEV)
    ZAVADV(JLEV) = ZVERVEL(JLEV) * ( PAT0(JLEV+1) - PAT0(JLEV) ) / PDELP9(JLEV)
    ZLVADV(JLEV) = ZVERVEL(JLEV) * ( PLT0(JLEV+1) - PLT0(JLEV) ) / PDELP9(JLEV)
    ZIVADV(JLEV) = ZVERVEL(JLEV) * ( PIT0(JLEV+1) - PIT0(JLEV) ) / PDELP9(JLEV)
   ELSE
    ZQVADV(JLEV) = 0.0
    ZTVADV(JLEV) = 0.0
    ZUVADV(JLEV) = 0.0
    ZVVADV(JLEV) = 0.0
    ZAVADV(JLEV) = 0.0
    ZLVADV(JLEV) = 0.0
    ZIVADV(JLEV) = 0.0
   ENDIF
  ENDIF

 ENDDO

ELSE

! centered in space differencing
! (unstable in combination with forward in time scheme
!  - can be fixed with Lax method)

  IF (LETADOT) THEN

! centered in space differencing - as in Eulerian IFS (Ritchie et al 1995)
! for etadotdpdeta input (on half levels)

    DO JLEV=2,NFLEVG-1
      ZQVADV(JLEV)=PETADOTDPDETA(JLEV)*(PQT0(JLEV+1)-PQT0(JLEV))
      ZQVADV(JLEV)=ZQVADV(JLEV)+PETADOTDPDETA(JLEV-1)*(PQT0(JLEV)-PQT0(JLEV-1))
      ZQVADV(JLEV)=ZQVADV(JLEV)/(2.0*PDELP9(JLEV))

      ZTVADV(JLEV)=PETADOTDPDETA(JLEV)*(PTT0(JLEV+1)-PTT0(JLEV))
      ZTVADV(JLEV)=ZTVADV(JLEV)+PETADOTDPDETA(JLEV-1)*(PTT0(JLEV)-PTT0(JLEV-1))
      ZTVADV(JLEV)=ZTVADV(JLEV)/(2.0*PDELP9(JLEV))

      ZUVADV(JLEV)=PETADOTDPDETA(JLEV)*(PUT0(JLEV+1)-PUT0(JLEV))
      ZUVADV(JLEV)=ZUVADV(JLEV)+PETADOTDPDETA(JLEV-1)*(PUT0(JLEV)-PUT0(JLEV-1))
      ZUVADV(JLEV)=ZUVADV(JLEV)/(2.0*PDELP9(JLEV))

      ZVVADV(JLEV)=PETADOTDPDETA(JLEV)*(PVT0(JLEV+1)-PVT0(JLEV))
      ZVVADV(JLEV)=ZVVADV(JLEV)+PETADOTDPDETA(JLEV-1)*(PVT0(JLEV)-PVT0(JLEV-1))
      ZVVADV(JLEV)=ZVVADV(JLEV)/(2.0*PDELP9(JLEV))

      ZAVADV(JLEV)=PETADOTDPDETA(JLEV)*(PAT0(JLEV+1)-PAT0(JLEV))
      ZAVADV(JLEV)=ZAVADV(JLEV)+PETADOTDPDETA(JLEV-1)*(PAT0(JLEV)-PAT0(JLEV-1))
      ZAVADV(JLEV)=ZAVADV(JLEV)/(2.0*PDELP9(JLEV))

      ZLVADV(JLEV)=PETADOTDPDETA(JLEV)*(PLT0(JLEV+1)-PLT0(JLEV))
      ZLVADV(JLEV)=ZLVADV(JLEV)+PETADOTDPDETA(JLEV-1)*(PLT0(JLEV)-PLT0(JLEV-1))
      ZLVADV(JLEV)=ZLVADV(JLEV)/(2.0*PDELP9(JLEV))

      ZIVADV(JLEV)=PETADOTDPDETA(JLEV)*(PIT0(JLEV+1)-PIT0(JLEV))
      ZIVADV(JLEV)=ZIVADV(JLEV)+PETADOTDPDETA(JLEV-1)*(PIT0(JLEV)-PIT0(JLEV-1))
      ZIVADV(JLEV)=ZIVADV(JLEV)/(2.0*PDELP9(JLEV))

    ENDDO

    ZQVADV(1)=0.5*(PETADOTDPDETA(1)*(PQT0(2)-PQT0(1)))/PDELP9(1)
    ZTVADV(1)=0.5*(PETADOTDPDETA(1)*(PTT0(2)-PTT0(1)))/PDELP9(1)
    ZUVADV(1)=0.5*(PETADOTDPDETA(1)*(PUT0(2)-PUT0(1)))/PDELP9(1)
    ZVVADV(1)=0.5*(PETADOTDPDETA(1)*(PVT0(2)-PVT0(1)))/PDELP9(1)

    ZAVADV(1)=0.5*(PETADOTDPDETA(1)*(PAT0(2)-PAT0(1)))/PDELP9(1)
    ZLVADV(1)=0.5*(PETADOTDPDETA(1)*(PLT0(2)-PLT0(1)))/PDELP9(1)
    ZIVADV(1)=0.5*(PETADOTDPDETA(1)*(PIT0(2)-PIT0(1)))/PDELP9(1)

    ZQVADV(NFLEVG)=0.5*(PETADOTDPDETA(NFLEVG-1)*(PQT0(NFLEVG)-PQT0(NFLEVG-1)))/PDELP9(NFLEVG)
    ZTVADV(NFLEVG)=0.5*(PETADOTDPDETA(NFLEVG-1)*(PTT0(NFLEVG)-PTT0(NFLEVG-1)))/PDELP9(NFLEVG)
    ZUVADV(NFLEVG)=0.5*(PETADOTDPDETA(NFLEVG-1)*(PUT0(NFLEVG)-PUT0(NFLEVG-1)))/PDELP9(NFLEVG)
    ZVVADV(NFLEVG)=0.5*(PETADOTDPDETA(NFLEVG-1)*(PVT0(NFLEVG)-PVT0(NFLEVG-1)))/PDELP9(NFLEVG)      

    ZAVADV(NFLEVG)=0.5*(PETADOTDPDETA(NFLEVG-1)*(PAT0(NFLEVG)-PAT0(NFLEVG-1)))/PDELP9(NFLEVG)
    ZLVADV(NFLEVG)=0.5*(PETADOTDPDETA(NFLEVG-1)*(PLT0(NFLEVG)-PLT0(NFLEVG-1)))/PDELP9(NFLEVG)
    ZIVADV(NFLEVG)=0.5*(PETADOTDPDETA(NFLEVG-1)*(PIT0(NFLEVG)-PIT0(NFLEVG-1)))/PDELP9(NFLEVG)

  ELSE

! pure centered space differencing - for omega input on full levels

    DO JLEV=2,NFLEVG-1
      ZQVADV(JLEV) = PVVEL0(JLEV) * ( PQT0(JLEV+1) - PQT0(JLEV-1) ) &
                   & / ( 2.0 * PDELP9(JLEV) )      
      ZTVADV(JLEV) = PVVEL0(JLEV) * ( PTT0(JLEV+1) - PTT0(JLEV-1) ) &
                   & / ( 2.0 * PDELP9(JLEV) )      
      ZUVADV(JLEV) = PVVEL0(JLEV) * ( PUT0(JLEV+1) - PUT0(JLEV-1) ) &
                   & / ( 2.0 * PDELP9(JLEV) )      
      ZVVADV(JLEV) = PVVEL0(JLEV) * ( PVT0(JLEV+1) - PVT0(JLEV-1) ) &
                   & / ( 2.0 * PDELP9(JLEV) )      
      ZAVADV(JLEV) = PVVEL0(JLEV) * ( PAT0(JLEV+1) - PAT0(JLEV-1) ) &
                   & / ( 2.0 * PDELP9(JLEV) )      
      ZLVADV(JLEV) = PVVEL0(JLEV) * ( PLT0(JLEV+1) - PLT0(JLEV-1) ) &
                   & / ( 2.0 * PDELP9(JLEV) )      
      ZIVADV(JLEV) = PVVEL0(JLEV) * ( PIT0(JLEV+1) - PIT0(JLEV-1) ) &
                   & / ( 2.0 * PDELP9(JLEV) )      
    ENDDO
    ZQVADV(1)      = PVVEL0(1) * ( PQT0(2) - PQT0(1) ) / ( 2.0 * PDELP9(1) )
    ZTVADV(1)      = PVVEL0(1) * ( PTT0(2) - PTT0(1) ) / ( 2.0 * PDELP9(1) )
    ZUVADV(1)      = PVVEL0(1) * ( PUT0(2) - PUT0(1) ) / ( 2.0 * PDELP9(1) )
    ZVVADV(1)      = PVVEL0(1) * ( PVT0(2) - PVT0(1) ) / ( 2.0 * PDELP9(1) )
    ZAVADV(1)      = PVVEL0(1) * ( PAT0(2) - PAT0(1) ) / ( 2.0 * PDELP9(1) )
    ZLVADV(1)      = PVVEL0(1) * ( PLT0(2) - PLT0(1) ) / ( 2.0 * PDELP9(1) )
    ZIVADV(1)      = PVVEL0(1) * ( PIT0(2) - PIT0(1) ) / ( 2.0 * PDELP9(1) )
    ZQVADV(NFLEVG) =  PVVEL0(NFLEVG) * ( PQT0(NFLEVG) - PQT0(NFLEVG-1) ) &
                   & / ( 2.0 * PDELP9(NFLEVG) )
    ZTVADV(NFLEVG) =  PVVEL0(NFLEVG) * ( PTT0(NFLEVG) - PTT0(NFLEVG-1) ) &
                   & / ( 2.0 * PDELP9(NFLEVG) )
    ZUVADV(NFLEVG) =  PVVEL0(NFLEVG) * ( PUT0(NFLEVG) - PUT0(NFLEVG-1) ) &
                   & / ( 2.0 * PDELP9(NFLEVG) )
    ZVVADV(NFLEVG) =  PVVEL0(NFLEVG) * ( PVT0(NFLEVG) - PVT0(NFLEVG-1) ) &
                   & / ( 2.0 * PDELP9(NFLEVG) )
    ZAVADV(NFLEVG) =  PVVEL0(NFLEVG) * ( PAT0(NFLEVG) - PAT0(NFLEVG-1) ) &
                   & / ( 2.0 * PDELP9(NFLEVG) )
    ZLVADV(NFLEVG) =  PVVEL0(NFLEVG) * ( PLT0(NFLEVG) - PLT0(NFLEVG-1) ) &
                   & / ( 2.0 * PDELP9(NFLEVG) )
    ZIVADV(NFLEVG) =  PVVEL0(NFLEVG) * ( PIT0(NFLEVG) - PIT0(NFLEVG-1) ) &
                   & / ( 2.0 * PDELP9(NFLEVG) )

  ENDIF


! Lax forward in time, centered in space differencing
! (this method adds a diffusion term to the centered space differencing)

  IF (LLAXMETHOD) THEN

    DO JLEV=2,NFLEVG-1
      ZQVADV(JLEV) = ZQVADV(JLEV) - &
                 & ( PQT0(JLEV+1) - 2 * PQT0(JLEV) + PQT0(JLEV-1) ) / 2 / PDT
      ZTVADV(JLEV) = ZTVADV(JLEV) - &
                 & ( PTT0(JLEV+1) - 2 * PTT0(JLEV) + PTT0(JLEV-1) ) / 2 / PDT
      ZUVADV(JLEV) = ZUVADV(JLEV) - &
                 & ( PUT0(JLEV+1) - 2 * PUT0(JLEV) + PUT0(JLEV-1) ) / 2 / PDT
      ZVVADV(JLEV) = ZVVADV(JLEV) - &
                 & ( PVT0(JLEV+1) - 2 * PVT0(JLEV) + PVT0(JLEV-1) ) / 2 / PDT
      ZAVADV(JLEV) = ZAVADV(JLEV) - &
                 & ( PAT0(JLEV+1) - 2 * PAT0(JLEV) + PAT0(JLEV-1) ) / 2 / PDT
      ZLVADV(JLEV) = ZLVADV(JLEV) - &
                 & ( PLT0(JLEV+1) - 2 * PLT0(JLEV) + PLT0(JLEV-1) ) / 2 / PDT
      ZIVADV(JLEV) = ZIVADV(JLEV) - &
                 & ( PIT0(JLEV+1) - 2 * PIT0(JLEV) + PIT0(JLEV-1) ) / 2 / PDT
    ENDDO
    ZQVADV(1)      = ZQVADV(1)      - ( PQT0(2)       - PQT0(1) )       / 2 / PDT
    ZQVADV(NFLEVG) = ZQVADV(NFLEVG) - ( PQT0(NFLEVG-1) - PQT0(NFLEVG) ) / 2 / PDT
    ZTVADV(1)      = ZTVADV(1)      - ( PTT0(2)       - PTT0(1) )       / 2 / PDT
    ZTVADV(NFLEVG) = ZTVADV(NFLEVG) - ( PTT0(NFLEVG-1) - PTT0(NFLEVG) ) / 2 / PDT
    ZUVADV(1)      = ZUVADV(1)      - ( PUT0(2)       - PUT0(1) )       / 2 / PDT
    ZUVADV(NFLEVG) = ZUVADV(NFLEVG) - ( PUT0(NFLEVG-1) - PUT0(NFLEVG) ) / 2 / PDT
    ZVVADV(1)      = ZVVADV(1)      - ( PVT0(2)       - PVT0(1) )       / 2 / PDT
    ZVVADV(NFLEVG) = ZVVADV(NFLEVG) - ( PVT0(NFLEVG-1) - PVT0(NFLEVG) ) / 2 / PDT
    ZAVADV(1)      = ZAVADV(1)      - ( PAT0(2)       - PAT0(1) )       / 2 / PDT
    ZAVADV(NFLEVG) = ZAVADV(NFLEVG) - ( PAT0(NFLEVG-1) - PAT0(NFLEVG) ) / 2 / PDT
    ZLVADV(1)      = ZLVADV(1)      - ( PLT0(2)       - PLT0(1) )       / 2 / PDT
    ZLVADV(NFLEVG) = ZLVADV(NFLEVG) - ( PLT0(NFLEVG-1) - PLT0(NFLEVG) ) / 2 / PDT
    ZIVADV(1)      = ZIVADV(1)      - ( PIT0(2)       - PIT0(1) )       / 2 / PDT
    ZIVADV(NFLEVG) = ZIVADV(NFLEVG) - ( PIT0(NFLEVG-1) - PIT0(NFLEVG) ) / 2 / PDT

  ENDIF      

ENDIF


!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF T TENDENCY.
!              --------------------------

IF (LDTWOTL) THEN

  !ZBETA=1.0_JPRB  ! fully implicit scheme
  ZBETA=0.5_JPRB  ! Crank-Nicolson scheme
  DO JLEV=1,NFLEVG
    ZCONST1=(PDT*(1._JPRB-ZBETA)*PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
    ZCONST2=(PDT*ZBETA          *PKAP(JLEV)*PVVEL0(JLEV))/PRESF(JLEV)
    ZCONST2=(1.0+ZCONST1)/(1.0-ZCONST2)
    ZTT1(JLEV)=ZCONST2*PTT0(JLEV)
    PTTEND(JLEV)=PTADV(JLEV)-ZTVADV(JLEV)
    PTTEND(JLEV)=PTTEND(JLEV)+(ZTT1(JLEV)-PTT0(JLEV))/PDT
  ENDDO

ELSE

  DO JLEV=1,NFLEVG
    PTTEND(JLEV)=PTADV(JLEV)-ZTVADV(JLEV)
    PTTEND(JLEV)=PTTEND(JLEV)+&
     &(PVVEL0(JLEV)*PKAP(JLEV)*PTT0(JLEV))/PRESF(JLEV)
  ENDDO

ENDIF


!     ------------------------------------------------------------------

!*       2.    MOMENTUM EQUATIONS.
!              -------------------


IF (LDTWOTL) THEN

  ZDTCORI=PDT*PRCORI
  ZBETA=0.5_JPRB ! Implicit factor (0 = explicit scheme, 0.5 = C-N scheme, 1 = fully implicit scheme)
  ZA   =ZDTCORI*(1._JPRB-ZBETA)
  ZB   =ZDTCORI*ZBETA
  ZALFA=1._JPRB/(1._JPRB + ZB*ZB)

  DO JLEV=1,NFLEVG

!      The Coriolis term is approximated in an implicit way
!      To solve: write U,V as complex variable in the Ekman equation,
!      eliminate and write components. 
    ZUT1(JLEV)=ZALFA*( (1._JPRB-ZA*ZB)*PUT0(JLEV) + ZDTCORI*PVT0(JLEV) &
             & + ZDTCORI*(ZB*PUG0(JLEV)-PVG0(JLEV)) )
  ! ZVT1(JLEV)=ZALFA*( (1._JPRB-ZA*ZB)*PVT0(JLEV) - ZDTCORI*PUT0(JLEV) &
  !          & + ZDTCORI*(ZB*PVG0(JLEV)+PUG0(JLEV)))
    ! Faster way substituting the previous result
    ZVT1(JLEV)=PVT0(JLEV) - ZA*PUT0(JLEV) - ZB*ZUT1(JLEV) &
             & + ZDTCORI*PUG0(JLEV)
  
    PUTEND(JLEV)=PUADV(JLEV)-ZUVADV(JLEV)
    PUTEND(JLEV)=PUTEND(JLEV)+(ZUT1(JLEV)-PUT0(JLEV))/PDT

    PVTEND(JLEV)=PVADV(JLEV)-ZVVADV(JLEV)
    PVTEND(JLEV)=PVTEND(JLEV)+(ZVT1(JLEV)-PVT0(JLEV))/PDT

  ENDDO

ELSE

  DO JLEV=1,NFLEVG
    PUTEND(JLEV)=PUADV(JLEV)-ZUVADV(JLEV)
    PUTEND(JLEV)=PUTEND(JLEV)+(PVT0(JLEV)-PVG0(JLEV))*PRCORI
    PVTEND(JLEV)=PVADV(JLEV)-ZVVADV(JLEV)
    PVTEND(JLEV)=PVTEND(JLEV)+(PUG0(JLEV)-PUT0(JLEV))*PRCORI
  ENDDO

ENDIF


!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF Q TENDENCY.
!        --------------------------------

DO JLEV=1,NFLEVG
  PQTEND(JLEV)=PQADV(JLEV)-ZQVADV(JLEV)
ENDDO


!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF A,L AND I TENDENCY.
!        ----------------------------------------

IF (LWADVCLD) THEN

  DO JLEV=1,NFLEVG
    PATEND(JLEV)=-ZAVADV(JLEV)
    PLTEND(JLEV)=-ZLVADV(JLEV)
    PITEND(JLEV)=-ZIVADV(JLEV)
!   PATEND(JLEV)=PAADV(JLEV)-ZAVADV(JLEV)
!   PLTEND(JLEV)=PLADV(JLEV)-ZLVADV(JLEV) 
!   PITEND(JLEV)=PIADV(JLEV)-ZIVADV(JLEV)
  ENDDO

ELSE

  DO JLEV=1,NFLEVG
    PATEND(JLEV)=0.0
    PLTEND(JLEV)=0.0
    PITEND(JLEV)=0.0
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPDYN1C',1,ZHOOK_HANDLE)


END SUBROUTINE CPDYN1C
