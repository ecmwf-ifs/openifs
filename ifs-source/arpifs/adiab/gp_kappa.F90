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

SUBROUTINE GP_KAPPA(&
 ! --- INPUT --------------------------------------------------
 & YDVAB, YDDYN,LDSLHD_STATIC,KPROMA,KSTART,KPROF,KFLEV,&
 & PDTS2,PUT0L,PVT0L,PVORT0,PDIVT0,&
 & PLHU0,PLHV0,PRES0,PSPT0L,PSPT0M,PRDELP,PSLHDA,PSLHDD0,&
 ! --- OUTPUT -------------------------------------------------
 & PKAPPA)

!**** *GP_KAPPA* input for SLHD diffusion.
!                Computation of the kappa function ("coefficient
!                of diffusion") based on rescaled horizontal
!                deformation of the flow.
!                Used by SLHD scheme.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *GP_KAPPA(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA      - horizontal dimension.
!          KSTART      - first element of work.
!          KPROF       - depth of work.
!          KFLEV       - number of layers.
!          PDTS2       - 0.5*"pdt", where "pdt" =
!                        time step for the first time-integration step of
!                        a leap-frog scheme or all time-integration steps of
!                        a two-time level scheme; 2*time step for the following
!                        time-integration steps of a leap-frog scheme.
!          PUT0L       - zonal derivative of U-wind at time t.
!          PVT0L       - zonal derivative of V-wind at time t.
!          PVORT0      - vorticity at time t.
!          PDIVT0      - divergence at t.
!          PLHU0       - zonal wind time t at half levels.
!          PLHV0       - meridian wind time t at half levels.
!          PRES0       - prehyds at t.
!          PSPT0L      - zonal derivative of "ln(prehyds)" at t.
!          PSPT0M      - merid derivative of "ln(prehyds)" at t.
!          PRDELP      - 1/(pressure depth of layers) at t.
!          PSLHDA      - Scaling factor of the deformation in f(d) function
!                        (including the model resolution correction)
!          PSLHDD0     - Treshold for deformation tensor enhancement

!        OUTPUT:
!          PKAPPA     - Kappa function ("coefficient of SLHD")
!                       based on the rescaled horizontal deformation
!                       of the flow evaluated at instant "t".

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LATTE_KAPPA.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!  Author.
!  -------
!    Original (F. VANA and K. YESSAD): AUGUST 2003.

!  Modifications.
!  --------------
!   19-Aug-2003 F. Vana     horizontal flow deformation for SLHD
!   01-Oct-2003 M. Hamrud   CY28 Cleaning
!   25-Oct-2004 F. Vana     New content of output + DFI fix -> renamed
!   07-Jun-2006 F. Vana     SLHDA,SLHDD0 changed to arays P...
!   30-Jun-2008 J. Masek    New SLHD interpolators (removed SLHDKMAX)
!   27-Aug-2008 F. Vana     LSLHD_STATIC case
!   27-Aug-2008 K. Yessad   rationalisation of dummy argument + renaming
!   09-Sep-2008 J. Masek    Flow deformation along pressure levels
!   21-May-2019 F. Vana     Computation restricted to levels where required.
! End Modifications
!     ------------------------------------------------------------------

USE YOMVERT  , ONLY : TVAB
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMDYN   , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
LOGICAL           ,INTENT(IN)    :: LDSLHD_STATIC
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0L(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0L(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVORT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLHU0(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLHV0(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0L(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0M(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKAPPA(KPROMA,KFLEV) 

!     -------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV,JROF
REAL(KIND=JPRB) :: ZABSDTS2
REAL(KIND=JPRB) :: ZC,ZAUX,ZDU,ZDV,ZUX,ZUY,ZVX,ZVY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GP_KAPPA',0,ZHOOK_HANDLE)
ASSOCIATE(SLHDB=>YDDYN%SLHDB, SLHDDIV=>YDDYN%SLHDDIV, SLHDHOR=>YDDYN%SLHDHOR, &
 & SLHDRATDDIV=>YDDYN%SLHDRATDDIV,NLEV_SPONGE=>YDDYN%NLEV_SPONGE)
!    --------------------------------------------------------

!*       0.   Basic settings

!  Positive time step (required in DFI)
ZABSDTS2=ABS(PDTS2)

!*       1.   Kappa function computation

IF (LDSLHD_STATIC) THEN
  PKAPPA(KSTART:KPROF,1:KFLEV)=1.0_JPRB
ELSE

  DO JLEV=1,NLEV_SPONGE
    ZC=0.5_JPRB*SLHDHOR*(YDVAB%VBH(JLEV)+YDVAB%VBH(JLEV-1))
!DEC$ IVDEP
    DO JROF=KSTART,KPROF

      ! diagnose wind derivatives along pressure levels (quasi horizontal)
      ZAUX=ZC*PRDELP(JROF,JLEV)*PRES0(JROF)
      ZDU=PLHU0(JROF,JLEV)-PLHU0(JROF,JLEV-1)
      ZDV=PLHV0(JROF,JLEV)-PLHV0(JROF,JLEV-1)
      ZUX=PUT0L(JROF,JLEV)-ZAUX*PSPT0L(JROF)*ZDU
      ZVX=PVT0L(JROF,JLEV)-ZAUX*PSPT0L(JROF)*ZDV
      ZUY=PVT0L(JROF,JLEV)-PVORT0(JROF,JLEV)-ZAUX*PSPT0M(JROF)*ZDU
      ZVY=PDIVT0(JROF,JLEV)-PUT0L(JROF,JLEV)-ZAUX*PSPT0M(JROF)*ZDV

      ! compute kappa function
      PKAPPA(JROF,JLEV)=(1.0_JPRB-SLHDDIV)*SQRT((ZUX-ZVY)**2+(ZVX+ZUY)**2)+&
       & SLHDDIV*SLHDRATDDIV*ABS(ZUX+ZVY)
      PKAPPA(JROF,JLEV)=2.0_JPRB*ZABSDTS2*&
        & PSLHDA(JROF)*PKAPPA(JROF,JLEV)*&
        & (MAX(1.0_JPRB,PKAPPA(JROF,JLEV)/PSLHDD0(JROF)))**SLHDB
      PKAPPA(JROF,JLEV)=PKAPPA(JROF,JLEV)/(PKAPPA(JROF,JLEV)+1.0_JPRB)

    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_KAPPA',1,ZHOOK_HANDLE)
END SUBROUTINE GP_KAPPA

