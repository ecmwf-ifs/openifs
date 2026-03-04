! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE VDFDIFMS (KIDIA, KFDIA, KLON, KLEV, KTOP,&
 & PTMST ,PVDIFTS,PUM1 ,PVM1  ,PAPHM1 , PCFM, PTOFDC,&
 & PVOM ,PVOL ,PUDIF ,PVDIF)
!     ------------------------------------------------------------------

!**   *VDFDIFMS* - DOES THE IMPLICIT CALCULATION FOR MOMENTUM DIFFUSION
!                  (Nonlinear version for trajectory in adjoint)

!     PURPOSE
!     -------

!     SOLVE TRIDIAGONAL MATRICES FOR MOMENTUM DIFFUSION

!     INTERFACE
!     ---------

!     *VDFDIFMS* IS CALLED BY *VDFMAINS*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KTOP*         INDEX FOR BOUNDARY LAYER TOP

!     INPUT PARAMETERS (REAL):

!     *PUM1*        X-VELOCITY COMPONENT AT T-1 (Trajectory)
!     *PVM1*        Y-VELOCITY COMPONENT AT T-1 (Trajectory)
!     *PAPHM1*      PRESSURE AT T-1             (Trajectory)
!     *PCFM*        PROP.TO EXCH.COEFF. FOR MOMENTUM (C,K-STAR IN DOC.)
!                                                (Trajectory)
!     *PTOFDC*      TURB. OROGR. DRAG (IMPLICIT) COEFFICIENTS (Trajectory)
!     *PVOM*        U-TENDENCY                  (Trajectory)
!     *PVOL*        V-TENDENCY                  (Trajectory)

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)

!     OUTPUT PARAMETERS (REAL):

!     *PUDIF*       U-DOUBLE-TILDE DEVIDED BY ALFA (Trajectory)
!     *PVDIF*       V-DOUBLE-TILDE DEVIDED BY ALFA (Trajectory)

!     METHOD
!     ------

!     *LU*-DECOMPOSITION AND BACK SUBSTITUTION IN ONE DOWNWARD SCAN
!     AND ONE UPWARD SCAN.

!     AUTHOR.
!     -------
!      J.F. MAHFOUF          E.C.M.W.F.    02/10/95

!     MODIFICATIONS.
!     --------------
!      P. LOPEZ    E.C.M.W.F.  25/02/05 (nonlinear version for trajectory in adjoint)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOFDC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIF(KLON,KLEV) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) ::  ZTCOE(KLON), ZEBSM(KLON,KLEV)

INTEGER(KIND=JPIM) :: ILEVM1, ITOPP1, JK, JL

REAL(KIND=JPRB) :: ZDISC, ZDISCI, ZFAC, ZQDP, ZTPFAC2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VDFDIFMS',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

ZTPFAC2 = 1.0_JPRB/PVDIFTS
ILEVM1 = -1+KLEV
ITOPP1 = 1+KTOP

!*         1.1     SETTING OF RIGHT HAND SIDES.

DO JK = KTOP, KLEV
  DO JL = KIDIA, KFDIA
    PUDIF(JL,JK) = ZTPFAC2*PUM1(JL,JK)+PTMST*PVOM(JL,JK)
    PVDIF(JL,JK) = ZTPFAC2*PVM1(JL,JK)+PTMST*PVOL(JL,JK)
  ENDDO
ENDDO

!*         1.2     TOP LAYER ELIMINATION.

DO JL = KIDIA, KFDIA
  ZTCOE(JL) = PCFM(JL,KTOP)
  ZQDP = 1.0_JPRB/(PAPHM1(JL,ITOPP1)-PAPHM1(JL,KTOP))
  ZDISCI=1.0_JPRB+PCFM(JL,KTOP)*ZQDP+PTOFDC(JL,KTOP)
  ZDISC=1.0_JPRB/ZDISCI
  ZEBSM(JL,KTOP) = ZDISC*(PCFM(JL,KTOP)*ZQDP)
  PUDIF(JL,KTOP) = ZDISC*PUDIF(JL,KTOP)
  PVDIF(JL,KTOP) = ZDISC*PVDIF(JL,KTOP)
ENDDO

!*         1.3     ELIMINATION FOR MIDDLE LAYERS.

DO JK = ITOPP1, ILEVM1
  DO JL = KIDIA, KFDIA
    ZQDP = 1.0_JPRB/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
    ZFAC = ZTCOE(JL)*ZQDP
    ZTCOE(JL) = PCFM(JL,JK)
    ZDISCI=1.0_JPRB+ZFAC*(1.0_JPRB-ZEBSM(JL,JK-1)) &
     &    +PCFM(JL,JK)*ZQDP+PTOFDC(JL,JK)
    ZDISC=1.0_JPRB/ZDISCI
    ZEBSM(JL,JK) = ZDISC*(PCFM(JL,JK)*ZQDP)
    PUDIF(JL,JK) =  ZDISC*(PUDIF(JL,JK)+ZFAC*PUDIF(JL,JK-1))
    PVDIF(JL,JK) =  ZDISC*(PVDIF(JL,JK)+ZFAC*PVDIF(JL,JK-1))
  ENDDO
ENDDO

!*         1.4     BOTTOM LAYER ELIMINATION.

DO JL = KIDIA, KFDIA
  ZQDP = 1.0_JPRB/(PAPHM1(JL,KLEV+1)-PAPHM1(JL,KLEV))
  ZFAC = ZTCOE(JL)*ZQDP
  ZDISCI=1.0_JPRB+ZFAC*(1.0_JPRB-ZEBSM(JL,ILEVM1)) &
   &    +PCFM(JL,KLEV)*ZQDP+PTOFDC(JL,KLEV)
  ZDISC=1.0_JPRB/ZDISCI
  PUDIF(JL,KLEV) =  ZDISC*(PUDIF(JL,KLEV)+ZFAC*PUDIF(JL,ILEVM1))
  PVDIF(JL,KLEV) =  ZDISC*(PVDIF(JL,KLEV)+ZFAC*PVDIF(JL,ILEVM1))
ENDDO

!*         1.5     BACK-SUBSTITUTION.

DO JK = ILEVM1, KTOP, -1
  DO JL = KIDIA, KFDIA
    PUDIF(JL,JK) = PUDIF(JL,JK)+ZEBSM(JL,JK)*PUDIF(JL,JK+1)
    PVDIF(JL,JK) = PVDIF(JL,JK)+ZEBSM(JL,JK)*PVDIF(JL,JK+1)
  ENDDO
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VDFDIFMS',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDIFMS

