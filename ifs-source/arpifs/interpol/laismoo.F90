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

SUBROUTINE LAISMOO(KPROMA,KPROMB,KSTART,KPROF,KFLEV,&
 & KFLDN,KFLDX,PDLAT,PDLO,KL0,PDVER,&
 & PXSL,PXF)  

!**** *LAISMOO  -  semi-LAgrangian scheme:
!                 Smoothing interpolation by linear least-square fit

!     Purpose.
!     --------
!       Performs smooting in the horizontal and linear in the vertical interpolation

!**   Interface.
!     ----------
!        *CALL* *LAISMOO(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA  - horizontal dimension for grid-point quantities.
!          KPROMB  - horizontal dimension for interpolation point
!                    quantities.
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          KFLEV   - vertical dimension.
!          KFLDN   - number of the first field.
!          KFLDX   - number of the last field.
!          PDLAT   - distance for horizontal interpolation
!          PDLO    - distances for horizontal interpolations
!          KL0     - index of the four western points
!                    of the 16 points interpolation grid.
!          PDVER   - weights for linear vertical interpolation
!          PXSL    - quantity to be interpolated.

!        OUTPUT:
!          PXF     - interpolated variable.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!        No external.

!     Reference.
!     ----------

!     Author.
!     -------
!       M.Hortal   ECMWF

!     Modifications.
!     --------------
!       Original : Dec 2003
!       K. Yessad (Sep 2008): update comments + cleanings.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPIA
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLAT(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(KPROMB,KFLEV,0:3) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMB,KFLEV,0:3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVER(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KPROMA*(KFLDX-KFLDN+1)) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXF(KPROMB,KFLEV) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIA) :: IV1L0, IV1L1, IV1L2, IV1L3,&
 & IV2L0, IV2L1, IV2L2, IV2L3  
INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZPOI(0:3,0:3,2)

REAL(KIND=JPRB) :: ZDVER, ZIS, ZISLO, ZISLO1, ZISLO2, ZISLO3, &
 & ZSI, ZSILO, ZSILO1, ZSILO2, ZSILO3  

REAL(KIND=JPRB) :: PD
REAL(KIND=JPRB) :: FSLO1
REAL(KIND=JPRB) :: FSLO2
REAL(KIND=JPRB) :: FSLO3
REAL(KIND=JPRB) :: FSLO4

REAL(KIND=JPRB),PARAMETER :: Z10_R=1.0_JPRB/10._JPRB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

FSLO1(PD)=(4._JPRB-3._JPRB*PD)*Z10_R
FSLO2(PD)=(3._JPRB-PD)*Z10_R
FSLO3(PD)=(2.0_JPRB+PD)*Z10_R
FSLO4(PD)=(1.0_JPRB+3._JPRB*PD)*Z10_R

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAISMOO',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    INTERPOLATIONS.
!              ---------------

IV1L0=  KPROMA
IV1L1=1+KPROMA
IV1L2=2+KPROMA
IV1L3=3+KPROMA
IV2L0=IV1L0+KPROMA
IV2L1=IV1L1+KPROMA
IV2L2=IV1L2+KPROMA
IV2L3=IV1L3+KPROMA

DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF

    ZPOI(1,1,1)=PXSL(KL0(JROF,JLEV,1)+IV1L1)
    ZPOI(2,1,1)=PXSL(KL0(JROF,JLEV,1)+IV1L2)
    ZPOI(1,2,1)=PXSL(KL0(JROF,JLEV,2)+IV1L1)
    ZPOI(2,2,1)=PXSL(KL0(JROF,JLEV,2)+IV1L2)
    ZPOI(1,1,2)=PXSL(KL0(JROF,JLEV,1)+IV2L1)
    ZPOI(2,1,2)=PXSL(KL0(JROF,JLEV,1)+IV2L2)
    ZPOI(1,2,2)=PXSL(KL0(JROF,JLEV,2)+IV2L1)
    ZPOI(2,2,2)=PXSL(KL0(JROF,JLEV,2)+IV2L2)

    ZDVER =PDVER(JROF,JLEV)

!     32 points smoothing interpolation.

    ZSILO =PXSL(KL0(JROF,JLEV,0)+IV1L0)*FSLO1(PDLO(JROF,JLEV,0)) + &
     & PXSL(KL0(JROF,JLEV,0)+IV1L1)*FSLO2(PDLO(JROF,JLEV,0)) + &
     & PXSL(KL0(JROF,JLEV,0)+IV1L2)*FSLO3(PDLO(JROF,JLEV,0)) + &
     & PXSL(KL0(JROF,JLEV,0)+IV1L3)*FSLO4(PDLO(JROF,JLEV,0))  
    ZSILO1=PXSL(KL0(JROF,JLEV,1)+IV1L0)*FSLO1(PDLO(JROF,JLEV,1)) + &
     & PXSL(KL0(JROF,JLEV,1)+IV1L1)*FSLO2(PDLO(JROF,JLEV,1)) + &
     & PXSL(KL0(JROF,JLEV,1)+IV1L2)*FSLO3(PDLO(JROF,JLEV,1)) + &
     & PXSL(KL0(JROF,JLEV,1)+IV1L3)*FSLO4(PDLO(JROF,JLEV,1))   
    ZSILO2=PXSL(KL0(JROF,JLEV,2)+IV1L0)*FSLO1(PDLO(JROF,JLEV,2)) + &
     & PXSL(KL0(JROF,JLEV,2)+IV1L1)*FSLO2(PDLO(JROF,JLEV,2)) + &
     & PXSL(KL0(JROF,JLEV,2)+IV1L2)*FSLO3(PDLO(JROF,JLEV,2)) + &
     & PXSL(KL0(JROF,JLEV,2)+IV1L3)*FSLO4(PDLO(JROF,JLEV,2))   
    ZSILO3=PXSL(KL0(JROF,JLEV,3)+IV1L0)*FSLO1(PDLO(JROF,JLEV,3)) + &
     & PXSL(KL0(JROF,JLEV,3)+IV1L1)*FSLO2(PDLO(JROF,JLEV,3)) + &
     & PXSL(KL0(JROF,JLEV,3)+IV1L2)*FSLO3(PDLO(JROF,JLEV,3)) + &
     & PXSL(KL0(JROF,JLEV,3)+IV1L3)*FSLO4(PDLO(JROF,JLEV,3))  

    ZISLO =PXSL(KL0(JROF,JLEV,0)+IV2L0)*FSLO1(PDLO(JROF,JLEV,0)) + &
     & PXSL(KL0(JROF,JLEV,0)+IV2L1)*FSLO2(PDLO(JROF,JLEV,0)) + &
     & PXSL(KL0(JROF,JLEV,0)+IV2L2)*FSLO3(PDLO(JROF,JLEV,0)) + &
     & PXSL(KL0(JROF,JLEV,0)+IV2L3)*FSLO4(PDLO(JROF,JLEV,0))  
    ZISLO1=PXSL(KL0(JROF,JLEV,1)+IV2L0)*FSLO1(PDLO(JROF,JLEV,1)) + &
     & PXSL(KL0(JROF,JLEV,1)+IV2L1)*FSLO2(PDLO(JROF,JLEV,1)) + &
     & PXSL(KL0(JROF,JLEV,1)+IV2L2)*FSLO3(PDLO(JROF,JLEV,1)) + &
     & PXSL(KL0(JROF,JLEV,1)+IV2L3)*FSLO4(PDLO(JROF,JLEV,1))  
    ZISLO2=PXSL(KL0(JROF,JLEV,2)+IV2L0)*FSLO1(PDLO(JROF,JLEV,2)) + &
     & PXSL(KL0(JROF,JLEV,2)+IV2L1)*FSLO2(PDLO(JROF,JLEV,2)) + &
     & PXSL(KL0(JROF,JLEV,2)+IV2L2)*FSLO3(PDLO(JROF,JLEV,2)) + &
     & PXSL(KL0(JROF,JLEV,2)+IV2L3)*FSLO4(PDLO(JROF,JLEV,2))  
    ZISLO3=PXSL(KL0(JROF,JLEV,3)+IV2L0)*FSLO1(PDLO(JROF,JLEV,3)) + &
     & PXSL(KL0(JROF,JLEV,3)+IV2L1)*FSLO2(PDLO(JROF,JLEV,3)) + &
     & PXSL(KL0(JROF,JLEV,3)+IV2L2)*FSLO3(PDLO(JROF,JLEV,3)) + &
     & PXSL(KL0(JROF,JLEV,3)+IV2L3)*FSLO4(PDLO(JROF,JLEV,3))  

    ZSI   = ZSILO *FSLO1(PDLAT(JROF,JLEV)) + &
     & ZSILO1*FSLO2(PDLAT(JROF,JLEV)) + &
     & ZSILO2*FSLO3(PDLAT(JROF,JLEV)) + &
     & ZSILO3*FSLO4(PDLAT(JROF,JLEV))  
    ZIS   = ZISLO *FSLO1(PDLAT(JROF,JLEV)) + &
     & ZISLO1*FSLO2(PDLAT(JROF,JLEV)) + &
     & ZISLO2*FSLO3(PDLAT(JROF,JLEV)) + &
     & ZISLO3*FSLO4(PDLAT(JROF,JLEV))  

    PXF(JROF,JLEV)=  ZSI + ZDVER*(ZIS - ZSI)

  ENDDO
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LAISMOO',1,ZHOOK_HANDLE)
END SUBROUTINE LAISMOO

