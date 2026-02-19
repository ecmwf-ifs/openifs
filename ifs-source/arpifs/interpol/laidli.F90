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

SUBROUTINE LAIDLI(KPROMA,KPROMB,KSTART,KPROF,KFLEV,&
 & KFLDN,KFLDX,&
 & PDLAT,PDLO,KL0,&
 & PXSL,PXF)  

!**** *LAIDLI  -  semi-LAgrangian scheme:
!                 Bilinear horizontal interpolations for one variable.

!     Purpose.
!     --------
!       Performs bilinear horizontal interpolations for one variable.

!**   Interface.
!     ----------
!        *CALL* *LAIDLI(...)

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
!          PDLAT   - weight (distance) for horizontal linear interpolation
!                    on a same latitude.
!          PDLO    - weights (distances) for horizontal linear interpolation
!                    on a same longitude.
!          KL0     - indices of the four western points
!                    of the 16 points interpolation grid.
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
!        K. YESSAD, after the subroutines LAGINL3
!        written by Maurice IMBARD, Alain CRAPLET and Michel ROCHAS
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!        Original : FEBRUARY 1992.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad (Sep 2008): update comments + cleanings.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(KPROMB,KFLEV,1:2) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMB,KFLEV,1:2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KPROMA*(KFLDX-KFLDN+1)) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXF(KPROMB,KFLEV) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZDLAT, ZDLO1, ZDLO2, ZVALLO1, ZVALLO2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAIDLI',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    INTERPOLATIONS.
!              ---------------

DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF

!     Computation of coordinates and distances.

    ZDLAT =PDLAT(JROF,JLEV)

    ZDLO1 =PDLO(JROF,JLEV,1)
    ZDLO2 =PDLO(JROF,JLEV,2)

!     Interpolation.

    ZVALLO1=PXSL(KL0(JROF,JLEV,1)+1) + ZDLO1*&
     & ( PXSL(KL0(JROF,JLEV,1)+2)-PXSL(KL0(JROF,JLEV,1)+1) )  
    ZVALLO2=PXSL(KL0(JROF,JLEV,2)+1) + ZDLO2*&
     & ( PXSL(KL0(JROF,JLEV,2)+2)-PXSL(KL0(JROF,JLEV,2)+1) )  
    PXF(JROF,JLEV)= ZVALLO1 + ZDLAT*(ZVALLO2-ZVALLO1)

  ENDDO
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LAIDLI',1,ZHOOK_HANDLE)
END SUBROUTINE LAIDLI

