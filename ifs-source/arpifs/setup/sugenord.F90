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

SUBROUTINE SUGENORD(KSTTYP,KFPRGP,PSTRET,PMUCEN,PLOCEN,&
 & PSINLAT,PSINLAM,PCOSLAM,PFPNORX,PFPNORY)  

!-----------------------------------------------------------------------
!**** *SUGENORD*  - INITIALIZE DIRECTION OF GEOGRAPHICAL NORTH

!     PURPOSE.
!     --------
!        COMPUTES THE ROTATION MATRIX TO FIND THE GEOGRAPHICAL NORTH, 
!        FROM GRID POINTS ON A TRANSFORMED SHERE.

!**   INTERFACE. *CALL* *SUGENORD*
!     ----------
!        EXPLICIT ARGUMENTS :
!        ------------------
!         KSTTYP           : type of transformation (1 : no rotation ; 2 : rotation) 
!         KFPRGP           : number of points 
!         PSTRET           : stretching factor
!         PMUCEN           : sine of geographical latitude of pole of interest
!         PLOCEN           : geographical longitude of pole of interest
!         PSINLAT          : sine of pseudo-latitude of points
!         PSINLAM, PCOSLAM : sine and cosine of pseudo-longitudes of points
!         PFPNORX          : rotation matrix (sin(alfa))        
!         PFPNORY          : rotation matrix (cos(alfa))       

!        IMPLICIT ARGUMENTS   : none
!        ------------------

!     EXTERNALS.  None
!     ----------
!     AUTHOR.    VINCENT CASSE, RYAD EL KHATIB *METEO-FRANCE*
!     -------
!     MODIFICATIONS.
!     --------------
!        970304 ORIGINAL from previous SUFPG
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KFPRGP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTTYP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRET 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMUCEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLAT(KFPRGP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLAM(KFPRGP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOSLAM(KFPRGP) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPNORX(KFPRGP) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPNORY(KFPRGP) 
INTEGER(KIND=JPIM) :: J

REAL(KIND=JPRB) :: ZD, ZL2C, ZLA, ZLAP, ZLB, ZLC, ZLC2M1, ZLC2P1,&
 & ZLCDELA, ZLCLAC, ZLCLAP, ZLCLOC, ZLCLOP, &
 & ZLD, ZLSDELA, ZLSLAC, ZLSLAP, ZLSLOC, ZLSLOP, &
 & ZLXX, ZX, ZY  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

!*       1. COMPUTE ROTATION
!           ----------------

IF (LHOOK) CALL DR_HOOK('SUGENORD',0,ZHOOK_HANDLE)
IF (KSTTYP == 1) THEN
  DO J = 1, KFPRGP
    PFPNORX(J) = 0.0_JPRB
    PFPNORY(J) = 1.0_JPRB
  ENDDO
ELSEIF(KSTTYP == 2) THEN
  ZLC    = PSTRET
  ZL2C   = 2.0_JPRB * PSTRET
  ZLC2P1 = PSTRET * PSTRET + 1.0_JPRB
  ZLC2M1 = PSTRET * PSTRET - 1.0_JPRB
  ZLSLAP = PMUCEN
  ZLCLAP = SQRT(1.0_JPRB - ZLSLAP*ZLSLAP)
  DO J = 1, KFPRGP
    ZLSLAC = PSINLAT(J)
    ZLA    = ZLC2P1 + ZLC2M1 * ZLSLAC
    ZLB    = ZLC2M1 + ZLC2P1 * ZLSLAC
    ZLCLAC = SQRT(1.0_JPRB - ZLSLAC* ZLSLAC)
    ZLSLOC = PSINLAM(J)
    ZLCLOC = PCOSLAM(J)
    ZLXX   = ZL2C*ZLCLAP*ZLCLAC*ZLCLOC + ZLB*ZLSLAP
    IF (ABS(ZLXX)  >=  ABS(ZLA)) THEN
      ZLAP = ZLC2P1 + ZLC2M1 * ZLSLAP
      ZLD  = - ZL2C * ZLCLAP / (ZLAP * ZLCLAC)
      ZLSLOP = SIN(PLOCEN)
      ZLCLOP = COS(PLOCEN)
      PFPNORX(J) = - ZLD * (ZLSLOC*ZLCLOP-ZLCLOC*ZLSLOP)
      PFPNORY(J) = ZLD * (ZLCLOC*ZLCLOP+ZLSLOC*ZLSLOP)
    ELSE
      ZLD = 1.0_JPRB / SQRT(ZLA*ZLA-ZLXX*ZLXX)
      PFPNORX(J) = - ZLCLAP * ZLSLOC * ZLA * ZLD
      PFPNORY(J) = (ZL2C*ZLSLAP*ZLCLAC-ZLB*ZLCLAP*ZLCLOC)*ZLD
    ENDIF
  ENDDO
ELSEIF(KSTTYP == 3) THEN
!       Computation with respect to a given reference (real or pseudo-geography) 
  ZLSLAP = PMUCEN
  ZLCLAP = SQRT(1.0_JPRB - ZLSLAP*ZLSLAP)
  ZLSLOP = SIN(PLOCEN)
  ZLCLOP = COS(PLOCEN)
  DO J = 1, KFPRGP
    ZLSLOC = PSINLAM(J)
    ZLCLOC = PCOSLAM(J)
    ZLSLAC = PSINLAT(J)
    ZLCLAC = SQRT(1.0_JPRB - ZLSLAC* ZLSLAC)
    ZLSDELA = ZLSLOC*ZLCLOP-ZLCLOC*ZLSLOP
    ZLCDELA = ZLCLOC*ZLCLOP+ZLSLOC*ZLSLOP
    ZX = ZLSLAC*ZLCLAP - ZLCLAC*ZLSLAP*ZLCDELA
    ZY = ZLCLAC*ZLSDELA
    ZD=1.0_JPRB / SQRT(ZX*ZX + ZY*ZY)
    PFPNORX(J) = ZLCLAP*ZLSDELA*ZD
    PFPNORY(J) = (ZLSLAP*ZLCLAC - ZLSLAC*ZLCLAP*ZLCDELA)*ZD
  ENDDO
ENDIF

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGENORD',1,ZHOOK_HANDLE)
END SUBROUTINE SUGENORD

