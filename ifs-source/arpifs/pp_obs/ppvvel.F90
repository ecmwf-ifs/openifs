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

SUBROUTINE PPVVEL(YDVETA,YDCVER,KPROMA,KST,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,KLEVB,&
 & PRPRES,LDBELO,LDBLOW,PRXP,PRXPD,PRDELP,PEVEL,PVVEL,&
 & PVVPP,LDETADOT)  

!**** *PPVVEL* - POST-PROCESS VERTICAL VELOCITY OMEGA OR ETADOT.

!     PURPOSE.
!     --------
!           COMPUTES VERTICAL VELOCITY OMEGA OR ETADOT AT FULL LEVELS AND
!       PERFORMS THE VERTICAL INTERPOLATION TO A GIVEN MODEL CO-ORDINATE LEVEL.
!       SETS OMEGA OR ETADOT TO 0 IF PRESSURE IS BIGGER
!       THAN THE PRESSURE OF THE LOWER FULL LEVEL .

!**   INTERFACE.
!     ----------
!       *CALL*  *PPVVEL(...)

!        EXPLICIT ARGUMENTS
!        --------------------

!        * INPUT:
!         KPROMA    - horizontal dimension.
!         KST       - start of work.
!         KPROF     - depth of work.
!         KFLEV     - number of input pressure levels.
!         KLEVP     - number of output pressure levels.
!         KLOLEV    - beginning for the interpolation.
!         KPPM     - Number of interpolation methods in post-processing (INPUT)
!         KLEVB     - input level below PRPRES (see routine PPFLEV).

!         PRPRES    - post-processing level pressures.
!         LDBELO    - .TRUE. if pressure is under lowest model level.
!         LDBLOW    - .TRUE. if LOBELO(J) is containing at least one .TRUE.
!         PRXP      - contains "prehyd" and ln(prehyd) at half and full levels.
!                     (see routine PPINIT).
!         PRXPD     - contains 1/[Delta prehyd] and  1/[Delta log(prehyd)].
!         PRDELP    - contains  1/[Delta prehyd] at full levels.
!         PEVEL     - "etadot (d prehyd / d eta)" at half levels.
!         PVVEL     - "(omega / prehyd)" at full levels.

!        * OUTPUT:
!         PVVPP     - post-processed vertical velocity.

!        * OPTIONAL INPUT:
!         LDETADOT  - compute "etadot" instead of "omega".

!        IMPLICIT ARGUMENTS :
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!    ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*
!      ORIGINAL : 89-01-26

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      N. Wedi      : add etadot as possible output
!      K. Yessad (Dec 2008) : remove useless dummy arguments.
!      K. Yessad (Feb 2009) : remove call to GPCTY.
!     ------------------------------------------------------------------

USE YOMVERT  , ONLY : TVETA
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCVER  , ONLY : TCVER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVETA)       ,INTENT(IN)    :: YDVETA
TYPE(TCVER)       ,INTENT(IN)    :: YDCVER
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVB(KPROMA,KLEVP,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRES(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBELO(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBLOW(KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXP(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXPD(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVVEL(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVVPP(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN),OPTIONAL    :: LDETADOT

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZVVEL(KPROMA,0:KFLEV)

INTEGER(KIND=JPIM) :: ISLCT, JL, JLEV, JLEVP

LOGICAL :: LL_ETAD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "ppintp.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPVVEL',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    POST-PROCESS VERTICAL VELOCITY OMEGA
!              ------------------------------------

!*       1.1   LL_ETAD.

LL_ETAD = .FALSE.
IF(PRESENT(LDETADOT)) THEN
  LL_ETAD = LDETADOT
ENDIF

!*       1.2   COMPUTE OMEGA OR ETADOT ON FULL LEVELS

IF( LL_ETAD ) THEN
  IF(YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEV
      DO JL=KST,KPROF
        ZVVEL(JL,JLEV)=PEVEL(JL,JLEV)&
         & *(YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))*PRDELP(JL,JLEV)  
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,KFLEV
      DO JL=KST,KPROF
        ZVVEL(JL,JLEV)=0.5_JPRB*(PEVEL(JL,JLEV)+PEVEL(JL,JLEV-1))&
         & *(YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))*PRDELP(JL,JLEV)  
      ENDDO
    ENDDO
  ENDIF
ELSE
  DO JLEV=1,KFLEV
    DO JL=KST,KPROF
      ZVVEL(JL,JLEV)=PRXP(JL,JLEV,2)*PVVEL(JL,JLEV)
    ENDDO
  ENDDO
ENDIF

!*       1.3   PREPARE FOR EXTRAPOLATION ABOVE TOP.

ZVVEL(KST:KPROF,0)=0.0_JPRB

!*       1.4   INTERPOLATE.

ISLCT=2
CALL PPINTP(KPROMA,KST,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,KLEVB,ISLCT,&
 & LDBELO,LDBLOW,PRPRES,PRXP,PRXPD,ZVVEL,PVVPP)  

!*       1.5   CONSTANT EXTRAPOLATION UNDER LOWEST MODEL SURFACE

DO JLEVP=KLOLEV,KLEVP
  IF (LDBLOW(JLEVP)) THEN
    DO JL=KST,KPROF
      IF (LDBELO(JL,JLEVP)) THEN
        PVVPP(JL,JLEVP) = ZVVEL(JL,KFLEV)
      ENDIF
    ENDDO
  ENDIF
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPVVEL',1,ZHOOK_HANDLE)
END SUBROUTINE PPVVEL
