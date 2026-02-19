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

SUBROUTINE PPQ(KPROMA,KST,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,&
 & KLEVB,PRPRES,LDBELO,LDBLOW,&
 & PRXP,PRXPD,PQF,PQPP,LDBOB,PBOB0HACK)

!**** *PPQ* - POST-PROCESS VARIABLES LIKE SPECIFIC HUMIDITY.

!     PURPOSE.
!     --------
!           PERFORMS THE VERTICAL INTERPOLATION OF SPECIFIC HUMIDITY GIVEN
!       MODEL CO-ORDINATE LEVEL.

!**   INTERFACE.
!     ----------
!        *CALL* *PPQ*

!        EXPLICIT ARGUMENTS
!        --------------------

!        * INPUT:
!         KPROMA    - horizontal dimension.
!         KST       - start of work.
!         KPROF     - depth of work.
!         KFLEV     - number of input pressure levels.
!         KLEVP     - number of output pressure levels.
!         KLOLEV    - beginning for the interpolation.
!        KPPM     - Number of interpolation methods in post-processing (INPUT-C)

!         KLEVB     - input level below PRPRES (see routine PPFLEV).
!         PRPRES    - post-processing level pressures.
!         LDBELO    - .TRUE. if pressure is under lowest model level.
!         LDBLOW    - .TRUE. if LOBELO(J) is containing at least one .TRUE.
!         PRXP      - contains "prehyd" and ln(prehyd) at half and full levels.
!                     (see routine PPINIT).
!         PRXPD     - contains 1/[Delta prehyd] and  1/[Delta log(prehyd)].
!         PQF       - full level input quantity to be interpolated.

!        * OUTPUT:
!         PQPP      - post-processed quantity.

!        * OPTIONAL INPUT:
!         LDBOB     - Simulates former routine BOB if T.

!        IMPLICIT ARGUMENTS :  NONE.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.  PPINTP - LINEAR INTERPOLATION
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-01-26
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K.Yessad (Apr 2009): merge BOB in PPQ.
!        A.Geer        24-Jul-2015 Pre-OOPS: 0:NFLEVG deprecated
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQPP(KPROMA,KLEVP) 
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LDBOB
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PBOB0HACK

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISLCT, JL, JLEVP
LOGICAL :: LLBOB
REAL(KIND=JPRB) :: ZS(KPROMA) 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZQF(KPROMA,0:KFLEV)

!     ------------------------------------------------------------------

#include "ppintp.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PPQ',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    VERTICAL INTERPOLATION.
!              -----------------------

ZQF(KST:KPROF,1:KFLEV) = PQF(KST:KPROF,1:KFLEV)

!*       1.1   PRELIMINARY CALCULATIONS.

IF (PRESENT(LDBOB)) THEN
  LLBOB=LDBOB
ELSE
  LLBOB=.FALSE.
ENDIF

IF (LLBOB) THEN
  ZQF(KST:KPROF,0)=PBOB0HACK
  ZS(KST:KPROF)=1.0_JPRB
  ISLCT=1
ELSE
  ZQF(KST:KPROF,0)=ZQF(KST:KPROF,1)
  ZS(KST:KPROF)=ZQF(KST:KPROF,KFLEV)
  ISLCT=2
ENDIF

!*       1.2   INTERPOLATE SPECIFIC HUMIDITY.

CALL PPINTP(KPROMA,KST,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,&
 & KLEVB,ISLCT,&
 & LDBELO,LDBLOW,PRPRES,PRXP,PRXPD,&
 & ZQF,PQPP)  

!*       1.3   EXTRAPOLATE  BELOW SURFACE.

DO JLEVP=KLOLEV,KLEVP
  IF (LDBLOW(JLEVP)) THEN
    DO JL=KST,KPROF
      IF(LDBELO(JL,JLEVP)) THEN
        PQPP(JL,JLEVP)=ZS(JL)
      ENDIF
    ENDDO
  ENDIF
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPQ',1,ZHOOK_HANDLE)
END SUBROUTINE PPQ
