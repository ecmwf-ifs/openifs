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

SUBROUTINE GPPECHOT(KPROMA,KSTART,KPROF,KFLEV,PPRESF,PSIMRFCDB,PECHOT)

USE PARKIND1,ONLY : JPRD, JPIM,  JPRB
USE YOMHOOK ,ONLY : LHOOK, DR_HOOK, JPHOOK


!**** *GPPRECHOT*  - Compute Simulated Echo Top -

!     PURPOSE.
!     --------
!        Compute summit in (hPa) of reflectivities superior to RECHOT 
!**   INTERFACE.
!     ----------
!       *CALL* *GPPECHOT*

!        EXPLICIT ARGUMENTS
!        --------------------
!            INPUT :
!        KPROMA    : Horizontal dimension
!        KSTART    : start of work
!        KPROF     : depth of work
!        KFLEV     : number of vertical levels
!        PPRESF    : pressure (KRPOMA,KFLEV)
!        PSIMRFCDB  : simulated reflectivies in dBZ (KRPOMA,KFLEV)
!            OUTPUT:
!        PECHOT    : pressure of echotop (KRPOMA)

!        IMPLICIT ARGUMENTS
!        --------------------
!           NONE

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
   
!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        Olivier Jaron - Nicolas Merlet *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 18-12-06
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA,KSTART,KPROF,KFLEV
REAL(KIND=JPRB),INTENT(IN)  :: PPRESF(KPROMA,KFLEV),PSIMRFCDB(KPROMA,KFLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PECHOT(KPROMA)

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JROF,JLEV
INTEGER(KIND=JPIM) :: IECHOTOP(KPROMA)
REAL(KIND=JPRB)    :: ZALPHA,ZBETA,ZMISVAL
REAL(KIND=JPRB) :: RECHOT

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPECHOT',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!* 1. Initialization  

ZMISVAL = 0._JPRB
RECHOT = 18._JPRB

!* 2. Compute reflectivities on each model vertical level (0D operator)*
!     Search upwards the last (ie uppest) level where PSIMRFCDB > RECHOT

IECHOTOP(KSTART:KPROF) = KFLEV
DO JLEV=KFLEV,2,-1
  DO JROF=KSTART,KPROF
    IECHOTOP(JROF) = KFLEV-MAX((KFLEV-JLEV)*INT(SIGN(1._JPRB,PSIMRFCDB(JROF,JLEV)-RECHOT)),KFLEV-IECHOTOP(JROF))
  ENDDO
ENDDO

!* 3. Interpolation and reset

DO JROF=KSTART,KPROF
   ! Interpolation along log(Pressure) and convert into hPa
   ZALPHA = ABS(PSIMRFCDB(JROF,IECHOTOP(JROF))-RECHOT)
   ZBETA  = ABS(RECHOT-PSIMRFCDB(JROF,IECHOTOP(JROF)-1_JPIM))
   PECHOT(JROF)= EXP((ZBETA*LOG(PPRESF(JROF,IECHOTOP(JROF)))+ZALPHA*LOG(PPRESF(JROF,IECHOTOP(JROF)-1_JPIM)))  &
           & / MAX(ZALPHA+ZBETA,1.E-08_JPRB))  * 0.01_JPRB
   ! Reset to zero where there is no PSIMRFCDB > RECHOT
   PECHOT(JROF)=PECHOT(JROF) * FLOAT(MAX(SIGN(1_JPIM,KFLEV-IECHOTOP(JROF)-1_JPIM),0_JPIM))&
    &         + ZMISVAL      * FLOAT(MAX(SIGN(1_JPIM,IECHOTOP(JROF)-KFLEV),0_JPIM))
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPPECHOT',1,ZHOOK_HANDLE)
END SUBROUTINE  GPPECHOT
