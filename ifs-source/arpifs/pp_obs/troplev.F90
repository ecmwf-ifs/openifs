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

SUBROUTINE TROPLEV (KPROMA,KSTART,KPROF,KFLEV,LDTROPHUM,PT,PQ,PRESF,KLTROP)

!      TROPLEV - Calculate model level closest to the tropopause.
!                Original code by Adrian Simmons.

!     Arguments :
!     -----------
!     Input:
!       KPROMA              - horizontal dimension
!       KSTART              - start of work
!       KPROF               - depth of work
!       KFLEV               - number of model levels
!       LDTROPHUM            - if true, calculate humidity based tropopause
!       PT(KPROMA,KFLEV)    - temperature on full model levels
!       PQ(KPROMA,KFLEV)    - humidity on full model levels
!       PRESF(KPROMA,KFLEV) - model full level pressures

!     Input/Output:
!       KLTROP(KPROMA)      - tropopause model levels

!     Author:  Adrian Simmons   *ECMWF*
!     -------

!     Method:  last verification and derivation by A. Simmons & P. Bechtold 26/03/2021
!              use Theta coordinates, but do not start comput from bottom and use test on given T gradient (WMO Def),
!              and use given buoyancy frequency N criterion and start computation from stratosphere
!              d ln theta = d ln T - K d ln P =>  P/theta dtheta/dp = P/T dT/dp - K
!              with 1/theta dtheta/dp = 1/theta dtheta/dz dz/dp = R T/P 1/g N^2

!     Modifications:
!     --------------
!     E.Holm       05-10-20   Moved code here from symtransin and jgnr/i/ad.
!     E.Holm       07-12-14   Alternative humidity tropopause definition
!     E.Holm       12-11-01   Move logical for humidity tropopause to call
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RD, RG, RKAPPA
IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
LOGICAL           ,INTENT(IN)    :: LDTROPHUM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESF(KPROMA,KFLEV)
INTEGER(KIND=JPIM),INTENT(OUT) :: KLTROP(KPROMA)

INTEGER(KIND=JPIM) :: ILEVM2, JROF, JLEV
REAL(KIND=JPRB) :: ZSTAB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TROPLEV',0,ZHOOK_HANDLE)

!*         1. COMPUTE TROPOPAUSE HEIGHT FROM LAPSE RATE DEFINITION
!*            OR HUMIDITY DEFINITION
!             ----------------------------------------------------

ILEVM2=KFLEV-2
DO JROF=KSTART,KPROF
  KLTROP(JROF)=ILEVM2
ENDDO

DO JLEV=1,ILEVM2
  DO JROF=KSTART,KPROF
    IF (PRESF(JROF,JLEV)>=7.E3_JPRB .AND. PRESF(JROF,JLEV)<=5.E4_JPRB) THEN
      IF (LDTROPHUM) THEN
!*         1.1 Humidity tropopause
        IF (KLTROP(JROF)==ILEVM2 .AND. PQ(JROF,JLEV)  >3.E-6_JPRB &
           &                     .AND. PQ(JROF,JLEV+2)>5.E-6_JPRB) THEN
          KLTROP(JROF)=JLEV
        ENDIF
      ELSE
!*         1.2 Lapse rate tropopause
         ZSTAB =  PRESF(JROF,JLEV+1)*(   PT(JROF,JLEV+2)-   PT(JROF,JLEV))&
               &  /(PT(JROF,JLEV+1)*(PRESF(JROF,JLEV+2)-PRESF(JROF,JLEV)))&
               &+ RD*PT(JROF,JLEV+1)*2.5E-4_JPRB/(RG*RG)
         IF (KLTROP(JROF)==ILEVM2 .AND. ZSTAB>RKAPPA) THEN
           KLTROP(JROF)=JLEV
         ENDIF
      ENDIF
    ELSEIF (PRESF(JROF,JLEV)>5.E4_JPRB .AND. KLTROP(JROF)==ILEVM2)THEN
      KLTROP(JROF)=JLEV-1
    ENDIF
  ENDDO
ENDDO


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TROPLEV',1,ZHOOK_HANDLE)
END SUBROUTINE TROPLEV

