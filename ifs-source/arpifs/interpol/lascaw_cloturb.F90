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

SUBROUTINE LASCAW_CLOTURB(KFLEV,KPROM,KST,KPROF,PKHTURB,PCLO,PCLOSLD)

! ------------------------------------------------------------------
! Purpose:
!   Modify PCLO and PCLOSLD when 3D turbulence is active
!
! INPUT:
!   KFLEV    - Vertical dimension
!   KPROM    - horizontal dimension.
!   KST      - first element of arrays where computations are performed.
!   KPROF    - depth of work.
!   PKHTURB  - horizontal exchange coefficients for 3D turbulence
!
! INPUT/OUTPUT:
!   PCLO     - weights for horizontal cubic interpolations in longitude.
!   PCLOSLD  - cf. PCLO, SLHD case.
!
! Author:
!   H Petithomme (Dec 2020): after lascaw_cla (original from K Yessad, 2009)
! ------------------------------------------------------------------

USE PARKIND1,ONLY: JPIA,JPIM,JPRB
USE YOMHOOK, ONLY: LHOOK,DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
REAL(KIND=JPRB)   , INTENT(IN)  :: PKHTURB(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PCLO(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PCLOSLD(KPROM,KFLEV,3)

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB)    :: ZKH,ZWDS1
REAL(KIND=JPHOOK)  :: ZHOOK

IF (LHOOK) CALL DR_HOOK('LASCAW_CLOTURB',0,ZHOOK)

! apply the horizontal Laplacian to both PCLO and PCLOSLD:
DO JLEV=1,KFLEV
  DO JROF=KST,KPROF
    ZKH = 1._JPRB-2._JPRB*PKHTURB(JROF,JLEV)
    PCLO(JROF,JLEV,3)=PKHTURB(JROF,JLEV)*PCLO(JROF,JLEV,2)+PCLO(JROF,JLEV,3)
    ZWDS1=ZKH*PCLO(JROF,JLEV,1)+PKHTURB(JROF,JLEV)*PCLO(JROF,JLEV,2)
    PCLO(JROF,JLEV,2)=ZKH*PCLO(JROF,JLEV,2)+PKHTURB(JROF,JLEV)*PCLO(JROF,JLEV,1)
    PCLO(JROF,JLEV,1)=ZWDS1

    PCLOSLD(JROF,JLEV,3)=PKHTURB(JROF,JLEV)*PCLOSLD(JROF,JLEV,2)+PCLOSLD(JROF,JLEV,3)
    ZWDS1=ZKH*PCLOSLD(JROF,JLEV,1)+PKHTURB(JROF,JLEV)*PCLOSLD(JROF,JLEV,2)
    PCLOSLD(JROF,JLEV,2)=ZKH*PCLOSLD(JROF,JLEV,2)+PKHTURB(JROF,JLEV)*PCLOSLD(JROF,JLEV,1)
    PCLOSLD(JROF,JLEV,1)=ZWDS1
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('LASCAW_CLOTURB',1,ZHOOK)
END SUBROUTINE LASCAW_CLOTURB
