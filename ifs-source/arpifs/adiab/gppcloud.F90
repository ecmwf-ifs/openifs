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

SUBROUTINE GPPCLOUD(KPROMA,KSTART,KPROF,KFLEV,PXCLDFRAC,PPRESF,PCLDFC,PGEOF,POROG,PPTOP,PPBASE,PHBASE)

USE PARKIND1, ONLY : JPIM,  JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST  , ONLY : RG
USE YOMLUN  , ONLY : NULOUT


!**** *GPPCLOUD*  - Compute Top and Base of clouds -

!     PURPOSE.
!     --------

!     Return the pressure (Pa) of base and top of clouds according to the threshold on cloudness fraction
!     and the height (m) above ground of base.

!**   INTERFACE.
!     ----------
!       *CALL* *GPPCLOUD*

!        EXPLICIT ARGUMENTS
!        --------------------
!            INPUT :
!        KPROMA    : Horizontal dimension
!        KSTART    : start of work
!        KPROF     : depth of work
!        KFLEV     : number of vertical levels
!        PXCLDFRAC : Threshlod for cloud fraction
!        PPRESF    : Pressure (KRPOMA,KFLEV)
!        PCLDFC    : Cloud Fraction
!        PGEOF     : Geopotentiel 
!        POROG     : Orography geopotentiel
!            OUTPUT:
!        PPTOP      : pressure of top of clouds (KRPOMA)
!        PPBASE     : pressure of base of clouds (KRPOMA)
!        PHBASE     : height   of base of clouds (KRPOMA)

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
!        Olivier Jaron  *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : june 2020
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA,KSTART,KPROF,KFLEV
REAL(KIND=JPRB),INTENT(IN)    :: PXCLDFRAC
REAL(KIND=JPRB),INTENT(IN)    :: PPRESF(KPROMA,KFLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PCLDFC(KPROMA,KFLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PGEOF(KPROMA,KFLEV)
REAL(KIND=JPRB),INTENT(IN)    :: POROG(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PPTOP(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PPBASE(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PHBASE(KPROMA)

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JROF,JLEV
INTEGER(KIND=JPIM) :: ITOP(KPROMA), IBASE(KPROMA), IDXCLD(KPROMA)
REAL(KIND=JPRB)    :: ZALPHA,ZBETA
REAL(KIND=JPRB)    :: RMISSVALH, RMISSVALP, Z1SRG ,RMINPRES

LOGICAL ::   LPTOP     ! indicate if top of cloud is required
LOGICAL ::   LPBASE    ! indicate if base (in pressure) of cloud is required
LOGICAL ::   LHBASE    ! indicate if base (in height)  of cloud is required

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPCLOUD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!* 1. Parameters  

RMISSVALH=20000._JPRB  ! Miss value for height coordinate
RMISSVALP=0._JPRB      ! Miss value for pressure coordinate
Z1SRG=1.0_JPRB/RG      ! 1/g
RMINPRES = 1000._JPRB  ! Minimum of pressure to detect cloud (avoid stratospheric clouds)

!* 2. Intialisation of comptutanional keys

LPTOP=PRESENT(PPTOP)
LPBASE=PRESENT(PPBASE)
LHBASE=PRESENT(PHBASE)


!* 2. Compute full level top and bottom of clouds

ITOP (KSTART:KPROF) =KFLEV
IBASE(KSTART:KPROF) =1_JPIM
DO JLEV=KFLEV,2,-1
  DO JROF=KSTART,KPROF
    IDXCLD(KPROMA)=INT(SIGN(1._JPRB,PCLDFC(JROF,JLEV)-PXCLDFRAC))&
     & * MAX(0_JPIM, INT(SIGN(1._JPRB,PPRESF(JROF,JLEV)-RMINPRES)))
    
    ITOP(JROF)=KFLEV-MAX(IDXCLD(KPROMA)*(KFLEV-JLEV),KFLEV-ITOP(JROF))
    IBASE(JROF)=     MAX(IDXCLD(KPROMA)*       JLEV ,     IBASE(JROF))
  ENDDO
ENDDO

!* 3. Interpolation verticale lineaire suivant log(P) en fonction de CLDFC
IF (LPTOP) THEN
  DO JROF=KSTART,KPROF
    IF (ITOP(JROF) < KFLEV) THEN
      ZALPHA = PCLDFC(JROF,ITOP(JROF))-PXCLDFRAC
      ZBETA  = PXCLDFRAC-PCLDFC(JROF,ITOP(JROF)-1_JPIM)
      PPTOP(JROF)=EXP((ZBETA*LOG(PPRESF(JROF,ITOP(JROF))) &
             & + ZALPHA*LOG(PPRESF(JROF,ITOP(JROF)-1_JPIM)))  &
             & / MAX(ZALPHA+ZBETA,1.E-12_JPRB))
    ELSE
      PPTOP(JROF)=RMISSVALP
    ENDIF  
  ENDDO
ENDIF  

IF (LPBASE) THEN
  DO JROF=KSTART,KPROF
    IF ((IBASE(JROF) > 1_JPIM) .AND. (IBASE(JROF) < KFLEV)) THEN
     ZALPHA = PXCLDFRAC-PCLDFC(JROF,IBASE(JROF)+1_JPIM)
     ZBETA  = PCLDFC(JROF,IBASE(JROF))-PXCLDFRAC
     PPBASE(JROF)=EXP((ZBETA*LOG(PPRESF(JROF,IBASE(JROF)+1_JPIM)) &
             & + ZALPHA*LOG(PPRESF(JROF,IBASE(JROF))))  &
             & / MAX(ZALPHA+ZBETA,1.E-12_JPRB))
    ELSEIF (IBASE(JROF) == KFLEV) THEN 
      PPBASE(JROF)=PPRESF(JROF,KFLEV)
    ELSE   
      PPBASE(JROF)=RMISSVALP
    ENDIF  
  ENDDO
ENDIF

!* 4. Interpolation verticale lineaire suivant la hauteur du geopotentiel en fonction de CLDFC

IF (LHBASE) THEN
  DO JROF=KSTART,KPROF
    IF ((IBASE(JROF) > 1_JPIM) .AND. (IBASE(JROF) < KFLEV)) THEN
      ZALPHA = PXCLDFRAC-PCLDFC(JROF,IBASE(JROF)+1_JPIM)
      ZBETA  = PCLDFC(JROF,IBASE(JROF))-PXCLDFRAC
      PHBASE(JROF)=(((ZBETA*PGEOF(JROF,IBASE(JROF)+1_JPIM) &
              & + ZALPHA*PGEOF(JROF,IBASE(JROF)))  &
              & / MAX(ZALPHA+ZBETA,1.E-12_JPRB)-POROG(JROF))) * Z1SRG
    ELSEIF (IBASE(JROF) == KFLEV) THEN
      PHBASE(JROF)=0._JPRB
    ELSE  
      PHBASE(JROF)=RMISSVALH
    ENDIF 
    IF (PHBASE(JROF) < 0._JPRB) PHBASE(JROF)=0._JPRB 
  ENDDO
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPPCLOUD',1,ZHOOK_HANDLE)
END SUBROUTINE  GPPCLOUD
