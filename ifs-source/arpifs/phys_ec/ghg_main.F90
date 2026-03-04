! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE GHG_MAIN( YDML_CHEM,YGFL,KIDIA, KFDIA, KLON, KLEV, KFLDX, KLEVX, KTRAC, PTSTEP, PRS1, &  
    & PLRCH4, PCFLX, PNEEFLX, PGELAT, PGFL, PTENGFL, PEXTRA, PCEN, PTENC ) 

!     PURPOSE.
!     --------
!     Apply surface flux sources/sinks for GHG 
!       GHG CH4 Tendencies 
!       GHG CO2 aircraft emissions
!
!**   INTERFACE.
!     ----------
!          *GHG_MAIN *CHEM_MAIN_LAYER .


! INPUTS:
!  -------

! INPUTS/OUTPUTS:

!-----------------------------------------------------------------------

!     Externals.
!     ---------


!     Author
!    --------
!         2016-09-12, J. Flemming 

!    
!-----------------------------------------------------------------------

USE MODEL_CHEM_MOD , ONLY : MODEL_CHEM_TYPE
USE PARKIND1       , ONLY : JPIM, JPRB
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCHEM        , ONLY : IEXTR_CH,IEXTR_EM
USE YOM_YGFL       , ONLY : TYPE_GFLD
USE YOM_GRIB_CODES , ONLY : NGRBGHG

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_CHEM_TYPE),INTENT(INOUT):: YDML_CHEM
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(KLON,KLEV,YGFL%NDIM), PTENGFL(KLON,KLEV,YGFL%NDIM1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLRCH4(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFLX(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNEEFLX(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KLON)

!-----------------------------------------------------------------------
! local variables  
INTEGER(KIND=JPIM) :: JL, JK, JT, IPCO2, IPCH4
INTEGER(KIND=JPIM), PARAMETER :: ITMODE=1  ! use MMR for loss tendency before emissions update 
!INTEGER(KIND=JPIM), PARAMETER :: ITMODE=2 ! use MMR                    after emissions update 
REAL(KIND=JPRB)    :: ZDELP(KLON,KLEV) 
REAL(KIND=JPRB)    :: ZTENC1(KLON,KLEV),  ZTENC2(KLON,KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "chem_inext.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GHG_MAIN',0,ZHOOK_HANDLE)
ASSOCIATE(LCHEM_DIA=>YDML_CHEM%YRCOMPO%LCHEM_DIA, YGHG=>YGFL%YGHG, & 
 & NCHEM=>YGFL%NCHEM, NACTAERO=>YGFL%NACTAERO, YLRCH4=>YGFL%YLRCH4, &
 & NGHG=>YGFL%NGHG)
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
     ZDELP(JL,JK) = PRS1(JL,JK)- PRS1(JL,JK-1)
   ENDDO
ENDDO

IPCO2=-999_JPIM
IPCH4=-999_JPIM

DO JT=1,NGHG
   IF (YGHG(JT)%IGRBCODE == NGRBGHG(1) ) IPCO2=JT
   IF (YGHG(JT)%IGRBCODE == NGRBGHG(2) ) IPCH4=JT
ENDDO

! Apply CH4 tendencies 
! GHG always first position in PTENC  - see gems_init.f90
! CH4 
IF(IPCH4 >= 0 ) THEN
  ZTENC1(KIDIA:KFDIA,1:KLEV) = PTENC(KIDIA:KFDIA,1:KLEV,IPCH4) 
  IF (ITMODE == 1 ) THEN 
    PTENC(KIDIA:KFDIA,1:KLEV,IPCH4) = PTENC(KIDIA:KFDIA,1:KLEV,IPCH4) - &
   &     PLRCH4(KIDIA:KFDIA,1:KLEV)*PGFL(KIDIA:KFDIA,1:KLEV,YGHG(IPCH4)%MP9_PH)
  ENDIF  
  IF (ITMODE == 2 ) THEN 
     PTENC(KIDIA:KFDIA,1:KLEV,IPCH4) = PTENC(KIDIA:KFDIA,1:KLEV,IPCH4) - &
  &     PLRCH4(KIDIA:KFDIA,1:KLEV)*PCEN(KIDIA:KFDIA,1:KLEV,IPCH4)
  ENDIF 
!  IF (ANY ( PGFL(KIDIA:KFDIA,1:KLEV,YLRCH4%MP9_PH)- PLRCH4(KIDIA:KFDIA,1:KLEV)  .NE. 0.0_JPRB )) THEN 
!   print*, ' CH4 tendencies wrong ? ' 
!  ENDIF 

! budget CH4 tendencies 
  IF (LCHEM_DIA) THEN
  ZTENC2(KIDIA:KFDIA,1:KLEV) = PTENC(KIDIA:KFDIA,1:KLEV,IPCH4) 
   ! CH4 chemistry tendency in diagnostics arrays 
    CALL CHEM_INEXT(KIDIA, KFDIA, KLON, KLEV,  1, 1,  ZDELP, &
   & PTSTEP, ZTENC2, ZTENC1, PEXTRA(:,NCHEM+NACTAERO+IPCH4, IEXTR_CH))
  ENDIF
ENDIF


IF(IPCO2 >= 0 ) THEN
  IF (LCHEM_DIA) THEN
! store nee flux as chemical tendency (Note that the sign of NEE has to be changed to follow atmospheric budget convention)
    PEXTRA(KIDIA:KFDIA,NCHEM+NACTAERO+IPCO2, IEXTR_CH) = PEXTRA(KIDIA:KFDIA,NCHEM+NACTAERO+IPCO2, IEXTR_CH) - & 
 &        PNEEFLX(KIDIA:KFDIA)  * PTSTEP
  ENDIF
ENDIF


! budget CO2 and CH4 surface emissions (other than NEE)
IF (LCHEM_DIA) THEN
  DO JT=1,NGHG
    DO JL=KIDIA,KFDIA
     ! GHG always in position 1 and 2 
      PEXTRA(JL,JT+NCHEM+NACTAERO, IEXTR_EM) = PEXTRA(JL,JT+NCHEM+NACTAERO, IEXTR_EM )- PCFLX(JL,JT) * PTSTEP
    ENDDO
  ENDDO
ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GHG_MAIN',1,ZHOOK_HANDLE)
END SUBROUTINE GHG_MAIN
