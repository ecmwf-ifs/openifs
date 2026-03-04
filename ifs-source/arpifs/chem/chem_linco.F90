! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_LINCO &
 &    (YDCHEM,YGFL,KIDIA  , KFDIA , KLON , KLEV , KGPLAT, &
 &     PTSTEP , PTP, PCEN ,&
 &     PTENC) 

!**   DESCRIPTION 
!     ----------
!
!   Linear CO chemical scheme routine for IFS chemistry 
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_LINCO* IS CALLED FROM *CHEM_MAIN*.

! INPUTS:
! -------
!
! - DIMENSIONS ETC.
!
! KIDIA   :  Start of Array  
! KFDIA   :  End  of Array 
! KLON    :  Length of Arrays 
! KLEV    :  Number of Levels
! KGPLAT  : DM-global number of the latitude of point jrof=KSTART
! PTSTEP  :  Time step in seconds 
!
! - 2D and 3D
!
! PTP(KLON,KLEV)                : FULL-LEVEL TEMPERATURE (W. DYN.TEND.) (K)
! PCEN(KLON,KLEV,1)             : CONCENTRATION OF TRACERS           (kg/kg)
!
! OUTPUTS:
! -------
!
! PTENC  (KLON,KLEV,1)          : TENDENCY OF CONCENTRATION OF TRACERS including chemistry(kg/kg s-1)
!
! LOCAL:
! -------
!
! ZKCO(KLON,KLEV,NCHEM_LCOCOEF) : COEFICIENTS OF THE LINEAR SCHEME (volume mixing ratio)
! ZCEN(KLON,KLEV,1)             : CONCENTRATION OF TRACERS           (kg/kg)
!
!     AUTHORS.
!     -------
!        JOHANNES FLEMMING  *ECMWF*
!        SEBASTIEN MASSART  *CERFACS/ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2009-07-22
!        LINEAR SCHEME : 2011-03-10 (SM)




USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL  ,ONLY : TYPE_GFLD
USE YOMCST    ,ONLY : RMD
USE YOMCHEM  , ONLY : TCHEM


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TCHEM)       ,INTENT(INOUT):: YDCHEM
TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KGPLAT(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV)   
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV,1) 
REAL(KIND=JPRB),INTENT(OUT)   :: PTENC(KLON,KLEV,1)

! * LOCAL 

!
REAL(KIND=JPRB),ALLOCATABLE    :: ZKCO(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE    :: ZTOP(:)
REAL(KIND=JPRB),ALLOCATABLE    :: ZCEN(:,:,:)
!
INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: JLEV
INTEGER(KIND=JPIM) :: JFLEV
!INTEGER(KIND=JPIM) :: ICHEM
!
REAL(KIND=JPRB) :: ZANEX
REAL(KIND=JPRB) :: ZRAPP
REAL(KIND=JPRB) :: ZEXPCO
!REAL(KIND=JPRB) :: ZTAUTOPCO
!REAL(KIND=JPRB) :: ZCLIMTOPCO
!
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "gpincoch.intfb.h"


IF (LHOOK) CALL DR_HOOK('CHEM_LINCO',0,ZHOOK_HANDLE)
ASSOCIATE(YCHEM=>YGFL%YCHEM,&
 & NCHEM_LCOCOEF=>YDCHEM%NCHEM_LCOCOEF, LCHEM_LCOMESO=>YDCHEM%LCHEM_LCOMESO, &
 & LCHEM_LCOCSTCLIM=>YDCHEM%LCHEM_LCOCSTCLIM, RCHEM_LCOTAUTOP=>YDCHEM%RCHEM_LCOTAUTOP, &
 & RCHEM_LCOCLIMTOP=>YDCHEM%RCHEM_LCOCLIMTOP, LCHEM_LCOLIMIT=>YDCHEM%LCHEM_LCOLIMIT, &
 & RCHEM_LCOCOEFA1=>YDCHEM%RCHEM_LCOCOEFA1)

!
! Allocate local arrays
!
ALLOCATE(ZKCO(KLON,KLEV,NCHEM_LCOCOEF))
ALLOCATE(ZTOP(KLON))
ALLOCATE(ZCEN(KLON,KLEV,1))

!
! Compute tendency
!
!
!     Update the coefficients
!
CALL GPINCOCH(YDCHEM,KIDIA,KFDIA,KLEV,KLON,KGPLAT,ZKCO,ZTOP)
!
! Set the number of level to apply the scheme
!
IF (LCHEM_LCOMESO) THEN
  JFLEV=2
ELSE
  JFLEV=1
ENDIF
!
! Loop over the levels and the grid points
!
DO JLEV=JFLEV,KLEV
  DO JL=KIDIA,KFDIA
      !
      ! Compute tendancy
      !
      ZANEX = ( ZKCO(JL,JLEV,2) * RCHEM_LCOCOEFA1 +&
 &              ZKCO(JL,JLEV,5)* (PTP(JL,JLEV)-ZKCO(JL,JLEV,4)) -&
 &              ZKCO(JL,JLEV,3)*ZKCO(JL,JLEV,1) )* PTSTEP
      ZRAPP  = 1._JPRB - PTSTEP * ZKCO(JL,JLEV,3)
      IF (LCHEM_LCOLIMIT) ZRAPP  = MAX(ZRAPP, 1.E-06_JPRB)
      ZANEX = ZANEX/RMD*YCHEM(1)%RMOLMASS
      ZANEX =  (ZANEX+PCEN(JL,JLEV,1)) / ZRAPP
      IF (LCHEM_LCOLIMIT) THEN
        ZANEX = ZANEX-PCEN(JL,JLEV,1)
        ZCEN(JL,JLEV,1)=MAX(1.E-30_JPRB, PCEN(JL,JLEV,1)+ZANEX)
        PTENC(JL,JLEV,1) = (ZCEN(JL,JLEV,1) - PCEN(JL,JLEV,1)) / PTSTEP
      ELSE 
        PTENC(JL,JLEV,1) = (ZANEX-PCEN(JL,JLEV,1)) / PTSTEP
      ENDIF  
      !
  ENDDO
ENDDO
!
! Relaxation toward constant top boundary climatology
! if we want a mesopsheric CO
!
IF (LCHEM_LCOMESO) THEN 
    IF (LCHEM_LCOCSTCLIM) THEN 
      ZEXPCO  = EXP(-PTSTEP/(RCHEM_LCOTAUTOP*86400._JPRB))
      !
      DO JL=KIDIA,KFDIA
        PTENC(JL,1,1) = (1._JPRB - ZEXPCO)&
 &                    * (RCHEM_LCOCLIMTOP - PCEN(JL,1,1)) / PTSTEP
      ENDDO
    !
    ! Relaxation toward time-dependent top boundary climatology
    !
    ELSE  ! LCHEM_LCOCSTCLIM
      ZEXPCO  = EXP(-PTSTEP/(RCHEM_LCOTAUTOP*86400._JPRB))
      !
      DO JL=KIDIA,KFDIA
        PTENC(JL,1,1) = (1._JPRB - ZEXPCO)&
 &                    * (ZTOP(JL) - PCEN(JL,1,1)) / PTSTEP
      ENDDO
    ENDIF   ! LCHEM_LCOCSTCLIM
ENDIF ! LCHEM_LCOMESO
!

!
! Deallocate local arrays
!
DEALLOCATE(ZKCO)
DEALLOCATE(ZTOP)
DEALLOCATE(ZCEN)


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_LINCO',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_LINCO
