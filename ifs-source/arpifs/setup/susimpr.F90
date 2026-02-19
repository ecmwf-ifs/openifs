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

SUBROUTINE SUSIMPR(YDGEOMETRY,YDMDDH)

!**** *SUSIMPR*  - SET-UP FOR HORIZONTAL MEANS IF DDH IS NOT USED

!     Purpose.
!     --------
!       - ALLOCATION AND COMPUTATION OF THE WEIGHTS FOR EACH POINT
!         WARNING: PREPARED ONLY FOR ARPEGE !!!!
!         Computes HDSF = weight of each point. Currently HDSF is set to GAW
!         (Gaussian weights).

!**   Interface.
!     ----------
!        *CALL* *SUSIMPR

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about simplified physics

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        Documentation about simplified physics

!     Author.
!     -------
!        Original : 99-02-15 by M. Janiskova,
!         adopted the computation of the weight of each point in the case
!         when DDH is not called.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LALLOPR
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : MY_REGION_NS, MY_REGION_EW, NPRINTLEV
USE YOMMDDH  , ONLY : TMDDH

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TMDDH)    ,INTENT(INOUT) :: YDMDDH
INTEGER(KIND=JPIM) :: ILOEN(YDGEOMETRY%YRDIM%NDGSAG:YDGEOMETRY%YRDIM%NDGENG)

INTEGER(KIND=JPIM) :: IGL, IGL1, IGL2, IGLOFF, IGUX, IIGL,&
 & ILAT, ILON, ILONDEP, IU, JGL  

LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSIMPR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA,  &
 & YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGENL=>YDDIM%NDGENL, NDGSAG=>YDDIM%NDGSAG, &
 & NGPTOT=>YDGEM%NGPTOT, NSTAGP=>YDGEM%NSTAGP, &
 & NFRSTLAT=>YDMP%NFRSTLAT, NLSTLAT=>YDMP%NLSTLAT, NONL=>YDMP%NONL, &
 & NPTRFRSTLAT=>YDMP%NPTRFRSTLAT)
!     ------------------------------------------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

!     1. Allocation of hdsf

ALLOCATE(YDMDDH%HDSF(NGPTOT))
IF(LLP)WRITE(IU,9) 'HDSF     ',SIZE(YDMDDH%HDSF  ),SHAPE(YDMDDH%HDSF  )

!     2. Initialization of the dimension

IGUX=NDGENL
ILONDEP=1
IIGL=0
IGLOFF = NPTRFRSTLAT(MY_REGION_NS)
IGL1 = NFRSTLAT(MY_REGION_NS)
IGL2 = NLSTLAT(MY_REGION_NS)
DO JGL=IGL1,IGL2
  IGL=IGLOFF+JGL-IGL1
  IIGL=IIGL+1
  ILOEN(IIGL)=NONL(IGL,MY_REGION_EW)
ENDDO

!     3. Weigths for each grid point

YDMDDH%HDSF(:) = YDGSGEOM_NB%GAW(:)

ILAT = IGUX/2
ILON = ILOEN(ILAT)/2
IF (NPRINTLEV >= 1) THEN
  WRITE (NULOUT,*) ' '
  WRITE (NULOUT,*) ' JGL = ',ILAT,' JLON = ',ILON
  WRITE (NULOUT,*) ' RW = ',YDCSGLEG%RW(ILAT),' ILOEN = ',ILOEN(ILAT)
  WRITE (NULOUT,*) ' YDMDDH%HDSF(',ILON+NSTAGP(ILAT)-ILONDEP,')= ',&
   & YDMDDH%HDSF(ILON+NSTAGP(ILAT)-ILONDEP)  
ENDIF

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSIMPR',1,ZHOOK_HANDLE)
END SUBROUTINE SUSIMPR
