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

MODULE YOMGWDIAG
! Diagnose gravity wave noise by surface pressure tendency
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMDIM       , ONLY : TDIM
USE YOMCT3       , ONLY : NSTEP
USE YOMLUN       , ONLY : NULOUT,NULNAM
USE IOSTREAM_MIX , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,&
 &                        TYPE_IOSTREAM , TYPE_IOREQUEST, CLOSE_IOREQUEST
USE YOMHOOK      , ONLY : LHOOK    ,DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

REAL(KIND=JPRB),ALLOCATABLE :: SURF_PRESSSURE(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: GWDIAGS(:,:,:)
LOGICAL :: LGWDIAGS_ON
LOGICAL :: LGWDIAGS_FIELD_OUT
TYPE(TYPE_IOSTREAM) :: YIOSTREAM_GWD
CONTAINS
!======================================================================
SUBROUTINE SETUP_GWDIAG(YDDIM)

!**** *SUCT0*   - Routine to setup gravity wave diagnostics

!     Purpose.
!     --------
!           Initialize gravity wave diagnostics
!**   Interface.
!     ----------
!        *CALL* *SETUP_GWDIAG*

!        Explicit arguments :
!        --------------------
!        
!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 2013-01-01
!      ----------------------------------------------------------------
IMPLICIT NONE
TYPE(TDIM) , INTENT(IN) :: YDDIM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
NAMELIST /NAMGWDIAG/ LGWDIAGS_ON,LGWDIAGS_FIELD_OUT
#include "posnam.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOMGWDIAG:SETUP_GWDIAG',0,ZHOOK_HANDLE)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA)

LGWDIAGS_ON = .FALSE.
LGWDIAGS_FIELD_OUT = .FALSE.
CALL POSNAM(NULNAM,'NAMGWDIAG')
READ(NULNAM,NAMGWDIAG)
IF(LGWDIAGS_ON) THEN
  IF(ALLOCATED(SURF_PRESSSURE)) DEALLOCATE(SURF_PRESSSURE)
  IF(ALLOCATED(GWDIAGS)) DEALLOCATE(GWDIAGS)
  ALLOCATE(SURF_PRESSSURE(NPROMA,1,NGPBLKS))
  ALLOCATE(GWDIAGS(NPROMA,4,NGPBLKS))
  GWDIAGS(:,:,:) = 0.0
  IF(LGWDIAGS_FIELD_OUT) THEN
    CALL SETUP_IOSTREAM(YIOSTREAM_GWD,'CIO','GWDIAGS_FIELDS',CDMODE='w',KIOMASTER=1)
  ENDIF
ENDIF
WRITE(NULOUT,*) 'SETUP_GWDIAG: LGWDIAGS_ON = ',LGWDIAGS_ON,' LGWDIAGS_FIELD_OUT = ',LGWDIAGS_FIELD_OUT

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOMGWDIAG:SETUP_GWDIAG',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_GWDIAG
!======================================================================
SUBROUTINE UPDATE_GWDIAG(YDGEOMETRY,YDRIP,KPROMA,KFLEV,KPROF,KSTEP,KSTGLO,PAPRS,PDIVDP)

!**** *SUCT0*   - Routine to update gravity wave diagnostics

!     Purpose.
!     --------
!           Initialize gravity wave diagnostics
!**   Interface.
!     ----------
!        *CALL* *UPDATE_GWDIAG*

!        Explicit arguments :
!        --------------------
!        KPROMA               - HORIZONTAL DIMENSIONS.                 (INPUT)
!        KFLEV                - NUMBER OF MODEL LEVELS                 (INPUT)
!        KPROF                - DEPTH OF WORK                          (INPUT)
!        KSTEP                - TIME STEP                              (INPUT)
!        KSTGLO               - GLOBAL OFFSET                          (INPUT)
!        PAPRS                - HALF LEVEL PRESSURE                    (INPUT)
!        PDIVDP               - grad(vec(V) * (Delta prehyd))          (INPUT)
!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 2013-01-01
!    IMPLICIT NONE
USE YOMRIP       , ONLY : TRIP
USE GEOMETRY_MOD , ONLY : GEOMETRY
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVDP(KPROMA,KFLEV)

INTEGER(KIND=JPIM) :: IBLK,JLEV
REAL(KIND=JPRB) :: ZSDIV(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB) :: ZPSDIV(KPROMA,KFLEV+1)
REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "verdisint.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOMGWDIAG:UPDATE_GWDIAG',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA, &
 & YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA,YDCVER=>YDGEOMETRY%YRCVER, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(TSTEP=>YDRIP%TSTEP)
ZEPS = TINY(ZEPS)
IBLK = (KSTGLO-1)/KPROMA+1
  
IF(KSTEP > 0) THEN
  GWDIAGS(1:KPROF,1,IBLK) = (PAPRS(1:KPROF,KFLEV)- SURF_PRESSSURE(1:KPROF,1,IBLK))/TSTEP
ENDIF
SURF_PRESSSURE(1:KPROF,1,IBLK) = PAPRS(1:KPROF,KFLEV)
IF(YDCVER%LVERTFE) THEN

  DO JLEV=1,KFLEV
    ZSDIV(1:KPROF,JLEV)=PDIVDP(1:KPROF,JLEV)*YDVETA%VFE_RDETAH(JLEV)
  ENDDO
  ZSDIV(1:KPROF,0) = 0.0_JPRB
  ZSDIV(1:KPROF,KFLEV+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',KPROMA,1,KPROF,KFLEV,ZSDIV,ZPSDIV)
  GWDIAGS(1:KPROF,2,IBLK) = ZPSDIV(1:KPROF,KFLEV+1)

  ZSDIV(1:KPROF,:) = ABS(ZSDIV(1:KPROF,:))
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',KPROMA,1,KPROF,KFLEV,ZSDIV,ZPSDIV)
  GWDIAGS(1:KPROF,3,IBLK) = ZPSDIV(1:KPROF,KFLEV+1)
  
ELSE
  GWDIAGS(1:KPROF,2,IBLK) = 0.0_JPRB
  GWDIAGS(1:KPROF,3,IBLK) = 0.0_JPRB
  DO JLEV=KFLEV,1,-1
    GWDIAGS(1:KPROF,2,IBLK) = GWDIAGS(1:KPROF,2,IBLK)-PDIVDP(1:KPROF,JLEV)
    GWDIAGS(1:KPROF,3,IBLK) = GWDIAGS(1:KPROF,3,IBLK)+ABS(PDIVDP(1:KPROF,JLEV))
  ENDDO
ENDIF
GWDIAGS(1:KPROF,2,IBLK) = 10800.0*ABS(GWDIAGS(1:KPROF,2,IBLK))
GWDIAGS(1:KPROF,3,IBLK) = 10800.0*GWDIAGS(1:KPROF,3,IBLK)
GWDIAGS(1:KPROF,4,IBLK) = 100.0*GWDIAGS(1:KPROF,2,IBLK)/MAX(GWDIAGS(1:KPROF,3,IBLK),ZEPS)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOMGWDIAG:UPDATE_GWDIAG',1,ZHOOK_HANDLE)
END SUBROUTINE UPDATE_GWDIAG
!========================================================================
SUBROUTINE PRINT_GWDIAG(YDGEOMETRY,YDRIP)
!**** *SUCT0*   - Routine to output gravity wave diagnostics

!     Purpose.
!     --------
!           Output gravity wave diagnostics
!**   Interface.
!     ----------
!        *CALL* *PRINT_GWDIAG*

!        Explicit arguments :
!        --------------------
!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 2013-01-01
USE YOMRIP       , ONLY : TRIP
USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE
#include "gpnorm1.intfb.h"
TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
TYPE(TRIP)     ,INTENT(INOUT):: YDRIP
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST
REAL(KIND=JPRB) :: ZNORMS(3,4),ZTIME1,ZTIME2
REAL(KIND=JPRB) :: ZGWDIAGS(YDGEOMETRY%YRDIM%NPROMA,4,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMGWDIAG:PRINT_GWDIAG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & TSTEP=>YDRIP%TSTEP)
ZTIME1 = (NSTEP-2)*TSTEP/3600.0
ZTIME2 = (NSTEP-1)*TSTEP/3600.0

IF(NSTEP >= 1) THEN
  ZGWDIAGS(:,:,:) = GWDIAGS(:,:,:)
  ZGWDIAGS(:,1,:) =  10800.0*ABS(ZGWDIAGS(:,1,:))
  CALL GPNORM1(YDGEOMETRY,ZGWDIAGS,4,.TRUE.,PNORMS=ZNORMS)
ENDIF
IF(NSTEP >= 2) THEN
  WRITE(NULOUT,'(1X,A,1X,F6.2,3(1X,E11.5))') 'NORM GWDIAGS1',ZTIME1,ZNORMS(:,1)
ENDIF
IF(NSTEP >= 1) THEN
  WRITE(NULOUT,'(1X,A,1X,F6.2,3(1X,E11.5))') 'NORM GWDIAGS2',ZTIME2,ZNORMS(:,2)
  WRITE(NULOUT,'(1X,A,1X,F6.2,3(1X,E11.5))') 'NORM GWDIAGS3',ZTIME2,ZNORMS(:,3)
  WRITE(NULOUT,'(1X,A,1X,F6.2,3(1X,E11.5))') 'NORM GWDIAGS4',ZTIME2,ZNORMS(:,4)
ENDIF
IF(LGWDIAGS_FIELD_OUT) THEN
  CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,&
   & KGRIB2D=(/158/),KLEVS2D=(/1/),CDLEVTYPE='ML',KPROMA=NPROMA,PTSTEP=YDRIP%TSTEP)
  CALL IO_PUT(YIOSTREAM_GWD,YL_IOREQUEST,PR3=GWDIAGS(:,1:1,:))
  CALL CLOSE_IOREQUEST(YL_IOREQUEST)
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOMGWDIAG:PRINT_GWDIAG',1,ZHOOK_HANDLE)
END SUBROUTINE PRINT_GWDIAG

END MODULE YOMGWDIAG
