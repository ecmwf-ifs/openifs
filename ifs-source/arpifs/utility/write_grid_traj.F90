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

SUBROUTINE WRITE_GRID_TRAJ(YDGEOMETRY,YDRIP,CDFILE,KGRIB_HD,KTYPE3D,KTYPE2D,K3D,K2D,&
                        & PGRID,PSPSP,KTYPE_INTERP,LDNEW,LDCLOSE,KGRIB_HD2)


!     Purpose.
!     --------
!     Gather grid point fields and write them out to GRIB file.

!     Arguments (all input)
!     ---------------------
!       CDFILE       : Filename
!       KGRIB_HD     : GRIB header for output
!       KGRIB_HD     : GRIB header for output additional for grib2 sfc
!       KTYPE3D      : List of 3D field variable GRIB definitions
!       KTYPE2D      : List of 2D field variable GRIB definitions
!       K3D          : Number of 3D field variables
!       K2D          : Number of 2D field variables
!       PGRID        : Grid point fields to write
!       PSPSP        : LnPs spectral field (required only if KTYPE_INTERP=3)
!       KTYPE_INTERP : Interpolation type
!                      1 - bilinear
!                      2 - bicubic
!                      3 - conserving bilinear/bicubic
!       LDNEW        : true : write to new file
!                      false: append to existing file
!       LDCLOSE      : close files

!     Author.
!     -------
!       L. Isaksen and Y. Tremolet

!     Modifications.
!     --------------
!       Original : 09-Jun-2004
!        L.Isaksen     18-Jan-2005 Introduce conserving interpolation
!        Y.Tremolet    18-Feb-2005 Split pre_grid_biconserv
!        Y.Tremolet    30-Mar-2005 PREGRBENC -> PRESET_GRIB_HEAD
!        E.Holm        05-Feb-2008 Allocate surface pressure for biconserv
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
! ----------------------------------------------------------------------

USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE PARDIM             , ONLY : JPMXLE
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMMP0             , ONLY : NPROC, MYPROC
USE IOSTREAM_MIX       , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,&
 &                              CLOSE_IOSTREAM, TYPE_IOSTREAM, TYPE_IOREQUEST,&
 &                              CLOSE_IOREQUEST, INQUIRE_GRIB
USE YOM_GRID_BICONSERV , ONLY : RGPPRS_HR, RGPPRS_LR
USE YOM_GRIB_CODES     , ONLY : NGRBSD,NGRBWSN, NGRBRSN, NGRBTSN, NGRBASN

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDFILE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIB_HD
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN) :: KGRIB_HD2
INTEGER(KIND=JPIM),INTENT(IN)    :: K3D, K2D
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE3D(K3D),KTYPE2D(K2D)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGRID(YDGEOMETRY%YRDIM%NPROMA,K3D*YDGEOMETRY%YRDIMV%NFLEVG+K2D,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPSP(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE_INTERP
LOGICAL           ,INTENT(IN)    :: LDNEW
LOGICAL           ,INTENT(IN)    :: LDCLOSE

INTEGER(KIND=JPIM) :: JF,ILEV,JLEV,IOPROC,JBLK,IOFF,IOFFG
INTEGER(KIND=JPIM) :: IGPOUT
INTEGER(KIND=JPIM) :: IOPROCS(NPROC),ILEVLIST(YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: IGRIB2D(YDGEOMETRY%YRDIMV%NFLEVG*K3D),ILEVS2D(YDGEOMETRY%YRDIMV%NFLEVG*K3D),&
 & IPROC2D(YDGEOMETRY%YRDIMV%NFLEVG*K3D)
INTEGER(KIND=JPIM) :: ILEVSSFC(K2D)
INTEGER(KIND = JPIM) :: IONPROCS,ICHUNKS,JROC,JCH,ISTLEV,IENLEV,ILEVS
REAL(KIND=JPRB),ALLOCATABLE :: ZGRID(:,:,:)
LOGICAL :: LLINTERP
LOGICAL :: LLFIRST=.TRUE.
CHARACTER(LEN=128),SAVE :: CLFILELEV(JPMXLE)
CHARACTER(LEN=128),SAVE :: CLFILESURF
CHARACTER(LEN=128) :: CLFILE
CHARACTER(LEN=3) :: CLEV
CHARACTER(LEN=1) :: CLMODE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
TYPE(TYPE_IOSTREAM),SAVE :: YL_IOSTREAM_SURF
TYPE(TYPE_IOSTREAM),SAVE :: YL_IOSTREAM(JPMXLE)
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST
INTEGER(KIND = JPIM):: ISDML,  IRSNML, ITSNML, IWSNML, IASNML 
INTEGER(KIND = JPIM):: IGRIBCD

#include "abor1.intfb.h"
#include "grid_psglobal.intfb.h"
! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRITE_GRID_TRAJ',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, &
 & NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOTG=>YDGEM%NGPTOTG)
IF(LLFIRST) THEN
  CLFILESURF='CLOSED'
  CLFILELEV(:)='CLOSED'
  LLFIRST=.FALSE.
ENDIF
IF (LDNEW) THEN
  CLMODE='w'
ELSE
  CLMODE='a'
ENDIF

CALL INQUIRE_GRIB(KGRIB_HD,KPTS=IGPOUT)

LLINTERP=(NGPTOTG /= IGPOUT)
IF (LLINTERP) THEN
  SELECT CASE (KTYPE_INTERP)
  CASE(1)
    WRITE(NULOUT,*)'WRITE_GRID_TRAJ: doing linear interpolation'
  CASE(2)
    WRITE(NULOUT,*)'WRITE_GRID_TRAJ: doing cubic interpolation'
  CASE(3)
    WRITE(NULOUT,*)'WRITE_GRID_TRAJ: doing conserving interpolation'
  CASE(4)
    WRITE(NULOUT,*)'WRITE_GRID_TRAJ: doing conserving interpolation w/o pressure'
  CASE DEFAULT
    CALL ABOR1('WRITE_GRID_TRAJ: unkown interpolation type')
  END SELECT
ENDIF

!-----------------------------------------------------------!
! To be cleaner, the definition of the grid should be       !
! passed as an argument, not a grid header.              YT !
!-----------------------------------------------------------!
! There should be a test in case the two grids are the same !
! set linterp to false in that case.                        !
! This test is a little crude...                         YT !
!-----------------------------------------------------------!

! Processor NPROC handles the surface fields

IF (K2D > 0) THEN
  IOPROC = NPROC
  IF (K3D>0) THEN
    ILEV = 0
    WRITE(CLEV,'(I3.3)')ILEV
    CLFILE = CDFILE//'L'//CLEV
  ELSE
    CLFILE = CDFILE
  ENDIF
  IF(LDNEW .OR. CLFILE /= CLFILESURF ) THEN
    IF(TRIM(CLFILESURF) /= 'CLOSED') THEN
      CALL CLOSE_IOSTREAM(YL_IOSTREAM_SURF)
      CLFILESURF='CLOSED'
    ENDIF
    CLFILESURF = CLFILE
    CALL SETUP_IOSTREAM(YL_IOSTREAM_SURF,'CIO',TRIM(CLFILESURF),CDMODE=CLMODE,&
   & KIOMASTER=IOPROC)
  ENDIF
! set level list of surface single/multi-layer
  ISDML  = 0
  IRSNML = 0
  ITSNML = 0
  IWSNML = 0
  IASNML = 0
  ILEVSSFC(1:K2D) = 0
  DO JF = 1, K2D
    IGRIBCD = KTYPE2D(JF)
    IF (IGRIBCD == NGRBSD) THEN
      ISDML = ISDML+1
      ILEVSSFC(JF) = ISDML
    ENDIF
    IF (IGRIBCD == NGRBRSN) THEN
      IRSNML = IRSNML+1
      ILEVSSFC(JF) = IRSNML
    ENDIF
    IF (IGRIBCD == NGRBTSN) THEN
      ITSNML = ITSNML+1
      ILEVSSFC(JF) = ITSNML
    ENDIF
    IF (IGRIBCD == NGRBWSN) THEN
      IWSNML = IWSNML+1
      ILEVSSFC(JF) = IWSNML
    ENDIF
    IF (IGRIBCD == NGRBASN) THEN
      IASNML = IASNML+1
      ILEVSSFC(JF) = IASNML
    ENDIF
  ENDDO
  CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,KRESOL=NRESOL,&
   & KGRIB2D=KTYPE2D(1:K2D), KLEVS2D = ILEVSSFC(1:K2D),CDLEVTYPE='SFC', KPROMA=NPROMA,&
   & LDINTERP=LLINTERP,KTYPE_INTERP=KTYPE_INTERP,KGRIB_HANDLE=KGRIB_HD,KGRIB_HANDLE2=KGRIB_HD2,&
   & KCHUNKSIZE=K2D,PTSTEP=YDRIP%TSTEP)

  CALL IO_PUT(YL_IOSTREAM_SURF,YL_IOREQUEST,&
   & PR3=PGRID(:,K3D*NFLEVG+1:K3D*NFLEVG+K2D,:))
  CALL CLOSE_IOREQUEST(YL_IOREQUEST)
  IF(LDCLOSE) THEN
    CALL CLOSE_IOSTREAM(YL_IOSTREAM_SURF)
    CLFILESURF='CLOSED'
  ENDIF
ENDIF

IF (K3D>0) THEN
  IF (KTYPE_INTERP==3) THEN
!   Create global gridpoint surface pressure on each processor
    IF (.NOT.ALLOCATED(RGPPRS_HR)) ALLOCATE(RGPPRS_HR(NGPTOTG))
    CALL GRID_PSGLOBAL(YDGEOMETRY,PSPSP)
  ENDIF
  IONPROCS = MIN(NPROC,NFLEVG)
  ICHUNKS  = (NFLEVG-1)/IONPROCS+1
  DO JROC=1,IONPROCS
    IOPROCS(JROC) = JROC
  ENDDO
  DO JROC=IONPROCS+1,NPROC
    CLFILE='XXXXXXX'
  ENDDO
  DO JCH=1,ICHUNKS
    ISTLEV = (JCH-1)*IONPROCS+1
    IENLEV = MIN(NFLEVG,JCH*IONPROCS)
    ILEVS  = IENLEV-ISTLEV+1
    DO JLEV=ISTLEV,IENLEV
      
      IOPROC = JLEV-ISTLEV+1
      ILEVLIST(IOPROC) = JLEV

      IF (IOPROC == MYPROC) THEN
        WRITE(CLEV,'(I3.3)')JLEV
        CLFILE=CDFILE//'L'//CLEV
      ENDIF
    ENDDO
    IF(LDNEW .OR. CLFILE /= CLFILELEV(JCH)) THEN
      IF(TRIM(CLFILELEV(JCH)) /= 'CLOSED') THEN
        CALL CLOSE_IOSTREAM(YL_IOSTREAM(JCH))
        CLFILELEV(JCH)='CLOSED'
      ENDIF
      CLFILELEV(JCH)=CLFILE   
      CALL SETUP_IOSTREAM(YL_IOSTREAM(JCH),'CIO',TRIM(CLFILELEV(JCH)),CDMODE=CLMODE,&
     & KIOPROCS=IOPROCS(1:ILEVS))
    ENDIF
    ALLOCATE(ZGRID(NPROMA,K3D*ILEVS,NGPBLKS))
    DO JF=1,K3D
      IOFF = (JF-1)*ILEVS
      IGRIB2D(IOFF+1:IOFF+ILEVS) = KTYPE3D(JF)
      ILEVS2D(IOFF+1:IOFF+ILEVS) = ILEVLIST(1:ILEVS)
      IPROC2D(IOFF+1:IOFF+ILEVS) = IOPROCS(1:ILEVS)
    ENDDO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JBLK,JF,IOFF,IOFFG)
    DO JBLK=1,NGPBLKS
      DO JF=1,K3D
        IOFF = (JF-1)*ILEVS
        IOFFG = (JF-1)*NFLEVG
        ZGRID(:,IOFF+1:IOFF+ILEVS,JBLK)=  PGRID(:,IOFFG+ISTLEV:IOFFG+IENLEV,JBLK)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

    CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,KRESOL=NRESOL,&
     & KGRIB2D=IGRIB2D(1:K3D*ILEVS),KLEVS2D=ILEVS2D(1:K3D*ILEVS),CDLEVTYPE='ML',&
     & KPROMA=NPROMA,LDINTERP=LLINTERP,KTYPE_INTERP=KTYPE_INTERP,KGRIB_HANDLE=KGRIB_HD,PTSTEP=YDRIP%TSTEP)
    CALL IO_PUT(YL_IOSTREAM(JCH),YL_IOREQUEST,&
     & PR3=ZGRID,KFLDPROC2D=IPROC2D(1:K3D*ILEVS))
    CALL CLOSE_IOREQUEST(YL_IOREQUEST)
    IF(LDCLOSE) THEN
      CALL CLOSE_IOSTREAM(YL_IOSTREAM(JCH))
      CLFILELEV(JCH)='CLOSED'
    ENDIF
    DEALLOCATE(ZGRID)
  ENDDO
ENDIF

IF (ALLOCATED(RGPPRS_HR)) DEALLOCATE(RGPPRS_HR)
IF (ALLOCATED(RGPPRS_LR)) DEALLOCATE(RGPPRS_LR)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRITE_GRID_TRAJ',1,ZHOOK_HANDLE)
END SUBROUTINE WRITE_GRID_TRAJ
