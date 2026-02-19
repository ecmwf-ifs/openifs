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

SUBROUTINE SUGRIDA_FIX_TOZ(YDGEOMETRY,YDOZO,YDMCC,YDDPHY,YDPHY,KFIELDS,PFIELD,YDFLDSC)

!**** *SUGRIDA_FIX_TOZ*  - Fix ozone fields after reading them

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 11-09-2012 from sugrida.F90

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!-----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHY       , ONLY : TPHY
USE YOMDPHY      , ONLY : TDPHY
USE YOMOZO       , ONLY : TOZO
USE YOMMCC       , ONLY : TMCC
USE YOMLUN       , ONLY : NULOUT
USE YOMMP0       , ONLY : MY_REGION_NS, MY_REGION_EW
USE IOFLDDESC_MOD, ONLY : IOFLDDESC

IMPLICIT NONE

TYPE(GEOMETRY)     , INTENT (IN)    :: YDGEOMETRY
TYPE(TDPHY)        , INTENT (IN)    :: YDDPHY
TYPE(TMCC)         , INTENT (IN)    :: YDMCC
TYPE(TOZO)         , INTENT (INOUT) :: YDOZO
TYPE(TPHY)         , INTENT (IN)    :: YDPHY
INTEGER (KIND=JPIM), INTENT (IN)    :: KFIELDS
REAL (KIND=JPRB)   , INTENT (IN)    :: PFIELD (YDGEOMETRY%YRGEM%NGPTOT,KFIELDS)
TYPE (IOFLDDESC)   , INTENT (IN)    :: YDFLDSC (KFIELDS) 

#include "abor1.intfb.h"

INTEGER (KIND=JPIM) :: JFLEV
INTEGER (KIND=JPIM) :: ICH, IROF, JROF, JGL, JV, IFLD, JFLD
INTEGER (KIND=JPIM) :: INUMLAT
REAL (KIND=JPRB) :: ZMAX, ZMIN, ZMOY
CHARACTER (LEN=16)  :: CLSUO3
CHARACTER (LEN=*), PARAMETER :: CLFMT = "(1X,'ARRAY ',A10,' ALLOCATED ',8I8)"

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


IF (LHOOK) CALL DR_HOOK ('SUGRIDA_FIX_TOZ',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NVCLIS=>YDDPHY%NVCLIS, NTOZ2D=>YDDPHY%NTOZ2D, NTOZ3D=>YDDPHY%NTOZ3D, &
 & NDGENL=>YDDIM%NDGENL, NDGSAL=>YDDIM%NDGSAL, &
 & NPTRFLOFF=>YDMP%NPTRFLOFF, NFRSTLAT=>YDMP%NFRSTLAT, NONL=>YDMP%NONL, &
 & NLSTLAT=>YDMP%NLSTLAT, &
 & NSTTYP=>YDGEM%NSTTYP, NGPTOT=>YDGEM%NGPTOT, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LOZONE=>YDPHY%LOZONE, LMCC01=>YDMCC%LMCC01)
IF ((LOZONE .OR. (NVCLIS == 1)) .AND. (.NOT.LMCC01)) THEN

  INUMLAT = NLSTLAT(MY_REGION_NS)-NFRSTLAT(MY_REGION_NS)+1


!*  There should be 7 coefficients
  IF (NVCLIS /= 7.AND.LOZONE) THEN
    WRITE (NULOUT,*) ' WRONG NVCLIS (it must be 7) : ', NVCLIS
    CALL ABOR1('SUGRIDA : ABOR1 CALLED (wrong NVCLIS)')
  ENDIF

!* Open 3D files if they are not allocated
  IF (NVCLIS * NTOZ3D > 0) THEN
    IF (.NOT. ALLOCATED (YDOZO%TOZ3DBL)) THEN
      ALLOCATE (YDOZO%TOZ3DBL (NGPTOT, NFLEVG*NVCLIS*NTOZ3D))
      WRITE(NULOUT,CLFMT) 'TOZ3DBL  ', SIZE(YDOZO%TOZ3DBL), SHAPE(YDOZO%TOZ3DBL)
    ENDIF
  ENDIF

!* Open 2D files if they are not allocated
  IF (NVCLIS*NTOZ2D > 0) THEN
    IF (.NOT. ALLOCATED (YDOZO%TOZ2DL)) THEN
      ALLOCATE (YDOZO%TOZ2DL (NDGSAL:NDGENL, NFLEVG*NVCLIS*NTOZ2D))
      WRITE(NULOUT,CLFMT) 'TOZ2DL   ', SIZE(YDOZO%TOZ2DL), SHAPE(YDOZO%TOZ2DL)
    ENDIF
  ENDIF

!* Start the vertical loop for the NFLEVG levels  
  DO JFLEV = 1, NFLEVG

    DO JV = 1, NVCLIS

      WRITE (UNIT=CLSUO3, FMT='(''OZONE.A'',I1)') JV

      IFLD = -1
      DO JFLD = 1, KFIELDS
        IF ((YDFLDSC (JFLD)%CPREF == 'S') .AND. (YDFLDSC (JFLD)%ILEVG == JFLEV)&
          & .AND. (YDFLDSC (JFLD)%CSUFF == CLSUO3)) THEN
          IFLD = JFLD
          EXIT
        ENDIF
      ENDDO

      IF (IFLD < 0) CALL ABOR1 ('SUGRIDA_FIX_TOZ: FIELD WAS NOT FOUND :'//CLSUO3)

      ICH = (JV-1)*NFLEVG + JFLEV

! FULL 3D COMMON (Longitude x Latitude x Level)
      IF (NTOZ3D == 1) THEN 

        DO JROF = 1, NGPTOT
          YDOZO%TOZ3DBL (JROF,ICH) = PFIELD (JROF, IFLD)
        ENDDO

        ZMIN = YDOZO%TOZ3DBL (1, ICH)
        ZMAX = YDOZO%TOZ3DBL (1, ICH)
        ZMOY = 0.0_JPRB

        DO JROF = 1, NGPTOT
          ZMOY = ZMOY+YDOZO%TOZ3DBL(JROF,ICH)/REAL(NGPTOT)
          ZMIN = MIN (ZMIN,YDOZO%TOZ3DBL(JROF,ICH))
          ZMAX = MAX (ZMAX,YDOZO%TOZ3DBL(JROF,ICH))
        ENDDO

        WRITE(NULOUT,'('' OZONE-3D A'',I1,'' LEVEL '',I3,&
         & '' MIN/MOY/MAX : '', 3(1x,E12.5) )')&
         & JV,JFLEV,ZMIN,ZMOY,ZMAX

      ENDIF   

! THE 2D COMMON (Latitude x Level)
      IF (NTOZ2D == 1) THEN

        IF (NSTTYP /= 1) THEN
          WRITE(NULOUT,&
           & FMT='('' SUGRIDA&
           & ! OZONE INCONSISTENCY : NTOZ2D='',I4 ,'' NSTTYP='',I4)') NTOZ2D, NSTTYP
          CALL ABOR1(' SUGRIDA : ABOR1 CALLED')
        ENDIF

        IROF = 1
        DO JGL = 1, INUMLAT
          YDOZO%TOZ2DL(JGL,ICH) = PFIELD (IROF, IFLD)
          IROF = IROF+NONL(NPTRFLOFF+JGL,MY_REGION_EW)
        ENDDO

        ZMIN = YDOZO%TOZ2DL (1,ICH)
        ZMAX = YDOZO%TOZ2DL (1,ICH)
        ZMOY = 0.0_JPRB

        DO JROF = 1, INUMLAT
          ZMOY = ZMOY+YDOZO%TOZ2DL(JROF,ICH)/REAL(INUMLAT)
          ZMIN = MIN (ZMIN,YDOZO%TOZ2DL(JROF,ICH))
          ZMAX = MAX (ZMAX,YDOZO%TOZ2DL(JROF,ICH))
        ENDDO

        WRITE(NULOUT,'('' OZONE-2D A'',I1,'' LEVEL '',I3,&
         & '' MIN/MOY/MAX : '', 3(1x,E12.5) )')&
         & JV,JFLEV, ZMIN,ZMOY,ZMAX

      ENDIF

    ENDDO   ! JV=1,NVCLIS

  ENDDO   ! JFLEV=1,NFLEVG

ENDIF   ! LOZONE=.T.

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK ('SUGRIDA_FIX_TOZ',1,ZHOOK_HANDLE)

END SUBROUTINE SUGRIDA_FIX_TOZ 
