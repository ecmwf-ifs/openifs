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

MODULE YOMJBECV

!     Purpose.
!     --------
!       Data and controls for extended control variable

!     Author.
!     -------
!       S. Massart

!     Modifications.
!     --------------
!       Original    April-2018
! ------------------------------------------------------------------


USE PARKIND1, ONLY: JPIM, JPRB
USE TYPE_ECV, ONLY: TECVGRIB, TECV_CONFIG, TECVDIM, ECV_CONTAINER

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!   LJB_ALPHA_CV    True if hybrid B using Alpha Control Variable
!   LSKTECV         True if skin temperature in ECV
!   LECPHYSPARECV   True if optimisation of parameters
!   LINVERACV       True if emission inversion (CAMS)
!
!   YRECVGRIB       Grib number of the ECV fields
!   YRECVCONFIG     ECV configuration
!   YRECV0          Background of ECV_CONTAINER
!   YRECV5          Trajectory of ECV_CONTAINER

!   NDIAECV         Level of output diagnostics
!   NINTERPECV      Interpolation type
!   CFNECVHRSPEC    Filename of the spectral ECV 'trajectory'
!   CFNECVHRGRID    Filename of the gridpoint ECV 'trajectory'

!     ------------------------------------------------------------------

LOGICAL                         :: LJB_ALPHA_CV
LOGICAL                         :: LSKTECV
LOGICAL                         :: LINVERACV
LOGICAL                         :: LECPHYSPARECV

TYPE(TECVGRIB)                  :: YRECVGRIB
TYPE(TECV_CONFIG)               :: YRECVCONFIG
TYPE(ECV_CONTAINER)             :: YRECV0
TYPE(ECV_CONTAINER)             :: YRECV5

INTEGER(KIND=JPIM)             :: NDIAECV
INTEGER(KIND=JPIM), PARAMETER  :: NINTERPECV = 4
CHARACTER(LEN=16),  PARAMETER  :: CFNECVHRSPEC = 'TRAJHR00/ecvspec'
CHARACTER(LEN=16),  PARAMETER  :: CFNECVHRGRID = 'TRAJHR00/ecvgrid'

CHARACTER(LEN =20), PARAMETER  :: CKNOWNECV(5) = (/  &
  &   'SOLAR_CONSTANT      ', &
  &   'SKTECV              ', &
  &   'JB_HYBRID_ALPHACV   ', &
  &   'FLUX_INVER          ', &
  &   'EC_PHYS             '/)

!-----------------------------------------------------------------------

INTERFACE ALLOCATE_ECV
MODULE PROCEDURE CREATE_ECV, COPY_ECV
END INTERFACE

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

SUBROUTINE CREATE_ECV(YDGEOMETRY,YDDIMECV,YDECV)

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR      , ONLY : LECV

IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN)     :: YDGEOMETRY
TYPE(TECVDIM)      ,INTENT(IN)     :: YDDIMECV
TYPE(ECV_CONTAINER),INTENT(INOUT)  :: YDECV
call abor1("oifs/fc-only - CREATE_ECV should never be called")


END SUBROUTINE CREATE_ECV

!-----------------------------------------------------------------------

SUBROUTINE COPY_ECV(YDGEOMETRY,YDDIMECV,YDECV,YDOTHER)

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR      , ONLY : LECV

IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN)     :: YDGEOMETRY
TYPE(TECVDIM)      ,INTENT(IN)     :: YDDIMECV
TYPE(ECV_CONTAINER),INTENT(OUT)    :: YDECV
TYPE(ECV_CONTAINER),INTENT(IN)     :: YDOTHER
call abor1("oifs/fc-only - COPY_ECV should never be called")

END SUBROUTINE COPY_ECV

!-----------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ECV(YDECV)

USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(ECV_CONTAINER),INTENT(INOUT)  :: YDECV

REAL(KIND=JPHOOK)              :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMJBECV:DEALLOCATE_ECV',0,ZHOOK_HANDLE)

IF (ALLOCATED(YDECV%RECV1D)) DEALLOCATE(YDECV%RECV1D)
IF (ALLOCATED(YDECV%RSPECV2D)) DEALLOCATE(YDECV%RSPECV2D)
IF (ALLOCATED(YDECV%RGPECV2D)) DEALLOCATE(YDECV%RGPECV2D)
IF (ALLOCATED(YDECV%RSPECV3D)) DEALLOCATE(YDECV%RSPECV3D)
IF (ALLOCATED(YDECV%RGPECV3D)) DEALLOCATE(YDECV%RGPECV3D)

IF (LHOOK) CALL DR_HOOK('YOMJBECV:DEALLOCATE_ECV',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLOCATE_ECV

!-----------------------------------------------------------------------

SUBROUTINE SPNORM_ECV(YDGEOMETRY,YDDIMECV,YDECV,CDLABEL)

USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMVAR      , ONLY : LECV
USE YOMLUN      , ONLY : NULOUT
USE YOMMP0      , ONLY : NPROC

IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN) :: YDGEOMETRY
TYPE(TECVDIM)      ,INTENT(IN) :: YDDIMECV
TYPE(ECV_CONTAINER),INTENT(IN) :: YDECV
CHARACTER(LEN=16)  ,INTENT(IN) :: CDLABEL
call abor1("oifs/fc-only - SPNORM_ECV should never be called")

END SUBROUTINE SPNORM_ECV

!-----------------------------------------------------------------------

SUBROUTINE ZERO_ECV(YDGEOMETRY,YDDIMECV,YDECV)


USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR      , ONLY : LECV
USE GEOMETRY_MOD, ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN)     :: YDGEOMETRY
TYPE(TECVDIM)      ,INTENT(IN)     :: YDDIMECV
TYPE(ECV_CONTAINER),INTENT(INOUT)  :: YDECV
call abor1("oifs/fc-only - ZERO_ECV should never be called")

END SUBROUTINE ZERO_ECV

!-----------------------------------------------------------------------
SUBROUTINE READ_FG_ECV(YDGEOMETRY,YDDIMECV,YDSURF,LDBCK)

USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR             , ONLY : LECV
USE YOMARG             , ONLY : CNMEXP
USE YOMLUN             , ONLY : NULERR
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ALLOCATE_SPEC, DEALLOCATE_SPEC
USE SURFACE_FIELDS_MIX , ONLY : TSURF

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN) :: YDGEOMETRY
TYPE(TSURF),OPTIONAL,INTENT(IN) :: YDSURF
TYPE(TECVDIM)     ,INTENT(IN) :: YDDIMECV
LOGICAL, OPTIONAL    ,INTENT(IN) :: LDBCK
call abor1("oifs/fc-only - READ_FG_ECV should never be called")


END SUBROUTINE READ_FG_ECV

!-----------------------------------------------------------------------
SUBROUTINE SAVE_FG_ECV(YDGEOMETRY,YDRIP,YDDIMECV,YDECV)

USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR             , ONLY : LECV, MUPTRA
USE YOMTRAJ            , ONLY : NSMAX_TRAJ, MAIN_GRIB
USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ALLOCATE_SPEC, DEALLOCATE_SPEC
USE YOMLUN             , ONLY : NULOUT


IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN) :: YDGEOMETRY
TYPE(TRIP)         ,INTENT(IN) :: YDRIP
TYPE(TECVDIM)      ,INTENT(IN) :: YDDIMECV
TYPE(ECV_CONTAINER),INTENT(IN) :: YDECV
call abor1("oifs/fc-only - SAVE_FG_ECV should never be called")


END SUBROUTINE SAVE_FG_ECV

!-----------------------------------------------------------------------

SUBROUTINE RANDOM_ECV(YDGEOMETRY,YDDIMECV,KSEED,YDECV)

USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR            , ONLY : LECV
USE YOMMP0            , ONLY : MYSETW
USE GEOMETRY_MOD      , ONLY : GEOMETRY
USE MPL_MODULE        , ONLY : MPL_ALLREDUCE
USE RANDOM_NUMBERS_MIX, ONLY : RANDOMNUMBERSTREAM,  &
                             & INITIALIZE_RANDOM_NUMBERS, GAUSSIAN_DISTRIBUTION

IMPLICIT NONE

TYPE(GEOMETRY)     , INTENT(IN)    :: YDGEOMETRY
TYPE(TECVDIM)      ,INTENT(IN)     :: YDDIMECV
INTEGER (KIND=JPIM), INTENT(INOUT) :: KSEED
TYPE(ECV_CONTAINER), INTENT(INOUT) :: YDECV
call abor1("oifs/fc-only - RANDOM_ECV should never be called")

END SUBROUTINE RANDOM_ECV

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE ECV_SPTOGP(YDGEOMETRY,YDECV,KECV,PGP2,PGP3)

!     Purpose.
!     --------
!       Convert 2D or 3D spectral ECV field to gridpoint


USE GEOMETRY_MOD, ONLY : GEOMETRY
USE TYPE_ECV,     ONLY : ECV_CONTAINER
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)  :: YDGEOMETRY
TYPE(ECV_CONTAINER)     ,INTENT(IN)  :: YDECV
INTEGER(KIND=JPIM)      ,INTENT(IN)  :: KECV
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP2(YDGEOMETRY%YRDIM%NPROMA,1,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP3(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS)
call abor1("oifs/fc-only - ECV_SPTOGP should never be called")


END SUBROUTINE ECV_SPTOGP

!-----------------------------------------------------------------------

SUBROUTINE ECV_SPTOGPAD(YDGEOMETRY,YDECV,KECV,PGP2,PGP3)

!     Purpose.
!     --------
!       Convert 2D or 3D spectral ECV field to gridpoint (adjoint version)


USE GEOMETRY_MOD, ONLY : GEOMETRY
USE TYPE_ECV,     ONLY : ECV_CONTAINER
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)    :: YDGEOMETRY
TYPE(ECV_CONTAINER)     ,INTENT(INOUT) :: YDECV
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KECV
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PGP2(YDGEOMETRY%YRDIM%NPROMA,1,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PGP3(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS)
call abor1("oifs/fc-only - ECV_SPTOGPAD should never be called")

END SUBROUTINE ECV_SPTOGPAD

!-----------------------------------------------------------------------

SUBROUTINE ECV_SPTOGP_INV(YDGEOMETRY,YDECV,KECV,PGP2,PGP3)

!     Purpose.
!     --------
!       Convert 2D or 3D spectral ECV field to gridpoint


USE GEOMETRY_MOD, ONLY : GEOMETRY
USE TYPE_ECV,     ONLY : ECV_CONTAINER
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)   :: YDGEOMETRY
TYPE(ECV_CONTAINER)     ,INTENT(INOUT):: YDECV
INTEGER(KIND=JPIM)      ,INTENT(IN)   :: KECV
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)   :: PGP2(YDGEOMETRY%YRDIM%NPROMA,1,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)   :: PGP3(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS)
call abor1("oifs/fc-only - ECV_SPTOGP_INV should never be called")

END SUBROUTINE ECV_SPTOGP_INV

!-----------------------------------------------------------------------

SUBROUTINE ECV_SPTOGP_INVAD(YDGEOMETRY,YDECV,KECV,PGP2,PGP3)

!     Purpose.
!     --------
!       Convert 2D or 3D spectral ECV field to gridpoint


USE GEOMETRY_MOD, ONLY : GEOMETRY
USE TYPE_ECV,     ONLY : ECV_CONTAINER
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)   :: YDGEOMETRY
TYPE(ECV_CONTAINER)     ,INTENT(INOUT):: YDECV
INTEGER(KIND=JPIM)      ,INTENT(IN)   :: KECV
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT):: PGP2(YDGEOMETRY%YRDIM%NPROMA,1,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT):: PGP3(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS)
call abor1("oifs/fc-only - ECV_SPTOGP_INVAD should never be called")


END SUBROUTINE ECV_SPTOGP_INVAD

!-----------------------------------------------------------------------

SUBROUTINE ECV_2DSP_TO_1D(YDGEOMETRY,KECV,PSP2D,P1D)

USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE MPL_MODULE  , ONLY : MPL_ALLREDUCE, MPL_BROADCAST
USE YOMMP0      , ONLY : MYSETW, MYPROC, NPROC

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)  :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)  :: KECV
REAL(KIND=JPRB),   INTENT(IN)  :: PSP2D(:,:)
REAL(KIND=JPRB),   INTENT(OUT) :: P1D
call abor1("oifs/fc-only - ECV_2DSP_TO_1D should never be called")

END SUBROUTINE ECV_2DSP_TO_1D

!-----------------------------------------------------------------------

SUBROUTINE ECV_2DSP_TO_1D_AD(YDGEOMETRY,KECV,PSP2D,P1D)

USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE GEOMETRY_MOD      , ONLY : GEOMETRY
USE MPL_MODULE        , ONLY: MPL_ALLREDUCE, MPL_BROADCAST
USE YOMMP0            , ONLY : MYSETW, MYPROC   ,NPROC

IMPLICIT NONE

TYPE(GEOMETRY) , INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)  :: KECV
REAL(KIND=JPRB), INTENT(INOUT) :: PSP2D(:,:)
REAL(KIND=JPRB), INTENT(INOUT)  :: P1D
call abor1("oifs/fc-only - ECV_2DSP_TO_1D_AD should never be called")

END SUBROUTINE ECV_2DSP_TO_1D_AD

!-----------------------------------------------------------------------

SUBROUTINE ECV_1D_TO_2DSP(YDGEOMETRY,KECV,P1D,PSP2D)

USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE GEOMETRY_MOD      , ONLY : GEOMETRY

IMPLICIT NONE

TYPE(GEOMETRY) , INTENT(IN)  :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN):: KECV
REAL(KIND=JPRB), INTENT(IN)  :: P1D
REAL(KIND=JPRB), INTENT(OUT) :: PSP2D(:,:)
call abor1("oifs/fc-only - ECV_2DSP_TO_1D_AD should never be called")

END SUBROUTINE ECV_1D_TO_2DSP

!-----------------------------------------------------------------------

END MODULE YOMJBECV
