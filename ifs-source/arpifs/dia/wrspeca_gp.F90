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

SUBROUTINE WRSPECA_GP(YDGEOMETRY,YDXFU,YDML_GCONF,YDDYN,YDDYNA,YDSPEC,YDGFL,YDSURF,YDFACTX,CDFIC,YDMCUF)

!**** *WRSPECA_GP*  - Write spectral fields in grid-point format

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 11-09-2012

!     Modifications:
!     --------------
!      P. Marguinaud : 10-10-2013 : Use FACTX
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, and YRSURF are now passed by argument
!      O. Marsden: Sept 2016 Removed use of SPA3, replaced by explicit spectral argument
!-----------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN                 , ONLY : TDYN
USE YOMDYNA                , ONLY : TDYNA
USE YOMXFU                 , ONLY : TXFU
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX     , ONLY : TSURF
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPRB, JPIM
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN                 , ONLY : NULOUT
USE YOMCT0                 , ONLY : LELAM, LWRSPECA_GP_UV
USE YOMMP0                 , ONLY : MYSETV, NPRTRV
USE YOMTAG                 , ONLY : MTAG_MFIO_WRSPECA_GP
USE IOSPECA_MOD            , ONLY :&
                                 & IOSPECA_SELECTD,&
                                 & IOSPECA_SELECTF,&
                                 & IOSPECA_CTX,&
                                 & IOSPECA_START,&
                                 & IOSPECA_COUNT,&
                                 & IOSPECA_FINISH,&
                                 & IOSPECA_VSETOFF,&
                                 & NSPECACT_WRITE_GP,&
                                 & IOSPECA_PFCF2UV_GP,&
                                 & IOSPECA_FLDDESC_SP2GP
USE IOFLDDESC_MOD          , ONLY : IOFLDDESC
USE FACTX_MOD              , ONLY : FACTX
USE SPECTRAL_FIELDS_MOD    , ONLY : SPECTRAL_FIELD
USE YOMMCUF                , ONLY : TMCUF

IMPLICIT NONE

TYPE(GEOMETRY)               , INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)                   , INTENT(INOUT) :: YDDYN
TYPE(TDYNA)                  , INTENT(INOUT) :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(INOUT) :: YDML_GCONF
TYPE(TXFU)                   , INTENT(INOUT) :: YDXFU
TYPE(SPECTRAL_FIELD)         , INTENT(INOUT) :: YDSPEC
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL
TYPE(TSURF)                  , INTENT(INOUT) :: YDSURF
TYPE (FACTX)                 , INTENT(INOUT) :: YDFACTX
CHARACTER(LEN=*)             , INTENT(IN)    :: CDFIC
TYPE(TMCUF)                  , INTENT(INOUT), OPTIONAL :: YDMCUF

#include "inv_trans.h"
#include "einv_trans_px.intfb.h"
#include "wrgp2fa.intfb.h"

TYPE (IOFLDDESC), ALLOCATABLE :: YLFLDSC_SP (:), YLFLDSC_GP (:)
REAL (KIND=JPRB), ALLOCATABLE :: ZSPBUFL (:,:) 
REAL (KIND=JPRB), ALLOCATABLE :: ZGPBUFL (:,:,:)  
INTEGER (KIND=JPIM) :: IFLDSCH (NPRTRV), IVSETOFF (NPRTRV)
INTEGER (KIND=JPIM) :: IFLDSPG  ! Total number of fields
INTEGER (KIND=JPIM) :: IFLDSPL  ! Number of fields in this V-set
INTEGER (KIND=JPIM) :: IFLGPOFF ! Offset of requested fields in grid-point array
INTEGER (KIND=JPIM) :: IUGPOFF 
INTEGER (KIND=JPIM) :: IVGPOFF 

LOGICAL :: LLUV ! U/V in output file instead of flow potential/current function
                ! This is for the global model only; AROME/ALADIN save U/V

TYPE (IOSPECA_CTX) :: YLIOSP

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('WRSPECA_GP',0,ZHOOK_HANDLE)


IF (.NOT. LELAM) THEN
  WRITE (NULOUT, *) "WRSPECA_GP: LWRSPECA_GP_UV = ", LWRSPECA_GP_UV,&
  & " (SAVE U/V INSTEAD OF PHI/KHI) "
  LLUV = LWRSPECA_GP_UV
ELSE
  LLUV = .FALSE.
ENDIF

CALL IOSPECA_START(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,NSPECACT_WRITE_GP,YLIOSP,YDSPEC)

CALL IOSPECA_COUNT(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_WRITE_GP,IFLDSPG,YDMCUF=YDMCUF)

IF (IFLDSPG > 0) THEN

  ALLOCATE (YLFLDSC_SP (IFLDSPG), YLFLDSC_GP (IFLDSPG))

  CALL IOSPECA_SELECTD(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_WRITE_GP,YLFLDSC_SP,YDMCUF=YDMCUF)

! Convert spectral descriptors to gridpoint

  CALL IOSPECA_FLDDESC_SP2GP (YDGEOMETRY, YLFLDSC_SP, YLFLDSC_GP)

  CALL IOSPECA_VSETOFF (YLFLDSC_SP%IVSET, IFLDSCH, IVSETOFF)

  IFLDSPL = IFLDSCH (MYSETV)

  IF (LLUV) THEN
    IUGPOFF  = 0
    IVGPOFF  = YDGEOMETRY%YRDIMV%NFLEVG
    IFLGPOFF = 2 * YDGEOMETRY%YRDIMV%NFLEVG ! NFLEVG U 2D fields + NFLEVG V 2D fields
  ELSE
    IUGPOFF  = -1
    IVGPOFF  = -1
    IFLGPOFF = 0
  ENDIF

  ALLOCATE (ZSPBUFL (IFLDSPL, YDGEOMETRY%YRDIM%NSPEC2), ZGPBUFL (YDGEOMETRY%YRGEM%NGPTOT, IFLGPOFF + IFLDSPG, 1))

  CALL IOSPECA_SELECTF(YDGEOMETRY,YDML_GCONF%YGFL,YDDYN, NSPECACT_WRITE_GP, YLIOSP, IFLDSPL,&
                      & IVSETOFF (MYSETV), 0_JPIM, YLFLDSC_SP, ZSPBUFL, YDSPEC,YDMCUF=YDMCUF)

! Spectral to grid-point

  IF (LELAM) THEN

    CALL EINV_TRANS_PX (PSPSCALAR=ZSPBUFL, KVSETSC=YLFLDSC_SP%IVSET, KRESOL=YDGEOMETRY%YRDIM%NRESOL,&
                      & KPROMA=YDGEOMETRY%YRGEM%NGPTOT, PGP=ZGPBUFL)
                      
  ELSEIF (LLUV) THEN

    CALL INV_TRANS (PSPVOR=YDSPEC%VOR, PSPDIV=YDSPEC%DIV, KVSETUV=YDGEOMETRY%YRMP%NBSETLEV,&
                  & PSPSCALAR=ZSPBUFL, KVSETSC=YLFLDSC_SP%IVSET, KRESOL=YDGEOMETRY%YRDIM%NRESOL,&
                  & KPROMA=YDGEOMETRY%YRGEM%NGPTOT, PGP=ZGPBUFL)

  ELSE

    CALL INV_TRANS (PSPSCALAR=ZSPBUFL, KVSETSC=YLFLDSC_SP%IVSET, KRESOL=YDGEOMETRY%YRDIM%NRESOL,&
                  & KPROMA=YDGEOMETRY%YRGEM%NGPTOT, PGP=ZGPBUFL)

  ENDIF


  IF (LLUV) THEN

! Replace gridpoint potential flow/current function with U/V

    CALL IOSPECA_PFCF2UV_GP (YLFLDSC_SP, YLFLDSC_GP,&
                           & PGPBUFLU=ZGPBUFL (1:YDGEOMETRY%YRGEM%NGPTOT,IUGPOFF+1:IUGPOFF+YDGEOMETRY%YRDIMV%NFLEVG,1),&
                           & PGPBUFLV=ZGPBUFL (1:YDGEOMETRY%YRGEM%NGPTOT,IVGPOFF+1:IVGPOFF+YDGEOMETRY%YRDIMV%NFLEVG,1),&
                           & PGPBUFL=ZGPBUFL (1:YDGEOMETRY%YRGEM%NGPTOT,IFLGPOFF+1:IFLGPOFF+IFLDSPG,1))

  ENDIF

! Write gridpoint fields

  CALL WRGP2FA(YDGEOMETRY,YDML_GCONF%YRRIP,IFLDSPG,ZGPBUFL(1:YDGEOMETRY%YRGEM%NGPTOT,IFLGPOFF+1:IFLGPOFF+IFLDSPG,1),&
 &             YLFLDSC_GP,YDFACTX,CDFIC,KTAG=MTAG_MFIO_WRSPECA_GP)

  DEALLOCATE (YLFLDSC_SP, YLFLDSC_GP, ZSPBUFL, ZGPBUFL)

ENDIF

CALL IOSPECA_FINISH (YDGEOMETRY, YDDYNA,NSPECACT_WRITE_GP, YLIOSP, YDSPEC)


IF (LHOOK) CALL DR_HOOK ('WRSPECA_GP',1,ZHOOK_HANDLE)

END SUBROUTINE WRSPECA_GP

