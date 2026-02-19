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

SUBROUTINE SUSPECA_GP(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,YDSP,KFILE,LDATA,LDINOR,YDMCUF)

!**** *SUSPECA_GP*  - Read dynamic fields in grid-point format

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 11-09-2012

!     Modifications.
!     --------------
!      P. Marguinaud : 10-10-2013 : Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      P. Marguinaud : 10-10-2014 : Use SUSPECA_FIXUP
!------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN                 , ONLY : TDYN
USE YOMDYNA                , ONLY : TDYNA
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPRB, JPIM
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0                 , ONLY : N3DINI, LELAM, LSUSPECA_GP_UV
USE YOMLUN                 , ONLY : NULOUT
USE YOMMP0                 , ONLY : NPRTRV, MYSETV

USE IOSPECA_MOD            , ONLY : IOSPECA_COUNT,&
                                  & NSPECACT_READ_GP,&
                                  & IOSPECA_CTX,&
                                  & IOSPECA_START,&
                                  & IOSPECA_FINISH,&
                                  & IOSPECA_SELECTD,&
                                  & IOSPECA_VSETOFF,&
                                  & IOSPECA_SELECTF,&
                                  & IOSPECA_FLDDESC_SP2GP,&
                                  & IOSPECA_UV2PFCF_GP,&
                                  & IOSPECA_PFCF2UV_GP 

USE IOFLDDESC_MOD          , ONLY : IOFLDDESC
USE YEMLBC_MODEL             , ONLY : TELBC_MODEL
USE SPECTRAL_FIELDS_MOD
USE YOMMCUF                , ONLY : TMCUF

!------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(IN)    :: YDML_GCONF
TYPE(TDYN)                   , INTENT(IN)    :: YDDYN
TYPE(TDYNA)                  , INTENT(IN)    :: YDDYNA
TYPE(TELBC_MODEL)            , INTENT(IN)    :: YDML_LBC
TYPE(SPECTRAL_FIELD)         , INTENT(INOUT) :: YDSP
INTEGER(KIND=JPIM)           , INTENT(IN)    :: KFILE 
LOGICAL                      , INTENT(OUT)   :: LDATA 
LOGICAL                      , INTENT(IN)    :: LDINOR 
TYPE(TMCUF)                  , INTENT(INOUT), OPTIONAL :: YDMCUF
!------------------------------------------------------------------------
#include "suspeca_fixup.intfb.h"
#include "dir_trans.h"
#include "edir_trans_px.intfb.h"
#include "rdfa2gp.intfb.h"

TYPE (IOFLDDESC), ALLOCATABLE :: YLFLDSC_SP (:), YLFLDSC_GP (:)
REAL (KIND=JPRB), ALLOCATABLE :: ZSPBUFL (:,:) 
REAL (KIND=JPRB), ALLOCATABLE :: ZGPBUFL (:,:,:) 
INTEGER (KIND=JPIM) :: IFLDSCH (NPRTRV), IVSETOFF (NPRTRV)
INTEGER (KIND=JPIM) :: IFLDSPG  ! Total number of fields
INTEGER (KIND=JPIM) :: IFLDSPL  ! Number of fields in this V-set
INTEGER (KIND=JPIM) :: IFLGPOFF ! Offset of requested fields in grid-point array
INTEGER (KIND=JPIM) :: IUGPOFF 
INTEGER (KIND=JPIM) :: IVGPOFF 

LOGICAL :: LLUV

TYPE (IOSPECA_CTX) :: YLIOSP

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK ('SUSPECA_GP',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------

LDATA = .TRUE.

IF (.NOT. LELAM) THEN
  WRITE (NULOUT, *) "SUSPECA_GP: LSUSPECA_GP_UV = ", LSUSPECA_GP_UV,&
  & " (READ U/V INSTEAD OF PHI/KHI) "
  LLUV = LSUSPECA_GP_UV
ELSE
  LLUV = .FALSE.
ENDIF

CALL IOSPECA_COUNT(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_READ_GP,IFLDSPG,LDINOR=LDINOR,K3DINI=N3DINI,YDMCUF=YDMCUF)

IF (IFLDSPG > 0) THEN

  CALL IOSPECA_START(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,NSPECACT_READ_GP,YLIOSP,YDSP)

  ALLOCATE (YLFLDSC_SP (IFLDSPG), YLFLDSC_GP (IFLDSPG))

  CALL IOSPECA_SELECTD(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_READ_GP,YLFLDSC_SP,LDINOR=LDINOR,K3DINI=N3DINI,YDMCUF=YDMCUF)

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

  CALL IOSPECA_FLDDESC_SP2GP (YDGEOMETRY, YLFLDSC_SP, YLFLDSC_GP)

  IF (LLUV) CALL IOSPECA_PFCF2UV_GP (YLFLDSC_SP, YLFLDSC_GP)
  
! Read grid-point data

  CALL RDFA2GP (YDGEOMETRY, YDML_GCONF%YRRIP, IFLDSPG, ZGPBUFL(:,IFLGPOFF+1:IFLGPOFF+IFLDSPG,1), &
              & YLFLDSC_GP, KFILE=KFILE, YDML_LBC=YDML_LBC)

  IF (LLUV) THEN

! If U/V was read, then fill the beginning of ZGPBUFL to go back to VOR/DIV

    CALL IOSPECA_UV2PFCF_GP (YLFLDSC_GP,&
                           & PGPBUFLU=ZGPBUFL (1:YDGEOMETRY%YRGEM%NGPTOT,IUGPOFF+1:IUGPOFF+YDGEOMETRY%YRDIMV%NFLEVG,1),&
                           & PGPBUFLV=ZGPBUFL (1:YDGEOMETRY%YRGEM%NGPTOT,IVGPOFF+1:IVGPOFF+YDGEOMETRY%YRDIMV%NFLEVG,1),&
                           & PGPBUFL=ZGPBUFL (1:YDGEOMETRY%YRGEM%NGPTOT,IFLGPOFF+1:IFLGPOFF+IFLDSPG,1))

  ENDIF


! Grid-point to spectral

  IF (LELAM) THEN

    CALL EDIR_TRANS_PX (PSPSCALAR=ZSPBUFL, KVSETSC=YLFLDSC_SP%IVSET,&
                      & KRESOL=YDGEOMETRY%YRDIM%NRESOL,&
                      & KPROMA=YDGEOMETRY%YRGEM%NGPTOT, PGP=ZGPBUFL)

  ELSEIF (LLUV) THEN

    CALL DIR_TRANS (PSPVOR=YDSP%VOR, PSPDIV=YDSP%DIV, KVSETUV=YDGEOMETRY%YRMP%NBSETLEV,&
                  & PSPSCALAR=ZSPBUFL, KVSETSC=YLFLDSC_SP%IVSET,&
                  & KRESOL=YDGEOMETRY%YRDIM%NRESOL,&
                  & KPROMA=YDGEOMETRY%YRGEM%NGPTOT, PGP=ZGPBUFL)

  ELSE

    CALL DIR_TRANS (PSPSCALAR=ZSPBUFL, KVSETSC=YLFLDSC_SP%IVSET,&
                  & KRESOL=YDGEOMETRY%YRDIM%NRESOL,&
                  & KPROMA=YDGEOMETRY%YRGEM%NGPTOT, PGP=ZGPBUFL)

  ENDIF

! Fill model variables
  CALL IOSPECA_SELECTF(YDGEOMETRY,YDML_GCONF%YGFL,YDDYN, NSPECACT_READ_GP, YLIOSP, IFLDSPL,&
                      & IVSETOFF (MYSETV), 0_JPIM, YLFLDSC_SP, ZSPBUFL, YDSP,YDMCUF=YDMCUF)


  DEALLOCATE (ZGPBUFL, ZSPBUFL, YLFLDSC_SP, YLFLDSC_GP)


  CALL IOSPECA_FINISH (YDGEOMETRY, YDDYNA,NSPECACT_READ_GP, YLIOSP, YDSP)

ENDIF

CALL SUSPECA_FIXUP(YDGEOMETRY,YDML_GCONF%YGFL,KFILE,YDSP%GFL,YDSP%SP)


IF (LHOOK) CALL DR_HOOK ('SUSPECA_GP',1,ZHOOK_HANDLE)

END SUBROUTINE SUSPECA_GP

