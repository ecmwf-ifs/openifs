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

SUBROUTINE SUOROG (YDGEOMETRY, PSPOR)

!**** *SUOROG*  - Routine to initialize the grid-point orography fields.

!     Purpose.
!     --------
!           Initialize the grid-point orography fields of the model.

!**   Interface.
!     ----------
!        *CALL* *SUOROG(.....)*

!        Explicit arguments :
!        --------------------
!          INPUT:
!            PSPOR   - spectral orography field

!        Implicit arguments :
!        --------------------
!        The grid-point fields of the model.
!        The boundary condition fields of the model.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      David Dent *ECMWF*
!      Original : 92-05-19

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      R. El Khatib 14-May-2018 move allocations to sugeometry and merge initializations with sueorog
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LELAM
#ifdef WITH_ATLAS
USE ATLAS_MODULE, ONLY : ATLAS_REAL, ATLAS_FIELD, ATLAS_MESH_NODES, &
 & ATLAS_FUNCTIONSPACE_NODECOLUMNS
#endif
USE TYPE_GEOMETRY , ONLY : GEOMETRY
USE YOMDYNA_STATIC  , ONLY : LGRADGP

IMPLICIT NONE

TYPE(GEOMETRY)  , INTENT(INOUT)    :: YDGEOMETRY
REAL(KIND=JPRB) , INTENT(INOUT) :: PSPOR(YDGEOMETRY%YRDIM%NSPEC2) 

REAL(KIND=JPRB) , ALLOCATABLE :: ZU(:), ZSWORK(:)

REAL(KIND=JPRB) :: ZDIV      (YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZZOROG    (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGL   (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGM   (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGLL  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGMM  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGLM  (YDGEOMETRY%YRGEM%NGPTOT)

#ifdef WITH_ATLAS
INTEGER(KIND=JPIM),ALLOCATABLE :: IPA(:)
REAL(KIND=JPRB), POINTER :: ZVAR(:,:)
REAL(KIND=JPRB), POINTER :: ZGRAD(:,:,:)
LOGICAL   , POINTER :: LLIS_GHOST(:)                                  
TYPE(ATLAS_FIELD) :: GHOSTFIELD
TYPE(ATLAS_FIELD) :: VARFIELD
TYPE(ATLAS_FIELD) :: GRADFIELD
TYPE(ATLAS_MESH_NODES) :: NODES
TYPE(ATLAS_FUNCTIONSPACE_NODECOLUMNS) :: NODE_COLUMNS
#endif
INTEGER(KIND=JPIM) :: JL,IP,INODE_SIZE
INTEGER(KIND=JPIM) :: JKGLO, IEND, IBL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fv_gradient.intfb.h"
#include "reespe.intfb.h"
#include "speree.intfb.h"
#include "speuv.intfb.h"
#include "sueorog.intfb.h"

!     ------------------------------------------------------------------

!*       1.    INITIALIZE GRIDPOINT OROGRAPHY.
!              -------------------------------

IF (LHOOK) CALL DR_HOOK('SUOROG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NSPEC2=>YDDIM%NSPEC2, &
 & NGPTOT=>YDGEM%NGPTOT)


IF (LELAM) THEN

!*       1.1   PLANE GEOMETRY.
  
  CALL SUEOROG (YDGEOMETRY, PSPOR, ZZOROG, ZZOROGL, ZZOROGM, ZZOROGLL, ZZOROGLM, ZZOROGMM)

ELSE

!*       1.2   SPHERICAL GEOMETRY.

!*       1.2.1   ALLOCATIONS.


!*       1.2.2   INITIALIZE GRIDPOINT OROGRAPHY

  CALL SPEREE(YDGEOMETRY,1,1,PSPOR,ZZOROG)

!*       1.2.3   COMPUTE FIRST ORDER DERIVATIVES OF OROGRAPHY.

  IF(.NOT.LGRADGP) THEN

    ZDIV(:)=0.0_JPRB
    CALL SPEUV(YDGEOMETRY,ZDIV,PSPOR,ZZOROGL,ZZOROGM,1,1,2)
    
  ELSE
#ifdef WITH_ATLAS
    NODE_COLUMNS = YDGEOMETRY%YRATLAS%FVM%NODE_COLUMNS()
    NODES        = YDGEOMETRY%YRATLAS%MESH%NODES()
    INODE_SIZE   = NODES%SIZE()
    GHOSTFIELD   = NODES%GHOST()
    CALL GHOSTFIELD%DATA(LLIS_GHOST)
    ALLOCATE(IPA(INODE_SIZE))
    IP = 0
    DO JL=1,INODE_SIZE
      IF(.NOT. LLIS_GHOST(JL)) THEN
        IP = IP+1
        IPA(JL) = IP
      ENDIF
    ENDDO
    VARFIELD  = NODE_COLUMNS%CREATE_FIELD(NAME="var_orog", KIND=ATLAS_REAL(JPRB),LEVELS=1)
    GRADFIELD = NODE_COLUMNS%CREATE_FIELD(NAME="grad_orog",KIND=ATLAS_REAL(JPRB),LEVELS=1,VARS=[2])
    CALL VARFIELD%DATA(ZVAR)
    CALL GRADFIELD%DATA(ZGRAD)
    DO JL=1,INODE_SIZE
      IF(.NOT. LLIS_GHOST(JL)) THEN
        IP = IPA(JL)
        ZVAR(1,JL) = ZZOROG(IP)
      ENDIF
    ENDDO
    CALL NODE_COLUMNS%HALO_EXCHANGE(VARFIELD)
    CALL FV_GRADIENT(YDGEOMETRY,ZVAR,ZGRAD)
    DO JL=1,INODE_SIZE
      IF(.NOT. LLIS_GHOST(JL)) THEN
        IP = IPA(JL)
        ZZOROGL(IP) = ZGRAD(1,1,JL)
        ZZOROGM(IP) = ZGRAD(2,1,JL)
      ENDIF
    ENDDO
    CALL VARFIELD%FINAL()
    CALL GRADFIELD%FINAL()
    CALL GHOSTFIELD%FINAL()
    CALL NODES%FINAL()
    CALL NODE_COLUMNS%FINAL()
#endif
  ENDIF

!*       1.2.4   COMPUTE SOME SECOND ORDER DERIVATIVES OF OROGRAPHY FOR NH.

  IF(YDGEOMETRY%LNONHYD_GEOM) THEN
    ALLOCATE(ZSWORK(NSPEC2))
    ALLOCATE(ZU(NGPTOT))
!       * orogll, oroglm.
    CALL REESPE(YDGEOMETRY,1,1,ZSWORK,ZZOROGL)
    CALL SPEUV(YDGEOMETRY,ZDIV,ZSWORK,ZZOROGLL,ZZOROGLM,1,1,2)
!       * orogmm:
    CALL REESPE(YDGEOMETRY,1,1,ZSWORK,ZZOROGM)
    CALL SPEUV(YDGEOMETRY,ZDIV,ZSWORK,ZU,ZZOROGMM,1,1,2)
    DEALLOCATE(ZU)
    DEALLOCATE(ZSWORK)
  ENDIF

ENDIF

!*       1.3    FILL YROROG
!               -----------

! Allocate and copy native NPROMA orography data structure
IF (.NOT. ALLOCATED(YDGEOMETRY%YROROG_B%OROG))  ALLOCATE(YDGEOMETRY%YROROG_B%OROG (NPROMA,NGPBLKS))
IF (.NOT. ALLOCATED(YDGEOMETRY%YROROG_B%OROGL)) ALLOCATE(YDGEOMETRY%YROROG_B%OROGL(NPROMA,NGPBLKS))
IF (.NOT. ALLOCATED(YDGEOMETRY%YROROG_B%OROGM)) ALLOCATE(YDGEOMETRY%YROROG_B%OROGM(NPROMA,NGPBLKS))

! * orography and horizontal gradient of orography:
DO JKGLO=1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1
  YDGEOMETRY%YROROG(IBL)%OROG (1:IEND)=ZZOROG (JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG(IBL)%OROGL(1:IEND)=ZZOROGL(JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG(IBL)%OROGM(1:IEND)=ZZOROGM(JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG_B%OROG (1:IEND,IBL)=ZZOROG (JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG_B%OROGL(1:IEND,IBL)=ZZOROGL(JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG_B%OROGM(1:IEND,IBL)=ZZOROGM(JKGLO:JKGLO+IEND-1)
ENDDO
  
! * horizontal second-order derivatives of orography:
IF (YDGEOMETRY%LNONHYD_GEOM) THEN
  DO JKGLO=1,NGPTOT,NPROMA
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    YDGEOMETRY%YROROG(IBL)%OROGLL(1:IEND)=ZZOROGLL(JKGLO:JKGLO+IEND-1)
    YDGEOMETRY%YROROG(IBL)%OROGLM(1:IEND)=ZZOROGLM(JKGLO:JKGLO+IEND-1)
    YDGEOMETRY%YROROG(IBL)%OROGMM(1:IEND)=ZZOROGMM(JKGLO:JKGLO+IEND-1)
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUOROG',1,ZHOOK_HANDLE)
END SUBROUTINE SUOROG
