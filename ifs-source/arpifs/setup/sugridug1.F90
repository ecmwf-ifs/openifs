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

SUBROUTINE SUGRIDUG1(YDGEOMETRY,YDGFL,YGFL)

!**** *SUGRIDUG1*  - Initialize the upper air grid-point fields using
!                    artificial data

!     Purpose.
!     --------
!           Initialize the upper air gridpoint fields of the model
!           using artificial data

!**   Interface.
!     ----------
!        *CALL* *SUGRIDUG1

!        Explicit arguments :
!        ------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      George Mozdzynski *ECMWF*
!      Original : 97-09-22 

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOM_YGFL , ONLY : TYPE_GFLD

IMPLICIT NONE

TYPE(GEOMETRY) , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)     , INTENT(INOUT) :: YDGFL
TYPE(TYPE_GFLD) ,INTENT(INOUT) :: YGFL

REAL(KIND=JPRB), ALLOCATABLE :: ZZQT0(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPQ(:)

INTEGER(KIND=JPIM) :: ICEND, IENDC, IM, IPCT, IPT, ISP, ISTC, J, JLEV, JMLOC, JN, JSTGLO&
 & , IBL, IFLD,JGFL  

LOGICAL :: LLINIT0
REAL(KIND=JPRB) :: ZQMAX, ZQMIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gpnorm3.intfb.h"
#include "speree.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGRIDUG1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & YQ=>YGFL%YQ, &
 & NPROMA=>YDDIM%NPROMA, NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, &
 & NUMP=>YDDIM%NUMP, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GFL=>YDGFL%GFL, &
 & MYMS=>YDLAP%MYMS, NASM0=>YDLAP%NASM0)
!     ------------------------------------------------------------------

LLINIT0=.TRUE.
IPT=0
! Set Gridpoint part of GFL to ZERO
DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LGP) THEN
    GFL(:,:,JGFL,:) =0.0_JPRB
  ENDIF
ENDDO

IF( YQ%LGP )THEN
  ALLOCATE (ZZQT0(NGPTOT,NFLEVG))
  ALLOCATE (ZSPQ(NSPEC2))
  ZSPQ(:)=0.0_JPRB
  ZZQT0(:,:)=0.0_JPRB
  ZQMIN= 0.25E-5_JPRB
  ZQMAX= 0.01_JPRB
  DO JLEV=1,NFLEVG
    DO JMLOC=1,NUMP
      IM = MYMS(JMLOC)
      DO JN=IM,NSMAX
        ISP = NASM0(IM)+(JN-IM)*2
        ZSPQ(ISP)=1.E-09_JPRB * IM * JN * JLEV&
         & / ( REAL(NSMAX*NSMAX,JPRB) * REAL(NFLEVG,JPRB) )  
        IF (IM  /=  0) THEN
          ZSPQ(ISP+1)=0.0_JPRB
        ENDIF
      ENDDO
    ENDDO
    IF (NUMP > 0.AND. MYMS(1) == 0) THEN
      ZSPQ(NASM0(0))=ZQMIN+ (ZQMAX-ZQMIN)*&
       & MAX(0.0_JPRB,REAL(2*JLEV,JPRB)/REAL(NFLEVG,JPRB)-1.0_JPRB)  
    ENDIF
    CALL SPEREE(YDGEOMETRY,1,1,ZSPQ,ZZQT0(1,JLEV))
  ENDDO
  DEALLOCATE(ZSPQ)
ENDIF

!     ------------------------------------------------------------------

!*       ADD TO UPPER AIR GRIDPOINT WORKFILE.
!        ------------------------------------

IF(YQ%LGP) THEN
  DO JSTGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
    IPCT=1
    ISTC=1
    IENDC=ICEND
    IBL=(JSTGLO-1)/NPROMA+1
    IFLD=0
    DO JLEV=1,NFLEVG
      DO J=ISTC,IENDC
        GFL(J,JLEV,YQ%MP,IBL)= ZZQT0(JSTGLO+J-1,JLEV)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF (ALLOCATED(ZZQT0)) DEALLOCATE (ZZQT0)

WRITE(NULOUT,'(A,A)') ' SUGRIDUG1 : STATISTICS FOR ALL ','GFL FIELDS'
DO JGFL=1,NUMFLDS
  WRITE(NULOUT,'(A)') YCOMP(JGFL)%CNAME
  CALL GPNORM3(YDGEOMETRY,GFL,NDIM,JGFL)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGRIDUG1',1,ZHOOK_HANDLE)
END SUBROUTINE SUGRIDUG1
