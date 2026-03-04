! (C) Copyright 2014- ECMWF.
! (C) Copyright 2014- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE CDMAP_MOD
CONTAINS
SUBROUTINE CDMAP(KM,KMLOC,KSL,KSLO,PEPSNM, KDIR, KDGNH, KDGNHD,&
& KFIELDS, PCOEFA, PCOEFS)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_FLT
USE TPM_GEOMETRY
USE TPM_DISTR       ,ONLY : D
USE TPM_TRANS       ,ONLY : FOUBUF_IN, FOUBUF
USE SEEFMM_MIX

!**** *CDMAP* - REMAP ROOTS
!
!     Purpose.
!     --------
! remap from one set of roots to another using Christoffel-Darboux formula, see Chien + Alpert, 1997.

!**   Interface.
!     ----------
!        *CALL* *CDMAP(...)

!        Explicit arguments :
!        --------------------
!          KM        - zonal wavenumber
!          KMLOC     - local zonal wavenumber
!
!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!      Chien + Alpert, 1997.

!     Author.
!     -------
!        Nils Wedi  *ECMWF*

!     Modifications.
!     --------------
!        Original : 14-05-14 
!     ------------------------------------------------------------------

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KSL
INTEGER(KIND=JPIM), INTENT(IN) :: KSLO
REAL(KIND=JPRB), INTENT(IN) :: PEPSNM
INTEGER(KIND=JPIM), INTENT(IN) :: KDIR ! direction of map
INTEGER(KIND=JPIM), INTENT(IN) :: KDGNH
INTEGER(KIND=JPIM), INTENT(IN) :: KDGNHD
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
REAL(KIND=JPRB), INTENT(INOUT) :: PCOEFA(:,:)
REAL(KIND=JPRB), INTENT(INOUT) :: PCOEFS(:,:)

INTEGER(KIND=JPIM) :: JGL, IGL, JF
REAL(KIND=JPRB), ALLOCATABLE :: ZALL(:,:), ZQX(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZALL1(:,:), ZQY(:,:)
INTEGER(KIND=JPIM) :: ISTN(KDGNH), ISTS(KDGNH)

INTEGER(KIND=JPIM) :: IGLS, IPROC, IPROCS, IEND, IENDO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('CDMAP_MOD',0,ZHOOK_HANDLE)

IF( KDIR == -1 ) THEN
  ! inverse map from internal (gg) roots to post-processing roots

  IENDO = 2*KDGNHD -  KSLO + 1
  IEND = 2*KDGNH -  KSL + 1

  !!!!! fourier buffer setup in output latitudes, may not work if different from input !!!!
  DO IGL=KSLO, KDGNHD
    IPROC = D%NPROCL(IGL)
    ISTN(IGL) = (D%NSTAGT0B(IPROC) + D%NPNTGTB1(KMLOC,IGL))*KFIELDS
    IGLS = 2*KDGNH+1-IGL
    IPROCS = D%NPROCL(IGLS)
    ISTS(IGL) = (D%NSTAGT0B(IPROCS) + D%NPNTGTB1(KMLOC,IGLS))*KFIELDS
  ENDDO

  ALLOCATE(ZALL(KFIELDS, 2*KDGNHD))
  ALLOCATE(ZALL1(KFIELDS, 2*KDGNHD))
  ALLOCATE(ZQX(KFIELDS, 2*KDGNH))
  ALLOCATE(ZQY(KFIELDS, 2*KDGNH))
  ZQX(:,1:KSL) = 0._JPRB
  ZQX(:,IEND:2*KDGNH) = 0._JPRB
  ZQY(:,1:KSL) = 0._JPRB
  ZQY(:,IEND:2*KDGNH) = 0._JPRB
  DO JGL=KSL, IEND
    ZQX(1:KFIELDS,JGL)=S%FA(KMLOC)%RPNMWI(JGL-KSL+1,1)*PCOEFA(1:KFIELDS,JGL)
    ZQY(1:KFIELDS,JGL)=S%FA(KMLOC)%RPNMWI(JGL-KSL+1,2)*PCOEFA(1:KFIELDS,JGL)
  ENDDO
  CALL SEEFMM_MULM(S%FMM_INTI,KFIELDS,1_JPIM,.TRUE.,ZQX,ZALL1)
  CALL SEEFMM_MULM(S%FMM_INTI,KFIELDS,1_JPIM,.TRUE.,ZQY,ZALL)
  DEALLOCATE(ZQX)
  DEALLOCATE(ZQY)
  ! minus sign comes from FMM ?!
  ! fill buffer
  DO IGL=KSLO,KDGNHD
    IGLS = 2*KDGNHD+1-IGL
    DO JF=1,KFIELDS
      FOUBUF_IN(ISTN(IGL)+JF) = S%FA(KMLOC)%RPNMWO(IGL-KSLO+1,1)*ZALL1(JF,IGL) & 
       & - S%FA(KMLOC)%RPNMWO(IGL-KSLO+1,2)*ZALL(JF,IGL)
      FOUBUF_IN(ISTS(IGL)+JF) = S%FA(KMLOC)%RPNMWO(IGLS-KSLO+1,1)*ZALL1(JF,IGLS) & 
       & - S%FA(KMLOC)%RPNMWO(IGLS-KSLO+1,2)*ZALL(JF,IGLS)
    ENDDO
  ENDDO
  DEALLOCATE(ZALL1)
  DEALLOCATE(ZALL)

ELSE
! direct map from post-processing/input field roots to internal (gg) roots
! this assumes essentially a nearest neighbour interpolation in latitude
! a more accurate approach may be 
! a local gridpoint interpolation of the input field to the target latitudes prior to the transforms

  IENDO = 2*KDGNHD -  KSLO + 1
  IEND   = 2*KDGNH -  KSL + 1

  !!!!! fourier buffer setup in input data latitudes, may not work if different from output !!!!
  DO JGL=KSLO, KDGNHD
    IPROC = D%NPROCL(JGL)
    ISTN(JGL) = (D%NSTAGT1B(IPROC) + D%NPNTGTB1(KMLOC,JGL))*KFIELDS
    IGLS = 2*KDGNHD+1-JGL
    IPROCS = D%NPROCL(IGLS)
    ISTS(JGL) = (D%NSTAGT1B(IPROCS) + D%NPNTGTB1(KMLOC,IGLS))*KFIELDS
  ENDDO

  ALLOCATE( ZQX( KFIELDS, 2*KDGNHD))
  ZQX(:,1:KSLO) = 0._JPRB
  ZQX(:,IENDO:2*KDGNHD) = 0._JPRB
  DO JGL=KSLO, KDGNHD
    IGLS = 2*KDGNHD+1-JGL
    DO JF=1,KFIELDS
      ZQX(JF,JGL)=FOUBUF(ISTN(JGL)+JF)
      ZQX(JF,IGLS)=FOUBUF(ISTS(JGL)+JF)
    ENDDO
  ENDDO

  ! split into symmetric / antisymmetric
  DO IGL=KSL,KDGNH
    IGLS = 2*KDGNH+1-IGL
    PCOEFS(1:KFIELDS,IGL) = ZQX(1:KFIELDS,IGL) + ZQX(1:KFIELDS,IGLS)
    PCOEFA(1:KFIELDS,IGL) = ZQX(1:KFIELDS,IGL) - ZQX(1:KFIELDS,IGLS)
  ENDDO

  DEALLOCATE(ZQX)
  
ENDIF

IF (LHOOK) CALL DR_HOOK('CDMAP_MOD',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE CDMAP
END MODULE CDMAP_MOD
