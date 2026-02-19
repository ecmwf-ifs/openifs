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

SUBROUTINE RELAXGP(YDGEOMETRY,YDGFL,YDGMV,YGFL,YDSPEC)

!**** *RELAXGP * - Perform relaxation in grib point space

!     Purpose.  
!     --------

!**   Interface.  CALL RELAXGP
!     ---------- 

!     Explicit arguments : none
!     --------------------

!     Externals.
!     ----------
!     INV_TRANS

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Linus Magnusson  *ECMWF*
!        The subroutine was splitted from TRANSINH
!        Original: march 2013
!
!     Modifications.
!     --------------
!      T. Stockdale (Nov 2015) Control on which spectral wavenumbers relaxed
!      J. Flemming (Jan 2014) Relaxation for ozone and fix for surface pressure relaxation
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL   , ONLY : TGFL
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCT3   , ONLY : NSTEP
USE YOMRLX   , ONLY : LRLXG, NFRLXG, LRLXVO, LRLXDI, LRLXTE, LRLXQ, LRLXLP,&
  & LRLXQL, LRLXQI, LRLXQC, LRLXO3, XRLXVO, XRLXTE, XRLXQ, XRLXLP, XRLXO3, ALATRLX1, ALATRLX2,&
  & ALONRLX1, ALONRLX2, NRLXLMIN, NRLXLMAX, NRLXLMINU, NRLXLMAXU, AXRLX, AYRLX, AZRLX, NRLXSMAX
USE YOMSRLX  , ONLY : TRLXVO, TRLXDI, TRLXTE, TRLXQ, TRLXLP,&
 & XPRLXG, TRLXQI, TRLXQL, TRLXQC, TRLXO3
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RPI
USE SPECTRAL_FIELDS_MOD,  ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE YOMMP0   , ONLY : MYSETV

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV
TYPE(TYPE_GFLD)     ,INTENT(INOUT) :: YGFL
TYPE(SPECTRAL_FIELD),INTENT(IN)    :: YDSPEC
INTEGER(KIND=JPIM) :: JSTEP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZPRLXVO(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2),&
  & ZPRLXDI(YDGEOMETRY%YRDIMV%NFLSUR,YDGEOMETRY%YRDIM%NSPEC2) ,&
  & ZPRLXTE(YDGEOMETRY%YRDIMV%NFLSUR,YDGEOMETRY%YRDIM%NSPEC2,1), ZPRLXLP(1,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZUSRLX

REAL(KIND=JPRB) :: ZGPRLXUV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,4,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB) :: ZGPRLXTE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,1,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB) :: ZGPRLXLP(YDGEOMETRY%YRDIM%NPROMA,1,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB) :: ZMASK
REAL(KIND=JPRB) :: ZALPHA
REAL(KIND=JPRB) :: ZPRLXQ,ZPRLXTV,ZPRLXO3
REAL(KIND=JPRB) :: ZAX, ZAY, ZAZ
REAL(KIND=JPRB) :: ZFLON1,ZFLON2,ZFLAT1,ZFLAT2
REAL(KIND=JPRB) :: ZDELTA_MIN,ZDELTA_MAX,ZDELTA_MINU,ZDELTA_MAXU
REAL(KIND=JPRB) :: ZFAC1,ZFAC2,ZFAC1U,ZFAC2U
REAL(KIND=JPRB) :: ZRR
LOGICAL :: LLSCDERS,LLVORGP,LLDIVGP,LLUVDER
INTEGER(KIND=JPIM) :: IVSETSC(1)
INTEGER(KIND=JPIM) :: JKGLO, ICEND, IBL, IOFF, JLEV, JROF
INTEGER(KIND=JPIM) :: JMLOC, IM, JN, IJSE

!     ------------------------------------------------------------------

#include "inv_trans.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RELAXGP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDMP=>YDGEOMETRY%YRMP, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB)
ASSOCIATE(YA=>YGFL%YA, YI=>YGFL%YI, YL=>YGFL%YL, YO3=>YGFL%YO3, YQ=>YGFL%YQ, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, &
 & NSPEC2=>YDDIM%NSPEC2, NUMP=>YDDIM%NUMP, NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, NFLSUR=>YDDIMV%NFLSUR, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GFL=>YDGFL%GFL, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, YT0=>YDGMV%YT0, &
 & MYMS=>YDLAP%MYMS, &
 & NBSETLEV=>YDMP%NBSETLEV, NBSETSP=>YDMP%NBSETSP)
!     ------------------------------------------------------------------


ZRR=1.61_JPRB

! Should we use TRANS_INQ to define MYMS, instead of taking values from module??


! Relaxation
IF (LRLXG .AND. NSTEP > 1) THEN
  WRITE(NULOUT,*)' Relaxation in relaxgp at step ',NSTEP,XPRLXG(1),XPRLXG(2)

  ! Time interpolation for spectral field before transform
  ZPRLXVO(:,:)=0.0_JPRB
  ZPRLXDI(:,:)=0.0_JPRB
  ZPRLXTE(:,:,:)=0.0_JPRB
  ZPRLXLP(:,:)=0.0_JPRB

  ! reference temperature
  IF (LRLXTE) THEN
    DO JSTEP = 1, NFRLXG
      ZPRLXTE(:,:,1)=ZPRLXTE(:,:,1)+TRLXTE(:,:,JSTEP)*XPRLXG(JSTEP)
    ENDDO
    DO JMLOC=1,NUMP
      IM=MYMS(JMLOC)
      DO JN=IM,NSMAX
        IF(JN>NRLXSMAX) THEN
          ! reset relaxation array to existing model value
          IJSE=YDSPEC%NASM0(IM)+2*(JN-IM)
          IF(IJSE>NSPEC2) CALL ABOR1('RELAXATION OOB')
          ZPRLXTE(:,IJSE,1)=YDSPEC%T(:,IJSE)
          ZPRLXTE(:,IJSE+1,1)=YDSPEC%T(:,IJSE+1)
        ENDIF
      ENDDO
    ENDDO
  ELSE
    ZPRLXTE(:,:,1)=YDSPEC%T(:,:)
  ENDIF
  ! reference vorticity
  IF (LRLXVO) THEN
    DO JSTEP = 1, NFRLXG
      ZPRLXVO(:,:)=ZPRLXVO(:,:)+TRLXVO(:,:,JSTEP)*XPRLXG(JSTEP)
    ENDDO
    DO JMLOC=1,NUMP
      IM=MYMS(JMLOC)
      DO JN=IM,NSMAX
        IF(JN>NRLXSMAX) THEN
          IJSE=YDSPEC%NASM0(IM)+2*(JN-IM)
          ZPRLXVO(:,IJSE)=YDSPEC%VOR(:,IJSE)
          ZPRLXVO(:,IJSE+1)=YDSPEC%VOR(:,IJSE+1)
        ENDIF
      ENDDO
    ENDDO
  ELSE
    ZPRLXVO(:,:)=YDSPEC%VOR(:,:)
  ENDIF
  ! reference divergence
  IF (LRLXDI) THEN
    DO JSTEP = 1, NFRLXG
      ZPRLXDI(:,:)=ZPRLXDI(:,:)+TRLXDI(:,:,JSTEP)*XPRLXG(JSTEP)
    ENDDO
    DO JMLOC=1,NUMP
      IM=MYMS(JMLOC)
      DO JN=IM,NSMAX
        IF(JN>NRLXSMAX) THEN
          IJSE=YDSPEC%NASM0(IM)+2*(JN-IM)
          ZPRLXDI(:,IJSE)=YDSPEC%DIV(:,IJSE)
          ZPRLXDI(:,IJSE+1)=YDSPEC%DIV(:,IJSE+1)
        ENDIF
      ENDDO
    ENDDO
  ELSE
    ZPRLXDI(:,:)=YDSPEC%DIV(:,:)
  ENDIF

  ! reference surface pressure
  IF ( MYSETV == NBSETSP ) THEN 
    IF (LRLXLP) THEN
      DO JSTEP = 1, NFRLXG
        ZPRLXLP(1,:)=ZPRLXLP(1,:)+TRLXLP(:,JSTEP)*XPRLXG(JSTEP)
      ENDDO
    ELSE
      ZPRLXLP(1,:)=YDSPEC%SP(:)
    ENDIF
  ENDIF

  LLSCDERS=.FALSE.
  LLVORGP=.FALSE.
  LLDIVGP=.FALSE.
  LLUVDER=.TRUE.
  IVSETSC(1) = NBSETSP

  CALL INV_TRANS(PSPVOR=ZPRLXVO,PSPDIV=ZPRLXDI,PSPSC2=ZPRLXLP,&
    & PSPSC3A=ZPRLXTE,&
    & LDSCDERS=LLSCDERS,LDVORGP=LLVORGP,LDDIVGP=LLDIVGP,LDUVDER=LLUVDER,&
    & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
    & KVSETSC3A=NBSETLEV,&
    & PGPUV=ZGPRLXUV,PGP2=ZGPRLXLP,PGP3A=ZGPRLXTE)

  ! actual relaxation in grid point space
  DO JKGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    IOFF=JKGLO

    ! smoothing parameters
    ZAX=AXRLX*180._JPRB/RPI
    ZAY=AYRLX*180._JPRB/RPI
    ZAZ=AZRLX

    ! upper air variables first
    DO JLEV=1,NFLEVG
      DO JROF=1,ICEND
        !Construct the mask
        ! longitudes
        IF (ALONRLX1*180._JPRB/RPI <=0.1 .AND.&
          & ALONRLX2*180._JPRB/RPI >= 359.9) THEN
 
         ZFLON1=1._JPRB
         ZFLON2=1._JPRB
       ELSE
         ZFLON1=1.0_JPRB/(1.0_JPRB+EXP(ZAX*(YDGSGEOM_NB%GELAM(IOFF+JROF-1)-ALONRLX1)))
         ZFLON2=1.0_JPRB-1.0_JPRB/(1.0_JPRB+EXP(ZAX*(YDGSGEOM_NB%GELAM(IOFF+JROF-1)-ALONRLX2)))
       ENDIF
       ! latitudes
       IF (ALATRLX1*180._JPRB/RPI >= 89.9) THEN
         ! here we go all the way to the north pole
         ZFLAT1=1._JPRB
       ELSE
         ZFLAT1=1.0_JPRB-1.0_JPRB/(1.0_JPRB+EXP(ZAY*(YDGSGEOM_NB%GELAT(IOFF+JROF-1)-ALATRLX1)))
       ENDIF
       IF (ALATRLX2*180._JPRB/RPI <= -89.9) THEN
         ! here we go all the way to the south pole
         ZFLAT2=1._JPRB
       ELSE
         ZFLAT2=1.0_JPRB/(1.0_JPRB+EXP(ZAY*(YDGSGEOM_NB%GELAT(IOFF+JROF-1)-ALATRLX2)))
       ENDIF
       ! contruct final mask
       IF(ALONRLX1 < ALONRLX2) THEN
         ZMASK=ZFLON1*ZFLON2*ZFLAT1*ZFLAT2
       ELSE
         ZMASK=(ZFLON1+ZFLON2)*ZFLAT1*ZFLAT2
       ENDIF

       ! the relaxation coefficient is changed by the mask
       ZDELTA_MIN=FLOAT(JLEV-NRLXLMIN)
       ZDELTA_MAX=FLOAT(NRLXLMAX-JLEV)
       IF (NRLXLMIN > 1.AND.NRLXLMAX < NFLEVG) THEN
         ZFAC1=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MIN))
         ZFAC2=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MAX))
       ELSEIF (NRLXLMIN == 1.AND.NRLXLMAX < NFLEVG) THEN
         ZFAC1=1.0
         ZFAC2=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MAX))
       ELSEIF (NRLXLMIN > 1.AND.NRLXLMAX == NFLEVG) THEN
         ZFAC1=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MIN))
         ZFAC2=1.0
       ELSE
         ZFAC1=1.0
         ZFAC2=1.0
       ENDIF

       ZDELTA_MINU=FLOAT(JLEV-NRLXLMINU)
       ZDELTA_MAXU=FLOAT(NRLXLMAXU-JLEV)
       IF (NRLXLMINU > 1.AND.NRLXLMAXU < NFLEVG) THEN
         ZFAC1U=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MINU))
         ZFAC2U=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MAXU))
       ELSEIF (NRLXLMINU == 1.AND.NRLXLMAXU < NFLEVG) THEN
         ZFAC1U=1.0
         ZFAC2U=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MAXU))
       ELSEIF (NRLXLMINU > 1.AND.NRLXLMAXU == NFLEVG) THEN
         ZFAC1U=1.0/(1.0+EXP(-1.0*ZAZ*ZDELTA_MINU))
         ZFAC2U=1.0
       ELSE
         ZFAC1U=1.0
         ZFAC2U=1.0
       ENDIF

       !temperature relaxation
       ZALPHA=XRLXTE*ZMASK*ZFAC1*ZFAC2
       ZUSRLX=1.0_JPRB/(1.0_JPRB+ZALPHA)
       IF (LRLXTE) THEN
         !ZPRLXQ=0.0_JPRB
         !DO JSTEP = 1, NFRLXG
         !  ZPRLXQ=ZPRLXQ+TRLXQ(JROF,JLEV,JSTEP,IBL)*XPRLXG(JSTEP)
         !ENDDO
         !ZPRLXTV=ZGZPRLXTE(JROF,JLEV,1,IBL)*(1_JPRB+(ZRR-1)*ZPRLXQ)
         ZPRLXTV=ZGPRLXTE(JROF,JLEV,1,IBL)*(1_JPRB+(ZRR-1)*GFL(JROF,JLEV,YQ%MP,IBL))
         GMV(JROF,JLEV,YT0%MT,IBL)=(GMV(JROF,JLEV,YT0%MT,IBL) +&
           & ZALPHA*ZPRLXTV)*ZUSRLX
       ENDIF
       
       !  u,v relaxation
       ZALPHA=XRLXVO*ZMASK*ZFAC1U*ZFAC2U
       ZUSRLX=1.0_JPRB/(1.0_JPRB+ZALPHA)
       IF (LRLXVO .OR. LRLXDI) THEN
         GMV(JROF,JLEV,YT0%MU,IBL)=(GMV(JROF,JLEV,YT0%MU,IBL) +&
           & ZALPHA*ZGPRLXUV(JROF,JLEV,1,IBL))*ZUSRLX
         GMV(JROF,JLEV,YT0%MV,IBL)=(GMV(JROF,JLEV,YT0%MV,IBL) +&
           & ZALPHA*ZGPRLXUV(JROF,JLEV,2,IBL))*ZUSRLX
       ENDIF

       ! Humidity related fields - all using the same time scale
       ZALPHA=XRLXQ*ZMASK*ZFAC1*ZFAC2
       ZUSRLX=1.0_JPRB/(1.0_JPRB+ZALPHA)
       !Specific humidity
       IF (LRLXQ) THEN
         ZPRLXQ=0.0_JPRB
         DO JSTEP = 1, NFRLXG
           ZPRLXQ=ZPRLXQ+TRLXQ(JROF,JLEV,JSTEP,IBL)*XPRLXG(JSTEP)
         ENDDO
         GFL(JROF,JLEV,YQ%MP,IBL)=(GFL(JROF,JLEV,YQ%MP,IBL) +&
            & ZALPHA*ZPRLXQ)*ZUSRLX
       ENDIF
       !Specific cloud liquid water content
       IF (LRLXQL) THEN
         ZPRLXQ=0.0_JPRB
         DO JSTEP = 1, NFRLXG
           ZPRLXQ=ZPRLXQ+TRLXQL(JROF,JLEV,JSTEP,IBL)*XPRLXG(JSTEP)
         ENDDO
         GFL(JROF,JLEV,YL%MP,IBL)=(GFL(JROF,JLEV,YL%MP,IBL) +&
           & ZALPHA*ZPRLXQ)*ZUSRLX
       ENDIF
       !Specific cloud ice water content
       IF (LRLXQI) THEN
         ZPRLXQ=0.0_JPRB
         DO JSTEP = 1, NFRLXG
           ZPRLXQ=ZPRLXQ+TRLXQI(JROF,JLEV,JSTEP,IBL)*XPRLXG(JSTEP)
         ENDDO
         GFL(JROF,JLEV,YI%MP,IBL)=(GFL(JROF,JLEV,YI%MP,IBL) +&
           & ZALPHA*ZPRLXQ)*ZUSRLX
       ENDIF
       !Cloud cover
       IF (LRLXQC) THEN
         ZPRLXQ=0.0_JPRB
         DO JSTEP = 1, NFRLXG
           ZPRLXQ=ZPRLXQ+TRLXQC(JROF,JLEV,JSTEP,IBL)*XPRLXG(JSTEP)
         ENDDO
         GFL(JROF,JLEV,YA%MP,IBL)=(GFL(JROF,JLEV,YA%MP,IBL) +&
           & ZALPHA*ZPRLXQ)*ZUSRLX
       ENDIF
       !ozone
       IF (LRLXO3) THEN
         ZALPHA=XRLXO3*ZMASK*ZFAC1*ZFAC2
         ZUSRLX=1.0_JPRB/(1.0_JPRB+ZALPHA)
         ZPRLXO3=0.0_JPRB
         DO JSTEP = 1, NFRLXG
           ZPRLXO3=ZPRLXO3+TRLXO3(JROF,JLEV,JSTEP,IBL)*XPRLXG(JSTEP)
         ENDDO
         GFL(JROF,JLEV,YO3%MP,IBL)=(GFL(JROF,JLEV,YO3%MP,IBL) +&
           & ZALPHA*ZPRLXO3)*ZUSRLX
       ENDIF
     ENDDO
   ENDDO

   ! now surface fields (LNSP)

   DO JROF=1,ICEND

     ! Construct the mask
     ! longitudes
     IF (ALONRLX1*180._JPRB/RPI <=0.1 .AND.&
       & ALONRLX2*180._JPRB/RPI >= 359.9) THEN
       ZFLON1=1._JPRB
       ZFLON2=1._JPRB
     ELSE
       ZFLON1=1.0_JPRB/(1.0_JPRB+EXP(ZAX*(YDGSGEOM_NB%GELAM(IOFF+JROF-1)-ALONRLX1)))
       ZFLON2=1.0_JPRB-1.0_JPRB/(1.0_JPRB+EXP(ZAX*(YDGSGEOM_NB%GELAM(IOFF+JROF-1)-ALONRLX2)))
     ENDIF

     ! latitudes
     IF (ALATRLX1*180._JPRB/RPI >= 89.9) THEN
       ! here we go all the way to the north pole
       ZFLAT1=1._JPRB
     ELSE
       ZFLAT1=1.0_JPRB-1.0_JPRB/(1.0_JPRB+EXP(ZAY*(YDGSGEOM_NB%GELAT(IOFF+JROF-1)-ALATRLX1)))
     ENDIF
     IF (ALATRLX2*180._JPRB/RPI <= -89.9) THEN
       ! here we go all the way to the south pole
       ZFLAT2=1._JPRB
     ELSE
       ZFLAT2=1.0_JPRB/(1.0_JPRB+EXP(ZAY*(YDGSGEOM_NB%GELAT(IOFF+JROF-1)-ALATRLX2)))
     ENDIF
     ! contruct final mask
     IF(ALONRLX1 < ALONRLX2) THEN
       ZMASK=ZFLON1*ZFLON2*ZFLAT1*ZFLAT2
     ELSE
       ZMASK=(ZFLON1+ZFLON2)*ZFLAT1*ZFLAT2
     ENDIF

     ZALPHA=XRLXLP*ZMASK
     ZUSRLX=1.0_JPRB/(1.0_JPRB+ZALPHA)
     IF (LRLXLP) THEN
       GMVS(JROF,YT0%MSP,IBL)=(GMVS(JROF,YT0%MSP,IBL) +&
         & ZALPHA*ZGPRLXLP(JROF,1,IBL))*ZUSRLX
     ENDIF
   ENDDO
 ENDDO
ENDIF
WRITE(NULOUT,*) 'Relaxation finished'


!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RELAXGP',1,ZHOOK_HANDLE)
END SUBROUTINE RELAXGP
