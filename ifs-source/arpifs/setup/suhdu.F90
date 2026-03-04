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

SUBROUTINE SUHDU(YDGEOMETRY,YDML_GCONF,YDDYN,YDDYNA)

!**** *SUHDU*   - Initialize the tridiagonal horizontal diffusion operator
!                 (stretched version of ARPEGE).

!     Purpose.    Initialize the tridiagonal horizontal diffusion operator,
!     --------    (stretched version of ARPEGE).
!                 The order in RDHI is assumed to be the following one:
!                 - 2D model: VOR,DIV,equivalent height.
!                 - 3D hyd model: VOR,DIV,T,NH variables (SPD,SVD,NHX),
!                   spectral GFL,ln(prehyds). 
!                 It has to be consistent with the one of arrays ZSPX and ZSPY
!                 in SPCHOR/SPCHORAD.

!**   Interface.
!     ----------
!        *CALL* *SUHDU

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------
!        COMMON YOMDYN

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
!      K. YESSAD (APRIL 1994)

!     Modifications.
!     --------------
!      M.Hamrud    : 03-08-01 GFL introduction
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.YESSAD (Dec 2003): cleaning of horizontal diffusion.
!      F. Vana and K. Yessad: modify diffusion when LSLHD=T.
!      F. Vana     Feb 2005: Split of the global LSLHD key.
!      K.YESSAD (Dec 2005): Diffusion for NH variables.
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      F. Vana  13-Jan-2009 : Removed argument + case LSLHD_STATIC
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYNA  , ONLY : TDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : MYSETV

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT):: YDDYN
TYPE(TDYNA)    ,INTENT(INOUT):: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF

REAL(KIND=JPRB),ALLOCATABLE :: ZD(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZE(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZES(:,:)

INTEGER(KIND=JPIM) :: IL, ILEVG, IS0, ISE, IVTH, IVTHS, JSPGFL, JL, JN
INTEGER(KIND=JPIM) :: IM, JMLOC

REAL(KIND=JPRB) :: ZDT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "abor1.intfb.h"
#include "suhert.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHDU',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
  & YDMP=>YDGEOMETRY%YRMP, &
  & YDLAP=>YDGEOMETRY%YRLAP, YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDRIP=>YDML_GCONF%YRRIP,YGFL=>YDML_GCONF%YGFL, &
  & YDDIMF=>YDML_GCONF%YRDIMF)

ASSOCIATE(NUMSPFLDS=>YGFL%NUMSPFLDS, &
 & NSMAX=>YDDIM%NSMAX, NUMP=>YDDIM%NUMP, &
 & LSPT=>YDDIMF%LSPT, NFTHER=>YDDIMF%NFTHER, NS3D=>YDDIMF%NS3D, &
 & NFLEVL=>YDDIMV%NFLEVL, &
 & HDTIME_STRHD=>YDDYN%HDTIME_STRHD, LSTRHD=>YDDYN%LSTRHD, RDHI=>YDDYN%RDHI, &
 & RDHS=>YDDYN%RDHS, RDIDIV=>YDDYN%RDIDIV, RDIGFL=>YDDYN%RDIGFL, &
 & RDIPD=>YDDYN%RDIPD, RDISP=>YDDYN%RDISP, RDITG=>YDDYN%RDITG, &
 & RDIVD=>YDDYN%RDIVD, RDIVOR=>YDDYN%RDIVOR, RDSDIV=>YDDYN%RDSDIV, &
 & RDSVD=>YDDYN%RDSVD, RDSVOR=>YDDYN%RDSVOR, &
 & LSLHD_W=>YDDYNA%LSLHD_W,LSLHD_SVD=>YDDYNA%LSLHD_SVD,LNHX=>YDDYNA%LNHX,LSLHD_STATIC=>YDDYNA%LSLHD_STATIC, &
 & MYMS=>YDLAP%MYMS, NSE0L=>YDLAP%NSE0L, &
 & NPTRLL=>YDMP%NPTRLL, &
 & TDT=>YDRIP%TDT, &
 & GMR=>YDSPGEOM%GMR)
!      ----------------------------------------------------------------

!*       1.    INITIALIZE UNIFIED HORIZONTAL DIFFUSION OPERATOR.
!              -------------------------------------------------

ZDT=ABS(TDT)
IF(HDTIME_STRHD /= ZDT.AND.LSTRHD)THEN

  WRITE(NULOUT,'(A)')&
   & ' SUHDU: COMPUTES MAIN HOR DIFFUSION OPERATOR FOR STRETCHED ARPEGE.'  
  HDTIME_STRHD=ZDT

  !     1.1   ALL VARIABLES FOR "HDI".

  ALLOCATE(ZD(NFLEVL*NS3D+1,0:NSMAX))
  ALLOCATE(ZE(NFLEVL*NS3D+1,0:NSMAX))
  ALLOCATE(ZES(NFLEVL*NS3D+1,0:NSMAX))

  ! * Preliminary set ZD to 1., ZE,ZES to zero.

  ZD=1._JPRB
  ZE=0._JPRB
  ZES=0._JPRB

  ! * Definite computation of ZD,ZE,ZES, then call to SUHERT.

  WAVE_LOOP_HDI : DO JMLOC=1,NUMP
    IM=MYMS(JMLOC)
    IS0=NSE0L(JMLOC)

    !     1.1.1   VORTICITY, DIVERGENCE.

    DO JL=1,NFLEVL

      DO JN=IM,NSMAX
        ISE=IS0+JN+1-IM
        ZD(JL       ,JN)=1.0_JPRB+GMR(ISE,1)*RDIVOR(JL,JN)*ZDT
        ZD(JL+NFLEVL,JN)=1.0_JPRB+GMR(ISE,1)*RDIDIV(JL,JN)*ZDT
      ENDDO

      DO JN=IM,NSMAX-1
        ISE=IS0+JN+1-IM
        ZES(JL       ,JN)=GMR(ISE,2)*RDIVOR(JL,JN+1)*ZDT
        ZES(JL+NFLEVL,JN)=GMR(ISE,2)*RDIDIV(JL,JN+1)*ZDT
        ZE(JL       ,JN)=GMR(ISE,2)*RDIVOR(JL,JN)*ZDT
        ZE(JL+NFLEVL,JN)=GMR(ISE,2)*RDIDIV(JL,JN)*ZDT
      ENDDO

    ENDDO

    !     1.1.2   THERMODYNAMIC VARIABLES OTHER THAN GFL (3D).

    IVTH=0

    IF (LSPT) THEN

      ! * Temperature.

      IVTH=IVTH+1
      IL=NFLEVL*(2+IVTH-1)
      DO JL=1,NFLEVL
        ILEVG=NPTRLL(MYSETV)-1+JL
        DO JN=IM,NSMAX
          ISE=IS0+JN+1-IM
          ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDITG(ILEVG,JN)*ZDT
        ENDDO
        DO JN=IM,NSMAX-1
          ISE=IS0+JN+1-IM
          ZES(JL+IL,JN)=GMR(ISE,2)*RDITG(ILEVG,JN+1)*ZDT
          ZE(JL+IL,JN)=GMR(ISE,2)*RDITG(ILEVG,JN)*ZDT
        ENDDO
      ENDDO

    ENDIF

    IF (YDDYNA%LNHDYN) THEN

      ! * Pressure departure.

      IVTH=IVTH+1
      IL=NFLEVL*(2+IVTH-1)
      DO JL=1,NFLEVL
        DO JN=IM,NSMAX
          ISE=IS0+JN+1-IM
          ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDIPD(JL,JN)*ZDT
        ENDDO
        DO JN=IM,NSMAX-1
          ISE=IS0+JN+1-IM
          ZES(JL+IL,JN)=GMR(ISE,2)*RDIPD(JL,JN+1)*ZDT
          ZE(JL+IL,JN)=GMR(ISE,2)*RDIPD(JL,JN)*ZDT
        ENDDO
      ENDDO

    ENDIF

    IF (YDDYNA%LNHDYN) THEN

      ! * Vertical divergence variable.

      IVTH=IVTH+1
      IL=NFLEVL*(2+IVTH-1)
      DO JL=1,NFLEVL
        DO JN=IM,NSMAX
          ISE=IS0+JN+1-IM
          ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDIVD(JL,JN)*ZDT
        ENDDO
        DO JN=IM,NSMAX-1
          ISE=IS0+JN+1-IM
          ZES(JL+IL,JN)=GMR(ISE,2)*RDIVD(JL,JN+1)*ZDT
          ZE(JL+IL,JN)=GMR(ISE,2)*RDIVD(JL,JN)*ZDT
        ENDDO
      ENDDO

    ENDIF

    IF (LNHX) THEN

      ! * "NHX" variable.

      IVTH=IVTH+1
      IL=NFLEVL*(2+IVTH-1)
      DO JL=1,NFLEVL
        DO JN=IM,NSMAX
          ISE=IS0+JN+1-IM
          ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDIVD(JL,JN)*ZDT
        ENDDO
        DO JN=IM,NSMAX-1
          ISE=IS0+JN+1-IM
          ZES(JL+IL,JN)=GMR(ISE,2)*RDIVD(JL,JN+1)*ZDT
          ZE(JL+IL,JN)=GMR(ISE,2)*RDIVD(JL,JN)*ZDT
        ENDDO
      ENDDO

    ENDIF

    IF (IVTH /= NFTHER) THEN
      WRITE (NULOUT,*)
      WRITE (NULOUT,*)' SUHDU: ERROR NUMBER OF NON-GFL THERMODYNAMIC VARIABLES.'
      WRITE (NULOUT,*)
      CALL ABOR1('SUHDU: ABOR1 CALLED')
    ENDIF

    !     1.1.3   GFL VARIABLES (3D).

    ! * Diffusion only applies to spectral GFL.
    !   RDIGFL is assumed to have been filled by zeros in SUHDF
    !   for spectral GFL which are not diffused.

    DO JSPGFL=1,NUMSPFLDS
      IVTH=IVTH+1
      IL=NFLEVL*(2+IVTH-1)
      DO JL=1,NFLEVL
        DO JN=IM,NSMAX
          ISE=IS0+JN+1-IM
          ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDIGFL(JL,JN,JSPGFL)*ZDT
        ENDDO
        DO JN=IM,NSMAX-1
          ISE=IS0+JN+1-IM
          ZES(JL+IL,JN)=GMR(ISE,2)*RDIGFL(JL,JN+1,JSPGFL)*ZDT
          ZE(JL+IL,JN)=GMR(ISE,2)*RDIGFL(JL,JN,JSPGFL)*ZDT
        ENDDO
      ENDDO
    ENDDO

    !     1.1.4   SURFACE PRESSURE (3D) OR EQUIVALENT HEIGHT (2D).

    IL=NFLEVL*NS3D+1

    DO JN=IM,NSMAX
      ISE=IS0+JN+1-IM
      ZD(IL,JN)=1.0_JPRB+GMR(ISE,1)*RDISP(JN)*ZDT
    ENDDO

    DO JN=IM,NSMAX-1
      ISE=IS0+JN+1-IM
      ZES(IL,JN)=GMR(ISE,2)*RDISP(JN+1)*ZDT
      ZE(IL,JN)=GMR(ISE,2)*RDISP(JN)*ZDT
    ENDDO

    !     1.1.5   LU FACTORISATION OF PENTADIAGONAL MATRIX.

    IL=NFLEVL*NS3D+1

    CALL SUHERT(NSMAX+1-IM,IL,IL,ZD(1,IM),ZE(1,IM),ZES(1,IM),&
     & RDHI(1,IS0+1,1),RDHI(1,IS0+1,2),RDHI(1,IS0+1,3))  

  ENDDO WAVE_LOOP_HDI

  DEALLOCATE(ZD)
  DEALLOCATE(ZE)
  DEALLOCATE(ZES)

  !     1.2   VORTICITY, DIVERGENCE, VERTICAL DIVERGENCE FOR "HDS".

  IF ((.NOT.LSLHD_STATIC).AND.(LSLHD_W.OR.LSLHD_SVD)) THEN

    ! * First compute the number of prognostic variables on which the HDS
    !   diffusion is applied (currently IVTHS=(IL/NFLEVL) is between 0 and 4).
    IVTHS=0
    IF (LSLHD_W) IVTHS=IVTHS+2
    IF (YDDYNA%LNHDYN.AND.LSLHD_SVD) IVTHS=IVTHS+1
    IF (LNHX.AND.LSLHD_SVD) IVTHS=IVTHS+1
    IL=IVTHS*NFLEVL

    ! * Allocate ZD,ZE,ZES.
    ALLOCATE(ZD(IL,0:NSMAX))
    ALLOCATE(ZE(IL,0:NSMAX))
    ALLOCATE(ZES(IL,0:NSMAX))

    ! * Preliminary set ZD to 1., ZE,ZES to zero.
    ZD=1._JPRB
    ZE=0._JPRB
    ZES=0._JPRB

    ! * Definite computation of ZD,ZE,ZES, then call to SUHERT.

    WAVE_LOOP_HDS : DO JMLOC=1,NUMP
      IM=MYMS(JMLOC)
      IS0=NSE0L(JMLOC)
      IVTHS=0

      !     1.2.1   VORTICITY, DIVERGENCE.

      IF (LSLHD_W) THEN

        ! * Vorticity.

        IVTHS=IVTHS+1
        IL=NFLEVL*(IVTHS-1)
        DO JL=1,NFLEVL
          DO JN=IM,NSMAX
            ISE=IS0+JN+1-IM
            ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDSVOR(JL,JN)*ZDT
          ENDDO
          DO JN=IM,NSMAX-1
            ISE=IS0+JN+1-IM
            ZES(JL+IL,JN)=GMR(ISE,2)*RDSVOR(JL,JN+1)*ZDT
            ZE(JL+IL,JN)=GMR(ISE,2)*RDSVOR(JL,JN)*ZDT
          ENDDO
        ENDDO

        ! * Divergence.

        IVTHS=IVTHS+1
        IL=NFLEVL*(IVTHS-1)
        DO JL=1,NFLEVL
          DO JN=IM,NSMAX
            ISE=IS0+JN+1-IM
            ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDSDIV(JL,JN)*ZDT
          ENDDO
          DO JN=IM,NSMAX-1
            ISE=IS0+JN+1-IM
            ZES(JL+IL,JN)=GMR(ISE,2)*RDSDIV(JL,JN+1)*ZDT
            ZE(JL+IL,JN)=GMR(ISE,2)*RDSDIV(JL,JN)*ZDT
          ENDDO
        ENDDO

      ENDIF

      !     1.2.2   VERTICAL DIVERGENCE AND "NHX" FOR NH MODEL.

      IF (YDDYNA%LNHDYN.AND.LSLHD_SVD) THEN

        ! * Vertical divergence variable.

        IVTHS=IVTHS+1
        IL=NFLEVL*(IVTHS-1)
        DO JL=1,NFLEVL
          DO JN=IM,NSMAX
            ISE=IS0+JN+1-IM
            ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDSVD(JL,JN)*ZDT
          ENDDO
          DO JN=IM,NSMAX-1
            ISE=IS0+JN+1-IM
            ZES(JL+IL,JN)=GMR(ISE,2)*RDSVD(JL,JN+1)*ZDT
            ZE(JL+IL,JN)=GMR(ISE,2)*RDSVD(JL,JN)*ZDT
          ENDDO
        ENDDO

      ENDIF

      IF (LNHX.AND.LSLHD_SVD) THEN

        ! * "NHX" variable.

        IVTHS=IVTHS+1
        IL=NFLEVL*(IVTHS-1)
        DO JL=1,NFLEVL
          DO JN=IM,NSMAX
            ISE=IS0+JN+1-IM
            ZD(JL+IL,JN)=1.0_JPRB+GMR(ISE,1)*RDSVD(JL,JN)*ZDT
          ENDDO
          DO JN=IM,NSMAX-1
            ISE=IS0+JN+1-IM
            ZES(JL+IL,JN)=GMR(ISE,2)*RDSVD(JL,JN+1)*ZDT
            ZE(JL+IL,JN)=GMR(ISE,2)*RDSVD(JL,JN)*ZDT
          ENDDO
        ENDDO

      ENDIF

      !     1.2.3   LU FACTORISATION OF PENTADIAGONAL MATRIX.

      IL=IVTHS*NFLEVL

      CALL SUHERT(NSMAX+1-IM,IL,IL,ZD(1,IM),ZE(1,IM),ZES(1,IM),&
       & RDHS(1,IS0+1,1),RDHS(1,IS0+1,2),RDHS(1,IS0+1,3))  

    ENDDO WAVE_LOOP_HDS

    DEALLOCATE(ZD)
    DEALLOCATE(ZE)
    DEALLOCATE(ZES)

  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHDU',1,ZHOOK_HANDLE)
END SUBROUTINE SUHDU
