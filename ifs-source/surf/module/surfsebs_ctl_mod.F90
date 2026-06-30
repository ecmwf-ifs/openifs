MODULE SURFSEBS_CTL_MOD
CONTAINS
SUBROUTINE SURFSEBS_CTL(KIDIA,KFDIA,KLON,KTILES,KTVL,KTVH,&
 & PTMST,PSSKM1M,PTSKM1M,PQSKM1M,PDQSDT,PRHOCHU,PRHOCQU,&
 & PALPHAL,PALPHAS,PSSRFL,PFRTI,PTSRF,&
 & PSNS,PRSN,PHLICE, & 
 & PSLRFL,PTSKRAD,PEMIS,PASL,PBSL,PAQL,PBQL,&
 & PTHKICE,PSNTICE,&
 & YDCST,YDEXC,YDVEG,YDFLAKE,YDSOIL,&
 !out
 & PJS,PJQ,PSSK,PTSK,PSSH,PSLH,PSTR,PG0,&
 & PSL,PQL, &
 & LNEMOLIMTHK)

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RVTMP2
USE YOS_CST   , ONLY : TCST
USE YOS_EXC   , ONLY : TEXC
USE YOS_VEG   , ONLY : TVEG
USE YOS_FLAKE , ONLY : TFLAKE
USE YOS_SOIL  , ONLY : TSOIL
! (C) Copyright 2003- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFSEB computes surface energy balance and skin temperature 
!    for each tile. 

!  SURFSEB is called by VDFDIFH

!  METHOD:
!    A linear relation between lowest model level dry static 
!    energy and moisture and their fluxes is specified as input. 
!    The surface energy balance equation is used to eliminate 
!    the skin temperature as in the derivation of the 
!    Penmann-Monteith equation. 

!    The routine can also be used in stand alone simulations by
!    putting PASL and PAQL to zero and by specifying for PBSL and PBQL 
!    the forcing with dry static energy and specific humidity. 

!  AUTHOR:
!    A. Beljaars       ECMWF April 2003   

!  REVISION HISTORY:
!    J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!    E. Dutra/G.Balsamo   01-05-08 add lake tile
!    Linus Magnusson      10-09-28 Sea-ice
!    G. Balsamo/A. Beljaars 09-08-2013 snow scheme stability fix
!    I. Sandu              24-02-2014  Lambda skin values by vegetation type instead of tile
!    E. Dutra              10/10/2014  net longwave tiled 

!  INTERFACE: 

!    Integers (In):
!      KIDIA   :    Begin point in arrays
!      KFDIA   :    End point in arrays
!      KLON    :    Length of arrays
!      KTILES  :    Number of tiles
!      KTVL    :    Dominant low vegetation type 
!      KTVH    :    Dominant high vegetation type

!    Reals with tile index (In): 
!      PTMST    :    Time Step
!      PSSKM1M :    Dry static energy of skin at T-1           (J/kg)
!      PTSKM1M :    Skin temperature at T-1                    (K)
!      PQSKM1M :    Saturation specific humidity at PTSKM1M    (kg/kg)
!      PDQSDT  :    dqsat/dT at PTSKM1M                        (kg/kg K)
!      PRHOCHU :    Rho*Ch*|U|                                 (kg/m2s)
!      PRHOCQU :    Rho*Cq*|U|                                 (kg/m2s)
!      PALPHAL :    multiplier of ql in moisture flux eq.      (-)
!      PALPHAS :    multiplier of qs in moisture flux eq.      (-)
!      PSSRFL  :    Net short wave radiation at the surface    (W/m2)
!      PFRTI   :    Fraction of surface area of each tile      (-)
!      PTSRF   :    Surface temp. below skin (e.g. Soil or SST)(K) 
!      PSNS    :    Snow mass per unit area                    (kg/m2)
!      PSNS    :    Snow density                               (kg/m3)
!      PHLICE  :    Lake ice thickness                         (m) 

!    Reals independent of tiles (In):
!      PSLRFL  :    Net long wave radiation at the surface     (W/m2) 
!      PTSKRAD :    Mean skin temp. at radiation time level    (K)
!      PEMIS   :    Surface emissivity                         (-)
!      PASL    :    Asl in Sl=Asl*Js+Bsl                       (m2s/kg)
!      PBSL    :    Bsl in Sl=Asl*Js+Bsl                       (J/kg)
!      PAQL    :    Aql in Ql=Aql*Jq+Bql                       (m2s/kg)
!      PBQL    :    Bql in Ql=Aql*Jq+Bql                       (kg/kg)
!      PTHKICE :    Sea ice thickness                          (m)
!      PSNTICE :    Thickness of snow layer on sea ice         (m)

!    Reals with tile index (Out):
!      PJS     :    Flux of dry static energy                  (W/m2)
!      PQS     :    Moisture flux                              (kg/m2s)
!      PSSK    :    New dry static energy of skin              (J/kg)
!      PTSK    :    New skin temperature                       (K)
!      PSSH    :    Surface sensible heat flux                 (W/m2)
!      PSLH    :    Surface latent heat flux                   (W/m2)
!      PSTR    :    Surface net thermal radiation              (W/m2)
!      PG0     :    Surface ground heat flux (solar radiation  (W/m2)
!                   leakage is not included in this term)

!    Reals independent of tiles (Out):
!      PSL     :    New lowest model level dry static energy   (J/kg)
!      PQL     :    New lowest model level specific humidity   (kg/kg)

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation
!    This routine uses the method suggested by Polcher and Best
!    (the basic idea is to start with a linear relation between 
!     the lowest model level varibles and their fluxes, which is 
!     obtained after the downward elimination of the vertical 
!     diffusion tridiagonal matrix). 

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)  :: KLON
INTEGER(KIND=JPIM), INTENT(IN)  :: KTILES
INTEGER(KIND=JPIM), INTENT(IN)  :: KTVL(KLON) 
INTEGER(KIND=JPIM), INTENT(IN)  :: KTVH(KLON)
REAL(KIND=JPRB),    INTENT(IN)  :: PTMST

REAL(KIND=JPRB),    INTENT(IN)  :: PSSKM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PTSKM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PQSKM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PDQSDT(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PRHOCHU(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PRHOCQU(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PALPHAL(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PALPHAS(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PSSRFL(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PTSRF(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PSLRFL(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PTSKRAD(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PEMIS(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PASL(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PBSL(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PAQL(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PBQL(:) 
REAL(KIND=JPRB),    INTENT(IN)  :: PSNS(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PRSN(:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PHLICE(:)  
REAL(KIND=JPRB),    INTENT(IN)  :: PTHKICE(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PSNTICE(:)
TYPE(TCST),         INTENT(IN)  :: YDCST
TYPE(TEXC),         INTENT(IN)  :: YDEXC
TYPE(TVEG),         INTENT(IN)  :: YDVEG
TYPE(TFLAKE),       INTENT(IN)  :: YDFLAKE
TYPE(TSOIL),        INTENT(IN)  :: YDSOIL

REAL(KIND=JPRB),    INTENT(OUT) :: PJS(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PJQ(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSSK(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PTSK(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSSH(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSLH(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSTR(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PG0(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSL(:)
REAL(KIND=JPRB),    INTENT(OUT) :: PQL(:)
LOGICAL,            INTENT(IN)  :: LNEMOLIMTHK

!        Local variables

REAL(KIND=JPRB) ::    ZEJS1(KLON),ZEJS2(KLON),ZEJS4(KLON),&
 & ZEJQ1(KLON),ZEJQ2(KLON),ZEJQ4(KLON),&
 & ZDLWDT(KLON),ZSTABEXSN(KLON)  

REAL(KIND=JPRB) ::     ZCJS1(KLON,KTILES),ZCJS3(KLON,KTILES),&
 & ZCJQ2(KLON,KTILES),ZCJQ3(KLON,KTILES),&
 & ZCJQ4(KLON,KTILES),ZDSS1(KLON,KTILES),&
 & ZDSS2(KLON,KTILES),ZDSS4(KLON,KTILES),&
 & ZDJS1(KLON,KTILES),&
 & ZDJS2(KLON,KTILES),ZDJS4(KLON,KTILES),&
 & ZDJQ1(KLON,KTILES),ZDJQ2(KLON,KTILES),&
 & ZDJQ4(KLON,KTILES),ZCPTM1(KLON,KTILES),&
 & ZICPTM1(KLON,KTILES),ZLAMSK(KLON,KTILES)  

INTEGER(KIND=JPIM) :: JL,JT
REAL(KIND=JPRB) ::    ZDELTA,ZLAM,&
 & ZZ,ZZ1,ZZ2,ZZ3,ZY1,ZY2,ZY3,ZJS,ZJQ,&
 & ZCOEF1,ZLARGE,ZLARGESN,ZFRSR,ZRTTMEPS,ZIZZ,ZLARGEWT  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZLICE(KLON),ZLWAT(KLON),ZSNOW,ZSNOWHVEG,ZSNOW_GLACIER
REAL(KIND=JPRB) :: ZTSKMEAN,ZDLWDTSTAR

!      1. Initialize constants

IF (LHOOK) CALL DR_HOOK('SURFSEB_CTL_MOD:SURFSEB_CTL',0,ZHOOK_HANDLE)
ASSOCIATE(RCPD=>YDCST%RCPD, RLSTT=>YDCST%RLSTT, RLVTT=>YDCST%RLVTT, &
 & RSIGMA=>YDCST%RSIGMA, RTT=>YDCST%RTT, &
 & LELWDD=>YDEXC%LELWDD, LELWTL=>YDEXC%LELWTL, &
 & RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & RHOCI=>YDSOIL%RHOCI, RHOICE=>YDSOIL%RHOICE, RQSNCR=>YDSOIL%RQSNCR, &
 & RVLAMSK=>YDVEG%RVLAMSK, RVLAMSKS=>YDVEG%RVLAMSKS, RVTRSR=>YDVEG%RVTRSR)
ZDELTA=RVTMP2              ! moisture coeff. in cp  
ZLARGE=1.E10_JPRB          ! large number to impose Tsk=SST
ZLARGESN=50._JPRB          ! large number to constrain Tsk variations in case
                           !   of melting snow 
ZLARGEWT=20._JPRB          ! 1/lamdaSK(w)+1/lamdaSK(tvh)                           
ZRTTMEPS=RTT-0.2_JPRB      ! slightly below zero to start snow melt
ZSNOW=7._JPRB
ZSNOW_GLACIER=8._JPRB
ZSNOWHVEG=20._JPRB

!* FIND LAKE POINTS WITH ICE COVER 
DO JL=KIDIA,KFDIA
  IF(PHLICE(JL) > RH_ICE_MIN_FLK ) THEN
    ZLICE(JL)=1._JPRB
    ZLWAT(JL)=0._JPRB
  ELSE
    ZLICE(JL)=0._JPRB
    ZLWAT(JL)=1._JPRB
  ENDIF
ENDDO 

!* STABILITY FACTOR FOR EXPLICIT SNOW SCHEME AND FOREST-SNOW MIX
!      3. Compute a stability factor for the explicit snow scheme and
!         forest-snow mix to prevent numerical instability for "thin-rough" snow 
!         when running with long time-step (e.g. 1-hour)
 
DO JL=KIDIA,KFDIA
  ZSTABEXSN(JL)=PTMST/(MAX(RQSNCR,(PSNS(JL,1)/PRSN(JL,1)))*RHOCI*PRSN(JL,1)/RHOICE)
ENDDO 
!      2. Prepare tile independent arrays

IF (LELWDD) THEN
  DO JL=KIDIA,KFDIA
    ZDLWDT(JL)=-4._JPRB*PEMIS(JL)*RSIGMA*PTSKRAD(JL)**3
  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    ZDLWDT(JL)=4._JPRB*(PSLRFL(JL)/PTSKRAD(JL))
  ENDDO
ENDIF

!      3. Compute coefficients for dry static energy flux Js and 
!         moisture flux Jq, expressed in Sl,Ql and Ssk 

DO JT=1,KTILES
  DO JL=KIDIA,KFDIA
!or    ZCPTM1(JL,JT)=PSSKM1M(JL,JT)/PTSKM1M(JL,JT)
    ZICPTM1(JL,JT)=PTSKM1M(JL,JT)/PSSKM1M(JL,JT)
    ZCJS1(JL,JT)=PRHOCHU(JL,JT)
    ZCJS3(JL,JT)=-PRHOCHU(JL,JT)
    ZCJQ2(JL,JT)=PRHOCQU(JL,JT)*PALPHAL(JL,JT)
    ZCJQ3(JL,JT)=-PRHOCQU(JL,JT)*PALPHAS(JL,JT)*PDQSDT(JL,JT)*ZICPTM1(JL,JT)
    ZCJQ4(JL,JT)=-PRHOCQU(JL,JT)*PALPHAS(JL,JT)*&
     & (PQSKM1M(JL,JT)-PDQSDT(JL,JT)*PTSKM1M(JL,JT))  
  ENDDO
ENDDO

!      4. Compute coefficients for dry static energy flux Js and 
!         moisture flux Jq, expressed in Sl and Ql (Ssk has been 
!         eliminated using surface energy balance. 

DO JT=1,KTILES

  IF (JT == 2 .OR. JT == 5) THEN
    ZLAM=RLSTT
  ELSE
    ZLAM=RLVTT
  ENDIF
  ZFRSR=1.0_JPRB

  SELECT CASE(JT)
  CASE(1)
    DO JL=KIDIA,KFDIA
      ZLAMSK(JL,JT)=ZLARGE
    ENDDO
  CASE(2) ! Sea ice with possibly a snow layer on top of it
     IF (LNEMOLIMTHK) THEN
        DO JL=KIDIA,KFDIA
           IF (PSNTICE(JL) > 0.0_JPRB) THEN
              ! For now use same conductivity as snow on land
              ! Possible refinement: take fractional snow cover for thin layers of snow into account
              IF (PTSKM1M(JL,JT) >= PTSRF(JL,JT).AND.PTSKM1M(JL,JT) > ZRTTMEPS) THEN
                 ZLAMSK(JL,JT)=ZLARGESN
              ELSE
                 ZLAMSK(JL,JT)=ZSNOW ! Snow tile!!
              ENDIF
           ELSE
              IF (PTSKM1M(JL,JT) > PTSRF(JL,JT)) THEN
                 ZLAMSK(JL,JT)=RVLAMSKS(12)
              ELSE
                 ZLAMSK(JL,JT)=RVLAMSK(12)
              ENDIF
           ENDIF
        ENDDO
     ELSE
        DO JL=KIDIA,KFDIA
           IF (PTSKM1M(JL,JT) > PTSRF(JL,JT)) THEN
              ZLAMSK(JL,JT)=RVLAMSKS(12)
           ELSE
              ZLAMSK(JL,JT)=RVLAMSK(12)
           ENDIF
        ENDDO
     ENDIF
  CASE(3)
    DO JL=KIDIA,KFDIA
      !IF (PTSKM1M(JL,JT) > PTSRF(JL,JT)) THEN
      !  ZLAMSK(JL,JT)=RVLAMSKS(KTVL(JL))
      !ELSE
      !  ZLAMSK(JL,JT)=RVLAMSK(KTVL(JL))
      !ENDIF
      ZLAMSK(JL,JT)=ZLARGEWT      
    ENDDO
  ZFRSR=1.0_JPRB-RVTRSR(1)
  CASE(4)
    DO JL=KIDIA,KFDIA
      IF (PTSKM1M(JL,JT) > PTSRF(JL,JT)) THEN
        ZLAMSK(JL,JT)=RVLAMSKS(KTVL(JL))
      ELSE
        ZLAMSK(JL,JT)=RVLAMSK(KTVL(JL))
      ENDIF
    ENDDO
  ZFRSR=1.0_JPRB-RVTRSR(1)
  CASE(5)
    DO JL=KIDIA,KFDIA
      IF (PTSKM1M(JL,JT) >= PTSRF(JL,JT).AND.PTSKM1M(JL,JT) > ZRTTMEPS) THEN
        ZLAMSK(JL,JT)=ZLARGESN
      ELSEIF (SUM(PSNS(JL,:),DIM=1)>9000.0_JPRB) THEN
        ZLAMSK(JL,JT)=ZSNOW_GLACIER
      ELSE
        ZLAMSK(JL,JT)=ZSNOW
      ENDIF
    ENDDO
  CASE(6)
     DO JL=KIDIA,KFDIA
      IF (PTSKM1M(JL,JT) > PTSRF(JL,JT)) THEN
        ZLAMSK(JL,JT)=RVLAMSKS(KTVH(JL))
      ELSE
        ZLAMSK(JL,JT)=RVLAMSK(KTVH(JL))
      ENDIF
    ENDDO
    ZFRSR=1.0_JPRB-RVTRSR(3)
  CASE(7)
     DO JL=KIDIA,KFDIA
      IF (PTSKM1M(JL,JT) > PTSRF(JL,JT)) THEN
        ZLAMSK(JL,JT)=RVLAMSKS(KTVH(JL))/(1+ZSTABEXSN(JL)*RVLAMSKS(KTVH(JL)))
      ELSE
        ZLAMSK(JL,JT)=ZSNOWHVEG/(1+ZSTABEXSN(JL)*ZSNOWHVEG)
      ENDIF
    ENDDO
    ZFRSR=1.0_JPRB-RVTRSR(3)
  CASE(8)
    DO JL=KIDIA,KFDIA
      IF (PTSKM1M(JL,JT) > PTSRF(JL,JT)) THEN
        ZLAMSK(JL,JT)=RVLAMSKS(8)
      ELSE
        ZLAMSK(JL,JT)=RVLAMSK(8)
      ENDIF
    ENDDO
   CASE(9)
    DO JL=KIDIA,KFDIA
      ZLAMSK(JL,JT)=ZLARGE
    ENDDO
  END SELECT

  DO JL=KIDIA,KFDIA
    IF (JT == 9 )  THEN                     
      ZLAM=RLVTT*ZLWAT(JL)+RLSTT*ZLICE(JL)  
    ENDIF                                   
    ZCOEF1=ZLAM-RCPD*PTSKM1M(JL,JT)*ZDELTA
    ZZ=(ZDLWDT(JL)-ZLAMSK(JL,JT))*ZICPTM1(JL,JT)+ZCJS3(JL,JT)&
     & +ZCOEF1*ZCJQ3(JL,JT)  
    ZIZZ = 1.0_JPRB/ZZ

    ZDSS1(JL,JT)=-ZCJS1(JL,JT)*ZIZZ
    ZDSS2(JL,JT)=-ZCOEF1*ZCJQ2(JL,JT)*ZIZZ
    ZDSS4(JL,JT)=(-PSSRFL(JL,JT)*ZFRSR-PSLRFL(JL)+ZDLWDT(JL)*PTSKRAD(JL)&
           &-ZCOEF1*ZCJQ4(JL,JT)-ZLAMSK(JL,JT)*PTSRF(JL,JT))*ZIZZ
    ZDJS1(JL,JT)=ZCJS1(JL,JT)+ZCJS3(JL,JT)*ZDSS1(JL,JT)
    ZDJS2(JL,JT)=             ZCJS3(JL,JT)*ZDSS2(JL,JT)
    ZDJS4(JL,JT)=ZCJS3(JL,JT)*ZDSS4(JL,JT)
    ZDJQ1(JL,JT)=             ZCJQ3(JL,JT)*ZDSS1(JL,JT)
    ZDJQ2(JL,JT)=ZCJQ2(JL,JT)+ZCJQ3(JL,JT)*ZDSS2(JL,JT)
    ZDJQ4(JL,JT)=ZCJQ3(JL,JT)*ZDSS4(JL,JT)+ZCJQ4(JL,JT)
  ENDDO
ENDDO

!      5.  Average coefficients over tiles

DO JL=KIDIA,KFDIA
  ZEJS1(JL)=0.0_JPRB
  ZEJS2(JL)=0.0_JPRB
  ZEJS4(JL)=0.0_JPRB
  ZEJQ1(JL)=0.0_JPRB
  ZEJQ2(JL)=0.0_JPRB
  ZEJQ4(JL)=0.0_JPRB
ENDDO

DO JT=1,KTILES
  DO JL=KIDIA,KFDIA
    ZEJS1(JL)=ZEJS1(JL)+PFRTI(JL,JT)*ZDJS1(JL,JT)
    ZEJS2(JL)=ZEJS2(JL)+PFRTI(JL,JT)*ZDJS2(JL,JT)
    ZEJS4(JL)=ZEJS4(JL)+PFRTI(JL,JT)*ZDJS4(JL,JT)
    ZEJQ1(JL)=ZEJQ1(JL)+PFRTI(JL,JT)*ZDJQ1(JL,JT)
    ZEJQ2(JL)=ZEJQ2(JL)+PFRTI(JL,JT)*ZDJQ2(JL,JT)
    ZEJQ4(JL)=ZEJQ4(JL)+PFRTI(JL,JT)*ZDJQ4(JL,JT)
  ENDDO
ENDDO

!      6.  Eliminate mean fluxes to find Sl and Ql

DO JL=KIDIA,KFDIA
  ZZ1=1.-PASL(JL)*ZEJS1(JL)
  ZZ2=PAQL(JL)*ZEJS2(JL)
  ZZ3=PBSL(JL)*ZEJS1(JL)+PBQL(JL)*ZEJS2(JL)+ZEJS4(JL)
  ZY1=1.0_JPRB-PAQL(JL)*ZEJQ2(JL)
  ZY2=PASL(JL)*ZEJQ1(JL)
  ZY3=PBSL(JL)*ZEJQ1(JL)+PBQL(JL)*ZEJQ2(JL)+ZEJQ4(JL)

  ZZ=ZZ1*ZY1-ZZ2*ZY2
  ZIZZ=1.0_JPRB/ZZ
  ZJQ=(ZZ3*ZY2+ZZ1*ZY3)*ZIZZ
  ZJS=(ZZ2*ZY3+ZZ3*ZY1)*ZIZZ

  PSL(JL)=PASL(JL)*ZJS+PBSL(JL)
  PQL(JL)=PAQL(JL)*ZJQ+PBQL(JL)

ENDDO

!      7.  Compute tile dependent fluxes and skin values

DO JT=1,KTILES
  IF (JT == 2 .OR. JT == 5) THEN
    ZLAM=RLSTT
  ELSE
    ZLAM=RLVTT
  ENDIF
  DO JL=KIDIA,KFDIA
    IF (JT == 9 ) THEN                      
      ZLAM=RLVTT*ZLWAT(JL)+RLSTT*ZLICE(JL)  
    ENDIF                                   
    PSSK(JL,JT)=ZDSS1(JL,JT)*PSL(JL)+ZDSS2(JL,JT)*PQL(JL)+ZDSS4(JL,JT)
    PJS(JL,JT) =ZDJS1(JL,JT)*PSL(JL)+ZDJS2(JL,JT)*PQL(JL)+ZDJS4(JL,JT)
    PJQ(JL,JT) =ZDJQ1(JL,JT)*PSL(JL)+ZDJQ2(JL,JT)*PQL(JL)+ZDJQ4(JL,JT)
    PTSK(JL,JT)=PSSK(JL,JT)*ZICPTM1(JL,JT)

!         Surface heat fluxes

    PSSH(JL,JT)=PJS(JL,JT)-RCPD*PTSKM1M(JL,JT)*ZDELTA*PJQ(JL,JT)
    PSLH(JL,JT)=ZLAM*PJQ(JL,JT)
    PSTR(JL,JT)=PSLRFL(JL)+ZDLWDT(JL)*(PTSK(JL,JT)-PTSKRAD(JL))
    PG0 (JL,JT)=ZLAMSK(JL,JT)*(PTSK(JL,JT)-PTSRF(JL,JT))
  ENDDO
ENDDO


!            8. Re-compute tiled net longwave centered in the Mean grib-box new skin temperature 
!               This allow energy conservation

IF ( LELWTL ) THEN
  DO JL=KIDIA,KFDIA

    ZTSKMEAN=0._JPRB
    ! calculate mean grid-box skin temperature 
    DO JT=1,KTILES
      ZTSKMEAN=ZTSKMEAN+PFRTI(JL,JT)*PTSK(JL,JT)
    ENDDO

    ! Linearized lonwave emission - centered in ZTSKMEAN
    ZDLWDTSTAR=-4._JPRB*PEMIS(JL)*RSIGMA*ZTSKMEAN**3
    ! re-compute PSTR (net longwave) 
    DO JT=1,KTILES
      PSTR(JL,JT)=PSLRFL(JL)+ZDLWDTSTAR*(PTSK(JL,JT)-ZTSKMEAN)
    ENDDO

  ENDDO
ENDIF 

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFSEB_CTL_MOD:SURFSEB_CTL',1,ZHOOK_HANDLE)

!      7.  Wrap up

END SUBROUTINE SURFSEBS_CTL
END MODULE SURFSEBS_CTL_MOD
