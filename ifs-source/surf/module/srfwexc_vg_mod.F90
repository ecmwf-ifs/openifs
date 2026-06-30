MODULE SRFWEXC_VG_MOD
CONTAINS
SUBROUTINE SRFWEXC_VG(KIDIA,KFDIA,KLON,KLEVS,KTILES,&
 & PTMST,KTVL,KTVH,KSOTY,PSDOR,PFRTI,PEVAPTI,&
 & PWSAM1M,PTSAM1M,PCUR,&
 & PTSFC,PTSFL,PMSN,PEMSSN,PEINTTI,PEVAPSNW,&
 & YDSOIL,YDVEG,YDURB,&
 & PROS,PCFW,PRHSW,&
 & PSAWGFL,PFWEL1,PFWE234,&
 & LDLAND,PDHWLS)  

USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RHOH2O
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_VEG   , ONLY : TVEG
USE YOS_URB   , ONLY : TURB

! (C) Copyright 1993- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.


!**** *SRFWEXC_VG* -  COMPUTES THE FLUXES BETWEEN THE SOIL LAYERS AND
!                  THE RIGHT-HAND SIDE OF THE SOIL WATER EQUATIONS.
!
!     PURPOSE.
!     --------
!          THIS ROUTINE COMPUTES THE DIFFERENT COEFFICIENTS IN THE
!     SOIL MOISTURE EQUATIONS (BEFORE SNOW MELTS). THE AIM IS TO COMPUTE
!     THE MODIFIED DIFFUSIVITIES AND THE RIGHT-HAND SIDE OF THE EQUATIONS.
!     IT SHOULD BE FOLLOWED BY A CALL TO *SRFWDIF* AND *SRFWINC*.

!**   INTERFACE.
!     ----------
!          *SRFWEXC_VG* IS CALLED FROM *SURFTSTP*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*      NUMBER OF SURFACE LAYERS
!    *KTILES*     NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                 OF SURFACE BOUNDARY CONDITION)
!    *KDHVWLS*    Number of variables for soil water budget
!    *KDHFWLS*    Number of fluxes for soil water budget

!    *KTVL*       VEGETATION TYPE FOR LOW VEGETATION FRACTION
!    *KTVH*       VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!    *KSOTY*      SOIL TYPE                                   (1-7)

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S
!    *PSDOR*      OROGRAPHIC PARAMETER                           m
!    *PFRTI*      TILE FRACTIONS                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10: URBAN
!    *PEVAPTI*      SURFACE MOISTURE FLUX                      KG/M**2/S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PWSAM1M*    MULTI-LAYER SOIL MOISTURE                    M**3/M**3
!    *PTSAM1M*    SOIL TEMPERATURE ALL LAYERS                    K
!    *PCUR*       URBAN COVER                                   (0-1)
!    *PTSFC*      CONVECTIVE THROUGHFALL                       KG/M**2/S
!    *PTSFL*      LARGE SCALE THROUGHFALL                      KG/M**2/S
!    *PMSN*       SNOW MELTING                                 KG/M**2/S
!    *PEMSSN*     EVAPORATIVE MISMATCH RESULTING FROM
!                  CLIPPING THE SNOW TO ZERO (AFTER P-E)       KG/M**2/S
!    *PEINTTI*    TILE EVAPORATION SEEN BY THE INTERCEPTION
!                 LAYER (INCLUDES NUMERICAL EVAPORATION
!                 MISMATCHES, FOR TILE 3, AND DEW DEPOSITION
!                 FOR TILES 3,4,6,7,8)                         KG/M**2/S
!    *PEVAPSNW*   EVAPORATION FROM SNOW UNDER FOREST           KG/M**2/S

!     OUTPUT PARAMETERS (REAL):
!    *PCFW*       MODIFIED DIFFUSIVITIES                         M
!    *PRHSW*      RIGHT-HAND SIDE OF SOIL MOISTURE EQUATIONS   m**3/m**3
!    *PROS*       RUN-OFF FOR THE SURFACE LAYER                kg/m**2
!    *PFWEL1*     BARE GROUND AND TOP LAYER EXTRACTION
!                 CONTRIBUTION TO EVAPORATION FROM
!                 THE SKIN AND TOP LAYER                       KG/M**2/S
!    *PFWE234*    ROOT EXTRACTION FROM LAYERS 2+3+4            KG/M**2/S

!     INSTANTANEOUS DIAGNOSTIC OUTPUT PARAMETERS (REAL):
!    *PSAWGFL*    GRAVITY PART OF WATER FLUX                   KG/M**2/S
!                   (positive downwards, at layer bottom)

!     OUTPUT PARAMETERS (DIAGNOSTIC):
!    *PDHWLS*     Diagnostic array for soil water (see module yomcdh)

!     METHOD.
!     -------
!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION.

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     ORIGINAL :
!     P.VITERBO      E.C.M.W.F.      9/02/93
!     P.VITERBO      E.C.M.W.F.      26-3-99
!        (Interface to tiling)
!     D.SALMOND      E.C.M.W.F.      000515
!     P.VITERBO      E.C.M.W.F.      17-05-2000
!        (Surface DDH for TILES)
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!     P. Viterbo    ECMWF    24-05-2004      Change surface units
!     G. Balsamo    ECMWF    08-01-2006      Include Van Genuchten Hydro.
!     G. Balsamo    ECMWF    11-01-2006      Include sub-grid surface runoff
!     G. Balsamo    ECMWF    03-07-2006      Add soil type
!     E. Dutra               07-07-2008      clean number of tiles dependence 
!     G. Balsamo    ECMWF    15-09-2009      protect surf-runoff for occasional overshooting
!     G. Balsamo    ECMWF    12-05-2010      cleaning and unit fix in VIC
!     F. Vana                05-Mar-2015     Support for single precision
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS
INTEGER(KIND=JPIM), INTENT(IN)   :: KTILES
INTEGER(KIND=JPIM), INTENT(IN)   :: KTVL(:)
INTEGER(KIND=JPIM), INTENT(IN)   :: KTVH(:)
INTEGER(KIND=JPIM), INTENT(IN)   :: KSOTY(:)

REAL(KIND=JPRB),    INTENT(IN)   :: PTMST
REAL(KIND=JPRB),    INTENT(IN)   :: PSDOR(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PTSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCUR(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PTSFC(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PTSFL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PMSN(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEMSSN(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEINTTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPSNW(:)
LOGICAL,            INTENT(IN)   :: LDLAND(:)
TYPE(TSOIL),        INTENT(IN)   :: YDSOIL
TYPE(TVEG),         INTENT(IN)   :: YDVEG
TYPE(TURB),         INTENT(IN)   :: YDURB

REAL(KIND=JPRB),    INTENT(INOUT):: PDHWLS(:,:,:)
REAL(KIND=JPRB),    INTENT(INOUT):: PFWEL1(:)
REAL(KIND=JPRB),    INTENT(INOUT):: PFWE234(:)

REAL(KIND=JPRB),    INTENT(OUT)  :: PROS(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PCFW(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PRHSW(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PSAWGFL(:,:)

!*         0.2    DECLARATION OF LOCAL VARIABLES.
!                 ----------- -- ----- ----------

REAL(KIND=JPRB) :: ZSAWEXT(KLON,KLEVS)
REAL(KIND=JPRB) :: ZROOTW(KLON,KLEVS,KTILES),ZEXT(KLON,KTILES)
REAL(KIND=JPRB) :: ZLIQ(KLON,KLEVS),ZF(KLON,KLEVS),ZCONDS(KLON)
LOGICAL :: LLFREEZ

INTEGER(KIND=JPIM) :: ITYP(KLON), JK, JL, JTILE, JS

REAL(KIND=JPRB) :: ZEPSILON,ZEPSILONADD
REAL(KIND=JPRB) :: Z_RHOH20, ZD, &
 & ZEXTK, ZINFMAX, ZK, ZKM, ZPSFR, ZROC, ZROL, ZDSURF, ZKSURF, &
 & ZROT, ZINVTMST, ZW, ZWM, ZWSFL,ZEVAP, ZFF, ZFFM,ZHOH2O, ZWS, &
 & ZDMAX, ZDMIN, ZALPHA, ZWFAC, ZLAM, ZMFAC, ZRMFAC, &
 & ZWCONS, ZKMD, &
 & ZRSFL, ZROEFF, ZSIGOR, ZBWS, ZB1, ZBM, ZWMAX, &
 & ZCONW1, ZLYEPS, ZLYSIC, ZVOL, ZROS, ZSUM, ZLIMRS
REAL(KIND=JPRD) :: ZDD, ZKD, ZSE

REAL(KIND=JPRB) :: ZFRK(KLEVS),ZWK(KLEVS),ZWMK(KLEVS)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
!*         1.    SET UP SOME CONSTANTS.
!                --- -- ---- ----------
!    SECURITY PARAMETERS

IF (LHOOK) CALL DR_HOOK('SRFWEXC_VG_MOD:SRFWEXC_VG',0,ZHOOK_HANDLE)
ASSOCIATE(LESSRO=>YDSOIL%LESSRO, RDAW=>YDSOIL%RDAW, RDMAXM=>YDSOIL%RDMAXM, &
 & RDMINM=>YDSOIL%RDMINM, RLAMBDAM=>YDSOIL%RLAMBDAM, RMFACM=>YDSOIL%RMFACM, &
 & RMVGALPHA=>YDSOIL%RMVGALPHA, RPSFR=>YDSOIL%RPSFR, &
 & RSIGORMAX=>YDSOIL%RSIGORMAX, RSIGORMIN=>YDSOIL%RSIGORMIN, &
 & RSIMP=>YDSOIL%RSIMP, RSRDEP=>YDSOIL%RSRDEP, RTF1=>YDSOIL%RTF1, &
 & RTF2=>YDSOIL%RTF2, RTF3=>YDSOIL%RTF3, RTF4=>YDSOIL%RTF4, &
 & RWCAPM=>YDSOIL%RWCAPM, RWCONSM=>YDSOIL%RWCONSM, RWPWPM=>YDSOIL%RWPWPM, &
 & RWRESTM=>YDSOIL%RWRESTM, RWSATM=>YDSOIL%RWSATM, &
 & RVROOTSA=>YDVEG%RVROOTSA, RURBALP=>YDURB%RURBALP, RURBCON=>YDURB%RURBCON,&
 & RURBLAM=>YDURB%RURBLAM, RURBSAT=>YDURB%RURBSAT, RURBSRES=>YDURB%RURBSRES)

ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

!    COMPUTATIONAL CONSTANTS.

ZHOH2O=1.0_JPRB/RHOH2O
ZPSFR=1.0_JPRB/RPSFR
ZINVTMST=1.0_JPRB/PTMST

LLFREEZ=.TRUE.

!  Water balance closure fix 
ZEPSILONADD=0._JPRB
IF ( YDSOIL%LEWBSOILFIX ) THEN
  ZEPSILONADD=ZEPSILON
ENDIF 
!     ------------------------------------------------------------------
!*          2. COMPUTATION OF THE MODIFIED DIFFUSIVITY COEFFICIENTS
!              ----------- -- --- -------- ----------- ------------
!*          2.1 Preliminary quantities related to root extraction
!               Compute first liquid fraction of soil water to 
!               be used later in stress functions.
DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    IF(PTSAM1M(JL,JK) < RTF1.AND.PTSAM1M(JL,JK) > RTF2) THEN
      ZF(JL,JK)=0.5_JPRB*(1.0_JPRB-SIN(RTF4*(PTSAM1M(JL,JK)-RTF3)))
    ELSEIF (PTSAM1M(JL,JK) <= RTF2) THEN
      ZF(JL,JK)=1.0_JPRB
    ELSE
      ZF(JL,JK)=0.0_JPRB
    ENDIF
    JS=KSOTY(JL)
    ZLIQ(JL,JK)=MAX(RWPWPM(JS),MIN(RWCAPM(JS),PWSAM1M(JL,JK)*(1._JPRB-ZF(JL,JK))))
  ENDDO
ENDDO

!*          2.2 Preliminary quantities related to root extraction
DO JL=KIDIA,KFDIA
  ZCONDS(JL)=0.0_JPRB
  DO JTILE=1,KTILES
    ZEXT(JL,JTILE)=0.0_JPRB
  ENDDO
ENDDO
DO JTILE=1,KTILES
  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      ZROOTW(JL,JK,JTILE)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

TILES: DO JTILE=1,KTILES
! no liquid evaporation contribution to soil from tile 1, 2, 5
! no root extraction from tile 3 and 8
    IF ( JTILE /= 4 .AND. JTILE /= 6 .AND. JTILE /= 7 ) CYCLE TILES 
    IF (JTILE == 4) THEN
      DO JL=KIDIA,KFDIA
        ITYP(JL)=KTVL(JL)
      ENDDO
    ELSEIF (JTILE == 6.OR. JTILE == 7) THEN
      DO JL=KIDIA,KFDIA
        ITYP(JL)=KTVH(JL)
      ENDDO
    ENDIF
! Layers 1-klevs
    LAYERS: DO JK=1,KLEVS
      DO JL=KIDIA,KFDIA
        JS=KSOTY(JL)
        IF (JTILE == 7) THEN
          ZEVAP=PFRTI(JL,JTILE)*(PEVAPTI(JL,JTILE)-PEVAPSNW(JL))
        ELSE
          ZEVAP=PFRTI(JL,JTILE)*PEVAPTI(JL,JTILE)
        ENDIF
        IF (ZEVAP > 0.0_JPRB) THEN
          IF (JK == 1) THEN
            ZROOTW(JL,JK,JTILE)=1.0_JPRB
          ENDIF
        ELSE
          ZROOTW(JL,JK,JTILE)=RVROOTSA(JK,ITYP(JL))*(ZLIQ(JL,JK)-RWPWPM(JS)+ZEPSILONADD)
        ENDIF
        ZEXT(JL,JTILE)=ZEXT(JL,JTILE)+ZROOTW(JL,JK,JTILE)
      ENDDO
    ENDDO LAYERS
    DO JL=KIDIA,KFDIA
      ZEXT(JL,JTILE)=MAX(ZEPSILON,ZEXT(JL,JTILE))
      !ZEXT(JL,JTILE)= SIGN(MAX(ABS(ZEXT(JL,JTILE)),ZEPSILON),ZEXT(JL,JTILE))
    ENDDO
ENDDO TILES

!*          2.3 COEFFICIENTS FOR TOP LAYER

DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN

!           HYDRAULIC PROPERTIES.

!           VAN GENUCHTEN

    JS=KSOTY(JL)
    ZDMAX=RDMAXM(JS)
    ZDMIN=RDMINM(JS)
    ZW=MAX(MAX(PWSAM1M(JL,1),PWSAM1M(JL,2)),RWRESTM(JS))
    ZSE=(ZW-RWRESTM(JS))/(RWSATM(JS)-RWRESTM(JS))
    ZWFAC=RWSATM(JS)-RWRESTM(JS)
    ZLAM=RLAMBDAM(JS)
    ZMFAC=RMFACM(JS)
    ZRMFAC=1./ZMFAC
    ZWCONS=RWCONSM(JS)
    ZALPHA=RMVGALPHA(JS)
    IF (ZW.LE.(1.001_JPRB*RWRESTM(JS))) THEN
      ZK=0.0_JPRB
      ZD=ZDMIN
    ELSEIF ((ZW.LE.(0.999_JPRB*RWSATM(JS))).AND.(ZW.GT.(1.001_JPRB*RWRESTM(JS)))) THEN
      ZKD=ZWCONS*ZSE**ZLAM* &
 &    (1.0_JPRB-((1.0_JPRB-(ZSE**ZRMFAC))**ZMFAC))**2.0_JPRB
      ! perhaps better to be rewritten
      ZDD=((1.0_JPRB-ZMFAC)*ZWCONS)/ &
 &    SIGN(MAX(ABS(ZALPHA*ZMFAC*ZWFAC),ZEPSILON),ZALPHA*ZMFAC*ZWFAC)
      ZDD= ZDD  *(ZSE**(ZLAM-ZRMFAC)) &
 &    *(((1.0_JPRB-(ZSE**ZRMFAC))**ZMFAC) &
 &    +((1.0_JPRB-(ZSE**ZRMFAC))**(-ZMFAC))-2.0_JPRB)
      ZK=ZKD
      ZD=ZDD
    ELSE
      ZK=RWCONSM(JS)
      ZD=ZDMAX
    ENDIF
    IF (LLFREEZ) THEN
      ZFF=MIN(ZF(JL,1),ZF(JL,2))
      ZD=ZFF*ZDMIN+(1.0_JPRB-ZFF)*ZD
! NOTE ZK = 0 for frozen soil
      ZK=ZFF*0.0_JPRB+(1.0_JPRB-ZFF)*ZK
    ENDIF
!
    IF (LESSRO) THEN

!          SURFACE RUNOFF DUE TO VARIABLE INFILTRATION CAPACITY (VIC)
!          ----------------------------------------------------------

!              Relative soil saturation is defined using 
!              liquid water in upper RSRDEP layer
      DO JK=1,KLEVS
        IF (JK > 1) THEN
          ZSUM=SUM(RDAW(1:JK-1))
        ELSE
          ZSUM=0.
        ENDIF
        ZFRK(JK)=MAX(0.0_JPRB,(MIN(RSRDEP,SUM(RDAW(1:JK)))-ZSUM)/RDAW(JK))
      ENDDO
      ZRSFL=PTSFL(JL)+PMSN(JL)+PTSFC(JL)                            !Units kg/m2/s
      IF (ZRSFL.GT.0.0_JPRB)THEN
!          SUBGRID SATURATION COEFFICIENTS
        ZROEFF=MAX(0.0_JPRB,(PSDOR(JL)-RSIGORMIN))/(PSDOR(JL)+RSIGORMAX)
        ZBWS=MAX(MIN(ZROEFF,0.5_JPRB),0.01_JPRB)
        ZB1=1.0_JPRB+ZBWS
        ZBM=1.0_JPRB/ZB1
        ZWK(:)=MAX(0.0_JPRB,(PWSAM1M(JL,:)-RWPWPM(JS)))*(1.0_JPRB-ZF(JL,:))+RWPWPM(JS)
        ZW=SUM(ZFRK(:)*ZWK(:)*RDAW(:))                              !Units m

        ZWMK(:)=((RWSATM(JS)-RWPWPM(JS))*(1.-ZF(JL,:))+RWPWPM(JS))
        IF ( KTILES .GT. 9 ) THEN
         ZWMK(:)=(((1.0_JPRB-PCUR(JL))*RWSATM(JS) + PCUR(JL)*RURBSAT &
           & -RWPWPM(JS))*(1.-ZF(JL,:))+RWPWPM(JS))
        ENDIF
        ZWMAX=SUM(ZFRK(:)*ZWMK(:)*RDAW(:))                          !Units m

        ZCONW1=ZWMAX*ZB1
        ZLYEPS=MAX(0.0_JPRB,ZW-ZWMAX)                               !Units m
        ZLIMRS=-1.0_JPRB*ZLYEPS
        IF (ZLYEPS.GT.0.1_JPRB*ZWMAX) THEN
           ZLIMRS=0.0_JPRB
        ENDIF
!       VIC Saturated area calculation
        ZLYSIC=MIN(1.0_JPRB,MAX(0.0_JPRB,(ZW-ZLYEPS)/ZWMAX))                     !Units -
        ZVOL=MAX(0.0_JPRB,(1.0_JPRB-ZLYSIC)**ZBM-ZRSFL/(ZINVTMST*RHOH2O*ZCONW1)) !Units -
!       Saturation excess runoff
        ZROS=ZRSFL/(ZINVTMST*RHOH2O)-MAX(ZWMAX-ZW,ZLIMRS)           !Units m
        IF (ZVOL.GT.0.0_JPRB) THEN
!       VIC runoff
          ZROS=ZROS+ZWMAX*ZVOL**ZB1                                 !Units m
        END IF
        ZROS=MAX(ZROS,0.0_JPRB)                                     !Units m
        ZROT=ZROS*RHOH2O*ZINVTMST                                   !Units kg/m2/s
      ELSE
        ZROT=0.
      ENDIF

    ELSE

!          SURFACE RUNOFF DUE TO INFILTRATION RATE LIMIT.
!          ----------------------------------------------

      ZDSURF=ZDMAX
      ZKSURF=RWCONSM(JS)
      IF (LLFREEZ) THEN
        ZFF=ZF(JL,1)
        ZDSURF=ZFF*ZDMIN+(1.-ZFF)*ZDMAX
        ZKSURF=ZFF*0.+(1.-ZFF)*RWCONSM(JS)
      ENDIF
      ZINFMAX=(ZDSURF*(RWSATM(JS)-PWSAM1M(JL,1))/(0.5_JPRB*RDAW(1))+ZKSURF)*RHOH2O
      IF ( KTILES .GT. 9 ) THEN
         ZINFMAX=(ZDSURF*(((1.0_JPRB-PCUR(JL))*RWSATM(JS) + PCUR(JL)*RURBSAT)-PWSAM1M(JL,1))&
          & /(0.5_JPRB*RDAW(1))+ZKSURF)*RHOH2O
      ENDIF
!          LARGE SCALE PRECIPITATION
      ZROL=MAX(0.0_JPRB,PTSFL(JL)+PMSN(JL)-ZINFMAX)

!          CONVECTIVE PRECIPITATION
      ZROC=MAX(0.0_JPRB,RPSFR*PTSFC(JL)-ZINFMAX)*ZPSFR
      ZROT=ZROL+ZROC
    ENDIF

!  Contribution of throughfall, melting, and runoff to the r.h.s.
    ZWSFL=PTSFL(JL)+PMSN(JL)+PTSFC(JL)-ZROT

!*             TILE CONTRIBUTIONS
!              ------------------
!  Tile by tile contribution of evaporation to the r.h.s.
!  Tile 3
    ZWSFL=ZWSFL+PFRTI(JL,3)*PEVAPTI(JL,3)-PEINTTI(JL,3)
!  Tile 4
    ZEXTK=(PFRTI(JL,4)*PEVAPTI(JL,4)-PEINTTI(JL,4))*ZROOTW(JL,1,4)/ZEXT(JL,4)
    ZSAWEXT(JL,1)=ZEXTK
    ZCONDS(JL)=ZCONDS(JL)+MAX((PFRTI(JL,4)*PEVAPTI(JL,4)-PEINTTI(JL,4)),0.0_JPRB)
    ZWSFL=ZWSFL+ZEXTK
!  Tile 6
    ZEXTK=(PFRTI(JL,6)*PEVAPTI(JL,6)-PEINTTI(JL,6))*ZROOTW(JL,1,6)/ZEXT(JL,6)
    ZSAWEXT(JL,1)=ZSAWEXT(JL,1)+ZEXTK
    ZCONDS(JL)=ZCONDS(JL)+MAX((PFRTI(JL,6)*PEVAPTI(JL,6)-PEINTTI(JL,6)),0.0_JPRB)
    ZWSFL=ZWSFL+ZEXTK
!  Tile 7
    ZEXTK=(PFRTI(JL,7)*(PEVAPTI(JL,7)-PEVAPSNW(JL))-PEINTTI(JL,7))&
     & *ZROOTW(JL,1,7)/ZEXT(JL,7)
    ZSAWEXT(JL,1)=ZSAWEXT(JL,1)+ZEXTK
    ZCONDS(JL)=ZCONDS(JL)+MAX((PFRTI(JL,7)*(PEVAPTI(JL,7)-PEVAPSNW(JL))-PEINTTI(JL,7)),0.0_JPRB)
    ZWSFL=ZWSFL+ZEXTK
!  Tile 8
    ZEXTK=PFRTI(JL,8)*PEVAPTI(JL,8)+PEMSSN(JL)-PEINTTI(JL,8)
    ZCONDS(JL)=ZCONDS(JL)+MAX((PFRTI(JL,8)*PEVAPTI(JL,8)+PEMSSN(JL)-PEINTTI(JL,8)),0.0_JPRB)
    ZWSFL=ZWSFL+ZEXTK
!  Tile 10
    IF ( KTILES .GT. 9 ) THEN
     ZEXTK=PFRTI(JL,10)*PEVAPTI(JL,10)+PEMSSN(JL)-PEINTTI(JL,10)
     ZCONDS(JL)=ZCONDS(JL)+MAX((PFRTI(JL,10)*PEVAPTI(JL,10)+PEMSSN(JL)-PEINTTI(JL,10)),0.0_JPRB)
     ZWSFL=ZWSFL+ZEXTK
    ENDIF

!           SOIL WATER DIFFUSIVITY
    PCFW(JL,1)=ZD*PTMST*RSIMP/(0.5_JPRB*(RDAW(1)+RDAW(2)))

!           RIGHT-HAND SIDE
    PRHSW(JL,1)=PTMST*(-ZK+ZWSFL*ZHOH2O)/RDAW(1)

!          BUDGETS AND RUN-OFF PROPERLY SCALED.
    PSAWGFL(JL,1)=RHOH2O*ZK
    PROS(JL)=PTMST*ZROT                                 !Units kg/m2
    PFWEL1(JL)=PFWEL1(JL)+&
     & (PFRTI(JL,3)*PEVAPTI(JL,3)-PEINTTI(JL,3)+ZSAWEXT(JL,1)+&
     & PFRTI(JL,8)*PEVAPTI(JL,8)-PEINTTI(JL,8))
    IF ( KTILES .GT. 9 ) THEN
     PFWEL1(JL)=PFWEL1(JL)+(PFRTI(JL,10)*PEVAPTI(JL,10)-PEINTTI(JL,10))
    ENDIF
  ELSE
!          SEA POINTS.
    PCFW(JL,1)=0.0_JPRB
    PRHSW(JL,1)=0.0_JPRB
    PROS(JL)=0.0_JPRB
    PSAWGFL(JL,1)=0.0_JPRB
    ZSAWEXT(JL,1)=0.0_JPRB
  ENDIF
ENDDO

!*          2.4 COEFFICIENTS FOR OTHER LAYERS

DO JK=2,KLEVS
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN

!           HYDRAULIC PROPERTIES.

!           VAN GENUCHTEN
      IF (JK < KLEVS) THEN
        JS=KSOTY(JL)
        ZW=MAX(MAX(PWSAM1M(JL,JK),PWSAM1M(JL,JK+1)),RWRESTM(JS))
      ELSE
        JS=KSOTY(JL)
        ZW=MAX(PWSAM1M(JL,JK),RWRESTM(JS))
      ENDIF
      ZDMAX=RDMAXM(JS)
      ZDMIN=RDMINM(JS)
      ZSE=(ZW-RWRESTM(JS))/(RWSATM(JS)-RWRESTM(JS))
      ZWFAC=RWSATM(JS)-RWRESTM(JS)
      ZLAM=RLAMBDAM(JS)
      ZMFAC=RMFACM(JS)
      ZRMFAC=1./ZMFAC
      ZWCONS=RWCONSM(JS)
      ZALPHA=RMVGALPHA(JS)
      IF (ZW.LE.(1.001*RWRESTM(JS))) THEN
      ZD=ZDMIN
      IF ( KTILES .GT. 9 ) THEN
       ZD=ZDMIN    * (1.0_JPRB-PCUR(JL))  + (PCUR(JL)*1.e-4_JPRD) !Test value for Urban
      ENDIF
        ZK=0.0_JPRB
      ELSEIF ((ZW.GT.(1.001_JPRB*RWRESTM(JS))).AND.(ZW.LE.(0.999_JPRB*RWSATM(JS)))) THEN
        ZKD=ZWCONS*ZSE**ZLAM* &
 &      (1.0_JPRB-((1.0_JPRB-(ZSE**ZRMFAC))**ZMFAC))**2.0_JPRB
        ! perhaps better to be rewritten
        ZDD=((1.0_JPRD-ZMFAC)*ZWCONS)/ &
 &      SIGN(MAX(ABS(ZALPHA*ZMFAC*ZWFAC),ZEPSILON),ZALPHA*ZMFAC*ZWFAC)
        ZDD= ZDD *REAL(ZSE**(ZLAM-ZRMFAC),JPRD) &
 &      *(((1.0_JPRD-(ZSE**ZRMFAC))**ZMFAC) &
 &      +((1.0_JPRD-(ZSE**ZRMFAC))**(-ZMFAC))-2.0_JPRD)
        ZK=ZKD
        ZD=ZDD
      ELSE
        ZK=RWCONSM(JS)
        ZD=ZDMAX
      ENDIF
      
      ZWM=MAX(MAX(PWSAM1M(JL,JK-1),PWSAM1M(JL,JK)),RWRESTM(JS))
      ZSE=(ZWM-RWRESTM(JS))/(RWSATM(JS)-RWRESTM(JS))
      IF (ZWM.LE.RWRESTM(JS)) THEN
        ZKM=0.0_JPRB
      ELSEIF ((ZWM.GT.RWRESTM(JS)).AND.(ZWM.LT.(0.999_JPRB*RWSATM(JS)))) THEN
        ZKMD=ZWCONS*ZSE**ZLAM*(1.0_JPRB-((1.0_JPRB-(ZSE**ZRMFAC))**ZMFAC))**2.0_JPRB
        ZKM=ZKMD
      ELSE 
        ZKM=RWCONSM(JS)
      ENDIF
          
      IF (JK == KLEVS) ZD=0.0_JPRB 

      IF (LLFREEZ) THEN
        IF (JK < KLEVS) THEN
          ZFF=MIN(ZF(JL,JK),ZF(JL,JK+1))
          ZD=ZFF*ZDMIN+(1.0_JPRB-ZFF)*ZD
        ELSE
          ZFF=ZF(JL,JK)
        ENDIF
        ZFFM=MIN(ZF(JL,JK-1),ZF(JL,JK))
        ZK=ZFF*0.0_JPRB+(1.0_JPRB-ZFF)*ZK
        ZKM=ZFFM*0.0_JPRB+(1.0_JPRB-ZFFM)*ZKM
      ENDIF

!*             TILE CONTRIBUTIONS
!              ------------------
!  Tile by tile contribution of root extraction to the r.h.s.
!  Tile 4
      ZEXTK=(PFRTI(JL,4)*PEVAPTI(JL,4)-PEINTTI(JL,4))*&
       & ZROOTW(JL,JK,4)/ZEXT(JL,4)  
      ZSAWEXT(JL,JK)=ZEXTK
      ZWSFL=ZEXTK
!  Tile 6
      ZEXTK=(PFRTI(JL,6)*PEVAPTI(JL,6)-PEINTTI(JL,6))*&
       & ZROOTW(JL,JK,6)/ZEXT(JL,6)  
      ZSAWEXT(JL,JK)=ZSAWEXT(JL,JK)+ZEXTK
      ZWSFL=ZWSFL+ZEXTK
!  Tile 7
      ZEXTK=(PFRTI(JL,7)*(PEVAPTI(JL,7)-PEVAPSNW(JL))-PEINTTI(JL,7))*&
       & ZROOTW(JL,JK,7)/ZEXT(JL,7)  
      ZSAWEXT(JL,JK)=ZSAWEXT(JL,JK)+ZEXTK
      ZWSFL=ZWSFL+ZEXTK

!           SOIL WATER DIFFUSIVITY
      IF (JK < KLEVS) THEN
        PCFW(JL,JK)=ZD*PTMST*RSIMP/(0.5_JPRB*(RDAW(JK)+RDAW(JK+1)))
      ELSE
        PCFW(JL,JK)=0.0_JPRB
      ENDIF

!           RIGHT-HAND SIDE
      PRHSW(JL,JK)=PTMST*(ZKM-ZK+ZWSFL*ZHOH2O)/RDAW(JK)

!          BUDGETS AND RUN-OFF PROPERLY SCALED.
      PSAWGFL(JL,JK)=RHOH2O*ZK
      PFWE234(JL)=PFWE234(JL)+ZSAWEXT(JL,JK)

    ELSE
!          SEA POINTS.
      PCFW(JL,JK)=0.0_JPRB
      PRHSW(JL,JK)=0.0_JPRB
      PSAWGFL(JL,JK)=0.0_JPRB
      ZSAWEXT(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!*          3. DDH diagnostics
!              ---------------
IF (SIZE(PDHWLS) > 0) THEN
! Soil water ans soil ice water
  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      PDHWLS(JL,JK,1)=RHOH2O*RDAW(JK)*PWSAM1M(JL,JK)
      PDHWLS(JL,JK,2)=ZF(JL,JK)*PDHWLS(JL,JK,1)
    ENDDO
  ENDDO
! Large-scale throughfall
  DO JL=KIDIA,KFDIA
    PDHWLS(JL,1,3)=PTSFL(JL)
! Convective throughfall
    PDHWLS(JL,1,4)=PTSFC(JL)
! Melt throughfall
    PDHWLS(JL,1,5)=PMSN(JL)
  ENDDO

! zero out fluxes
  DO JK=2,KLEVS
    DO JL=KIDIA,KFDIA
      PDHWLS(JL,JK,3)=0.0_JPRB
      PDHWLS(JL,JK,4)=0.0_JPRB
      PDHWLS(JL,JK,5)=0.0_JPRB
      PDHWLS(JL,JK,10)=0.0_JPRB
    ENDDO
  ENDDO

! Runoff top layer (negative values mean water lost by the layer)
  DO JK=1,KLEVS
    DO JL=KIDIA,KFDIA
      PDHWLS(JL,JK,6)=0.0_JPRB
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    PDHWLS(JL,1,6)=PDHWLS(JL,1,6)-ZINVTMST*PROS(JL)  
  ENDDO

! Root extraction (<0) without condensation on tile 4 6 7 8
  DO JL=KIDIA,KFDIA
    PDHWLS(JL,1,7)=ZSAWEXT(JL,1)-ZCONDS(JL)
    DO JK=2,KLEVS
      PDHWLS(JL,JK,7)=ZSAWEXT(JL,JK)
    ENDDO
  ENDDO
! Bare ground evaporation including mismatches from snow and interception layer
  DO JL=KIDIA,KFDIA
    PDHWLS(JL,1,9)=PFRTI(JL,8)*PEVAPTI(JL,8)-&
      & PEINTTI(JL,8)+&
      & PEMSSN(JL)+&
      & PFRTI(JL,3)*PEVAPTI(JL,3)-&
      & PEINTTI(JL,3)
    DO JK=2,KLEVS
      PDHWLS(JL,JK,9)=0.0_JPRB
    ENDDO
  ENDDO
! Condensation (>0) due to excess dew deposition on tile 3 (clipped to WLmax)
  DO JL=KIDIA,KFDIA
    PDHWLS(JL,1,10)=ZCONDS(JL)
  ENDDO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFWEXC_VG_MOD:SRFWEXC_VG',1,ZHOOK_HANDLE)

END SUBROUTINE SRFWEXC_VG
END MODULE SRFWEXC_VG_MOD
