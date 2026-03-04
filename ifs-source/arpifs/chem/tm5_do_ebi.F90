! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_DO_EBI(YGFL,KIDIA,KFDIA,KLON,KMAXIT,PDT,PRR,PRJ,PY0,PY,PVD,PRES)

!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!--------------------------------------------------------------------------
!   Eulerian backward Iteration
!   Chemistry solver for the CBM4 scheme 
!--------------------------------------------------------------------------
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5_ebi* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! PDT   :  Time step in seconds 
! KMAXIT :  maximum number of iterations, depending on model level 
! PRR  (KLON,NREAC)      : reaction rates
! PRJ  (KLON,NPHOTO)     : photolysis rates
! PY0(KLON,NCHEM)        : initial volume ratios OF TRACERS           (mol/mol)
! PVD  (KLON,NCHEM)      : dry deposition velocities (may be equal to zero)
! PRES (KLON)            : FULL-LEVEL PRESSURE           (Pa)
!
!
! OUTPUTS:
! -------
! PY (KLON,NCHEM+3)        : final   volume ratios OF TRACERS           (mol/mol)
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!        TM5-community    
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2009-09-08



USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE TM5_CHEM_MODULE, ONLY :  IO3 , IH2O2, ICH4, ICO,&
   &   IHNO3,  ICH3O2H,  ICH2O,  IPAR, IETH,  IOLE,  IALD2  , IPAN,&
   &   IROOH,  IORGNTR, IISOP,   ISO2,   IDMS,   INH3,   ISO4,  INH4,&
   &   IMSA,   IMGLY, IRN222,  IPB210,   INO,   IHO2,   ICH3O2,&
   &   IOH,    INO2,   INO3,   IN2O5,   IHNO4,  IC2O3,   IROR,   IRXPAR,&
   &   IXO2,   IXO2N,   INH2,   ICH3OH,   IHCOOH,  IMCOOH,   IC2H6,&
   &   IETHOH,   IC3H8,   IC3H6,  ITERP,   IISPD,  INO3_A, IACO2,IACET,&
   &   IIC3H7O2,IHYPROPO2,&
!* reaction rates
   &   KNOO3,  KHO2NO,  KMO2NO,  KNO2OH,&
   &   KOHHNO3,  KNO2O3,  KNONO3 ,  KNO2NO3 ,  KN2O5,  KHNO4OH ,  KNO2HO2,  KHNO4M  ,&
   &   KO3HO2,  KCOOH,  KO3OH,  KHPOH,  KFRMOH,  KCH4OH  ,&
   &   KOHMPER,  KOHROOH,  KMO2HO2A,KMO2HO2B ,  KMO2MO2 ,  KHO2OH  ,  KHO2HO2,  KN2O5L,&
   &   KN2O5_AER,KHO2_AER,KNO3_AER, KHO2L,&
   &   KH2OH,  KC41,  KC43,  KC44,  KC46,  KC47,  KC48,  KC49,&
   &   KC50A, KC50B,  KC52,  KC53,  KC54,  KC57,  KC58,  KC59,  KC61,&
   &   KC62,  KC73,  KC76,  KC77,  KC78,  KC79,  KC80,  KC81,&
   &   KC82,  KC83,  KC84,  KC85,  KDMSOHA ,  KDMSOHB,  KDMSNO3 ,  KSO2OH  ,&
   &   KNH3SO4 ,  KNH3OH ,  KNH2NO  ,  KNH2NO2 ,  KNH2HO2 ,  KNH2O2 ,  KOHCH3OH,&
   &   KOHHCOOH,  KNO3HO2 ,  KNO3MO2 ,  KNO3C2O3,  KNO3XO2 ,  KOHMCOOH,  KOHC2H6 ,&
   &   KOHETHOH,  KOHC3H8 ,  KOHC3H6 ,  KO3C3H6 ,  KNO3C3H6,&
   &   KOHTERP ,  KO3TERP ,  KNO3TERP,  KRN222  , KO3PO3 ,&
   &   KOHACET,KACO2HO2,KACO2MO2,KACO2NO,KACO2XO2,&
   &   KXO2XO2N,KXO2N,&
   &   KOHISPD ,  KO3ISPD ,  KNO3ISPD,&
  &  KNOHYPROPO2,KHO2HYPROPO2,KNOIC3H7O2,KHO2IC3H7O2, NREAC
USE TM5_PHOTOLYSIS , ONLY : NPHOTO,&
  !* photolysis rates
   & JO3D  ,JNO2  ,JH2O2 ,JHNO3 ,JHNO4 ,&
   & JN2O5 ,JACH2O,JBCH2O,JMEPE ,JANO3 ,JBNO3 ,JPANA, JORGN ,J45   ,J74   ,&
   & JROOH ,JO2   ,JISPD ,JA_ACET,JB_ACET  

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KMAXIT
REAL(KIND=JPRB),INTENT(IN)    :: PDT
REAL(KIND=JPRB),INTENT(IN)    :: PRR(KLON,NREAC)
REAL(KIND=JPRB),INTENT(IN)    :: PRJ(KLON,NPHOTO)   
REAL(KIND=JPRB),INTENT(IN)    :: PY0(KLON,YGFL%NCHEM)   ! initial concentrations
REAL(KIND=JPRB),INTENT(IN)    :: PVD(KLON,YGFL%NCHEM)   ! deposition velocities ... in case we apply them here...
REAL(KIND=JPRB),INTENT(INOUT) :: PY(KLON,YGFL%NCHEM+3)  ! final concentrations
REAL(KIND=JPRB),INTENT(IN)    :: PRES(KLON)        ! full level pressure

! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

INTEGER(KIND=JPIM)    :: ITER
REAL (KIND=JPRB)      ::  ZR57, ZR89, ZP1, ZR12, ZR21, ZXL1, ZP2, ZXL2, ZP3, ZP32,&
                      &   ZXL3, ZX1, ZX2, ZX3, ZC1, ZC2, ZC3, ZY2, ZXJT, ZR21T,   &
                      &   ZR12T, ZACUB, ZBCUB, ZCCUB, ZCUBDET, ZDNO2, ZR56,    &
                      &   ZR65, ZR75, ZP5, ZXL5, ZR66, ZX5, ZP6, ZXL6, ZX6, ZC6,  &
                      &   ZXL7, ZY1, ZC7, ZR98, ZP8, ZXL8, ZX4, ZC5, ZXL9,       &
                      &   ZR1920, ZR1919, ZP19, ZXL19, ZR2019, ZXL20, ZXLHNO3, &
                      &   ZPH2O2, ZXLH2O2, ZPCH2O, ZPCO, ZPHNO3, ZXLCH2O,     &
                      &   ZPCH3O2, ZXLCH3O2, ZPCH3O2H, ZXLCH3O2H, ZPALD2,    &
                      &   ZXLALD2, ZPMGLY, ZXLMGLY, ZPOLE, ZXLETH, ZXLOLE,    &
                      &   ZXLISOP, ZPRXPAR, ZXLRXPAR, ZPPAR, ZXLPAR, ZPROR,   &
                      &   ZXLROR, ZPXO2, ZXLXO2, ZPXO2N, ZXLXO2N, ZPROOH, ZXLTERP, &
                      &   ZPHCOOH,ZPMCOOH,ZXLC3H6, &
                      &   ZXLROOH, ZPORGNTR, ZXLORGNTR, ZXLCO, ZQDMS, ZPSO2, &
                      &   ZQSO2, ZQSO2D, ZQNH3, ZPNH2, ZQNH2, ZQDMS1, ZQDMS2,  &
                      &   ZPMSA, ZPISPD,ZXLISPD,ZXLACET,ZPACET,ZPACO2,ZXLACO2, &
                      &   ZPCH3OH, ZPIC3H7O2,ZXLIC3H7O2,ZPHYPROPO2,ZXLHYPROPO2
REAL(KIND=JPRB)       :: ZDT2,ZY0_IACID,ZY_IACID
! * counters
INTEGER(KIND=JPIM) :: JL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_DO_EBI',0,ZHOOK_HANDLE )



ZDT2=PDT*PDT


      DO ITER=1        ,KMAXIT


         DO JL=KIDIA,KFDIA
            ! --- Short living compounds & groups
            ! --- First group: NO NO2 O3
            ZP1=PRJ(JL,JBNO3)*PY(JL,INO3) ! NOx (NO) emissions are removed here:  ... +emino(JL) ...
            ZR12=0._JPRB
            ZR21=PRR(JL,KHO2NO)*PY(JL,IHO2)+PRR(JL,KMO2NO)*PY(JL,ICH3O2)&
              &   +PRR(JL,KC79)*PY(JL,IXO2)+PRR(JL,KC46)*PY(JL,IC2O3)&
              &   +PRR(JL,KACO2NO)*PY(JL,IACO2)&
              &   +PRR(JL,KNOIC3H7O2)*PY(JL,IIC3H7O2)&
              &   +PRR(JL,KNOHYPROPO2)*PY(JL,IHYPROPO2)
            ZXL1=PRR(JL,KNONO3)*PY(JL,INO3)+PRR(JL,KC81)*PY(JL,IXO2N)
            ZXL1 = ZXL1 + PVD(JL,INO)
            ZP2=PRJ(JL,JHNO3)*PY(JL,IHNO3)+PRJ(JL,JN2O5)*PY(JL,IN2O5)&
              &   +PRR(JL,KN2O5)*PY(JL,IN2O5)+PRJ(JL,JANO3)*PY(JL,INO3)&
              &   +PY(JL,IHNO4)*(PRJ(JL,JHNO4)+PRR(JL,KHNO4M)+PRR(JL,KHNO4OH)&
              &   *PY(JL,IOH))+2.*PRR(JL,KNONO3)*PY(JL,INO3)*PY(JL,INO)&
              &   +PRR(JL,KC48)*PY(JL,IPAN)+PRR(JL,KC59)*PY(JL,IOLE)*PY(JL,INO3)&
              &   +PRJ(JL,JORGN)*PY(JL,IORGNTR)&
              &   +PRJ(JL,JPANA)*PY(JL,IPAN)&
              &   +0.2*PRR(JL,KC78) *PY(JL,IISOP)*PY(JL,INO3)&
              &   +PRR(JL,KNO3MO2)*PY(JL,ICH3O2)*PY(JL,INO3)&
              &   +PRR(JL,KNO3C2O3)*PY(JL,IC2O3)*PY(JL,INO3)&
              &   +PRR(JL,KNO3XO2)*PY(JL,IXO2)*PY(JL,INO3)&
              &   +0.47*PRR(JL,KNO3TERP)*PY(JL,INO3)*PY(JL,ITERP)
            ZXL2=PRR(JL,KNO2OH)*PY(JL,IOH)+PRR(JL,KNO2NO3)*PY(JL,INO3)&
              &   +PRR(JL,KNO2HO2)*PY(JL,IHO2)+PRR(JL,KNO2O3)*PY(JL,IO3)&
              &   +PRR(JL,KC47)*PY(JL,IC2O3)
            ZXL2 = ZXL2 + PVD(JL,INO2)
            ZP3=PRJ(JL,JANO3)*PY(JL,INO3)+PRJ(JL,JO2) ! * O3P + O2 = > O3
            ZXL3=PRR(JL,KO3HO2)*PY(JL,IHO2)+PRR(JL,KO3OH)*PY(JL,IOH)&
              &  +PRR(JL,KNO2O3)*PY(JL,INO2)+PRJ(JL,JO3D)+PRR(JL,KO3PO3)&
              &  +PRR(JL,KC58)*PY(JL,IOLE)&
              &  +PRR(JL,KC62)*PY(JL,IETH)&
              &  +PRR(JL,KC77)*PY(JL,IISOP)&
              &  +PRR(JL,KO3C3H6)*PY(JL,IC3H6)&
              &  +PRR(JL,KO3TERP)*PY(JL,ITERP)&
              &  +PRR(JL,KO3ISPD)*PY(JL,IISPD)
            ZXL3 = ZXL3 + PVD(JL,IO3)

            ZX1=PY0(JL,INO)+ZP1*PDT
            ZX2=PY0(JL,INO2)+ZP2*PDT
            ZX3=PY0(JL,IO3)+ZP3*PDT
            ZC1=1._JPRB+ZXL1*PDT
            ZC2=1._JPRB+ZXL2*PDT
            ZC3=1._JPRB+ZXL3*PDT
            ZY1=PRR(JL,KNOO3)*PDT
            ZR21T=ZR21*PDT
            ZR12T=ZR12*PDT
            ZXJT=PRJ(JL,JNO2)*PDT
            ZP32=PRR(JL,KC50B)*PY(JL,IC2O3)*PY(JL,IHO2)
            ! --- Solve unknown x
            ZACUB=-2._JPRB*ZY1*(ZC2+ZR12T+ZC2*ZR21T/ZC1)
            ZBCUB=2._JPRB*ZC1*ZC2*ZC3+2._JPRB*ZC1*ZC3*(ZR12T+ZXJT)+2._JPRB*ZC2*ZC3*ZR21T+&
               &  2._JPRB*ZY1*(ZR12T*(ZX1-ZX2)+2.*ZC2*ZR21T*ZX1/ZC1+ZC2*(ZX1+ZX3))
            ZCCUB=2._JPRB*ZC1*ZC3*ZX2*(ZR12T+ZXJT)-2._JPRB*ZC2*ZC3*ZX1*ZR21T+2._JPRB*ZY1*ZX1*&
               &  (ZX2*ZR12T-ZC2*ZX3-ZC2*ZR21T*ZX1/ZC1)
            ZCUBDET=MAX(0.0_JPRB,ZBCUB*ZBCUB-4.*ZACUB*ZCCUB)
            ZDNO2=(-1._JPRB*ZBCUB+SQRT(ZCUBDET))/(2._JPRB*ZACUB)
            ZDNO2=MIN(ZX1,ZDNO2)
            PY(JL,INO2)=(ZX2+ZDNO2)/ZC2
            PY(JL,INO)=(ZX1-ZDNO2)/ZC1
            PY(JL,IO3)=(ZX3+(ZP32*PDT)+ZXJT*PY(JL,INO2))/(ZC3+ZY1*PY(JL,INO))
            ! --- Second group: PY(JL,iho2) PY(JL,ioh) PY(JL,ihno4)
            ZR57=PRJ(JL,JHNO4)+PRR(JL,KHNO4M)
            ZR56=PRR(JL,KCOOH)*PY(JL,ICO)+PRR(JL,KO3OH)*PY(JL,IO3)+PRR(JL,KHPOH)&
               &  *PY(JL,IH2O2)+PRR(JL,KFRMOH)*PY(JL,ICH2O)+PRR(JL,KH2OH)&
               &  +PRR(JL,KSO2OH)*PY(JL,ISO2) 

            ZP5=2._JPRB*PRJ(JL,JBCH2O)*PY(JL,ICH2O)&
              &   +PRR(JL,KMO2NO)*PY(JL,ICH3O2)*PY(JL,INO)&
              &   +PRJ(JL,JMEPE)*PY(JL,ICH3O2H)&
              &   +2.0_JPRB*PRJ(JL,J45)*PY(JL,IALD2)&
              &   +PRJ(JL,J74)*PY(JL,IMGLY)+0.11_JPRB*PRR(JL,KC52)*PY(JL,IPAR)*PY(JL,IOH)&
              &   +0.94_JPRB*PRR(JL,KC53)*PY(JL,IROR)+PRR(JL,KC54)*PY(JL,IROR)&
              &   +1.57_JPRB*PRR(JL,KC57)*PY(JL,IOLE)*PY(JL,IOH)&
              &   +0.76_JPRB*PRR(JL,KC58)*PY(JL,IO3)*PY(JL,IOLE)&
              &   +0.56_JPRB*PRR(JL,KC59)*PY(JL,INO3)*PY(JL,IOLE)&
              &   +PRR(JL,KC61)*PY(JL,IETH)*PY(JL,IOH)+0.22_JPRB*PRR(JL,KC62)&
              &   *PY(JL,IETH)*PY(JL,IO3)&
              &   +0.066_JPRB*PRR(JL,KC77)*PY(JL,IISOP)*PY(JL,IO3)&
              &   +PRR(JL,KC41)*PY(JL,ICH2O)*PY(JL,INO3)&
              &   +0.9_JPRB*PRR(JL,KC84)*PY(JL,IOH)*PY(JL,IORGNTR)&
              &   +0.9_JPRB*PRJ(JL,JORGN)*PY(JL,IORGNTR)&
              &   +0.8_JPRB*PRR(JL,KC78)*PY(JL,IISOP)*PY(JL,INO3)&
              &   +0.74_JPRB*PRR(JL,KMO2MO2)*PY(JL,ICH3O2)*PY(JL,ICH3O2)&
              &   +0.9_JPRB*PRJ(JL,JROOH)*PY(JL,IROOH)&
              &   +PRR(JL,KNO3MO2)*PY(JL,ICH3O2)*PY(JL,INO3)&
              &   +0.19_JPRB*PRR(JL,KO3C3H6)*PY(JL,IO3)*PY(JL,IC3H6)&
              &   +0.28_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
              &   +0.75_JPRB*PRR(JL,KNO3TERP)*PY(JL,INO3)*PY(JL,ITERP)&
              &   +0.154_JPRB*PRR(JL,KO3ISPD)*PY(JL,IO3 ) *PY(JL,IISPD)&
              &   +0.925_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3 )*PY(JL,IISPD)&
              &   +1.033_JPRB*PRJ(JL,JISPD)*PY(JL,IISPD)&
              &   +PRR(JL,KOHCH3OH)*PY(JL,ICH3OH)&
              &   +PRR(JL,KOHHCOOH)*PY(JL,IHCOOH)&
              &   +PRR(JL,KOHETHOH)*PY(JL,IETHOH)&
              &   +1.22*PRR(JL,KOHTERP)*PY(JL,ITERP)&
              &   +PRR(JL,KOHC2H6)*PY(JL,IC2H6)&
              &   +0.5*PRR(JL,KC76)*PY(JL,IOH)*PY(JL,IISOP)&! reduced according to archibald et al, AE, 2011
              &   +0.503*PRR(JL,KOHISPD)*PY(JL,IOH)*PY(JL,IISPD)&
              &   +0.5*PRR(JL,KACO2MO2)*PY(JL,IACO2)*PY(JL,ICH3O2)&
              &   +PRR(JL,KACO2NO) *PY(JL,IACO2)*PY(JL,INO)&
              &   +PRR(JL,KNOIC3H7O2)*PY(JL,IIC3H7O2)*PY(JL,INO)&
              &   +PRR(JL,KNOHYPROPO2)*PY(JL,IHYPROPO2)*PY(JL,INO)
            ZXL5=PRR(JL,KHO2NO)   *PY(JL,INO)&
              &  +PRR(JL,KNO2HO2)*PY(JL,INO2)&
              &  +PRR(JL,KO3HO2) *PY(JL,IO3)&
              &  +(PRR(JL,KMO2HO2A)+PRR(JL,KMO2HO2B))*PY(JL,ICH3O2)&
              &  +(PRR(JL,KC50A)+PRR(JL,KC50B)) *PY(JL,IC2O3)&
              &  +PRR(JL,KHO2OH) *PY(JL,IOH)&
              &  +PRR(JL,KC82)   *PY(JL,IXO2)&
              &  +PRR(JL,KC85)   *PY(JL,IXO2N)&
              &  +PRR(JL,KNO3HO2)*PY(JL,INO3)&
              &  +PRR(JL,KACO2HO2)*PY(JL,IACO2)&
              &  +PRR(JL,KHO2IC3H7O2)*PY(JL,IIC3H7O2)&
              &  +PRR(JL,KHO2HYPROPO2)*PY(JL,IHYPROPO2)&
              &  +PRR(JL,KHO2_AER)&
              &  +PRR(JL,KHO2L)
            ZR66=2._JPRB*PRR(JL,KHO2HO2)
            ZX5=PY0(JL,IHO2)+ZP5*PDT
            ZR65=PRR(JL,KHO2NO)*PY(JL,INO)+PRR(JL,KO3HO2)*PY(JL,IO3)


            ZP6=PRJ(JL,JHNO3)*PY(JL,IHNO3)&
             &  +2._JPRB*PRJ(JL,JO3D)*PY(JL,IO3)&
             &  +2._JPRB*PRJ(JL,JH2O2)*PY(JL,IH2O2)&
             &  +PRJ(JL,JMEPE)*PY(JL,ICH3O2H)&
             &  +0.1_JPRB*PRR(JL,KC58)*PY(JL,IO3)*PY(JL,IOLE)&
             &  +0.266_JPRB*PRR(JL,KC77)*PY(JL,IISOP)*PY(JL,IO3)&
             &  +PRJ(JL,JROOH)*PY(JL,IROOH)&
             &  +0.12_JPRB*PRR(JL,KC62)*PY(JL,IETH)*PY(JL,IO3)&
             &  +0.33_JPRB*PRR(JL,KO3C3H6)*PY(JL,IO3)*PY(JL,IC3H6)&
             &  +0.57_JPRB*PRR(JL,KO3TERP)*PY(JL,ITERP)*PY(JL,IO3)&
             &  +0.268_JPRB*PRR(JL,KO3ISPD)*PY(JL,IISPD)*PY(JL,IO3)&
             &  +PRR(JL,KC50A)*PY(JL,IC2O3)*PY(JL,IHO2)

            ZXL6=PRR(JL,KHNO4OH)*PY(JL,IHNO4)&
             &  +PRR(JL,KHO2OH)*PY(JL,IHO2)&
             &  +PRR(JL,KNO2OH)*PY(JL,INO2)&
             &  +PRR(JL,KOHHNO3)*PY(JL,IHNO3)&
             &  +PRR(JL,KCOOH)*PY(JL,ICO)&
             &  +PRR(JL,KO3OH)*PY(JL,IO3)&
             &  +PRR(JL,KHPOH)*PY(JL,IH2O2)&
             &  +PRR(JL,KFRMOH)*PY(JL,ICH2O)&
             &  +PRR(JL,KCH4OH)*PY(JL,ICH4)&
             &  +0.6_JPRB*PRR(JL,KOHMPER)*PY(JL,ICH3O2H)&
             &  +PRR(JL,KC43)*PY(JL,IALD2)&
             &  +PRR(JL,KC73)*PY(JL,IMGLY)&
             &  +PRR(JL,KC52)*PY(JL,IPAR)&
             &  +PRR(JL,KC57)*PY(JL,IOLE)&
             &  +PRR(JL,KC61)*PY(JL,IETH)&
             &  +PRR(JL,KC76)*PY(JL,IISOP)&
             &  +0.77_JPRB*PRR(JL,KOHROOH)*PY(JL,IROOH)&! note: change from '0.7' to '0.77'
             &  +PRR(JL,KC84)*PY(JL,IORGNTR)&
             &  +PRR(JL,KH2OH)&
             &  +PRR(JL,KSO2OH)*PY(JL,ISO2)&! bug found by Jason 01/2008
             &  +(PRR(JL,KDMSOHA)+PRR(JL,KDMSOHB)) *PY(JL,IDMS)&
             &  +PRR(JL,KNH3OH)*PY(JL,INH3)&!sulfur
             &  +PRR(JL,KOHCH3OH)*PY(JL,ICH3OH)&
             &  +PRR(JL,KOHHCOOH)*PY(JL,IHCOOH)&
             &  +PRR(JL,KOHETHOH)*PY(JL,IETHOH)&
             &  +PRR(JL,KOHTERP)*PY(JL,ITERP)&
             &  +PRR(JL,KOHISPD)*PY(JL,IISPD)&
             &  +PRR(JL,KOHMCOOH)*PY(JL,IMCOOH)&
             &  +PRR(JL,KOHC2H6)*PY(JL,IC2H6)&
             &  +PRR(JL,KOHC3H8)*PY(JL,IC3H8)&
             &  +PRR(JL,KOHC3H6)*PY(JL,IC3H6)&
             &  +PRR(JL,KOHACET)*PY(JL,IACET)  
            ZX6=PY0(JL,IOH)+ZP6*PDT
            ZC6=1._JPRB+ZXL6*PDT
            ZR75=PRR(JL,KNO2HO2)*PY(JL,INO2)
            ZXL7=PRJ(JL,JHNO4)+PRR(JL,KHNO4OH)*PY(JL,IOH)+PRR(JL,KHNO4M)
            ZXL7 = ZXL7 + PVD(JL,IHNO4)
            ZC7=1._JPRB+ZXL7*PDT
            ZY1=ZR57/ZC7
            ZY2=ZR56/ZC6
            ZACUB=ZR66*PDT
            ZBCUB=1._JPRB+ZXL5*PDT-ZDT2*(ZY1*ZR75+ZY2*ZR65)
            ZCCUB=-1._JPRB*ZX5-PDT*(ZY1*PY0(JL,IHNO4)+ZY2*ZX6)
            ZCUBDET=ZBCUB*ZBCUB-4._JPRB*ZACUB*ZCCUB
            ZCUBDET=MAX(ZCUBDET,1.E-20_JPRB)
            PY(JL,IHO2)=MAX(0.1_JPRB,(-1._JPRB*ZBCUB+SQRT(ZCUBDET))/(2._JPRB*ZACUB))
            PY(JL,IOH)=(ZX6+ZR65*PY(JL,IHO2)*PDT)/ZC6
            PY(JL,IHNO4)=(PY0(JL,IHNO4)+ZR75*PDT*PY(JL,IHO2))/ZC7
            !
            ! --- Third group: NO3 N2O5
            !
            ZR89=PRJ(JL,JN2O5)+PRR(JL,KN2O5)
            ZP8=PRR(JL,KOHHNO3)*PY(JL,IHNO3)*PY(JL,IOH)+PRR(JL,KNO2O3)*PY(JL,INO2)*PY(JL,IO3)
            ZXL8=PRJ(JL,JBNO3)+PRJ(JL,JANO3)&
             &    +PRR(JL,KNONO3)*PY(JL,INO)&
             &    +PRR(JL,KNO2NO3)*PY(JL,INO2)&
             &    +PRR(JL,KC44)*PY(JL,IALD2)&
             &    +PRR(JL,KC59)*PY(JL,IOLE)&
             &    +PRR(JL,KC78)*PY(JL,IISOP)&
             &    +PRR(JL,KC41)*PY(JL,ICH2O)&
             &    +PRR(JL,KDMSNO3)*PY(JL,IDMS)&
             &    +PRR(JL,KNO3HO2)*PY(JL,IHO2)&
             &    +PRR(JL,KNO3MO2)*PY(JL,ICH3O2)&
             &    +PRR(JL,KNO3C2O3)*PY(JL,IC2O3)&
             &    +PRR(JL,KNO3XO2)*PY(JL,IXO2)&
             &    +PRR(JL,KNO3C3H6)*PY(JL,IC3H6)&
             &    +PRR(JL,KNO3TERP)*PY(JL,ITERP)&
             &    +PRR(JL,KNO3ISPD)*PY(JL,IISPD)&
             &    +PRR(JL,KNO3_AER)
            ZXL8 = ZXL8 + PVD(JL,INO3) 
            ZX4=PY0(JL,INO3)+ZP8*PDT
            ZC5=1._JPRB+ZXL8*PDT
            ZR98=PRR(JL,KNO2NO3)*PY(JL,INO2)
            ! VH: TM5- precalculate as first order N2O5 loss, this is now done directly in the ebi-solver
            ! PRR(i,j,kn2o5aq) = PRR(i,j,kn2o5aq)*(PY0(JL,iso4) + PY0(JL,imsa)  + PY0(JL,INO3_A) )   
            ZXL9=PRJ(JL,JN2O5)+PRR(JL,KN2O5)&
              &   +PRR(JL,KN2O5L)&!cmk rates now idependent from y
              &   +PRR(JL,KN2O5_AER) !VH - inclusion of heterogeneous reaction on aerosol
            ZXL9 = ZXL9 + PVD(JL,IN2O5)
            ZC6=1._JPRB+ZXL9*PDT
            ZC7=(ZC5*ZC6-ZR89*ZR98*ZDT2)
            PY(JL,IN2O5)=(ZC5*PY0(JL,IN2O5)+ZR98*PDT*ZX4)/ZC7
            PY(JL,INO3)=(ZC6*ZX4+ZR89*PDT*PY0(JL,IN2O5))/ZC7
            !
            ! --- Fourth group: C2O3 PAN
            !
            ZR1920=PRR(JL,KC48)+PRJ(JL,JPANA)
            ZR1919=PRR(JL,KC49)

            ZP19=PRR(JL,KC43)*PY(JL,IALD2)*PY(JL,IOH)&
              &   +PRR(JL,KC44)*PY(JL,IALD2)*PY(JL,INO3)&
              &   +PRR(JL,KC73)*PY(JL,IMGLY)*PY(JL,IOH)&
              &   +PRJ(JL,J74)*PY(JL,IMGLY)&
              &   +0.2_JPRB*PRR(JL,KC77)*PY(JL,IISOP)*PY(JL,IO3)&
              &   +0.74_JPRB*PRJ(JL,JROOH)*PY(JL,IROOH)&
              &   +0.74_JPRB*PRJ(JL,JORGN)*PY(JL,IORGNTR)&
              &   +0.74_JPRB*PRR(JL,KC84)*PY(JL,IOH)*PY(JL,IORGNTR)&
              &   +0.39_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
              &   +0.498_JPRB*PRR(JL,KOHISPD)*PY(JL,IOH)*PY(JL,IISPD)&
              &   +0.114_JPRB*PRR(JL,KO3ISPD)*PY(JL,IO3)*PY(JL,IISPD)&
              &   +0.075_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)&
              &   +0.967_JPRB*PRJ(JL,JISPD)*PY(JL,IISPD)&
              &   +0.3*PRR(JL,KACO2MO2)*PY(JL,IACO2)*PY(JL,ICH3O2)&
              &   +PRR(JL,KACO2NO)*PY(JL,IACO2)*PY(JL,INO)&
              &   +PRJ(JL,JB_ACET)*PY(JL,IACET)      
            ZXL19=PRR(JL,KC46)*PY(JL,INO)+(PRR(JL,KC50A)+PRR(JL,KC50B))*PY(JL,IHO2)&
              &   +PRR(JL,KC47)*PY(JL,INO2)&
              &   +PRR(JL,KNO3C2O3)*PY(JL,INO3)
            ZXL19 = ZXL19 + PVD(JL,IC2O3)
            ZR2019=PRR(JL,KC47)*PY(JL,INO2)
            ZXL20=PRR(JL,KC48)+PRJ(JL,JPANA)
            ZXL20 = ZXL20 + PVD(JL,IPAN)
            ZACUB=2._JPRB*ZR1919*PDT*(1._JPRB+ZXL20*PDT)
            ZBCUB=(1._JPRB+ZXL20*PDT)*(1+ZXL19*PDT)-ZR1920*PDT*ZR2019*PDT
            ZCCUB=(1._JPRB+ZXL20*PDT)*(PY0(JL,IC2O3)+ZP19*PDT)+ZR1920*PDT*PY0(JL,IPAN)
            ZCUBDET=ZBCUB*ZBCUB+4._JPRB*ZACUB*ZCCUB
            PY(JL,IC2O3)=MAX(1E-8_JPRB,(-1._JPRB*ZBCUB+SQRT(MAX(0.0_JPRB,ZCUBDET)))/(2._JPRB*ZACUB))    !cmk  put max here....
            PY(JL,IPAN)=(PY0(JL,IPAN)+ZR2019*PY(JL,IC2O3)*PDT)/(1._JPRB+ZXL20*PDT)
            !
            ! --- CH4 chemistry (short living radicals)
            !
            ZPCH3O2=PRR(JL,KCH4OH)*PY(JL,ICH4)*PY(JL,IOH)&
             &   +0.6_JPRB*PRR(JL,KOHMPER)*PY(JL,IOH)*PY(JL,ICH3O2H)&
             &   +PRR(JL,KNO3C2O3)*PY(JL,INO3)*PY(JL,IC2O3)&
             &   +PRR(JL,KC46)*PY(JL,INO)*PY(JL,IC2O3)&
             &   +2.0_JPRB*PRR(JL,KC49)*PY(JL,IC2O3)*PY(JL,IC2O3)&
             &   +PRJ(JL,J45)*PY(JL,IALD2)&
             &   +0.74_JPRB*PRJ(JL,JROOH)*PY(JL,IROOH)&
             &   +0.74_JPRB*PRR(JL,KC84)*PY(JL,IORGNTR)*PY(JL,IOH)&
             &   +0.74_JPRB*PRJ(JL,JORGN)*PY(JL,IORGNTR)&
             &   +PRR(JL,KOHMCOOH)*PY(JL,IOH)*PY(JL,IMCOOH)&
             &   +0.31_JPRB*PRR(JL,KO3C3H6)*PY(JL,IO3)*PY(JL,IC3H6)&
             &   +0.39_JPRB*PRR(JL,KO3TERP)*PY(JL,ITERP)*PY(JL,IO3)&
             &   +2.0*PRJ(JL,JA_ACET)*PY(JL,IACET)+PRJ(JL,JB_ACET)*PY(JL,IACET)&
             &   +PRR(JL,KC50A)*PY(JL,IC2O3)*PY(JL,IHO2) 
            ZXLCH3O2=PRR(JL,KMO2NO)*PY(JL,INO)+(PRR(JL,KMO2HO2A)+PRR(JL,KMO2HO2B))*PY(JL,IHO2)&
             &    +2.0_JPRB*PRR(JL,KMO2MO2)*PY(JL,ICH3O2)+PRR(JL,KNO3MO2)*PY(JL,INO3)&
             &    +PRR(JL,KACO2MO2)*PY(JL,IACO2)
            PY(JL,ICH3O2)=(PY0(JL,ICH3O2)+ZPCH3O2*PDT)/(1._JPRB+ZXLCH3O2*PDT)
            !
            ! -------- ISPD chemistry 
            !
            ZPISPD=0.4*PRR(JL,KC76)*PY(JL,IOH)*PY(JL,IISOP)&
              &  +0.65*PRR(JL,KC77)*PY(JL,IO3)*PY(JL,IISOP)&
              &  +0.2*PRR(JL,KC78)*PY(JL,INO3)*PY(JL,IISOP)
            ZXLISPD=PRR(JL,KOHISPD)*PY(JL,IOH)&
             &  +PRR(JL,KO3ISPD)*PY(JL,IO3)&
             &  +PRR(JL,KNO3ISPD)*PY(JL,INO3)&
             &  +PRJ(JL,JISPD)
            PY(JL,IISPD)=(PY0(JL,IISPD)+ZPISPD*PDT)/(1._JPRB+ZXLISPD*PDT)
            !
            ! -------- ACO2 chemistry 
            !
            ZPACO2=PRR(JL,KOHACET)*PY(JL,IACET)*PY(JL,IOH)
          
            ZXLACO2=PRR(JL,KACO2HO2)*PY(JL,IHO2)+PRR(JL,KACO2MO2)*PY(JL,ICH3O2)&
              &   +PRR(JL,KACO2NO)*PY(JL,INO)+PRR(JL,KACO2XO2)*PY(JL,IXO2)
          
            PY(JL,IACO2)=(PY0(JL,IACO2)+ZPACO2*PDT)/(1._JPRB+ZXLACO2*PDT)
            !
            ! --- CBM4 chem.(short living compounds & operators)
            !
            ZPRXPAR=0.11_JPRB*PRR(JL,KC52)*PY(JL,IOH)*PY(JL,IPAR)&
             &  +2.1_JPRB*PRR(JL,KC53)*PY(JL,IROR)&
             &  +0.7_JPRB*PRR(JL,KC57)*PY(JL,IOLE)*PY(JL,IOH)&
             &  +PRR(JL,KC58)*PY(JL,IO3)*PY(JL,IOLE)&
             &  +PRR(JL,KC59)*PY(JL,IOLE)*PY(JL,INO3)&
             &  +PRR(JL,KOHROOH)*PY(JL,IOH)*PY(JL,IROOH)&
             &  +1.98_JPRB*PRJ(JL,JROOH)*PY(JL,IROOH)&
             &  +1.98_JPRB*PRR(JL,KC84)*PY(JL,IOH)*PY(JL,IORGNTR)&
             &  +1.98_JPRB*PRJ(JL,JORGN)*PY(JL,IORGNTR)

            ZXLRXPAR=PRR(JL,KC83)*PY(JL,IPAR)
            PY(JL,IRXPAR)=(PY0(JL,IRXPAR)+ZPRXPAR*PDT)/(1._JPRB+ZXLRXPAR*PDT)
            ZXLISOP=PRR(JL,KC76)*PY(JL,IOH)+PRR(JL,KC77)*PY(JL,IO3)+PRR(JL,KC78)*PY(JL,INO3)
            PY(JL,IISOP)=PY0(JL,IISOP)/(1._JPRB+ZXLISOP*PDT)
            ZPROR=0.76_JPRB*PRR(JL,KC52)*PY(JL,IPAR)*PY(JL,IOH)+0.02_JPRB*PRR(JL,KC53)*PY(JL,IROR)
            ZXLROR=PRR(JL,KC53)+PRR(JL,KC54)
            PY(JL,IROR)=(PY0(JL,IROR)+ZPROR*PDT)/(1._JPRB+ZXLROR*PDT)

            ZXLTERP=PRR(JL,KOHTERP)*PY(JL,IOH)+PRR(JL,KO3TERP)*PY(JL,IO3)&
             & +PRR(JL,KNO3TERP)*PY(JL,INO3)     
            PY(JL,ITERP)=PY0(JL,ITERP)/(1.+ZXLTERP*PDT)

            ZPXO2=PRR(JL,KC73)*PY(JL,IMGLY)*PY(JL,IOH)+0.87_JPRB*PRR(JL,KC52)*PY(JL,IPAR)*PY(JL,IOH)&
             &    +0.96_JPRB*PRR(JL,KC53)*PY(JL,IROR)+0.8_JPRB*PRR(JL,KC57)*PY(JL,IOLE)*PY(JL,IOH)&
             &    +0.22_JPRB*PRR(JL,KC58)*PY(JL,IO3)*PY(JL,IOLE)+0.91*PRR(JL,KC59)&
             &    *PY(JL,IOLE)*PY(JL,INO3)+PRR(JL,KC61)*PY(JL,IETH)*PY(JL,IOH)&
             &    +0.7_JPRB*PRR(JL,KC76)*PY(JL,IISOP)*PY(JL,IOH)&
             &    +PRR(JL,KC78)*PY(JL,IISOP)*PY(JL,INO3)&
             &    +0.77_JPRB*PRR(JL,KOHROOH)*PY(JL,IROOH)*PY(JL,IOH)&
             &    +0.5_JPRB*PRJ(JL,JROOH)*PY(JL,IROOH)&
             &    +0.51_JPRB*PRR(JL,KC84)*PY(JL,IOH)*PY(JL,IORGNTR)&
             &    +0.51_JPRB*PRJ(JL,JORGN)*PY(JL,IORGNTR)&
             &    +0.2_JPRB*PRR(JL,KC77)*PY(JL,IISOP)*PY(JL,IO3)&
             &    +0.1_JPRB*PRR(JL,KOHETHOH)*PY(JL,IOH)*PY(JL,IETHOH)&
             &    +0.991_JPRB*PRR(JL,KOHC2H6)*PY(JL,IOH)*PY(JL,IC2H6)&
             !VH &    +PRR(JL,kohc3h8)*PY(JL,ic3h8)*PY(JL,ioh)&
             !VH &    +PRR(JL,kohc3h6)*PY(JL,ic3h6)*PY(JL,ioh)&
             &    +1.25_JPRB*PRR(JL,KOHTERP)*PY(JL,IOH)*PY(JL,ITERP)&
             &    +0.76_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
             &    +1.03_JPRB*PRR(JL,KNO3TERP)*PY(JL,INO3)*PY(JL,ITERP)&
             &    +0.713_JPRB*PRR(JL,KOHISPD)*PY(JL,IOH)*PY(JL,IISPD)&
             &    +0.064_JPRB*PRR(JL,KO3ISPD)*PY(JL,IO3)*PY(JL,IISPD)&
             &    +0.075_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)&
             &    +0.7_JPRB*PRJ(JL,JISPD)           
            ZXLXO2=PRR(JL,KC79)*PY(JL,INO)+2.*PRR(JL,KC80)*PY(JL,IXO2)&
             &    +PRR(JL,KC82)*PY(JL,IHO2)+PRR(JL,KNO3XO2)*PY(JL,INO3)&
             &    +PRR(JL,KACO2XO2)*PY(JL,IACO2)+PRR(JL,KXO2XO2N)*PY(JL,IXO2N)

            PY(JL,IXO2)=(PY0(JL,IXO2)+ZPXO2*PDT)/(1.+ZXLXO2*PDT)

            ZPXO2N=0.13_JPRB*PRR(JL,KC52)*PY(JL,IPAR)*PY(JL,IOH)&
             &    +0.04*PRR(JL,KC53)*PY(JL,IROR)&
             &    +0.09_JPRB*PRR(JL,KC59)*PY(JL,IOLE)*PY(JL,INO3)&
             &    +0.009_JPRB*PRR(JL,KOHC2H6)*PY(JL,IOH)*PY(JL,IC2H6)&
             &    +0.088_JPRB*PRR(JL,KC76)*PY(JL,IISOP)*PY(JL,IOH)&
             &    +0.25_JPRB*PRR(JL,KOHTERP)*PY(JL,IOH)*PY(JL,ITERP)&
             &    +0.18_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
             &    +0.25_JPRB*PRR(JL,KNO3TERP)*PY(JL,INO3)*PY(JL,ITERP)
            ZXLXO2N=PRR(JL,KC81)*PY(JL,INO)+PRR(JL,KC85)*PY(JL,IHO2)&
             &    +PRR(JL,KXO2XO2N)*PY(JL,IXO2N)+PRR(JL,KXO2N)*PY(JL,IXO2N)

            PY(JL,IXO2N)=(PY0(JL,IXO2N)+ZPXO2N*PDT)/(1._JPRB+ZXLXO2N*PDT)

            ZPIC3H7O2 =PRR(JL,KOHC3H8)*PY(JL,IOH)*PY(JL,IC3H8)
            ZXLIC3H7O2=PRR(JL,KNOIC3H7O2)*PY(JL,INO)+PRR(JL,KHO2IC3H7O2)*PY(JL,IHO2)
            PY(JL,IIC3H7O2)=(PY0(JL,IIC3H7O2)+ZPIC3H7O2*PDT)/(1._JPRB+ZXLIC3H7O2*PDT)

            ZPHYPROPO2 =PRR(JL,KOHC3H6)*PY(JL,IOH)*PY(JL,IC3H6)
            ZXLHYPROPO2=PRR(JL,KNOHYPROPO2)*PY(JL,INO)+PRR(JL,KHO2HYPROPO2)*PY(JL,IHO2)
            PY(JL,IHYPROPO2)=(PY0(JL,IHYPROPO2)+ZPHYPROPO2*PDT)/(1._JPRB+ZXLHYPROPO2*PDT)

         ENDDO !JL

         IF ( MOD(ITER,2) == 0 ) THEN

            DO JL=KIDIA,KFDIA
               ! --- Species with intermediate lifetimes
               ! --- Inorganic compounds (HNO3 H2O2)
               ! 
               ! VH: TM5-precalculate as first order N2O5 loss, this is now done directly in the ebi-solver
               ! PRR(i,j,kn2o5aq) = PRR(i,j,kn2o5aq)*(PY0(JL,iso4) + PY0(JL,imsa)  )   
               ZPHNO3=PRR(JL,KNO2OH)*PY(JL,INO2)*PY(JL,IOH)&
                &    +2._JPRB*(&
                &    +PRR(JL,KN2O5L)&
                &    +PRR(JL,KN2O5_AER)&
                &    )*PY(JL,IN2O5)&
                &    +PRR(JL,KNO3_AER)*PY(JL,INO3)&
                &    +PRR(JL,KC44)*PY(JL,IALD2)*PY(JL,INO3)&
                &    +PRR(JL,KC41)*PY(JL,ICH2O)*PY(JL,INO3)&
                &    +PRR(JL,KC84)*PY(JL,IOH)*PY(JL,IORGNTR)&
                &    +PRR(JL,KNO3HO2)*PY(JL,INO3)*PY(JL,IHO2)&
                &    +0.15_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)
               ZXLHNO3=PRJ(JL,JHNO3)+PRR(JL,KOHHNO3)*PY(JL,IOH)
               ZXLHNO3=ZXLHNO3 + PVD(JL,IHNO3) 
               PY(JL,IHNO3)=(PY0(JL,IHNO3)+ZPHNO3*PDT)/(1._JPRB+ZXLHNO3*PDT)
               ZPH2O2=PRR(JL,KHO2HO2)*PY(JL,IHO2)*PY(JL,IHO2)&
               ! VH - accoring to Emmons et al., GMD 2011, include H2O2, but
               ! VH - according to Mao et al, ACP 2013 don't include H2O2 formation
               ! VH - for now include it again 
                   &  +0.5_JPRB*PRR(JL,KHO2_AER)*PY(JL,IHO2)&
               !VH but do include contribution from cloud
                   &  +0.5_JPRB*PRR(JL,KHO2L)*PY(JL,IHO2)
               ZXLH2O2=PRJ(JL,JH2O2)+PRR(JL,KHPOH)*PY(JL,IOH)
               ZXLH2O2=ZXLH2O2 + PVD(JL,IH2O2)
               PY(JL,IH2O2)=(PY0(JL,IH2O2)+ZPH2O2*PDT)/(1._JPRB+ZXLH2O2*PDT)
               ! --- CH4-chemistry (methyl peroxide formaldehyde)
               ZPCH3O2H=(PRR(JL,KMO2HO2A)+PRR(JL,KMO2HO2B))*PY(JL,ICH3O2)*PY(JL,IHO2)
               ZXLCH3O2H=PRR(JL,KOHMPER)*PY(JL,IOH)+PRJ(JL,JMEPE)
               ZXLCH3O2H=ZXLCH3O2H + PVD(JL,ICH3O2H)
               PY(JL,ICH3O2H)=(PY0(JL,ICH3O2H)+ZPCH3O2H*PDT)/(1._JPRB+ZXLCH3O2H*PDT)

               ZPCH2O=0.4*PRR(JL,KOHMPER)*PY(JL,ICH3O2H)*PY(JL,IOH)&
                &  +PRR(JL,KMO2NO)*PY(JL,ICH3O2)*PY(JL,INO)&
                &  +1.37_JPRB*PRR(JL,KMO2MO2)*PY(JL,ICH3O2)*PY(JL,ICH3O2)&
                &  +PRJ(JL,JMEPE)*PY(JL,ICH3O2H)&
                &  +0.8_JPRB*PRR(JL,KC57)*PY(JL,IOLE)*PY(JL,IOH)&
                &  +0.74_JPRB*PRR(JL,KC58)*PY(JL,IOLE)*PY(JL,IO3)&
                &  +PRR(JL,KC59)*PY(JL,IOLE)*PY(JL,INO3)&
                &  +1.56_JPRB*PRR(JL,KC61)*PY(JL,IETH)*PY(JL,IOH)&
                &  +PRR(JL,KC62)*PY(JL,IETH)*PY(JL,IO3)&
                &  +0.629_JPRB*PRR(JL,KC76)*PY(JL,IISOP)*PY(JL,IOH)&
                &  +0.6_JPRB*PRR(JL,KC77)*PY(JL,IISOP)*PY(JL,IO3)&
                &  +0.03_JPRB*PRR(JL,KC78)*PY(JL,IISOP)*PY(JL,INO3)&
                &  +PRR(JL,KOHCH3OH)*PY(JL,IOH)*PY(JL,ICH3OH)&
                &  +PRR(JL,KNO3MO2)*PY(JL,INO3)*PY(JL,ICH3O2)&
                &  +0.1_JPRB*PRR(JL,KOHETHOH)*PY(JL,IOH)*PY(JL,IETHOH)&
                &  +0.54_JPRB*PRR(JL,KO3C3H6)*PY(JL,IO3)*PY(JL,IC3H6)&
                &  +1.22_JPRB*PRR(JL,KOHTERP)*PY(JL,IOH)*PY(JL,ITERP)&
                &  +1.8_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
                &  +0.167_JPRB*PRR(JL,KOHISPD)*PY(JL,IOH)*PY(JL,IISPD)&
                &  +0.15_JPRB*PRR(JL,KO3ISPD)*PY(JL,IO3)*PY(JL,IISPD)&
                &  +0.282_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)&
                &  +0.9_JPRB*PRJ(JL,JISPD)*PY(JL,IISPD)&
                &  +PRR(JL,KACO2NO)*PY(JL,IACO2)*PY(JL,INO)&
                &  +PRR(JL,KNOHYPROPO2)*PY(JL,IHYPROPO2)*PY(JL,INO)


               ZXLCH2O=PRJ(JL,JACH2O)+PRJ(JL,JBCH2O)+PY(JL,IOH)*PRR(JL,KFRMOH)&
                &  +PRR(JL,KC41)*PY(JL,INO3)
               ZXLCH2O=ZXLCH2O + PVD(JL,ICH2O)
               PY(JL,ICH2O)=(PY0(JL,ICH2O)+ZPCH2O*PDT)/(1._JPRB+ZXLCH2O*PDT)
               !
               ! --- CBIV-elements for higher HC-chemistry: ALD2 MGLY
               ! --- ETH OLE ISOP ROOH ORGNTR
               !
               ZPALD2=0.11*PRR(JL,KC52)*PY(JL,IPAR)*PY(JL,IOH)&
                & +1.1_JPRB*PRR(JL,KC53)*PY(JL,IROR)&
                & +0.95_JPRB*PRR(JL,KC57)*PY(JL,IOLE)*PY(JL,IOH)&
                & +0.5_JPRB*PRR(JL,KC58)*PY(JL,IOLE)*PY(JL,IO3)&
                & +0.91_JPRB*PRR(JL,KC59)*PY(JL,IOLE)*PY(JL,INO3)&
                & +0.22_JPRB*PRR(JL,KC61)*PY(JL,IETH)*PY(JL,IOH)&
                & +0.8_JPRB*PRR(JL,KC78)*PY(JL,IISOP)*PY(JL,INO3)&
                & +0.04_JPRB*PRR(JL,KOHROOH)*PY(JL,IOH)*PY(JL,IROOH)&
                & +0.991*PRR(JL,KOHC2H6)*PY(JL,IOH)*PY(JL,IC2H6)&
                & +0.3_JPRB*PRJ(JL,JROOH)*PY(JL,IROOH)&
                & +0.3*PRR(JL,KC84)*PY(JL,IOH)*PY(JL,IORGNTR)&
                & +0.3_JPRB*PRJ(JL,JORGN)*PY(JL,IORGNTR)&
                & +PRR(JL,KOHETHOH)*PY(JL,IOH)*PY(JL,IETHOH)&
                & +0.5_JPRB*PRR(JL,KO3C3H6)*PY(JL,IO3)*PY(JL,IC3H6)&
                & +0.47_JPRB*PRR(JL,KOHTERP)*PY(JL,IOH)*PY(JL,ITERP)&
                & +0.21_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
                & +0.47_JPRB*PRR(JL,KNO3TERP)*PY(JL,INO3)*PY(JL,ITERP)&
                & +0.273_JPRB*PRR(JL,KOHISPD)*PY(JL,IOH)*PY(JL,IISPD)&
                & +0.02_JPRB*PRR(JL,KO3ISPD)*PY(JL,IO3)*PY(JL,IISPD)&
                & +0.357_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)&
                & +0.067_JPRB*PRJ(JL,JISPD)*PY(JL,IISPD)&
                & +0.7*PRR(JL,KACO2MO2)*PY(JL,IACO2)*PY(JL,ICH3O2)&
                & +0.27_JPRB*PRR(JL,KNOIC3H7O2)*PY(JL,IIC3H7O2)*PY(JL,INO)&
                & +PRR(JL,KNOHYPROPO2)*PY(JL,IHYPROPO2)*PY(JL,INO)

               ZXLALD2=PRR(JL,KC43)*PY(JL,IOH)+PRR(JL,KC44)*PY(JL,INO3)+PRJ(JL,J45)
               ZXLALD2=ZXLALD2 + PVD(JL,IALD2)
               PY(JL,IALD2)=(PY0(JL,IALD2)+ZPALD2*PDT)/(1._JPRB+ZXLALD2*PDT)
               ZPMGLY=0.19_JPRB*PRR(JL,KOHROOH)*PY(JL,IOH)*PY(JL,IROOH)&
                & +0.168_JPRB*PRR(JL,KOHISPD)*PY(JL,IISPD)*PY(JL,IOH)&
                & +0.85_JPRB*PRR(JL,KO3ISPD)*PY(JL,IISPD)*PY(JL,IO3)&
                & +0.5*PRR(JL,KACO2MO2)*PY(JL,IACO2)*PY(JL,ICH3O2)
               ZXLMGLY=PRR(JL,KC73)*PY(JL,IOH)+PRJ(JL,J74)
               PY(JL,IMGLY)=(PY0(JL,IMGLY)+ZPMGLY*PDT)/(1._JPRB+ZXLMGLY*PDT)
               ZXLETH=PRR(JL,KC61)*PY(JL,IOH)+PRR(JL,KC62)*PY(JL,IO3)
               PY(JL,IETH)=PY0(JL,IETH)/(1._JPRB+ZXLETH*PDT)
               ZPOLE=0._JPRB
               ZXLOLE=PRR(JL,KC57)*PY(JL,IOH)+PRR(JL,KC58)*PY(JL,IO3)+PRR(JL,KC59)*PY(JL,INO3)
               PY(JL,IOLE)=(PY0(JL,IOLE)+ZPOLE*PDT)/(1._JPRB+ZXLOLE*PDT)
               ZPROOH=PRR(JL,KC82)*PY(JL,IXO2)*PY(JL,IHO2)&
                 &   +PRR(JL,KC85)*PY(JL,IHO2)*PY(JL,IXO2N)&
                 &   +PRR(JL,KACO2HO2)*PY(JL,IACO2)*PY(JL,IHO2)&
                 &   +PRR(JL,KHO2IC3H7O2)*PY(JL,IIC3H7O2)*PY(JL,IHO2)&
                 &   +PRR(JL,KHO2HYPROPO2)*PY(JL,IHYPROPO2)*PY(JL,IHO2)
               ZXLROOH=PRJ(JL,JROOH)+PRR(JL,KOHROOH)*PY(JL,IOH)
               ZXLROOH = ZXLROOH + PVD(JL,IROOH)
               PY(JL,IROOH)=(PY0(JL,IROOH)+ZPROOH*PDT)/(1._JPRB+ZXLROOH*PDT)

               ZPORGNTR=PRR(JL,KC81)*PY(JL,INO)*PY(JL,IXO2N)&
                & +0.8_JPRB*PRR(JL,KC78)*PY(JL,IISOP)*PY(JL,INO3)&
                & +PRR(JL,KNO3C3H6)*PY(JL,IC3H6)*PY(JL,INO3)&
                & +0.53_JPRB*PRR(JL,KNO3TERP)*PY(JL,INO3)*PY(JL,ITERP)&
                & +0.85_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)

               ZXLORGNTR=PRR(JL,KC84)*PY(JL,IOH)+PRJ(JL,JORGN)
               ZXLORGNTR=ZXLORGNTR+PVD(JL,IORGNTR)

               PY(JL,IORGNTR)=(PY0(JL,IORGNTR)+ZPORGNTR*PDT)/(1._JPRB+ZXLORGNTR*PDT)
 
               ZXLACET=PRJ(JL,JA_ACET)+PRJ(JL,JB_ACET)+PRR(JL,KOHACET)*PY(JL,IOH)
               !
               ! Extend with source from IC3H7O2...
               ZPACET = 0.82_JPRB*PRR(JL,KNOIC3H7O2)*PY(JL,IIC3H7O2)*PY(JL,INO)
               PY(JL,IACET)=(PY0(JL,IACET)+ZPACET*PDT)/(1.+ZXLACET*PDT)       

               ! PY(JL,IACET)=PY0(JL,IACET)/(1.+ZXLACET*PDT)       


               ! gas phase sulfur   & ammonia

               ZQDMS1=PRR(JL,KDMSOHA)*PY(JL,IOH)+PRR(JL,KDMSNO3)*PY(JL,INO3)
               ZQDMS2=PRR(JL,KDMSOHB)*PY(JL,IOH)
               ZQDMS=ZQDMS1+ZQDMS2
               PY(JL,IDMS)=PY0(JL,IDMS)/(1._JPRB+ZQDMS*PDT)
               ZPSO2=PY(JL,IDMS)*(ZQDMS1+0.75_JPRB*ZQDMS2)
               ZPMSA=PY(JL,IDMS)*0.25_JPRB*ZQDMS2
               ZQSO2=PRR(JL,KSO2OH)*PY(JL,IOH)
               ZQSO2D=ZQSO2 + PVD(JL,ISO2)
               PY(JL,ISO2)=(PY0(JL,ISO2)+ZPSO2*PDT) /(1._JPRB+ZQSO2D*PDT)  !ZQSO2d includes deposition
               PY(JL,IMSA)=(PY0(JL,IMSA)+ZPMSA*PDT) /(1._JPRB+PVD(JL,IMSA)*PDT)  
               PY(JL,ISO4)=(PY0(JL,ISO4)+ZQSO2*PY(JL,ISO2)*PDT) /(1._JPRB + PVD(JL,ISO4)*PDT)   

               ! VH/FD include dry dep for no3_a? (should be copy of so4)
               PY(JL,INO3_A)=PY0(JL,INO3_A) /(1._JPRB+PVD(JL,INO3_A)*PDT) 
               ZY0_IACID=MAX(2._JPRB*PY0(JL,ISO4)+PY0(JL,IMSA)-PY0(JL,INH4),0._JPRB)
               ZY_IACID =(ZY0_IACID+(ZPMSA+2.*ZQSO2*PY(JL,ISO2))*PDT)/&
                  &  (1._JPRB+PRR(JL,KNH3SO4)*PY(JL,INH3)*PDT)
               PY(JL,INH4)=(PY0(JL,INH4))/(1._JPRB+PVD(JL,INH4)*PDT)
               ZPNH2=PY(JL,IOH)*PRR(JL,KNH3OH)
               ZQNH3 = ZPNH2+ PVD(JL,INH3)
               PY(JL,INH3)=PY0(JL,INH3)/(1._JPRB+ZQNH3*PDT)
               ZQNH2= PRR(JL,KNH2NO)*PY(JL,INO)+PRR(JL,KNH2NO2)*PY(JL,INO2)&
                  &  +PRR(JL,KNH2HO2)*PY(JL,IHO2) +PRR(JL,KNH2O2) !VH +PRR(JL,KNH2O3)*PY(JL,IO3)
               PY(JL,INH2)=(PY0(JL,INH2)+PY(JL,INH3)*ZPNH2*PDT)/(1._JPRB+ZQNH2*PDT)
            ENDDO  !JL

         ENDIF

         IF ( MOD(ITER,KMAXIT) == 0 ) THEN

            ! --- Long living compounds
            DO JL=KIDIA,KFDIA
               PY(JL,ICH4)=PY0(JL,ICH4)/(1._JPRB+PRR(JL,KCH4OH)*PY(JL,IOH)*PDT)

               ZPCO=PY(JL,ICH2O)*(PRJ(JL,JACH2O)+PRJ(JL,JBCH2O)&
               & +PY(JL,IOH)*PRR(JL,KFRMOH))&
               & +PRJ(JL,J45)*PY(JL,IALD2)&
               & +PRJ(JL,J74)*PY(JL,IMGLY)&
               & +0.62_JPRB*PRR(JL,KC57)*PY(JL,IOLE)*PY(JL,IOH)&
               & +0.65_JPRB*PRR(JL,KC58)*PY(JL,IOLE)*PY(JL,IO3)&
               & +0.56_JPRB*PRR(JL,KC59)*PY(JL,IOLE)*PY(JL,INO3)&
               & +0.24_JPRB*PRR(JL,KC62)*PY(JL,IETH)*PY(JL,IO3)&
               & +0.066_JPRB*PRR(JL,KC77)*PY(JL,IISOP)*PY(JL,IO3)&
               & +PRR(JL,KC41)*PY(JL,ICH2O)*PY(JL,INO3)&
               & +0.56_JPRB*PRR(JL,KO3C3H6)*PY(JL,IO3)*PY(JL,IC3H6)&
               & +0.47_JPRB*PRR(JL,KOHTERP)*PY(JL,IOH)*PY(JL,ITERP)&
               & +0.211_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
               & +0.47_JPRB*PRR(JL,KNO3TERP)*PY(JL,INO3)*PY(JL,ITERP)&
               & +0.334_JPRB*PRR(JL,KOHISPD)*PY(JL,IOH)*PY(JL,IISPD)&
               & +0.225_JPRB*PRR(JL,KO3ISPD)*PY(JL,IO3)*PY(JL,IISPD)&
               & +0.643*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)&
               & +0.333*PRJ(JL,JISPD)*PY(JL,IISPD)&
               & +PRJ(JL,JA_ACET)*PY(JL,IACET) 
               ZXLCO = PRR(JL,KCOOH)*PY(JL,IOH)
               ZXLCO = ZXLCO + PVD(JL,ICO)
               PY(JL,ICO)=(PY0(JL,ICO)+ZPCO*PDT)/(1._JPRB+ZXLCO*PDT)

               ZPCH3OH = 0.63_JPRB*PRR(JL,KMO2MO2)*PY(JL,ICH3O2)*PY(JL,ICH3O2)&
               & + 0.5*PRR(JL,KACO2MO2)*PY(JL,ICH3O2)*PY(JL,IACO2)
               PY(JL,ICH3OH)=(PY0(JL,ICH3OH)+(ZPCH3OH*PDT))/&
                &     (1._JPRB+(PVD(JL,ICH3OH)+PRR(JL,KOHCH3OH)*PY(JL,IOH))*PDT)
               ZPHCOOH=0.52_JPRB*PRR(JL,KC62)*PY(JL,IO3)*PY(JL,IETH)&
               & + 0.25*PRR(JL,KO3C3H6)*PY(JL,IO3)*PY(JL,IC3H6)      
               PY(JL,IHCOOH)=(PY0(JL,IHCOOH)+(ZPHCOOH*PDT))/&
                &     (1._JPRB+PRR(JL,KOHHCOOH)*PY(JL,IOH)*PDT)
                     
               ZPMCOOH=PRR(JL,KC50B)*PY(JL,IC2O3)*PY(JL,IHO2)
                             
               PY(JL,IMCOOH)=(PY0(JL,IMCOOH)+(ZPMCOOH*PDT))/&
                &     (1._JPRB+(PVD(JL,IMCOOH)+PRR(JL,KOHMCOOH)*PY(JL,IOH))*PDT)
                     
                     
               PY(JL,IC2H6)=(PY0(JL,IC2H6))/&
                &    (1._JPRB+PRR(JL,KOHC2H6)*PY(JL,IOH)*PDT)
           
               PY(JL,IETHOH)=(PY0(JL,IETHOH)/(1._JPRB+(PVD(JL,IETHOH)+PRR(JL,KOHETHOH)*PY(JL,IOH))*PDT))
           
               PY(JL,IC3H8)=(PY0(JL,IC3H8)/(1._JPRB+PRR(JL,KOHC3H8)*PY(JL,IOH)*PDT))
          
               ZXLC3H6=PRR(JL,KOHC3H6)*PY(JL,IOH)&
                &  +PRR(JL,KO3C3H6)*PY(JL,IO3)&
                &  +PRR(JL,KNO3C3H6)*PY(JL,INO3)
               PY(JL,IC3H6)=(PY0(JL,IC3H6)/(1._JPRB+ZXLC3H6*PDT))
                             
               ZPPAR=0.35_JPRB*PRR(JL,KC77)*PY(JL,IISOP)*PY(JL,IO3)&
                & +2.4_JPRB*PRR(JL,KC78)*PY(JL,IISOP)*PY(JL,INO3)&
                & +5.0_JPRB*PRR(JL,KOHTERP)*PY(JL,IOH)*PY(JL,ITERP)&
                & +6.0_JPRB*PRR(JL,KO3TERP)*PY(JL,IO3)*PY(JL,ITERP)&
                & +1.565_JPRB*PRR(JL,KOHISPD)*PY(JL,IOH)*PY(JL,IISPD)&
                & +0.36_JPRB*PRR(JL,KO3ISPD)*PY(JL,IO3)*PY(JL,IISPD)&
                & +1.282_JPRB*PRR(JL,KNO3ISPD)*PY(JL,INO3)*PY(JL,IISPD)&
                & +0.832_JPRB*PRJ(JL,JISPD)*PY(JL,IISPD)&
                & +0.6*PRR(JL,KACO2MO2)*PY(JL,ICH3O2)*PY(JL,IACO2)  
                
               ZXLPAR=PRR(JL,KC52)*PY(JL,IOH)+PRR(JL,KC83)*PY(JL,IRXPAR)
               PY(JL,IPAR)=(PY0(JL,IPAR)+ZPPAR*PDT)/(1._JPRB+ZXLPAR*PDT)

               !
               !cmk ____added rn222 chemistry in EBI language
               !
               PY(JL,IRN222) = PY0(JL,IRN222)/(1._JPRB+PRR(JL,KRN222)*PDT)
               PY(JL,IPB210) = PY0(JL,IPB210)+PY0(JL,IRN222)-PY(JL,IRN222)
               !IF(y(JL,ipb210) < 0.0 .or. PY(JL,irn222) < 0.0) THEN
               !   print *, 'Negatives .....rn222, pb210', PY0(JL,irn222), PY(JL,irn222) , PY0(JL,ipb210),  PY(JL,ipb210)
               !ENDIF
            ENDDO   !JL

         ENDIF

      ENDDO !ITER

IF (LHOOK) CALL DR_HOOK('TM5_DO_EBI',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_DO_EBI
