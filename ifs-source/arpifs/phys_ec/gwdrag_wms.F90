! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GWDRAG_WMS(YDSTA,YDEGWD,YDEGWWMS,KIDIA,  KFDIA,   KLON,  KLEV,   KLAUNCH, PTSTEP,&
                     &PTM1 ,  PUM1,    PVM1,  PAPM1,  PAPHM1, PGEO1 ,&
                     &PGELAT, PGAW,    PPRECIP,&
                     &PTENU,  PTENV,   PFLUXU, PFLUXV)

!**** *GWDRAG_WMS*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME

!     Original Fortran Code by   J. SCINOCCIA
!     Rewritten in IFS format by A. ORR          E.C.M.W.F.     August 2008

!     PURPOSE
!     -------

!          THIS ROUTINE COMPUTES NON-OROGRAPHIC GRAVITY WAVE DRAG
!     AFTER SCINOCCA (2003) AND Mc LANDRESS AND SCINOCCIA (JAS 2005)
!     HYDROSTATIC NON-ROTATIONAL SIMPLIFIED VERSION OF THE 
!     WARNER AND MCINTYRE (1996) NON-OROGRAPHIC GRAVITY WAVE PARAMETERIZATION 
!     CONSTANTS HAVE BEEN OPTIMIZED FOLLOWING M. ERN ET AL. (ATMOS. CHEM. PHYS. 2006)

!     REFERENCE: Orr, A., P. Bechtold, J. Scinoccia, M. Ern, M. Janiskova, 2010: 
!                Improved middle atmosphere climate and analysis in the ECMWF forecasting system 
!                through a non-orographic gravity wave parametrization. J.  Climate., 23, 5905-5926.

!     LAUNCH SPECTRUM - GENERALIZED DESAUBIES 
!     INCLUDES A CRITICAL-LEVEL CORRECTION THAT PREVENTS THE 
!     MOMEMTUM DEPOSITION IN EACH AZIMUTH FROM DRIVING THE FLOW TO SPEEDS FASTER 
!     THAN THE PHASE SPEED OF THE WAVES, I.E. WHEN WAVES BREAK THEY DRAG THE MEAN 
!     FLOW TOWARDS THEIR PHASE SPEED - NOT PAST IT.  

!**   INTERFACE.
!     ----------
!          *GWDRAG_WMS* IS CALLED FROM *CALLPAR*

!     MODIFICATIONS
!     -------------
!           October 2008 : Cleaning and full adaptation   P. Bechtold/JJMorcrette
!                          Optimisation+bug corrections   P. Bechtold
!           November 2012: possibility of different ZGAM transform T. Stockdale
!           November 2015: resolution scaling for octahedral grid  P. Bechtold
! ---------------------------------------------------------------------------------


USE YOMSTA     , ONLY : TSTA
USE PARKIND1,    ONLY : JPIM, JPRB
USE YOEGWWMS,    ONLY : TEGWWMS
USE YOMCST,      ONLY : RG, RD, RCPD, RPI, RA
USE YOMHOOK,     ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOEGWD,      ONLY : TEGWD

IMPLICIT NONE

!in
TYPE(TSTA)        ,INTENT(IN) :: YDSTA
TYPE(TEGWD)       ,INTENT(INOUT):: YDEGWD
TYPE(TEGWWMS)     ,INTENT(INOUT):: YDEGWWMS
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON                 ! horizontal dimension
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV                 ! vertical levels
INTEGER(KIND=JPIM),INTENT(IN) :: KLAUNCH              ! index for launch level 
REAL(KIND=JPRB)   ,INTENT(IN) :: PTSTEP               ! model time step
REAL(KIND=JPRB)   ,INTENT(IN) :: PUM1(KLON,KLEV)      ! full model level zonal velocity (t-dt)
REAL(KIND=JPRB)   ,INTENT(IN) :: PVM1(KLON,KLEV)      ! full model level meridional velocity (t-dt)
REAL(KIND=JPRB)   ,INTENT(IN) :: PTM1(KLON,KLEV)      ! full model level temperature (t-dt)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPM1(KLON,KLEV)     ! full model level pressure (t-dt)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPHM1(KLON,KLEV+1)  ! half-model level pressure (t-dt)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEO1(KLON,KLEV)     ! full model level geopotential
REAL(KIND=JPRB)   ,INTENT(IN) :: PGELAT(KLON)         ! latitude
REAL(KIND=JPRB)   ,INTENT(IN) :: PGAW(KLON)           ! normalised gaussian quadrature weight/number of longitude points
                                                      ! local sub-area == 4*RPI*RA**2 * PGAW
REAL(KIND=JPRB)   ,INTENT(IN) :: PPRECIP(KLON)        ! total surface precipitation

!out
REAL(KIND=JPRB)   ,INTENT(OUT):: PTENU(KLON,KLEV)     ! full-model level zonal momentum tendency
REAL(KIND=JPRB)   ,INTENT(OUT):: PTENV(KLON,KLEV)     ! full-model level meridional momentum tendency
REAL(KIND=JPRB)   ,INTENT(OUT):: PFLUXU(KLON,KLEV+1)  ! = zonal component of vertical momentum flux (Pa)  
REAL(KIND=JPRB)   ,INTENT(OUT):: PFLUXV(KLON,KLEV+1)  ! = meridional component of vertical momentum flux (Pa)
!work
INTEGER(KIND=JPIM), PARAMETER  :: IAZIDIM=4     !number of azimuths
INTEGER(KIND=JPIM), PARAMETER  :: INCDIM=20     !number of discretized c spectral elements in launch spectrum  

REAL(KIND=JPRB) :: ZUHM1(KLON,KLEV)             !half-model level zonal velocity
REAL(KIND=JPRB) :: ZVHM1(KLON,KLEV)             !half-model level meridional velocity
REAL(KIND=JPRB) :: ZBVFHM1(KLON,KLEV)           !half-model level Brunt-Vaisalla frequency
REAL(KIND=JPRB) :: ZRHOHM1(KLON,KLEV)           !half-model level density
REAL(KIND=JPRB) :: ZX(INCDIM)                   !coordinate transformation
REAL(KIND=JPRB) :: ZCI(INCDIM)                  !phase speed element
REAL(KIND=JPRB) :: ZDCI(INCDIM)
REAL(KIND=JPRB) :: ZUI(KLON,KLEV,IAZIDIM)       !intrinsic velocity
REAL(KIND=JPRB) :: ZUL(KLON,IAZIDIM)            !velocity in azimuthal direction at launch level
REAL(KIND=JPRB) :: ZBVFL(KLON)                  !buoyancy at launch level
REAL(KIND=JPRB) :: ZCOSANG(IAZIDIM)             !cos of azimuth angle
REAL(KIND=JPRB) :: ZSINANG(IAZIDIM)             !sin of azimuth angle
REAL(KIND=JPRB) :: ZFCT(KLON,KLEV)
REAL(KIND=JPRB) :: ZFNORM(KLON)                 !normalisation factor (A)
REAL(KIND=JPRB) :: ZCI_MIN(KLON,IAZIDIM)
REAL(KIND=JPRB) :: ZTHM1(KLON,KLEV)             !temperature on half-model levels
REAL(KIND=JPRB) :: ZFLUX(KLON,INCDIM,IAZIDIM)   !momentum flux at each vertical level and azimuth 
REAL(KIND=JPRB) :: ZPU(KLON,KLEV,IAZIDIM)       !momentum flux
REAL(KIND=JPRB) :: ZDFL(KLON,KLEV,IAZIDIM)
REAL(KIND=JPRB) :: ZACT(KLON,INCDIM,IAZIDIM)    !if =1 then critical level encountered  
REAL(KIND=JPRB) :: ZACC(KLON,INCDIM,IAZIDIM)
REAL(KIND=JPRB) :: ZCRT(KLON,KLEV,IAZIDIM)

INTEGER(KIND=JPIM) :: ILAUNCH                   !model level from which GW spectrum is launched
INTEGER(KIND=JPIM) :: INC, JK, JL, IAZI
REAL(KIND=JPRB) :: ZRADTODEG, ZGELATDEG
REAL(KIND=JPRB) :: ZCIMIN, ZCIMAX
REAL(KIND=JPRB) :: ZGAM, ZPEXP, ZXMAX, ZXMIN, ZXRAN, ZDX, ZX1, ZX2, ZDXA, ZDXB, ZDXS
REAL(KIND=JPRB) :: ZANG, ZAZ_FCT, ZNORM, ZANG1, ZTX
REAL(KIND=JPRB) :: ZU, ZCIN, ZCPEAK
REAL(KIND=JPRB) :: ZCIN4, ZBVFL4, ZCIN2, ZBVFL2, ZCIN3, ZBVFL3, ZCINC
REAL(KIND=JPRB) :: ZATMP, ZFLUXS, ZDEP, ZFLUXSQ, ZULM, ZDFT, ZE1, ZE2
REAL(KIND=JPRB) :: ZMS_L,ZMS, Z0P5, Z0P0, Z50S
REAL(KIND=JPRB) :: ZGAUSS(KLON), ZFLUXLAUN(KLON), ZCNGL(KLON)
REAL(KIND=JPRB) :: ZCONS1,ZCONS2,ZDELP,ZRGPTS
REAL(KIND=JPRB) :: ZTHSTD,ZRHOSTD,ZBVFSTD
REAL(KIND=JPRB) :: ZGAUSSB,ZFLUXGLOB

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GWDRAG_WMS',0,ZHOOK_HANDLE)
ASSOCIATE(GSSEC=>YDEGWD%GSSEC, &
 & GCOEFF=>YDEGWWMS%GCOEFF, GCSTAR=>YDEGWWMS%GCSTAR, &
 & GFLUXLAUNL=>YDEGWWMS%GFLUXLAUNL, GGAUSSA=>YDEGWWMS%GGAUSSA, &
 & GGAUSSB=>YDEGWWMS%GGAUSSB, GMSTAR_L=>YDEGWWMS%GMSTAR_L, &
 & GPTWO=>YDEGWWMS%GPTWO, LGACALC=>YDEGWWMS%LGACALC, LGINDL=>YDEGWWMS%LGINDL, &
 & LGSATL=>YDEGWWMS%LGSATL, LOZPR=>YDEGWWMS%LOZPR, NGAUSS=>YDEGWWMS%NGAUSS, &
 & NLAUNCHL=>YDEGWWMS%NLAUNCHL, NSLOPE=>YDEGWWMS%NSLOPE, &
 & STPHI=>YDSTA%STPHI, STPREH=>YDSTA%STPREH, STTEM=>YDSTA%STTEM)
!--------------------------------------------------------------------------

!*       INPUT PARAMETERS
!*       ---------------- 

ZRADTODEG=57.29577951_JPRB

! Set parameters which are a function of launch height
ILAUNCH=NLAUNCHL(KLAUNCH)
ZFLUXGLOB=GFLUXLAUNL(KLAUNCH)
ZGAUSSB=GGAUSSB(KLAUNCH)

ZMS_L=GMSTAR_L(KLAUNCH)
ZMS=2*RPI/ZMS_L    

!*       INITIALIZE FIELDS TO ZERO
!*       ------------------------- 

PTENU(:,:)=0.0_JPRB
PTENV(:,:)=0.0_JPRB

DO IAZI=1,IAZIDIM
   DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
         ZPU(JL,JK,IAZI)=0.0_JPRB
         ZCRT(JL,JK,IAZI)=0.0_JPRB
         ZDFL(JL,JK,IAZI)=0.0_JPRB
      ENDDO
   ENDDO      
ENDDO         

DO JK=1,KLEV+1
   DO JL=KIDIA,KFDIA
     PFLUXU(JL,JK)=0.0_JPRB
     PFLUXV(JL,JK)=0.0_JPRB
   ENDDO
ENDDO


!*       INITIALIZE PARAMETERS FOR COORDINATE TRANSFORM
!*       ----------------------------------------------

! ZCIMIN,ZCIMAX - min,max intrinsic launch-level phase speed (c-U_o) (m/s)
! ZGAM - half=width of coordinate stretch

ZCIMIN=0.50_JPRB
ZCIMAX=100.0_JPRB
ZGAM=0.25_JPRB

ZPEXP=GPTWO/2.0_JPRB

! set initial min ci in each column and azimuth (used for critical levels)

DO IAZI=1,IAZIDIM
   DO JL=KIDIA,KFDIA
      ZCI_MIN(JL,IAZI)=ZCIMIN
   ENDDO
ENDDO

   
!*       DEFINE HALF MODEL LEVEL WINDS AND TEMPERATURE
!*       -----------------------------------

DO JK=2,KLEV
  DO JL=KIDIA,KFDIA                                        
     ZTHM1(JL,JK) =0.5_JPRB*(PTM1(JL,JK-1)+PTM1(JL,JK)) 
     ZUHM1(JL,JK) =0.5_JPRB*(PUM1(JL,JK-1)+PUM1(JL,JK))
     ZVHM1(JL,JK) =0.5_JPRB*(PVM1(JL,JK-1)+PVM1(JL,JK))
  ENDDO
ENDDO
JK=1
DO JL=KIDIA,KFDIA                                        
   ZTHM1(JL,JK)=PTM1(JL,JK) 
   ZUHM1(JL,JK)=PUM1(JL,JK)
   ZVHM1(JL,JK)=PVM1(JL,JK)
ENDDO


!*       DEFINE STATIC STABILITY AND AIR DENSITY ON HALF MODEL LEVELS
!*       ------------------------------------------------------------

ZCONS1=1.0_JPRB/RD
ZCONS2=RG**2/RCPD
DO JK=KLEV,2,-1
  DO JL=KIDIA,KFDIA
!    ZDELP=PAPM1(JL,JK)-PAPM1(JL,JK-1)
    ZDELP=PGEO1(JL,JK)-PGEO1(JL,JK-1)
    ZRHOHM1(JL,JK)=PAPHM1(JL,JK)*ZCONS1/ZTHM1(JL,JK)
    ZBVFHM1(JL,JK)=ZCONS2/ZTHM1(JL,JK)*&
!     & (1.0_JPRB-RCPD*ZRHOHM1(JL,JK)*(PTM1(JL,JK)-PTM1(JL,JK-1))/ZDELP)  
     & (1.0_JPRB+RCPD*(PTM1(JL,JK)-PTM1(JL,JK-1))/ZDELP)
    ZBVFHM1(JL,JK)=MAX(ZBVFHM1(JL,JK),GSSEC)
    ZBVFHM1(JL,JK)=SQRT(ZBVFHM1(JL,JK))
  ENDDO
ENDDO

!*       SET UP AZIMUTH DIRECTIONS AND SOME TRIG FACTORS
!*       -----------------------------------------------     

ZANG=2*RPI/IAZIDIM
ZAZ_FCT=1.0_JPRB

! get normalization factor to ensure that the same amount of momentum
! flux is directed (n,s,e,w) no mater how many azimuths are selected.
! note, however, the code below assumes a symmetric distribution of
! of azimuthal directions (ie 4,8,16,32,...)

ZNORM=0.0_JPRB
DO IAZI=1,IAZIDIM
  ZANG1=(IAZI-1)*ZANG
  ZCOSANG(IAZI)=COS(ZANG1)
  ZSINANG(IAZI)=SIN(ZANG1)
  ZNORM=ZNORM+ABS(ZCOSANG(IAZI))
ENDDO
ZAZ_FCT=2._JPRB*ZAZ_FCT/ZNORM


!*       DEFINE COORDINATE TRANSFORM
!*       -----------------------------------------------     

! note that this is expresed in terms of the intrinsic phase speed
! at launch ci=c-u_o so that the transformation is identical at every
! launch site.
! See Eq. 28-30 of Scinocca 2003.
     
ZXMAX=1.0_JPRB/ZCIMIN
ZXMIN=1.0_JPRB/ZCIMAX

ZXRAN=ZXMAX-ZXMIN
ZDX=ZXRAN/REAL(INCDIM-1)
IF(LGACALC) ZGAM=(ZXMAX-ZXMIN)/LOG(ZXMAX/ZXMIN)
ZX1=ZXRAN/(EXP(ZXRAN/ZGAM)-1.0_JPRB)
ZX2=ZXMIN-ZX1

DO INC=1,INCDIM
   ZTX=REAL(INC-1)*ZDX+ZXMIN
   ZX(INC)=ZX1*EXP((ZTX-ZXMIN)/ZGAM)+ZX2                       !Eq. 29 of Scinocca 2003
   ZCI(INC)=1.0_JPRB/ZX(INC)                                   !Eq. 28 of Scinocca 2003
   ZDCI(INC)=ZCI(INC)**2*(ZX1/ZGAM)*EXP((ZTX-ZXMIN)/ZGAM)*ZDX  !Eq. 30 of Scinocca 2003
ENDDO


!*       DEFINE INTRINSIC VELOCITY (RELATIVE TO LAUNCH LEVEL VELOCITY) U(Z)-U(Zo), AND COEFFICINETS
!*       ------------------------------------------------------------------------------------------        

DO IAZI=1,IAZIDIM
   DO JL=KIDIA,KFDIA
      ZUL(JL,IAZI)=ZCOSANG(IAZI)*ZUHM1(JL,ILAUNCH)+ZSINANG(IAZI)*ZVHM1(JL,ILAUNCH)
   ENDDO
ENDDO
DO JL=KIDIA,KFDIA
   ZBVFL(JL)=ZBVFHM1(JL,ILAUNCH)
ENDDO

DO JK=2,ILAUNCH
   DO IAZI=1,IAZIDIM
      DO JL=KIDIA,KFDIA
         ZU=ZCOSANG(IAZI)*ZUHM1(JL,JK)+ZSINANG(IAZI)*ZVHM1(JL,JK)
         ZUI(JL,JK,IAZI)=ZU-ZUL(JL,IAZI)
      ENDDO
   ENDDO
ENDDO

!*       DEFINE RHO(Zo)/N(Zo)
!*       ------------------- 
DO JK=2,ILAUNCH
   DO JL=KIDIA,KFDIA
      ZFCT(JL,JK)=ZRHOHM1(JL,JK)/ZBVFHM1(JL,JK)
   ENDDO
ENDDO

! Optionally set ZFCT at launch level using standard atmos values, to ensure saturation is
! independent of location
IF (LGINDL) THEN
  ZCONS1=1.0_JPRB/RD
  ZCONS2=RG**2/RCPD
  ZDELP=STPHI(ILAUNCH)-STPHI(ILAUNCH-1)
  ZTHSTD=0.5_JPRB*(STTEM(ILAUNCH-1)+STTEM(ILAUNCH)) 
  ZRHOSTD=STPREH(ILAUNCH-1)*ZCONS1/ZTHSTD
  ZBVFSTD=ZCONS2/ZTHSTD*(1.0_JPRB+RCPD*(STTEM(ILAUNCH)-STTEM(ILAUNCH-1))/ZDELP)
  ZBVFSTD=MAX(ZBVFSTD,GSSEC)
  ZBVFSTD=SQRT(ZBVFSTD)
  DO JL=KIDIA,KFDIA
     ZFCT(JL,ILAUNCH)=ZRHOSTD/ZBVFSTD
  ENDDO
ENDIF

!*       SET LAUNCH MOMENTUM FLUX SPECTRAL DENSITY
!*       ----------------------------------------- 

! Eq. (25) of Scinocca 2003 (not including the 'A' component), and with U-Uo=0
! do this for only one azimuth since it is identical to all azimuths, and it will be renormalized
! Initial spectrum fully saturated if LGSATL

IF(NSLOPE==1) THEN
! s=1 case
   DO INC=1,INCDIM
      ZCIN=ZCI(INC)
      ZCIN4=(ZMS*ZCIN)**4
      DO JL=KIDIA,KFDIA
         ZBVFL4=ZBVFL(JL)**4
         IF(LGSATL) THEN
           ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*MIN(ZCIN/ZCIN4,ZCIN/ZBVFL4)
         ELSE
           ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*ZCIN/(ZBVFL4+ZCIN4)
         ENDIF
         ZACT(JL,INC,1)=1.0_JPRB
      ENDDO
   ENDDO
ELSEIF(NSLOPE==2) THEN
! s=2 case
   DO INC=1,INCDIM
      ZCIN=ZCI(INC)
      ZCIN4=(ZMS*ZCIN)**4
      DO JL=KIDIA,KFDIA
         ZBVFL4=ZBVFL(JL)**4
         ZCPEAK=ZBVFL(JL)/ZMS
         IF(LGSATL) THEN
           ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*MIN(ZCPEAK/ZCIN4,ZCIN/ZBVFL4)
         ELSE
           ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*ZCIN*ZCPEAK/(ZBVFL4*ZCPEAK+ZCIN4*ZCIN)
         ENDIF
         ZACT(JL,INC,1)=1.0_JPRB
      ENDDO
   ENDDO
ELSEIF(NSLOPE==-1) THEN
! s=-1 case
  DO INC=1,INCDIM
     ZCIN=ZCI(INC)
     ZCIN2=(ZMS*ZCIN)**2
     DO JL=KIDIA,KFDIA
        ZBVFL2=ZBVFL(JL)**2
        ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL2*ZCIN/(ZBVFL2+ZCIN2)
        ZACT(JL,INC,1)=1.0_JPRB
     ENDDO
  ENDDO
ELSEIF(NSLOPE==0) THEN
! s=0 case
  DO INC=1,INCDIM
     ZCIN=ZCI(INC)
     ZCIN3=(ZMS*ZCIN)**3
     DO JL=KIDIA,KFDIA
        ZBVFL3=ZBVFL(JL)**3
        ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL3*ZCIN/(ZBVFL3+ZCIN3)
        ZACT(JL,INC,1)=1.0_JPRB
        ZACC(JL,INC,1)=1.0_JPRB
     ENDDO
  ENDDO
ENDIF

!*       NORMALIZE LAUNCH MOMENTUM FLUX
!*       ------------------------------

! (rho x F^H = rho_o x F_p^total)

! integrate (ZFLUX x dX)
DO INC=1,INCDIM
   ZCINC=ZDCI(INC)
   DO JL=KIDIA,KFDIA
      ZPU(JL,ILAUNCH,1)=ZPU(JL,ILAUNCH,1)+ZFLUX(JL,INC,1)*ZCINC
   ENDDO
ENDDO

!*       NORMALIZE GFLUXLAUN TO INCLUDE SENSITIVITY TO PRECIPITATION
!*       -----------------------------------------------------------

! Also other options to alter tropical values

! A=ZFNORM in Scinocca 2003.  A is independent of height.
ZDXA=1.0_JPRB/29.E3_JPRB
ZDXB=1.0_JPRB/3.5E3_JPRB
DO JL=KIDIA,KFDIA       
  ZDX=MAX(1.E2_JPRB,2*RA*SQRT(RPI*PGAW(JL))) !grid resolution (m)
   !Scaling factor for launch flux depending on grid resolution
   ! smooth reduction below 30 km
  ZDXS=1.0_JPRB-MIN(1.0_JPRB,ATAN((MAX(1.0_JPRB/ZDX,ZDXA)-ZDXA)/(ZDXB-ZDXA)))
  ZFLUXLAUN(JL)=ZFLUXGLOB*ZDXS
  ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)
ENDDO

! If LOZPR=TRUR then vary EPLAUNCH over tropics      
IF (LOZPR) THEN
  IF (NGAUSS==1) THEN
    DO JL=KIDIA,KFDIA       
       ZFLUXLAUN(JL)=ZFLUXLAUN(JL)*(1.0_JPRB+MIN(0.5_JPRB,GCOEFF*PPRECIP(JL)))     !precip
       ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)
    ENDDO
  ELSEIF (NGAUSS==2) THEN           
    DO JL=KIDIA,KFDIA       
      ZGELATDEG=PGELAT(JL)*ZRADTODEG 
      ZGAUSS(JL)=ZGAUSSB*EXP((-ZGELATDEG*ZGELATDEG)/(2*GGAUSSA*GGAUSSA))
      ZFLUXLAUN(JL)=(1.0_JPRB+ZGAUSS(JL))*ZFLUXLAUN(JL)
      ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)
    ENDDO
  ELSEIF (NGAUSS==4) THEN           
! Set latitudinal dependence to optimize stratospheric winds for 36r1
    Z50S=-50.0_JPRB
    DO JL=KIDIA,KFDIA       
      ZGELATDEG=PGELAT(JL)*ZRADTODEG-Z50S
      ZGAUSS(JL)=ZGAUSSB*EXP((-ZGELATDEG*ZGELATDEG)/(2*GGAUSSA*GGAUSSA))
      ZFLUXLAUN(JL)=(1.0_JPRB+ZGAUSS(JL))*ZFLUXLAUN(JL)
      ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)
    ENDDO
  ENDIF
ENDIF 
 
DO IAZI=1,IAZIDIM
   DO JL=KIDIA,KFDIA
      ZPU(JL,ILAUNCH,IAZI)=ZFLUXLAUN(JL)
   ENDDO
ENDDO

!*       ADJUST CONSTANT ZFCT
!*       --------------------
DO JK=2,ILAUNCH
   DO JL=KIDIA,KFDIA
      ZFCT(JL,JK)=ZFNORM(JL)*ZFCT(JL,JK)
   ENDDO
ENDDO

!*       RENORMALIZE EACH SPECTRAL ELEMENT IN FIRST AZIMUTH
!*       --------------------------------------------------
DO INC=1,INCDIM
   DO JL=KIDIA,KFDIA
      ZFLUX(JL,INC,1)=ZFNORM(JL)*ZFLUX(JL,INC,1)
   ENDDO
ENDDO

!*       COPY ZFLUX INTO ALL OTHER AZIMUTHS
!*       --------------------------------

! ZACT=1 then no critical level
! ZACT=0 then critical level

DO IAZI=2,IAZIDIM
   DO INC=1,INCDIM
      DO JL=KIDIA,KFDIA
         ZFLUX(JL,INC,IAZI)=ZFLUX(JL,INC,1)
         ZACT(JL,INC,IAZI)=1.0_JPRB
         ZACC(JL,INC,IAZI)=1.0_JPRB
      ENDDO
   ENDDO
ENDDO

! -----------------------------------------------------------------------------

!*       BEGIN MAIN LOOP OVER LEVELS
!*       ---------------------------

!* begin IAZIDIM do-loop
!* --------------------

DO IAZI=1,IAZIDIM

!* begin JK do-loop
!* ----------------

   DO JK=ILAUNCH-1,2,-1


!* first do critical levels
!* ------------------------

         DO JL=KIDIA,KFDIA
            ZCI_MIN(JL,IAZI)=MAX(ZCI_MIN(JL,IAZI),ZUI(JL,JK,IAZI))               
         ENDDO

!* set ZACT to zero if critical level encountered
!* ----------------------------------------------

         Z0P5=0.5_JPRB
         DO INC=1,INCDIM
            ZCIN=ZCI(INC)
            DO JL=KIDIA,KFDIA
               ZATMP=Z0P5+SIGN(Z0P5,ZCIN-ZCI_MIN(JL,IAZI))
               ZACC(JL,INC,IAZI)=ZACT(JL,INC,IAZI)-ZATMP
               ZACT(JL,INC,IAZI)=ZATMP
            ENDDO
         ENDDO

!* integrate to get critical-level contribution to mom deposition on this level, i.e. ZACC=1
!* ----------------------------------------------------------------------------------------

         DO INC=1,INCDIM
            ZCINC=ZDCI(INC)
            DO JL=KIDIA,KFDIA
               ZDFL(JL,JK,IAZI)=ZDFL(JL,JK,IAZI)+&
     &                ZACC(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZCINC
            ENDDO
         ENDDO

!* get weighted average of phase speed in layer
!* --------------------------------------------

         DO JL=KIDIA,KFDIA
            IF(ZDFL(JL,JK,IAZI)>0.0_JPRB) THEN
               ZATMP=ZCRT(JL,JK,IAZI)
               DO INC=1,INCDIM
                  ZATMP=ZATMP+ZCI(INC)*&
     &                   ZACC(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZDCI(INC)
               ENDDO
               ZCRT(JL,JK,IAZI)=ZATMP/ZDFL(JL,JK,IAZI)
            ELSE
               ZCRT(JL,JK,IAZI)=ZCRT(JL,JK+1,IAZI)
            ENDIF
         ENDDO

!* do saturation (Eq. (26) and (27) of Scinocca 2003)
!* -------------------------------------------------

         IF(GPTWO==3.0_JPRB) THEN
             DO INC=1,INCDIM
                ZCIN=ZCI(INC)
                ZCINC=1.0_JPRB/ZCIN
                DO JL=KIDIA,KFDIA
                   ZE1=ZCIN-ZUI(JL,JK,IAZI)
                   ZE2=GCSTAR*ZFCT(JL,JK)*ZE1
                   ZFLUXSQ=ZE2*ZE2*ZE1*ZCINC
                !  ZFLUXSQ=ZE2*ZE2*ZE1/ZCIN
                   ZDEP=ZACT(JL,INC,IAZI)*(ZFLUX(JL,INC,IAZI)**2-ZFLUXSQ)
                   IF(ZDEP>0.0_JPRB) THEN
                      ZFLUX(JL,INC,IAZI)=SQRT(ZFLUXSQ)
                   ENDIF
                ENDDO
             ENDDO
         ELSEIF(GPTWO==2.0_JPRB) THEN
             DO INC=1,INCDIM
                ZCIN=ZCI(INC)
                ZCINC=1.0_JPRB/ZCIN
                DO JL=KIDIA,KFDIA
                   ZFLUXS=GCSTAR*ZFCT(JL,JK)*&
     &                    (ZCIN-ZUI(JL,JK,IAZI))**2*ZCINC
                !  ZFLUXS=GCSTAR*ZFCT(JL,JK)*(ZCIN-ZUI(JL,JK,IAZI))**2/ZCIN
                   ZDEP=ZACT(JL,INC,IAZI)*(ZFLUX(JL,INC,IAZI)-ZFLUXS)
                   IF(ZDEP>0.0_JPRB) THEN
                      ZFLUX(JL,INC,IAZI)=ZFLUXS
                   ENDIF
                ENDDO
             ENDDO
         ENDIF
            
!* integrate spectrum
!* ------------------            

         DO INC=1,INCDIM
            ZCINC=ZDCI(INC)
            DO JL=KIDIA,KFDIA
               ZPU(JL,JK,IAZI)=ZPU(JL,JK,IAZI)+&
     &               ZACT(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZCINC
            ENDDO
         ENDDO

!* end JK do-loop
!* --------------

     ENDDO

!* end IAZIDIM do-loop
!* ---------------

ENDDO

! -----------------------------------------------------------------------------

!*       MAKE CORRECTION FOR CRITICAL-LEVEL MOMENTUM DEPOSITION
!*       ------------------------------------------------------

Z0P0=0._JPRB
ZRGPTS=1.0_JPRB/(RG*PTSTEP)
DO IAZI=1,IAZIDIM
   DO JL=KIDIA,KFDIA
     ZCNGL(JL)=0.0_JPRB
   ENDDO
   DO JK=2,ILAUNCH
      DO JL=KIDIA,KFDIA
         ZULM=ZCOSANG(IAZI)*PUM1(JL,JK)+ZSINANG(IAZI)*PVM1(JL,JK)-ZUL(JL,IAZI)
         ZDFL(JL,JK-1,IAZI)=ZDFL(JL,JK-1,IAZI)+ZCNGL(JL)
         ZDFT=MIN(ZDFL(JL,JK-1,IAZI),2.0_JPRB*(PAPM1(JL,JK-1)-PAPM1(JL,JK))*&
&                 (ZCRT(JL,JK-1,IAZI)-ZULM)*ZRGPTS)
         ZDFT=MAX(ZDFT,Z0P0)
         ZCNGL(JL)=(ZDFL(JL,JK-1,IAZI)-ZDFT)
         ZPU(JL,JK,IAZI)=ZPU(JL,JK,IAZI)-ZCNGL(JL)
       ENDDO
   ENDDO
ENDDO


!*       SUM CONTRIBUTION FOR TOTAL ZONAL AND MERIDIONAL FLUX
!*       ---------------------------------------------------

DO IAZI=1,IAZIDIM
   DO JK=ILAUNCH,2,-1
      DO JL=KIDIA,KFDIA
         PFLUXU(JL,JK)=PFLUXU(JL,JK)+ZPU(JL,JK,IAZI)*ZAZ_FCT*ZCOSANG(IAZI)
         PFLUXV(JL,JK)=PFLUXV(JL,JK)+ZPU(JL,JK,IAZI)*ZAZ_FCT*ZSINANG(IAZI)
      ENDDO
   ENDDO
ENDDO


!*    UPDATE U AND V TENDENCIES
!*    ----------------------------   

ZCONS1=1.0_JPRB/RCPD
DO JK=1,ILAUNCH
   DO JL=KIDIA, KFDIA
      ZDELP= RG/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))  
      ZE1=(PFLUXU(JL,JK+1)-PFLUXU(JL,JK))*ZDELP
      ZE2=(PFLUXV(JL,JK+1)-PFLUXV(JL,JK))*ZDELP   
      PTENU(JL,JK)=ZE1
      PTENV(JL,JK)=ZE2
   ENDDO
ENDDO

!--------------------------------------------------------------------------- 

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GWDRAG_WMS',1,ZHOOK_HANDLE)

END SUBROUTINE GWDRAG_WMS
