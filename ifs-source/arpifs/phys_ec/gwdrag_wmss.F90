! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GWDRAG_WMSS(YDEGWD,YDEGWWMS,KIDIA, KFDIA, KLON  , KLEV  , PTSTEP,&
                     & PTM1 , PUM1 , PVM1  , PAPM1 , PAPHM1, PGEO1 ,&
                     & PGELAT,PGAW ,&
                     & PTENU, PTENV, PFLUXU, PFLUXV)

!**** *GWDRAG_WMSS*  MASTER ROUTINE FOR NON-OROGRAPHIC GRAVITY WAVE DRAG
!                    SIMPLIFIED VERSION

!     Original Fortran Code by     J. SCINOCCIA
!     Rewritten in IFS format by   A. ORR         E.C.M.W.F.  August 2008
!     Simplified&optimised version P. BECHTOLD    E.C.M.W.F.  February/July 2009

!     PURPOSE
!     -------

!          THIS ROUTINE COMPUTES NON-OROGRAPHIC GRAVITY WAVE DRAG
!     AFTER SCINOCCA (2003) AND Mc LANDRESS AND SCINOCCIA (JAS 2005)
!     HYDROSTATIC NON-ROTATIONAL SIMPLIFIED VERSION OF THE 
!     WARNER AND MCINTYRE (1996) NON-OROGRAPHIC GRAVITY WAVE PARAMETERIZATION 
!     CONSTANTS HAVE BEEN OPTIMIZED FOLLOWING M. ERN ET AL. (ATMOS. CHEM. PHYS. 
!     2006)
! 
!     LAUNCH SPECTRUM - GENERALIZED DESAUBIES 
!     INCLUDES A CRITICAL-LEVEL CORRECTION THAT PREVENTS THE 
!     MOMEMTUM DEPOSITION IN EACH AZIMUTH FROM DRIVING THE FLOW TO SPEEDS 
!     FASTER THAN THE PHASE SPEED OF THE WAVES, I.E. WHEN WAVES BREAK THEY DRAG
!      THE MEAN FLOW TOWARDS THEIR PHASE SPEED - NOT PAST IT.  

!**   INTERFACE.
!     ----------
!          *GWDRAG_WMSS* IS CALLED FROM *CALLPAR*

!     MODIFICATIONS
!     -------------
!         
!  
! ------------------------------------------------------------------------------


USE PARKIND1,    ONLY : JPIM, JPRB
USE YOEGWWMS,    ONLY : TEGWWMS
USE YOMCST,      ONLY : RG, RD, RCPD, RPI, RA
USE YOMHOOK,     ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOEGWD,      ONLY : TEGWD

IMPLICIT NONE

!in
TYPE(TEGWD)       ,INTENT(INOUT):: YDEGWD
TYPE(TEGWWMS)     ,INTENT(INOUT):: YDEGWWMS
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON              ! horizontal dimension
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV              ! vertical levels
!INTEGER(KIND=JPIM),INTENT(IN) :: KLAUNCH           ! index for launch level selection
REAL(KIND=JPRB),INTENT(IN) :: PTSTEP               ! model time step
REAL(KIND=JPRB),INTENT(IN) :: PUM1(KLON,KLEV)      ! full model level zonal 
                                                   ! velocity (t-dt)
REAL(KIND=JPRB),INTENT(IN) :: PVM1(KLON,KLEV)      ! full model level meridional
                                                   ! velocity (t-dt)
REAL(KIND=JPRB),INTENT(IN) :: PTM1(KLON,KLEV)      ! full model level 
                                                   ! temperature (t-dt)
REAL(KIND=JPRB),INTENT(IN) :: PAPM1(KLON,KLEV)     ! full model level pressure 
                                                   ! (t-dt)
REAL(KIND=JPRB),INTENT(IN) :: PAPHM1(KLON,KLEV+1)  ! half-model level pressure 
                                                   ! (t-dt)
REAL(KIND=JPRB),INTENT(IN) :: PGEO1(KLON,KLEV)     ! full mod.level geopotential
REAL(KIND=JPRB),INTENT(IN) :: PGELAT(KLON)         ! latitude
REAL(KIND=JPRB),INTENT(IN) :: PGAW(KLON)           !normalised gaussian quadrature weight/number of longitude points
                                                   !local sub-area == 4*RPI*RA**2 * PGAW

!inout
REAL(KIND=JPRB),INTENT(OUT):: PTENU(KLON,KLEV)     ! full-model level zonal 
                                                   ! momentum tendency
REAL(KIND=JPRB),INTENT(OUT):: PTENV(KLON,KLEV)     ! full-model level meridional
                                                   ! momentum tendency
REAL(KIND=JPRB),INTENT(OUT):: PFLUXU(KLON,KLEV+1)  ! zonal component of vertical
                                                   ! momentum flux (Pa)  
REAL(KIND=JPRB),INTENT(OUT):: PFLUXV(KLON,KLEV+1)  ! meridional component of 
                                                   ! vertical momentum flux (Pa)
!work
INTEGER(KIND=JPIM), PARAMETER  :: IAZIDIM=4   !number of azimuths
INTEGER(KIND=JPIM), PARAMETER  :: INCDIM=20   !number of discretized c
                                              !spectral elements in launch spectrum

REAL(KIND=JPRB) :: ZUHM1(KLON,KLEV)           !half-model level zonal velocity
REAL(KIND=JPRB) :: ZVHM1(KLON,KLEV)           !half-model level meridional velocity
REAL(KIND=JPRB) :: ZBVFHM1(KLON,KLEV)         !half-model level Brunt-Vaisalla 
                                              !frequency
REAL(KIND=JPRB) :: ZRHOHM1(KLON,KLEV)         !half-model level density
REAL(KIND=JPRB) :: ZX(INCDIM)                 !coordinate transformation
REAL(KIND=JPRB) :: ZCI(INCDIM)                !phase speed element
REAL(KIND=JPRB) :: ZDCI(INCDIM)
REAL(KIND=JPRB) :: ZUI(KLON,KLEV,IAZIDIM)     !intrinsic velocity
REAL(KIND=JPRB) :: ZUL(KLON,IAZIDIM)          !velocity in azimuthal direction 
                                              !at launch level
REAL(KIND=JPRB) :: ZBVFL(KLON)                !buoyancy at launch level
REAL(KIND=JPRB) :: ZCOSANG(IAZIDIM)           !cos of azimuth angle
REAL(KIND=JPRB) :: ZSINANG(IAZIDIM)           !sin of azimuth angle
REAL(KIND=JPRB) :: ZFCT(KLON,KLEV)
REAL(KIND=JPRB) :: ZFNORM(KLON)               !normalisation factor (A)
REAL(KIND=JPRB) :: ZCI_MIN(KLON,IAZIDIM)
REAL(KIND=JPRB) :: ZTHM1(KLON,KLEV)           !temperature on half-model levels
REAL(KIND=JPRB) :: ZFLUX(KLON,INCDIM,IAZIDIM) !momentum flux at each vertical 
                                              !level and azimuth 
REAL(KIND=JPRB) :: ZPU(KLON,KLEV,IAZIDIM)     !momentum flux
REAL(KIND=JPRB) :: ZDFL(KLON,KLEV,IAZIDIM)
REAL(KIND=JPRB) :: ZACT(KLON,INCDIM,IAZIDIM)  !if =1 then critical level 
                                              !encountered  
REAL(KIND=JPRB) :: ZACC(KLON,INCDIM,IAZIDIM)
REAL(KIND=JPRB) :: ZCRT(KLON,KLEV,IAZIDIM)
REAL(KIND=JPRB) :: ZPUI(KLON)

INTEGER(KIND=JPIM) :: INC, JK, JL, IAZI
REAL(KIND=JPRB) :: ZRADTODEG, ZPTWO, ZGELATDEG
REAL(KIND=JPRB) :: ZCIMIN, ZCIMAX
REAL(KIND=JPRB) :: ZGAM, ZPEXP, ZXMAX, ZXMIN, ZXRAN, ZDX, ZX1, ZX2, ZDXA, ZDXB, ZDXS
REAL(KIND=JPRB) :: ZANG, ZAZ_FCT, ZNORM, ZANG1, ZTX
REAL(KIND=JPRB) :: ZU, ZCIN
REAL(KIND=JPRB) :: ZCIN4, ZBVFL4, ZCINC
REAL(KIND=JPRB) :: ZATMP, ZFLUXS, ZDEP, ZULM, ZDFT
REAL(KIND=JPRB) :: ZMS_L,ZMS, Z0P5, Z0P0 
REAL(KIND=JPRB) :: ZGAUSS(KLON), ZFLUXLAUN(KLON), ZCNGL(KLON)
REAL(KIND=JPRB) :: ZCONS1,ZCONS2,ZDELP,ZRGPTS
REAL(KIND=JPRB) :: ZEXP1, ZDIV1, ZDFLI, ZPRULM
REAL(KIND=JPRB) :: ZTHMI1, ZTHMIB, ZDELPI, ZBVFHMI1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GWDRAG_WMSS',0,ZHOOK_HANDLE)
ASSOCIATE(GSSEC=>YDEGWD%GSSEC, &
 & GCSTAR=>YDEGWWMS%GCSTAR, GFLUXLAUN=>YDEGWWMS%GFLUXLAUN, &
 & GGAUSSA=>YDEGWWMS%GGAUSSA, GGAUSSB=>YDEGWWMS%GGAUSSB, LOZPR=>YDEGWWMS%LOZPR, &
 & NGAUSS=>YDEGWWMS%NGAUSS, NLAUNCH=>YDEGWWMS%NLAUNCH)
!--------------------------------------------------------------------------

!*       INPUT PARAMETERS
!*       ---------------- 

ZRADTODEG=57.29577951_JPRB

!m_star
ZMS_L=2000.0_JPRB
ZMS=2*RPI/ZMS_L    
ZPTWO=2.0_JPRB          !assume GPTWO=2 in sugwwms.F90
                        !2*p, p is the exponent of omega 

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

ZCIMIN=0.25_JPRB
ZCIMAX=100.0_JPRB
ZGAM=0.25_JPRB

ZPEXP=ZPTWO/2.0_JPRB

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
    ZDELPI = 1.0_JPRB/ZDELP
    ZTHMI1 = 1.0_JPRB/ZTHM1(JL,JK)
    ZTHMIB = (PTM1(JL,JK)-PTM1(JL,JK-1))*ZTHMI1
    ZRHOHM1(JL,JK)=PAPHM1(JL,JK)*ZCONS1*ZTHMI1
    ZBVFHM1(JL,JK) = ZCONS2*(ZTHMI1+RCPD*ZTHMIB*ZDELPI)
    ZBVFHM1(JL,JK)=MAX(ZBVFHM1(JL,JK),GSSEC)
    ZBVFHM1(JL,JK)=SQRT(ZBVFHM1(JL,JK))
  ENDDO
ENDDO

!*       SET UP AZIMUTH DIRECTIONS AND SOME TRIG FACTORS
!*       -----------------------------------------------     

ZANG=2*RPI/IAZIDIM
ZAZ_FCT=1.0

! get normalization factor to ensure that the same amount of momentum
! flux is directed (n,s,e,w) no mater how many azimuths are selected.
! note, however, the code below assumes a symmetric distribution of
! of azimuthal directions (ie 4,8,16,32,...)

ZNORM=0.0
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
ZX1=ZXRAN/(EXP(ZXRAN/ZGAM)-1.0_JPRB)
ZX2=ZXMIN-ZX1

DO INC=1,INCDIM
  ZTX=REAL(INC-1)*ZDX+ZXMIN
  ZEXP1 = EXP((ZTX-ZXMIN)/ZGAM)
  ZX(INC)=ZX1*ZEXP1+ZX2                           !Eq. 29 of Scinocca 2003
  ZCI(INC)=1.0_JPRB/ZX(INC)                       !Eq. 28 of Scinocca 2003
  ZDCI(INC)=ZCI(INC)**2*(ZX1/ZGAM)*ZEXP1*ZDX      !Eq. 30 of Scinocca 2003
ENDDO


!*   DEFINE INTRINSIC VELOCITY (RELATIVE TO LAUNCH LEVEL VELOCITY) U(Z)-U(Zo)
!    AND COEFFICINETS
!*   --------------------------------------------------------------------------

DO IAZI=1,IAZIDIM
  DO JL=KIDIA,KFDIA
    ZUL(JL,IAZI)=ZCOSANG(IAZI)*ZUHM1(JL,NLAUNCH)+ZSINANG(IAZI)*ZVHM1(JL,NLAUNCH)
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  ZBVFL(JL)=ZBVFHM1(JL,NLAUNCH)
ENDDO


DO IAZI=1,IAZIDIM
  DO JK=2,NLAUNCH
    DO JL=KIDIA,KFDIA
      ZU=ZCOSANG(IAZI)*ZUHM1(JL,JK)+ZSINANG(IAZI)*ZVHM1(JL,JK)
      ZUI(JL,JK,IAZI)=ZU-ZUL(JL,IAZI)
    ENDDO
  ENDDO
ENDDO

!*       DEFINE RHO(Zo)/N(Zo)
!*       ------------------- 

DO JK=2,NLAUNCH
  DO JL=KIDIA,KFDIA
    ZBVFHMI1 = 1.0_JPRB/ZBVFHM1(JL,JK)
    ZFCT(JL,JK)=ZRHOHM1(JL,JK)*ZBVFHMI1
  ENDDO
ENDDO

!*       SET LAUNCH MOMENTUM FLUX SPECTRAL DENSITY
!*       ----------------------------------------- 

! Eq. (25) of Scinocca 2003 (not including the 'A' component), and with U-Uo=0
! do this for only one azimuth since it is identical to all azimuths, and
! it will be renormalized

! s=1 case : standard, NSLOPE in sugwwms.F90 supposed to be 1
DO INC=1,INCDIM
  ZCIN=ZCI(INC)
  ZCIN4=(ZMS*ZCIN)**4
  DO JL=KIDIA,KFDIA
    ZBVFL4=ZBVFL(JL)**4
    ZDIV1 = 1.0_JPRB/(ZBVFL4+ZCIN4)
    ZFLUX(JL,INC,1)=ZFCT(JL,NLAUNCH)*ZBVFL4*ZCIN*ZDIV1
    ZACT(JL,INC,1)=1.0_JPRB
  ENDDO
ENDDO

!*       NORMALIZE LAUNCH MOMENTUM FLUX
!*       ------------------------------

! (rho x F^H = rho_o x F_p^total)

! integrate (ZFLUX x dX)
DO INC=1,INCDIM
  ZCINC=ZDCI(INC)
  DO JL=KIDIA,KFDIA
    ZPU(JL,NLAUNCH,1)=ZPU(JL,NLAUNCH,1)+ZFLUX(JL,INC,1)*ZCINC
  ENDDO
ENDDO

!*       NORMALIZE GFLUXLAUN TO INCLUDE SENSITIVITY TO PRECIPITATION
!*       -----------------------------------------------------------

! A=ZFNORM in Scinocca 2003.  A is independent of height.
ZDXA=1.0_JPRB/29.E3_JPRB
ZDXB=1.0_JPRB/3.5E3_JPRB
DO JL=KIDIA,KFDIA       
  ZDX=MAX(1.E2_JPRB,2*RA*SQRT(RPI*PGAW(JL))) !grid resolution (m)
   !Scaling factor for launch flux depending on grid resolution
  ZDXS=1.0_JPRB-MIN(1.0_JPRB,ATAN((MAX(1.0_JPRB/ZDX,ZDXA)-ZDXA)/(ZDXB-ZDXA)))
  ZFLUXLAUN(JL)=GFLUXLAUN*ZDXS
  ZPUI(JL) = 1.0_JPRB/ZPU(JL,NLAUNCH,1)
  ZFNORM(JL)=ZFLUXLAUN(JL)*ZPUI(JL)
ENDDO

IF (LOZPR) THEN
  IF (NGAUSS==2) THEN           
    DO JL=KIDIA,KFDIA       
      ZGELATDEG=PGELAT(JL)*ZRADTODEG 
      ZGAUSS(JL)=GGAUSSB(1)*EXP((-ZGELATDEG*ZGELATDEG)/(2*GGAUSSA*GGAUSSA))
      ZFLUXLAUN(JL)=(1.0_JPRB+ZGAUSS(JL))*ZFLUXLAUN(JL)
      ZFNORM(JL)=ZFLUXLAUN(JL)*ZPUI(JL)
    ENDDO
  ENDIF
ENDIF

DO IAZI=1,IAZIDIM
  DO JL=KIDIA,KFDIA
    ZPU(JL,NLAUNCH,IAZI)=ZFLUXLAUN(JL)
  ENDDO
ENDDO

!*       ADJUST CONSTANT ZFCT
!*       --------------------
DO JK=2,NLAUNCH
  DO JL=KIDIA,KFDIA
    ZFCT(JL,JK)=ZFNORM(JL)*ZFCT(JL,JK)
  ENDDO
ENDDO

!*       RENORMALIZE EACH SPECTRAL ELEMENT IN ONE AZIMUTH
!*       ------------------------------------------------
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

  DO JK=NLAUNCH-1,2,-1

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

!* integrate to get critical-level contribution to mom deposition on this level,
! i.e. ZACC=1
!* -----------------------------------------------------------------------------

    DO INC=1,INCDIM
      ZCINC=ZDCI(INC)
      DO JL=KIDIA,KFDIA
        ZDFL(JL,JK,IAZI)=ZDFL(JL,JK,IAZI)+&
         &  ZACC(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZCINC
      ENDDO
    ENDDO

!* get weighted average of phase speed in layer
!* --------------------------------------------

    DO JL=KIDIA,KFDIA
      IF (ZDFL(JL,JK,IAZI)>0.0_JPRB) THEN
        ZATMP=ZCRT(JL,JK,IAZI)
        DO INC=1,INCDIM
          ZATMP=ZATMP+ZCI(INC)*&
           & ZACC(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZDCI(INC)
        ENDDO
        ZDFLI = 1.0_JPRB/ZDFL(JL,JK,IAZI)
        ZCRT(JL,JK,IAZI)=ZATMP*ZDFLI
      ELSE
        ZCRT(JL,JK,IAZI)=ZCRT(JL,JK+1,IAZI)
      ENDIF
    ENDDO

!* do saturation (Eq. (26) and (27) of Scinocca 2003)
!* -------------------------------------------------
!  case GPTWO=2=ZPTWO in sugwwms.F90

    DO INC=1,INCDIM
      ZCIN=ZCI(INC)
      ZCINC=1.0_JPRB/ZCIN
      DO JL=KIDIA,KFDIA
        ZFLUXS=GCSTAR*ZFCT(JL,JK)*(ZCIN-ZUI(JL,JK,IAZI))**2*ZCINC
     !  ZFLUXS=GCSTAR*ZFCT(JL,JK)*(ZCIN-ZUI(JL,JK,IAZI))**2/ZCIN
        ZDEP=ZACT(JL,INC,IAZI)*(ZFLUX(JL,INC,IAZI)-ZFLUXS)
        IF (ZDEP>0.0_JPRB) THEN
          ZFLUX(JL,INC,IAZI)=ZFLUXS
        ENDIF
      ENDDO
    ENDDO
            
!* integrate spectrum
!* ------------------            

    DO INC=1,INCDIM
      ZCINC=ZDCI(INC)
      DO JL=KIDIA,KFDIA
        ZPU(JL,JK,IAZI)=ZPU(JL,JK,IAZI)+&
         & ZACT(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZCINC
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

Z0P0=0.0_JPRB
ZRGPTS=1.0_JPRB/(RG*PTSTEP)
DO IAZI=1,IAZIDIM
  DO JL=KIDIA,KFDIA
    ZCNGL(JL)=0.0_JPRB
  ENDDO
  DO JK=2,NLAUNCH
    DO JL=KIDIA,KFDIA
      ZULM=ZCOSANG(IAZI)*PUM1(JL,JK)+ZSINANG(IAZI)*PVM1(JL,JK)-ZUL(JL,IAZI)
      ZDFL(JL,JK-1,IAZI)=ZDFL(JL,JK-1,IAZI)+ZCNGL(JL)
      ZPRULM = 2.0_JPRB*(PAPM1(JL,JK-1)-PAPM1(JL,JK)) &
       & * (ZCRT(JL,JK-1,IAZI)-ZULM)*ZRGPTS
      ZDFT=MIN(ZDFL(JL,JK-1,IAZI),ZPRULM)
      ZDFT=MAX(ZDFT,Z0P0)
      ZCNGL(JL)=(ZDFL(JL,JK-1,IAZI)-ZDFT)
      ZPU(JL,JK,IAZI)=ZPU(JL,JK,IAZI)-ZCNGL(JL)
    ENDDO
  ENDDO
ENDDO


!*       SUM CONTRIBUTION FOR TOTAL ZONAL AND MERIDIONAL FLUX
!*       ---------------------------------------------------

DO IAZI=1,IAZIDIM
  DO JK=NLAUNCH,2,-1
    DO JL=KIDIA,KFDIA
      PFLUXU(JL,JK)=PFLUXU(JL,JK)+ZPU(JL,JK,IAZI)*ZAZ_FCT*ZCOSANG(IAZI)
      PFLUXV(JL,JK)=PFLUXV(JL,JK)+ZPU(JL,JK,IAZI)*ZAZ_FCT*ZSINANG(IAZI)
    ENDDO
  ENDDO
ENDDO


!*    UPDATE U AND V TENDENCIES
!*    ----------------------------   

DO JK=1,NLAUNCH
  DO JL=KIDIA, KFDIA
    ZDELP= 1.0_JPRB/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))  
    PTENU(JL,JK)=(PFLUXU(JL,JK+1)-PFLUXU(JL,JK))*RG*ZDELP
    PTENV(JL,JK)=(PFLUXV(JL,JK+1)-PFLUXV(JL,JK))*RG*ZDELP   
  ENDDO
ENDDO

!--------------------------------------------------------------------------- 

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GWDRAG_WMSS',1,ZHOOK_HANDLE)

END SUBROUTINE GWDRAG_WMSS
