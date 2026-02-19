! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DIAG_TURB&
 & ( YDEPHY, YDECUMF,KIDIA,  KFDIA,  KLON,   KLEV,  &
 &   LEGWWMS, KLAUNCH,KTYPE, KCTOP, &
 &   POROG,  POROGSTD,       PGAW,&
 &   PAP,    PAPH,   PGEO,   PGEOH,&
 &   PT,     PQ,     PU,     PV,     PVERVEL, &
 &   PTDX,   PTDY,   PUDX,   PVDX,   PVOR,   PDIV, PDISSCU, &
 &   PTEND_UDIF, PTEND_VDIF, PTEND_T,PRI,    PTENOGW_U, PTENOGW_V, PEDRP )


!    MARTINA BRAMBERGER (DLR)  & P.BECHTOLD     E.C.M.W.F.     11/2018

!    -------

!    TO PROVIDE TURBULENCE DIAGNOSTICS 

!    INTERFACE
!    ---------
!    THIS ROUTINE IS CALLED FROM *CALLPAR*.

!    METHOD.
!    --------
!    See IFS Tech. Memo No. 874
!    Eddy dissipation rates derived from climatological projections
!    following Sharman and Pearson, JAMC 2017, 317-337
!    Computes Model turbulent dissipation=DISS, ZEL1=Ellrod1, ZF3=3d Frontogen function
!    and convective gravity wave dissipation
!    Deduces CAT and mountain wave turbulence preferably using DISS and/or Ellrod1

!    NOTA: due to computational reasons we decided to define the operational CAT=DISS. 
!          Optionally we can compute a CAT based on Ellrod1 which is stored in 'MWT' Grib field.

!    PARAMETER     DESCRIPTION                                   UNITS 
!    ---------     -----------                                   ----- 
!    INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        CONVECTION TYPE
!    *KCTOP*        CONVECTION TOP LEVEL
!    *KLAUNCH*      LAUNCH-LEVEL OF NON-OROG GRAVITY WAVES

!    INPUT PARAMETERS (LOGICAL):
!    *LEGWWMS*      SWITCH IF NON-OROG GRAVITY WAVE INPUT AVAILABLE
 
!     INPUT PARAMETERS (REAL)
!    *POROG*        HEIGHT OF OROGRAPHY (M)
!    *POROGSTD*     STANDARD DEVIATION OF OROGRAPHY (M)
!    *PGAW*         NORMALISED GAUSSIAN QUADRATURE WEIGHT / NUMBER OF LONGITUDE POINTS
!                           LOCAL SUB-AREA == 4*RPI*RA**2 * PGAW
!    *PAP*          PRESSURE                                      PA
!    *PAPH*         PRESSURE HALF LEVEL                           PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PT*           TEMPERATURE                                    K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG
!    *PU*           ZONAL WIND                                    M/S 
!    *PV*           MERIDIONAL WIND                               M/S 
!    *PVERVEL*      OMEGA                                        PA/S 
!    *PTDX*         ZONAL  GRADIENT T                             K/M
!    *PTDY*         MERIOD GRADIENT T                             K/M
!    *PUDX*         ZONAL  GRADIENT U                             1/S 
!    *PVDX*         ZONAL  GRADIENT V                             1/S 
!    *PRI*          RICHARDSON NUMBER                  
!    *PVOR*         VORTICITY                                     1/S 
!    *PDIV*         DIVERGENCE                                    1/S 
!    *PDIV*         DIVERGENCE                                    1/S 
!    *PDISSCU*      DISIIPATION FROM CU CONVECTION            (M2/S3)^1/3
!    *PTEND_UDIF*   DIFFUSION TENDENCY U                          M/S2
!    *PTEND_VDIF*   DIFFUSION TENDENCY V                          M/S2
!    *PTEND_T*      CONVECTIVE TENDENCY T                         K/S

!    *PTENOGW_U*    U GW-MOMENTUM TENDENCIES (FLUXES)             M/S2
!    *PTENOGW_V*    V GW-MOMENTUM TENDENCIES (FLUXES)             M/S2
!    NOTA: ZUDY= PUDX-PVOR    ZVDY=PDIV-PVDX

!    OUTPUT PARAMETERS (REAL) 
!    *PEDRP*        CLEAR AIR TURBULENCE EDR                   (M2/S3)^1/3
!                   MOUNTAIN WAVE TURBULENCE EDR               (M2/S3)^1/3

!    NOTA: WE HAVE 2 GRIB VALUES AVAILABLE

!          MODIFICATIONS
!          -------------

!----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RD, RCPD, RPI, RA
USE YOEPHY   , ONLY : TEPHY
USE YOECUMF  , ONLY : TECUMF

IMPLICIT NONE

TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY
TYPE(TECUMF)      ,INTENT(IN)    :: YDECUMF
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAUNCH
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON)
LOGICAL           ,INTENT(IN)    :: LEGWWMS

REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGSTD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTDX(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTDY(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUDX(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDX(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOR(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDISSCU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEND_UDIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEND_VDIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEND_T(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENOGW_U(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENOGW_V(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEDRP(KLON,KLEV,2) 

INTEGER(KIND=JPIM) :: JL, JK, JP, IK
INTEGER(KIND=JPIM), PARAMETER :: NP=5

!Parameters to compute
LOGICAL :: LLCOMP_ELLR=.FALSE.,LLCOMP_F3D=.FALSE.
!Climatological values for Eddy dissipation mapping
REAL(KIND=JPRB)   , PARAMETER    :: ZC1=-2.572_JPRB, ZC2=0.5087_JPRB  
!Parameter for calculating convective gravity wave drag
REAL(KIND=JPRB)   , PARAMETER    :: ZC_EPS=1.0_JPRB !0.93_JPRB
!Lognormal distribution fit parameters: expected value and standard deviation
REAL(KIND=JPRB)                  :: ZEV(NP), ZSD(NP)
!Mapping parameters onto eddy dissipation
REAL(KIND=JPRB)                  :: ZA(NP), ZB(NP), ZEA(NP)

REAL(KIND=JPRB)    :: ZRG, ZDZ, ZRDOCP, ZSHR, ZDEF, ZUV950, ZUV850, ZEXN

!Indices
REAL(KIND=JPRB) :: ZEL1(KLON,KLEV) , ZF3(KLON,KLEV), ZTDISS(KLON,KLEV), ZCGWD(KLON,KLEV)

REAL(KIND=JPRB) :: ZDUDZ(KLON,KLEV), ZDVDZ(KLON,KLEV), ZDTHETADZ(KLON,KLEV)
REAL(KIND=JPRB) :: ZDUDY(KLON,KLEV), ZDVDY(KLON,KLEV)
REAL(KIND=JPRB) :: ZDTHETADX(KLON,KLEV), ZDTHETADY(KLON,KLEV)
REAL(KIND=JPRB) :: ZDWDX, ZDWDY, ZDWDZ, ZTKE
REAL(KIND=JPRB) :: ZTHETA(KLON,KLEV), ZDELTHETA(KLON,KLEV), ZF3D(KLON,KLEV) 
REAL(KIND=JPRB) :: ZDX(KLON), ZSMW(KLON), ZHEAT(KLON)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIAG_TURB',0,ZHOOK_HANDLE)
ASSOCIATE(NJKT8=>YDECUMF%NJKT8,NJKT4=>YDECUMF%NJKT4,NJKT3=>YDECUMF%NJKT3,NJKT5=>YDECUMF%NJKT5,&
&NJKT1=>YDECUMF%NJKT1, LDIAGTURBGRAD_EC=>YDEPHY%LDIAGTURBGRAD_EC)

! set fit values for Ellrod1, F3D, MWT3, GWD, TURB
!TCo639
ZEV = (/-15.4_JPRB, -16.3_JPRB, -6.15_JPRB,-2.2_JPRB,-3.3_JPRB/)
ZSD = (/1.25_JPRB, 1.80_JPRB, 2.50_JPRB,0.52_JPRB,0.6_JPRB/)
!TCo1279
ZEV = (/-15.4_JPRB, -16.3_JPRB, -5.80_JPRB,-2.4_JPRB,-3.3_JPRB/)
ZSD = (/1.40_JPRB, 1.80_JPRB, 2.50_JPRB,0.65_JPRB,0.6_JPRB/)
ZSD=SQRT(ZSD)
! compute mapping factors onto eddy dissipation
DO JP=1,NP
  ZB(JP)=ZC2/ZSD(JP)
  ZA(JP)=ZC1-ZB(JP)*ZEV(JP) 
  ZEA(JP)=EXP(ZA(JP))
ENDDO 

ZRG=1.0_JPRB/RG
ZRDOCP=RD/RCPD

! Disable computation of Ellrod and F3D if horizontal gradients are not enabled
IF(.NOT.LDIAGTURBGRAD_EC) THEN
  LLCOMP_ELLR=.FALSE.
  LLCOMP_F3D=.FALSE.
ENDIF

! Initialisations
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PEDRP(JL,JK,1) =0.0_JPRB
    PEDRP(JL,JK,2) =0.0_JPRB
    ZTDISS(JL,JK)=0.0_JPRB
    ZCGWD(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!------------------------------------------------------------------------------------------------
! calculate turbulent dissipation using diffusion tendencies (contain orographic effects) 
!------------------------------------------------------------------------------------------------

  DO JK=KLEV,NJKT8,-1
    DO JL=KIDIA,KFDIA
       ZTDISS(JL,JK)=ABS(PU(JL,JK)*PTEND_UDIF(JL,JK)+PV(JL,JK)*PTEND_VDIF(JL,JK))**0.3333
       ZTDISS(JL,JK)=ZTDISS(JL,JK)+PDISSCU(JL,JK)
       ZTDISS(JL,JK)=ZEA(5)*ZTDISS(JL,JK)**ZB(5)
       PEDRP(JL,JK,1)=0.66_JPRB*MIN(ZTDISS(JL,JK),1.0_JPRB)
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------------------------

IF(LLCOMP_ELLR.OR.LLCOMP_F3D) THEN
!-----------------------------------------------------------------------------------------
! calculate meridional wind gradients
!----------------------------------------------------------------------------------------
  DO JK=KLEV,NJKT8,-1
    DO JL=KIDIA,KFDIA
      ZDUDY(JL,JK) = PVDX(JL,JK)-PVOR(JL,JK)    
      ZDVDY(JL,JK) = PDIV(JL,JK)-PUDX(JL,JK)
    ENDDO
  ENDDO
!------------------------------------------------------------------------------------------
! calculate vertical gradients and static stability
!-----------------------------------------------------------------------------------------
  IK=NJKT8-1
  DO JL=KIDIA,KFDIA
    ZTHETA(JL,IK)=PT(JL,IK)*(1.E5_JPRB/PAP(JL,IK))**ZRDOCP
  ENDDO
  DO JK=NJKT8,KLEV
    IK=JK-1
    DO JL=KIDIA,KFDIA
      ZDZ = RG/(PGEO(JL,IK)-PGEO(JL,JK))
      ZDUDZ(JL,JK) = (PU(JL,IK)-PU(JL,JK))*ZDZ  
      ZDVDZ(JL,JK) = (PV(JL,IK)-PV(JL,JK))*ZDZ  
      ZEXN=(1.E5_JPRB/PAP(JL,JK))**ZRDOCP
      ZDTHETADX(JL,JK) = PTDX(JL,JK)*ZEXN
      ZDTHETADY(JL,JK) = PTDY(JL,JK)*ZEXN
      ZTHETA(JL,JK) = PT(JL,JK)*ZEXN
      ZDTHETADZ(JL,JK) = (ZTHETA(JL,IK)-ZTHETA(JL,JK))*ZDZ
    ENDDO
  ENDDO
ENDIF

IF(LLCOMP_ELLR) THEN
!------------------------------------------------------------------------------------------------
! calculate Ellrod Index 1
!------------------------------------------------------------------------------------------------
!     --- Computes Ellrod indices (Ellrod and Knapp, Wea. Forecasting, 7, 1992)
!     --- on the input grid.
!     --- Note should be evaluated on constant height (z) surfaces.
!     --- Ref: Ellrod, G. P., and D. L. Knapp, 1992: An objective clear-air
!     --- turbulence forecasting technique: Verification and operational
!     --- use.  Wea. Forecasting, 7 

! need ZDUDX,DUDY,DVDX,DVDY

  DO JK=KLEV,NJKT8,-1
    DO JL=KIDIA,KFDIA
! calculate vertical wind shear
      ZSHR = SQRT(ZDUDZ(JL,JK)**2 + ZDVDZ(JL,JK)**2)   ! 1/s ! calculate deformation 
      ZDEF = SQRT((PVDX(JL,JK) + ZDUDY(JL,JK))**2 + (PUDX(JL,JK) - ZDVDY(JL,JK))**2) 
! calculate Ellrod Index 1
      ZEL1(JL,JK) = ZSHR*ZDEF               ! 1/sec^2
! map index onto eddy dissipation rate
      ZEL1(JL,JK) =ZEA(1)*ZEL1(JL,JK)**ZB(1)
! define CAT from Ellrod only or Ellrod and Ri
    ! PEDRP(JL,JK,2)=0.70*ZEL1(JL,JK)
    ! scaling factor was 0.95 in raw output
      PEDRP(JL,JK,2)=0.86_JPRB*ZEL1(JL,JK)*MIN(1.0_JPRB,0.5_JPRB/MAX(1.E-6_JPRB,PRI(JL,JK)))
      PEDRP(JL,JK,2)=MIN(PEDRP(JL,JK,1),1.0_JPRB)
    ENDDO  ! JL loop
  ENDDO  ! JK loop

ENDIF

!   --------------------------------------------------------------------------------------
IF(LLCOMP_F3D) THEN
!-------------------------------------------------------------------------------------------------
!  calculate 3D Fronotgenesis Function
!-------------------------------------------------------------------------------------------------

! --- Computes 3D frontogenesis function F on a const z surface.
! --- F=(1/del|theta|)*{-(dtheta/dx)*[(dtheta/dx)*(du/dx)+(dtheta/dy)*(dv/dx)
!                                    +(dtheta/dz)*(dw/dx)*]
!                       -(dtheta/dy)*[(dtheta/dx)*(du/dy)+(dtheta/dy)*(dv/dy)
!                                    +(dtheta/dz)*(dw/dy)]
!                       -(dtheta/dz)*[(dtheta/dx)*(du/dz)+(dtheta/dy)*(dv/dz)
!                                    +(dtheta/dz)*(dw/dz)] }
! --- ref: Bluestein vol.2, (2.3.20,21)
! --- Note input U,V are grid relative.  Output is in ZF3.

! We need  DTHETADX, DTHETADY, 
! set dw/dx and dw/dy to 0 as their gradients are negligible on a 8km hor. grid
! or set them to w/dx=w/dy
  ZDWDX=0.0_JPRB
  ZDWDY=0.0_JPRB

  DO JL=KIDIA,KFDIA
! grid size
     ZDX(JL)=2*RA*SQRT(RPI*PGAW(JL))
     ZDX(JL)=MAX(ZDX(JL),1.E2_JPRB)
! wind and orography depndent scaling factor for MW
     ZUV950=SQRT(PU(JL,NJKT4)*PU(JL,NJKT4)+PV(JL,NJKT4)*PV(JL,NJKT4))
     ZUV850=SQRT(PU(JL,NJKT3)*PU(JL,NJKT3)+PV(JL,NJKT3)*PV(JL,NJKT3))
     ZSMW(JL)=MAX(ZUV950,ZUV850)*MIN(MAX(0.0_JPRB,POROG(JL)),2750._JPRB)
! use only when std of surface orography is significant (nota: GWD param active if POROGSTD>50)
     IF(POROGSTD(JL)<10.0_JPRB) ZSMW(JL)=0.0_JPRB
  ENDDO

  DO JK=KLEV,NJKT8,-1
    DO JL=KIDIA,KFDIA

! test using also w

!     ZRHO=PAP(JL,JK)/(RD*PT(JL,JK))
!     ZDWDX=-PVERVEL(JL,JK)*ZRG/(ZDX(JL)*ZRHO)
!     ZDWDY=0.

      ZDWDZ = -PDIV(JL,JK) ! by continuity
!     --- Compute |del(theta)|
      ZDELTHETA(JL,JK) = SQRT(ZDTHETADX(JL,JK)**2 + ZDTHETADY(JL,JK)**2 + ZDTHETADZ(JL,JK)**2)
!     --- Compute the 3d gradients
      ZF3D(JL,JK)  = -ZDTHETADX(JL,JK)*(ZDTHETADX(JL,JK)*PUDX(JL,JK)+ZDTHETADY(JL,JK)*PVDX(JL,JK)+ZDTHETADZ(JL,JK)*ZDWDX)&
        & -ZDTHETADY(JL,JK)*(ZDTHETADX(JL,JK)*ZDUDY(JL,JK)+ZDTHETADY(JL,JK)*ZDVDY(JL,JK)+ZDTHETADZ(JL,JK)*ZDWDY) &
        & -ZDTHETADZ(JL,JK)*(ZDTHETADX(JL,JK)*ZDUDZ(JL,JK)+ZDTHETADY(JL,JK)*ZDVDZ(JL,JK)+ZDTHETADZ(JL,JK)*ZDWDZ)
      IF (ZDELTHETA(JL,JK) < 1.0E-12_JPRB) THEN 
         ZF3D(JL,JK)=0._JPRB
      ELSE 
         ZF3D(JL,JK) = ZF3D(JL,JK)/ZDELTHETA(JL,JK)
      ENDIF
!      NOTE: Here using the absolute value of the frontogenesis fn.  I.e.,
!      frontogenesis F3D>0 and frontolysis F3d<0 are treated as equally likely
!      to produce turbulence.
      ZF3(JL,JK)=ABS(ZF3D(JL,JK)) 

! For mountain wave turb multiply ZF3 by "surface parameter=low-lev windspeed x orography"

      PEDRP(JL,JK,2)=ZF3(JL,JK)*ZSMW(JL)
! Mountain wave turbulence
      PEDRP(JL,JK,2)=ZEA(3)*PEDRP(JL,JK,2)**ZB(3)
      PEDRP(JL,JK,2)=MIN(PEDRP(JL,JK,2),1.0_JPRB)
      ZF3(JL,JK)  =ZEA(2)*ZF3(JL,JK)**ZB(2)
! arithmetic mean Ellrod F3D
!     PEDRP(JL,JK,2)=0.5_JPRB*(ZEL1(JL,JK)+ZF3(JL,JK))

    ENDDO
  ENDDO
ENDIF
!-------------------------------------------------------------------------------------------------

!  calculate convective gravity wave drag based on subgrid scaled non-orographic fluxes/tendencies
!-------------------------------------------------------------------------------------------------
! previously use formulation of Chun and Baik 1998

IF(LEGWWMS) THEN
  ! calculate wave stress above convective clouds based on non-orographic wave flux
  ! add to CAT

  ZHEAT(:)=0.0_JPRB
  DO JK=NJKT5,NJKT1,-1
    DO JL=KIDIA,KFDIA
      IF(KTYPE(JL) == 1) THEN
         ZHEAT(JL)=ZHEAT(JL)+PTEND_T(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG*RCPD
      ENDIF
    ENDDO
  ENDDO
  DO JK=KLEV-1,NJKT8,-1
    DO JL=KIDIA,KFDIA
      IF (JK <=KLAUNCH.AND.ZHEAT(JL)>2.0_JPRB) THEN
        ZTKE = ABS(PTENOGW_U(JL,JK)*PU(JL,JK)+PTENOGW_V(JL,JK)*PV(JL,JK))*MAX(0.0_JPRB,ZHEAT(JL))
        ZTKE=MIN(ZTKE,1.E5_JPRB)
        ! derive EDR from TKE (based on Lilly (1966) Eq. 47)
        ZCGWD(JL,JK)= (ZC_EPS*ZTKE)**0.33333_JPRB
        ZCGWD(JL,JK)= ZEA(4)*ZCGWD(JL,JK)**ZB(4)
        PEDRP(JL,JK,1)=PEDRP(JL,JK,1)+ZCGWD(JL,JK)
        PEDRP(JL,JK,2)=PEDRP(JL,JK,2)+ZCGWD(JL,JK)
        PEDRP(JL,JK,1)=MIN(PEDRP(JL,JK,1),1.0_JPRB)
        PEDRP(JL,JK,2)=MIN(PEDRP(JL,JK,2),1.0_JPRB)
      ENDIF
    ENDDO
  ENDDO

ENDIF

END ASSOCIATE 
IF (LHOOK) CALL DR_HOOK('DIAG_TURB',1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_TURB
