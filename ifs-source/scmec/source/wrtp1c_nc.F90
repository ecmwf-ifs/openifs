! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WRTP1C_NC(YDGEOMETRY,YDMODEL,YDSURF)

!**** *WRTP1C_NC*  - Write prognostic variables of the one-column model

!     Purpose.
!     --------
!     Write out prognostic variables in NetCDF format

!**   Interface.
!     ----------
!        *CALL* *WRTP1C_NC

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation 
!        of the single column model

!     Author.
!     -------
!        Martin Koehler  *ECMWF*

!     Modifications.
!     --------------
!        Original      94-01-11
!        J.Teixeira   Jan.-1995  new output files.
!        M.Koehler    Sep.-2000  converted to netCDF output.
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!        G. Carver    Aug/2012   Fixes for gfortran compiler
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE TYPE_MODEL  , ONLY : MODEL
USE YOMHOOK     , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE PARDIM1C
USE PARDIM
USE YOMPHYDS
USE YOETHF      , ONLY : RHOH2O
USE YOMCT3      , ONLY : NSTEP
USE YOMGT1C0    , ONLY : UT0      ,VT0      ,TT0      ,QT0      ,&
                        &WT0      ,ST0      ,AT0      ,SPT0     ,&
                        &RNT0     ,SNT0
USE YOMGP1C0    , ONLY : TSA0     ,WSA0     ,SNS0     ,TL0      ,&
                        &WL0      ,RSN0     ,ASN0     ,TSN0
USE YOMGPD1C    , ONLY : VQSAT
USE YOMLOG1C    , ONLY : NPOSPRG  ,NPOSASC
USE YOETHF
USE YOMCST
USE INTDYN_MOD  , ONLY : YYTXYB
USE SURFACE_FIELDS_MIX, ONLY : TSURF


IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
TYPE(TSURF),    INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM) :: ICSS

INTEGER(KIND=JPIM) :: JSLEV, IST, IEND, JK, JALEV

REAL(KIND=JPRB)    :: ZSNNU0, ZRENU0, ZSCALE, ZQS

REAL(KIND=JPRB)    :: ZLINU0(YDSURF%YSP_SBD%NLEVS), ZRH(YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB)    :: ZPRESH (0:YDGEOMETRY%YRDIMV%NFLEVG)   ! HALF LEVEL PRESSURE
REAL(KIND=JPRB)    :: ZPRESF (  YDGEOMETRY%YRDIMV%NFLEVG)   ! FULL LEVEL PRESSURE
REAL(KIND=JPRB)    :: ZXYB9  (1,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)

REAL(KIND=JPRB)    :: ZDUM   (  YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB)    :: ZR0    (  YDGEOMETRY%YRDIMV%NFLEVG)   ! R
REAL(KIND=JPRB)    :: ZCP0   (  YDGEOMETRY%YRDIMV%NFLEVG)   ! CP
REAL(KIND=JPRB)    :: ZKAP   (  YDGEOMETRY%YRDIMV%NFLEVG)   ! K=R/CP

REAL(KIND=JPRB)    :: ZPHI0  (0:YDGEOMETRY%YRDIMV%NFLEVG)   ! HALF LEVEL GEOPOTENTIAL
REAL(KIND=JPRB)    :: ZPHIF0 (  YDGEOMETRY%YRDIMV%NFLEVG)   ! FULL LEVEL GEOPOTENTIAL

REAL(KIND=JPRB)    :: ZRDAW(YDSURF%YSP_SBD%NLEVS)! SOIL LEVEL THICKNESS

REAL(KIND=JPRB)    :: ZTHETA(YDGEOMETRY%YRDIMV%NFLEVG)    ,ZTHETA_E(YDGEOMETRY%YRDIMV%NFLEVG) &
                    & ,ZDRY_ST(YDGEOMETRY%YRDIMV%NFLEVG) ,&
                    & ZMOIST_ST(YDGEOMETRY%YRDIMV%NFLEVG) ,CPD_MOIST(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: ISTATUS, IDIMID, IDIMLEN, INCID, IVARID
INTEGER(KIND=JPIM) :: ISTART1, ICOUNT1, ISTART2(2), ICOUNT2(2), ICOUNT3(2), ICOUNT4(2)
INTEGER(KIND=JPIM) :: ITIMESECS, IDAY, IHOUR, IMIN, ISEC

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "netcdf.inc"
#include "fcttre.func.h"

!     ------------------------------------------------------------------
#include "surf_inq.h"
#include "gpgeo.intfb.h"
#include "handle_err_nc.intfb.h"
#include "gphpre.intfb.h"
#include "varwrite1c_nc.intfb.h"
#include "gprcp.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRTP1C_NC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM,YDMP=>YDGEOMETRY%YRMP, &
 & YDVAB=>YDGEOMETRY%YRVAB,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, LEFLAKE=>YDEPHY%LEFLAKE, YSURF=>YDEPHY%YSURF, &
 & TSTEP=>YDMODEL%YRML_GCONF%YRRIP%TSTEP, &
 & YSD_VF=>YDSURF%YSD_VF, YSP_SBD=>YDSURF%YSP_SBD, YSP_SL=>YDSURF%YSP_SL,&
 & SD_VF=>YDSURF%SD_VF, SP_SL=>YDSURF%SP_SL)
ICSS=YSP_SBD%NLEVS

!        1.    SCALING SOME VARIABLES AND NEW VARIABLES.
!              -----------------------------------------

ZSNNU0 = SNS0(1) / RHOH2O
ZRENU0 = WL0 (1) / RHOH2O
CALL SURF_INQ(YSURF,PRDAW=ZRDAW)
ZLINU0(1:ICSS) = WSA0(1:ICSS) / RHOH2O / ZRDAW(1:ICSS)

IST  = 1
IEND = 1

ZPRESH(NFLEVG) = EXP(SPT0)


!        1.1   HEIGHT AND PRESSURE ON HALF AND FULL LEVELS.

!   computation of pressure at half and full model levels.
CALL GPHPRE(NPROMA,NFLEVG,IST,IEND,YDVAB,YDGEOMETRY%YRCVER,ZPRESH,PXYB=ZXYB9,PRESF=ZPRESF)

!   computation of r, cp and kappa.
CALL GPRCP(NPROMA,IST,IEND,NFLEVG,PQ=QT0,PQI=ST0,PQL=WT0,PQR=RNT0,PQS=SNT0,PCP=ZCP0,PR=ZR0,PKAP=ZKAP)

!   computation of hydrostatic equation.
ZPHI0(NFLEVG) = YDGEOMETRY%YROROG(1)%OROG(1)*RG
CALL GPGEO(NPROMA,IST,IEND,NFLEVG,ZPHI0,ZPHIF0,TT0,ZR0, &
         & ZXYB9(1,:,YYTXYB%M_LNPR),ZXYB9(1,:,YYTXYB%M_ALPH),YDGEOMETRY%YRVERT_GEOM)

!        1.2   computation of rh.

DO JK=1,NFLEVG
  ZQS=FOEEWM(TT0(JK))/ZPRESF(JK)
  ZQS=MIN(0.5_JPRB,ZQS)
  ZQS=ZQS/(1.0_JPRB-RETV*ZQS)
  ZRH(JK)=QT0(JK)/ZQS
ENDDO

!        1.3   thermodynamic conserved variables.
!              ...caution:  cpd_moist might be slightly misused

CPD_MOIST(1:NFLEVG) = RCPD * ( 1.0 + RVTMP2 * QT0(1:NFLEVG) )

! potential temperature
ZTHETA(1:NFLEVG)    = TT0(1:NFLEVG) * ( 1.0E5 / ZPRESF(1:NFLEVG) ) ** ( ZR0(1:NFLEVG) / CPD_MOIST(1:NFLEVG) )

! equivalent potential temperature (Emanuel, 1994, p120)
ZTHETA_E(1:NFLEVG)  = TT0(1:NFLEVG) * ( 1.0E5 / ZPRESF(1:NFLEVG) ) ** ( RD  / CPD_MOIST(1:NFLEVG) )   &  
     &   *  ZRH(1:NFLEVG) ** ( - QT0(1:NFLEVG) * RV / CPD_MOIST(1:NFLEVG) )                           &
     &   *  EXP( RLVTT * QT0(1:NFLEVG) / TT0(1:NFLEVG) / CPD_MOIST(1:NFLEVG) )

! dry static energy
ZDRY_ST(1:NFLEVG)   = ZPHIF0(1:NFLEVG) + TT0(1:NFLEVG) * CPD_MOIST(1:NFLEVG)

! moist static energy
ZMOIST_ST(1:NFLEVG) = ZDRY_ST(1:NFLEVG) + RLVTT * QT0(1:NFLEVG)


!        2.    WRITE TO DIAGNOSTIC TEXT FILE.
!              ------------------------------

!        2.1   WRITE TIME STEP.

WRITE(NPOSASC,*) 'TIME STEP'
WRITE(NPOSASC,'(I6)') NSTEP
ITIMESECS = NSTEP*TSTEP
IDAY  = INT(ITIMESECS/86400.)
IHOUR = INT(MOD(ITIMESECS,86400)/3600.)
IMIN  = INT(MOD(ITIMESECS,3600)/60.)
ISEC  = INT(MOD(ITIMESECS,60))
WRITE(*,'(A8,1X,I6,3X,A3,1X,I3,3X,A4,1X,I0.2,A1,I0.2,A1,I0.2)') &
     &  'Timestep',NSTEP,'Day',IDAY,'Time',IHOUR,':',IMIN,':',ISEC

WRITE(NPOSASC,*) 'PROGNOSTIC VARIABLES'

!*       2.3   WRITE ATMOSPHERIC VARIABLES.

WRITE(NPOSASC,3019) '  P  ','  U  ','  V  ','  T  ','  Q  '&
                 &,'  A  ','  L  ','  I  ','  RH  '
DO JALEV=1,NFLEVG
  WRITE(NPOSASC,3020) ZPRESF(JALEV),UT0(JALEV),VT0(JALEV)&
   &,TT0(JALEV),QT0(JALEV),AT0(JALEV),WT0(JALEV),ST0(JALEV),ZRH(JALEV)
ENDDO
3019 FORMAT(9(2X,A12))
3020 FORMAT(9(2X,E12.6))

WRITE(NPOSASC,*) 'SURFACE PRESSURE  -  LOG SURFACE PRESSURE'
WRITE(NPOSASC,'(2(2x,E12.6))') ZPRESH(NFLEVG),SPT0

!*       2.4   WRITE SOIL VARIABLES.

WRITE(NPOSASC,*) ' SOIL TEMPERATURE - SOIL MOISTURE '
DO JSLEV=1,ICSS
  WRITE(NPOSASC,4020) TSA0(JSLEV) , ZLINU0(JSLEV)
ENDDO
4020 FORMAT(2(2X,E12.6))

!*       2.5   WRITE SKIN VARIABLES.

WRITE(NPOSASC,*) ' SKIN TEMP. - SKIN. RES. CONT. '
WRITE(NPOSASC,5020) TL0 , ZRENU0
5020 FORMAT(2(2X,E12.6))

!*       2.6   WRITE SNOW DEPTH AND OTHER SNOW VARIABLES.

WRITE(NPOSASC,*) 'SNOW: DEPTH     TEMPERATURE   ALBEDO        DENSITY'
WRITE(NPOSASC,'(4(2x,E12.6))') ZSNNU0,TSN0,ASN0,RSN0

!*       2.7  WRITE LAKE VARIABLES.
IF (LEFLAKE) THEN
WRITE(NPOSASC,*) 'LAKE: DEPTH     MX-LAY TEMP   MX-LAY DPTH   BOTTOM TEMP'
WRITE(NPOSASC,'(4(2x,E12.6))') SD_VF(:,YSD_VF%YDL%MP,1),SP_SL(:,YSP_SL%YLMLT%MP9,1), &
 & SP_SL(:,YSP_SL%YLMLD%MP9,1), SP_SL(:,YSP_SL%YLBLT%MP9,1)
WRITE(NPOSASC,*) 'LAKE: TOT TEMP  SHAPE FACT    ICE TEMP      ICE DEPTH     L. COVER'   
WRITE(NPOSASC,'(5(2x,E12.6))') SP_SL(:,YSP_SL%YLTLT%MP9,1),SP_SL(:,YSP_SL%YLSHF%MP9,1), &
 & SP_SL(:,YSP_SL%YLICT%MP9,1), SP_SL(:,YSP_SL%YLICD%MP9,1),SD_VF(:,YSD_VF%YCLK%MP,1)

ENDIF


!        3.    WRITE netCDF OUTPUT.
!              --------------------

INCID = NPOSPRG         ! netCDF file unit number (not = fortran unit #)

!        3.1   set-up

ISTATUS = NF_INQ_DIMID   (INCID, 'time', IDIMID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_INQ_DIMLEN  (INCID, IDIMID, IDIMLEN)
CALL HANDLE_ERR_NC(ISTATUS)

ISTART1    = IDIMLEN+1  ! 1-d variables: starting index
ICOUNT1    = 1          !      -"-       written indices
ISTART2(1) = 1          ! 2-d variables - dim 1: starting index
ICOUNT2(1) = NFLEVG     !      -"-               written indices
ISTART2(2) = IDIMLEN+1  ! 2-d variables - dim 2: starting index
ICOUNT2(2) = 1          !      -"-               written indices
ICOUNT3(1) = NFLEVG+1   ! half level variables
ICOUNT3(2) = 1
ICOUNT4(1) = ICSS       ! soil/sea ice variables
ICOUNT4(2) = 1

!        3.2   write time & time step

CALL VARWRITE1C_NC (INCID, 1, (/ISTART1/),(/ICOUNT1/), 'time', (/NSTEP*TSTEP/) )

!istatus = NF_INQ_VARID   (incid, 'timestp', ivarid)
!call handle_err_nc(istatus)
!istatus = NF_PUT_VAR1_INT(incid,ivarid,idimlen+1, nstep)     
!call handle_err_nc(istatus)

!        3.3   write atmospheric variables.

!        pressure (full levels - T,q,u,v)
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'pressure_f',       ZPRESF)
!        pressure (half levels - w, fluxes)
CALL VARWRITE1C_NC (INCID, NFLEVG+1, ISTART2,ICOUNT3, 'pressure_h',       ZPRESH)
!        pressure (full levels - T,q,u,v)
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'height_f',         ZPHIF0/RG)
!        pressure (half levels - w, fluxes)
CALL VARWRITE1C_NC (INCID, NFLEVG+1, ISTART2,ICOUNT3, 'height_h',         ZPHI0/RG)
!        temperature
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   't',                TT0)
!        u-wind
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'u',                UT0)
!        v-wind
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'v',                VT0)
!        water vapor mixing ratio
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'q',                QT0)
!        relative humidity
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'relative_humidity',ZRH)
!        cloud fraction
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'cloud_fraction',   AT0)
!        liquid water mixing ratio
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'ql',               WT0 )
!        ice water mixing ratio
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'qi',               ST0 )
!        rain water mixing ratio
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'qr',               RNT0 )
!        snow water mixing ratio
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'qsn',               SNT0 )
!        potential temperature
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'pot_temperature',  ZTHETA)
!        equivalent potential temperature
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'pot_temp_e',       ZTHETA_E)
!        dry static energy
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'dry_st_energy',    ZDRY_ST)
!        moist static energy
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'moist_st_energy',  ZMOIST_ST)
!        saturation mixing ratio
CALL VARWRITE1C_NC (INCID, NFLEVG, ISTART2,ICOUNT2,   'q_sat',            VQSAT(1:NFLEVG))

!        3.4   write land/ocean/sea-ice variables.

!        snow depth (equivalent liquid water depth)
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'snow',         (/ZSNNU0/) )
!        snow temperature
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 't_snow',       TSN0)
!        snow albedo
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'albedo_snow',  ASN0)
!        snow density
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'density_snow', RSN0)
!        skin temperature
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 't_skin',       TL0)
!        skin reservoir water content
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'q_skin',       (/ZRENU0/) )
!        soil temperature
CALL VARWRITE1C_NC (INCID, ICSS,  ISTART2,ICOUNT4,   't_soil',           TSA0)
!        soil moisture
CALL VARWRITE1C_NC (INCID, ICSS,  ISTART2,ICOUNT4,   'q_soil',           ZLINU0)

!        3.5   write lake variables.

IF (LEFLAKE) THEN
!        Lake depth
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'dl',   SD_VF(:,YSD_VF%YDL%MP,1))
!        Lake mix-layer temperature
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'lmlt', SP_SL(:,YSP_SL%YLMLT%MP9,1))
!        Lake mix-layer depth
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'lmld', SP_SL(:,YSP_SL%YLMLD%MP9,1))
!        Lake bottom temperature
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'lblt', SP_SL(:,YSP_SL%YLBLT%MP9,1))
!        Lake total layer temperature
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'ltlt', SP_SL(:,YSP_SL%YLTLT%MP9,1))
!        Lake shape factor
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'lshf', SP_SL(:,YSP_SL%YLSHF%MP9,1))
!        Lake ice temperature
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'lict', SP_SL(:,YSP_SL%YLICT%MP9,1))
!        Lake ice depth
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'licd', SP_SL(:,YSP_SL%YLICD%MP9,1))
!        Lake cover
CALL VARWRITE1C_NC (INCID, 1,     (/ISTART1/),(/ICOUNT1/), 'cl',   SD_VF(:,YSD_VF%YCLK%MP,1))
ENDIF

!     --------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRTP1C_NC',1,ZHOOK_HANDLE)
END SUBROUTINE WRTP1C_NC
