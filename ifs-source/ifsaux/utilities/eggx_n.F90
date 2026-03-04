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

SUBROUTINE EGGX_N(PI,PRA,KROTEQ,PLONR,PLATR,PBETA,PLON1,PLAT1,PLON2,PLAT2,&
     & PLON0,PLAT0,PRPK,KULOUT,KSOTRP,KGIVO,&
     & PGELAM,PGELAT,PGM,PGNORX,PGNORY,KDLSA,&
     & KDLSUR,KDGSA,KDGEN,KDLUN,KDLUX,KDGUN,KDGUX,&
     & PDELX,PDELY,PLONC,PLATC)  

! Version 2006.1016 by JD GRIL

!** *EGGX_N*  - the interface to both old and new geographic package of ALADIN

!     Purpose.
!     --------
!      To provide an interface to both new and old geographic setup
!      routines MAKDO and EGGX.
!      Convert between the old EGGX domain definition and the new 
!      domain definition
!      The old definition uses corners, number of grid points and EGGX
!      projection definition parameters
!      The new definition uses the centre of domain, number of grid points
!      and the resolution in x and y.

!**   Interface.
!     ----------
!        *CALL* *EGGX_N

!     Explicit arguments :
!     --------------------

!     INPUT:
!      PI : pi (3.14ETC)
!      PRA  : radius of spherical planet
!      KROTEQ : previous rotation parameter
!               here it is a control of the direction of the conversion
!               since the options KROTEQ>0 are already no more supported
!      PLONR : geographic longitude of reference point of rotation
!      PLATR : geographic latitude of reference point of rotation
!      PLON0 : longitude of reference for the projection
!      PLAT0 : latitude of reference for the projection
!      PBETA : angle (in rd) between x-axis and rotated latitude circles
!              at the reference longitude
!              (usually, pbeta = 0. : gives pure projections)
!      KSOTRP : isotropy parameter under projection
!      KGIVO  : choice of reference point for projection
!      KDLSA:KDLSUR : lower and upper first dimensions of arrays (X)
!      KDGSA:KDGEN  : lower and upper second dimensions of arrays (Y)
!      KDLUN:KDLUX  : lower and upper first dimensions of
!                     the domain of interest, where arrays are initialized.
!      KDGUN:KDGUX  : lower and upper second dimensions of
!                     the domain of interest, where arrays are initialized.
!      KULOUT : unit of control prints file

!     INPUT/OUTPUT (depending on KROTEQ):
!      PLON1, PLAT1 : latitude of the south-west corner of useful domain
!      PLON2, PLAT2 : latitude of the north-east corner of useful domain
!      PLONC, PLATC : longitude and latitude of the centre of domain
!      PDELX, PDELY : horizontal resolution in x and y direction
!      PRPK  : projection parameter and definition in the old EGGX
!              PRPK = 10. projection type self determined
!                         by minimizing the variation of the map factor
!              PRPK = 1.  polar stereographic projection
!              0. < PRPK < 1.  lambert conformal projection with
!                              cone parameter prpk
!              PRPK = 0.  mercator conformal projection
!              PRPK < 0.  no projection
!             on output, PRPK contains the effective projection
!             parameter that has been used

!     OUTPUT:
!      PGELAM, PGELAT : longitude and latitude of the grid points
!      PGM            : map factor at the grid points
!      PGNORX, PGNORY : components of the vector pointing to the north pole
!                       at the grid points locations

!     Implicit arguments :
!     --------------------

!     Method.
!     -------
!      The parameter KROTEQ controls the direction in which the conversion
!      is performed:
!        KROTEQ<0: the new parameter set defining the domain and projection
!                  (PLON0,PLAT0,PLONC,PLATC,PDELX,PDELY) is converted
!                  to the old one by the call of MAKDO
!                     KROTEQ = -1 Normal mode
!                     KROTEQ = -2 Mercator Rotated-Tilted mode
!        KROTEQ=0: the old parameter set defining the domain and projection
!                  (PLONR,PLATR,PBETA,PLON1,PLAT1,PLON2,PLAT2,PLON0,PLAT0,
!                  PRPK,KSOTRP,KGIVO) is converted to the new one by the
!                  call of EGGX
!      All geographic coordinates must be in radians.
!      All latitudes in <-PI/2;PI/2>, all longitudes <-PI;+PI>.
!********* future : All latitudes in <-PI/2;PI/2>, all longitudes <0;2*PI>

!     Externals.
!     ----------
!      EGGX : old geographic setup routine
!      MAKDO: new geographic setup routine

!     Reference.
!     ----------

!     Author.
!     -------
!       Jean-Daniel GRIL, 2000-2001

!     Modifications.
!     --------------
!       Modified in April 2001 by M.Janousek:
!           all input/output angles are in radians
!           add check of unsupported old EGGX domains
!        C. Fischer & J.D. Gril 02-04-16 : Improve new EGGX calls
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        O.Spaniel     Oct-2004 cleaning AL29
!        JD Gril       17-Nov-2004 Mercator Rotated-Tilted case 
!        JD Gril       18-Nov-2005 add KGIVO=0 and KGIVO=0
!        JD Gril       03-Jui-2006 comment lines 328/329
!        JD Gril       15-Sep-2006 correct both previous case
!        JD Gril       16-Oct-2006 cleaning
!      F. Vana  05-Mar-2015  Support for single precision
!        P. Marguinaud 04-10-2016 Port to single precision
!        R. El Khatib 22-Mar-2017 disable verbosity if KULOUT < 0
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE EGGPACK   ,ONLY : LOLA,XY,NBPTS,PGN,DELTA,ERROR,DOMI,PARAM_PROJ,MAKDO,&
 & REF_DATAS,LATLON_TO_XY,XY_TO_LATLON
USE EGGANGLES ,ONLY : ANGLE_DOMAIN

!     ------------------------------------------------------------------

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(INOUT) :: KROTEQ
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT
INTEGER(KIND=JPIM),INTENT(INOUT) :: KSOTRP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KGIVO
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLSA
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLSUR
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUN
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUN
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRA 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLONR 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLATR 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBETA 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLON1 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT1 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLON2 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT2 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLON0 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT0 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRPK 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGELAM(KDLSA:KDLSUR,KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGELAT(KDLSA:KDLSUR,KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGM(KDLSA:KDLSUR,KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGNORX(KDLSA:KDLSUR,KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGNORY(KDLSA:KDLSUR,KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDELX 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDELY 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLONC 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLATC 

!     ------------------------------------------------------------------

TYPE (LOLA)              :: YL_TLKRES, YL_TLCENT, YL_TLSW_LOLA, YL_TLNE_LOLA
TYPE (LOLA), ALLOCATABLE :: YL_TLGRID_LOLA(:,:)
TYPE (XY)                :: YL_TLSW_XY , YL_TLNE_XY, YL_TLCENT_XY
TYPE (NBPTS)             :: YL_TLNB_PTS
TYPE (PGN), ALLOCATABLE  :: YL_TLGRID_PGN(:,:)
TYPE (DELTA)             :: YL_TLDEL
TYPE (ERROR)             :: YL_TLERR
TYPE (DOMI)              :: YL_TLGRID_INFO
TYPE (PARAM_PROJ)        :: YL_TLMODDOM
REAL(KIND=JPRB)          :: ZGRID_MF(KDLUX-KDLUN+1,KDGUX-KDGUN+1) 
REAL(KIND=JPRB)          :: ZRTD
REAL(KIND=JPRB)          :: ZPI, ZRA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB)          :: ZEPS
LOGICAL                  :: LLMRT

!     ------------------------------------------------------------------

#include "eggx.h"

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EGGX_N',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

! The routine can be sometimes called before constants are initialized
! Check if it is the case and then set defaults
IF (INT(PI*100._JPRB) == 314) THEN
  ZPI=REAL(PI,JPRB)
  ZRA=REAL(PRA,JPRB)
ELSE
  ZPI=ASIN(1.0_JPRB)*2.0_JPRB
  ZRA=6371229._JPRB
ENDIF
ZEPS=EPSILON(1.0_JPRB)*100.0_JPRB
ZRTD   = 180.0_JPRB/ZPI
PGELAM = 0.0_JPRB
PGELAT = 0.0_JPRB
PGM    = 0.0_JPRB
PGNORX = 0.0_JPRB
PGNORY = 0.0_JPRB

IF (KULOUT >= 0) WRITE(KULOUT,*) '********* INFO of Input data in EGGX_N **************'
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLON0 (rd) = ',PLON0,'PLON0 (dg) = ',PLON0*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLAT0 (rd) = ',PLAT0,'PLAT0 (dg) = ',PLAT0*ZRTD 
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLONC (rd) = ',PLONC,'PLONC (dg) = ',PLONC*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLATC (rd) = ',PLATC,'PLATC (dg) = ',PLATC*ZRTD 
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLON1 (rd) = ',PLON1,'PLON1 (dg) = ',PLON1*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLAT1 (rd) = ',PLAT1,'PLAT1 (dg) = ',PLAT1*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLON2 (rd) = ',PLON2,'PLON2 (dg) = ',PLON2*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLAT2 (rd) = ',PLAT2,'PLAT2 (dg) = ',PLAT2*ZRTD 
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PDELX      = ',PDELX
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PDELY      = ',PDELY
IF (KULOUT >= 0) WRITE(KULOUT,*) 'KROTEQ     = ',KROTEQ
IF (KULOUT >= 0) WRITE(KULOUT,*) '****************************************************'

IF (KROTEQ < 0) THEN
  IF (KULOUT >= 0) WRITE(KULOUT,*) 'KROTEQ < 0 : New Eggx domain'
  ! the input parameters are in the new style of the domain definition
  LLMRT = (KROTEQ == -2)
  IF (.NOT.LLMRT) KROTEQ = -1
  IF (KULOUT >= 0) WRITE(KULOUT,*) 'KROTEQ = ',KROTEQ,'LLMRT = ',LLMRT
  IF (LLMRT .AND. (ABS(PLAT0) >= ZEPS)) THEN
    WRITE(ABS(KULOUT),*) 'EGGX_N: PLAT0=',PLAT0,&
     & ' MUST BE EQUAL ZERO IF LLMRT IS TRUE!'
    CALL ABOR1('EGGX_N: LLMRT & PLAT0 INCONSISTENT')
  ENDIF
  KSOTRP              = 0_JPIM
  KGIVO               = 0_JPIM
  PLONR               = 0.0_JPRB
  PLATR               = 0.0_JPRB
  PBETA               = 0.0_JPRB
  YL_TLKRES%LON       = PLON0*ZRTD
  YL_TLKRES%LAT       = PLAT0*ZRTD
  YL_TLCENT%LON       = PLONC*ZRTD
  YL_TLCENT%LAT       = PLATC*ZRTD
  YL_TLKRES           = ANGLE_DOMAIN(YL_TLKRES,REAL(ZPI,JPRB),'-+','D')
  YL_TLCENT           = ANGLE_DOMAIN(YL_TLCENT,REAL(ZPI,JPRB),'-+','D')
  YL_TLDEL%ONX        = PDELX
  YL_TLDEL%ONY        = PDELY
  YL_TLNB_PTS%ONX     = KDLUX-KDLUN+1
  YL_TLNB_PTS%ONY     = KDGUX-KDGUN+1
  ALLOCATE(YL_TLGRID_LOLA(KDLUX-KDLUN+1,KDGUX-KDGUN+1))
  ALLOCATE(YL_TLGRID_PGN(KDLUX-KDLUN+1,KDGUX-KDGUN+1))
  CALL MAKDO(YL_TLKRES,YL_TLCENT,YL_TLDEL,YL_TLNB_PTS,YL_TLGRID_LOLA,&
   & ZGRID_MF,YL_TLGRID_PGN,YL_TLGRID_INFO,YL_TLERR,.TRUE.,.TRUE.,&
   & REAL(ZPI,JPRB),REAL(ZRA,JPRB),KULOUT,LLMRT)  
  PLON1                           = YL_TLGRID_LOLA(1,1)%LON
  PLAT1                           = YL_TLGRID_LOLA(1,1)%LAT
  PLON2                           = YL_TLGRID_LOLA(YL_TLNB_PTS%ONX,YL_TLNB_PTS%ONY)%LON
  PLAT2                           = YL_TLGRID_LOLA(YL_TLNB_PTS%ONX,YL_TLNB_PTS%ONY)%LAT
  PLON0                           = YL_TLKRES%LON
  PLAT0                           = YL_TLKRES%LAT 
  PLONC                           = YL_TLCENT%LON
  PLATC                           = YL_TLCENT%LAT
  PRPK                            = YL_TLGRID_INFO%INFO_PROJ%KL
  PGELAM(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_LOLA(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%LON  
  PGELAT(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_LOLA(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%LAT  
  PGM(KDLUN:KDLUX,KDGUN:KDGUX)    = ZGRID_MF(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)
  PGNORX(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_PGN(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%ONX  
  PGNORY(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_PGN(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%ONY  
  DEALLOCATE(YL_TLGRID_LOLA)
  DEALLOCATE(YL_TLGRID_PGN)
ELSE
  ! KROTEQ>0 => the input is in the old style
  ! Some old EGGX domains are no more supported in ALADIN
  ! Check if that is not the case of this domain
  IF (KROTEQ > 0) THEN
    WRITE(ABS(KULOUT),*) 'EGGX_N: NROTEQ=',KROTEQ,&
     & ' IS NOT VALID VALUE, IT MUST BE ZERO!'  
    CALL ABOR1('EGGX_N: UNSUPPORTED NROTEQ')
  ELSEIF (PBETA /= 0.0_JPRB .AND. PRPK == 0.0_JPRB) THEN
    WRITE(ABS(KULOUT),*) 'EGGX_N: ROTATED DOMAIN IN MERCATOR PROJECTION NOT&
     & SUPPORTED (EBETA HAS TO BE 0)'
    CALL ABOR1('EGGX_N: UNSUPPORTED EBETA')
  ELSEIF ( ABS(PRPK-ABS(SIN(PLAT0))) > 1.E-7 ) THEN
    WRITE(ABS(KULOUT),*) 'EGGX_N: YOU SEEM TO HAVE A SECANT CASE OF PROJECTION'
    WRITE(ABS(KULOUT),*) '       ERPK=',PRPK,'  SIN(ELAT0)=',SIN(PLAT0)
    CALL ABOR1('EGGX_N: UNSUPPORTED SECANT PROJECTION')
  ENDIF
  IF (KULOUT >= 0) WRITE(KULOUT,*) 'KROTEQ = 0 : Old Eggx domain'
  ! Call EGGX to handle cases when corners may change
  ! Fill in every case the arrays. Not needed in model
  ! (call echien) but outside. So, either with old eggx,
  ! either with makdo (see below)
  IF(KSOTRP/=0 .OR. KGIVO/=0 .OR. PRPK==10._JPRB) THEN
    IF (KULOUT >= 0) WRITE(KULOUT,*) 'Call old EGGX to handle cases when corners may change'
    IF (KULOUT >= 0) WRITE(KULOUT,*) 'KSOTRP = ',KSOTRP,'  KGIVO = ',KGIVO,'  PRPK = ',PRPK
    CALL EGGX(REAL(ZPI,JPRB),REAL(ZRA,JPRB),KROTEQ,PLONR,PLATR,PBETA,PLON1,PLAT1,PLON2,PLAT2,&
     & PLON0,PLAT0,PRPK,KULOUT,KSOTRP,KGIVO,&
     & PGELAM,PGELAT,PGM,PGNORX,PGNORY,KDLSA,KDLSUR,KDGSA,KDGEN,&
     & KDLUN,KDLUX,KDGUN,KDGUX,PDELX,PDELY)  
  ENDIF
  ! Now calculate x,y coordinates of the corners, compute the centre
  ! point and convert it to lat-lon by the EGGPACK functions
  IF (KULOUT >= 0) WRITE(KULOUT,*) 'COMPUTATION OF CENTER'
  YL_TLKRES%LON    = PLON0*ZRTD
  YL_TLKRES%LAT    = PLAT0*ZRTD
  YL_TLKRES        = ANGLE_DOMAIN(YL_TLKRES,ZPI,'-+','D')
  YL_TLMODDOM      = REF_DATAS(YL_TLKRES,ZRA)
  YL_TLSW_LOLA%LON = PLON1
  YL_TLSW_LOLA%LAT = PLAT1
  YL_TLNE_LOLA%LON = PLON2
  YL_TLNE_LOLA%LAT = PLAT2
  YL_TLSW_LOLA     = ANGLE_DOMAIN(YL_TLSW_LOLA,ZPI,'-+','R')
  YL_TLNE_LOLA     = ANGLE_DOMAIN(YL_TLNE_LOLA,ZPI,'-+','R')
  YL_TLSW_XY       = LATLON_TO_XY(YL_TLSW_LOLA,YL_TLMODDOM,ZPI)
  YL_TLNE_XY       = LATLON_TO_XY(YL_TLNE_LOLA,YL_TLMODDOM,ZPI)
  YL_TLCENT_XY%X   = (YL_TLSW_XY%X+YL_TLNE_XY%X)*0.5_JPRB
  YL_TLCENT_XY%Y   = (YL_TLSW_XY%Y+YL_TLNE_XY%Y)*0.5_JPRB
  YL_TLCENT        = ANGLE_DOMAIN(XY_TO_LATLON(YL_TLCENT_XY,YL_TLMODDOM,ZPI),ZPI,'0+','R')
  PLONC            = YL_TLCENT%LON
  PLATC            = YL_TLCENT%LAT 
  ! If KSOTRP=0 and KGIVO=0 then the values of SW,NE,REF are fixed. They come from
  ! - old eggx (rare) but can be computed by new eggx
  ! - new eggx but with NCADFORM=0 (old "cadre")
  ! In both cases we can use new eggx to recompute missing values, this way protects
  ! us from old eggx possible bugs in case number 2 (new eggx + NCADFORM=0) 
  ! We use Makdo to compute all arrays
  IF(KSOTRP==0 .AND. KGIVO==0 .AND. PRPK/=10._JPRB) THEN
    ! Computation of resolution to use Makdo
    ! We test the case "point" or "linear" wide
    ! Protect from divided by zero
    IF (KULOUT >= 0) WRITE(KULOUT,*) 'KSOTRP==0 .AND. KGIVO==0 .AND. PRPK/=10'
    IF (KULOUT >= 0) WRITE(KULOUT,*) 'COMPUTATION OF RESOLUTION AND USE OF MAKDO'
    IF (KULOUT >= 0) WRITE(KULOUT,*) 'because cadre is in old style but domain may be created'
    IF (KULOUT >= 0) WRITE(KULOUT,*) 'by new eggx (may be not supported by old eggx)'  
    IF ((KDLUX-KDLUN) == 0) THEN
      PDELX = 0.0_JPRB
    ELSE
      PDELX = ABS(YL_TLNE_XY%X-YL_TLSW_XY%X)/REAL(KDLUX-KDLUN,JPRB)
    ENDIF
    IF ((KDGUX-KDGUN) == 0) THEN
      PDELY = 0.0_JPRB
    ELSE
      PDELY = ABS(YL_TLNE_XY%Y-YL_TLSW_XY%Y)/REAL(KDGUX-KDGUN,JPRB)
    ENDIF
    YL_TLCENT%LON       = PLONC*ZRTD
    YL_TLCENT%LAT       = PLATC*ZRTD
    YL_TLKRES           = ANGLE_DOMAIN(YL_TLKRES,ZPI,'-+','D')
    YL_TLCENT           = ANGLE_DOMAIN(YL_TLCENT,ZPI,'-+','D')
    YL_TLDEL%ONX        = PDELX
    YL_TLDEL%ONY        = PDELY
    YL_TLNB_PTS%ONX     = KDLUX-KDLUN+1
    YL_TLNB_PTS%ONY     = KDGUX-KDGUN+1
    ALLOCATE(YL_TLGRID_LOLA(KDLUX-KDLUN+1,KDGUX-KDGUN+1))
    ALLOCATE(YL_TLGRID_PGN(KDLUX-KDLUN+1,KDGUX-KDGUN+1))
    CALL MAKDO(YL_TLKRES,YL_TLCENT,YL_TLDEL,YL_TLNB_PTS,YL_TLGRID_LOLA,&
     & ZGRID_MF,YL_TLGRID_PGN,YL_TLGRID_INFO,YL_TLERR,.TRUE.,.TRUE.,&
     & REAL(ZPI,JPRB),REAL(ZRA,JPRB),KULOUT,.FALSE.)  
    PLON1                           = YL_TLGRID_LOLA(1,1)%LON
    PLAT1                           = YL_TLGRID_LOLA(1,1)%LAT
    PLON2                           = YL_TLGRID_LOLA(YL_TLNB_PTS%ONX,YL_TLNB_PTS%ONY)%LON
    PLAT2                           = YL_TLGRID_LOLA(YL_TLNB_PTS%ONX,YL_TLNB_PTS%ONY)%LAT
    PLON0                           = YL_TLKRES%LON
    PLAT0                           = YL_TLKRES%LAT 
    PLONC                           = YL_TLCENT%LON
    PLATC                           = YL_TLCENT%LAT
    PRPK                            = YL_TLGRID_INFO%INFO_PROJ%KL
    PGELAM(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_LOLA(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%LON  
    PGELAT(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_LOLA(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%LAT  
    PGM(KDLUN:KDLUX,KDGUN:KDGUX)    = ZGRID_MF(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)
    PGNORX(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_PGN(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%ONX  
    PGNORY(KDLUN:KDLUX,KDGUN:KDGUX) = YL_TLGRID_PGN(1:YL_TLNB_PTS%ONX,1:YL_TLNB_PTS%ONY)%ONY  
    DEALLOCATE(YL_TLGRID_LOLA)
    DEALLOCATE(YL_TLGRID_PGN)
  ENDIF
  IF (KULOUT >= 0) WRITE(KULOUT,*) 'SWX = ',YL_TLSW_XY%X,'NEX = ',YL_TLNE_XY%X,'CEX = ',YL_TLCENT_XY%X
  IF (KULOUT >= 0) WRITE(KULOUT,*) 'SWY = ',YL_TLSW_XY%Y,'NEY = ',YL_TLNE_XY%Y,'CEY = ',YL_TLCENT_XY%Y   
ENDIF

IF (KULOUT >= 0) WRITE(KULOUT,*) '********* INFO before Return out of EGGX_N *********'
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLON0 (rd) = ',PLON0,'PLON0 (dg) = ',PLON0*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLAT0 (rd) = ',PLAT0,'PLAT0 (dg) = ',PLAT0*ZRTD 
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLONC (rd) = ',PLONC,'PLONC (dg) = ',PLONC*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLATC (rd) = ',PLATC,'PLATC (dg) = ',PLATC*ZRTD 
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLON1 (rd) = ',PLON1,'PLON1 (dg) = ',PLON1*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLAT1 (rd) = ',PLAT1,'PLAT1 (dg) = ',PLAT1*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAM(KDLUN,KDGUN):SW (rd) = ',PGELAM(KDLUN,KDGUN)
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAM(KDLUN,KDGUN):SW (dg) = ',PGELAM(KDLUN,KDGUN)*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAT(KDLUN,KDGUN):SW (rd) = ',PGELAT(KDLUN,KDGUN)
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAT(KDLUN,KDGUN):SW (dg) = ',PGELAT(KDLUN,KDGUN)*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLON2 (rd) = ',PLON2,'PLON2 (dg) = ',PLON2*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PLAT2 (rd) = ',PLAT2,'PLAT2 (dg) = ',PLAT2*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAM(KDLUX,KDGUX):NE (rd) = ',PGELAM(KDLUX,KDGUX)
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAM(KDLUX,KDGUX):NE (dg) = ',PGELAM(KDLUX,KDGUX)*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAT(KDLUX,KDGUX):NE (rd) = ',PGELAT(KDLUX,KDGUX)
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGELAT(KDLUX,KDGUX):NE (dg) = ',PGELAT(KDLUX,KDGUX)*ZRTD
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PRPK = ',PRPK
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGM(KDLUN,KDGUN)    (SW) = ',PGM(KDLUN,KDGUN)
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGNORX(KDLUN,KDGUN) (SW) = ',PGNORX(KDLUN,KDGUN)
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PGNORY(KDLUN,KDGUN) (SW) = ',PGNORY(KDLUN,KDGUN)
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PDELX = ',PDELX
IF (KULOUT >= 0) WRITE(KULOUT,*) 'PDELY = ',PDELY
IF (KULOUT >= 0) WRITE(KULOUT,*) '****************************************************'

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EGGX_N',1,ZHOOK_HANDLE)
END SUBROUTINE EGGX_N
