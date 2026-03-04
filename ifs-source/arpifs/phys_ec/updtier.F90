! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE UPDTIER(YDML_PHY_RAD,YDDYNA,YDRIP,KSTEP,KRADFR,PTSTEP)

!**** *UPDTIER* - UPDATE TIME FOR ECMWF FULL RADIATION COMPUTATIONS

!     Purpose.
!     --------
!     UPDATE TIME OF THE MODEL + 1/2 TIME BETWEEN 2 FULL RADIATION STEPS

!**   Interface.
!     ----------
!        *CALL* *UPDTIER(...)      from *ECRADFR*

!        Explicit arguments :
!        --------------------
!        KSTEP : TIME STEP INDEX
!        KRADFR: FREQUENCY OF FULL RADIATION COMPUTATIONS
!        PTSTEP: TIME STEP LENGTH

!        Implicit arguments :
!        --------------------
!        YOERIP

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
!      Jean-Jacques Morcrette             *ECMWF*
!      Original : 93-02-22

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      P Bechtold 18/05/2012   Introduce invariant RDAYI for ZWSOVRM
!                              otherwise SW wrong for small planet
!      N.Semane+P.Bechtold 04-10-2012 replace RDAYI by RDAY consistently with utility/updtim.F90
!      K. Yessad (July 2014): Move some variables.
!      R Hogan       Dec 2014  Correct time by half a timestep
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOMDYNA    , ONLY : TDYNA
USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMRIP0   , ONLY : NSSSSS, RTIMST
USE YOMRIP    , ONLY : TRIP
USE YOMCST    , ONLY : RPI, RDAY, REA, REPSM, RI0, RV, RCPV, RETV, RCW, RCS, RTT, &
 &                     RLVTT, RLSTT, RALPW, RBETW, RGAMW, RALPS, RBETS, RGAMS, RALPD, RBETD, RGAMD  
USE YOMDYNCORE, ONLY : LAPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_PHYSICS_RADIATION_TYPE),INTENT(INOUT):: YDML_PHY_RAD
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KRADFR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISTADD, ISTASS, ITIME, IZT, ISEC

REAL(KIND=JPRB) :: ZANGOZC, ZCOTHOZ, ZDEASOM, ZDECLIM, ZEQTIMM,&
 & ZHGMT, ZI0, ZSITHOZ, ZSOVRM, ZSTATI, ZTETA, ZTHETOZ, &
 & ZTIMTR, ZWSOVRM  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "suecozo.intfb.h"

#include "fctast.func.h"
#include "fcttim.func.h"
#include "fcttrm.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UPDTIER',0,ZHOOK_HANDLE)
ASSOCIATE(LCENTREDTIMESZA=>YDML_PHY_RAD%YRERAD%LCENTREDTIMESZA, LPERPET=>YDML_PHY_RAD%YRERAD%LPERPET, &
 & NHINCSOL=>YDML_PHY_RAD%YRERAD%NHINCSOL, &
 & RSOLINC=>YDML_PHY_RAD%YRERDI%RSOLINC, &
 & RCODECM=>YDML_PHY_RAD%YRERIP%RCODECM, RSIDECM=>YDML_PHY_RAD%YRERIP%RSIDECM, RCOVSRM=>YDML_PHY_RAD%YRERIP%RCOVSRM, &
 & RSIVSRM=>YDML_PHY_RAD%YRERIP%RSIVSRM, &
 & RTIMTR=>YDRIP%RTIMTR)
!     ------------------------------------------------------------------

!          1.  TIME EQUATION
!              -------------

ITIME=NINT(PTSTEP)
IF (YDDYNA%LTWOTL) THEN
  ! In the two-level timestepping scheme, the solar zenith angle
  ! should be computed at a time half-way between the current and
  ! future timestep, hence the 0.5 here.  IZT is the number of seconds
  ! since the start of the forecast.
  IZT=NINT(PTSTEP*(REAL(KSTEP,JPRB)+0.5_JPRB))
ELSE
  IZT=ITIME*KSTEP
ENDIF
IZT=ITIME*KSTEP

!--
IF (LPERPET) THEN
  ISEC=IZT/NINT(RDAY)
  IZT=IZT-ISEC*NINT(RDAY)
ENDIF
!--


IF (LCENTREDTIMESZA) THEN
  ! IZT is the number of seconds into the forecast and is already half
  ! a model timestep ahead.  For radiation every timestep (KRADFR=1),
  ! we don't want to modify this, hence the -1 below.
  ZSTATI=REAL(IZT,JPRB)+0.5_JPRB*(KRADFR-1)*ITIME
ELSE
  ! The older scheme adds half a radiation timestep to a time that
  ! is already half a model timestep ahead.
  ZSTATI=REAL(IZT,JPRB)+0.5_JPRB*KRADFR*ITIME
ENDIF
ISTADD=IZT/NINT(RDAY)
ISTASS=MOD(IZT,NINT(RDAY))
ZTIMTR=RTIMST+ZSTATI
ZHGMT=REAL(MOD(NINT(ZSTATI)+NSSSSS,NINT(RDAY)),JPRB)

ZTETA=RTETA(ZTIMTR)
IF( LAPE ) THEN
  ZDEASOM=RRSAQUA(ZTETA)
  ZDECLIM=RDSAQUA(ZTETA)
  ZEQTIMM=RETAQUA(ZTETA)
ELSE
  ZDEASOM=RRS(ZTETA)
  ZDECLIM=RDS(ZTETA)
  ZEQTIMM=RET(ZTETA)
ENDIF
ZSOVRM =ZEQTIMM+ZHGMT
ZWSOVRM=ZSOVRM*2.0_JPRB*RPI/RDAY

IF (NHINCSOL /= 0) THEN
  ZI0=RSOLINC
ELSE
  ZI0=RI0
ENDIF

RCODECM=COS(ZDECLIM)
RSIDECM=SIN(ZDECLIM)

RCOVSRM=COS(ZWSOVRM)
RSIVSRM=SIN(ZWSOVRM)
ZTHETOZ=RTETA(RTIMTR)
ZANGOZC=REL(ZTHETOZ)-1.7535_JPRB
ZCOTHOZ=COS(ZANGOZC)
ZSITHOZ=SIN(ZANGOZC)

!     ------------------------------------------------------------------

!          2.  COMPUTES TIME-DEPENDENT PARAMETERS FOR OZONE CLIMATOLOGY
!              --------------------------------------------------------

CALL SUECOZO (ZANGOZC,YDRIP%YREOZOC)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDTIER',1,ZHOOK_HANDLE)
END SUBROUTINE UPDTIER
