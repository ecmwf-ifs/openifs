MODULE SUVEXC_MOD
CONTAINS
SUBROUTINE SUVEXC(LD_LEOCWA,LD_LEOCCO,LD_LEOCSA,LD_LEOCLA,&
          & LD_LWCOU, LD_LWCOU2W, LD_LWCOUHMF,&
          & LD_LSCMEC,LD_LROUGH,PEXTZ0M,PEXTZ0H,PRPLRG,&
          & YDEXC)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_EXC  , ONLY : TEXC
! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *SUVEXC* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOS_EXC*

!     A.C.M. BELJAARS         E.C.M.W.F.       2/11/89

!     PURPOSE
!     -------

!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOS_EXC*

!     INTERFACE.
!     ----------

!    Logicals (In):

!      LD_LEOCWA : .T. if WARM OCEAN LAYER PARAMETRIZATION active
!      LD_LEOCCO : .T. if COOL OCEAN SKIN PARAMETRIZATION active
!      LD_LEOCSA : .T. if SALINTY EFFECT ON SATURATION AT OCEAN SURFACE active
!      LD_LEOCLA : .T. if LANGMUIR CIURCULATION EFFECT IN VOSKIN active

!      LD_LWCOU    : .T. COUPLING TO WAVE MODEL
!      LD_LWCOU2W  : .T. COUPLING TO WAVE MODEL WITH FEEDBACK TO ATMOSPHERE
!      LD_LWCOUHMF : .T. SEA STATE DEPENDENT HEAT AND MOISTURE FLUXES IF COUPLED TO WAVE MODEL 



!     CALL *SUVEXC* FROM *SUSURF*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     NONE.

!     REFERENCE.
!     ----------

!     MODIFICATIONS
!     -------------
!     J.-J. MORCRETTE         E.C.M.W.F.      91/07/14
!     M.Hamrud                01-Oct-2003     CY28 Cleaning
!     P. Viterbo              09/06/2005      Externalise surf
!     N.Semane+P.Bechtold 04-10-2012 Add PRPLRG factor for small planet
!     E. Dutra                10/10/2014      add LELWTL
!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL           ,INTENT(IN)    :: LD_LEOCWA
LOGICAL           ,INTENT(IN)    :: LD_LEOCCO
LOGICAL           ,INTENT(IN)    :: LD_LEOCSA
LOGICAL           ,INTENT(IN)    :: LD_LEOCLA
LOGICAL           ,INTENT(IN)    :: LD_LWCOU 
LOGICAL           ,INTENT(IN)    :: LD_LWCOU2W
LOGICAL           ,INTENT(IN)    :: LD_LWCOUHMF
LOGICAL           ,INTENT(IN)    :: LD_LSCMEC
LOGICAL           ,INTENT(IN)    :: LD_LROUGH
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTZ0M
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTZ0H
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPLRG
TYPE(TEXC)        ,INTENT(INOUT) :: YDEXC

REAL(KIND=JPRB) :: Z_CEPZ0O
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUVEXC_MOD:SUVEXC',0,ZHOOK_HANDLE)
ASSOCIATE(LELWDD=>YDEXC%LELWDD, LELWTL=>YDEXC%LELWTL, LEOCCO=>YDEXC%LEOCCO, &
 & LEOCLA=>YDEXC%LEOCLA, LEOCSA=>YDEXC%LEOCSA, LEOCWA=>YDEXC%LEOCWA, &
 & LWCOU=>YDEXC%LWCOU, LWCOU2W=>YDEXC%LWCOU2W, LWCOUHMF=>YDEXC%LWCOUHMF, &
 & LROUGH=>YDEXC%LROUGH, LSCMEC=>YDEXC%LSCMEC, RCHAR=>YDEXC%RCHAR, &
 & REPDU2=>YDEXC%REPDU2, REPUST=>YDEXC%REPUST, REXTZ0H=>YDEXC%REXTZ0H, &
 & REXTZ0M=>YDEXC%REXTZ0M, RKAP=>YDEXC%RKAP, RNU=>YDEXC%RNU, RNUH=>YDEXC%RNUH, RNUM=>YDEXC%RNUM, &
 & RNUQ=>YDEXC%RNUQ, RPARZI=>YDEXC%RPARZI, RZ0ICE=>YDEXC%RZ0ICE)

! First set of constants
RKAP   =0.4_JPRB
RZ0ICE =0.001_JPRB/PRPLRG

! Simplified physics constant
RCHAR  =0.0155_JPRB

! Other constants

REPDU2 =(0.1_JPRB)**2
REPUST=0.0001_JPRB  

!     KINEMATIC VISCOSITY OF AIR
RNU   =1.5E-5_JPRB/PRPLRG
RNUM  =0.11_JPRB*RNU
RNUH  =0.40_JPRB*RNU
RNUQ  =0.62_JPRB*RNU

!     ENTRAINMENT PARAMETRIZATION

RPARZI=1000._JPRB/PRPLRG

!     COMPUTATION OF SKIN TEMPERATURE

LELWDD=.TRUE.

!     COMPUTATION OF Tiled net longwave 

LELWTL=.TRUE.   ! if true use the Tiled net longwave in the surface scheme

!     Ocean skin effects and ocean salinity effect on esat

LEOCWA=LD_LEOCWA
LEOCCO=LD_LEOCCO
LEOCSA=LD_LEOCSA
LEOCLA=LD_LEOCLA

!     Coupling to wave model

LWCOU=LD_LWCOU
LWCOU2W=LD_LWCOU2W
LWCOUHMF=LD_LWCOUHMF

!     SCM parameters for surface roughness length

LSCMEC =LD_LSCMEC
LROUGH =LD_LROUGH
REXTZ0M=PEXTZ0M
REXTZ0H=PEXTZ0H


!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVEXC_MOD:SUVEXC',1,ZHOOK_HANDLE)
END SUBROUTINE SUVEXC
END MODULE SUVEXC_MOD

