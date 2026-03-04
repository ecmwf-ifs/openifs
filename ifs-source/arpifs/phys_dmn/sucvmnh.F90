! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUCVMNH(YDML_PHY_MF,KULOUT)

! Purpose:
! --------
! Setup of KFB convection

! Interface:
! ----------
! KULOUT: logical unit for the output

! Externals:
! ----------
! None

! Method:
! -------

! Reference:
! ----------

! Author:
! ------
! 04-02-05 G. Hello

! Modifications:
! --------------
! End Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        E. BAZILE     03-August-2006 Read NAMCVMNH
!        E. BAZILE     21-Feb.-2008 Correction for XWTRIG
!----------------------------------------------------------------
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULNAM
IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT), TARGET :: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

REAL(KIND=JPRB) , POINTER ::  XCDEPTH_D
REAL(KIND=JPRB) , POINTER ::  XATPERT
REAL(KIND=JPRB) , POINTER ::  OTADJD
REAL(KIND=JPRB) , POINTER ::  XTFRZ1
REAL(KIND=JPRB) , POINTER ::  OTADJS
REAL(KIND=JPRB) , POINTER ::  XZLCL
REAL(KIND=JPRB) , POINTER ::  XCRAD
REAL(KIND=JPRB) , POINTER ::  XSTABT
LOGICAL         , POINTER ::  LSETTADJ, LSMOOTH
REAL(KIND=JPRB) , POINTER ::  XZPBL
REAL(KIND=JPRB) , POINTER ::  XWTRIG
REAL(KIND=JPRB) , POINTER ::  XTFRZ2
REAL(KIND=JPRB) , POINTER ::  XBTPERT
REAL(KIND=JPRB) , POINTER ::  XDTPERT
REAL(KIND=JPRB) , POINTER ::  XNHGAM
REAL(KIND=JPRB) , POINTER ::  XENTR
REAL(KIND=JPRB) , POINTER ::  XCDEPTH
REAL(KIND=JPRB) , POINTER ::  XBW
REAL(KIND=JPRB) , POINTER ::  XSTABC
REAL(KIND=JPRB) , POINTER ::  XA25
REAL(KIND=JPRB) , POINTER ::  XAW

#include "namcvmnh.nam.h"

! 1. Set default values
! ---------------------

IF (LHOOK) CALL DR_HOOK('SUCVMNH',0,ZHOOK_HANDLE)
ASSOCIATE(LKFBD=>YDML_PHY_MF%YRARPHY%LKFBD, LKFBS=>YDML_PHY_MF%YRARPHY%LKFBS, &
 & LDOWN=>YDML_PHY_MF%YRCVMNH%LDOWN, NIICE=>YDML_PHY_MF%YRCVMNH%NIICE, LDEEP=>YDML_PHY_MF%YRCVMNH%LDEEP, &
 & NSETENS=>YDML_PHY_MF%YRCVMNH%NSETENS, LDIAGCONV=>YDML_PHY_MF%YRCVMNH%LDIAGCONV, &
 & LREFRESH_ALL=>YDML_PHY_MF%YRCVMNH%LREFRESH_ALL, LSHALLOW=>YDML_PHY_MF%YRCVMNH%LSHALLOW, &
 & LCVPPKF=>YDML_PHY_MF%YRPHY%LCVPPKF, YDCVMNH=>YDML_PHY_MF%YRCVMNH)
XNHGAM => YDCVMNH%XNHGAM
XA25 => YDCVMNH%XA25
XTFRZ1 => YDCVMNH%XTFRZ1
XTFRZ2 => YDCVMNH%XTFRZ2
XENTR => YDCVMNH%XENTR
OTADJD => YDCVMNH%OTADJD
LSMOOTH => YDCVMNH%LSMOOTH
XZLCL => YDCVMNH%XZLCL
OTADJS => YDCVMNH%OTADJS
XBW => YDCVMNH%XBW
XDTPERT => YDCVMNH%XDTPERT
XCRAD => YDCVMNH%XCRAD
XBTPERT => YDCVMNH%XBTPERT
XAW => YDCVMNH%XAW
XCDEPTH => YDCVMNH%XCDEPTH
XZPBL => YDCVMNH%XZPBL
XATPERT => YDCVMNH%XATPERT
XSTABT => YDCVMNH%XSTABT
LSETTADJ => YDCVMNH%LSETTADJ
XCDEPTH_D => YDCVMNH%XCDEPTH_D
XWTRIG => YDCVMNH%XWTRIG
XSTABC => YDCVMNH%XSTABC

LDEEP=LKFBD
LSHALLOW=LKFBS
LDIAGCONV=.TRUE.
LSETTADJ=.TRUE.
LREFRESH_ALL=.TRUE.
LDOWN=.TRUE.
LSMOOTH=.TRUE.

OTADJD=3600._JPRB
OTADJS=3600._JPRB

XA25     = 625.E6_JPRB    ! 25 km x 25 km reference grid area
XCRAD       = 50._JPRB    ! cloud radius
XCDEPTH     = 0.5E3_JPRB  ! minimum necessary shallow cloud depth
XCDEPTH_D   = 2.5E3_JPRB  ! maximum allowed shallow cloud depth
XDTPERT     = .2_JPRB     ! add small Temp perturbation at LCL
XATPERT     = 0._JPRB     ! 0.=original scheme , recommended = 1000. 
XBTPERT     = 1._JPRB     ! 1.=original scheme , recommended = 0.
XENTR    = 0.02_JPRB      ! entrainment constant (m/Pa) = 0.2 (m)
XZLCL    = 0.5E3_JPRB   ! maximum allowed allowed height
                     ! difference between the DPL and the surface
XZPBL    = 40.E2_JPRB     ! minimum mixed layer depth to sustain convection
XNHGAM   = 1.3333_JPRB    ! accounts for non-hydrost. pressure
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XWTRIG   = 6._JPRB   ! constant in vertical velocity trigger
XTFRZ1   = 268.16_JPRB    ! begin of freezing interval
XTFRZ2   = 248.16_JPRB    ! end of freezing interval
XSTABT   = 0.75_JPRB      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95_JPRB      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
XAW      = 0._JPRB   ! 0.= Original scheme , 1 = recommended 
XBW      = 1._JPRB   ! 1.= Original scheme , 0 = recommended



!*       2.    Modify default values.
!              ----------------------

CALL POSNAM(NULNAM,'NAMCVMNH')
READ(NULNAM,NAMCVMNH)

NSETENS=0.0_JPRB
NIICE=1.0_JPRB

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMCVMNH '')')

WRITE(UNIT=KULOUT,FMT='('' LDEEP= '',L5&
 & ,'' LSHALLOW= '',L5&
 & ,'' LDIAGCONV= '',L5&
 & ,'' LSETTADJ= '',L5&
 & ,'' LREFRESH_ALL= '',L5&
 & ,'' LDOWN= '',L5&
 & ,'' LSMOOTH= '',L5&
 & )')&
 & LDEEP,LSHALLOW,LDIAGCONV,LSETTADJ,LREFRESH_ALL,LDOWN,LSMOOTH

WRITE(UNIT=KULOUT,FMT='('' NSETENS= '',I2&
 & ,'' NIICE= '',I2&
 & )')&
 & NSETENS,NIICE  
IF (LCVPPKF) THEN
   WRITE(UNIT=KULOUT,FMT='('' COMMON YOMCVMNH PARTIE SHALLOW KFB'')')
   WRITE(UNIT=KULOUT,FMT='('' OTADJD= '',E14.7&
     & ,'' OTADJS= '',E14.7, '' XA25= '',E14.7,'' XCRAD= '',E14.7&
     & ,'' XCDEPTH= '',E14.7, '' XCDEPTH_D= '',E14.7,'' XDTPERT= '',E14.7&
     & ,'' XATPERT= '',E14.7, '' XBTPERT= '',E14.7&
     & ,'' XENTR= '',E14.7, '' XZLCL= '',E14.7,'' XZPBL= '',E14.7&
     & ,'' XWTRIG= '',E14.7&
     & ,'' XNHGAM= '',E14.7, '' XTFRZ1= '',E14.7,'' XTFRZ2= '',E14.7&
     & ,'' XSTABT= '',E14.7, '' XSTABC= '',E14.7&
     & ,'' XAW= '',E14.7, '' XBW= '',E14.7&
     & )')&
     & OTADJD,OTADJS  ,&
     & XA25, XCRAD, XCDEPTH, XCDEPTH_D, XDTPERT, XATPERT, XBTPERT,XENTR,&
     & XZLCL, XZPBL, XWTRIG, XNHGAM, XTFRZ1, XTFRZ2, XSTABT, XSTABC,&
     & XAW, XBW
ENDIF
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCVMNH',1,ZHOOK_HANDLE)

END SUBROUTINE SUCVMNH
