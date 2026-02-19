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

SUBROUTINE SUPPVI

!**** *SUPPVI*  - Initialize variables used in the vertical interpolator
!                 (for example FULL-POS vertical interpolator or observation
!                 vertical interpolator).

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *SUPPVI

!        EXPLICIT ARGUMENTS
!        ------------------

!        IMPLICIT ARGUMENTS
!        ------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        K. Yessad (Aug 2009)

!     MODIFICATIONS.
!     --------------
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El Khatib 13-Sep-2016 Interoperability EC GRIB2 vs FA
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT,NULNAM  
USE YOMPPVI  , ONLY : LESCALE  ,LESCALE_T, LESCALE_Q, LESCALE_U, LESCALE_PD, &
  & LESCALE_GFL, LRPPUV_CSTEXT, LRPPUV_CALLITPQ, LPPVIVX, RPPVIVX, RPPVIVP,  &
  & LNOTS_T
USE YOMCT0   , ONLY : LOLDPP, LECMWF, LARPEGEF

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "namppvi.nam.h"

#include "posnam.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUPPVI',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1. SET DEFAULT VALUES
!           ------------------

LESCALE=.NOT.LECMWF.OR.LARPEGEF
LESCALE_T=LESCALE
LESCALE_Q=LESCALE
LESCALE_U=LESCALE
LESCALE_PD=LESCALE
LESCALE_GFL=LESCALE
LNOTS_T=.FALSE.
LRPPUV_CSTEXT=.FALSE.
LRPPUV_CALLITPQ=LOLDPP

LPPVIVX=.FALSE.
RPPVIVX=120.0_JPRB
RPPVIVP=1000.0_JPRB

!*       2. READ NAMELIST
!           -------------

CALL POSNAM(NULNAM,'NAMPPVI')
READ(NULNAM,NAMPPVI)


!*       3. PRINT OUT FINAL VALUES
!           ----------------------

WRITE(UNIT=NULOUT,FMT='('' MODULE YOMPPVI'')')
WRITE(UNIT=NULOUT,FMT='('' LESCALE = '',L2,'' LESCALE_T = '',L2,&
 & '' LESCALE_Q = '',L2,'' LESCALE_U = '',L2,'' LESCALE_PD = '',L2,&
 & '' LESCALE_GFL = '',L2,'' LRPPUV_CSTEXT = '',L2,&
 & '' LRPPUV_CALLITPQ = '',L2,'' LNOTS_T = '',L2)')&
 & LESCALE,LESCALE_T,LESCALE_Q,LESCALE_U,LESCALE_PD,&
 & LESCALE_GFL,LRPPUV_CSTEXT,LRPPUV_CALLITPQ,LNOTS_T
WRITE(UNIT=NULOUT,FMT='('' LPPVIVX = '',L2 &
 & ,'' RPPVIVX = '',F9.2,'' RPPVIVP = '',F9.2)')&
 & LPPVIVX ,RPPVIVX ,RPPVIVP

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPPVI',1,ZHOOK_HANDLE)
END SUBROUTINE SUPPVI
