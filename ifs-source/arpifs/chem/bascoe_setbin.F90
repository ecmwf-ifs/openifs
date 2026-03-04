! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_SETBIN
!**   DESCRIPTION 
!     ----------
!
!   Part of BASCOE / TM5 routines for IFS chemistry: 
!     AUTHOR.
!     -------
!        Coded in C-IFS by VINCENT HUIJNEN    *KNMI*
!        Original code from BASCOE_CTM v4s09, simonc@oma.be, June 2008
!
!     Subroutine SETBIN must be called once before any call to
!     subroutines PSCBOX, or INITSP.
!
!     This subroutine calculates the radius, surface area, and volume of
!     each particle size bin on a geometrically increasing volume scale.
!     Further the lower and upper volume bin borders are calculated, and
!     the bin radius-width in um.
!     These values are stored in array PTSIZE to be used througout the
!     PSC box model, and posssibly also in the calling program unit. The
!     values in PTSIZE must not be changed by the calling program unit.
!     The content of array PTSIZE is as follows:
!         Radius (m):                            PTSIZE(I,1)
!         Particle surface area (m**2):          PTSIZE(I,2)
!         Particle volume (m**3):                PTSIZE(I,3)
!         Lower bin border (m**3):               PTSIZE(I,4)
!         Upper bin border (m**3):               PTSIZE(I,5)
!         Bin radius-width (microns):            PTSIZE(I,6)
!         ln(r):                                 PTSIZE(I,7)
!         (ln(r))**2:                            PTSIZE(I,8)
!     In array RADIUSCLASS is stored limiting radii in cumulated size distributions.
!
!
!
!     Input:
!     INTEGER NBINS
!     REAL(KIND=8)  RMIN,RMAX
!
!     Output:
!     REAL PTSIZE(NBINS,8)
!
!     Input:
!         NBINS:                  Number of radii groups
!         RMIN:    (microns)      Minimum particle radius
!         RMAX:    (microns)      Minimum particle radius  in initial state
!
!     Output:
!         PTSIZE:  (m,m**2,m**3)  Particle radii, surface, volume, etc.
!-----------------------------------------------------------------------

USE BASCOE_MODULE      , ONLY :  NBINS, RMIN, RMAX, PTSIZE, D1
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST             , ONLY : RPI


IMPLICIT NONE
!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
REAL(KIND=JPRB)                :: ZVRATIO
INTEGER(KIND=JPIM)             :: JK
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_SETBIN',0,ZHOOK_HANDLE )

ZVRATIO=EXP(LOG(RMAX/RMIN)/(REAL(NBINS-1)/3.0))
D1=((2.0/(ZVRATIO+1.0))**(1.0/3.0))*(ZVRATIO**(1.0/3.0)-1.0)
!VH LOGVR=LOG(ZVRATIO)
DO JK=1,NBINS
!         Radius (m):
          PTSIZE(JK,1)=RMIN*1.0E-6*ZVRATIO**((JK-1.0)/3.0)
!         Particle surface area (m**2):
          PTSIZE(JK,2)=4.0*RPI*PTSIZE(JK,1)**2
!         Particle volume (m**3):
          PTSIZE(JK,3)=4.0*RPI*(PTSIZE(JK,1)**3)/3.0
!         Lower bin border (m**3):
          PTSIZE(JK,4)=PTSIZE(JK,3)*2.0/(ZVRATIO+1.0)
!         Upper bin border (m**3):
          PTSIZE(JK,5)=PTSIZE(JK,4)*ZVRATIO
!         Bin radius-width (microns):
          PTSIZE(JK,6)=PTSIZE(JK,1)*D1*1.0E6
!         ln(r):
          PTSIZE(JK,7)=LOG(PTSIZE(JK,1))
!         (ln(r))**2:
          PTSIZE(JK,8)=PTSIZE(JK,7)**2
ENDDO
!----------------------------------------------------------------------------
!      VR=VRATIO
!      ZERO=NDMIN
      D1=D1/SQRT(2.0*RPI)
!      SEDI=.TRUE.
!----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('BASCOE_SETBIN',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_SETBIN
