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

SUBROUTINE GPGEO_SACC(KPROMA,KSTART,KPROF,KFLEV,PHI,PHIF,PT,PU,PR,PLNPR,PALPH,PVGEOM)

!**** *GPGEO_SACC* - Computes half and full level geopotential height "gz".
!                    Shallow-atmosphere complete-Coriolis formulation of Tort & Dubos (2013)
!     Purpose.
!     --------

!      Computes half and full level geopotential height "gz".

!      Laplace relation writes:

!       d (gz)/d prehyd =  RT/pre (g/(2*omega*U-g))= (RT/prehyd) * (prehyd/pre)*((g/(2*omega*U-g))

!      where:
!       - "gz" is the geopotential height.
!       - "prehyd" is the hydrostatic pressure.
!       - "pre" is the total pressure including non-hydrostatic effects.
!       - "R" is the air constant (including moisture effects).
!       - "T" is the temperature.
!       - "U" is zonal wind (u*cos phi)
!       - "omega" is Coriolis frequency
!       - "g" is gravitational acceleration

!      It is important to note that the geopotential height here is not the usual hydrostatic "gz",
!      but the one consistent with the vertical momentum equation in Tort & Dubos (2013). Hence the extra "U"-term 
!      at the surface the geop. height in both cases is the and defined by Phi_s = g z[surf].

!      Integrating the Laplace equations yields the following discretisation
!      for "gz".

!      * "gz" at interlayer "lbar":

!        g z[lbar] = g z[surf]
!        + sum[k=L to l+1] (prehyd/pre)[k] R[k] T[k] delta[k] *((g/(2*omega*U[k]-g))

!      * "gz" at layer "l":

!        g z[l] = g z[lbar] + (prehyd/pre)[l] R[l] T[l] alpha[l]*((g/(2*omega*U[k]-g))

!**   Interface.
!     ----------
!        *CALL* *GPGEO_SACC(...)

!        Explicit arguments :
!        --------------------
!          KPROMA : horizontal dimensioning                          (input)
!          KSTART : start of work                                    (input)
!          KPROF  : depth of work                                    (input)
!          KFLEV  : number of levels                                 (input)
!          PHI    : geopotential height "gz" at interlayers          (output)
!          PHIF   : geopotential height "gz" at layers               (output)
!          PT     : temperature at layers                            (input)
!          PU     : zonal wind at layers                             (input)
!          PR     : "R" at layers for hydrostatic model              (input)
!                    "(prehyd/pre) R" at layers for NH model
!          PLNPR  : term "delta" on layers                           (input)
!                   (= logarithm of ratio of pressure if "ndlnpr=0")
!          PALPH  : term "alpha" on layers                           (input)
!          PVGEOM : vertical geometry from the model                 (input)

!        Implicit arguments :    None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!     Inna Polichtchouk  *ECMWF*: based on GPGEO, but adapted for LSACC model

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVERT  , ONLY : TVERTICAL_GEOM
USE YOMCST   , ONLY : RG, ROMEGA

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PHI(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PHIF(KPROMA,KFLEV) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PT(KPROMA,KFLEV) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PU(KPROMA,KFLEV)
REAL(KIND=JPRB)     ,INTENT(IN)    :: PR(KPROMA,KFLEV) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PLNPR(KPROMA,KFLEV) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PALPH(KPROMA,KFLEV) 
TYPE(TVERTICAL_GEOM),INTENT(IN)    :: PVGEOM

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZPHI(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB) :: ZOUT(KPROMA,0:KFLEV)

INTEGER(KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB) :: ZTWOOM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGEO_SACC',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    COMPUTES HALF AND FULL LEVEL GEOPOTENTIAL HEIGHT.
!              -------------------------------------------------
ZTWOOM=2.0_JPRB*ROMEGA
IF(PVGEOM%YRCVER%LVERTFE) THEN
  DO JLEV=1,KFLEV
    DO JLON=KSTART,KPROF
      ZPHI(JLON,JLEV)=PR(JLON,JLEV)*PT(JLON,JLEV)*RG/(ZTWOOM*PU(JLON,JLEV)-RG)*&
       & PLNPR(JLON,JLEV)*PVGEOM%YRVETA%VFE_RDETAH(JLEV)  
    ENDDO
  ENDDO
  ZPHI(KSTART:KPROF,0)=0.0_JPRB
  ZPHI(KSTART:KPROF,KFLEV+1)=0.0_JPRB
  CALL VERDISINT(PVGEOM%YRVFE,PVGEOM%YRCVER,'IBOT','11',KPROMA,KSTART,KPROF,KFLEV,ZPHI,ZOUT)

  DO JLEV=1,KFLEV
    DO JLON=KSTART,KPROF
      PHIF(JLON,JLEV) =ZOUT(JLON,JLEV-1)+PHI(JLON,KFLEV)
    ENDDO
  ENDDO
  DO JLEV=KFLEV,1,-1
    DO JLON=KSTART,KPROF
      PHI(JLON,JLEV-1)=PHI(JLON,JLEV)&
       & -PR(JLON,JLEV)*PT(JLON,JLEV)*RG/(ZTWOOM*PU(JLON,JLEV)-RG)*PLNPR(JLON,JLEV)  
    ENDDO
  ENDDO
ELSE
  DO JLEV=KFLEV,1,-1
    DO JLON=KSTART,KPROF
      PHI(JLON,JLEV-1)=PHI(JLON,JLEV)&
       & -PR(JLON,JLEV)*PT(JLON,JLEV)*RG/(ZTWOOM*PU(JLON,JLEV)-RG)*PLNPR(JLON,JLEV)  

      PHIF(JLON,JLEV) =PHI(JLON,JLEV)&
       & -PALPH(JLON,JLEV)*PR(JLON,JLEV)*PT(JLON,JLEV)*RG/(ZTWOOM*PU(JLON,JLEV)-RG)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGEO_SACC',1,ZHOOK_HANDLE)
END SUBROUTINE GPGEO_SACC
