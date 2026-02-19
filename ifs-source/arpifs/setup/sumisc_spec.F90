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

SUBROUTINE SUMISC_SPEC(YDGEOMETRY,KTESTCASE,PTEMP,PDIV,PVOR,PLNSP,POROG)

!     Purpose.
!     --------
!           Initialize spectral fields for idealized test cases

!        Explicit arguments :
!        --------------------


!     Author.
!     -------
!        S.Malardel *ECMWF*

!     Modifications.
!     --------------

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULERR
USE YOMCST    , ONLY : RPI, RA, RD, RG, ROMEGA, RCPD, RATM
USE YOMDYNCORE, ONLY : RDELTA_T, RHS_KAPPA, RPRESSURE_SCALE, NOISEVOR, &
 &                     RP00_DYN, RU00_DYN, RST0_DYN
USE INTDYN_MOD, ONLY : YYTXYB
IMPLICIT NONE

TYPE(GEOMETRY),     INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN)    :: KTESTCASE
REAL(KIND=JPRB) ,   INTENT(INOUT) :: PTEMP(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,   INTENT(INOUT) :: PDIV(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,   INTENT(INOUT) :: PVOR(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ,   INTENT(INOUT) :: PLNSP(YDGEOMETRY%YRDIM%NSPEC2),POROG(YDGEOMETRY%YRDIM%NSPEC2)

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JLEV, JWORD, JK

!HS
REAL(KIND=JPRB) :: ZVP00, ZU0, ZST
REAL(KIND=JPRB) :: ZU(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZV(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPSP(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZTVERT(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPRESHX(0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPRS
REAL(KIND=JPRB) :: ZT(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVORG(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZOROGG(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZGELAM(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ::    ZXYB   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) ::    ZPRESH(YDGEOMETRY%YRGEM%NGPTOT,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPRESF (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPHI(YDGEOMETRY%YRGEM%NGPTOT,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPHIF(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZRDGAZ(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)

!Klemp profiles
REAL(KIND=JPRB) ::    ZPREF_K(41)
REAL(KIND=JPRB) ::    ZT_K(41),ZU_K(41),ZV_K(41)

! SQUALL LINE
REAL(KIND=JPRB) ::    ZDEEP, ZTHETAM
REAL(KIND=JPRB) ::    ZRADX,ZRADY,ZLATC,ZLONC,ZA,ZB,ZZ,ZZZ

! noise generator
INTEGER(KIND=JPIM) :: im,ia,ic,iran,itop
!     ------------------------------------------------------------------

#include "abor1.intfb.h"

#include "reespe.intfb.h"
#include "uvspe.intfb.h"
#include "gphpre.intfb.h"
#include "gpgeo.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMISC_SPEC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,YDCVER=>YDGEOMETRY%YRCVER, &
 & YDVAB=>YDGEOMETRY%YRVAB,YDVERT_GEOM=>YDGEOMETRY%YRVERT_GEOM)
ASSOCIATE(NRESOL=>YDDIM%NRESOL, NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & NGPTOT=>YDGEM%NGPTOT)

!==================================================================
! Case 71: Held-Suarez
!==================================================================
  !  HS initial temperature prescribed as if at the pole
  ! this was necessary for comparison with Eulag
IF (KTESTCASE == 71) THEN

  ZVP00 = 1.E5_JPRB
  ZPRESHX(0)=YDVAB%VAH(0)+YDVAB%VBH(0)*ZVP00
  DO JLEV=NFLEVG,1,-1
    ZPRESHX(JLEV-1)=YDVAB%VAH(JLEV-1)+YDVAB%VBH(JLEV-1)*ZVP00
    ZPRESHX(JLEV)=YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*ZVP00
    ZPRS=0.5_JPRB*(ZPRESHX(JLEV)+ZPRESHX(JLEV-1))*RPRESSURE_SCALE
    ZTVERT(JLEV) = MAX ( 200._JPRB,&
     &(315._JPRB- RDELTA_T)* ZPRS**RHS_KAPPA )
  ENDDO

  DO JLEV=1,NFLEVG
    DO JWORD=1,NGPTOT
      ZT(JWORD,JLEV)= ZTVERT(JLEV)
    ENDDO
  ENDDO

! Prescribe surface pressure
! constant surface pressure
  ZVP00 = 1.E5_JPRB
  DO JWORD=1,NGPTOT
    ZPSP(JWORD) = LOG ( ZVP00 )
  ENDDO

!!!! From GP space to SP space
! spectral surface pressure
CALL REESPE(YDGEOMETRY,1,1,PLNSP,ZPSP)
! spectral temperature
CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PTEMP,ZT)
! DIV and OROG in spectral space
PDIV=0._JPRB
POROG=0._JPRB
! VORG in GP space
ZVORG=0._JPRB
IF (NOISEVOR == 1) THEN
    ! noise generator
    ! add some noise in vorticity
    im=86436
    ia=1093
    ic=18254
    iran=1   
    itop=30
    DO JLEV=1,NFLEVG
      DO JWORD=1,NGPTOT
        iran=mod(iran*ia+ic,im)
        ZVORG(JWORD,JLEV) = ZVORG(JWORD,JLEV) + 1.e-5*float(iran)/float(im)
      ENDDO
    ENDDO
ENDIF

  ! spectral
  CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PVOR,ZVORG)

!==================================================================
! Case 91: Squall line
!==================================================================
ELSEIF ((KTESTCASE == 91) .OR. (KTESTCASE == 92)) THEN

ZVP00 = RP00_DYN
ZST   = RST0_DYN
!!!!!!!!!!!!!!!!!!!!!!!!!!
! flat orography
DO JWORD=1,NGPTOT
  ZOROGG(JWORD) = 0.0_JPRB    
ENDDO

! set spectral orography
CALL REESPE(YDGEOMETRY,1,1,POROG,ZOROGG)

!!!!!!!!!!!!!!!!!!!!!!!!!!
! constant ln(surface pressure)
DO JWORD=1,NGPTOT
  ZPSP(JWORD) = LOG (ZVP00)
ENDDO
! spectral surface pressure
CALL REESPE(YDGEOMETRY,1,1,PLNSP,ZPSP)

! Compute pressure functions and pressures at half and full levels
DO JWORD=1,NGPTOT
  ZPRESH(JWORD,NFLEVG) = ZVP00
ENDDO

CALL GPHPRE (NGPTOT,NFLEVG,1,NGPTOT,YDVAB,YDCVER,ZPRESH,PXYB=ZXYB,PRESF=ZPRESF)

!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature profile (like all WKS papers)
! VERTICAL PROFIL FROM THE ORIGINAL 41 levels 
ZPREF_K(1)=  271.82816
ZPREF_K(2)=  539.6529
ZPREF_K(3)=  1282.4026
ZPREF_K(4)=  2319.1091
ZPREF_K(5)=  3633.8188
ZPREF_K(6)=  5210.409
ZPREF_K(7)=  7032.732
ZPREF_K(8)=  9084.624
ZPREF_K(9)=  11349.933
ZPREF_K(10)= 13812.486
ZPREF_K(11)= 16456.13
ZPREF_K(12)= 19264.697
ZPREF_K(13)= 22222.035
ZPREF_K(14)= 25311.957
ZPREF_K(15)= 28518.324
ZPREF_K(16)= 31824.938
ZPREF_K(17)= 35215.67
ZPREF_K(18)= 38674.36
ZPREF_K(19)= 42184.836
ZPREF_K(20)= 45730.86
ZPREF_K(21)= 49296.402
ZPREF_K(22)= 52865.195
ZPREF_K(23)= 56421.18
ZPREF_K(24)= 59947.996
ZPREF_K(25)= 63429.74
ZPREF_K(26)= 66849.99
ZPREF_K(27)= 70192.836
ZPREF_K(28)= 73441.836
ZPREF_K(29)= 76581.15
ZPREF_K(30)= 79594.445
ZPREF_K(31)= 82465.37
ZPREF_K(32)= 85177.92
ZPREF_K(33)= 87716.055
ZPREF_K(34)= 90063.38
ZPREF_K(35)= 92204.2
ZPREF_K(36)= 94121.6
ZPREF_K(37)= 95799.85
ZPREF_K(38)= 97223.336
ZPREF_K(39)= 98373.88
ZPREF_K(40)= 99238.28
ZPREF_K(41)= 99799.18

ZT_K(1)=  233.66
ZT_K(2)=  224.28
ZT_K(3)=  221.57
ZT_K(4)=  220.15
ZT_K(5)=  219.25
ZT_K(6)=  218.62
ZT_K(7)=  218.14
ZT_K(8)=  217.77
ZT_K(9)=  217.47
ZT_K(10)= 217.22
ZT_K(11)= 217.01
ZT_K(12)= 216.83
ZT_K(13)= 221.29
ZT_K(14)= 227.16
ZT_K(15)= 232.63
ZT_K(16)= 237.75
ZT_K(17)= 242.55
ZT_K(18)= 247.07
ZT_K(19)= 251.32
ZT_K(20)= 255.34
ZT_K(21)= 259.13
ZT_K(22)= 262.72
ZT_K(23)= 266.12
ZT_K(24)= 269.33
ZT_K(25)= 272.37
ZT_K(26)= 275.25
ZT_K(27)= 277.96
ZT_K(28)= 280.53
ZT_K(29)= 282.94
ZT_K(30)= 285.21
ZT_K(31)= 287.33
ZT_K(32)= 289.31
ZT_K(33)= 291.13
ZT_K(34)= 292.81
ZT_K(35)= 294.34
ZT_K(36)= 295.70
ZT_K(37)= 296.90
ZT_K(38)= 297.93
ZT_K(39)= 298.77
ZT_K(40)= 299.41
ZT_K(41)= 299.84

! vertical interpolation on the current full levels 

DO JWORD = 1,NGPTOT
DO JLEV = 1, NFLEVG
  DO JK=1,41
    IF (ZPRESF(JWORD,JLEV) < ZPREF_K(JK)) THEN
! JLEV is between JK-1 and JK (if JK/=1 !)
      IF (JK==1) THEN
! level above the first full pressure level of original profil
! extrapolation above
      ZT(JWORD,JLEV) = ZT_K(1)+ (ZT_K(2)-ZT_K(1))&
                     & * LOG(ZPRESF(JWORD,JLEV)/ZPREF_K(1))&
                     & / LOG(ZPREF_K(2)/ZPREF_K(1))
      EXIT
      ELSE
      ZT(JWORD,JLEV) = ZT_K(JK-1)+ (ZT_K(JK)-ZT_K(JK-1))&
                     & * LOG(ZPRESF(JWORD,JLEV)/ZPREF_K(JK-1))&
                     & / LOG(ZPREF_K(JK)/ZPREF_K(JK-1))
      EXIT
      ENDIF
    ELSE 
! level under the last full pressure level of original profil
! extrapolation below
      ZT(JWORD,JLEV) = ZT_K(41)+ (ZT_K(40)-ZT_K(41))&
                     & * LOG(ZPRESF(JWORD,JLEV)/ZPREF_K(41))&
                     & / LOG(ZPREF_K(40)/ZPREF_K(41))
    ENDIF
  ENDDO
ENDDO
ENDDO

!!!! Cold (density current) anomaly to trigger the squall line
! anomaly maxi at the surface (dtheta=-8K) and linearly decreases to 0 at 2.5~km
!!! Geopotential
  ZRDGAZ = RD
  ZPHI=0._JPRB
  ZPHIF=0._JPRB

  ZDEEP=2500._JPRB
  ZTHETAM= 8._JPRB

  ! 2 foci points N/S from equator, determine the width of the cold pool
  ! centre of ellipse at equator
  ZLATC =  0._JPRB
  ZLONC =  0._JPRB
  ! long distance of ellipse, x-direction
  ZA = 50000.0_JPRB
  ! long distance of ellipse, y-direction
  ZB =100000._JPRB

! compute geopotentiel of reference atmosphere
  CALL GPGEO (NGPTOT,1,NGPTOT,NFLEVG,&
              & ZPHI, ZPHIF, &
              & ZT,ZRDGAZ,ZXYB(1,1,YYTXYB%M_LNPR),&
              & ZXYB(1,1,YYTXYB%M_ALPH), YDVERT_GEOM)

DO JWORD = 1,NGPTOT
DO JLEV = 1, NFLEVG
ZZZ=ZPHIF(JWORD,JLEV)/RG
      IF (YDGSGEOM_NB%GELAM(JWORD)>RPI) THEN
           ZGELAM(JWORD)= YDGSGEOM_NB%GELAM(JWORD)-2._JPRB*RPI
      ELSE
           ZGELAM(JWORD)= YDGSGEOM_NB%GELAM(JWORD)
      ENDIF
    ! projection on x, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(ZLATC)+&
     &COS(ZLATC)*COS(ZLATC)*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
    ZRADX = RA*ACOS(ZZ)      
    ! projection on y, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
     &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(ZLONC-ZLONC)
    ZRADY = RA*ACOS(ZZ)
IF ((ZZZ<=ZDEEP) .AND. (ZGELAM(JWORD)<0._JPRB)) THEN
ZT(JWORD,JLEV) = ZT(JWORD,JLEV) + (ZTHETAM/ZDEEP*ZZZ-ZTHETAM )&
 & * (ZPRESF(JWORD,JLEV)/100000._JPRB)**(RD/RCPD)&
 & /(sqrt(1._JPRB+(ZRADX/ZA)**2+(ZRADY/ZB)**2))
ENDIF
ENDDO
ENDDO

! T to spectral space
    CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PTEMP,ZT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Wind profile at the equator (U=0 at surf, U=ZU0 at ZDEEP)
ZU0= RU00_DYN
ZU = ZU0

DO JWORD = 1,NGPTOT
DO JLEV = 1, NFLEVG
ZZZ=ZPHIF(JWORD,JLEV)/RG
IF (ZZZ<=ZDEEP) THEN
ZU(JWORD,JLEV) =  ZU0/ZDEEP*ZZZ
ENDIF
ENDDO
ENDDO

! Wind U=0 toward the poles (solid body rotation)
DO JWORD = 1,NGPTOT
DO JLEV = 1, NFLEVG
ZU(JWORD,JLEV) = (ZU(JWORD,JLEV)-ZST)*YDGSGEOM_NB%GSQM2(JWORD)
ZV(JWORD,JLEV) = 0._JPRB
ENDDO
ENDDO

CALL UVSPE(YDGEOMETRY,PVOR,PDIV,ZU,ZV,NFLEVL,NFLEVG,1)

ELSE
  WRITE(NULERR,'('' INVALID SETUP NTESTCASE FOR MISC FAMILY'')')
  CALL ABOR1(' INVALID SETUP NTESTCASE FOR MISC FAMILY')
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUMISC_SPEC',1,ZHOOK_HANDLE)
END SUBROUTINE SUMISC_SPEC
