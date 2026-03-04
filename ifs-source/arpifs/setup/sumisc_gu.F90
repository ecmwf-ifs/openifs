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

SUBROUTINE SUMISC_GU(YDGEOMETRY,YDDIMF,YDSPEC,KTESTCASE,PQ)

!     Purpose.
!     --------
!           uuper air grid point fields for MISC test cases

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
USE YOMLUN    , ONLY : NULOUT, NULERR
USE YOMDIMF   , ONLY : TDIMF
USE YOMCST    , ONLY : RPI, RA, RD, RG, ROMEGA, RCPD, RATM
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD
USE INTDYN_MOD, ONLY : YYTXYB
IMPLICIT NONE

TYPE(GEOMETRY),     INTENT(IN)     :: YDGEOMETRY
TYPE(TDIMF)        ,INTENT(INOUT)  :: YDDIMF
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC
INTEGER(KIND=JPIM), INTENT(IN)     :: KTESTCASE
REAL(KIND=JPRB),    INTENT(INOUT)  :: PQ(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)

! squall line
REAL(KIND=JPRB) ::    ZPRESH(YDGEOMETRY%YRGEM%NGPTOT,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPRESF (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZXYB(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)

REAL(KIND=JPRB) ::    ZTEMP(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ::    ZLNSP(YDGEOMETRY%YRDIM%NSPEC2),ZSPOROG(YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) ::    ZT(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPSP(YDGEOMETRY%YRGEM%NGPTOT),ZGPOROG(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ::    ZPHI(YDGEOMETRY%YRGEM%NGPTOT,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPHIF(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZRDGAZ(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZRADX,ZRADY,ZLATC,ZLONC,ZA,ZB,ZZ,ZZZ
REAL(KIND=JPRB) ::    ZGELAM(YDGEOMETRY%YRGEM%NGPTOT)
!Klemp profiles
REAL(KIND=JPRB) ::    ZPREF_K(41)
REAL(KIND=JPRB) ::    ZQ_K(41)

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JLEV, JWORD, JK

!     ------------------------------------------------------------------

!#include "abor1.intfb.h"
#include "gphpre.intfb.h"
#include "speree.intfb.h"
#include "gpgeo.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMISC_GU',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, &
 & YDCVER=>YDGEOMETRY%YRCVER,YDVERT_GEOM=>YDGEOMETRY%YRVERT_GEOM)
ASSOCIATE(NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & NGPTOT=>YDGEM%NGPTOT)

!==================================================================
! Case 91:  Squall line
!==================================================================
IF((KTESTCASE == 91) .OR. (KTESTCASE == 92)) THEN

!! compute the necessary pressure information before interpolation

! Get surface pressure and T from spectral space
  ZLNSP(:)=YDSPEC%SP(:)
  ZTEMP(:,:)=YDSPEC%T(:,:)
  ZSPOROG(:)=YDSPEC%OROG(:)
  CALL SPEREE(YDGEOMETRY,1,1,ZLNSP,ZPSP)
  CALL SPEREE(YDGEOMETRY,1,1,ZSPOROG,ZGPOROG)
  CALL SPEREE(YDGEOMETRY,NFLEVL,NFLEVG,ZTEMP,ZT)
  DO JWORD = 1,NGPTOT
    ZPRESH(JWORD,NFLEVG) = EXP( ZPSP(JWORD) )
    ZPHI(JWORD,NFLEVG) = ZGPOROG(JWORD)
  ENDDO

  CALL GPHPRE (NGPTOT,NFLEVG,1,NGPTOT,YDVAB,YDCVER,ZPRESH,PXYB=ZXYB,PRESF=ZPRESF)

! compute geopotentiel
  ZRDGAZ = RD
  CALL GPGEO (NGPTOT,1,NGPTOT,NFLEVG,&
              & ZPHI, ZPHIF,&
              & ZT,ZRDGAZ,ZXYB(1,1,YYTXYB%M_LNPR),&
              & ZXYB(1,1,YYTXYB%M_ALPH),YDVERT_GEOM)

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

ZQ_K(1)=   0.013807
ZQ_K(2)=   0.002353
ZQ_K(3)=   0.000671
ZQ_K(4)=   0.000308
ZQ_K(5)=   0.000175
ZQ_K(6)=   0.000112
ZQ_K(7)=   0.000078
ZQ_K(8)=   0.000058
ZQ_K(9)=   0.000045
ZQ_K(10)=  0.000036
ZQ_K(11)=  0.000029
ZQ_K(12)=  0.000024
ZQ_K(13)=  0.000043
ZQ_K(14)=  0.000090
ZQ_K(15)=  0.000167
ZQ_K(16)=  0.000285
ZQ_K(17)=  0.000455
ZQ_K(18)=  0.000686
ZQ_K(19)=  0.000989
ZQ_K(20)=  0.001372
ZQ_K(21)=  0.001842
ZQ_K(22)=  0.002406
ZQ_K(23)=  0.003066
ZQ_K(24)=  0.003824
ZQ_K(25)=  0.004679
ZQ_K(26)=  0.005628
ZQ_K(27)=  0.006666
ZQ_K(28)=  0.007786
ZQ_K(29)=  0.008977
ZQ_K(30)=  0.010227
ZQ_K(31)=  0.011522
ZQ_K(32)=  0.012845
ZQ_K(33)=  0.013807
ZQ_K(34)=  0.013807
ZQ_K(35)=  0.013807
ZQ_K(36)=  0.013807
ZQ_K(37)=  0.013807
ZQ_K(38)=  0.013807
ZQ_K(39)=  0.013807
ZQ_K(40)=  0.013807
ZQ_K(41)=  0.013807

! vertical interpolation on the current full levels 

DO JWORD = 1,NGPTOT

DO JLEV = 1, NFLEVG
  
  DO JK=1,41

    IF (ZPRESF(JWORD,JLEV) < ZPREF_K(JK)) THEN
! JLEV is between JK-1 and JK (if JK/=1 !)

      IF (JK==1) THEN
! level above the first full pressure level of original profil
! extrapolation above
      PQ(JWORD,JLEV) = ZQ_K(1)+ (ZQ_K(2)-ZQ_K(1))&
                     & * LOG(ZPRESF(JWORD,JLEV)/ZPREF_K(1))&
                     & / LOG(ZPREF_K(2)/ZPREF_K(1))

      EXIT
      ELSE
      PQ(JWORD,JLEV) = ZQ_K(JK-1)+ (ZQ_K(JK)-ZQ_K(JK-1))&
                     & * LOG(ZPRESF(JWORD,JLEV)/ZPREF_K(JK-1))&
                     & / LOG(ZPREF_K(JK)/ZPREF_K(JK-1))
      EXIT
      ENDIF

    ELSE 
! level under the last full pressure level of original profil
! extrapolation below
      PQ(JWORD,JLEV) = ZQ_K(41)+ (ZQ_K(40)-ZQ_K(41))&
                     & * LOG(ZPRESF(JWORD,JLEV)/ZPREF_K(41))&
                     & / LOG(ZPREF_K(40)/ZPREF_K(41))
    ENDIF

  ENDDO

ENDDO

ENDDO

IF (KTESTCASE==92) THEN
! reduce humidity outside the cold pool region

  ! 2 foci points N/S from equator, determine the width of the cold pool
  ! centre of ellipse at equator
  ZLATC =  0._JPRB
  ZLONC =  0._JPRB
  ! long distance of ellipse, x-direction
  ZA = 50000.0_JPRB
  ! long distance of ellipse, y-direction
  ZB =100000._JPRB
DO JWORD = 1,NGPTOT
DO JLEV = 1, NFLEVG
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

IF ((ZGELAM(JWORD)>0._JPRB)) THEN
    PQ(JWORD,JLEV)=PQ(JWORD,JLEV)*0.8_JPRB
ELSE
    PQ(JWORD,JLEV)= PQ(JWORD,JLEV)* &
 &    max(0.8_JPRB,1._JPRB/(sqrt(1._JPRB+(ZRADX/ZA)**2+(ZRADY/ZB)**2)) )
ENDIF
ENDDO
ENDDO
ENDIF


  WRITE(NULOUT,'(''END OF GFL setup in sugridug2 for NTESTCASE=91/92'')')
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUMISC_GU',1,ZHOOK_HANDLE)
END SUBROUTINE SUMISC_GU
