! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUPHYDS1C(YDMODEL)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULNAM
USE YOMGPD1C , ONLY : VEXTR2   ,VEXTRA

#ifdef DOC

!**** *SUPHYDS1C*  - Routine to setup physics fields desciptors

!     Purpose.
!     --------
!           Setup the discriptors for initializing the physics
!        fields including default values for the different fields,
!        ARPEGE field names and grib codes.

!**   Interface.
!     ----------
!        *CALL* *SUPHYDS1C*

!        Explicit arguments : None
!        --------------------

!        Implicit arguments : Common YOMPHYDS
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.     None
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-08-01
!        R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        R. EL Khatib : 94-03-04 auxilary climatological fields
!        M.Deque      : 94-10-12 4 layer soil temperature
!        M.Hamrud     : 94-06-27 GRIB code variables from YOMGRB
!        EB and DG    : 94-10-07 implementation of ISBA
!        H. Douville  : 95-01-13 Snow parameterization
!        J. Bidlot    : 97-06-12 wave field definition
!        E. Bazile    : 97-11-18 Soil freezing (LFGEL)
!        C. Madeira   : 97-07-23 ISBA Fields in FULL-POS
!        R. EL Khatib : 98-04-15 Move print out of NFLDPTP from SUGRIDASM/DM
!        R. EL Khatib : 98-08-31 Adiabatic pp.
!        E. Bazile    : 99-02-12 Surface Soil freezing (LFGELS).
!        J. Teixeira  : 99-03-26 Surface tiling
!        A.Beljaars   : 00-03-06 New fields:TSRC,TTRC,SSRC,STRC,ES,SMLT,10FG,LSPF
!        D.Dent       : 00-12-20 Define character names for 2D extra fields
!        D. Giard     : 00-10-25 Albedos for vegetation and bare-ground
!        R. El Khatib : 01-08-07 Pruning options
!        D.Dent       : 02-05-15 Checks for extra field codes
!        JJMorcrette  : 02-09-30 ECMWF PAR, UV, CAPE
!        Y. Bouteloup : 02-03-29 Climatological ozone profiles (xxO3ABC)
!        E. Bazile    : 02-09-04 New snow scheme (LVGSN).
!        R. El Khatib : 03-04-17 Cleanups
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Marquet   : 02-07-17 add VCLIA for aerosol files
!        G. Hello     : 05-04-25 Lground options
!        M. Ko"hler   : 6-6-2006 Single Column Model integration within IFS 
!     ------------------------------------------------------------------
#endif

USE TYPE_MODEL, ONLY : MODEL

IMPLICIT NONE

TYPE(MODEL), INTENT(INOUT) :: YDMODEL

INTEGER(KIND=JPIM) :: JF

character (len=2)  :: ntxt

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) , POINTER ::  NVEXTRAGB(:)
INTEGER(KIND=JPIM) , POINTER ::  NSFLUX
REAL(KIND=JPRB) , POINTER ::  VEXTRADF(:)
REAL(KIND=JPRB) , POINTER ::  EXTRPDF(:,:)
CHARACTER (LEN = 16) , POINTER ::   CVEXTR2(:)
INTEGER(KIND=JPIM) , POINTER ::  NVEXTRDIGB(:)
REAL(KIND=JPRB) , POINTER ::  VEXTRDI(:)
REAL(KIND=JPRB) , POINTER ::  VEXRADDF(:)
LOGICAL , POINTER ::  LREQIN_VEXTRDI
CHARACTER (LEN = 16) , POINTER ::   CVEXTRDI(:)
INTEGER(KIND=JPIM) , POINTER ::  NEXTRPGB(:,:)
CHARACTER (LEN = 16) , POINTER ::   CVEXTRA(:)
CHARACTER (LEN = 16) , POINTER ::   CVEXRAD(:)
INTEGER(KIND=JPIM) , POINTER ::  NSFORC
CHARACTER (LEN = 16) , POINTER ::   CEXTRP(:,:)
INTEGER(KIND=JPIM) , POINTER ::  NVEXRADGB(:)
LOGICAL , POINTER ::  LREQOUT_VEXTRDI
INTEGER(KIND=JPIM) , POINTER ::  NVEXTR2GB(:)
REAL(KIND=JPRB) , POINTER ::  VEXTR2DF(:)


#include "namphyds.nam.h"

!     ------------------------------------------------------------------
!*       1.    SET DEFAULTS.
!              ------------

!        1.1 Set implicit default values

IF (LHOOK) CALL DR_HOOK('SUPHYDS1C',0,ZHOOK_HANDLE)
ASSOCIATE(YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS)
ASSOCIATE(NXTRP2GB=>YDPHYDS%NXTRP2GB, JPVXP2=>YDPHYDS%JPVXP2, &
 & JPVEXTRDI=>YDPHYDS%JPVEXTRDI, JPVEXTR=>YDPHYDS%JPVEXTR, &
 & JPVXP=>YDPHYDS%JPVXP, JPVXTR2=>YDPHYDS%JPVXTR2, XTRP2DF=>YDPHYDS%XTRP2DF, &
 & CXTRP2=>YDPHYDS%CXTRP2, JPCXP=>YDPHYDS%JPCXP, &
 & CEXTRP => YDPHYDS%CEXTRP, &
 & VEXTRDI => YDPHYDS%VEXTRDI, &
 & NVEXRADGB => YDPHYDS%NVEXRADGB, &
 & NVEXTR2GB => YDPHYDS%NVEXTR2GB, &
 & EXTRPDF => YDPHYDS%EXTRPDF, &
 & LREQOUT_VEXTRDI => YDPHYDS%LREQOUT_VEXTRDI, &
 & CVEXTRDI => YDPHYDS%CVEXTRDI, &
 & NSFLUX => YDPHYDS%NSFLUX, &
 & VEXTRADF => YDPHYDS%VEXTRADF, &
 & CVEXTR2 => YDPHYDS%CVEXTR2, &
 & LREQIN_VEXTRDI => YDPHYDS%LREQIN_VEXTRDI, &
 & NVEXTRAGB => YDPHYDS%NVEXTRAGB, &
 & CVEXRAD => YDPHYDS%CVEXRAD, &
 & VEXTR2DF => YDPHYDS%VEXTR2DF, &
 & VEXRADDF => YDPHYDS%VEXRADDF, &
 & NSFORC => YDPHYDS%NSFORC, &
 & CVEXTRA => YDPHYDS%CVEXTRA, &
 & NVEXTRDIGB => YDPHYDS%NVEXTRDIGB, &
 & NEXTRPGB => YDPHYDS%NEXTRPGB)


DO JF=1,JPVEXTR
  write(ntxt,"(2i1)") int(jf/10), mod(jf,10)
  CVEXTRA  (JF) = 'Extra Column Variable #'//ntxt
  VEXTRADF (JF) = 0.0_JPRB
ENDDO

DO JF=1,JPVXTR2
  write(ntxt,"(2i1)") int(jf/10), mod(jf,10)
  CVEXTR2  (JF) = 'Extra Surface Variable #'//ntxt
  VEXTR2DF (JF) = 0.0_JPRB
ENDDO

!     ------------------------------------------------------------------
!*       2.    READ NAMELIST.
!              -------------

CALL POSNAM(NULNAM,'NAMPHYDS')
READ(NULNAM,NAMPHYDS)

END ASSOCIATE
END ASSOCIATE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUPHYDS1C',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHYDS1C
