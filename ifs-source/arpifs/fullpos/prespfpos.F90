! (C) Copyright 1989- Meteo-France.

SUBROUTINE PRESPFPOS(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,YDSPEC)

!**** *PRESPFPOS*  - Prepare before dynamic post-processing

!     Purpose.
!     --------
!        To pack the data before post-processing
!         - pack spectral arrays
!         - pack gridpoint upperair fields

!**   Interface.
!     ----------
!        *CALL* *PRESPFPOS

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    
!     ----------    
!        SPACONVERT - convert model fields spectral array into file fields
!                     spectral arrays and reverse
!        PKGRIDA  - pack gridpoint arrays/workfile - Arpege
!        PKSPECA  - pack spectral arrays/workfile - Arpege

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib  *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 98-01-28 from STEPO
!     R. El Khatib : 01-03-16 Control with NGRSP for aladin
!     R. El Khatib 02-03-19 lagged Fourier diffusion before packing
!     R. El Khatib : 01-08-07 Pruning options
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     K.Yessad (Dec 2003): cleaning in horizontal diffusion.
!     R. El Khatib : 04-12-08 remove broken backup of upper-air gp fields
!     R. El Khatib  25-Apr-2005 SPA* not saved if 1st part 927 (Only 1 scan)
!        M.Hamrud      01-Jul-2006  Revised surface fields
!      R. El Khatib 22-Mar-2012 Fix uninitialized variables
!      R. El Khatib 16-Jul-2012 Fullpos move away from STEPO.
!      R. El Khatib 31-Jul-2012 Split into PREDYNFPOS + PRESPFPOS
!      R. El khatib 16-May-2014 Optimization of in-line/off-line post-processing reproducibility
!      R. El khatib 24_Sep-2021 Bugfix : add gridpoint GFL packing
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYNA  , ONLY : TDYNA
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE YOMFA    , ONLY : YFAL, YFAI, YFAR, YFAS, YFALRAD, YFAIRAD, YFACLF, YFATKE, YFAG, YFAH, YFAO3, YFACPF, YFASPF

IMPLICIT NONE

!     ------------------------------------------------------------------

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(TDYNA)         ,INTENT(IN)    :: YDDYNA
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC


REAL(KIND=JPRB) :: ZSPEC(YDGEOMETRY%YRDIM%NSPEC2,YDGEOMETRY%YRDIMV%NFLEVL*YDML_GCONF%YRDIMF%NS3D+YDML_GCONF%YRDIMF%NS2D)
REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:)

INTEGER(KIND=JPIM) :: IFLD, JF, JLEV, JS, JKGLO, IBL, IEND, IGPGFL
INTEGER(KIND=JPIM), ALLOCATABLE :: IBITGFL(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "pkspeca.intfb.h"
#include "pkgrida.intfb.h"
#include "spaconvert.intfb.h"

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PRESPFPOS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YGFL=>YDML_GCONF%YGFL, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NSPEC2=>YDDIM%NSPEC2, NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NGPTOT=>YDGEM%NGPTOT, &
 & NS2D=>YDDIMF%NS2D, NS3D=>YDDIMF%NS3D, NUMFLDS=>YGFL%NUMFLDS, YDCOMP=>YGFL%YCOMP, &
 & NFLEVL=>YDDIMV%NFLEVL, NFLEVG=>YDDIMV%NFLEVG, &
 & NBSETSP=>YDMP%NBSETSP)
!      -----------------------------------------------------------

!*       1.    VARIOUS INITIALIZATION
!              ----------------------

!   Convert spectral field before packing
    CALL SPACONVERT(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,.TRUE.,YDSPEC)

!   Spectral packing
    IFLD=0
    DO JF=1,NS3D
      DO JLEV=1,NFLEVL
        IFLD=IFLD+1
        DO JS=1,NSPEC2
          ZSPEC(JS,IFLD)=YDSPEC%SP3D(JLEV,JS,JF)
        ENDDO
      ENDDO
    ENDDO
    IF (MYSETV==NBSETSP) THEN
      DO JF=1,NS2D
        IFLD=IFLD+1
        DO JS=1,NSPEC2
          ZSPEC(JS,IFLD)=YDSPEC%SP2D(JS,JF)
        ENDDO
      ENDDO
    ENDIF

    CALL PKSPECA(YDGEOMETRY,ZSPEC,IFLD)

    IFLD=0
    DO JF=1,NS3D
      DO JLEV=1,NFLEVL
        IFLD=IFLD+1
        DO JS=1,NSPEC2
          YDSPEC%SP3D(JLEV,JS,JF)=ZSPEC(JS,IFLD)
        ENDDO
      ENDDO
    ENDDO
!   Here below : orography is packed as the others. This is not safe, 
!   but it is correct since orography is not used in spectral space
!   in the post-processing and SUOROG is not called.
    IF (MYSETV==NBSETSP) THEN
      DO JF=1,NS2D
        IFLD=IFLD+1
        DO JS=1,NSPEC2
          YDSPEC%SP2D(JS,JF)=ZSPEC(JS,IFLD)
        ENDDO
      ENDDO
    ENDIF

!   Gridpoint packing for gridpoint GFL :
    IGPGFL=NFLEVG*COUNT(YDCOMP(1:NUMFLDS)%LGP.AND.YDCOMP(1:NUMFLDS)%LREQOUT)
    ALLOCATE(IBITGFL(IGPGFL))
    ALLOCATE(ZGP(NGPTOT,IGPGFL))
    IFLD=0
    DO JF=1,NUMFLDS
      IF (YDCOMP(JF)%LGP.AND.YDCOMP(JF)%LREQOUT) THEN
        DO JLEV=1,NFLEVG
          IFLD=IFLD+1
          IF (YDCOMP(JF)%CNAME == YFAL%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAL%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFAI%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAI%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFALRAD%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFALRAD%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFAIRAD%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAIRAD%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFACLF%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFACLF%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFAR%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAR%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFAS%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAS%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFAG%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAG%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFAH%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAH%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFATKE%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFATKE%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFAO3%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFAO3%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFACPF%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFACPF%NBITS
          ELSEIF (YDCOMP(JF)%CNAME == YFASPF%CLNAME(1:12)) THEN
            IBITGFL(IFLD) = YFASPF%NBITS
          ELSE
            IBITGFL(IFLD) = -1
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    DO JKGLO=1,NGPTOT,NPROMA
      IBL=(JKGLO-1)/NPROMA+1
      IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IFLD=0
      DO JF=1,NUMFLDS
        IF (YDCOMP(JF)%LGP.AND.YDCOMP(JF)%LREQOUT) THEN
          DO JLEV=1,NFLEVG
            IFLD=IFLD+1
            ZGP(JKGLO:JKGLO-1+IEND,IFLD)=YDGFL%GFL(1:IEND,JLEV,YDCOMP(JF)%MP,IBL)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CALL PKGRIDA(YDGEOMETRY,IGPGFL,IBITGFL,ZGP)
    DO JKGLO=1,NGPTOT,NPROMA
      IBL=(JKGLO-1)/NPROMA+1
      IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IFLD=0
      DO JF=1,NUMFLDS
        IF (YDCOMP(JF)%LGP.AND.YDCOMP(JF)%LREQOUT) THEN
          DO JLEV=1,NFLEVG
            IFLD=IFLD+1
            YDGFL%GFL(1:IEND,JLEV,YDCOMP(JF)%MP,IBL)=ZGP(JKGLO:JKGLO-1+IEND,IFLD)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(ZGP)
    
!   Revert conversion of spectral fields
    CALL SPACONVERT(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,.FALSE.,YDSPEC)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PRESPFPOS',1,ZHOOK_HANDLE)

END SUBROUTINE PRESPFPOS
