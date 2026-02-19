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

SUBROUTINE SUGRIDUG2(YDGEOMETRY,YDML_GCONF,YDSPEC,YDGFL)

!**** *SUGRIDUG2*  - Initialize the upper air grid-point fields for N3DINI=2

!     Purpose.
!     --------
!           Initialize the upper air gridpoint fields of the model for N3DINI=2

!**   Interface.
!     ----------
!        *CALL* *SUGRIDUG2

!        Explicit arguments :
!        ------------------

!        Implicit arguments :
!        --------------------

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
!        Nils Wedi *ECMWF*

!     Modifications.
!     --------------
!              adapted from sugridug.F90
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      S.Malardel: Jan. 2016 Cleaning + DCMIP12 cases
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULERR 
USE YOMVAR   , ONLY : LSPINT
USE YOM_GRIB_CODES  , ONLY : NGRBQ
USE YOMDYNCORE, ONLY : CTEST_FAMILY, NTESTCASE
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD

IMPLICIT NONE

! * INPUT:
TYPE(GEOMETRY), INTENT(IN)         :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC
TYPE(TGFL),     INTENT(INOUT)      :: YDGFL

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IGRIB3D(YDML_GCONF%YGFL%NUMFLDS)
INTEGER(KIND=JPIM) :: ICL, JGFL, ILEN, JSTGLO, ICEND, ISTC, IBL, ICL2
INTEGER(KIND=JPIM) :: JLEV, J, IQ
REAL(KIND=JPRB),ALLOCATABLE :: ZZGFL(:,:,:,:)
REAL(KIND=JPRB)    :: ZQ(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)

LOGICAL :: LLSPINT, LLINCR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "sudcmip12_gu.intfb.h"
#include "sudcmip16_gu.intfb.h"
#include "sumisc_gu.intfb.h"
#include "gpnorm_gfl.intfb.h"
#include "abor1.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGRIDUG2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GFL=>YDGFL%GFL)
!     ------------------------------------------------------------------

WRITE(UNIT=NULOUT,FMT='('' INITIALISING ARTIFICIAL 3D GRID-POINT DATA FOR: '')')
WRITE(UNIT=NULOUT,FMT='('' CTEST_FAMILY = '',A)') CTEST_FAMILY
WRITE(UNIT=NULOUT,FMT='('' NTESTCASE    = '',I3)') NTESTCASE

!==================================================================
LLINCR  = .FALSE.
LLSPINT = LSPINT.AND.LLINCR

ICL = 0
IQ = 0

IF(NUMFLDS > 0 ) THEN
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LGP .AND. YCOMP(JGFL)%NREQIN == 1 .AND. YCOMP(JGFL)%LGPINGP&
       & .AND. .NOT. LLSPINT) THEN  
      ICL = ICL+1
      IGRIB3D(ICL) = YCOMP(JGFL)%IGRBCODE
      IF( YCOMP(JGFL)%IGRBCODE == NGRBQ ) THEN
        IQ=ICL
      ENDIF
    ENDIF
  ENDDO
ENDIF

IF (ICL /= 0) THEN

  ALLOCATE (ZZGFL(NGPTOT,NFLEVG,ICL,1))

!==================================================================
! Initial setup
  ! set all gfl to 0.
  ZZGFL(:,:,1:ICL,:) = 0._JPRB
  ZQ=0.0_JPRB
!==================================================================
  IF      (CTEST_FAMILY == "DCMIP12     ") THEN
    CALL SUDCMIP12_GU(YDGEOMETRY,YDML_GCONF%YRDIMF,YDSPEC,NTESTCASE,ZQ)

  ELSEIF  (CTEST_FAMILY == "DCMIP16     ") THEN
    CALL SUDCMIP16_GU(YDGEOMETRY,YDML_GCONF%YRDIMF,YDSPEC,NTESTCASE,ZQ)

  ELSEIF  (CTEST_FAMILY == "MISC        ") THEN
    CALL SUMISC_GU(YDGEOMETRY,YDML_GCONF%YRDIMF,YDSPEC,NTESTCASE,ZQ)
  ELSE
    WRITE(NULERR,'(''  IN SUGRIDUG2'')')
    CALL ABOR1(' INVALID SETUP for CTEST_FAMILY IN SUGRIDUG2 ')
  ENDIF

! set q profile
  IF (IQ /= 0) THEN
    DO JLEV=1,NFLEVG
      ZZGFL(1:NGPTOT,JLEV,IQ,1) = ZQ(1:NGPTOT,JLEV)
    ENDDO
  ENDIF

!*       ADD TO GFL GRIDPOINT
!        --------------------

  CALL GSTATS(1411,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSTGLO,ICEND,ISTC,JLEV,J,&
!$OMP& JGFL,IBL,ICL2)
  DO JSTGLO=1,NGPTOT,NPROMA
    ICEND = MIN(NPROMA,NGPTOT-JSTGLO+1)
    ISTC  = 1
    IBL   = (JSTGLO-1)/NPROMA+1
    ICL2  = 0
    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LGP .AND. YCOMP(JGFL)%NREQIN == 1 .AND. YCOMP(JGFL)%LGPINGP&
     & .AND. .NOT. LLSPINT) THEN  
        ICL2 = ICL2+1
        DO JLEV=1,NFLEVG
          DO J=ISTC,ICEND
            GFL(J,JLEV,YCOMP(JGFL)%MP,IBL) = ZZGFL(JSTGLO+J-1,JLEV,ICL2,1)
          ENDDO
        ENDDO
      ELSEIF(YCOMP(JGFL)%LGP .AND. .NOT. YCOMP(JGFL)%NREQIN == 1) THEN  
      GFL(:,:,YCOMP(JGFL)%MP,IBL) = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1411,1)

  IF (ALLOCATED(ZZGFL)) DEALLOCATE (ZZGFL)

  WRITE(NULOUT,'(A)') ' SUGRIDUG2 : STATISTICS FOR ALL GFL FIELDS after SUGRIDUG2'
  CALL GPNORM_GFL(YDGEOMETRY,YDGFL)

! end test if some GFLs have to be setup (ICL /= 0)
ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGRIDUG2',1,ZHOOK_HANDLE)

END SUBROUTINE SUGRIDUG2
