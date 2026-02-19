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

SUBROUTINE GETNEMODIAG3D(KGPTOT,KPROMA,KSTEP,YDSURF,YDMCC,YDGEOMETRY)
!
!**** *GETNEMODIAG*  - Retrieve diagnostics (from IFS POV) from NEMO.
!
!     Purpose.
!     --------
!       Interpolate NEMO data onto the reduced gaussian
!       grid and to surface data structre
!
!**   Interface.
!     ----------
!       *CALL*  *GETNEMODIAG*
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!     Method:
!     ------
!       
!     Externals:
!     ---------
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!       K. Mogensen, ECMWF
!
!     Modifications.
!     --------------
!     -----------------------------------------------------------
   
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMMCC   , ONLY : TMCC
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE PARKIND_OCEAN, ONLY : JPRO
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE MPL_MODULE, ONLY : MPL_COMM
USE YOMMP0   , ONLY : MYPROC, NPROC

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KGPTOT
INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP
TYPE(TSURF)     ,INTENT(INOUT) :: YDSURF
TYPE(TMCC)         ,INTENT(IN) :: YDMCC
TYPE(GEOMETRY)    , INTENT(IN) :: YDGEOMETRY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRO), DIMENSION(KGPTOT,YDMCC%NNEMO3DLEVS) :: ZGT, ZGS, ZGU, ZGV
INTEGER(KIND=JPIM) :: JSTGLO, IEND, IST, IBL, JROF, JK, JV
REAL(KIND=JPRB), PARAMETER :: ZMDI=HUGE(ZMDI)

#include "gpnorm3.intfb.h"

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETNEMODIAG',0,ZHOOK_HANDLE)
ASSOCIATE(SP_OM=>YDSURF%SP_OM,YSP_OM=>YDSURF%YSP_OM,YSP_OMD=>YDSURF%YSP_OMD,&
  &       SD_VF=>YDSURF%SD_VF,YSD_VF=>YDSURF%YSD_VF,&
  &       RNEMOMASK=>YDMCC%RNEMOMASK)

! Regrid the NEMO data.

#ifdef WITH_NEMO
IF(KSTEP==0) THEN

  WRITE(NULOUT,*)
  WRITE(NULOUT,*)'INITIAL 3D DIAGNOSTIC FIELDS RETRIEVED FROM NEMO.'
  WRITE(NULOUT,*)

ELSE

  WRITE(NULOUT,*)
  WRITE(NULOUT,*)'3D DIAGNOSTIC FIELDS RETRIEVED FROM NEMO.'
  WRITE(NULOUT,*)

ENDIF

ZGT(:,:)=0.0_JPRO
ZGS(:,:)=0.0_JPRO
ZGU(:,:)=0.0_JPRO
ZGV(:,:)=0.0_JPRO
CALL NEMOGCMCOUP_MLFLDS_GET( MYPROC-1, NPROC, MPL_COMM, YDMCC%NNEMO3DLEVS, &
  & KGPTOT, ZGT, ZGS, ZGU, ZGV )

DO JK=1,YSP_OMD%NLEVS
  DO JSTGLO=1,KGPTOT,KPROMA
    IEND=MIN(KPROMA,KGPTOT-JSTGLO+1)
    IST=1
    IBL=(JSTGLO-1)/KPROMA+1
    DO JROF=1,IEND
      SP_OM(JROF,JK,YSP_OM%YTO%MP,IBL) = ZGT(JSTGLO+JROF-1,JK)
      SP_OM(JROF,JK,YSP_OM%YSO%MP,IBL) = ZGS(JSTGLO+JROF-1,JK)
      SP_OM(JROF,JK,YSP_OM%YUO%MP,IBL) = ZGU(JSTGLO+JROF-1,JK)
      SP_OM(JROF,JK,YSP_OM%YVO%MP,IBL) = ZGV(JSTGLO+JROF-1,JK)
    ENDDO
  ENDDO
ENDDO

WRITE(NULOUT,'(A)') ' GETNEMODIAG3D: STATISTICS FOR ALL NEMO OCEAN LEVEL FIELDS'
IF (YSP_OMD%NUMFLDS>0) THEN
  DO JV=1,YSP_OMD%NUMFLDS
    WRITE(NULOUT,'(A,I4)') ' NEMO PROGNOSTIC FIELDS ',JV
    CALL GPNORM3(YDGEOMETRY,SP_OM,YSP_OMD%NUMFLDS,JV,KLEVS=YSP_OMD%NLEVS,LDLEVELS=.TRUE.)
  ENDDO
ENDIF
WRITE(NULOUT,*)

DO JK=1,YSP_OMD%NLEVS
  DO JSTGLO=1,KGPTOT,KPROMA
    IEND=MIN(KPROMA,KGPTOT-JSTGLO+1)
    IST=1
    IBL=(JSTGLO-1)/KPROMA+1
    DO JROF=1,IEND
      IF ((SD_VF(JROF,YSD_VF%YCLK%MP,IBL)>0.5_JPRB) &
        & .OR.(RNEMOMASK(JSTGLO+JROF-1,JK)<0.99999_JPRB)) THEN
        SP_OM(JROF,JK,YSP_OM%YTO%MP,IBL) = ZMDI
        SP_OM(JROF,JK,YSP_OM%YSO%MP,IBL) = ZMDI
        SP_OM(JROF,JK,YSP_OM%YUO%MP,IBL) = ZMDI
        SP_OM(JROF,JK,YSP_OM%YVO%MP,IBL) = ZMDI
     ENDIF
    ENDDO
  ENDDO
ENDDO

#else

CALL ABOR1('ININEMO: COMPILED WITHOUT WITH_NEMO')

#endif

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GETNEMODIAG',1,ZHOOK_HANDLE)

!     -----------------------------------------------------------

END SUBROUTINE GETNEMODIAG3D
