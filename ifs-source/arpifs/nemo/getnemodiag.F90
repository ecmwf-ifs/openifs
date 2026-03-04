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

SUBROUTINE GETNEMODIAG(KGPTOT,KPROMA,KSTEP,YDSURF,YDMCC)
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
USE YOMMCC  , ONLY : TMCC
USE PARKIND1 , ONLY : JPRD, JPIM, JPRB
USE PARKIND_OCEAN, ONLY : JPRO
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE MPL_MODULE, ONLY : MPL_COMM
USE YOMMP0   , ONLY : MYPROC, NPROC

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KGPTOT
INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP
TYPE(TSURF), INTENT(INOUT)     :: YDSURF
TYPE(TMCC), INTENT(INOUT)      :: YDMCC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRO), DIMENSION(KGPTOT) :: ZGSSH, ZGMLD, ZG20D, ZGSSS,&
 & ZGTEM300, ZGSAL300
INTEGER(KIND=JPIM) :: JSTGLO, IEND, IST, IBL, JROF
REAL(KIND=JPRB), PARAMETER :: ZMDI=HUGE(ZMDI)

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETNEMODIAG',0,ZHOOK_HANDLE)

ASSOCIATE(SD_OC=>YDSURF%SD_OC, YSD_OC=>YDSURF%YSD_OC, &
  &       SD_VF=>YDSURF%SD_VF,YSD_VF=>YDSURF%YSD_VF, &
  &       CPLNG_FLD=>YDMCC%CPLNG_FLD)

! Regrid the NEMO data.

#ifdef WITH_NEMO

IF(KSTEP==0) THEN

   WRITE(NULOUT,*)
   WRITE(NULOUT,*)'INITIAL DIAGNOSTIC FIELDS RETRIEVED FROM NEMO.'
   WRITE(NULOUT,*)

ELSE

   WRITE(NULOUT,*)
   WRITE(NULOUT,*)'DIAGNOSTIC FIELDS RETRIEVED FROM NEMO.'
   WRITE(NULOUT,*)

ENDIF

CALL NEMOGCMCOUP_EXFLDS_GET( MYPROC-1, NPROC, MPL_COMM, KGPTOT,&
  & ZGSSH, ZGMLD, ZG20D, ZGSSS, ZGTEM300, ZGSAL300 )

DO JSTGLO=1,KGPTOT,KPROMA
  IEND=MIN(KPROMA,KGPTOT-JSTGLO+1)
  IST=1
  IBL=(JSTGLO-1)/KPROMA+1
  DO JROF =1,IEND
    IF ((SD_VF(JROF,YSD_VF%YLSM%MP,IBL)>0.5_JPRB).OR.&
       &(SD_VF(JROF,YSD_VF%YCLK%MP,IBL)>0.5_JPRB)) THEN
      SD_OC(JROF,YSD_OC%YSSH%MP,IBL)  = ZMDI
      SD_OC(JROF,YSD_OC%YMLD%MP,IBL)  = ZMDI
      SD_OC(JROF,YSD_OC%Y20D%MP,IBL)  = ZMDI
      SD_OC(JROF,YSD_OC%YSSS%MP,IBL)  = ZMDI
      SD_OC(JROF,YSD_OC%YTEM3%MP,IBL) = ZMDI
      SD_OC(JROF,YSD_OC%YSAL3%MP,IBL) = ZMDI 
      SD_OC(JROF,YSD_OC%YICTH%MP,IBL) = ZMDI
    ELSE
      SD_OC(JROF,YSD_OC%YSSH%MP,IBL) = ZGSSH(JSTGLO+JROF-1)
      SD_OC(JROF,YSD_OC%YMLD%MP,IBL) = ZGMLD(JSTGLO+JROF-1)
      SD_OC(JROF,YSD_OC%Y20D%MP,IBL) = ZG20D(JSTGLO+JROF-1)
      SD_OC(JROF,YSD_OC%YSSS%MP,IBL) = ZGSSS(JSTGLO+JROF-1)
      IF (ZGSAL300(JSTGLO+JROF-1)<0.0_JPRD) THEN
        SD_OC(JROF,YSD_OC%YTEM3%MP,IBL) = ZMDI
        SD_OC(JROF,YSD_OC%YSAL3%MP,IBL) = ZMDI
      ELSE
        SD_OC(JROF,YSD_OC%YTEM3%MP,IBL) = ZGTEM300(JSTGLO+JROF-1)
        SD_OC(JROF,YSD_OC%YSAL3%MP,IBL) = ZGSAL300(JSTGLO+JROF-1)
      ENDIF
      IF (YDMCC%LNEMOLIMGET) THEN
        SD_OC(JROF,YSD_OC%YICTH%MP,IBL) = &
          & CPLNG_FLD(YDMCC%IP_A_ICE_THICKNESS)%D(JSTGLO+JROF-1,1,1) 
      ELSE
        SD_OC(JROF,YSD_OC%YICTH%MP,IBL) = ZMDI
      ENDIF
    ENDIF
  ENDDO
ENDDO

#else

CALL ABOR1('ININEMO: COMPILED WITHOUT WITH_NEMO')

#endif

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GETNEMODIAG',1,ZHOOK_HANDLE)

!     -----------------------------------------------------------

END SUBROUTINE GETNEMODIAG
