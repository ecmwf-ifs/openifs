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

SUBROUTINE SPFILT(YDGEOMETRY,KSTA,KEND,PSPVOR,PSPDIV)
!     ------------------------------------------------------------------

!**** *SPFILT* - HORIZONTAL SPECTRAL SPACE FILTERING

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPFILT(..)
!     Author.
!     -------
!     Nils Wedi (2011)
!
!     Modifications
!     -------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE


TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVOR(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIV(:,:) 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZFAC, ZN
INTEGER(KIND=JPIM) :: JLEV, JSP, IN, IGG, ISHIFT
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPFILT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NFLEVL=>YDDIMV%NFLEVL)

! dealiasing adjustments
IGG=NINT(REAL(2*NSMAX,JPRB)/3._JPRB)-1
ISHIFT = NINT(0.83_JPRB*REAL(NSMAX,JPRB))-IGG

!$OMP PARALLEL DO SCHEDULE(STATIC,1) &
!$OMP& PRIVATE(JLEV,JSP,IN,ZFAC,ZN)
DO JSP=KSTA,KEND
  
  IN = YDLAP%NVALUE(JSP) + ISHIFT
  IF( IN > IGG ) THEN
    ! ZFAC = 8._JPRB*(EXP(-0.5_JPRB*REAL(IN-NSMAX,JPRB)**2/REAL(IN-IGG,JPRB)**2))
    ZFAC = 32._JPRB*(EXP(-0.5_JPRB*REAL(IN-NSMAX,JPRB)**2/REAL(IN-IGG,JPRB)**2))
    ! ZN = 1._JPRB/(1._JPRB+(REAL(IN,JPRB)/REAL(IGG,JPRB))**IEXP)
  ELSE
    ZFAC = 0._JPRB
  ENDIF
  ZN = 1._JPRB/(1._JPRB+ZFAC*ZFAC)
  
  DO JLEV=1,NFLEVL
    PSPVOR(JLEV,JSP)=ZN*PSPVOR(JLEV,JSP)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPFILT',1,ZHOOK_HANDLE)

END SUBROUTINE SPFILT
