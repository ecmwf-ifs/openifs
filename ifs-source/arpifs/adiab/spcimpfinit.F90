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

SUBROUTINE SPCIMPFINIT(YDGEOMETRY,YDRIP,YDDYN,KM,KSTA,KEND,PSPVOR,PSPAUX)

!**** *SPCIMPFINIT* - PREPARATION FOR IMPLICIT CORIOLIS CASE.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCIMFINIT(..)

!        Explicit arguments :
!        --------------------  KM     - Zonal wavenumber        (in)
!                              KSTA   - First column processed  (in)
!                              KEND   - Last column processed   (in)
!                              PSPVOR - Distributed vorticity   (inout)
!                              PSPAUX - Correction to ZSPDIV    (out)

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Tomas Wilhelmsson *ECMWF*

!     Modifications.
!     --------------
!        22-09-2009 Original with extracted code from spcsi.F90
!        K. Yessad (Feb 2012): test on LLDOSI in the caller; simplifications.
!        K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : ROMEGA
USE YOMDYN       , ONLY : TDYN
USE YOMRIP       , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KM
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVOR(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPAUX(:,:) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZALPHA (0:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZDENIM (0:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZEPSI  (0:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZFPLUS (0:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZFMINUS(0:YDGEOMETRY%YRDIM%NSMAX+1)

INTEGER(KIND=JPIM) :: ILO, IMSP, IN, IRSP, JL, JLEV, JN 

REAL(KIND=JPRB) :: ZAL, ZBDT, ZEM, ZEN, ZF, ZTEMP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCIMPFINIT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NFLEVL=>YDDIMV%NFLEVL, &
 & RBTS2=>YDDYN%RBTS2, &
 & TDT=>YDRIP%TDT)
!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

!*        2.1  Preliminary initialisations.

ZBDT=RBTS2*TDT

!*        2.2  Set up helper arrays for implicit Coriolis case.

ZEM=REAL(KM,JPRB)
ZAL=2.0_JPRB*ZBDT*ROMEGA*ZEM
ILO=KM
IF (KM == 0) THEN
  ZALPHA(0)=0.0_JPRB
  ZDENIM(0)=0.0_JPRB
  ZEPSI(0)=0.0_JPRB
  ILO=1
ENDIF
DO JL=ILO,NSMAX
  ZEN=REAL(JL,JPRB)
  ZALPHA(JL)=ZAL/(ZEN*(ZEN+1.0_JPRB))
  ZDENIM(JL)=1.0_JPRB/(1.0_JPRB+ZALPHA(JL)**2)
  ZEPSI(JL)=SQRT((ZEN*ZEN-ZEM*ZEM)/(4.0_JPRB*ZEN*ZEN-1.0_JPRB))
ENDDO
ZALPHA(NSMAX+1)=0.0_JPRB
ZDENIM(NSMAX+1)=0.0_JPRB

IF (KM == 0) THEN
  ZFPLUS(0)=0.0_JPRB
  ZFMINUS(0)=0.0_JPRB
ENDIF
ZF=2.0_JPRB*ZBDT*ROMEGA
DO JL=ILO,NSMAX-1
  ZEN=REAL(JL,JPRB)
  ZFPLUS(JL)=ZF*ZEN*ZEPSI(JL+1)/(ZEN+1.0_JPRB)
  ZFMINUS(JL)=ZF*(ZEN+1.0_JPRB)*ZEPSI(JL)/ZEN
ENDDO
ZEN=REAL(NSMAX,JPRB)
ZFPLUS(NSMAX)=0.0_JPRB
ZFMINUS(NSMAX)=ZF*(ZEN+1.0_JPRB)*ZEPSI(NSMAX)/ZEN
ZFPLUS(NSMAX+1)=0.0_JPRB
ZFMINUS(NSMAX+1)=0.0_JPRB

!      Multiply rhs of vorticity equn by INV[J]

IF (KM > 0) THEN
  DO JN=KM,NSMAX
    DO JLEV=1,NFLEVL
      IRSP=KSTA+(JN-KM)*2
      IMSP=IRSP+1
      ZTEMP=ZDENIM(JN)*(PSPVOR(JLEV,IMSP)+ZALPHA(JN)*PSPVOR(JLEV,IRSP))  
      PSPVOR(JLEV,IRSP)=ZDENIM(JN)*(PSPVOR(JLEV,IRSP)-ZALPHA(JN)*PSPVOR(JLEV,IMSP))  
      PSPVOR(JLEV,IMSP)=ZTEMP
    ENDDO
  ENDDO
ENDIF
  
DO JN=KM,NSMAX
  DO JLEV=1,NFLEVL
    IRSP=KSTA+(JN-KM)*2
    IMSP=IRSP+1
    PSPAUX(JLEV,IRSP) = 0.0_JPRB
    PSPAUX(JLEV,IMSP) = 0.0_JPRB
  ENDDO
ENDDO
  
!      Add [F] * result to rhs of Helmholtz equation
  
DO JN=KM+1,NSMAX
  DO JLEV=1,NFLEVL
    IRSP=KSTA+(JN-KM)*2
    IMSP=IRSP+1
    PSPAUX(JLEV,IRSP)=PSPAUX(JLEV,IRSP)+ZFMINUS(JN)*PSPVOR(JLEV,IRSP-2)  
    PSPAUX(JLEV,IMSP)=PSPAUX(JLEV,IMSP)+ZFMINUS(JN)*PSPVOR(JLEV,IMSP-2)  
  ENDDO
ENDDO

DO JN=KM,NSMAX-1
  DO JLEV=1,NFLEVL
    IRSP=KSTA+(JN-KM)*2
    IMSP=IRSP+1
    PSPAUX(JLEV,IRSP)=PSPAUX(JLEV,IRSP)+ZFPLUS(JN)*PSPVOR(JLEV,IRSP+2)
    PSPAUX(JLEV,IMSP)=PSPAUX(JLEV,IMSP)+ZFPLUS(JN)*PSPVOR(JLEV,IMSP+2)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCIMPFINIT',1,ZHOOK_HANDLE)
END SUBROUTINE SPCIMPFINIT
