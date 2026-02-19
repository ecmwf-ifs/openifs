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

SUBROUTINE SUHEG(YDGEOMETRY,YDRIP,YDDYN,YDDYNA)

!**** *SUHEG*   - Initialize the solver of the Helmholz equation

!     Purpose.    Initialize the solver of the Helmholz equation
!     --------    in the semi-implicit, in case of stretching and
!                 use of geographical divergence (LSIDG=.T.).
!                 Cases where the unknown is DIV in the Helmholz equation.

!                 This code is valid only for a triangular truncation.
!                 (the condition "n is in [m,NSMAX]" is explicitly used).

!**   Interface.
!     ----------
!        *CALL* *SUHEG

!        Explicit arguments :
!        --------------------

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
!      K. YESSAD (OCTOBER 1993)

!     Modifications.
!     --------------
!      K. Yessad and A. Mary (Feb 2012): new SPGEOM_MOD; optimisations.
!      K. Yessad (July 2014): Move some variables.
!      K. Yessad (Dec 2016): Prune obsolete options.
!      K. Yessad (June 2017): Introduce NHQE model.
!      K. Yessad (June 2018): Alternate NHEE SI scheme elimination.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : TDYNA
USE YOMLUN   , ONLY : NULOUT
USE YOMRIP   , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT):: YDDYN
TYPE(TDYNA)    ,INTENT(INOUT):: YDDYNA
TYPE(TRIP)     ,INTENT(INOUT):: YDRIP

REAL(KIND=JPRB) :: ZD(YDGEOMETRY%YRDIMV%NFLEVG,0:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZE(YDGEOMETRY%YRDIMV%NFLEVG,0:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZF(YDGEOMETRY%YRDIMV%NFLEVG,0:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZES(YDGEOMETRY%YRDIMV%NFLEVG,0:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZFS(YDGEOMETRY%YRDIMV%NFLEVG,0:YDGEOMETRY%YRDIM%NSMAX)

INTEGER(KIND=JPIM) :: ILUX, IS0, ISE, JL, JN, JSE, IM, JMLOC
LOGICAL :: LLCOMPSIHEG

REAL(KIND=JPRB) :: ZBDT2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "suher.h"
#include "suhes.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHEG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDLAP=>YDGEOMETRY%YRLAP, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NSPEC=>YDDIM%NSPEC, NUMP=>YDDIM%NUMP, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LSIDG=>YDDYN%LSIDG, RBTS2=>YDDYN%RBTS2, SIHEG=>YDDYN%SIHEG, &
 & SIHEG2=>YDDYN%SIHEG2, SITIME=>YDDYN%SITIME, SIVP=>YDDYN%SIVP, &
 & MYMS=>YDLAP%MYMS, NSE0L=>YDLAP%NSE0L, RLAPDI=>YDLAP%RLAPDI, &
 & RLAPIN=>YDLAP%RLAPIN, &
 & TDT=>YDRIP%TDT, &
 & SCGMAP=>YDSPGEOM%SCGMAP)
!      ----------------------------------------------------------------

IF (.NOT.YDDYNA%LNHDYN) THEN
  IF (SITIME /= TDT) THEN
    SITIME=TDT
    LLCOMPSIHEG=.TRUE.
    WRITE(NULOUT,*) ' SUHEG: COMPUTES HELMHOLTZ EQUATION OPERATOR'&
     & ,' IN CASE OF SEMI-IMPLICIT SCHEME IN NOT'&
     & ,' REDUCED DIVERGENCE (LSIDG=.T.): HYDROSTATIC MODEL'  
  ELSE
    LLCOMPSIHEG=.FALSE.
  ENDIF
ELSEIF (YDDYNA%LNHEE.AND.(.NOT.YDDYNA%LSI_NHEE)) THEN
  ! Contrary to the hydrostatic case, the test on SITIME
  ! and the resetting of SITIME is already done in SUNHSI, so it has not to be redone there.
  ! For the time being, SUHEG should not be called in this case, but that may change in the future.
  LLCOMPSIHEG=.FALSE.
ELSEIF (YDDYNA%LNHEE.AND.YDDYNA%LSI_NHEE) THEN
  ! Contrary to the hydrostatic case, the test on SITIME
  ! and the resetting of SITIME is already done in SUNHEESI, so it has not to be redone there.
  LLCOMPSIHEG=.TRUE.
  WRITE(NULOUT,*) ' SUHEG: COMPUTES HELMHOLTZ EQUATION OPERATOR'&
   & ,' IN CASE OF SEMI-IMPLICIT SCHEME IN NOT'&
   & ,' REDUCED DIVERGENCE (LSIDG=.T.): NHEE MODEL'
ELSEIF (YDDYNA%LNHQE) THEN
  ! Contrary to the hydrostatic case, the test on SITIME
  ! and the resetting of SITIME is already done in SUNHQESI, so it has not to be redone there.
  LLCOMPSIHEG=.TRUE.
  WRITE(NULOUT,*) ' SUHEG: COMPUTES HELMHOLTZ EQUATION OPERATOR'&
   & ,' IN CASE OF SEMI-IMPLICIT SCHEME IN NOT'&
   & ,' REDUCED DIVERGENCE (LSIDG=.T.): NHQE MODEL'  
ENDIF

IF(LLCOMPSIHEG)THEN

!*       1.    INITIALIZE HELMHOLTZ EQUATION SOLVER.
!              -------------------------------------

  ILUX=NSPEC
  ZBDT2=(RBTS2*SITIME)**2

  DO JL=1,NFLEVG
    DO JSE=1,ILUX
      SIHEG(JL,JSE,1)=0.0_JPRB
      SIHEG(JL,JSE,2)=0.0_JPRB
      SIHEG(JL,JSE,3)=0.0_JPRB
    ENDDO
  ENDDO

  DO JMLOC=1,NUMP

    IM=MYMS(JMLOC)

    IF (IM == 0) THEN

      DO JL=1,NFLEVG
        DO JSE=1,NSMAX+1
          SIHEG2(JL,JSE,2)=0.0_JPRB
          SIHEG2(JL,JSE,3)=0.0_JPRB
        ENDDO
      ENDDO

      IS0=NSE0L(JMLOC)
      DO JL=1,NFLEVG
        DO JN=0,NSMAX
          ISE=IS0+JN+1
          ZD(JL,JN)=1.0_JPRB-ZBDT2*SIVP(JL)*RLAPDI(JN)*SCGMAP(ISE,1)
        ENDDO
        DO JN=0,NSMAX-1
          ISE=IS0+JN+1
          ZE(JL,JN)=-ZBDT2*SIVP(JL)*RLAPDI(JN+1)*SCGMAP(ISE,2)
          ZES(JL,JN)=-ZBDT2*SIVP(JL)*RLAPDI(JN)*SCGMAP(ISE,2)
        ENDDO
        DO JN=0,NSMAX-2
          ISE=IS0+JN+1
          ZF(JL,JN)=-ZBDT2*SIVP(JL)*RLAPDI(JN+2)*SCGMAP(ISE,3)
          ZFS(JL,JN)=-ZBDT2*SIVP(JL)*RLAPDI(JN)*SCGMAP(ISE,3)
        ENDDO
      ENDDO
      CALL SUHER(NSMAX+1,NFLEVG,NFLEVG,ZD(1,0),ZE(1,0),ZES(1,0),&
       & ZF(1,0),ZFS(1,0),SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),&
       & SIHEG(1,IS0+1,3),SIHEG2(1,1,2),SIHEG2(1,1,3))  

    ELSE

      IS0=NSE0L(JMLOC)
      DO JL=1,NFLEVG
        DO JN=IM,NSMAX
          ISE=IS0+JN+1-IM
          ZD(JL,JN)=RLAPIN(JN)-ZBDT2*SIVP(JL)*SCGMAP(ISE,1)
          ZE(JL,JN)=-ZBDT2*SIVP(JL)*SCGMAP(ISE,2)
          ZF(JL,JN)=-ZBDT2*SIVP(JL)*SCGMAP(ISE,3)
        ENDDO
      ENDDO
      CALL SUHES(NSMAX+1-IM,NFLEVG,NFLEVG,ZD(1,IM),ZE(1,IM),ZF(1,IM),&
       & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3))  

    ENDIF

  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHEG',1,ZHOOK_HANDLE)
END SUBROUTINE SUHEG
