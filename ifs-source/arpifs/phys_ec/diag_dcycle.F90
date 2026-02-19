! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DIAG_DCYCLE &
 & ( KIDIA,  KFDIA,  KLON,   KLEV,  &
 &   KLEVX,  KFLDX,  PTSTEP,&
 &   PRS1,   PFPLCL, PFPLCN, PFPLSL, PFPLSN,&
 &   PDIFTS, PDIFTQ, PFRTH,  PFRTHC,&
 &   PL,     PI,     PTL,    PTCFL,  PQCFL,&
 &   PEXTRA)

!          A.BELJAARS+P.BECHTOLD     E.C.M.W.F.     05/2009

!          PURPOSE
!          -------

!          TO PROVIDE DIAGNOSTICS FOR A CERTAIN NUMBER OF FIELDS
!          IN HOURLY RESOLUTION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CALLPAR*.

!          METHOD.
!          --------
!          DATA IS STORED IN MODEL-LEVEL EXTRA FIELD
!          LEVELS CORRESPOND TO HOURS

!     PARAMETER     DESCRIPTION                                   UNITS 
!     ---------     -----------                                   ----- 
!     INPUT PARAMETERS (INTEGER): 

!    see decsription in code below

!    OUTPUT PARAMETERS (INTEGER):

!    *PEXTRA*      3 3DEXTRA FIELDS 

!          MODIFICATIONS
!          -------------

!----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG
USE YOMLUN   , ONLY : NULOUT
USE YOMCT3   , ONLY : NSTEP

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVX
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX

REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRS1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCN(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSN(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIFTS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIFTQ(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTHC(KLON,0:1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTCFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCFL(KLON) 

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX) 

REAL(KIND=JPRB) ::     ZFIELD2D(KLON)
REAL(KIND=JPRB) ::     ZDP, ZRG, ZHSTEP, ZH

INTEGER(KIND=JPIM) :: JL, JK, JK1, JK2, JK3, ID, IH

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIAG_DCYCLE',0,ZHOOK_HANDLE)
!----------------------------------------------------------------------

IF (KLEV<74) THEN
  WRITE(NULOUT,*)'WARNING!! ROUTINE DIAG_CYCLE ONLY WORKS for KLEV>=74 LEVELS'
  WRITE(NULOUT,*)'IT WILL NOT BE EXECUTED, SO NO DIURNAL CYCLE FIELDS'
ELSE
  ZRG=1.0_JPRB/RG
  ZHSTEP=NSTEP*PTSTEP/3600.
  
  ID=ZHSTEP/24
  ZH=ZHSTEP-ID*24
  IH=NINT(ZH)
  JK1= 1+IH
  JK2=25+IH
  JK3=50+IH
  
  IF(NSTEP==0) THEN
    PEXTRA(:,:,:)=0.0_JPRB
  ENDIF
  ZFIELD2D(:)=0.0_JPRB
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZDP=(PRS1(JL,JK)-PRS1(JL,JK-1))*ZRG
      ZFIELD2D(JL)=ZFIELD2D(JL)+(PL(JL,JK)+PI(JL,JK))*ZDP
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
  ! surface precip
    PEXTRA(JL,JK1,1)=PEXTRA(JL,JK1,1)&
                    &+( PFPLCL(JL,KLEV)+PFPLCN(JL,KLEV)&
                    &  +PFPLSL(JL,KLEV)+PFPLSN(JL,KLEV) )*PTSTEP
  ! surface sensible heat flux
    PEXTRA(JL,JK2,1)=PEXTRA(JL,JK2,1)+PDIFTS(JL,KLEV)*PTSTEP
  ! surface moisture flux
    PEXTRA(JL,JK3,1)=PEXTRA(JL,JK3,1)+PDIFTQ(JL,KLEV)*PTSTEP
  
  ! OLR
    PEXTRA(JL,JK1,2)=PEXTRA(JL,JK1,2)+PFRTH(JL,0)*PTSTEP
  ! cloudy OLR=OLR-clear sky OLR
    PEXTRA(JL,JK2,2)=PEXTRA(JL,JK2,2)+(PFRTH(JL,0)-PFRTHC(JL,0))*PTSTEP
  ! total column water
    PEXTRA(JL,JK3,2)=PEXTRA(JL,JK3,2)+ZFIELD2D(JL)*PTSTEP
  
  ! skin temperature
    PEXTRA(JL,JK1,3)=PEXTRA(JL,JK1,3)+PTL  (JL)*PTSTEP
  ! 2m temperature
    PEXTRA(JL,JK2,3)=PEXTRA(JL,JK2,3)+PTCFL(JL)*PTSTEP
  ! 2m specific humidity
    PEXTRA(JL,JK3,3)=PEXTRA(JL,JK3,3)+PQCFL(JL)*PTSTEP
  ENDDO
ENDIF
!-----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIAG_DCYCLE',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------------------
END SUBROUTINE DIAG_DCYCLE
