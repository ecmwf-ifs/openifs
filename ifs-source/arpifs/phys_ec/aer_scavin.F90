! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SCAVIN &
 & ( YDEAERSNK,KIDIA, KFDIA, KLON , KLEV , PSCAVFR,KSTART, KSTEP, &
 &   PDP   , PDZ    , PRHO , PTP, PFLXR, PFLXS, PCLCOV, PCLWAT, PWATER, PICE, PAERO, PTAERI, &
 &   PTSPHY, &
 &   PTAERO, PFAERO &
 & )

!*** * AER_SCAVIN* - IN-CLOUD SCAVENGING OF TRACERS


! N.B.: only SEA-SALT and DUST are presently catered for.
! =======================================================


!**   INTERFACE.
!     ----------
!          *AER_SCAVIN* IS CALLED FROM *CALLPAR*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        FROM O. BOUCHER (LOA, 1998-03) in_scav

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2004-05-07
!        P.BECHTOLD 2009-04-22 Clean/optimise and rewrite to externalise 
!                              scavenging coefficient

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    ,ONLY : RG

USE YOEAERSNK ,ONLY : TEAERSNK

!USE YOEDBUG   ,ONLY : KSTPDBG, NSTPDBG

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

TYPE(TEAERSNK)    ,INTENT(INOUT):: YDEAERSNK
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV, KSTART, KSTEP

REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV)    , PDZ(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PCLCOV(KLON,KLEV) , PCLWAT(KLON,KLEV)  , PRHO(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PWATER(KLON,KLEV), PICE(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PFLXR(KLON,KLEV+1), PFLXS(KLON,KLEV+1)
REAL(KIND=JPRB),INTENT(IN)    :: PAERO(KLON,KLEV)  , PTAERI(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB),INTENT(IN)    :: PSCAVFR

REAL(KIND=JPRB),INTENT(OUT)   :: PFAERO(KLON), PTAERO(KLON,KLEV)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZAERI(KLON,KLEV), ZFAERO(KLON)
REAL(KIND=JPRB) :: ZBETAI_W, ZBETAK_W, ZBETAR, ZTMP
REAL(KIND=JPRB) :: ZBETAI_I, ZBETAK_I, ZBETA
REAL(KIND=JPRB) :: ZDP, ZDX, ZDZ, ZFRAC, ZFUNC, ZRHO
REAL(KIND=JPRB) :: ZEPSQLIQ, ZEPSFLX, ZTSPHY, ZRG, ZSCAVFR
LOGICAL :: LLPRINT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SCAVIN',0,ZHOOK_HANDLE)
ASSOCIATE(RFRAER=>YDEAERSNK%RFRAER)
!- PAERO  in kg cm-3
!  PDZ    in m
!  PFLXx  in kg m-2 s-1
!  PCLWAT  in kg kg-1
!  PRHO   in kg m-3
!  PSCAVFR a fraction N.D.
!  PTAERI in kg m-3 s-1
!  PTAERO in KG m-3 s-1


! NB: PCLWAT is the in-cloud water mixing ratio

ZRG=1.0_JPRB/RG
ZTSPHY=1.0_JPRB/PTSPHY

LLPRINT=.FALSE.
!DO JL=1,NSTPDBG
!  IF (KSTEP == KSTPDBG(JL)) THEN
!    LLPRINT=.TRUE.
!  ENDIF
!ENDDO

ZEPSFLX =1.E-18_JPRB
ZEPSQLIQ=1.E-18_JPRB
ZFRAC=RFRAER

DO JL=KIDIA,KFDIA
  PFAERO(JL)=0._JPRB
  ZFAERO(JL)=0._JPRB
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
!-- initialisation

    ZAERI(JL,JK) = PAERO(JL,JK) + PTSPHY * PTAERI(JL,JK)
     
    ZRHO=PRHO(JL,JK)                          ! in kg m-3
    ZDZ =PDZ(JL,JK)                           ! in m
    ZDP =PDP(JL,JK)                           ! in Pa (kg m-1 s-2) 

!-- scavenging
    ZBETAI_W=PFLXR(JL,JK+1)-PFLXR(JL,JK)   ! in kg m-2 s-1 (w/o cloud fraction effect)
    ZBETAI_I=PFLXS(JL,JK+1)-PFLXS(JL,JK)   ! in kg m-2 s-1 (w/o cloud fraction effect)

    ZBETAK_W=( PFLXR(JL,JK+1)-PFLXR(JL,JK) ) / MAX (0.01_JPRB, PCLCOV(JL,JK) )  
    ZBETAK_I=( PFLXS(JL,JK+1)-PFLXS(JL,JK) ) / MAX (0.01_JPRB, PCLCOV(JL,JK) )  
        
!  liquid contribution                                              ! in kg m-2 s-1
    ZBETA=ZBETAK_W

!-- security added
    IF (PWATER(JL,JK)+PICE(JL,JK) > ZEPSQLIQ .AND. PSCAVFR > 0._JPRB) THEN
      ZSCAVFR=0.06_JPRB*(PICE(JL,JK)/(PICE(JL,JK)+PWATER(JL,JK)))+PSCAVFR*(1._JPRB-(PICE(JL,JK)/(PICE(JL,JK)+PWATER(JL,JK))))
    ELSE
      ZSCAVFR=0._JPRB
    ENDIF

    ZSCAVFR=PSCAVFR

    IF (PWATER(JL,JK) > ZEPSQLIQ ) THEN
      ZBETA=ZBETA/(ZDP*PWATER(JL,JK)*ZRG) ! in s-1
    ELSE
      ZBETA=0._JPRB
    ENDIF

    ZBETA=MIN( 200._JPRB, MAX( ZBETA, 0._JPRB ))
    ZFUNC=EXP(-ZSCAVFR*ZBETA*PTSPHY)                           ! N.D.       (N.D.)

    ZDX = ZAERI(JL,JK)*(ZFUNC - 1._JPRB) * PCLCOV(JL,JK)      ! in kg kg-1

    PTAERO(JL,JK) = PTAERI(JL,JK) + ZDX*ZTSPHY                ! in kg kg-1 s-1

!-- factor ZMASSE/RNAVO =1.

    PFAERO(JL) = PFAERO(JL) - ZDX*ZTSPHY * ZDP*ZRG             ! in kg m-2 s-1
    
!-- reevaporation  (NB: this uses ZBETAI, i.e., independent of the cloud fraction)
!-- security on fluxes added
    IF (ZBETAI_W < 0._JPRB .AND. PFLXR(JL,JK) > ZEPSFLX) THEN
      ZBETAR=ZBETAI_W/(PFLXR(JL,JK))
    ELSE
      ZBETAR=0._JPRB
    ENDIF

    IF ( PFLXR(JL,JK+1) <= ZEPSFLX) THEN
!-- total reevaporation
      ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)              ! N.D.
    ELSE
!-- non total reevaporation for the aerosols
      ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)        ! N.D.
      ZTMP=ZBETAR
      ZBETAR=ZBETAR*(1._JPRB-EXP(-2._JPRB*ZBETAR**0.5_JPRB)* &
             & (1._JPRB + 2._JPRB*ZBETAR**0.5_JPRB + 2._JPRB*ZBETAR + 1.33_JPRB*ZBETAR**1.5_JPRB))* &
             & (1._JPRB-ZBETAR) + ZBETAR**2._JPRB
      !write(*,*) "REEVAP",JL,JK,ZBETAR,ZTMP
      
    ENDIF
!-- factor ZMASSE/RNAVO =1.

    ZDX = ZBETAR*PFAERO(JL) * RG/ZDP 
    PTAERO(JL,JK) = PTAERO(JL,JK) + ZDX
    PFAERO(JL) = (1._JPRB-ZBETAR) * PFAERO(JL)



!  solid contribution                                              ! in kg m-2 s-1
    ZBETA=ZBETAK_I
    ZSCAVFR=0.06_JPRB

!-- security added
    IF (PICE(JL,JK) > ZEPSQLIQ .AND. PSCAVFR > 0._JPRB) THEN
      ZBETA=ZBETA/(ZDP*PICE(JL,JK)*ZRG) ! in s-1
    ELSE
      ZBETA=0._JPRB
    ENDIF

    ZBETA=MIN( 200._JPRB, MAX( ZBETA, 0._JPRB ))
! modulate scav efficiency following Bourgeois and Bey (2011)
    ZFUNC=EXP(-ZSCAVFR*ZBETA*PTSPHY)                           ! N.D.       (N.D.)

    ZDX = ZAERI(JL,JK)*(ZFUNC - 1._JPRB) * PCLCOV(JL,JK)      ! in kg kg-1

    PTAERO(JL,JK) = PTAERO(JL,JK) + ZDX*ZTSPHY                ! in kg kg-1 s-1

!-- factor ZMASSE/RNAVO =1.

    ZFAERO(JL) = ZFAERO(JL) - ZDX*ZTSPHY * ZDP*ZRG             ! in kg m-2 s-1
    
!-- reevaporation  (NB: this uses ZBETAI, i.e., independent of the cloud fraction)
!-- security on fluxes added
    IF (ZBETAI_I < 0._JPRB .AND. PFLXS(JL,JK) > ZEPSFLX) THEN
      ZBETAR=ZBETAI_I/(PFLXS(JL,JK))
    ELSE
      ZBETAR=0._JPRB
    ENDIF

    IF ( PFLXS(JL,JK+1) <= ZEPSFLX) THEN
!-- total reevaporation
      ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)              ! N.D.
    ELSE
!-- non total reevaporation for the aerosols
      ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)        ! N.D.
      ZTMP=ZBETAR
      !ZBETAR=ZBETAR*RFRAER
      ZBETAR=ZBETAR*(1._JPRB-EXP(-2._JPRB*ZBETAR**0.5_JPRB)* &
             & (1._JPRB + 2._JPRB*ZBETAR**0.5_JPRB + 2._JPRB*ZBETAR + 1.33_JPRB*ZBETAR**1.5_JPRB))* &
             & (1._JPRB-ZBETAR) + ZBETAR**2._JPRB
      !write(*,*) "REEVAP",JL,JK,ZBETAR,ZTMP
      
    ENDIF
!-- factor ZMASSE/RNAVO =1.

    ZDX = ZBETAR*ZFAERO(JL) * RG/ZDP 
    PTAERO(JL,JK) = PTAERO(JL,JK) + ZDX
    ZFAERO(JL) = (1._JPRB-ZBETAR) * ZFAERO(JL)

  ENDDO
ENDDO
PFAERO(KIDIA:KFDIA)=PFAERO(KIDIA:KFDIA)+ZFAERO(KIDIA:KFDIA)

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SCAVIN',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SCAVIN
