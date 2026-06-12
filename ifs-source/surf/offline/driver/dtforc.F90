SUBROUTINE DTFORC

USE PARKIND1  ,ONLY : JPIM     ,JPRB , JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMFORC1S, ONLY : UFI      ,VFI      ,TFI      ,QFI      ,&
            &PSFI     ,SRFFI    ,TRFFI    ,R30FI    ,S30FI    ,&
            &R30FI_C  ,S30FI_C  ,&
            &DTIMFC   ,RTSTFC   ,NSTPFC, CO2FI
USE YOMRIP   , ONLY : RTIMTR   ,RTIMST, NSSSSS
USE YOMLUN1S , ONLY : NULOUT
USE YOMDYN1S , ONLY : NSTEP    ,TSTEP    ,NACCTYPE,LPREINT,LSWINT,LFLXINT
USE YOMGF1S  , ONLY : UNLEV0   ,VNLEV0   ,TNLEV0   ,QNLEV0   ,&
           &PNLP0    ,UNLEV1   ,VNLEV1   ,TNLEV1   ,QNLEV1   , PNLP1,&
           &FSSRD    ,FSTRD    ,FLSRF    ,FCRF     ,FLSSF    ,FCSF, &
           &CNLEV0, CNLEV1
USE YOMLOG1S , ONLY : LDBGS1,IDBGS1
USE YOMGC1S  , ONLY : GELAM   ,GELAT, GEMU
USE YOERIP   , ONLY : RCODECM  ,RSIDECM  ,RCOVSRM  ,RSIVSRM
USE YOMCST   , ONLY : RDAY, REA , REPSM, RMD, RMCO2, RANCO2, RCO2REFYEAR, RNCO2YEARS
USE YOMDPHY  , ONLY : NPOI, NTRAC
USE YOEPHY   , ONLY : LEAIRCO2COUP


#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *DTFORC *  - TIME INTERPOLATION OF THE ATMOSPHERIC FORCING DATA

!     PURPOSE.
!     --------
!        Computes the atmospheric forcing for the current time step

!**   INTERFACE.
!     ----------
!        *CALL* *DTFORC*

!-----------------------------------------------------------------------

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------
!        IYMD2C

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE 
!        ONE COLUMN SURFACE MODEL

!     AUTHOR.
!     -------
!        JEAN-FRANCOIS MAHFOUF AND PEDRO VITERBO  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 95-03-13
!        BART VD HURK: MULTIPLE GRID POINTS (2000-07-13)
!        S.Boussetta: interpolation following the solar zenith angle for SW radiation (Jan 2013)
!        A. Agusti-Panareda: Atmospheric CO2 forcing (2020-11-17)
        
!     --------------------------

!
!     A FEW IMPORTANT PARAMETERS
!     --------------------------
!        RTIMST    Time of start of forecast (Julian seconds)
!        TSTEP     Time step of the forecast model (seconds)
!        NSTEP     Current model step (0 for first step)
!        RTIMTR    Midpoint time of current time step (Julian seconds)
!
!        RTSTFC    Time of start of forcing (Julian seconds)
!                  (For fluxes: This can be start, midpoint or end of 
!                   the first fluxinterval with lenth DTIMFC) 
!        DTIMFC    Time step of forcing (seconds)
!
!
!        NACCTYPE  Accum. type: 0  =  Centred flux
!                               1  =  Time stamp at start of flux interval
!                               2  =  Time stamp at end of flux interval
!
! Solar related parameters
!     RCODECM: COSINE OF THE DECLINATION
!     RSIDECM:   SINE OF THE DECLINATION
!     RCOVSRM: COSINE OF TRUE SOLAR TIME
!     RSIVSRM:   SINE OF TRUE SOLAR TIME
!
! Grid point related parameters
!     GEMU    -  SIN of latitude on the real earth
!     GELAM   -  longitude on the real earth
!     GELAT   -  latitude on the real earth



#endif
IMPLICIT NONE

!* LOCAL VARIABLES
REAL(KIND=JPRD) :: ZTIMCUR,ZW,ZWP1,ZWF,ZWFP1,ZWP,ZWPP1,TP1,TP2,TP3
REAL(KIND=JPRD) :: ZWa(6,3),ZWTMP(NPOI,6),ZWSUM(NPOI)
INTEGER(KIND=JPIM) :: IF,IFP1,IFF,IFF1,IFF2,IFF3,IFFP1,JL,IFPREC,IFP,IFPP1,IFSOLAR
REAL(KIND=JPRD) :: Z_SINLAT, Z_COSLAT, Z_PI, Z_TWOPI,Z_RADCON, Z_CONRAD, Z_RLHH, Z_COSLHH
REAL(KIND=JPRD) :: Z_SINE, Z_DECL, Z_ANGLE
REAL(KIND=JPRD) :: ZMU0M(NPOI),ZANG(NPOI),SUMZANG(NPOI,7),ZWS(NPOI)
REAL(KIND=JPRD) :: Z_ZENITH,Z_S_ELEV
REAL(KIND=JPRD) :: ZT1, ZT2, ZTETA1, ZTETA2, ZDECLIM1,ZDECLIM2, ZCODECM1,ZCODECM2,ZSIDECM1, ZSIDECM2
REAL(KIND=JPRD) :: ZHGMT1, ZHGMT2, ZEQTIMM1, ZEQTIMM2, ZSOVRM1, ZSOVRM2, ZWSOVRM1, ZWSOVRM2
REAL(KIND=JPRD) :: ZCOVSRM1, ZCOVSRM2, ZSIVSRM1, ZSIVSRM2, ZANG1, ZANG2, ZZT1, ZZT2

REAL(KIND=JPRD) :: ZT, ZTETA, ZDECLIM, ZCODECM,ZSIDECM
REAL(KIND=JPRD) :: ZHGMT, ZEQTIMM, ZSOVRM, ZWSOVRM
REAL(KIND=JPRD) :: ZCOVSRM, ZSIVSRM, ZZT, TT,  ZWSS
REAL(KIND=JPRD) :: Z_ZENITH1,Z_S_ELEV1, Z_ZENITH2,Z_S_ELEV2,Z_ANGLE1,Z_ANGLE2
REAL(KIND=JPRD) :: ZJUL,ZCO2


INTEGER(KIND=JPIM) :: NFORC, JK, JKK, STACT, IFFF,JJ
INTEGER(KIND=JPIM) :: IYMD,IHM,IYYYY,YYYY

INTEGER(KIND=JPIM) :: ITMP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fctast.h"

IF (LHOOK) CALL DR_HOOK('DTFORC',0,ZHOOK_HANDLE)

IF ( LSWINT ) THEN
    Z_PI=2.0_JPRD*ASIN(1.0_JPRD)
    Z_TWOPI=2.0_JPRD*Z_PI
    Z_RADCON=Z_TWOPI/360.0_JPRD
    Z_CONRAD=360.0_JPRD/Z_TWOPI
    ! Solar zenith angle computation
    ZMU0M(:)=MAX( RSIDECM*GEMU(:) &
        & -RCODECM*RCOVSRM*SQRT(1._JPRD-GEMU(:)**2)*COS(GELAM(:)) &
        & +RCODECM*RSIVSRM*SQRT(1._JPRD-GEMU(:)**2)*SIN(GELAM(:)) &
        & ,0._JPRD)
ENDIF 




! ztimcur is the seconds at the beginning of the time step 
ztimcur=RTIMST+real(NSTEP,KIND=JPRD)*TSTEP
IF=(ztimcur-RTSTFC)/DTIMFC+1
IFP1=IF+1
ZW=1._JPRD-(ztimcur-RTSTFC-REAL(IF-1,KIND=JPRD)*DTIMFC)/DTIMFC
ZWP1=1._JPRD-ZW
IF (IFP1 == NSTPFC+1) THEN
   IFP1=IF
   ZW=1._JPRD
   ZWP1=0._JPRD
ELSEIF (IFP1 > NSTPFC+1) THEN
   WRITE(*,*) 'out of bounds IFP1:',IFP1,NSTPFC+1
   CALl ABORT()
ENDIF



IF (NACCTYPE.eq.0) THEN
   ! Fluxed centred on the timestamp |---X---| CENTRED ACCUMULATION

   !! Linear interpolation for SWdown and LWdown
   IFF=(rtimtr-RTSTFC)/DTIMFC+1
   ZWf=1._JPRD-(rtimtr-RTSTFC-REAL(IFF-1,KIND=JPRD)*DTIMFC)/DTIMFC
   IFFP1=IFF+1
   ZWfP1=1._JPRD-ZWf

   !! no flx interpolatin 
   IF (.NOT. LFLXINT) THEN
       IF ( ZWf > 0.5 ) THEN 
          ZWf=1._JPRD
       ELSE
          ZWf=0._JPRD
       ENDIF
       ZWfP1=1._JPRD-ZWf
    ENDIF 

   
   IF (IFFP1 > NSTPFC) THEN
      IFFP1=IFF
      ZWf=1._JPRD
      ZWfP1=0._JPRD
   ENDIF

ENDIF
IF (NACCTYPE.eq.1) THEN
   ! Fluxed beginning at the timestamp |X------| FORWARD ACCUMULATION
   IFF=IF
   IFFP1=IFF+1
   ZWf=1._JPRD
   ZWfP1=0._JPRD
   IF (IFFP1 > NSTPFC) THEN
      IFFP1=IFF
   ENDIF
ENDIF
IF (NACCTYPE.eq.2) THEN
   ! Fluxed ending at the timestamp |------X| BACKWARD ACCUMULATION
   
   ! FLX Interpolation
   IFF=(rtimtr-RTSTFC+DTIMFC/2)/DTIMFC+1
   ZWf=1._JPRD-(rtimtr-RTSTFC+DTIMFC/2-REAL(IFF-1,KIND=JPRD)*DTIMFC)/DTIMFC
   ! NO flx INTERPOLATION 
   IF (.NOT. LFLXINT) THEN
      IFF=IF
      ZWf=0._JPRD
   ENDIF 
   IFFP1=IFF+1
   ZWfP1=1._JPRD-ZWf


    IF ( LSWINT ) THEN 
    
        !=====================================================================
        !Compute the interpoloation coefficients following the solar zenith angle
        !=====================================================================
        !start of the interpolation interval
        IFFF=(rtimtr-RTSTFC)/DTIMFC+1
        ZT1=(IFFF-1)*DTIMFC+RTSTFC

        ! number of steps within the forcing interval
        NFORC=NINT(DTIMFC/TSTEP)

        ! Determine which of the NFORC FRACTIONS is the ACTUAL ONE in the forecast period
        JK=1
        DO WHILE ( (ZT1+REAL((JK)*TSTEP,KIND=JPRD)) .LT. ZTIMCUR+TSTEP)
            JK=JK+1
        ENDDO
        STACT=JK
        IF ( IDBGS1 > 1 ) THEN
          WRITE(NULOUT,*) "STACT= ",STACT
        ENDIF

        ! Loop to compute all Solar zenith angle within the forcing interval
        ! We loop through all TSTEPs in the period
        ! For each TSTEP we compute one zenith angle average based on 4 values (every 1/3 including start and end)
        ! The actual zenith angle related radiation value is determined by
        ! (actual) TSTEP-Average/SUM-of-TSTEP-averages-in-the-forcing-period
        SUMZANG(:,:)=0._JPRD
        DO JK=1,NFORC
            ZANG(:)=0._JPRD
            DO JKK=1,4
        !         ZT=ZT1+REAL((JK)*TSTEP)
                ZT=ZT1+REAL((JK-1)*TSTEP+(JKK-1)*TSTEP/3,KIND=JPRD)
                ZTETA=RTETA(ZT)
                ZDECLIM=RDS(ZTETA)
                ZCODECM=COS(ZDECLIM)
                ZSIDECM=SIN(ZDECLIM)
        !         TT=((NSTEP-1)+JK)*TSTEP
        !         print *,"TT= ",TT
        !         ZZT=MOD(NINT(((IFFF-1)*NFORC+JK)*TSTEP),NINT(RDAY))
!                  ZZT=MOD(NINT((IFFF-1)*DTIMFC+(JK-1)*TSTEP+(JKK-1)*TSTEP/3),NINT(RDAY))
!                ZTT calculation assumes forcing starts at 00:00, if not we should add here that !!!
                  ZZT=MOD(NINT(ABS(RTSTFC)+(IFFF-1)*DTIMFC+(JK-1)*TSTEP+(JKK-1)*TSTEP/3),NINT(RDAY))
                 WRITE(NULOUT,"('ZT=',F14.2,' ZTT=',F14.2,' jk',I6,' jkk=',I6)") ZT,ZZT,JK,JKK
!                 print *,"ZT=", ZT,"ZZT=",ZZT,JK,JKK
                ZHGMT=REAL(ZZT,KIND=JPRD)
                ZEQTIMM=RET(ZTETA)
                ZSOVRM =ZEQTIMM+ZHGMT                                 
                ZWSOVRM=ZSOVRM*2._JPRD*Z_PI/RDAY
                ZCOVSRM=COS(ZWSOVRM)
                ZSIVSRM=SIN(ZWSOVRM)
                ZANG(:)=ZANG(:)+MAX( ZSIDECM*GEMU(:) &
                    & -ZCODECM*ZCOVSRM*SQRT(1._JPRD-GEMU(:)**2)*COS(GELAM(:)) &
                    & +ZCODECM*ZSIVSRM*SQRT(1._JPRD-GEMU(:)**2)*SIN(GELAM(:)) &
                    & ,0._JPRD)
            END DO
            ZANG(:)=ZANG(:)/4

            SUMZANG(:,JK+1)=ZANG(:)
            SUMZANG(:,1)=SUMZANG(:,1)+ZANG(:)
            !print *, "ZANG(JL)= ", ZANG, "ZT=", ZT

            !WRITE(*,*) NSTEP,JK
            !DO JJ=1000,20000,1000
            !   WRITE(*,'(2I3,I10,F15.0,14F12.3,3I10)') NSTEP,JK,JJ,ZT,ZT-ZT1,RSIDECM,RCOVSRM,RSIVSRM,ZTETA,&
            !                ZDECLIM,ZCODECM,ZSIDECM,ZHGMT,ZEQTIMM,ZSOVRM,ZWSOVRM,ZCOVSRM,ZSIVSRM,&
            !                NINT((IFFF-1)*DTIMFC+JK*TSTEP),NINT((IFFF-1)*DTIMFC+(JK-1)*TSTEP+(JKK-1)*TSTEP/3),NINT(RDAY)
            !END DO
        ENDDO

        DO JJ=1,NPOI
            IF (SUMZANG(JJ,1).GT.0._JPRD) THEN
                !ZWS(JJ)=ZMU0M(JJ)/SUMZANG(JJ,1)
                ZWS(JJ)=SUMZANG(JJ,STACT+1)/SUMZANG(JJ,1)
            ELSE
                ZWS(JJ)=0._JPRD
            ENDIF 
            IF ( SUMZANG(JJ,1) .LE. 0._JPRD .AND. SRFFI(JJ,MIN(IFFF+1, NSTPFC)) .GT. 0._JPRD ) THEN
                !print*, 'in dtfor,jj,SUMZANG(JJ,1),SRFFI(JJ,IFFF+1)',jj,SUMZANG(JJ,1),SRFFI(JJ,IFFF+1)
                ZWS(JJ) = (TSTEP/DTIMFC)
            ENDIF
        ENDDO
    ENDIF ! SOLAR ANGLE 


   IF (IFFP1 > NSTPFC) THEN
      IFFP1=IFF
      ZWf=1._JPRD
      ZWfP1=0._JPRD
   ENDIF
   
ENDIF ! NACCTYPE == 2

! print*, 'DBG'
! print*,'RTIMST',RTIMST
! print*,'ztimcur-RTIMST',ztimcur-RTIMST
! print*,'rtimtr-RTIMST',rtimtr-RTIMST
! print*,'rtimtr-RTIMST',rtimtr-RTIMST
! print*,'RTSTFC-ztimcur',RTSTFC-ztimcur

IF ( NACCTYPE .eq.2 .AND. LPREINT ) THEN 
   !# COMPUTE THE WEIGHTS FOR PRECIP TYPE PARAMETERS TO DISTRIBUTE INTO SUBSTEPS
   !# THREE CONSECUTIVE FORCING TIME STEPS ARE USED.
   !# ONE BEFORE AND ONE AFTER THE ACTUAL ONE THAT WE DEVIDE INTO SUBSTEPS
   !# THE DIVISION IS DONE IN A WAY TO PRESERVE THE ORIGINAL FORECAST VALUE (AS SUM OF SUBVALUES)
   ZWa(:,:)=0.
   IFF1=IFFF !actually = (rtimtr-RTSTFC)/DTIMFC+1
   IFF2=IFF1+1
   IFF3=IFF2+1
   IF (IFF1.LE.1) THEN
      IFF1=-99
   ENDIF
   IF (IFF2.EQ.NSTPFC) THEN
      IFF3=-99
   ENDIF
   IF (IFF1.LT.0 .AND. IFF3.LT.0) then
      WRITE(*,*) "ERROR - The precip values in the 1st and 3rd time steps are both empty"
      WRITE(*,*) "BUT only one can be empty at the beginning or at the end of the period"
      STOP 2
   ELSE IF (IFF1.GT.0 .AND. IFF3.LT.0) THEN
      IF (NFORC.EQ.2) THEN
         ZWa(1,1)=0.6
         ZWa(1,2)=0.4
         ZWa(2,1)=0.2
         ZWa(2,2)=0.8
      ELSE IF (NFORC.EQ.3) THEN
         ZWa(1,1)=0.6
         ZWa(1,2)=0.4
         ZWa(2,1)=0.2
         ZWa(2,2)=0.8
         ZWa(3,2)=1.0
      ELSE IF (NFORC.EQ.6) THEN
         ZWa(1,1)=0.8
         ZWa(1,2)=0.2
         ZWa(2,1)=0.6
         ZWa(2,2)=0.4
         ZWa(3,1)=0.4
         ZWa(3,2)=0.6
         ZWa(4,1)=0.2
         ZWa(4,2)=0.8
         ZWa(5,1)=0.1
         ZWa(5,2)=0.9
         ZWa(6,2)=1.0
      ENDIF
   ELSE IF (IFF1.LT.0 .AND. IFF3.GT.0) THEN
      IF (NFORC.EQ.2) THEN
         ZWa(1,2)=1.0
         ZWa(2,2)=0.4
         ZWa(2,3)=0.6
      ELSE IF (NFORC.EQ.3) THEN
         ZWa(1,2)=1.0
         ZWa(2,2)=0.8
         ZWa(2,3)=0.2
         ZWa(3,2)=0.4
         ZWa(3,3)=0.6
      ELSE IF (NFORC.EQ.6) THEN
         ZWa(1,2)=1.0
         ZWa(2,2)=0.9
         ZWa(2,3)=0.1
         ZWa(3,2)=0.8
         ZWa(3,3)=0.2
         ZWa(4,2)=0.6
         ZWa(4,3)=0.4
         ZWa(5,2)=0.4
         ZWa(5,3)=0.6
         ZWa(6,2)=0.2
         ZWa(6,3)=0.8
      ENDIF
   ELSE IF (IFF1.GT.0 .AND. IFF3.GT.0) THEN
      IF (NFORC.EQ.2) THEN
         ZWa(1,1)=0.6
         ZWa(1,2)=0.4
         ZWa(2,2)=0.4
         ZWa(2,3)=0.6
      ELSE IF (NFORC.EQ.3) THEN
         ZWa(1,1)=0.6
         ZWa(1,2)=0.4
         ZWa(2,2)=1.0
         ZWa(3,2)=0.4
         ZWa(3,3)=0.6
      ELSE IF (NFORC.EQ.6) THEN
         ZWa(1,1)=0.8
         ZWa(1,2)=0.2
         ZWa(2,1)=0.5
         ZWa(2,2)=0.5
         ZWa(3,1)=0.2
         ZWa(3,2)=0.8
         ZWa(4,2)=0.8
         ZWa(4,3)=0.2
         ZWa(5,2)=0.5
         ZWa(5,3)=0.5
         ZWa(6,2)=0.2
         ZWa(6,3)=0.8
      ENDIF
   ENDIF
ENDIF 



! Parameters with default linear interpolation (instantaneous parameters)
UNLEV0(:)=ZW*UFI(:,IF)+ZWP1*UFI(:,IFP1)
VNLEV0(:)=ZW*VFI(:,IF)+ZWP1*VFI(:,IFP1)
TNLEV0(:)=ZW*TFI(:,IF)+ZWP1*TFI(:,IFP1)
QNLEV0(:)=ZW*QFI(:,IF)+ZWP1*QFI(:,IFP1)
PNLP0(:)=ZW*PSFI(:,IF)+ZWP1*PSFI(:,IFP1)

! FLUXES INTERPLOATION 
FSTRD(:)=ZWf*TRFFI(:,IFF)+ZWfP1*TRFFI(:,IFFP1)
! SW radiation Interpolated based on Solar zenith angle
IF ( LSWINT ) THEN
    IFSOLAR=IFFF+1
    ! for last timestep
    IF (IFSOLAR == NSTPFC+1) THEN
        IFSOLAR=IFSOLAR-1
    ELSEIF (IFSOLAR > NSTPFC+1) THEN
        WRITE(NULOUT,*) " STOP IN ROUTINE DTFORC"
        WRITE(NULOUT,*) " IFSOLAR > NSTPFC+1"
        STOP
    ENDIF
    FSSRD(:)=(DTIMFC/TSTEP)*ZWS(:)*SRFFI(:,IFSOLAR)
ELSE
    FSSRD(:)=ZWf*SRFFI(:,IFF)+ZWfP1*SRFFI(:,IFFP1)
ENDIF




!PV
IF (LDBGS1) THEN
  JL=1
  WRITE(NULOUT,*) ' JL= ',JL
  WRITE(NULOUT,*) ' IF,IFp1 = ',IF,IFp1
  WRITE(NULOUT,*) ' IFF,IFFp1 = ',IFF,IFFp1
  WRITE(NULOUT,*) ' ZW,ZWP1= ',ZW,ZWP1
  WRITE(NULOUT,*) ' ZWf,ZWfP1= ',ZWf,ZWfP1
  WRITE(NULOUT,*) ' UFI(JL,IF),VFI(JL,IF),TFI(JL,IF),QFI(JL,IF),&
   &PSFI(JL,IF),SRFFI(JL,IFF),TRFFI(JL,IFF)= '
  WRITE(NULOUT,*) UFI(JL,IF),VFI(JL,IF),TFI(JL,IF),QFI(JL,IF),&
   &PSFI(JL,IF),&
   &SRFFI(JL,IFF),TRFFI(JL,IFF)
  WRITE(NULOUT,*) ' UFI(JL,IFP1),VFI(JL,IFP1),TFI(JL,IFP1)&
   &QFI(JL,IFP1),PSFI(JL,IFP1),SRFFI(JL,IFFP1),TRFFI(JL,IFFP1)= '
  WRITE(NULOUT,*) UFI(JL,IFP1),VFI(JL,IFP1),TFI(JL,IFP1),&
   &QFI(JL,IFP1),PSFI(JL,IFP1),&
   &SRFFI(JL,IFFP1),TRFFI(JL,IFFP1)
  WRITE(NULOUT,*) ' R30FI(JL,IFF)= ',R30FI(JL,IFF)
  WRITE(NULOUT,*) ' R30FI(JL,IFFP1)= ',R30FI(JL,IFFP1)
  WRITE(NULOUT,*) ' S30FI(JL,IFF)= ',S30FI(JL,IFF)
  WRITE(NULOUT,*) ' S30FI(JL,IFFP1)= ',S30FI(JL,IFFP1)
  WRITE(NULOUT,*) ' FSSRD(JL): ',FSSRD(JL)
ENDIF


! for precipitation(s) there are TWO options:
! 1. take rate in current forcing step (IFPREC)
! 2. forcing values are distributed in the period according to the relation to the neighburing periods

!# (1st option) SIMPLY CURRENT FORCING STEP VALUE IS TAKEN
!##########
IF (NACCTYPE.eq.0) THEN
   ! Fluxed assumed centred on the timestamp |---X---|
   IF(ZWf > 0.5)THEN
      IFPREC=IFF
   ELSE
      IFPREC=IFFP1
   ENDIF
ENDIF
IF (NACCTYPE.eq.1) THEN
   ! Fluxed beginning at the timestamp |X------| FORWARD ACCUMULATION
   IFPREC=IFF
ENDIF
IF (NACCTYPE.eq.2) THEN
   ! Fluxed ending at the timestamp |------X| BACKWARD ACCUMULATION
   IFPREC=IF+1 !IFFP1
ENDIF

! for last timestep
IF (IFPREC == NSTPFC+1) THEN
      IFPREC=IFPREC-1
ELSEIF (IFPREC > NSTPFC+1) THEN
    WRITE(NULOUT,*) " STOP IN ROUTINE DTFORC"
    WRITE(NULOUT,*) " IFPREC > NSTPFC+1"
    STOP
ENDIF

FLSRF(:)=R30FI(:,IFPREC)
FLSSF(:)=S30FI(:,IFPREC)
FCRF(:)=R30FI_C(:,IFPREC)
FCSF(:)=S30FI_C(:,IFPREC)


!# (3rd option) MORE REALISTIC LOOKING DISTRIBUTION IS APPLIED (WEIGHTING BY THE NEIGHBURING STEPS)
!# Applied ONLY for backward accumulation and ONLY IF we devide the forcing period into 2, 3 or 6 parts
!##########
IF (LPREINT .AND. NACCTYPE.eq.2 .and. (NFORC.eq.2 .or. NFORC.eq.3 .or. NFORC.eq.6)) THEN
   ! Fluxed ending at the timestamp |------X| BACKWARD ACCUMULATION
   ! This interpolation works only for division into 2, 3 or 6 parts
   ! (depending on the relation between the length of the forcing and the model time step)

   !# COMPUTE THE WEIGHTS (ZWTMP and ZWSUM)
   DO JJ=1,NPOI
      DO JK=1,NFORC
         TP2=DTIMFC*(R30FI(JJ,IFF2)+S30FI(JJ,IFF2))
         IF (TP2.LE.0.001) THEN
            ZWTMP(JJ,JK)=0.
         ELSE
            IF (IFF1.LT.0) THEN
               TP3=DTIMFC*(R30FI(JJ,IFF3)+S30FI(JJ,IFF3))
               ZWTMP(JJ,JK)=(ZWa(JK,2)*TP2+ZWa(JK,3)*TP3)/TP2
            ELSE IF (IFF3.LT.0) THEN
               TP1=DTIMFC*(R30FI(JJ,IFF1)+S30FI(JJ,IFF1))
               ZWTMP(JJ,JK)=(ZWa(JK,1)*TP1+ZWa(JK,2)*TP2)/TP2
            ELSE
               TP1=DTIMFC*(R30FI(JJ,IFF1)+S30FI(JJ,IFF1))
               TP3=DTIMFC*(R30FI(JJ,IFF3)+S30FI(JJ,IFF3))
               ZWTMP(JJ,JK)=(ZWa(JK,1)*TP1+ZWa(JK,2)*TP2+ZWa(JK,3)*TP3)/TP2
            END IF
         END IF 
      END DO
   END DO

   ZWSUM(:)=0.
   DO JK=1,NFORC
      ZWSUM(:)=ZWSUM(:)+ZWTMP(:,JK)
   ENDDO

   DO JJ=1,NPOI
      IF (ZWSUM(JJ).GT.0.001) THEN
         FLSRF(JJ)=NFORC*R30FI(JJ,IFF2)*(ZWTMP(JJ,STACT)/ZWSUM(JJ))
         FLSSF(JJ)=NFORC*S30FI(JJ,IFF2)*(ZWTMP(JJ,STACT)/ZWSUM(JJ))
      ELSE
         FLSRF(JJ)=0.
         FLSSF(JJ)=0.
      END IF
   ENDDO
ENDIF



IFP=(ztimcur+TSTEP-RTSTFC)/DTIMFC+1
IFPP1=IFP+1
ZWP=1._JPRD-(ztimcur+TSTEP-RTSTFC-REAL(IFP-1,KIND=JPRD)*DTIMFC)/DTIMFC
ZWPP1=1._JPRD-ZWP
IF (IFP == NSTPFC + 1 .OR. IFPP1 == NSTPFC + 1 ) THEN
  ! LAST interpolation for future timesetep
  IFP=NSTPFC
  ZWP=1._JPRD
  IFPP1=NSTPFC
  ZWPP1=0._JPRD
ELSEIF(IFP > NSTPFC + 1 ) THEN
  WRITE(*,*) 'out of bounds IFP',IFP,NSTPFC+1
  CALL ABORT()
ENDIF


IF (LDBGS1) THEN
  JL=1
  WRITE(NULOUT,*) ' JL= ',JL
  WRITE(NULOUT,*) ' ZWP,ZWPP1= ',ZWP,ZWPP1
  WRITE(NULOUT,*)' UFI(IFP),VFI(IFP),TFI(IFP),QFI(IFP),PSFI(IFP)='
  WRITE(NULOUT,*) UFI(JL,IFP),VFI(JL,IFP),TFI(JL,IFP)
  WRITE(NULOUT,*) ' UFI(IFPP1),VFI(IFPP1),TFI(IFPP1),QFI(IFPP1),'&
   &,'PSFI(IFPP1)= '
  WRITE(NULOUT,*) UFI(JL,IFPP1),VFI(JL,IFPP1),TFI(JL,IFPP1)&
   &,QFI(JL,IFPP1)&
   &,PSFI(JL,IFPP1)
ENDIF

UNLEV1(:)=ZWP*UFI(:,IFP)+ZWPP1*UFI(:,IFPP1)
VNLEV1(:)=ZWP*VFI(:,IFP)+ZWPP1*VFI(:,IFPP1)
TNLEV1(:)=ZWP*TFI(:,IFP)+ZWPP1*TFI(:,IFPP1)
QNLEV1(:)=ZWP*QFI(:,IFP)+ZWPP1*QFI(:,IFPP1)
PNLP1(:) =ZWP*PSFI(:,IFP)+ZWPP1*PSFI(:,IFPP1)


IF ( IDBGS1 > 1 ) THEN
  WRITE(NULOUT,"('dtfor:UVTQ(t)   =Z1*U(t1)+Z2*U(t2),Z1=',F5.2,' Z2=',F5.2,' t1=',I6,' t2=',I6)") ZW,ZWP1,IF,IFP1
  WRITE(NULOUT,"('dtfor:RAD(t/t+1)=Z1*U(t1)+Z2*U(t2),Z1=',F5.2,' Z2=',F5.2,' t1=',I6,' t2=',I6)") ZWf,ZWfP1,IFF,IFFP1
  IF ( LSWINT ) THEN
      WRITE(NULOUT,"('dtfor:SOL(t/t+1)=Z1*SWdown(t1),Z1=',F5.2,' t1=',I6)") (DTIMFC/TSTEP)*ZWS(1),IFSOLAR
  ENDIF 
  WRITE(NULOUT,"('dtfor:PRE(t/t+1)=Rainf(t1),t1=',I6,' NSTPFC=',I6)") IFPREC,NSTPFC
  WRITE(NULOUT,"('dtfor:UVTQ(t+1) =Z1*U(t1)+Z2*U(t2),Z1=',F5.2,' Z2=',F5.2,' t1=',I6,' t2=',I6)") ZWP,ZWPP1,IFP,IFPP1
ENDIF
UNLEV0(:)=ZW*UFI(:,IF)+ZWP1*UFI(:,IFP1)
VNLEV0(:)=ZW*VFI(:,IF)+ZWP1*VFI(:,IFP1)
TNLEV0(:)=ZW*TFI(:,IF)+ZWP1*TFI(:,IFP1)
QNLEV0(:)=ZW*QFI(:,IF)+ZWP1*QFI(:,IFP1)
PNLP0(:)=ZW*PSFI(:,IF)+ZWP1*PSFI(:,IFP1)

! Use variable atmospheric CO2 from the forcing fields

IF (LEAIRCO2COUP) THEN

  CNLEV0(:,1)=ZW*CO2FI(:,IF)+ZWP1*CO2FI(:,IFP1)
  CNLEV1(:,1)=ZWP*CO2FI(:,IFP)+ZWPP1*CO2FI(:,IFPP1)

ELSE

! Include atmospheric CO2 in the forcing fields (from global annual mean)
! Using values from timeseries used in IFS radiation code from CMIP6 file (CY48)

  call dattim(ZJUL,IYMD,IHM)
  YYYY=INT(IYMD/10000._JPRB)
  IYYYY=MIN(MAX(1,YYYY-RCO2REFYEAR+1),RNCO2YEARS)

  ZCO2=RANCO2(IYYYY)
  WRITE(*,*) ' ***** ATMOSPHERIC CO2:',ZCO2

  CNLEV0(:,:)=ZCO2*RMCO2/(RMD*1000000._JPRB)
  CNLEV1(:,:)=ZCO2*RMCO2/(RMD*1000000._JPRB)

ENDIF
IF (LHOOK) CALL DR_HOOK('DTFORC',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE DTFORC
