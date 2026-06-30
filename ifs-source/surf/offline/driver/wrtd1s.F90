      SUBROUTINE WRTD1S
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *WRTP1S*  - Writing diagnostic variables of the one-column
!                 surface model

!     Purpose.
!     --------
!     Write out diagnostic variables

!**   Interface.
!     ----------
!        *CALL* *WRTD1S

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation 
!        of the one-column surface model

!     Author.
!     -------
!        Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-04-05

!     ------------------------------------------------------------------


      USE YOMGF1S  , ONLY : UNLEV0   ,VNLEV0   ,TNLEV0   ,QNLEV0   , &
     &            FSSRD    ,FSTRD    ,FLSRF    ,FCRF     ,FLSSF    ,FCSF
      USE YOMLUN1S , ONLY : NPOSDFO  ,NPOSDBD  ,&
     &            NPOSDTI1 ,NPOSDTI2 ,NPOSDTI3 ,NPOSDTI4 ,&
     &            NPOSDTI5 ,NPOSDTI6 ,NPOSDTI7 ,NPOSDTI8 ,&
     &            NPOSDT0  ,NPOSDST  ,&
     &            NPOSDT1  ,NPOSDT2  ,NPOSDT3  ,NPOSDT4  ,&
     &            NPOSDSW  ,NPOSDW0  ,&
     &            NPOSDW1  ,NPOSDW2  ,NPOSDW3  ,NPOSDW4,NPOSRC,&
     &            NPOSDTI9 ,& ! ENDUTRA 
     &            NPOSVEG ,NPOSCO2  

      USE YOMDYN1S , ONLY : NSTEP    ,TSTEP
      USE YOMDPHY, ONLY : NTILES, NLON
      USE YOMGDI1S , ONLY : &
     &            D1STISRD2,D1STISRU2,D1STITRD2,D1STITRU2,&
     &            D1STIH2  ,D1STILE2 ,D1STIGFL2,D1STII2  ,D1STIEVAP2,&
     &            D1SRFLD2 ,D1SRFLU2 ,D1TRFLD2 ,D1TRFLU2 ,&
     &            D1AHFS2  ,D1AHFL2  ,&
     &            D1STNSRD2,D1STNSRU2,D1STNTRD2,D1STNTRU2,&
     &            D1STNGFL2,D1STNM2  ,&
     &            D1ST1SRD2,D1ST1SRU2,D1ST1TRD2,D1ST1TRU2,&
     &            D1ST1H2  ,D1ST1LE2 ,D1STAGFL2,D1STASF2 , &
     &            D1SSFL2  ,D1SSFC2  ,D1SWNJQ2 ,D1SWNM2  ,&
     &            D1SWLIT2 ,D1SWLJQ2 ,&
     &            D1SW1JBG2,D1SW1TF2 ,D1SWAGFL2,D1SW1RI2 ,&
     &            D1SWARS2 ,D1SWAEXT2,D1SW1M2  ,D1SWAC2  ,&
     &            D1STIFR  ,D1STITK  ,D1STIALB ,D1SNDEPTH,D1SWAFR,&
     &            D1SVTRC2,  D1SVTRA2,D1T2M2   ,D1D2M2 &
! CO2	   
           &,D1SAN2,D1SAG2,D1SRD2,D1SRSOIL_STR2,D1SRECO,D1SCO2FLUX2 &
	   &,D1SVTAN2,D1SVTAG2,D1SVTRD2,D1SVTRSOIL_STR2,D1SVTRECO2,D1SVTCO2FLUX2 &
! LAI, biomass
           &,D1SLAI2,D1SBIOM2,D1SBLOSS2,D1SBGAIN2,D1SBIOMSTR2,D1SBIOMSTRB2 &
	   &,D1SVTLAI2,D1SVTBIOM2,D1SVTBLOSS2,D1SVTBGAIN2 &
	   &,D1SVTBIOMSTR2,D1SVTBIOMSTRB2,D1SVTRECO2 &
! vegetation variables
	   &,D1SVTGC2 ,D1SVTGA2 ,D1SVTF22 &
	   &,D1SVTLE2 ,D1SVTFR2 ,D1SVTDS2 ,D1SVTDMAX2 &
! fractions
           &,D1SWFR2 &
! tiled ouput
           &,D1STIALB2 &
! extra
	   &,D1SWDSSL2


      USE YOMLOG1S , ONLY : LACCUMW
      USE YOMCST   , ONLY : RDAY
      USE PARKIND1  , ONLY : JPIM,JPRB, JPRD
      IMPLICIT LOGICAL (L)

      INTEGER (KIND=JPIM) :: IPOS(NTILES),IYYMD,IHM
      REAL(KIND=JPRB) ::FOEEWMO
      REAL(KIND=JPRD) ::ZJUL
!     ------------------------------------------------------------------

!*       1.   WRITE OUT FORCING.
!             ------------------


!   time stuff
      CALL DATTIM(ZJUL,IYYMD,IHM)
      write(*,*) 'ZJUL=',zjul,IYYMD,IHM
      JL=1
      IF (NSTEP == 0) THEN
        WRITE(NPOSDFO,'(a)') '#    juld     date time      step'&
     &  //'            u            v            t            q'&
     &  //'        sradd        tradd        cnvrf        cnvsf'&
     &  //'         lsrf         lssf'&
     &  //'          T2m          D2m         RH2m'

        WRITE(NPOSDFO,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         ms-1         ms-1            K       kgkg-1'&
     &  //'         Wm-2         Wm-2     kgm-2s-1     kgm-2s-1'&
     &  //'     kgm-2s-1     kgm-2s-1'&
     &  //'          K            K            % '
      ENDIF
      WRITE(NPOSDFO,'(f10.3,1X,I8,1X,I4,1X,I8,13(1X,E13.6))')&
     &      zjul,IYYMD,IHM,NSTEP &
     &      ,UNLEV0(JL),VNLEV0(JL),TNLEV0(JL),QNLEV0(JL)&
     &     ,FSSRD(JL),FSTRD(JL) &
     &     ,FCRF(JL),FCSF(JL),FLSRF(JL),FLSSF(JL) &
     &     ,D1T2M2(JL,1),D1D2M2(JL,1)&
     &     ,100.0_JPRB*FOEEWMO(D1D2M2(JL,1))/FOEEWMO(D1T2M2(JL,1))

!     ------------------------------------------------------------------

!*       2.   WRITE OUT OVERALL BUDGETS.
!             --------------------------

      IF (LACCUMW) THEN
!     ZWA - FACTOR TO TRANSFORM W/M**2*S INTO W/M**2
!     ZMM - FACTOR TO TRANSFORM KG/M**2 INTO MM
        ZWA=1./((NSTEP+1)*TSTEP)
        ZMM=1.
        IA=2
      ELSE
!     ZWA - FACTOR TO TRANSFORM W/M**2 INTO W/M**2
!     ZMM - FACTOR TO TRANSFORM KG/M**2/S INTO MM/DAY
        ZWA=1.
        ZMM=RDAY
        IA=1
      ENDIF

      IF (NSTEP == 0) THEN
        WRITE(NPOSDBD,'(a)') '#    juld     date time      step'&
     &  //'         srad         trad          thf          nhf'&
     &  //'          e'
        if (laccumw) then
          WRITE(NPOSDBD,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &  //'        kgm-2'
        else
          WRITE(NPOSDBD,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &  //'     kgm-2d-1'
        endif
      ENDIF
      WRITE(NPOSDBD,'(f10.3,1X,I8,1X,I4,1X,I8,8(1X,E13.6))')&
     &      zjul,IYYMD,IHM,NSTEP,(D1SRFLD2(JL,IA)+D1SRFLU2(JL,IA))*ZWA&
     &     ,(D1TRFLD2(JL,IA)+D1TRFLU2(JL,IA))*ZWA&
     &     ,(D1AHFS2(JL,IA)+D1AHFL2(JL,IA))*ZWA&
     &     ,(D1SRFLD2(JL,IA)+D1SRFLU2(JL,IA)+D1TRFLD2(JL,IA)+&
     &      D1TRFLU2(JL,IA)+D1AHFS2(JL,IA)+D1AHFL2(JL,IA))*ZWA&
     &     ,(D1SWNJQ2(JL,IA)+D1SWLJQ2(JL,IA)+D1SW1JBG2(JL,IA)+&
     &      D1SWAEXT2(JL,1,IA)+D1SWAEXT2(JL,2,IA)+D1SWAEXT2(JL,3,IA)+&
     &      D1SWAEXT2(JL,4,IA))*ZMM

!     ------------------------------------------------------------------

!*       3.   WRITE OUT TILES BUDGET.
!             -----------------------

      IPOS(:)=(/NPOSDTI1, NPOSDTI2, NPOSDTI3, NPOSDTI4,&
     &          NPOSDTI5, NPOSDTI6, NPOSDTI7, NPOSDTI8,NPOSDTI9/) ! ENDUTRA 
      LDSEA=(D1STIFR(JL,1)+D1STIFR(JL,2)) >= 0.5
      LDLAND=.NOT. LDSEA
      TILES: DO JTI=1,NTILES
        IF ((JTI == 1 .OR. JTI == 2) .AND. LDLAND) CYCLE TILES
        IF (JTI >= 3 .AND. LDSEA) CYCLE TILES

        IF (NSTEP == 0) THEN
          WRITE(IPOS(JTI),'(a)') &
     &                         '#    juld     date time      step'&
     &    //'           fr            t          alb            e'&
     &    //'        sradd        sradu        tradd        tradu'&
     &    //'          shf          lhf         ghf0          imp'
          if (laccumw) then 
            WRITE(IPOS(JTI),'(a)') &
     &                           '#   jjj.j yyyymmdd hhmm          '&
     &      //'            -            K            -        kgm-2'&
     &      //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &      //'         Wm-2         Wm-2         Wm-2         Wm-2'
          else
            WRITE(IPOS(JTI),'(a)')&
     &                           '#   jjj.j yyyymmdd hhmm          '&
     &      //'            -            K            -     kgm-2d-1'&
     &      //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &      //'         Wm-2         Wm-2         Wm-2         Wm-2'
          endif
        ENDIF

        WRITE(IPOS(JTI),'(f10.3,1X,I8,1X,I4,1X,I8,12(1X,E13.6E3))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1STIFR(JL,JTI),D1STITK(JL,JTI)&
     &     ,D1STIALB(JL,JTI),D1STIEVAP2(JL,JTI,IA)*ZMM&
     &     ,D1STISRD2(JL,JTI,IA)*ZWA,D1STISRU2(JL,JTI,IA)*ZWA&
     &     ,D1STITRD2(JL,JTI,IA)*ZWA,D1STITRU2(JL,JTI,IA)*ZWA&
     &     ,D1STIH2(JL,JTI,IA)*ZWA,D1STILE2(JL,JTI,IA)*ZWA&
     &     ,D1STIGFL2(JL,JTI,IA)*ZWA,D1STII2(JL,JTI,IA)*ZWA
      ENDDO TILES
      
!     ------------------------------------------------------------------

!*       4.   WRITE OUT SKIN TEMPERATURE BUDGET.
!             ----------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDT0,'(a)') '#    juld     date time      step'&
     &  //'        sradd        sradu        tradd        tradu'&
     &  //'          shf          lhf         ghf0          imp'
        WRITE(NPOSDT0,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'
      ENDIF

!     WRITE(NPOSDT0,'(f10.3,1X,I8,1X,I4,1X,I8,8(1X,E13.6))')&
!    &      zjul,IYYMD,IHM,NSTEP&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STISRD2(JL,:,IA))*ZWA&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STISRU2(JL,:,IA))*ZWA&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STITRD2(JL,:,IA))*ZWA&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STITRU2(JL,:,IA))*ZWA&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STIH2(JL,:,IA))*ZWA&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STILE2(JL,:,IA))*ZWA&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STIGFL2(JL,:,IA))*ZWA&
!    &     ,DOT_PRODUCT(D1STIFR(JL,:),D1STII2(JL,:,IA))*ZWA

! All fluxes are already tile-area weighted therefore the sum gives the 
! grid-point flux

      WRITE(NPOSDT0,'(f10.3,1X,I8,1X,I4,1X,I8,8(1X,E13.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,SUM(D1STISRD2(JL,:,IA))*ZWA&
     &     ,SUM(D1STISRU2(JL,:,IA))*ZWA&
     &     ,SUM(D1STITRD2(JL,:,IA))*ZWA&
     &     ,SUM(D1STITRU2(JL,:,IA))*ZWA&
     &     ,SUM(D1STIH2(JL,:,IA))*ZWA&
     &     ,SUM(D1STILE2(JL,:,IA))*ZWA&
     &     ,SUM(D1STIGFL2(JL,:,IA))*ZWA&
     &     ,SUM(D1STII2(JL,:,IA))*ZWA

!     ------------------------------------------------------------------

!*       5.   WRITE OUT SNOW TEMPERATURE BUDGET.
!             ----------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDST,'(a)') '#    juld     date time      step'&
     &  //'           fr          alb            d'&
     &  //'        sradd        sradu        tradd        tradu'&
     &  //'          shf          lhf         ghfb         melt'
        WRITE(NPOSDST,'(a)') '#   jjj.j yyyymmdd hhmm&          '&
     &  //'            -            -            m'&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'
      ENDIF
      WRITE(NPOSDST,'(f10.3,1X,I8,1X,I4,1X,I8,11(1X,E13.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1STIFR(JL,5)+D1STIFR(JL,7)&
     &     ,D1STIFR(JL,5)*D1STIALB(JL,5)+D1STIFR(JL,7)*D1STIALB(JL,7)&
     &     ,D1SNDEPTH(JL)&
     &     ,D1STNSRD2(JL,IA)*ZWA,D1STNSRU2(JL,IA)*ZWA&
     &     ,D1STNTRD2(JL,IA)*ZWA,D1STNTRU2(JL,IA)*ZWA&
     &     ,(D1STIFR(JL,5)*D1STIH2(JL,5,IA)+&
     &      D1STIFR(JL,7)*D1STIH2(JL,7,IA))*ZWA&
     &     ,(D1STIFR(JL,5)*D1STILE2(JL,5,IA)+&
     &      D1STIFR(JL,7)*D1STILE2(JL,7,IA))*ZWA&
     &     ,D1STNGFL2(JL,IA)*ZWA,D1STNM2(JL,IA)*ZWA

!     ------------------------------------------------------------------

!*       6.   WRITE OUT SOIL TEMPERATURE LAYER 1 BUDGET.
!             ------------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDT1,'(a)') '#    juld     date time      step'&
     &  //'        sradd        sradu        tradd        tradu'&
     &  //'          shf          lhf        gsnow         ghfb'&
     &  //'       wfreez'
        WRITE(NPOSDT1,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &  //'         Wm-2         Wm-2         Wm-2         Wm-2'&
     &  //'         Wm-2'
      ENDIF
      WRITE(NPOSDT1,'(F10.3,1X,I8,1X,I4,1X,I8,9(1X,E13.6))')&
     &    zjul,IYYMD,IHM,NSTEP&
     &   ,D1ST1SRD2(JL,IA)*ZWA,D1ST1SRU2(JL,IA)*ZWA&
     &   ,D1ST1TRD2(JL,IA)*ZWA,D1ST1TRU2(JL,IA)*ZWA&
     &   ,D1ST1H2(JL,IA)*ZWA,D1ST1LE2(JL,IA)*ZWA&
     &   ,D1STNGFL2(JL,IA),D1STAGFL2(JL,1,IA)*ZWA&
     &   ,D1STASF2(JL,1,IA)*ZWA

!     ------------------------------------------------------------------

!*       7.   WRITE OUT SOIL TEMPERATURE LAYER 2 BUDGET.
!             ------------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDT2,'(a)') '#    juld     date time      step'&
     &  //'         ghft         ghfb        freez'          
        WRITE(NPOSDT2,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         Wm-2         Wm-2         Wm-2'
      ENDIF
      WRITE(NPOSDT2,'(F10.3,1X,I8,1X,I4,1X,I8,3(1X,E12.6E3))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1STAGFL2(JL,1,IA)*ZWA,D1STAGFL2(JL,2,IA)*ZWA&
     &     ,D1STASF2(JL,2,IA)*ZWA

!     ------------------------------------------------------------------

!*       8.   WRITE OUT SOIL TEMPERATURE LAYER 3 BUDGET.
!             ------------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDT3,'(a)') '#    juld     date time      step'&
     &  //'         ghft         ghfb       wfreez'
        WRITE(NPOSDT3,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         Wm-2         Wm-2         Wm-2'
      ENDIF
      WRITE(NPOSDT3,'(F10.3,1X,I8,1X,I4,1X,I8,3(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1STAGFL2(JL,2,IA)*ZWA,D1STAGFL2(JL,3,IA)*ZWA&
     &     ,D1STASF2(JL,3,IA)*ZWA

!     ------------------------------------------------------------------

!*       9.   WRITE OUT SOIL TEMPERATURE LAYER 4 BUDGET.
!             ------------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDT4,'(a)') '#    juld     date time      step'&
     &  //'         ghft       wfreez'
        WRITE(NPOSDT4,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'         Wm-2         Wm-2'
      ENDIF
      WRITE(NPOSDT4,'(F10.3,1X,I8,1X,I4,1X,I8,3(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1STAGFL2(JL,3,IA)*ZWA,D1STASF2(JL,4,IA)*ZWA

!     ------------------------------------------------------------------

!*      10.   WRITE OUT SNOW WATER BUDGET.
!             ----------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDSW,'(a)') '#    juld     date time      step'&
     &  //'            e           sf         melt'
        if (laccumw) then
          WRITE(NPOSDSW,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        kgm-2        kgm-2        kgm-2'
        else
          WRITE(NPOSDSW,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'     kgm-2d-1     kgm-2d-1     kgm-2d-1'
        endif
      ENDIF
      WRITE(NPOSDSW,'(F10.3,1X,I8,1X,I4,1X,I8,3(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1SWNJQ2(JL,IA)*ZMM,(D1SSFL2(JL,IA)+D1SSFC2(JL,IA))*ZMM&
     &     ,D1SWNM2(JL,IA)*ZMM

!     ------------------------------------------------------------------

!*      11.   WRITE OUT INTERCEPTION RESERVOIR WATER BUDGET.
!             ----------------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDW0,'(a)') '#    juld     date time      step'&
     &  //'          int            e'
        if (laccumw) then
          WRITE(NPOSDW0,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        kgm-2        kgm-2'
        else
          WRITE(NPOSDW0,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'     kgm-2d-1     kgm-2d-1'
        endif
      ENDIF
      WRITE(NPOSDW0,'(F10.3,1X,I8,1X,I4,1X,I8,2(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1SWLIT2(JL,IA)*ZMM,D1SWLJQ2(JL,IA)*ZMM

!     ------------------------------------------------------------------

!*       12.   WRITE OUT SOIL LAYER 1 WATER BUDGET.
!              ------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDW1,'(a)') '#    juld     date time      step'&
     &  //'         w_fr          tfa         melt          ebg'&
     &  //'          ere         gwfb         roil         rose'&
     &  //'         clip'
        if (laccumw) then
          WRITE(NPOSDW1,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3        kgm-2        kgm-2        kgm-2'&
     &  //'        kgm-2        kgm-2        kgm-2        kgm-2'&
     &  //'        kgm-2'
        else
          WRITE(NPOSDW1,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3     kgm-2d-1     kgm-2d-1     kgm-2d-1'&
     &  //'     kgm-2d-1     kgm-2d-1     kgm-2d-1     kgm-2d-1'&
     &  //'     kgm-2d-1'
        endif
      ENDIF
      WRITE(NPOSDW1,'(F10.3,1X,I8,1X,I4,1X,I8,9(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1SWAFR(JL,1),D1SW1TF2(JL,IA)*ZMM&
     &     ,D1SW1M2(JL,IA)*ZMM,D1SW1JBG2(JL,IA)*ZMM&
     &     ,D1SWAEXT2(JL,1,IA)*ZMM,D1SWAGFL2(JL,1,IA)*ZMM&
     &     ,D1SW1RI2(JL,IA)*ZMM,D1SWARS2(JL,1,IA)*ZMM&
     &     ,D1SWAC2(JL,1,IA)*ZMM

!     ------------------------------------------------------------------

!*       13.   WRITE OUT SOIL LAYER 2 WATER BUDGET.
!              ------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDW2,'(a)') '#    juld     date time      step'&
     &  //'         w_fr         gwft         gwfb         rose'&
     &  //'          ere         clip'
        if (laccumw) then
          WRITE(NPOSDW2,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3        kgm-2        kgm-2        kgm-2'&
     &  //'        kgm-2        kgm-2'
        else
          WRITE(NPOSDW2,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3     kgm-2d-1     kgm-2d-1     kgm-2d-1'&
     &  //'     kgm-2d-1     kgm-2d-1'
        endif
      ENDIF
      WRITE(NPOSDW2,'(F10.3,1X,I8,1X,I4,1X,I8,6(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1SWAFR(JL,2),D1SWAGFL2(JL,1,IA)*ZMM&
     &     ,D1SWAGFL2(JL,2,IA)*ZMM,D1SWARS2(JL,2,IA)*ZMM&
     &     ,D1SWAEXT2(JL,2,IA)*ZMM,D1SWAC2(JL,2,IA)*ZMM

!     ------------------------------------------------------------------

!*       14.   WRITE OUT SOIL LAYER 3 WATER BUDGET.
!             -------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDW3,'(a)') '#    juld     date time      step'&
     &  //'         w_fr         gwft         gwfb         rose'&
     &  //'          ere         clip'
        if (laccumw) then
          WRITE(NPOSDW3,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3        kgm-2        kgm-2        kgm-2'&
     &  //'        kgm-2        kgm-2'
        else
          WRITE(NPOSDW3,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3     kgm-2d-1     kgm-2d-1     kgm-2d-1'&
     &  //'     kgm-2d-1     kgm-2d-1'
        endif
      ENDIF
      WRITE(NPOSDW3,'(F10.3,1X,I8,1X,I4,1X,I8,6(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1SWAFR(JL,3),D1SWAGFL2(JL,2,IA)*ZMM&
     &     ,D1SWAGFL2(JL,3,IA)*ZMM,D1SWARS2(JL,3,IA)*ZMM&
     &     ,D1SWAEXT2(JL,3,IA)*ZMM,D1SWAC2(JL,3,IA)*ZMM

!     ------------------------------------------------------------------

!*       15.   WRITE OUT SOIL LAYER 4 WATER BUDGET.
!             -------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSDW4,'(a)') '#    juld     date time      step'&
     &  //'         w_fr         gwft         gwfb         rose'&
     &  //'          ere         clip'
        if (laccumw) then
          WRITE(NPOSDW4,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3        kgm-2        kgm-2        kgm-2'&
     &  //'        kgm-2        kgm-2'
        else
          WRITE(NPOSDW4,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        m3m-3     kgm-2d-1     kgm-2d-1     kgm-2d-1'&
     &  //'     kgm-2d-1     kgm-2d-1'

        endif
      ENDIF
      WRITE(NPOSDW4,'(F10.3,1X,I8,1X,I4,1X,I8,6(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &     ,D1SWAFR(JL,4),D1SWAGFL2(JL,3,IA)*ZMM&
     &     ,D1SWAGFL2(JL,4,IA)*ZMM,D1SWARS2(JL,4,IA)*ZMM&
     &     ,D1SWAEXT2(JL,4,IA)*ZMM,D1SWAC2(JL,4,IA)*ZMM

      
      WRITE(NPOSRC,'(F10.3,1X,I8,1X,I4,1X,I8,4(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP,&
     &      D1SVTRC2(JL,1,IA)*ZWA, D1SVTRC2(JL,2,IA)*ZWA,&
     &      D1SVTRA2(JL,1,IA)*ZWA, D1SVTRA2(JL,2,IA)*ZWA



!*       16.   WRITE OUT CO2,LAI, BIOMASS, 
!             ------------------------------------------------
!            (unit of CO2 fluxes kg_CO2/m2/s)
      IF (NSTEP == 0) THEN
        WRITE(NPOSCO2,'(a)') '#    juld     date time      step'&
     &  //'         An         Ag         Rd         Rsoil_str'&
     &  //'         CO2flux    Anvt       Anvt       Agvt'&
     &  //'         Agvt       Rdvt       Rdvt       Rsoil_strvt'&
     &  //'         Rsoil_strvt        CO2fluxvt       CO2fluxvt'&
     &  //'         lai       biomass       Bloss       Bgain'&
     &  //'         Biomstr       Biomstr2       laivt       laivt'&
     &  //'         biomassvt       biomassvt       Blossvt       Blossvt'&
     &  //'         Bgainvt       Bgainvt       Biomstrvt       Biomstrvt'&
     &  //'         Biomstr2vt       Biomstr2vt        Recovt       Recovt'

          WRITE(NPOSCO2,'(a)') '#   jjj.j yyyymmdd hhmm          '&
     &  //'        kgCO2m-2s-1     kgCO2m-2s-1     kgCO2m-2s-1     kgCO2m-2s-1'&
     &  //'     kgCO2m-2s-1     kgCO2m-2s-1     kgCO2m-2s-1     kgCO2m-2s-1'&
     &  //'     kgCO2m-2s-1     kgCO2m-2s-1     kgCO2m-2s-1     kgCO2m-2s-1'&
     &  //'     kgCO2m-2s-1     kgCO2m-2s-1     kgCO2m-2s-1'&
     &  //'     m2m-2     kgm-2     kgm-2     kgm-2'&
     &  //'     kgm-2     kgm-2     m2m-2     m2m-2'&
     &  //'     kgm-2     kgm-2     kgm-2     kgm-2'&
     &  //'     kgm-2     kgm-2     kgm-2     kgm-2'&
     &  //'     kgm-2     kgm-2     kgCO2m-2s-1     kgCO2m-2s-1'
      ENDIF
!
! WRITE(*,*)   'ZWA=',ZWA, '  ZMM=',ZMM,' IA=',IA

 WRITE(NPOSCO2,'(F10.3,1X,I8,1X,I4,1X,I8,35(1X,E14.6E3))')&
     &       zjul,IYYMD,IHM,NSTEP&
     &      ,D1SAN2(JL,IA)*ZWA,D1SAG2(JL,IA)*ZWA&
     &      ,D1SRD2(JL,IA)*ZWA,D1SRSOIL_STR2(JL,IA)*ZWA,D1SCO2FLUX2(JL,IA)*ZWA &
     &      ,D1SVTAN2(JL,1,IA)*ZWA,     D1SVTAN2(JL,2,IA)*ZWA      &
     &      ,D1SVTAG2(JL,1,IA)*ZWA,     D1SVTAG2(JL,2,IA)*ZWA      &
     &      ,D1SVTRD2(JL,1,IA)*ZWA,     D1SVTRD2(JL,2,IA)*ZWA      &
     &      ,D1SVTRSOIL_STR2(JL,1,IA)*ZWA,  D1SVTRSOIL_STR2(JL,2,IA)*ZWA   &     
     &      ,D1SVTCO2FLUX2(JL,1,IA)*ZWA,D1SVTCO2FLUX2(JL,2,IA)*ZWA &
     &      ,D1SLAI2(JL,IA)*ZWA,        D1SBIOM2(JL,IA)*ZWA        &
     &      ,D1SBLOSS2(JL,IA)*ZWA,      D1SBGAIN2(JL,IA)*ZWA       &
     &      ,D1SBIOMSTR2(JL,IA)*ZWA,    D1SBIOMSTRB2(JL,IA)*ZWA    &
     &      ,D1SVTLAI2(JL,1,IA)*ZWA,    D1SVTLAI2(JL,2,IA)*ZWA     &
     &      ,D1SVTBIOM2(JL,1,IA)*ZWA,   D1SVTBIOM2(JL,2,IA)*ZWA    &
     &      ,D1SVTBLOSS2(JL,1,IA)*ZWA,  D1SVTBLOSS2(JL,2,IA)*ZWA   &
     &      ,D1SVTBGAIN2(JL,1,IA)*ZWA,  D1SVTBGAIN2(JL,2,IA)*ZWA   &
     &      ,D1SVTBIOMSTR2(JL,1,IA)*ZWA,D1SVTBIOMSTR2(JL,2,IA)*ZWA &
     &      ,D1SVTBIOMSTRB2(JL,1,IA)*ZWA,D1SVTBIOMSTRB2(JL,2,IA)*ZWA &
     &      ,D1SVTRECO2(JL,1,IA)*ZWA,   D1SVTRECO2(JL,2,IA)*ZWA 

!*       17.   WRITE OUT VEGETATION VARIABLE 
!             ------------------------------------------------

      IF (NSTEP == 0) THEN
        WRITE(NPOSVEG,'(a)') '#    juld     date time      step'&
     &  //'         rc_h         rc_l         ra_h         ra_l'&
     &  //'         f2vt_h         f2vt_l         LEvt_h       LEvt_l'&
     &  //'         vtfr_h         vtfr_l         dsvt_h       dsvt_l'&
     &  //'         dmaxvt_h       dmaxvt_l       intfr'
      ENDIF

 WRITE(NPOSVEG,'(F10.3,1X,I8,1X,I4,1X,I8,15(1X,E12.6))')&
     &      zjul,IYYMD,IHM,NSTEP&
     &      ,D1SVTGC2(JL,1,IA)*ZWA,D1SVTGC2(JL,2,IA)*ZWA & 
     &      ,D1SVTGA2(JL,1,IA)*ZWA,D1SVTGA2(JL,2,IA)*ZWA &
     &      ,D1SVTF22(JL,1,IA)*ZWA,D1SVTF22(JL,2,IA)*ZWA &
     &      ,D1SVTLE2(JL,1,IA)*ZWA,D1SVTLE2(JL,2,IA)*ZWA &
     &      ,D1SVTFR2(JL,1,IA)*ZWA,D1SVTFR2(JL,2,IA)*ZWA &
     &      ,D1SVTDS2(JL,1,IA)*ZWA,D1SVTDS2(JL,2,IA)*ZWA &
     &      ,D1SVTDMAX2(JL,1,IA)*ZWA,D1SVTDMAX2(JL,2,IA)*ZWA &
! fractions
           &,D1SWFR2(JL,IA)*ZWA 

!write(*,*),'end wrtd1s ZJUL=',zjul,IYYMD,IHM


      RETURN
      END SUBROUTINE WRTD1S

FUNCTION FOEEWMO(T)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
!---------------------------------------------------------------------------

REAL(KIND=JPRB), PARAMETER :: R=1.380658*6.0221367
REAL(KIND=JPRB), PARAMETER :: RD=1000.*R/28.9644
REAL(KIND=JPRB), PARAMETER :: RV=1000.*R/18.0153
REAL(KIND=JPRB), PARAMETER :: RTT=273.16
REAL(KIND=JPRB), PARAMETER :: R2ES=611.21*RD/RV
REAL(KIND=JPRB), PARAMETER :: R3LES=17.502
REAL(KIND=JPRB), PARAMETER :: R4LES=32.19
REAL(KIND=JPRB), INTENT(IN) :: T
REAL(KIND=JPRB) :: FOEEWMO

FOEEWMO=R2ES*EXP(R3LES*(T-RTT)/(T-R4LES))
END FUNCTION FOEEWMO

