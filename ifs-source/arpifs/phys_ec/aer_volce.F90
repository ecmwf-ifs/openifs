! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_VOLCE&
 &(YDGEOMETRY, YDEAERATM,YDEAERVOL,YDML_GCONF,KIDIA, KFDIA, KLON  , KTDIA, KLEV , KSTART, KSTEP, KSTGLO, &
 &KTRAC, KAERO, KNBAER,&
 &PALTH, PAPHI, PCEN  , PGLAT, PGLON, PRHO  ,&
 &PCFLX, PTENC&
 &)

!*** * AER_VOLCE* - SOURCE TERMS OF AEROSOLS FOR VOLCANIC AEROSOLS 

!**   INTERFACE.
!     ----------
!          *AER_VOLCE* IS CALLED FROM *AER_SRC*

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        ORIGINAL : 2011-07-19

!     MODIFICATIONS.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,  DR_HOOK, JPHOOK
USE YOMCST    ,ONLY : RA, RPI, RDAY, RG
USE YOMLUN    ,ONLY : NULOUT
USE YOMRIP0   ,ONLY : NINDAT, NSSSSS
USE YOEAERATM ,ONLY : TEAERATM
!! USE YOEAERSRC ,ONLY : NTYPAER
USE YOEAERVOL ,ONLY : TEAERVOL
USE YOE_AERVOLE,ONLY: RVMASSVI, RVHGHTEM

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS IN
!              ------------

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TEAERATM)    ,INTENT(INOUT) :: YDEAERATM
TYPE(TEAERVOL)    ,INTENT(INOUT) :: YDEAERVOL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON  , KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV  , KTDIA, KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART, KSTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO(YDML_GCONF%YGFL%NAERO)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PALTH(KLON,0:KLEV), PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGLAT(KLON), PGLON(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHO(KLON,KLEV)


!*       0.2   ARGUMENTS OUT or INOUT
!              ----------------------

INTEGER(KIND=JPIM),INTENT(INOUT) :: KNBAER

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFLX(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCEN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: IBASPL, ITOPPL, IGVOLC, IGVOLE, IGLGLO
INTEGER(KIND=JPIM) :: JVOLC, JVOLE, JERUP, JAERVOLC
INTEGER(KIND=JPIM) :: JK, JL

INTEGER(KIND=JPIM) :: INDLAT(KLON)
INTEGER(KIND=JPIM) :: IVDATES, IVDAY, IVSECND
INTEGER(KIND=JPIM) :: IYY, IMM, IDD, IMDATE 
INTEGER(KIND=JPIM) :: IY0, IM0, ID0, INC, IMON(12)

REAL(KIND=JPRB)    :: ZEQUATOR, ZDLAT, ZDLON, ZGRDLAT, ZGRDLAT2
REAL(KIND=JPRB)    :: ZAREA, ZGELAV, ZGELAA, ZGELAB, ZGELOV, ZCIRCLE , ZGELOA , ZGELOB
REAL(KIND=JPRB)    :: ZDEEPLUME, ZVOLUME, ZLATSQ, ZLONSQ, ZPI2, ZPIH, Z1GP
REAL(KIND=JPRB)    :: ZADDMASS, ZZ0, ZZ1, ZDZ, ZVKGKG

!! REAL(KIND=JPRB)    :: ZGELAT(KLON)  , ZGRDLON(KLON), ZGRDLON2
REAL(KIND=JPRB)    :: ZGDLAT(KLON), ZGDLON(KLON)
REAL(KIND=JPRB)    :: ZGRDLON(KLON)
REAL(KIND=JPRB)    :: ZVOLCASH(KLON), ZVOLCSO2(KLON)
REAL(KIND=JPRB)    :: ZVOLEASH(KLON), ZVOLESO2(KLON)

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "updcal.intfb.h"
#include "fcttim.func.h"


IF (LHOOK) CALL DR_HOOK('AER_VOLCE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDRIP=>YDML_GCONF%YRRIP,YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NAERO=>YGFL%NAERO, &
 & NDGLG=>YDDIM%NDGLG, &
 & LAERVOL=>YDEAERATM%LAERVOL, &
 & NAERVOLC=>YDEAERVOL%NAERVOLC, NAERVOLE=>YDEAERVOL%NAERVOLE, &
 & NVOLDATS=>YDEAERVOL%NVOLDATS, NVOLERUZ=>YDEAERVOL%NVOLERUZ, &
 & RAERVOLC=>YDEAERVOL%RAERVOLC, RAERVOLE=>YDEAERVOL%RAERVOLE, &
 & RVOLERUZ=>YDEAERVOL%RVOLERUZ, &
 & NLOENG=>YDGEM%NLOENG, &
 & NGLOBALAT=>YDMP%NGLOBALAT, &
 & NSTASS=>YDRIP%NSTASS, RSTATI=>YDRIP%RSTATI)
!-----------------------------------------------------------------------

!*       1.0   TIME AND DATE OF THE MODEL
!              --------------------------
 
IY0=NCCAA(NINDAT)
IM0=NMM(NINDAT)
ID0=NDD(NINDAT)
INC=(NSSSSS + NINT(RSTATI))/NINT(RDAY)
CALL UPDCAL (ID0, IM0, IY0, INC,  IDD, IMM, IYY, IMON, -1)
IMDATE=IYY*10000+IMM*100+IDD

!-----------------------------------------------------------------------

!*       6.0   FLY ASH FROM VOLCANIC ERUPTIONS
!              -------------------------------
!IF (NTYPAER(6) /= 0) THEN
IF (LAERVOL) THEN
!
  ZDLAT  = 180._JPRB / NDGLG

  DO JL=KIDIA,KFDIA
!-- surface flux is set to zero as tendency is directly updated over the volcano.
    PCFLX(JL,KAERO(KNBAER+1))=0._JPRB
    IGLGLO=NGLOBALAT(KSTGLO+JL-1)
!    IGLGLO=MYLATS(JL)
    INDLAT(JL)=IGLGLO
    Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
    ZDLON=Z1GP*2.0_JPRB*RPI
    ZGDLAT(JL)=ZDLAT    ! latitude increment in degrees
    ZGDLON(JL)=360._JPRB*Z1GP                       ! longitude increment in degrees
  ENDDO

  ZPI2 = 2 * RPI
  ZPIH = RPI/2._JPRB
  ZDLAT  = 180._JPRB / NDGLG
  ZEQUATOR = ZPI2 * RA
  ZLATSQ = ZEQUATOR / (2*NDGLG)  ! length (m) of the latitude side of the grid
  ZGRDLAT= 180._JPRB / NDGLG
  ZGRDLAT2= ZGRDLAT*0.55_JPRB

! use information in YOEAERVOL to set the fluxes over the volcanoes
  
  JAERVOLC=0
  IF (JAERVOLC /= 0) THEN
!-- Continous volcanoes
    DO JVOLC=1,NAERVOLC

      ZGELAV=RAERVOLC(JVOLC,1)
      ZGELOV=RAERVOLC(JVOLC,2)
!-- find closest latitude points on the model grid
      ZCIRCLE=ZEQUATOR*COS(ZGELAV*RPI/180._JPRB)

      DO JL=KIDIA,KFDIA
        IGLGLO=NGLOBALAT(KSTGLO+JL-1)
!        IGLGLO=MYLATS(JL)
        INDLAT(JL)=IGLGLO
        Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
        ZDLON=Z1GP*2.0_JPRB*RPI                      ! longitude increment in radians
        ZVOLCASH(JL)=0._JPRB
        ZVOLCSO2(JL)=0._JPRB
        IBASPL=KLEV
        ITOPPL=1
        IGVOLC=0
        Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
        ZGRDLON(JL)=ZDLON
        ZGDLAT(JL)=ZDLAT
        ZGDLON(JL)=360._JPRB*Z1GP 

        ZGELAA = MIN( 90._JPRB, RAERVOLC(JVOLC,1)+ZGDLAT(JL)*0.55_JPRB)
        ZGELAB = MAX(-90._JPRB, RAERVOLC(JVOLC,1)-ZGDLAT(JL)*0.55_JPRB)
        ZGELOA = MAX(  0._JPRB, RAERVOLC(JVOLC,2)-ZGDLON(JL)*0.55_JPRB)
        ZGELOB = MIN(360._JPRB, RAERVOLC(JVOLC,2)+ZGDLON(JL)*0.55_JPRB)

        IF (PGLON(JL) >= ZGELOA .AND. PGLON(JL) < ZGELOB .AND.&
          & PGLAT(JL) <= ZGELAA .AND. PGLAT(JL) > ZGELAB ) THEN

          IGVOLC=JL
          ZLONSQ=ZEQUATOR * ZGRDLON(JL) / (2._JPRB*RPI) ! length (m) of the longitude side of the grid 
          ZAREA=ZLATSQ*ZLONSQ
          DO JK=1,KLEV
            IF (RAERVOLC(JVOLC,5) < PALTH(JL,JK) .AND. RAERVOLC(JVOLC,5) >= PALTH(JL,JK+1)) THEN
              IBASPL=JK
            ENDIF
            IF (RAERVOLC(JVOLC,6) < PALTH(JL,JK) .AND. RAERVOLC(JVOLC,6) >= PALTH(JL,JK+1)) THEN
              ITOPPL=JK
            ENDIF
          ENDDO
          IF (IBASPL == ITOPPL) THEN
            ITOPPL=ITOPPL-1
          ENDIF
          ZDEEPLUME=(PAPHI(JL,ITOPPL)-PAPHI(JL,IBASPL))/RG
          ZVOLUME=ZAREA*ZDEEPLUME
          ZVOLCASH(JL)=RAERVOLC(JVOLC,3)*RAERVOLC(JVOLC,7)/ZVOLUME
          ZVOLCSO2(JL)=RAERVOLC(JVOLC,4)*RAERVOLC(JVOLC,8)/ZVOLUME

          DO JK=ITOPPL,IBASPL
!- updating fly ash tendency (should be in variable 13)
            PTENC(JL,JK,KAERO(KNBAER+1))=PTENC(JL,JK,KAERO(KNBAER+1))+&
             & ZVOLCASH(JL)*(PAPHI(JL,JK-1)-PAPHI(JL,JK))/RG
!- updating volcanic SO2 tendency (should be in variable 15)
            PTENC(JL,JK,KAERO(KNBAER+2))=PTENC(JL,JK,KAERO(KNBAER+3))+&
             & ZVOLCSO2(JL)*(PAPHI(JL,JK-1)-PAPHI(JL,JK))/RG
          ENDDO

        ENDIF
      ENDDO

    ENDDO
  ENDIF

  IF (NAERVOLE /= 0) THEN    
!-- Explosive volcanoes
    DO JVOLE=1,NAERVOLE
      ZGELAV=RAERVOLE(JVOLE,1)
      ZGELOV=RAERVOLE(JVOLE,2)
!-- find closest latitude points on the model grid
      ZCIRCLE=ZEQUATOR*COS(ZGELAV*RPI/180._JPRB)
      IVDATES=NVOLDATS(JVOLE)
      IVDAY=IVDATES/100
      IVSECND=(IVDATES-IVDAY*100)*3600
      IF (KSTEP <= KSTART + 10) THEN
        WRITE(NULOUT,FMT='(" IVDATES IVDAY IVSECND:",2x,I10,2x,I8,2x,I6," NINDAT=",I10,&
     &                     " NSSSSS=",I6," IMDATE=",I10," NSTASS=",I10)')&
     &                     IVDATES,IVDAY,IVSECND,NINDAT,NSSSSS,IMDATE,NSTASS
      ENDIF

      IF (IMDATE > IVDAY .OR. (IMDATE == IVDAY .AND. NSTASS >= IVSECND)) THEN 

        DO JL=KIDIA,KFDIA
          IGLGLO=NGLOBALAT(KSTGLO+JL-1)
!          IGLGLO=MYLATS(JL)
          INDLAT(JL)=IGLGLO
          Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
          ZDLON=Z1GP*2.0_JPRB*RPI                      ! longitude increment in radians
          ZGRDLON(JL)=ZDLON
          ZVOLEASH(JL)=0._JPRB
          ZVOLESO2(JL)=0._JPRB
          IBASPL=KLEV
          ITOPPL=1
          IGVOLE=0

          ZGDLAT(JL)=ZDLAT
          ZGDLON(JL)=360._JPRB*Z1GP 

          ZGELAA = MIN( 90._JPRB, RAERVOLE(JVOLE,1)+ZGDLAT(JL)*0.55_JPRB)
          ZGELAB = MAX(-90._JPRB, RAERVOLE(JVOLE,1)-ZGDLAT(JL)*0.55_JPRB)
          ZGELOA = MAX(  0._JPRB, RAERVOLE(JVOLE,2)-ZGDLON(JL)*0.55_JPRB)
          ZGELOB = MIN(360._JPRB, RAERVOLE(JVOLE,2)+ZGDLON(JL)*0.55_JPRB)

          IF (PGLON(JL) >= ZGELOA .AND. PGLON(JL) < ZGELOB .AND.&
            & PGLAT(JL) <= ZGELAA .AND. PGLAT(JL) > ZGELAB ) THEN

            IGVOLE=JL
            ZLONSQ = ZEQUATOR * ZGRDLON(JL) / (2._JPRB*RPI)  ! length (m) of the longitude side of the grid 
            ZAREA = ZLATSQ*ZLONSQ

!-- distribution of the eruption mass over the vertical
!
!------ homogeneously between bottom and top of plumes (as given by RAERVOLE(5 .. and 6) 
            IF (NVOLERUZ(JVOLE) == 0) THEN
              DO JK=1,KLEV
                IF (RAERVOLE(JVOLE,5) < PALTH(JL,JK) .AND. RAERVOLE(JVOLE,5) >= PALTH(JL,JK+1)) THEN
                  IBASPL=JK
                ENDIF
                IF (RAERVOLE(JVOLE,6) < PALTH(JL,JK) .AND. RAERVOLE(JVOLE,6) >= PALTH(JL,JK+1)) THEN
                  ITOPPL=JK
                ENDIF
              ENDDO
              IF (IBASPL == ITOPPL) THEN
                ITOPPL=ITOPPL-1
              ENDIF
              ZDEEPLUME=(PAPHI(JL,ITOPPL)-PAPHI(JL,IBASPL))/RG
!-- NB: ZVOLUME refers to a "pile of grid-boxes" over the depth of the plume
              ZVOLUME=ZAREA*ZDEEPLUME
              ZVOLEASH(JL)=RAERVOLE(JVOLE,3)*RAERVOLE(JVOLE,7)/ZVOLUME
              ZVOLESO2(JL)=RAERVOLE(JVOLE,4)*RAERVOLE(JVOLE,8)/ZVOLUME

              DO JK=ITOPPL,IBASPL
!- updating fly ash tendency (should be in variable 13)
                PTENC(JL,JK,KAERO(KNBAER+1))=PTENC(JL,JK,KAERO(KNBAER+1))+&
                 & ZVOLEASH(JL)/PRHO(JL,JK)
!- updating volcanic SO2 tendency (should be in variable 15)
                PTENC(JL,JK,KAERO(KNBAER+2))=PTENC(JL,JK,KAERO(KNBAER+3))+&
                 & ZVOLESO2(JL)/PRHO(JL,JK)
              ENDDO

!------ as defined by profile information entered as an additional file or routine 
            ELSEIF (NVOLERUZ(JVOLE) == 1) THEN
              DO JK=1,KLEV
                ZZ0 = PAPHI(JL,JK-1)/RG
                ZZ1 = PAPHI(JL,JK)/RG
                ZDZ = ZZ0 - ZZ1
!-- NB: here ZVOLUME refers to a single grid box
                ZVOLUME=ZAREA*ZDZ

                DO JERUP=1,19
                  IF (ZZ1 < RVHGHTEM(JERUP-1) .AND. ZZ1 >= RVHGHTEM(JERUP)) THEN
!-- mass added from volcanic emission: kg s-1
                    ZADDMASS = (RVMASSVI(JERUP)-RVMASSVI(JERUP-1))&
                      & / (RVHGHTEM(JERUP-1)-RVHGHTEM(JERUP)) * ZDZ
!-- corresponding increment in mass mixing ratio  
!                 [ kg kg^-1 s^-1 = mass (kg s-1) / (density (kg m^-3) * volume (m^3)) ]
                    ZVKGKG = ZADDMASS / (PRHO(JL,JK) * ZVOLUME) * RVOLERUZ(JVOLE)

! NB: In this very specific configuration, volcanic matter is divided equally between 
!     varaibales 13 and 15
                    PTENC(JL,JK,KAERO(KNBAER+1))=PTENC(JL,JK,KAERO(KNBAER+1)) + ZVKGKG*0.5_JPRB
                    PTENC(JL,JK,KAERO(KNBAER+3))=PTENC(JL,JK,KAERO(KNBAER+3)) + ZVKGKG*0.5_JPRB
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_VOLCE',1,ZHOOK_HANDLE)
END SUBROUTINE AER_VOLCE
