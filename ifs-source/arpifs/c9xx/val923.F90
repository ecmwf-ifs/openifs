! (C) Copyright 1989- Meteo-France.

SUBROUTINE VAL923(LDNEW)

!**** *GEO923*

!     PURPOSE.
!     --------
!      Compute the constants (YOMCLI) which are used by configuration 923.

!     INTERFACE.
!     ----------
!      CALL VAL923(...)
!        LDNEW = .FALSE. if old fields required
!        Results in YOMCLI.

!     AUTHORS.
!     --------
!      D. Giard 97-05-06

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCLI   , ONLY : NTPMER   ,NTPGLA   ,NTPDES   ,NTPLAC   ,&
 & SMASK    ,SMANQ    ,STHER    ,SALBN    ,SALBX    ,&
 & SALBM    ,SALBG    ,SALBB    ,SALBD    ,SEMIN    ,&
 & SEMIX    ,SEMIM    ,SEMIG    ,SEMIB    ,SEMID    ,&
 & SDEPN    ,SDEPX    ,SDEPD    ,SARGN    ,SARGX    ,&
 & SARGD    ,SSABN    ,SSABX    ,SSABD    ,SRSMN    ,&
 & SRSMX    ,SRSMD    ,SZZ0N    ,SZZ0M    ,SZZ0B    ,&
 & SZZ0U    ,SZZ0D  
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

LOGICAL           ,INTENT(IN)    :: LDNEW 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('VAL923',0,ZHOOK_HANDLE)

!  Threshold defining the mask
SMASK= 0.5_JPRB
!  Value for missing data + 1
SMANQ=-9998._JPRB
!  Land-use types for sea, ice-cap, desert, lakes
NTPMER= 1
NTPGLA= 2
NTPDES= 3
NTPLAC= 5
!  Roughness length : minimum, sea, sea-ice, urban areas, desert
SZZ0N= 0.001_JPRB
SZZ0M= 0.001_JPRB
SZZ0B= 0.001_JPRB
SZZ0U= 2.500_JPRB
SZZ0D= 0.001_JPRB
!  Ration of thermal to kinetic roughness length
STHER= 0.10_JPRB
!  Albedo : minimum, maximum, sea, ice-cap, sea-ice, desert
IF (LDNEW) THEN
  SALBN= 0.05_JPRB
  SALBX= 0.80_JPRB
ELSE
  SALBN= 0.07_JPRB
  SALBX= 0.70_JPRB
ENDIF
SALBM= 0.07_JPRB
SALBG= 0.75_JPRB
SALBB= 0.65_JPRB
SALBD= 0.10_JPRB
!  Emissivity : minimum, maximum, sea, ice-cap, sea-ice, desert
SEMIN= 0.90_JPRB
SEMIX= 1.00_JPRB
SEMIM= 0.96_JPRB
SEMIG= 0.98_JPRB
SEMIB= 0.97_JPRB
SEMID= 0.943_JPRB
!  Soil depth : minimum, maximum, desert
SDEPN= 0.10_JPRB
SDEPX= 8.00_JPRB
SDEPD= 0.10_JPRB
!  Percentage of clay : minimum, maximum, desert
SARGN=  3._JPRB
SARGX= 58._JPRB
SARGD=  3._JPRB
!  Percentage of sand : minimum, maximum, desert
SSABN=  6._JPRB
SSABX= 92._JPRB
SSABD= 92._JPRB
!  Minimum surface resistance : minimum, maximum, desert
SRSMX=5000._JPRB
SRSMN=   1.0_JPRB
SRSMD=5000._JPRB

WRITE(UNIT=NULOUT,FMT=111) SMASK,SMANQ,STHER,&
 & NTPMER,NTPGLA,NTPDES,NTPLAC  
WRITE(UNIT=NULOUT,FMT=112) SZZ0N,SZZ0M,SZZ0B,SZZ0U,SZZ0D
WRITE(UNIT=NULOUT,FMT=113) SALBN,SALBX,SALBM,SALBG,SALBB,SALBD,&
 & SEMIN,SEMIX,SEMIM,SEMIG,SEMIB,SEMID  
WRITE(UNIT=NULOUT,FMT=114) SDEPN,SDEPX,SDEPD,SARGN,SARGX,SARGD,&
 & SSABN,SSABX,SSABD,SRSMN,SRSMX,SRSMD  
111 FORMAT(' COMMON YOMCLI',/,&
 & ' SMASK=',F4.2,' SMANQ=',F6.0,' STHER=',F4.2,/&
 & ' NTPMER=',I2,' NTPGLA=',I2,' NTPDES=',I2,' NTPLAC=',I2)  
112 FORMAT(' LONGUEUR DE RUGOSITE :',/,&
 & ' minimum    mer    banquise  villes   desert ',&
 & /,5F9.3)  
113 FORMAT(' ALBEDO ET EMISSIVITE :',/,&
 & ' minimum  maximum    mer    glacier  banquise  desert ',&
 & 2(/,6F9.3))  
114 FORMAT(' PROFONDEUR, % ARGILE, % SABLE, RESIS. MIN. :',/,&
 & ' minimum  maximum   desert ',&
 & 4(/,3F9.3))  

IF (LHOOK) CALL DR_HOOK('VAL923',1,ZHOOK_HANDLE)
END SUBROUTINE VAL923
