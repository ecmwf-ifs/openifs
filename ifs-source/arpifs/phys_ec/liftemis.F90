! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LIFTEMIS (YDML_GCONF,KIDIA  , KFDIA , KLON, KLEV, KTRAC, KCHEM, KSTGLO, PDELP, PAPHIF, &
 &      PFLUX ,PTENC, KFLDX, KLEVX, PEXTRA) 

!* DESCRIPTION 
!     ----------
!   Injection of raised emissions, i.e. for volcanoes  
!   change tendency at prescribed level, set flux to zero 
!   for volcanos, introduce time dependent emissions 
!
!* INTERFACE.
!  ----------
!     *LIFTEMIS* is called from gems_init *.
! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV  :  Number of Levels
! KTRAC  :  Number of GFL tracer fields
! KSTGLO  :  BLOCK Index
! PAPHIF(KLON,KLEV)            : Geopotential on full levels 
! PCEN(KLON,KLEV,KTRAC)        : Mass mixing ratio   (kg/kg)
! PFLUX(KLON)                  : Surface emissions in kg/m2 

! OUTPUTS:
! -------
! PTENC  (KLON,KLEV,NCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS BECAUSE OF EMISSION(kg/kg s-1)
! PFLUX(KLON)                  : Surface emissions in kg/m2 
!
!
!     AUTHOR.
!     -------
!        JOHANNES FLEMMING  *ECMWF*
!        ORIGINAL : 2013-01-22

!     MODIFICATIONS.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE PARKIND1               , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMRIP0                , ONLY : NINDAT, NSSSSS
USE YOMCST                 , ONLY : RG, RDAY  
USE YOMLUN                 , ONLY : NULOUT
!USE YOMPHY2               , ONLY : TSPHY
USE YOMVOLCANO             , ONLY : SEMIVOCFLX, SEMIVOC, SEMIVOCFLXENS, LVOCMP, IVOCGP, &
&                                   IVOCBLK, SLVOC1, SLVOC2, NVOCDATES, IVOCSTART, SLVOCES1, SLVOCES2, LVOCENS
!USE YOMCHEM               , ONLY : YRCHEM, IEXTR_EM
!USE YOMCOMPO              , ONLY : YRCOMPO

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KLEV, KTRAC, KLEVX,KFLDX 
INTEGER(KIND=JPIM),INTENT(IN) :: KCHEM(YDML_GCONF%YGFL%NCHEM), KSTGLO
REAL(KIND=JPRB) ,INTENT(IN)   :: PDELP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB),INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB),INTENT(INOUT) :: PFLUX(KLON,KTRAC) 
REAL(KIND=JPRB),INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX)

!
!*       0.5   LOCAL VARIABLES
!              ---------------
INTEGER(KIND=JPIM) :: ID0, IM0, IY0, IDINCR,IDD,IMM,IYY,ILMON(12),IH0
INTEGER(KIND=JPIM) :: JK, JL, JT, JT1, IT, IVOL, IDATE 
INTEGER(KIND=JPIM) :: JTENS(10),ILVOC1, ILVOC2, IX(1)
REAL(KIND=JPRB) :: ZDELPVOC  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE  
CHARACTER(LEN=8) :: CLTEST     

!-----------------------------------------------------------------------

#include "fcttim.func.h"
#include "updcal.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LIFTEMIS',0,ZHOOK_HANDLE)
ASSOCIATE(NCHEM=>YDML_GCONF%YGFL%NCHEM, YCHEM=>YDML_GCONF%YGFL%YCHEM, &
 & RSTATI=>YDML_GCONF%YRRIP%RSTATI)
!-----------------------------------------------------------------------

! find index of volcanic SO2 in YCHEM
IVOL = -999_JPIM
JTENS(:) =  -999_JPIM
DO JT = 1, NCHEM
  IF (TRIM(YCHEM(JT)%CNAME) == 'SO2' ) IVOL = JT
  ! ensemble 
  IF (LVOCENS) THEN
    DO JT1 = 1, 10
      WRITE(CLTEST,'(A6,I2.2)') 'SO2_VE' , JT1  
      IF (TRIM(YCHEM(JT)%CNAME) == CLTEST ) JTENS(JT1) = JT
    ENDDO
  ENDIF
ENDDO

IF ( IVOL == -999_JPIM ) THEN 
  IF (LHOOK) CALL DR_HOOK('LIFTEMIS',1,ZHOOK_HANDLE )
  RETURN
ENDIF

!*       2.     CALCULATE TENDECIES  
!              -------------------

! volcanoe emissions on this process
IF ( LVOCMP .AND. NCHEM > 0  ) THEN
  IY0=NCCAA(NINDAT)
  IM0=NMM(NINDAT)
  ID0=NDD(NINDAT)
  IDINCR=(NSSSSS+NINT(RSTATI))/NINT(RDAY)
  IH0=INT(MOD(NSSSSS+NINT(RSTATI),NINT(RDAY))/3600_JPRB)
  CALL UPDCAL(ID0,IM0,IY0,IDINCR,IDD,IMM,IYY,ILMON,-1)
  IDATE=IYY*1000000 + IMM*10000 + IDD*100 + IH0  
  IT=-1

  ! find block with volcano 
  IF ( KSTGLO == IVOCBLK ) THEN 
  ! jlon of volcano
    JL =  IVOCGP

    ! Volcano Emissions with time dependency 
    ! find date in list
    DO JK=1,NVOCDATES
       IF (IDATE >= IVOCSTART(JK) ) THEN     
         IT=JK
       ENDIF 
    ENDDO
    ! time is right  
    IF (IT > 0) THEN  

      IX=MINLOC( ABS( PAPHIF(JL,1:KLEV)/RG - SLVOC1(IT)*1000.0_JPRB))
      ILVOC1=IX(1) 
      IX=MINLOC( ABS( PAPHIF(JL,1:KLEV)/RG - SLVOC2(IT)*1000.0_JPRB))
      ILVOC2=IX(1) 

      WRITE(NULOUT,*) ' VOLCANO EMISS HEIGHT LEVELS ', ILVOC1, ILVOC2, PAPHIF(JL, ILVOC1)/RG , PAPHIF(JL, ILVOC2)/RG 
      WRITE(NULOUT,*) ' VOLCANO EMISS FLUX and MASS  kg/ms ',  SEMIVOCFLX(IT) , SEMIVOC(IT),SLVOC1(IT),SLVOC2(IT) 
     ! calculate total detltap over injected levels 
      ZDELPVOC=0.0_JPRB
      DO JK = ILVOC1, ILVOC2
        ZDELPVOC = ZDELPVOC + PDELP(JL,JK)  
      ENDDO
      DO JK = ILVOC1, ILVOC2
         PTENC(JL,JK,KCHEM(IVOL)) = PTENC(JL,JK,KCHEM(IVOL)) + & 
      &  SEMIVOCFLX(IT) * RG / ZDELPVOC  
      ENDDO
     ! set surface fluc to zero 
      PFLUX(JL,KCHEM(IVOL)) = 0.0_JPRB       
!     IF (YRCOMPO%LCHEM_DIA) THEN
!        PEXTRA(JL,IVOL,IEXTR_EM)= PEXTRA(JL,IVOL,IEXTR_EM) + SEMIVOCFLX(IT) * TSPHY 
!     ENDIF
 
    ENDIF  ! emission time  

!   volcano ensemble always on   
    IF ( LVOCENS ) THEN
!    forecast time for puffs 
!    IZT=NINT(TSTEP*(REAL(KSTEP,JPRB)+0.5_JPRB))
      DO JT1 = 1, 10 
        IX=MINLOC( ABS( PAPHIF(JL,1:KLEV)/RG - SLVOCES1(JT1)*1000.0_JPRB))
        ILVOC1=IX(1) 
        IX=MINLOC( ABS( PAPHIF(JL,1:KLEV)/RG - SLVOCES2(JT1)*1000.0_JPRB))
        ILVOC2=IX(1) 
        WRITE(NULOUT,*) ' VOLCANO EMISS ENSEMBLE  LEVELS ', ILVOC1, ILVOC2, PAPHIF(JL, ILVOC1)/RG , PAPHIF(JL, ILVOC2)/RG 
        JT =  JTENS(JT1) 
        IF ( JT > 0 ) THEN 
          ZDELPVOC=0.0_JPRB
          DO JK = ILVOC1, ILVOC2
            ZDELPVOC = ZDELPVOC + PDELP(JL,JK)
          ENDDO
          DO JK = ILVOC1, ILVOC2
            PTENC(JL,JK,KCHEM(JT)) = PTENC(JL,JK,KCHEM(JT)) + &
      &     SEMIVOCFLXENS * RG / ZDELPVOC
          ENDDO
! the emissions appear as flux err, since no surface flux data 
!         IF (YRCOMPO%LCHEM_DIA) THEN
!           PEXTRA(JL,JT,IEXTR_EM)= PEXTRA(JL,JT,IEXTR_EM) + SEMIVOCFLXENS * TSPHY 
!        ENDIF
        ENDIF  
      ENDDO
    ENDIF ! LVOCENS  

  ENDIF
ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LIFTEMIS',1,ZHOOK_HANDLE )
END SUBROUTINE LIFTEMIS 
