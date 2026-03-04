! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_RAD &
 & ( YDEAERATM,YDML_PHY_AER,YDCOMPO,YGFL,KIDIA  , KFDIA  , KLON   , KTDIA, KLEV , KACTAERO, KNWAVL, KTWAVL, KSTART, KSTEP, &
 &   PAP    , PAPH   , PQ     , PT   , PAERO, &
 &   PTAUAER, POMGAER, PASYAER &
 & )

!*** * AER_RAD* - COMPUTATION OF AEROSOL OPTICAL THICKNESSES FOR RADIATION

!**   INTERFACE.
!     ----------
!          *AER_RAD* IS CALLED FROM *CALLPAR*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2008-05-27
!        SRemy    : 2017-04-21  add sea-salt mmr switch (80% RH or dry)

!     PURPOSE.
!     --------
!     - computes absorption optical thicknesses at 550 nm for the 7 aerosol
!       types recognized by the UV-radiation processor. 
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_AEROSOL_MOD , ONLY : MODEL_PHYSICS_AEROSOL_TYPE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE YOEAEROP , ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA, &
                  &   ASY_SU, ASY_OM, ASY_DD, ASY_SS, ASY_BC, ASY_NI, ASY_AM, ASY_SOA, &
                  &   OMG_SU, OMG_OM, OMG_DD, OMG_SS, OMG_BC, OMG_NI, OMG_AM, OMG_SOA 
USE YOEAERATM, ONLY : TEAERATM
USE YOMCOMPO , ONLY : TCOMPO
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOEPHLI  , ONLY : YREPHLI
IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERATM)    ,INTENT(INOUT):: YDEAERATM
TYPE(MODEL_PHYSICS_AEROSOL_TYPE), INTENT(INOUT) :: YDML_PHY_AER
TYPE(TCOMPO)      ,INTENT(INOUT):: YDCOMPO
TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)  :: KTDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)  :: KACTAERO
INTEGER(KIND=JPIM),INTENT(IN)  :: KSTART, KSTEP
INTEGER(KIND=JPIM),INTENT(IN)  :: KNWAVL, KTWAVL(19)

REAL(KIND=JPRB)   ,INTENT(IN)  :: PAP(KLON,KLEV), PAPH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PQ(KLON,KLEV) , PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PAERO(KLON,KLEV,YGFL%NAERO)

REAL(KIND=JPRB)   ,INTENT(OUT) :: PTAUAER(KLON,8,KLEV), POMGAER(KLON,8,KLEV), PASYAER(KLON,8,KLEV)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: IBIN, IFLAG, ITYP, IWAVL,IAE8,IEFRH,IIRH
INTEGER(KIND=JPIM) :: JAER, JK, JL, JTAB, JAE7
INTEGER(KIND=JPIM) :: IRH(KLON,KLEV)
REAL(KIND=JPRB) :: ZALF(KACTAERO), &
  & ZASY(KACTAERO), ZOMG(KACTAERO), ZQSAT(KLON,KLEV), ZRHCL(KLON,KLEV)
REAL(KIND=JPRB) :: ZAERMSS(KLON,KLEV,KACTAERO), ZAERTAU(KLON,KLEV,KACTAERO), &
  & ZAEROMG(KLON,KLEV,KACTAERO), ZAERASY(KLON,KLEV,KACTAERO)
REAL(KIND=JPRB) :: ZEPSTAU, ZFAC

LOGICAL :: LLPHYLIN
!-----------------------------------------------------------------------

#include "satur.intfb.h"

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_RAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDEAERSNK=>YDML_PHY_AER%YREAERSNK, &
 & YDEAERVOL=>YDML_PHY_AER%YREAERVOL)
ASSOCIATE(NAERO=>YGFL%NAERO, &
 & RRHTAB=>YDEAERSNK%RRHTAB, &
 & RSS_RH80_MASSFAC=>YDEAERATM%RSS_RH80_MASSFAC, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME, &
 & LAERNITRATE => YDCOMPO%LAERNITRATE, &
 & NVOLOPTP=>YDEAERVOL%NVOLOPTP, &
 & YAERO_DESC=>YDEAERATM%YAERO_DESC)

ZEPSTAU=1.E-20_JPRB

LLPHYLIN = YREPHLI%LPHYLIN

!-- units:
!   ------
! ZALF is the extinction coefficient           m2 g-1
! PAERO is the aerosol mass mixing ratio      kg kg-1 

!-- the effective relative hunidity is the low value (20%) assumed for
!hydrophobic component of OM
IEFRH=3


IRH(KIDIA:KFDIA,1:KLEV)=1
PTAUAER(KIDIA:KFDIA,:,1:KLEV)=0._JPRB
POMGAER(KIDIA:KFDIA,:,1:KLEV)=0._JPRB
PASYAER(KIDIA:KFDIA,:,1:KLEV)=0._JPRB

IFLAG=2
CALL SATUR (KIDIA, KFDIA, KLON  , KTDIA , KLEV, LLPHYLIN, &
  & PAP  , PT    , ZQSAT, IFLAG)  

!-- define RH index from "clear-sky" (not yet!) relative humidity
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHCL(JL,JK)= PQ(JL,JK) / ZQSAT(JL,JK)
    DO JTAB=1,12
      IF (ZRHCL(JL,JK)*100._JPRB > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO

!-- mass of aerosols in each layer

DO JAER=1,KACTAERO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAERMSS(JL,JK,JAER) = MAX(ZEPSTAU, PAERO(JL,JK,JAER))*(PAPH(JL,JK)-PAPH(JL,JK-1))/RG 
    ENDDO
  ENDDO
ENDDO
      
!-- IWAVL = 9 corresponds to the reference wavelength (0.55 um) for this way of handling 
!  radiation transfer for aerosols
  IWAVL=9

  DO JAER=1,KACTAERO
    SELECT CASE (TRIM(AERO_SCHEME))

      CASE ("glomap")

        ! TODO: GLOMAP needs proper radiative code, not just mapping onto
        !   the AER types!

        IF (LAERNITRATE) THEN
          IF (JAER == 16) THEN
            ITYP=1
            IBIN=1
          ELSEIF (JAER == 24) THEN
            ITYP=1
            IBIN=3
          ELSEIF (JAER == 17 .OR. JAER == 30) THEN
            ITYP=2
            IBIN=1
          ELSEIF (JAER == 25 .OR. JAER == 32) THEN
            ITYP=2
            IBIN=3
          ELSEIF (JAER == 9 .OR. JAER == 15 .OR. JAER==23 .OR. JAER==28) THEN
            ITYP=3
            IBIN=1
          ELSEIF (JAER == 8 .OR. JAER == 14 .OR. JAER == 22 .OR. JAER == 27) THEN
            ITYP=4
            IBIN=1
          ELSEIF (JAER == 5 .OR. JAER == 11 .OR. JAER==19) THEN
            ITYP=5
            IBIN=1
          ELSEIF (JAER == 6 .OR. JAER == 12 .OR. JAER==20) THEN
            ITYP=5
            IBIN=1
          ELSEIF (JAER == 7 .OR. JAER == 13 .OR. JAER==21) THEN
            ITYP=5
            IBIN=1
          ELSEIF (JAER==2 .OR. JAER==3) THEN
          ! Nucleation mode
            EXIT
          ELSEIF (JAER==1 .OR. JAER==4 .OR. JAER==10 .OR. JAER==18 .OR. JAER==26 .OR. JAER==29 .OR. JAER==31) THEN
            ! Number tracer
            EXIT
          ELSEIF (JAER>=33 .AND. JAER<=36) THEN
            ! Precursor
            EXIT
          ELSEIF (JAER>=38 .AND. JAER<=41) THEN
            ! Water content
            EXIT
          ENDIF
        ELSE
          IF (JAER == 12) THEN
            ITYP=1
            IBIN=1
          ELSEIF (JAER == 18) THEN
            ITYP=1
            IBIN=3
          ELSEIF (JAER == 13 .OR. JAER == 24) THEN
            ITYP=2
            IBIN=1
          ELSEIF (JAER == 19 .OR. JAER == 26) THEN
            ITYP=2
            IBIN=3
          ELSEIF (JAER == 7 .OR. JAER == 11 .OR. JAER==17 .OR. JAER==22) THEN
            ITYP=3
            IBIN=1
          ELSEIF (JAER == 6 .OR. JAER == 10 .OR. JAER == 16 .OR. JAER == 21) THEN
            ITYP=4
            IBIN=1
          ELSEIF (JAER == 5 .OR. JAER == 9 .OR. JAER==15) THEN
            ITYP=5
            IBIN=1
          ELSEIF (JAER==2 .OR. JAER==3) THEN
          ! Nucleation mode
            EXIT
          ELSEIF (JAER==1 .OR. JAER==4 .OR. JAER==8 .OR. JAER==14 .OR. JAER==20 .OR. JAER==23 .OR. JAER==25) THEN
            ! Number tracer
            EXIT
          ELSEIF (JAER>=27 .AND. JAER<=30) THEN
            ! Precursor
            EXIT
          ELSEIF (JAER>=32 .AND. JAER<=35) THEN
            ! Water content
            EXIT
          ENDIF

        ENDIF

      CASE ("aer")

        ITYP=YAERO_DESC(JAER)%NTYP
        IBIN=YAERO_DESC(JAER)%NBIN

    END SELECT

    IF (ITYP <= 8) THEN
      IAE8 = ITYP
    ELSEIF (ITYP >= 9) THEN
      IAE8 = 4
    ELSEIF (ITYP == 10) THEN
      IAE8 = 5
    ENDIF

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:NI, !7:AM, 8: SOA, 9: VAsh, 10: Vsu
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1-2:NI,   1:AM, 1-3: SOA, 1: VAsh, 1-2: Vsu
!   N.B.: extinction coefficients are in m2 g-1

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        IIRH=IRH(JL,JK)
        ZFAC = 1.0_JPRB
        IF (ITYP == 1) THEN
          ZALF(JAER)=ALF_SS(IIRH,IWAVL,IBIN)
          ZOMG(JAER)=OMG_SS(IIRH,IWAVL,IBIN)
          ZASY(JAER)=ASY_SS(IIRH,IWAVL,IBIN)
          ZFAC = RSS_RH80_MASSFAC
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALF_DD(IBIN,IWAVL)
          ZOMG(JAER)=OMG_DD(IBIN,IWAVL)
          ZASY(JAER)=ASY_DD(IBIN,IWAVL)
        ELSEIF (ITYP == 3) THEN
          IF (IBIN == 2) IIRH=IEFRH
          ZALF(JAER)=ALF_OM(IIRH,IWAVL)
          ZOMG(JAER)=OMG_OM(IIRH,IWAVL)
          ZASY(JAER)=ASY_OM(IIRH,IWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)= ALF_BC(IWAVL)
          ZOMG(JAER)= OMG_BC(IWAVL)
          ZASY(JAER)=ASY_BC(IWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
!          ZALF(JAER)=ALF_SU(IIRH,IWAVL)/RMSO4
          ZALF(JAER)= ALF_SU(IIRH,IWAVL)
          ZOMG(JAER)= OMG_SU(IIRH,IWAVL)
          ZASY(JAER)=ASY_SU(IIRH,IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
            ZOMG(JAER)=0._JPRB
            ZASY(JAER)=0._JPRB
          ENDIF
        ELSEIF (ITYP == 6) THEN
          ZALF(JAER)= ALF_NI(IIRH,IWAVL,IBIN)
          ZOMG(JAER)= OMG_NI(IIRH,IWAVL,IBIN)
          ZASY(JAER)= ASY_NI(IIRH,IWAVL,IBIN)
        ELSEIF (ITYP == 7) THEN
          ZALF(JAER)= ALF_AM(IIRH,IWAVL)
          ZOMG(JAER)= OMG_AM(IIRH,IWAVL)
          ZASY(JAER)= ASY_AM(IIRH,IWAVL)
        ELSEIF (ITYP == 8) THEN
          IF (IBIN < 3) THEN
            ZALF(JAER)= ALF_SOA(IIRH,IWAVL,IBIN)
            ZOMG(JAER)= OMG_SOA(IIRH,IWAVL,IBIN)
            ZASY(JAER)= ASY_SOA(IIRH,IWAVL,IBIN)
          ELSE
            ZALF(JAER)= 0._JPRB
            ZOMG(JAER)= 0._JPRB
            ZASY(JAER)= 0._JPRB
          ENDIF
        ELSEIF (ITYP ==9) THEN
          IF (NVOLOPTP == 1) THEN
            IIRH=IEFRH
            ZALF(JAER)= ALF_SU(IIRH, IWAVL)
            ZOMG(JAER)= OMG_SU(IIRH, IWAVL)
            ZASY(JAER)=ASY_SU(IIRH, IWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
            ZALF(JAER)= ALF_BC(IWAVL)
            ZOMG(JAER)= OMG_BC(IWAVL)
            ZASY(JAER)=ASY_BC(IWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
            ZALF(JAER)= ALF_DD(3,IWAVL)
            ZOMG(JAER)= OMG_DD(3,IWAVL)
            ZASY(JAER)=ASY_DD(3,IWAVL)
          ENDIF
        ENDIF

!- ZAERTAU (ND = g m-2 * m2 g-1)
        ZAERTAU(JL,JK,JAER) = ZAERMSS(JL,JK,JAER) * 1.E+03_JPRB * ZALF(JAER)*ZFAC
        ZAEROMG(JL,JK,JAER) = ZAERTAU(JL,JK,JAER)*ZOMG(JAER)
        ZAERASY(JL,JK,JAER) = ZAEROMG(JL,JK,JAER)*ZASY(JAER)
      ENDDO
    ENDDO
!--end loop on JK

!-- total optical thickness per aerosol type: as sum over all individual components
    ITYP=YAERO_DESC(JAER)%NTYP
    IBIN=YAERO_DESC(JAER)%NBIN
    IF (ITYP <= 8) THEN
      IAE8 = ITYP
    ELSEIF (ITYP == 10) THEN
      IAE8 = 5
    ELSEIF (ITYP == 9) THEN
      IAE8 = 4
    ENDIF

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PTAUAER(JL,IAE8,JK)=PTAUAER(JL,IAE8,JK)+ZAERTAU(JL,JK,JAER)
      POMGAER(JL,IAE8,JK)=POMGAER(JL,IAE8,JK)+ZAEROMG(JL,JK,JAER)
      PASYAER(JL,IAE8,JK)=PASYAER(JL,IAE8,JK)+ZAERASY(JL,JK,JAER)
    ENDDO
  ENDDO
ENDDO
!--end loop on JAER

DO JAE7=1,7
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF( PTAUAER(JL,JAE7,JK) > 0._JPRB .AND. POMGAER(JL,JAE7,JK) > 0._JPRB) THEN
        PASYAER(JL,JAE7,JK)=PASYAER(JL,JAE7,JK)/POMGAER(JL,JAE7,JK)
        POMGAER(JL,JAE7,JK)=POMGAER(JL,JAE7,JK)/PTAUAER(JL,JAE7,JK)
      ELSE
        PASYAER(JL,JAE7,JK)=0._JPRB
        POMGAER(JL,JAE7,JK)=0._JPRB
        PTAUAER(JL,JAE7,JK)=0._JPRB
      ENDIF
    ENDDO

  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_RAD',1,ZHOOK_HANDLE)
END SUBROUTINE AER_RAD
