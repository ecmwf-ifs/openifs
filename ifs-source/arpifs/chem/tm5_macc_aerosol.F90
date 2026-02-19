! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_MACC_AEROSOL ( KIDIA,KFDIA,KLON,KLEV, KACTAERO, &
    &   PAPH   , PAEROK  , PRHCL   ,   &
    &   PTAUS_AER,PTAUA_AER, PMAER )

 
!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!
!
! assignment of aerosol optical depths 
! To be used for TM5 phosolysis at specific wavelengths.
! Interpolated from MACC fields
!
!------------------------------------------------------------------
!
!
!**   INTERFACE.
!     ----------
!          *TM5_MACC_AEROSOL* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV  :  NMBER OF LEVELS         (INPUT)
!
! OUTPUTS:
! -------
!
! PMAER
! PTAUS_AER
! PTAUA_AER
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!
!        Vincent Huijnen     *KNMI*
!         
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL :       2012-07-25
!        Vincent Huijnen: 2017-08-30 : Merged with computation of aerosol optical prop's


USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TM5_PHOTOLYSIS , ONLY : NBANDS_TROP, NGRID, WL_EFF, WL_AER

USE YOMCST    ,ONLY : RG
USE YOEAERATM ,ONLY : YREAERATM
USE YOEAEROP  ,ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA, &
                  &  ASY_SU, ASY_OM, ASY_DD, ASY_SS, ASY_BC, ASY_NI, ASY_AM, ASY_SOA, &
                  &  OMG_SU, OMG_OM, OMG_DD, OMG_SS, OMG_BC, OMG_NI, OMG_AM, OMG_SOA
USE YOEAERSNK ,ONLY : YREAERSNK
USE YOEAERVOL ,ONLY : YREAERVOL

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA, KFDIA, KLON, KLEV
INTEGER(KIND=JPIM),INTENT(IN)  :: KACTAERO  ! Number of active aerosol species

REAL(KIND=JPRB)   ,INTENT(IN)  :: PAPH(KLON,0:KLEV)                        
REAL(KIND=JPRB)   ,INTENT(IN)  :: PAEROK(KLON,KLEV,KACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PRHCL(KLON,KLEV)                  

REAL(KIND=JPRB), INTENT(OUT)   :: PTAUS_AER(KLON,KLEV,NBANDS_TROP,NGRID)
REAL(KIND=JPRB), INTENT(OUT)   :: PTAUA_AER(KLON,KLEV,NBANDS_TROP,NGRID)
REAL(KIND=JPRB), INTENT(OUT)   :: PMAER(KLON,KLEV,NBANDS_TROP,NGRID)



! * LOCAL 
REAL(KIND=JPRB)   :: ZAERMSS(KLON,KLEV,KACTAERO)
REAL(KIND=JPRB)   :: ZAERTAU

INTEGER(KIND=JPIM) :: ITWAVL(2)

INTEGER(KIND=JPIM) :: IBIN, IEFRH, IIRH, ITYP, IWAVL
INTEGER(KIND=JPIM) :: JAER, JTAB, JWAVL, JB, JK, JL
INTEGER(KIND=JPIM) :: IRH(KLON,KLEV)

REAL(KIND=JPRB) :: ZALF(KACTAERO), ZASY(KACTAERO), ZOMG(KACTAERO), ZETA(KLON,KLEV)
REAL(KIND=JPRB) :: ZAEROMGLT(KLON,KLEV,2)

REAL(KIND=JPRB) :: ZFAC

REAL(KIND=JPRB)    :: ZF1,ZF2

REAL(KIND=JPRB)  :: ZTAUT(KLON,KLEV,2)   ! Aerosol optical depth at 2 wavelengths
REAL(KIND=JPRB)  :: ZTAUB(KLON,KLEV,2)   ! Aerosol absorption optical depth at 2 wavelengths
REAL(KIND=JPRB)  :: ZASYL(KLON,KLEV,2)   ! Aerosol assymetry factor at 2 wavelengths.

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_MACC_AEROSOL',0,ZHOOK_HANDLE )


ASSOCIATE( &
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & RRHTAB=>YREAERSNK%RRHTAB, RSSGROWTH_RHTAB=>YREAERSNK%RSSGROWTH_RHTAB, &
 & RSSDENS_RHTAB=>YREAERSNK%RSSDENS_RHTAB, &
 & YAERO_DESC=>YREAERATM%YAERO_DESC, &
 & NVOLOPTP=>YREAERVOL%NVOLOPTP, RMMD_SS=>YREAERSNK%RMMD_SS)


!-- units:
!   ------
! ZALF is the extinction coefficient           m2 g-1
! PAEROK is the aerosol mass mixing ratio      kg kg-1 

IRH(:,:)=1
ZAERMSS(KIDIA:KFDIA,1:KLEV,1:KACTAERO)= 0.0_JPRB
ZTAUT(KIDIA:KFDIA,1:KLEV,:)=0._JPRB
ZTAUB(KIDIA:KFDIA,1:KLEV,:)=0._JPRB
ZASYL(KIDIA:KFDIA,1:KLEV,:)=0._JPRB
ZAEROMGLT(KIDIA:KFDIA,1:KLEV,:)=0._JPRB

!-- the effective relative humidity is the low value (20%) assumed for hydrophobic component of OM
IEFRH=3

!-- define RH index from "clear-sky" (not yet!) relative humidity
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZETA(JL,JK)= (PAPH(JL,JK-1)+PAPH(JL,JK))/(2._JPRB*PAPH(JL,KLEV))
      DO JTAB=1,12
        IF (PRHCL(JL,JK)*100._JPRB > RRHTAB(JTAB)) THEN
          IRH(JL,JK)=JTAB
        ENDIF
      ENDDO
    ENDDO
  ENDDO

!-- mass of aerosols in each layer

  DO JAER=1,KACTAERO
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
!-- PAEROK in kg kg-1 ; ZAERMSS in kg m-2
        ZAERMSS(JL,JK,JAER) = PAEROK(JL,JK,JAER)*(PAPH(JL,JK)-PAPH(JL,JK-1))/RG 
      ENDDO
    ENDDO
  ENDDO




  ITWAVL(1)= 1  ! 340 nm  
  ITWAVL(2)= 10 ! 645 nm  
      
  DO JWAVL=1,2
    IWAVL=ITWAVL(JWAVL)

    DO JAER=1,KACTAERO
      ITYP=YAERO_DESC(JAER)%NTYP
      IBIN=YAERO_DESC(JAER)%NBIN

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:FA,   7:BS,   8=VS,
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1:FA,   1:SB    1=VS,
!   N.B.: extinction coefficients are in m2 g-1

      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          IIRH=IRH(JL,JK)
          ZFAC = 1.0_JPRB
          IF (ITYP == 1) THEN
            ZALF(JAER)=ALF_SS( IIRH ,IWAVL,IBIN)
            ZOMG(JAER)=OMG_SS( IIRH ,IWAVL,IBIN)
            ZASY(JAER)=ASY_SS( IIRH ,IWAVL,IBIN)
            ZFAC = RSS_RH80_MASSFAC
          ELSEIF (ITYP == 2) THEN
            ZALF(JAER)=ALF_DD(IBIN,IWAVL)
            ZOMG(JAER)=OMG_DD(IBIN,IWAVL)
            ZASY(JAER)=ASY_DD(IBIN,IWAVL)
          ELSEIF (ITYP == 3) THEN
!-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
            IF (IBIN == 2) IIRH=IEFRH
            ZALF(JAER)=ALF_OM( IIRH ,IWAVL)
            ZOMG(JAER)=OMG_OM( IIRH ,IWAVL)
            ZASY(JAER)=ASY_OM( IIRH ,IWAVL)
          ELSEIF (ITYP == 4) THEN
            ZALF(JAER)=ALF_BC(IWAVL)
            ZOMG(JAER)=OMG_BC(IWAVL)
            ZASY(JAER)=ASY_BC(IWAVL)
          ELSEIF (ITYP == 5 .OR. ITYP == 9) THEN
            ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
            ZOMG(JAER)=OMG_SU( IIRH ,IWAVL)
            ZASY(JAER)=ASY_SU( IIRH ,IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
            IF (IBIN == 2) THEN
              ZALF(JAER)=0._JPRB
              ZOMG(JAER)=0._JPRB
              ZASY(JAER)=0._JPRB
            ENDIF
          ELSEIF (ITYP == 6) THEN
            ZALF(JAER)=ALF_NI(IIRH ,IWAVL,IBIN)
            ZOMG(JAER)=OMG_NI(IIRH ,IWAVL,IBIN)
            ZASY(JAER)=ASY_NI(IIRH ,IWAVL,IBIN)
          ELSEIF (ITYP == 7) THEN
            ZALF(JAER)=ALF_AM(IIRH ,IWAVL)
            ZOMG(JAER)=OMG_AM(IIRH ,IWAVL)
            ZASY(JAER)=ASY_AM(IIRH ,IWAVL)
          ELSEIF (ITYP == 8) THEN
            IF (IBIN < 3) THEN
              ZALF(JAER)=ALF_SOA(IIRH,IWAVL,IBIN)
              ZOMG(JAER)=OMG_SOA(IIRH,IWAVL,IBIN)
              ZASY(JAER)=ASY_SOA(IIRH,IWAVL,IBIN)
            ELSE
              ZALF(JAER)=0._JPRB
              ZOMG(JAER)=0._JPRB
              ZASY(JAER)=0._JPRB
            ENDIF
          ELSEIF (ITYP == 9) THEN
            IF (NVOLOPTP == 1) THEN
!-- use sulphate optical properties at 20% RH
              IIRH=IEFRH
              ZALF(JAER)=ALF_SU( IIRH ,IWAVL)
              ZOMG(JAER)=OMG_SU( IIRH ,IWAVL)
              ZASY(JAER)=ASY_SU( IIRH ,IWAVL)
            ELSEIF (NVOLOPTP == 2) THEN
!-- use black carbon optical properties
              ZALF(JAER)=ALF_BC(IWAVL)
              ZOMG(JAER)=OMG_BC(IWAVL)
              ZASY(JAER)=ASY_BC(IWAVL)
            ELSEIF (NVOLOPTP == 3) THEN
!-- use dust for 0.9-20 um bin
              ZALF(JAER)=ALF_DD(3,IWAVL)
              ZOMG(JAER)=OMG_DD(3,IWAVL)
              ZASY(JAER)=ASY_DD(3,IWAVL)
            ENDIF
          ELSEIF (ITYP == 10) THEN
            ZALF(JAER)=0._JPRB
            ZOMG(JAER)=0._JPRB
            ZASY(JAER)=0._JPRB
          ENDIF


!- ZAERTAU (ND = g m-2 * m2 g-1)
!- For SS, vary depending on RH
          ZAERTAU  = ZAERMSS(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER)

          ZTAUT(JL,JK,JWAVL) = ZTAUT(JL,JK,JWAVL) + ZAERTAU ! AOD @ 645 nm
          ZTAUB(JL,JK,JWAVL) = ZTAUB(JL,JK,JWAVL) + ZAERTAU*(1._JPRB-ZOMG(JAER))*(1._JPRB-ZETA(JL,JK))
          ZAEROMGLT(JL,JK,JWAVL) = ZAEROMGLT(JL,JK,JWAVL) + ZAERTAU*ZOMG(JAER)
          ZASYL (JL,JK,JWAVL) = ZASYL (JL,JK,JWAVL) + ZAERTAU*ZOMG(JAER)*ZASY(JAER)
        ENDDO
      ENDDO



    ENDDO

  ENDDO

! Compute Asymmetry factor at selected wavelengths:
  DO JL=KIDIA,KFDIA
    DO JK =1,KLEV
      DO JWAVL=1,2 
        IF( ZTAUT(JL,JK,JWAVL) > 0._JPRB .AND. ZAEROMGLT(JL,JK,JWAVL) > 0._JPRB) THEN
          ZASYL(JL,JK,JWAVL)=ZASYL(JL,JK,JWAVL)/ZAEROMGLT(JL,JK,JWAVL)
        ELSE
          ZASYL(JL,JK,JWAVL)=0._JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDDO




DO JB = 1,5
  ! For these wavelengths don't interpolate the asymetry.
  !  just take its value at lowest wavelength.
  ! For AOD, assume AOD(lambda) ~ lambda^(-1), so:
  ! AOD(lambda) = AOD(lambda=340) *(340/lambda)
  DO JK=1,KLEV 
    DO JL=KIDIA,KFDIA 
     ! Fill Absorption aerosol layer
     PTAUA_AER(JL,JK,JB,1)=MAX(0._JPRB,ZTAUB(JL,JK,1)*WL_AER(1)/WL_EFF(JB))
     PTAUA_AER(JL,JK,JB,2)=PTAUA_AER(JL,JK,JB,1)

 
     ! Fill scattering aerosol layer: scattering = total - absorption
     PTAUS_AER(JL,JK,JB,1)=MAX(0._JPRB,(ZTAUT(JL,JK,1)-ZTAUB(JL,JK,1))*WL_AER(1)/WL_EFF(JB))
     PTAUS_AER(JL,JK,JB,2)= PTAUS_AER(JL,JK,JB,1)
    
     PMAER(JL,JK,JB,1)=MAX(MIN(ZASYL(JL,JK,1),1._JPRB),0._JPRB)
     PMAER(JL,JK,JB,2)=PMAER(JL,JK,JB,1)
    ENDDO
  ENDDO
ENDDO
        
! Interpolate for the two largest wavelengths.
DO JB = 6 ,7  
  ZF1=(WL_EFF(JB)-WL_AER(1)) / (WL_AER(2)-WL_AER(1))
  ! ensure validity
  ZF1=MAX(0._JPRB,MIN(1._JPRB,ZF1))
  ZF2=1.-ZF1
  DO JK=1,KLEV
    DO JL = KIDIA,KFDIA
      PTAUA_AER(JL,JK,JB,1)=MAX(0._JPRB,ZF2*ZTAUB(JL,JK,1)+ ZF1*ZTAUB(JL,JK,2))
      PTAUA_AER(JL,JK,JB,2)=PTAUA_AER(JL,JK,JB,1)
    
      PTAUS_AER(JL,JK,JB,1)=MAX(0._JPRB,ZF2*(ZTAUT(JL,JK,1)-ZTAUB(JL,JK,1))+ ZF1*(ZTAUT(JL,JK,2)-ZTAUB(JL,JK,2)))
      PTAUS_AER(JL,JK,JB,2)=PTAUS_AER(JL,JK,JB,1)
      IF (PTAUS_AER(JL,JK,JB,1) > 0._JPRB ) THEN
        PMAER(JL,JK,JB,1)=MAX(MIN(ZF2*ZASYL(JL,JK,1)+ZF1*ZASYL(JL,JK,2),1._JPRB),0._JPRB)
        PMAER(JL,JK,JB,2)=PMAER(JL,JK,JB,1)
      ELSE
        PMAER(JL,JK,JB,1)=0._JPRB
        PMAER(JL,JK,JB,2)=0._JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO
  

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TM5_MACC_AEROSOL',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_MACC_AEROSOL

