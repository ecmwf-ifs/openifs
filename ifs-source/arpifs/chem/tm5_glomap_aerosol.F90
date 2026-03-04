! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_GLOMAP_AEROSOL ( KIDIA,KFDIA,KLON,KLEV, &
    &   PAERAOT, PAERAAOT, PAERASY,  &
    &   PTAUS_AER,PTAUA_AER, PMAER )

 
!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!
!
! assignment of aerosol optical depths 
! To be used for TM5 phosolysis at specific wavelengths.
! Interpolated from GLOMAP fields
!
!------------------------------------------------------------------
!
!
!**   INTERFACE.
!     ----------
!          *TM5_GLOMAP_AEROSOL* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA                       : Start of Array  
! KFDIA                       : End  of Array 
! KLON                        : Length of Arrays 
! KLEV                        : NUMBER OF LEVELS     
! PAERAOT(KLON,KLEV,6)        : Glomap extinction AOD per model level at 6 wavelengths
! PAERAAOT(KLON,KLEV,6)       : Glomap absorption AOD per model level at 6 wavelengths
! PAERASY(KLON,KLEV,6)        : Glomap asymetry factor
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
!      Original : 2018-08-31 Vincent Huijnen (KNMI)
!
!
!     MODIFICATIONS.
!     --------------
!        
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TM5_PHOTOLYSIS , ONLY : NBANDS_TROP, NGRID, WL_EFF, WL_GLOMAP


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA, KFDIA, KLON, KLEV

REAL(KIND=JPRB)   ,INTENT(IN) :: PAERAOT(KLON,KLEV,6)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAERAAOT(KLON,KLEV,6)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAERASY(KLON,KLEV,6)

REAL(KIND=JPRB), INTENT(OUT)   :: PTAUS_AER(KLON,KLEV,NBANDS_TROP,NGRID)
REAL(KIND=JPRB), INTENT(OUT)   :: PTAUA_AER(KLON,KLEV,NBANDS_TROP,NGRID)
REAL(KIND=JPRB), INTENT(OUT)   :: PMAER(KLON,KLEV,NBANDS_TROP,NGRID)



! * LOCAL 
INTEGER(KIND=JPIM) ::   JB, JK, JL

! REAL(KIND=JPRB)    :: ZFAC, ZF1,ZF2

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_GLOMAP_AEROSOL',0,ZHOOK_HANDLE )

! 6 GLOMAP wavelengths correspond to:  380, 440, 550, 670, 870 and 1020 nm 

DO JB = 1,6
  ! For these wavelengths don't interpolate the asymetry.
  !  just take its value at lowest wavelength.
  ! For AOD, assume AOD(lambda) ~ lambda^(-1), so:
  ! AOD(lambda) = AOD(lambda=380) *(380/lambda)
  DO JK=1,KLEV 
    DO JL=KIDIA,KFDIA 
      ! Fill Absorption aerosol layer
      PTAUA_AER(JL,JK,JB,1)=PAERAAOT(JL,JK,1)*WL_GLOMAP(1)/WL_EFF(JB)
      PTAUA_AER(JL,JK,JB,2)=PTAUA_AER(JL,JK,JB,1)
 
  
      ! Fill scattering aerosol layer: (note: scattering = extinction - absorption)
      PTAUS_AER(JL,JK,JB,1)= MAX(0._JPRB, PAERAOT(JL,JK,1)-PAERAAOT(JL,JK,1))*WL_GLOMAP(1)/WL_EFF(JB)
      PTAUS_AER(JL,JK,JB,2)= PTAUS_AER(JL,JK,JB,1)
     
      PMAER(JL,JK,JB,1)=MAX(MIN(PAERASY(JL,JK,1),1._JPRB),0._JPRB)
      PMAER(JL,JK,JB,2)=PMAER(JL,JK,JB,1)
    ENDDO
  ENDDO
ENDDO

JB=7
! For 580 nm just take closest value which is at 550 nm, which is index 3
DO JK=1,KLEV 
  DO JL=KIDIA,KFDIA 
    ! Fill Absorption aerosol layer
    PTAUA_AER(JL,JK,JB,1)=PAERAAOT(JL,JK,3)
    PTAUA_AER(JL,JK,JB,2)=PTAUA_AER(JL,JK,JB,1)
 
 
    ! Fill scattering aerosol layer: (note: scattering = extinction - absorption)
    PTAUS_AER(JL,JK,JB,1)= MAX(0._JPRB, PAERAOT(JL,JK,3)-PAERAAOT(JL,JK,3))
    PTAUS_AER(JL,JK,JB,2)= PTAUS_AER(JL,JK,JB,1)
   
    PMAER(JL,JK,JB,1)=MAX(MIN(PAERASY(JL,JK,3),1._JPRB),0._JPRB)
    PMAER(JL,JK,JB,2)=PMAER(JL,JK,JB,1)
  ENDDO
ENDDO

        
! ! Code to interpolate for specified wavelengths . Not used
! DO JB = 7 ,7  
!   ZF1=(WL_EFF(JB)-WL_GLOMAP(1)) / (WL_GLOMAP(2)-WL_GLOMAP(1))
!   ! ensure validity
!   ZF1=MAX(0._JPRB,MIN(1._JPRB,ZF1))
!   ZF2=1.-ZF1
!   DO JK=1,KLEV
!     DO JL = KIDIA,KFDIA
!       PTAUA_AER(JL,JK,JB,1)=MAX(0._JPRB,ZF2*PAERAAOT(JL,JK,1)+ ZF1*PAERAAOT(JL,JK,2))
!       PTAUA_AER(JL,JK,JB,2)=PTAUA_AER(JL,JK,JB,1)
!     
!       PTAUS_AER(JL,JK,JB,1)=MAX(0._JPRB,ZF2*(PAERAOT(JL,JK,1)-PAERAAOT(JL,JK,1))+ ZF1*(PAERAOT(JL,JK,2)-PAERAAOT(JL,JK,2)))
!       PTAUS_AER(JL,JK,JB,2)=PTAUS_AER(JL,JK,JB,1)
!       IF (PTAUS_AER(JL,JK,JB,1) > 0._JPRB ) THEN
!         PMAER(JL,JK,JB,1)=MAX(MIN(ZF2*PAERASY(JL,JK,1)+ZF1*PAERASY(JL,JK,2),1._JPRB),0._JPRB)
!         PMAER(JL,JK,JB,2)=PMAER(JL,JK,JB,1)
!       ELSE
!         PMAER(JL,JK,JB,1)=0._JPRB
!         PMAER(JL,JK,JB,2)=0._JPRB
!       ENDIF
!     ENDDO
!   ENDDO
! ENDDO
  

IF (LHOOK) CALL DR_HOOK('TM5_GLOMAP_AEROSOL',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_GLOMAP_AEROSOL

