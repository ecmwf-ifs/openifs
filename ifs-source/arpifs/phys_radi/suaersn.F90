! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUAERSN (YDESWRT, KTSW, KSW)

!**** *SUAERS*   - INITIALIZE COMMON YOEAER

!     PURPOSE.
!     --------
!           INITIALIZE YOEAER, THE COMMON THAT CONTAINS THE
!           RADIATIVE CHARACTERISTICS OF THE AEROSOLS

!**   INTERFACE.
!     ----------
!              -----        -----

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOEAER

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "IFS MODEL"

!     AUTHOR.
!     -------
!      JEAN-JACQUES MORCRETTE *ECMWF*
!      ORIGINAL : 88-02-15

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOESW    , ONLY : TESWRT

!      ----------------------------------------------------------------

IMPLICIT NONE

TYPE(TESWRT)      ,INTENT(INOUT) :: YDESWRT
INTEGER(KIND=JPIM),INTENT(IN)    :: KTSW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW 
REAL(KIND=JPRB) :: ZTAUA2(2,6)  ,ZPIZA2(2,6)  ,ZCGA2(2,6)
REAL(KIND=JPRB) :: ZTAUA4(4,6)  ,ZPIZA4(4,6)  ,ZCGA4(4,6)
REAL(KIND=JPRB) :: ZTAUA6(6,6)  ,ZPIZA6(6,6)  ,ZCGA6(6,6)

INTEGER(KIND=JPIM) :: JAER, JNU
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

!*       1.    SHORTWAVE COEFFICIENTS
!              ----------------------
!=======================================================================
!-- The (old) five aerosol types were respectively:

!  1/ continental average (+desert)       2/ maritime
!  3/ urban                               4/ volcanic active
!  5/ stratospheric background

!-- old values were not spectrally defined:
! ZTAU2  = .730719, .912819, .725059, .745405, .682188
! ZPIZA2 = .872212, .982545, .623143, .944887, .997975
! ZCGA2  = .647596, .739002, .580845, .662657, .624246
!=======================================================================

!-- The six aerosol types are respectively:

!  1/ continental average                 2/ maritime
!  3/ desert                              4/ urban
!  5/ volcanic active                     6/ stratospheric background

! The quantities given are:
! TAU : ratio of average optical thickness in interval to that at 0.55 
!       micron
! PIZA: average single scattering albedo
! CGA : average asymmetry factor

! computed from Hess and Koepke (con, mar, des, urb)
!          from Bonnel et al.   (vol, str)

!        1.1   TWO SPECTRAL INTERVALS (0.25-0.69-4.00microns)

IF (LHOOK) CALL DR_HOOK('SUAERSN',0,ZHOOK_HANDLE)
ZTAUA2(1, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA2(2, :)= (/&
 & 0.40174_JPRB , 0.89383_JPRB , 0.89546_JPRB , 0.40741_JPRB , 0.51143_JPRB , 0.32646_JPRB /)  
 
ZPIZA2(1, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA2(2, :)= (/&
 & .8814597_JPRB, .9920407_JPRB, .9239428_JPRB, .7546879_JPRB, .9515548_JPRB, .9938563_JPRB/)  
 
ZCGA2(1, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA2(2, :)= (/&
 & 0.663224_JPRB, 0.793746_JPRB, 0.696315_JPRB, 0.652612_JPRB, .6608509_JPRB, .6318786_JPRB/)  

!        1.2   FOUR SPECTRAL INTERVALS (0.25-0.69-1.19-2.38-4.00microns)

ZTAUA4(1, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA4(2, :)= (/&
 & 0.52838_JPRB , 0.93285_JPRB , 0.93449_JPRB , 0.53078_JPRB , 0.67148_JPRB , 0.46608_JPRB /)  
ZTAUA4(3, :)= (/&
 & 0.20543_JPRB , 0.84642_JPRB , 0.84958_JPRB , 0.21673_JPRB , 0.28270_JPRB , 0.10915_JPRB /)  
ZTAUA4(4, :)= (/&
 & 0.10849_JPRB , 0.66699_JPRB , 0.65255_JPRB , 0.11600_JPRB , 0.06529_JPRB , 0.04468_JPRB /)  
 
ZPIZA4(1, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA4(2, :)= (/&
 & .8970131_JPRB, .9984940_JPRB, .9245594_JPRB, .7768385_JPRB, .9532763_JPRB, .9999999_JPRB/)  
ZPIZA4(3, :)= (/&
 & .8287144_JPRB, .9949396_JPRB, .9279543_JPRB, .6765051_JPRB, .9467578_JPRB, .9955938_JPRB/)  
ZPIZA4(4, :)= (/&
 & .5230504_JPRB, .7868518_JPRB, .8531531_JPRB, .4048149_JPRB, .8748231_JPRB, .2355667_JPRB/)  
 
ZCGA4(1, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA4(2, :)= (/&
 & 0.668431_JPRB, 0.788530_JPRB, 0.698682_JPRB, 0.657422_JPRB, .6735182_JPRB, .6519706_JPRB/)  
ZCGA4(3, :)= (/&
 & 0.636342_JPRB, 0.802467_JPRB, 0.691305_JPRB, 0.627497_JPRB, .6105750_JPRB, .4760794_JPRB/)  
ZCGA4(4, :)= (/&
 & 0.700610_JPRB, 0.818871_JPRB, 0.702399_JPRB, 0.689886_JPRB, .4629866_JPRB, .1907639_JPRB/)  

!        1.3   SIX SPECTRAL INTERVALS (0.185-0.25-0.44-0.69-1.19-2.38-4.00microns)

ZTAUA6(1, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA6(2, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA6(3, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA6(4, :)= (/&
 & 0.52838_JPRB , 0.93285_JPRB , 0.93449_JPRB , 0.53078_JPRB , 0.67148_JPRB , 0.46608_JPRB /)  
ZTAUA6(5, :)= (/&
 & 0.20543_JPRB , 0.84642_JPRB , 0.84958_JPRB , 0.21673_JPRB , 0.28270_JPRB , 0.10915_JPRB /)  
ZTAUA6(6, :)= (/&
 & 0.10849_JPRB , 0.66699_JPRB , 0.65255_JPRB , 0.11600_JPRB , 0.06529_JPRB , 0.04468_JPRB /)  
 
ZPIZA6(1, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA6(2, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA6(3, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA6(4, :)= (/&
 & .8970131_JPRB, .9984940_JPRB, .9245594_JPRB, .7768385_JPRB, .9532763_JPRB, .9999999_JPRB/)  
ZPIZA6(5, :)= (/&
 & .8287144_JPRB, .9949396_JPRB, .9279543_JPRB, .6765051_JPRB, .9467578_JPRB, .9955938_JPRB/)  
ZPIZA6(6, :)= (/&
 & .5230504_JPRB, .7868518_JPRB, .8531531_JPRB, .4048149_JPRB, .8748231_JPRB, .2355667_JPRB/)  
 
ZCGA6(1, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA6(2, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA6(3, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA6(4, :)= (/&
 & 0.668431_JPRB, 0.788530_JPRB, 0.698682_JPRB, 0.657422_JPRB, .6735182_JPRB, .6519706_JPRB/)  
ZCGA6(5, :)= (/&
 & 0.636342_JPRB, 0.802467_JPRB, 0.691305_JPRB, 0.627497_JPRB, .6105750_JPRB, .4760794_JPRB/)  
ZCGA6(6, :)= (/&
 & 0.700610_JPRB, 0.818871_JPRB, 0.702399_JPRB, 0.689886_JPRB, .4629866_JPRB, .1907639_JPRB/)  

!      ----------------------------------------------------------------

IF (KSW == 2) THEN
  DO JNU=1,KSW
    DO JAER=1,6
      YDESWRT%RTAUA(JNU,JAER)=ZTAUA2(JNU,JAER)
      YDESWRT%RPIZA(JNU,JAER)=ZPIZA2(JNU,JAER)
      YDESWRT%RCGA(JNU,JAER) =ZCGA2 (JNU,JAER)
    ENDDO
  ENDDO
ELSEIF (KSW == 4) THEN
  DO JNU=1,KSW
    DO JAER=1,6
      YDESWRT%RTAUA(JNU,JAER)=ZTAUA4(JNU,JAER)
      YDESWRT%RPIZA(JNU,JAER)=ZPIZA4(JNU,JAER)
      YDESWRT%RCGA(JNU,JAER) =ZCGA4 (JNU,JAER)
    ENDDO
  ENDDO
ELSEIF (KSW == 6) THEN
  DO JNU=1,KSW
    DO JAER=1,6
      YDESWRT%RTAUA(JNU,JAER)=ZTAUA6(JNU,JAER)
      YDESWRT%RPIZA(JNU,JAER)=ZPIZA6(JNU,JAER)
      YDESWRT%RCGA(JNU,JAER) =ZCGA6 (JNU,JAER)
    ENDDO
  ENDDO
ELSEIF (KSW == 14 .AND. KTSW == 14) THEN
  WRITE (NULOUT,*) 'SUAERSN: 14-SPECTRAL INTERVALS --> RRTM_SW'
ELSE
  STOP 'SUAERSN: WRONG NUMBER OF SPECTRAL INTERVALS'
ENDIF

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUAERSN',1,ZHOOK_HANDLE)
END SUBROUTINE SUAERSN
