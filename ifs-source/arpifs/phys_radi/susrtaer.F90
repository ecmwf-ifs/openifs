! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUSRTAER 

!**** *SUSRTAER*   - INITIALIZE COMMON YOESRTAER

!     PURPOSE.
!     --------
!           INITIALIZE YOESRTAER, THE COMMON THAT CONTAINS THE
!           RADIATIVE CHARACTERISTICS OF THE AEROSOLS

!**   INTERFACE.
!     ----------
!              -----        -----

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOESRTAER

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
!        JEAN-JACQUES MORCRETTE *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-03-06
!        03-03-06  JJMorcrette   14 spectral intervals for RRTM_SW
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOESRTAER   , ONLY : RSRTAUA     ,RSRPIZA    ,RSRASYA

!      ----------------------------------------------------------------

IMPLICIT NONE

!     DUMMY INTEGER REAL
REAL(KIND=JPRB) :: ZTAUA14(14,6),ZPIZA14(14,6),ZCGA14(14,6)

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

!       1.1    FOURTEEN SPECTRAL INTERVALS FOR RRTM_SW

! SW : 14 spectral intervals 
!  3.846 -  3.077
!  3.077 -  2.500
!  2.500 -  2.150
!  2.150 -  1.942
!  1.942 -  1.626
!  1.626 -  1.299
!  1.299 -  1.242
!  1.242 -  0.7782
!  0.7782-  0.6250
!  0.6250-  0.4415
!  0.4415-  0.3448
!  0.3448-  0.2632
!  0.2632-  0.2000
! 12.195 -  3.846

IF (LHOOK) CALL DR_HOOK('SUSRTAER',0,ZHOOK_HANDLE)
ZTAUA14( 1, :)= (/&
 & 0.10849_JPRB , 0.66699_JPRB , 0.65255_JPRB , 0.11600_JPRB , 0.06529_JPRB , 0.04468_JPRB /)  
ZTAUA14( 2, :)= (/&
 & 0.10849_JPRB , 0.66699_JPRB , 0.65255_JPRB , 0.11600_JPRB , 0.06529_JPRB , 0.04468_JPRB /)  
ZTAUA14( 3, :)= (/&
 & 0.20543_JPRB , 0.84642_JPRB , 0.84958_JPRB , 0.21673_JPRB , 0.28270_JPRB , 0.10915_JPRB /)  
ZTAUA14( 4, :)= (/&
 & 0.20543_JPRB , 0.84642_JPRB , 0.84958_JPRB , 0.21673_JPRB , 0.28270_JPRB , 0.10915_JPRB /)  
ZTAUA14( 5, :)= (/&
 & 0.20543_JPRB , 0.84642_JPRB , 0.84958_JPRB , 0.21673_JPRB , 0.28270_JPRB , 0.10915_JPRB /)  
ZTAUA14( 6, :)= (/&
 & 0.20543_JPRB , 0.84642_JPRB , 0.84958_JPRB , 0.21673_JPRB , 0.28270_JPRB , 0.10915_JPRB /)  
ZTAUA14( 7, :)= (/&
 & 0.20543_JPRB , 0.84642_JPRB , 0.84958_JPRB , 0.21673_JPRB , 0.28270_JPRB , 0.10915_JPRB /)  
ZTAUA14( 8, :)= (/&
 & 0.52838_JPRB , 0.93285_JPRB , 0.93449_JPRB , 0.53078_JPRB , 0.67148_JPRB , 0.46608_JPRB /)  
ZTAUA14( 9, :)= (/&
 & 0.52838_JPRB , 0.93285_JPRB , 0.93449_JPRB , 0.53078_JPRB , 0.67148_JPRB , 0.46608_JPRB /)  
ZTAUA14(10, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA14(11, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA14(12, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA14(13, :)= (/&
 & 1.69446_JPRB , 1.11855_JPRB , 1.09212_JPRB , 1.72145_JPRB , 1.03858_JPRB , 1.12044_JPRB /)  
ZTAUA14(14, :)= (/&
 & 0.10849_JPRB , 0.66699_JPRB , 0.65255_JPRB , 0.11600_JPRB , 0.06529_JPRB , 0.04468_JPRB /)  
 
ZPIZA14( 1, :)= (/&
 & .5230504_JPRB, .7868518_JPRB, .8531531_JPRB, .4048149_JPRB, .8748231_JPRB, .2355667_JPRB/)  
ZPIZA14( 2, :)= (/&
 & .5230504_JPRB, .7868518_JPRB, .8531531_JPRB, .4048149_JPRB, .8748231_JPRB, .2355667_JPRB/)  
ZPIZA14( 3, :)= (/&
 & .8287144_JPRB, .9949396_JPRB, .9279543_JPRB, .6765051_JPRB, .9467578_JPRB, .9955938_JPRB/)  
ZPIZA14( 4, :)= (/&
 & .8287144_JPRB, .9949396_JPRB, .9279543_JPRB, .6765051_JPRB, .9467578_JPRB, .9955938_JPRB/)  
ZPIZA14( 5, :)= (/&
 & .8287144_JPRB, .9949396_JPRB, .9279543_JPRB, .6765051_JPRB, .9467578_JPRB, .9955938_JPRB/)  
ZPIZA14( 6, :)= (/&
 & .8287144_JPRB, .9949396_JPRB, .9279543_JPRB, .6765051_JPRB, .9467578_JPRB, .9955938_JPRB/)  
ZPIZA14( 7, :)= (/&
 & .8287144_JPRB, .9949396_JPRB, .9279543_JPRB, .6765051_JPRB, .9467578_JPRB, .9955938_JPRB/)  
ZPIZA14( 8, :)= (/&
 & .8970131_JPRB, .9984940_JPRB, .9245594_JPRB, .7768385_JPRB, .9532763_JPRB, .9999999_JPRB/)  
ZPIZA14( 9, :)= (/&
 & .8970131_JPRB, .9984940_JPRB, .9245594_JPRB, .7768385_JPRB, .9532763_JPRB, .9999999_JPRB/)  
ZPIZA14(10, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA14(11, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA14(12, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA14(13, :)= (/&
 & .9148907_JPRB, .9956173_JPRB, .7504584_JPRB, .8131335_JPRB, .9401905_JPRB, .9999999_JPRB/)  
ZPIZA14(14, :)= (/&
 & .5230504_JPRB, .7868518_JPRB, .8531531_JPRB, .4048149_JPRB, .8748231_JPRB, .2355667_JPRB/)  

ZCGA14( 1, :)= (/&
 & 0.700610_JPRB, 0.818871_JPRB, 0.702399_JPRB, 0.689886_JPRB, .4629866_JPRB, .1907639_JPRB/)  
ZCGA14( 2, :)= (/&
 & 0.700610_JPRB, 0.818871_JPRB, 0.702399_JPRB, 0.689886_JPRB, .4629866_JPRB, .1907639_JPRB/)  
ZCGA14( 3, :)= (/&
 & 0.636342_JPRB, 0.802467_JPRB, 0.691305_JPRB, 0.627497_JPRB, .6105750_JPRB, .4760794_JPRB/)  
ZCGA14( 4, :)= (/&
 & 0.636342_JPRB, 0.802467_JPRB, 0.691305_JPRB, 0.627497_JPRB, .6105750_JPRB, .4760794_JPRB/)  
ZCGA14( 5, :)= (/&
 & 0.636342_JPRB, 0.802467_JPRB, 0.691305_JPRB, 0.627497_JPRB, .6105750_JPRB, .4760794_JPRB/)  
ZCGA14( 6, :)= (/&
 & 0.636342_JPRB, 0.802467_JPRB, 0.691305_JPRB, 0.627497_JPRB, .6105750_JPRB, .4760794_JPRB/)  
ZCGA14( 7, :)= (/&
 & 0.636342_JPRB, 0.802467_JPRB, 0.691305_JPRB, 0.627497_JPRB, .6105750_JPRB, .4760794_JPRB/)  
ZCGA14( 8, :)= (/&
 & 0.668431_JPRB, 0.788530_JPRB, 0.698682_JPRB, 0.657422_JPRB, .6735182_JPRB, .6519706_JPRB/)  
ZCGA14( 9, :)= (/&
 & 0.668431_JPRB, 0.788530_JPRB, 0.698682_JPRB, 0.657422_JPRB, .6735182_JPRB, .6519706_JPRB/)  
ZCGA14(10, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA14(11, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA14(12, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA14(13, :)= (/&
 & 0.729019_JPRB, 0.803129_JPRB, 0.784592_JPRB, 0.712208_JPRB, .7008249_JPRB, .7270548_JPRB/)  
ZCGA14(14, :)= (/&
 & 0.700610_JPRB, 0.818871_JPRB, 0.702399_JPRB, 0.689886_JPRB, .4629866_JPRB, .1907639_JPRB/)  

!      ----------------------------------------------------------------

DO JNU=1,14
  DO JAER=1,6
    RSRTAUA(JNU,JAER)=ZTAUA14(JNU,JAER)
    RSRPIZA(JNU,JAER)=ZPIZA14(JNU,JAER)
    RSRASYA(JNU,JAER)=ZCGA14 (JNU,JAER)
  ENDDO
ENDDO

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSRTAER',1,ZHOOK_HANDLE)
END SUBROUTINE SUSRTAER

