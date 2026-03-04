! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

FUNCTION FCGENERALIZED_GAMMA(PX)  

!****  *FCGENERALIZED_GAMMA * -  Gamma  function


!    PURPOSE
!    -------
!       The purpose of this function is to compute the Generalized gamma
!    function of its argument.


!**  METHOD
!    ------

!    EXTERNAL
!    --------
!      NONE

!    IMPLICIT ARGUMENTS
!    ------------------
!      None

!    REFERENCE
!    ---------
!      Press, Teukolsky, Vetterling and Flannery: Numerical Recipes, 206-207

!    AUTHOR
!    ------
!      Jean-Pierre Pinty *LA/OMP*

!    MODIFICATIONS
!    -------------
!      Original     7/11/95
!      S.Martinez   2008/03/11 - rename function GAMMA to GENERALIZED_GAMMA
!                                to avoid confusion with not always intrinsic
!                                GAMMA function...
!      E. Wattrelot 2009/01/14 - transform module yomgamma in function generalized_gamma
!      E. Wattrelot 2010/04/22 - vectorization of gamma function

!*       0. DECLARATIONS
!           ------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST  ,ONLY : RPI

IMPLICIT NONE

!*       0.1 declarations of arguments and result

REAL(KIND=JPRB), INTENT(IN)                     :: PX

!*       0.2 declarations of local variables

!INTEGER(KIND=JPIM)      :: JJ ! Loop index
!REAL(KIND=JPRB)         :: ZSER,ZSTP,ZTMP,ZX,ZY,ZCOEF(6)
!REAL(KIND=JPRB)         :: ZPI

REAL(KIND=JPRB)         :: ZSTP=2.5066282746310005
REAL(KIND=JPRB), DIMENSION(6)   :: ZCOEF=(/ 76.18009172947146,-86.50532032941677,24.01409824083091, &
       &          -1.231739572450155,0.1208650973866179E-2,-0.5395239384953E-5 /)
REAL(KIND=JPRB)         :: ZSERCST=1.000000000190015

REAL(KIND=JPRB)         :: ZGAMMA,  FCGENERALIZED_GAMMA
REAL(KIND=JPRB)         :: ZSER,ZTTMMP,ZTMP,ZX,ZY

REAL(KIND=JPHOOK)         :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('FCGENERALIZED_GAMMA',0,ZHOOK_HANDLE)

!ZCOEF(1) = 76.18009172947146
!ZCOEF(2) =-86.50532032941677
!ZCOEF(3) = 24.01409824083091
!ZCOEF(4) = -1.231739572450155
!ZCOEF(5) =  0.1208650973866179E-2
!ZCOEF(6) = -0.5395239384953E-5
!ZSTP     =  2.5066282746310005

!ZPI = 3.141592654
!IF (PX < 0.) THEN
!  ZX = 1.- PX
!ELSE 
!  ZX = PX
!ENDIF
!ZY = ZX
!ZTMP =  ZX + 5.5
!ZTMP = (ZX + 0.5)*LOG(ZTMP) - ZTMP
!ZSER = 1.000000000190015

!DO JJ = 1 , 6
!  ZY = ZY + 1.0
!  ZSER = ZSER + ZCOEF(JJ)/ZY
!ENDDO

!IF (PX < 0.) THEN
!  ZGAMMA = RPI/SIN(RPI*PX)/EXP( ZTMP + LOG( ZSTP*ZSER/ZX ) )
!ELSE
!  ZGAMMA = EXP( ZTMP + LOG( ZSTP*ZSER/ZX ) )
!ENDIF

   ZX = (1._JPRB-PX)*(1._JPRB-MAX(SIGN(1._JPRB,PX),0._JPRB)) & 
  &   + PX*MAX(SIGN(1._JPRB,PX),0._JPRB)   
ZY = ZX
ZTTMMP =  ZX + 5.5
ZTMP = (ZX + 0.5)*LOG(ZTTMMP) - ZTTMMP


    ZSER = ZSERCST + ZCOEF(1)/(1.+ZY) + ZCOEF(2)/(2.+ZY) &
   &     + ZCOEF(3)/(3.+ZY) + ZCOEF(4)/(4.+ZY) &
   &     + ZCOEF(5)/(5.+ZY) + ZCOEF(6)/(6.+ZY)

   ZGAMMA = EXP( ZTMP + LOG( ZSTP*ZSER/ZX ) )
   FCGENERALIZED_GAMMA = ZGAMMA * MAX(SIGN(1._JPRB,PX),0._JPRB) &
    &                    + (RPI/SIN(RPI*PX)/ZGAMMA)         &
    &                    *(1._JPRB-MAX(SIGN(1._JPRB,PX),0._JPRB))


IF (LHOOK) CALL DR_HOOK('FCGENERALIZED_GAMMA',1,ZHOOK_HANDLE)

END FUNCTION FCGENERALIZED_GAMMA
