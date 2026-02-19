! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_SOA(KNSOA,PT,PSOG,PSOA,PORGAEROP)


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!   Solution of SOG - SOA equilibrium
!
!
!**   INTERFACE.
!     ----------
!          *TM5_SOA* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KNSOA                  : Number of SOA tracers
! PT                     : Temperature [K]
! PSOG                   : initial SOA precursor gas (product of VOC oxidation) (ug/m3)
! PSOA                   : initial secondary aerosol concentration              (ug/m3)
! PORGAEROP              : CAMS Organic aerosol concentrations                  (ug/m3)
! 
!
!
! OUTPUTS:
! -------
! PSOG (KNSOA)               : final, re-distributed SOA precursor gas                ( ug/m3)
! PSOA (KNSOA)               : final, re-distributed secondary aerosol concentration  ( ug/m3)
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        Vincent Huijnen *KNMI*
!
!     MODIFICATIONS.
!     --------------
!        First implementation : 2017-04-21
!      
!      references
!      ---------
!      
!_________________________________________________________________________________________





USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! USE YOMLUN   , ONLY : NULOUT, NULERR 

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------


INTEGER(KIND=JPIM), INTENT(IN)    :: KNSOA
REAL(KIND=JPRB)   , INTENT(IN)    :: PT
REAL(KIND=JPRB)   , INTENT(INOUT) :: PSOG(KNSOA),PSOA(KNSOA)
REAL(KIND=JPRB)   , INTENT(IN)    :: PORGAEROP

! * LOCAL 

!* Currently two volatility classes are implemented, matching
!* with the isoprene and terpene oxidation schemes.
!* First  class has K coefficient of ~ 0.1 - 1.0 (terpenes)
!* Second class has K coefficient of ~ 1.0 - 10  (isop+OH)
!* For now, with limited reactions, just take the most representative  value of K.

!* Equilibrium gas-particle partition coefficients of 
!* semi-volatile compounts [ug-1 m**3]
REAL(KIND=JPRB), DIMENSION(2),PARAMETER :: ZKOM_REF = (/ &
  & 1.0_JPRB , &      !  OH / O3 + Terpene reaction product; OH + ISOP. Net K of unity might be on high side.
!VH  & 0.01_JPRB, &      !  Anthro SOA, Low K (high C*) -> No longer explicitly included.
  & 1.0_JPRB  /)      !  Aged Anthro SOA, high K (low C*)
REAL(KIND=JPRB), DIMENSION(2) :: ZKOM
REAL(KIND=JPRB), DIMENSION(2) :: ZORG_AER, ZORG_GAS
REAL(KIND=JPRB), DIMENSION(2) :: ZSOA_IN, ZSOG_IN

! Heat of vaporization over R (Takekawa, AE 2003)
! is also dependant on SOA type..
REAL(KIND=JPRB), DIMENSION(2), PARAMETER    :: ZHEAT_VAPORR = (/ &
  & 1.0E4_JPRB, 5.E3_JPRB /)   

REAL(KIND=JPRB), PARAMETER    :: ZRT310=1./310_JPRB
REAL(KIND=JPRB), PARAMETER    :: ZRT295=1./295_JPRB
REAL(KIND=JPRB), PARAMETER    :: ZRELAX = 0.2 ! Relaxation parameter to obtain more smooth solution

REAL(KIND=JPRB)               :: ZRPT, ZT_310,ZT_295, ZTMP2,ZMT_ORGAER,ZFRAC

! Number of iterations to assumed needed to reach convergence 
INTEGER(KIND=JPIM)            :: ITER,JK
INTEGER(KIND=JPIM), PARAMETER :: JNITER = 3

REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('TM5_SOA',0,ZHOOK_HANDLE )

IF ( KNSOA > 3 ) THEN
    CALL ABOR1('TM5_SOA: too many SOA sources specified')
ENDIF

!* Initialization of ZKOM:

!Reciprocal temperature [1/K]
ZRPT = 1.0_JPRB / PT

! Divide TEMP by 310K 
ZT_310 =  PT / 310_JPRB

! Compute the heat-of-vaporization exponential term 
ZTMP2 = EXP( ZHEAT_VAPORR(1) * ( ZRPT - ZRT310 ) )

! first for terpenes
ZKOM(1) = ZKOM_REF(1) * ZT_310 * ZTMP2
 
ZT_295 = PT / 295_JPRB
! Compute the heat-of-vaporization exponential term outside the DO loop
ZTMP2 = EXP( ZHEAT_VAPORR(2) * ( ZRPT - ZRT295 ) )

! next for VOC's - No longer included 
!VH ZKOM(2) = ZKOM_REF(2) * ZT_295 * ZTMP2

! next for VOC's
ZKOM(2) = ZKOM_REF(2) * ZT_295 * ZTMP2

! Initial SOA
DO JK=1,KNSOA
  ZSOA_IN(JK)=PSOA(JK)
  ZSOG_IN(JK)=PSOG(JK)
ENDDO



!-----------------------------------------------------------
! Compute SOG condensation onto OC aerosol

! do several iterations
DO ITER =1,JNITER

  ! First estimate of total organic aerosol available for condensation 
  ! of additional SOA
  ZMT_ORGAER=PORGAEROP
  DO JK=1, KNSOA
    ZMT_ORGAER  = ZMT_ORGAER+PSOA(JK)
  ENDDO


  DO JK=1,KNSOA
    ZFRAC=ZKOM(JK)*ZMT_ORGAER / &
      &  (1._JPRB + ZKOM(JK) * ZMT_ORGAER ) 
    ZFRAC=MAX(0.0_JPRB,MIN(ZFRAC,1.0_JPRB))
    

    ZORG_AER(JK) = ZFRAC * (PSOG(JK) + PSOA(JK))
    ! Apply some relaxation towards original concentration,
    ZORG_AER(JK) = ZRELAX*PSOA(JK) + (1.-ZRELAX)*ZORG_AER(JK)

    ! Limit conversion of aerosol back to gas conversion to 30%  at most for each time step.
    ZORG_AER(JK) =MAX(ZSOA_IN(JK)*0.7,ZORG_AER(JK))

    ! All that is not in aerosol phase is in gas-phase:
    ZORG_GAS(JK) = ( PSOG(JK)+PSOA(JK) ) - ZORG_AER(JK)  
  ENDDO
  
  ! update fields 
  DO JK=1,KNSOA
    PSOG(JK)=ZORG_GAS(JK)
    PSOA(JK)=ZORG_AER(JK)
  ENDDO
ENDDO

! Check on concentrations
DO JK=1,KNSOA
  ! Limit SOA to SOA+SOG
  PSOA(JK)=MIN(PSOA(JK),ZSOA_IN(JK)+ZSOG_IN(JK))
  ! Compute SOG to be consistent with SOA_in+SOG_in - SOA_out
  PSOG(JK)=ZSOA_IN(JK)+ZSOG_IN(JK) - PSOA(JK)
ENDDO


IF (LHOOK) CALL DR_HOOK('TM5_SOA',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_SOA

