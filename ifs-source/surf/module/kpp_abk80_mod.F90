MODULE KPP_ABK80_MOD
CONTAINS
SUBROUTINE KPP_ABK80 &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEV     ,LDKPPCAL ,&
  &   PSAL     ,PTEMP    ,PRES     ,PALPHA   ,PBETA    ,&
  &   PSIG0    ,PSIG     ,PRHO0    ,LDALPHA  ,LDBETA   ) 

! (C) Copyright 1987- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
! Subroutine to coordinate computation of the expansion coefficients  
! of seawater with respect to temperature (alpha), salinity (beta) 
! using the 1980 equation of state for seawater. 

! Interface :
! ---------
!   Call *ABK80* with PSAL(practical salinity)
!                     PTEMP(temperature (degC))
!                     PRES(pressure (dbar))

! Method :
! ------
! Analytic calculation of expansion coefficients with the polynomial 
! structure of the 1980 equation of state. 

! Externals :
! ---------

! Reference :
! ---------
!             UNESCO(1981): Background papers and supporting data on
!                           the international equation of state of
!                           seawater, 1980. UNESCO Tech. Pap. in Mar.
!                           Sci., No. 38, 192 pp.
!
!             FOFONOFF, N. P., and R. C. Millard, Jr. (1983): 
!                           Algorithms for computation of fundamental 
!                           properties of seawater.  UNESCO Tech. Pap. 
!                           in Mar. Sci., No. 44, 53 pp.
!
!             Lillibridge, J.L. (1988): Computing the seawater expansion
!                           coefficients directly from the 1980 equation
!                           of state.  J. Atm. Ocean. Tech., 59-66.

! Modifications :
! -------------
!     25-Mar-1987  John L. Lillibridge, URI/GSO   Original codes
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

IMPLICIT NONE 

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV

REAL(KIND=JPRB),INTENT(IN) :: PSAL(KIDIA:KFDIA,KLEV)  ! practical salinity
REAL(KIND=JPRB),INTENT(IN) :: PTEMP(KIDIA:KFDIA,KLEV) ! temperature (deg. C)
REAL(KIND=JPRB),INTENT(IN) :: PRES(KIDIA:KFDIA,KLEV)  ! pressure (dbar)
REAL(KIND=JPRB),INTENT(INOUT) :: PALPHA(KIDIA:KFDIA,KLEV) 
                                                   ! thermal expansion coef. 
                                                   ! alpha (degC^-1)
REAL(KIND=JPRB),INTENT(INOUT) :: PBETA(KIDIA:KFDIA,KLEV)
                                                   ! haline contraction coef.
                                                   ! beta (nondim.)
                                                   ! {divided by 10**-3 
                                                   ! for 1978 practical sal.}
REAL(KIND=JPRB),INTENT(INOUT) :: PSIG(KLON,KLEV)   ! sigma  (sigma-t kg/m**3) 
REAL(KIND=JPRB),INTENT(INOUT) :: PSIG0(KLON,KLEV)  ! sigma0 
REAL(KIND=JPRB),INTENT(OUT)   :: PRHO0(KLON,KLEV)  ! ...

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)       
LOGICAL,INTENT(IN) :: LDALPHA
LOGICAL,INTENT(IN) :: LDBETA

INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: JZ

REAL(KIND=JPRB) :: ZTEMP          ! temperature (>-2)
REAL(KIND=JPRB) :: ZR1            ! Polynomial coefficient
REAL(KIND=JPRB) :: ZR2            ! ...
REAL(KIND=JPRB) :: ZR3            ! ...
REAL(KIND=JPRB) :: ZR4            ! ...
REAL(KIND=JPRB) :: ZA             ! ...
REAL(KIND=JPRB) :: ZB             ! ...
REAL(KIND=JPRB) :: ZC             ! ...
REAL(KIND=JPRB) :: ZD             ! ... 
REAL(KIND=JPRB) :: ZE             ! ...
REAL(KIND=JPRB) :: ZA1            ! ...
REAL(KIND=JPRB) :: ZB1            ! ...   
REAL(KIND=JPRB) :: ZP0            ! ...
REAL(KIND=JPRB) :: ZSR            ! ...
REAL(KIND=JPRB) :: ZP0K           ! ...
REAL(KIND=JPRB) :: ZK             ! ...
REAL(KIND=JPRB) :: ZRHO           ! ...
REAL(KIND=JPRB) :: ZABFAC         ! ...   
REAL(KIND=JPRB) :: ZDA
REAL(KIND=JPRB) :: ZDB
REAL(KIND=JPRB) :: ZDK
REAL(KIND=JPRB) :: ZDK0
REAL(KIND=JPRB) :: ZDRHO
REAL(KIND=JPRB) :: ZSR5
REAL(KIND=JPRB) :: ZK0
REAL(KIND=JPRB) :: ZAW
REAL(KIND=JPRB) :: ZBW
REAL(KIND=JPRB) :: ZKW
REAL(KIND=JPRB) :: ZALPH0
REAL(KIND=JPRB) :: ZALPHAA
REAL(KIND=JPRB) :: ZALPHB
REAL(KIND=JPRB) :: ZALPHK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!LOGICAL         :: LLKAPFLG   
LOGICAL         :: LLABFLG
LOGICAL         :: LLPRES

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_ABK80_MOD:KPP_ABK80',0,ZHOOK_HANDLE)

!LLKAPFLG = .FALSE. 
ZR4 = 4.8314E-4_JPRB
ZD  = 1.91075E-4_JPRB

! 1. Compute the haline expansion coefficient
!    ----------------------------------------

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
    DO JZ = 1, KLEV
      ZTEMP = MAX( PTEMP(JL,JZ), -2.0_JPRB )
      LLABFLG  = .TRUE.

! subroutine KPP_SIG80 is inlined here

      ZP0 = PRES(JL,JZ) * 0.1_JPRB
      ZSR = SQRT( ABS(PSAL(JL,JZ)) )

      ! compute sigma pure water at atm press
      ZR1 = ( ( ( ( 6.536332E-9_JPRB * ZTEMP      &
          &         - 1.120083E-6_JPRB ) * ZTEMP  &
          &       + 1.001685E-4_JPRB ) * ZTEMP    &
          &     - 9.095290E-3_JPRB ) * ZTEMP      &
          &   + 6.793952E-2_JPRB ) * ZTEMP - 0.157406_JPRB

      ! Note constant for SIGMA vs. RHO

      ! seawater sigma at atm press
      ZR2 = ( ( ( 5.3875E-9_JPRB * ZTEMP          &
          &       - 8.2467E-7_JPRB ) * ZTEMP      &
          &     + 7.6438E-5_JPRB ) * ZTEMP        &
          &   - 4.0899E-3_JPRB ) * ZTEMP + 8.24493E-1_JPRB
      ZR3 = ( - 1.6546E-6_JPRB * ZTEMP            &
          &   + 1.0227E-4_JPRB ) * ZTEMP - 5.72466E-3_JPRB
!      ZR4 = 4.8314E-4_JPRB
      PSIG0(JL,JZ) = ( ZR4 * PSAL(JL,JZ) + ZR3 * ZSR &
                   &  + ZR2 ) * PSAL(JL,JZ) + ZR1
      PRHO0(JL,JZ) = 1000.0_JPRB + PSIG0(JL,JZ)

!     sigma at atm press
      IF( PRES(JL,JZ) == 0.0_JPRB )THEN
        PSIG(JL,JZ) = PSIG0(JL,JZ)
        ZRHO = PRHO0(JL,JZ)
      ENDIF

! compute compression terms for computation of kappa

! this is the entry point when the bulk modulus k is desired, which
! doesn't require knowledge of density. entry from kappa will have
! kapflg set .true. so that we exit before computing sig.

      ZA1 = ( ( -6.1670E-5_JPRB * ZTEMP        &
          &     + 1.09987E-2_JPRB ) * ZTEMP    &
          &   - 0.603459_JPRB ) * ZTEMP + 54.6746_JPRB
      ZB1 = ( - 5.3009E-4_JPRB * ZTEMP         &
          &   + 1.6483E-2_JPRB ) * ZTEMP       &
          & + 7.944E-2_JPRB
      ZKW = ( ( ( - 5.155288E-5_JPRB * ZTEMP   &
          &       + 1.360477E-2_JPRB ) * ZTEMP &
          &     - 2.327105_JPRB ) * ZTEMP      &
          &   + 148.4206_JPRB ) * ZTEMP        &
          & +19652.21_JPRB
      ZK0 = (ZB1*ZSR+ZA1)*PSAL(JL,JZ)+ZKW

      IF( PRES(JL,JZ) == 0.0_JPRB )THEN
        ZK = ZK0
      ELSE

!       pressure terms in bulk modulus
        ZE = ( 9.1697E-10_JPRB * ZTEMP + 2.0816E-8_JPRB )      &
           &  * ZTEMP -9.9348E-7_JPRB
        ZBW = ( 5.2787E-8_JPRB * ZTEMP - 6.12293E-6_JPRB )     &
            & * ZTEMP + 8.50935E-5_JPRB
        ZB  = ZBW + ZE * PSAL(JL,JZ)
        ZC  = ( -1.6078E-6_JPRB * ZTEMP - 1.0981E-5_JPRB )     &
            & * ZTEMP + 2.2838E-3_JPRB
!        ZD  = 1.91075E-4_JPRB
        ZAW = ( ( -5.77905E-7_JPRB * ZTEMP + 1.16092E-4_JPRB ) &
            &   * ZTEMP + 1.43713E-3_JPRB ) * ZTEMP            &
            & + 3.239908_JPRB
        ZA  = ( ZD * ZSR + ZC ) * PSAL(JL,JZ) + ZAW

!       evaluate press polynomial and return
        ZK = ( ZB * ZP0 + ZA ) * ZP0 + ZK0

        ZP0K = ZP0 / ZK
!        IF ( .NOT. LLKAPFLG ) THEN
          PSIG(JL,JZ) = ( 1000.0_JPRB * ZP0K + PSIG0(JL,JZ) ) &
                      & / ( 1.0_JPRB - ZP0K )
          ZRHO = 1000.0_JPRB + PSIG(JL,JZ)
!        ENDIF

      ENDIF

! KPP_SIG80 end

! subroutine KPP_BET80 inlined here

      ZSR5  = ZSR * 1.5_JPRB
      ZDRHO = ZR2 + ZSR5 * ZR3 &
                & + ( PSAL(JL,JZ) + PSAL(JL,JZ) ) * ZR4
      IF ( PRES(JL,JZ) /= 0.0_JPRB ) THEN
!       next compute the derivative terms for the bulk modulus
        ZDK0 = ZA1 + ZSR5 * ZB1
        ZDA  = ZC + ZSR5 * ZD
        ZDB  = ZE

!       derivative dk/ds
        ZDK = ( ZDB * ZP0 + ZDA ) * ZP0 + ZDK0

!       assemble beta from all the terms (rho0,k,pk from sig80)
        ZABFAC = PRHO0(JL,JZ) * ZP0  / ( ( ZK - ZP0 ) * ( ZK - ZP0 ) )
        PBETA(JL,JZ) = ZDRHO / ( 1.0_JPRB - ZP0K ) - ZABFAC * ZDK
        PBETA(JL,JZ) = PBETA(JL,JZ) / ZRHO
        LLABFLG = .FALSE.
      ELSE ! if PRES == 0.0
        PBETA(JL,JZ) = ZDRHO / ZRHO
      ENDIF

! 2. Compute the thermal expansion coefficient 
!    ----------------------------------------

! subroutine KPP_ALF80 inlined here

!     compute alpha pure water at atm press
      ZR1 = ( ( ( 0.3268166E-7_JPRB * ZTEMP                 &
          &       - 0.4480332E-5_JPRB ) * ZTEMP             &
          &     + 0.3005055E-3_JPRB ) * ZTEMP               &
          &   - 0.1819058E-1_JPRB ) * ZTEMP + 6.793952E-2_JPRB

!     seawater alpha at atm press
      ZR2 = ( ( 0.215500E-7_JPRB * ZTEMP             &
                 &     - 0.247401E-5_JPRB ) * ZTEMP         &
                 &   + 0.152876E-3_JPRB ) * ZTEMP - 4.0899E-3_JPRB

      ZR3 = -0.33092E-5_JPRB * ZTEMP + 1.0227E-4_JPRB
      ZALPH0 = ( ZR3 * ZSR + ZR2 ) * PSAL(JL,JZ)  + ZR1

!     alpha at atm press
      IF( PRES(JL,JZ) /= 0.0_JPRB ) THEN

!       compute compression terms

        ZB1 = -0.106018E-2_JPRB * ZTEMP + 1.6483E-2_JPRB
        ZA1 = ( -0.18501E-3_JPRB * ZTEMP + 0.219974E-1_JPRB ) &
            & * ZTEMP - 0.603459_JPRB
        ZKW = ( ( - 0.2062115E-3_JPRB * ZTEMP                 &
            &     + 0.4081431E-1_JPRB ) * ZTEMP               &
            &   - 0.4654210E+1_JPRB ) * ZTEMP + 148.4206_JPRB
        ZK0 = ( ZB1 * ZSR + ZA1 ) * PSAL(JL,JZ) + ZKW

!         pressure terms in bulk modulus

        ZE  = 0.183394E-8_JPRB * ZTEMP + 2.0816E-8_JPRB
        ZBW = 0.105574E-6_JPRB * ZTEMP - 6.12293E-6_JPRB
        ZALPHB = ZBW + ZE * PSAL(JL,JZ)

        ZC  = -0.32156E-5_JPRB * ZTEMP - 1.0981E-5
        ZAW = ( -0.1733715E-5_JPRB * ZTEMP + 0.232184E-3_JPRB ) &
            & * ZTEMP + 1.43713E-3_JPRB
        ZALPHAA = ZC * PSAL(JL,JZ) + ZAW

!       evaluate press polynomial and return

        ZALPHK = ( ZALPHB*ZP0 + ZALPHAA ) * ZP0 + ZK0

        IF ( LLABFLG ) THEN
          ZABFAC = PRHO0(JL,JZ) * ZP0  / ( (ZK-ZP0) * (ZK-ZP0) )
        ENDIF
        PALPHA(JL,JZ) = ZALPH0 / (1.0_JPRB - ZP0K ) &
                      & - ZABFAC * ZALPHK
        PALPHA(JL,JZ) = - PALPHA(JL,JZ) / ZRHO

      ELSE ! if PRES == 0.0_JPRB

        PALPHA(JL,JZ) = -ZALPH0 / ZRHO

      ENDIF

    ENDDO !JZ
  ENDIF
ENDDO !JL

IF (LHOOK) CALL DR_HOOK('KPP_ABK80_MOD:KPP_ABK80',1,ZHOOK_HANDLE)
 
END SUBROUTINE KPP_ABK80 
END MODULE KPP_ABK80_MOD
