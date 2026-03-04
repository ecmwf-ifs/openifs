! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_TROPOPAUSE(KIDIA,KFDIA,KLON,KLEV,PTP,PRSF1,PGEOH,PLAT,KTROPOP)

!**   DESCRIPTION 
!     ----------
! Calculate KTROPOP, the pressure index of tropopause
!                                                simonc, v2s35, Jul 2003
!
!       05/03/2014 - taken over to IFS, V. Huijnen
!-----------------------------------------------------------------------
USE PARKIND1      , ONLY : JPIM ,   JPRB
USE YOMHOOK       , ONLY : LHOOK,   DR_HOOK, JPHOOK
! USE YOMLUN        , ONLY : NULOUT
USE BASCOE_MODULE , ONLY : PTROP_SOC_B,NLAT_PTROPO
USE YOMCST        , ONLY : RG
IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA , KFDIA , KLON , KLEV
REAL(KIND=JPRB)   , INTENT(IN)  :: PTP(KLON,KLEV),PRSF1(KLON,KLEV),PGEOH(KLON,KLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PLAT(KLON)
INTEGER(KIND=JPIM), INTENT(OUT) :: KTROPOP(KLON)

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM)              :: JL,JK,JK2,JLAT_SOC,JTROPOP_SOC,JK_TROP
REAL(KIND=JPRB)                 :: ZTHKNESS,ZGRADT
REAL(KIND=JPRB)                 :: ZRGI,ZKM
REAL(KIND=JPHOOK)                 :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('BASCOE_TROPOPAUSE',0,ZHOOK_HANDLE )

ZRGI=1.0_JPRB/RG


   DO JL=KIDIA,KFDIA

    KTROPOP(JL) = -9
    DO JK = KLEV-1, 1, -1
!-----------------------------------------------------------------------
! Calculate KTROPOP using std test: going up from 400hPa, 1st level
! where dT/dz > -2K/km . For lat<-60, tropo not higher than 250hPa
!-----------------------------------------------------------------------
      IF( PRSF1(JL,JK) < 40000._JPRB .and. KTROPOP(JL) < 0_JPIM ) THEN
        Zthkness = 0.5*(PTP(JL,JK+1)+PTP(JL,JK))*287./9.806* &
&               LOG(PRSF1(JL,JK+1)/PRSF1(JL,JK))  ! km, from routine ALTITUDE
        ZGRADT = ( PTP(JL,JK) - PTP(JL,JK+1) ) / (1.e-3*Zthkness)   ! K/km
        IF( ZGRADT > -2._JPRB ) THEN
          KTROPOP(JL) = JK
          IF( PLAT(JL) <= -60._JPRB .and. PRSF1(JL,JK+1) < 25000._JPRB) THEN
            KTROPOP(JL) = -9_JPIM
            DO JK2 = KLEV-1, 1, -1
              IF( KTROPOP(JL)<0 .and. PRSF1(JL,JK2)<25000._JPRB .and. &
               &                       PRSF1(JL,JK2+1) >= 25000._JPRB ) THEN
                  KTROPOP(JL) = JK2
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO

!-----------------------------------------------------------------------
! If tropopause not found or not realistic, use lat-dep vector from SOCRATES
!   Keep diags: this happens *rarely* !! If it happens often, something is wrong
!-----------------------------------------------------------------------
    JK_TROP = KTROPOP(JL)
    ZKM = PGEOH(JL,JK_TROP)  * ZRGI *1e-3 ! height in km
    IF( JK_TROP < 0_JPIM .or. ZKM < 6.6_JPRB .or.  ZKM> 20.5_JPRB ) THEN

      JLAT_SOC = NINT((PLAT(JL) + 90.0_JPRB ) / (180./NLAT_PTROPO) + 1.0_JPRB)
      JLAT_SOC = MAX(1_JPIM,MIN(NLAT_PTROPO,JLAT_SOC))
      JTROPOP_SOC  = -9
      DO JK = KLEV-1, 1, -1
        IF( JTROPOP_SOC<0 .AND. PRSF1(JL,JK)   <  PTROP_SOC_B(JLAT_SOC) .and. &
&                               PRSF1(JL,JK+1) >= PTROP_SOC_B(JLAT_SOC) ) JTROPOP_SOC = JK
      ENDDO

!             print*,'TROPOPAUSE warning: not found/irrealistic'
!             WRITE(NULOUT,'(a,2f9.2,a,i3)') '   at (lat)= ',JLAT(JL) &
!     &          ,' : KTROPOP= ',JK_TROP
!             WRITE(NULOUT,'(3(a,f12.5))')'   ZKM= ',ZKM &
!     &          ,' ; PRSF1(KTROPOP)= ',0.01*PRSF1(JL,JK_TROP) &
!     &          ,' with ZGRADT= ',ZGRADT
!             WRITE(NULOUT,*)'  Using SOCRATES tropopause pressure instead: ' &
!     &                 ,' KTROPOP= ',JTROPOP_SOC

      KTROPOP(JL) = JTROPOP_SOC

    ENDIF
  ENDDO

IF (LHOOK) CALL DR_HOOK('BASCOE_TROPOPAUSE',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_TROPOPAUSE

