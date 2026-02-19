! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SOACONV &
  &( YDRIP, YDEAERSNK,KIDIA   , KFDIA    , KLON     , KLEV  , KLEVTROP, &
  &  PTP, PRHO, PGELAM, PGELAT, PTSPHY, POM, PSOAA, PSOAB, PAERGASA, PAERGASB,&
  & PTAERGASA, PTAERSOAA,PTAERGASB, PTAERSOAB)

!*** * AER_SOACONV Conversion of SOA precursor gas to SOA 

!**   INTERFACE.
!     ----------
!          *AER_SOACONV* IS CALLED FROM *AER_PHY3*.

!     AUTHOR.
!     -------
!        SAMUEL REMY HYGEOS

!     SOURCE.
!     -------

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 20190712
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEAERSNK ,ONLY : TEAERSNK
USE YOMRIP   , ONLY : TRIP

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
TYPE(TEAERSNK)    ,INTENT(INOUT) :: YDEAERSNK
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVTROP(KLON)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PTP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: POM(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOAA(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOAB(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERGASA(KLON,KLEV),PAERGASB(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAERGASA(KLON,KLEV), PTAERSOAA(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAERGASB(KLON,KLEV), PTAERSOAB(KLON,KLEV)
REAL(KIND=JPRB), DIMENSION(3),PARAMETER :: ZKOM_REF = (/ &
  & 1.0_JPRB , &      !  OH / O3 + Terpene reaction product; OH + ISOP. Net K of unity might be on high side.
  & 0.01_JPRB, &      !  Anthro SOA, Low K (high C*)
  & 1.0_JPRB  /)      !  Aged Anthro SOA, high K (low C*)
REAL(KIND=JPRB), PARAMETER    :: ZRT310=1._JPRB/310_JPRB
REAL(KIND=JPRB), PARAMETER    :: ZRELAX = 0.2 ! Relaxation parameter to obtainmore smooth solution
REAL(KIND=JPRB)               :: ZRPT, ZT_310,ZT_295, ZTMP2_A, ZTMP2_B,ZMT_ORGAER
! Heat of vaporization over R (Takekawa, AE 2003)
! is also dependant on SOA type..
REAL(KIND=JPRB), DIMENSION(3), PARAMETER    :: ZHEAT_VAPORR = (/ &
  & 1.0E4_JPRB, 5.E3_JPRB,5.E3_JPRB /)
REAL(KIND=JPRB)  :: ZSOAA(KLON,KLEV)
REAL(KIND=JPRB)  :: ZSOAB(KLON,KLEV)
REAL(KIND=JPRB)  :: ZSOGA(KLON,KLEV)
REAL(KIND=JPRB)  :: ZSOGB(KLON,KLEV)
REAL(KIND=JPRB)  :: ZORG_GAS_A, ZORG_GAS_B, ZORG_AER_A, ZORG_AER_B



INTEGER(KIND=JPIM), PARAMETER ::IITER = 3


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JK, ITER
REAL(KIND=JPRB) ::  ZTSCALI, ZKOM_A, ZKOM_B
! Diurnal cycle

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SOACONV',0,ZHOOK_HANDLE)
ASSOCIATE(RSO2CV1=>YDEAERSNK%RSO2CV1,RSO2CV2=>YDEAERSNK%RSO2CV2)

ZTSCALI =  1./(1._JPRB+2.315E-6_JPRB*PTSPHY)
ZSOAA(KIDIA:KFDIA,1:KLEV)=PSOAA(KIDIA:KFDIA,1:KLEV)*PRHO(KIDIA:KFDIA,1:KLEV)*1E9_JPRB
ZSOAB(KIDIA:KFDIA,1:KLEV)=PSOAB(KIDIA:KFDIA,1:KLEV)*PRHO(KIDIA:KFDIA,1:KLEV)*1E9_JPRB
ZSOGA(KIDIA:KFDIA,1:KLEV)=PAERGASA(KIDIA:KFDIA,1:KLEV)*PRHO(KIDIA:KFDIA,1:KLEV)*1E9_JPRB
ZSOGB(KIDIA:KFDIA,1:KLEV)=PAERGASB(KIDIA:KFDIA,1:KLEV)*PRHO(KIDIA:KFDIA,1:KLEV)*1E9_JPRB
DO JL=KIDIA,KFDIA
  DO JK=KLEVTROP(JL),KLEV
    !Reciprocal temperature [1/K]
    ZRPT = 1.0_JPRB / PTP(JL,JK)

    ZTMP2_A = EXP( ZHEAT_VAPORR(3) * ( ZRPT - ZRT310 ) )
    ZTMP2_B = EXP( ZHEAT_VAPORR(1) * ( ZRPT - ZRT310 ) )
    ZT_295 = PTP(JL,JK) / 295_JPRB
    ZT_310 = PTP(JL,JK) / 310_JPRB
    ZKOM_A = ZKOM_REF(3) * ZT_295 * ZTMP2_A 
    ZKOM_B = ZKOM_REF(1) * ZT_310 * ZTMP2_B 
    ! First estimate of total organic aerosol available for condensation
    ! of additional SOA
    ZMT_ORGAER=(POM(JL,JK)+PSOAA(JL,JK)+PSOAB(JL,JK))*PRHO(JL,JK)*1E9_JPRB
    ! do several iterations
    DO ITER =1,IITER

     ! Compute the heat-of-vaporization exponential term
     ZORG_AER_A = ZKOM_A*ZMT_ORGAER / &
      &  (1._JPRB + ZKOM_A * ZMT_ORGAER ) * (ZSOGA(JL,JK) + ZSOAA(JL,JK))
    ! Apply some relaxation towards original concentration,
     ZORG_AER_A = ZRELAX*ZSOAA(JL,JK) + (1._JPRB-ZRELAX)*ZORG_AER_A

    ! Limit conversion of aerosol back to gas conversion to 1%  at most for
    ! each time step.
     ZORG_AER_A =MAX(PSOAA(JL,JK)*PRHO(JL,JK)*1E9_JPRB*0.99_JPRB,ZORG_AER_A)
    ! All that is not in aerosol phase is in gas-phase:
     ZORG_GAS_A = ( ZSOGA(JL,JK)+ZSOAA(JL,JK) ) - ZORG_AER_A
    
      ZORG_AER_B = ZKOM_B*ZMT_ORGAER / &
      &  (1._JPRB + ZKOM_B * ZMT_ORGAER ) * (ZSOGB(JL,JK) + ZSOAB(JL,JK))
    ! Apply some relaxation towards original concentration,
     ZORG_AER_B = ZRELAX*ZSOAB(JL,JK) + (1._JPRB-ZRELAX)*ZORG_AER_B

    ! Limit conversion of aerosol back to gas conversion to 1%  at most for
    ! each time step.
     ZORG_AER_B =MAX(PSOAB(JL,JK)*PRHO(JL,JK)*1E9_JPRB*0.99,ZORG_AER_B)
    ! All that is not in aerosol phase is in gas-phase:
     ZORG_GAS_B = ( ZSOGB(JL,JK)+ZSOAB(JL,JK) ) - ZORG_AER_B

     ZSOGA(JL,JK)=ZORG_GAS_A
     ZSOAA(JL,JK)=ZORG_AER_A
     ZSOGB(JL,JK)=ZORG_GAS_B
     ZSOAB(JL,JK)=ZORG_AER_B
    ENDDO

    ! Limit SOA to SOA+SOG
    ZSOAA(JL,JK)=MIN(ZSOAA(JL,JK),PSOAA(JL,JK)*PRHO(JL,JK)*1E9_JPRB+PAERGASA(JL,JK)*PRHO(JL,JK)*1E9_JPRB)
    ZSOAB(JL,JK)=MIN(ZSOAB(JL,JK),PSOAB(JL,JK)*PRHO(JL,JK)*1E9_JPRB+PAERGASB(JL,JK)*PRHO(JL,JK)*1E9_JPRB)
    ! Compute SOG to be consistent with SOA_in+SOG_in - SOA_out
    ZSOGA(JL,JK)=PSOAA(JL,JK)*PRHO(JL,JK)*1E9_JPRB+PAERGASA(JL,JK)*PRHO(JL,JK)*1E9_JPRB- ZSOAA(JL,JK)
    ZSOGB(JL,JK)=PSOAB(JL,JK)*PRHO(JL,JK)*1E9_JPRB+PAERGASB(JL,JK)*PRHO(JL,JK)*1E9_JPRB - ZSOAB(JL,JK)

    PTAERSOAA(JL,JK)= (ZTSCALI*ZSOAA(JL,JK)*1E-9/PRHO(JL,JK) - PSOAA(JL,JK))/PTSPHY
    PTAERSOAB(JL,JK)= (ZTSCALI*ZSOAB(JL,JK)*1E-9/PRHO(JL,JK) - PSOAB(JL,JK))/PTSPHY
    PTAERGASA(JL,JK)= (ZTSCALI*ZSOGA(JL,JK)*1E-9/PRHO(JL,JK) - PAERGASA(JL,JK))/PTSPHY
    PTAERGASB(JL,JK)= (ZTSCALI*ZSOGB(JL,JK)*1E-9/PRHO(JL,JK) - PAERGASB(JL,JK))/PTSPHY
  ENDDO
  DO JK=1,KLEVTROP(JL)-1
    PTAERSOAA(JL,JK)= (ZTSCALI*ZSOAA(JL,JK)*1E-9/PRHO(JL,JK) - PSOAA(JL,JK))/PTSPHY
    PTAERSOAB(JL,JK)= (ZTSCALI*ZSOAB(JL,JK)*1E-9/PRHO(JL,JK) - PSOAB(JL,JK))/PTSPHY
    PTAERGASA(JL,JK)= (ZTSCALI*ZSOGA(JL,JK)*1E-9/PRHO(JL,JK) - PAERGASA(JL,JK))/PTSPHY
    PTAERGASB(JL,JK)= (ZTSCALI*ZSOGB(JL,JK)*1E-9/PRHO(JL,JK) - PAERGASB(JL,JK))/PTSPHY
  ENDDO
ENDDO
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SOACONV',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SOACONV
