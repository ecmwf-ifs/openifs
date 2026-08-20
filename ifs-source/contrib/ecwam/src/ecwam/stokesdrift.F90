! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES)
 
!
!***  *STOKESDRIFT*   DETERMINES THE STOKES DRIFT
!
!     PETER JANSSEN MARCH 2009
!
!     PURPOSE.
!     --------
!
!              DETERMINATION OF STOKES DRIFT VECTOR
!
!     INTERFACE.
!     ----------
!              *CALL*  *STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE,WDWAVE,CICOVER,USTOKES,VSTOKES)*
!
!                       INPUT:
!                            *KIJS*   - FIRST GRIDPOINT
!                            *KIJL*   - LAST GRIDPOINT
!                            *FL1*    - 2-D SPECTRUM
!                            *STOKFAC*- FACTOR TO COMPUTE THE STOKES DRIFT
!                            Auxilliary fields to specify Stokes when model sea ice cover the blocking threshold
!                            as 0.016*WSWAVE, aligned in the wind direction
!                            *WSWAVE* - WIND SPEED IN M/S.
!                            *WDWAVE* - WIND DIRECTION IN RADIANS.
!                            *CICOVER*- SEA ICE COVER.
!
!                       OUTPUT: 
!                            *USTOKES*   - U-COMPONENT STOKES DRIFT
!                            *VSTOKES*   - V-COMPONENT STOKES DRIFT
!
!     METHOD.
!     -------
!              DETERMINE U- AND V-COMPONENT OF STOKES DRIFT FOLLOWING
!              K.E. KENYON, J.G.R., 74, 6991-6994
!
!     EXTERNALS.
!     ----------
!              NONE
!
!
!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,    &
     &                      DFIM_SIM ,FRATIO   ,COSTH    ,SINTH
      USE YOWICE   , ONLY : LICERUN  ,LWAMRSETCI, CITHRSH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
       
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: STOKFAC 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, WDWAVE, CICOVER
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: USTOKES, VSTOKES


      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB), PARAMETER :: STMAX=1.5_JWRB   ! maximum magnitude (this is for safety when coupled)
      REAL(KIND=JWRB) :: CONST, FAC, FAC1, FAC2, FAC3
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: STFAC

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STOKESDRIFT',0,ZHOOK_HANDLE)

 
!***  1. DETERMINE STOKE DRIFT VECTOR.
!     --------------------------------

      CONST = 2.0_JWRB*DELTH*ZPI**3/G*FR(NFRE_ODD)**4

!***  1.1 PERFORM INTEGRATION.
!     ------------------------
 
      DO IJ = KIJS,KIJL
         USTOKES(IJ) = 0.0_JWRB
         VSTOKES(IJ) = 0.0_JWRB
      ENDDO

      DO M=1,NFRE_ODD
         DO IJ = KIJS,KIJL
           STFAC(IJ) = STOKFAC(IJ,M)*DFIM_SIM(M)
         ENDDO
         DO K=1,NANG
            DO IJ = KIJS,KIJL
               FAC3 = STFAC(IJ)*FL1(IJ,K,M)
               USTOKES(IJ) = USTOKES(IJ)+FAC3*SINTH(K)
               VSTOKES(IJ) = VSTOKES(IJ)+FAC3*COSTH(K)
            ENDDO
         ENDDO
      ENDDO
 
!***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
!     -----------------------------------------
 
      DO K=1,NANG
         FAC1 = CONST*SINTH(K)
         FAC2 = CONST*COSTH(K)
         DO IJ = KIJS,KIJL
            USTOKES(IJ) = USTOKES(IJ)+FAC1*FL1(IJ,K,NFRE_ODD)
            VSTOKES(IJ) = VSTOKES(IJ)+FAC2*FL1(IJ,K,NFRE_ODD)
         ENDDO
      ENDDO


!***  1.3 Sea Ice exception
!     ---------------------
      IF (LICERUN .AND. LWAMRSETCI) THEN
       DO IJ=KIJS,KIJL
         IF (CICOVER(IJ) > CITHRSH) THEN
           USTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*SIN(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
           VSTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*COS(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
         ENDIF
       ENDDO
     ENDIF

!***  1.4 Protection
!     --------------

      DO IJ = KIJS,KIJL
         USTOKES(IJ) = MIN(MAX(USTOKES(IJ),-STMAX),STMAX)
         VSTOKES(IJ) = MIN(MAX(VSTOKES(IJ),-STMAX),STMAX)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('STOKESDRIFT',1,ZHOOK_HANDLE)

      END SUBROUTINE STOKESDRIFT
