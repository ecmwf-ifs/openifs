! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

       SUBROUTINE SMOOTHSARSPEC(LSPEC,NSPEC,NANG,NFRE,SPEC) 
 
!---------------------------------------------------------------------

!     PURPOSE
!     -------
!     TO SMOOTH SPECTRA THAT ARE USED IN THE PARTITIONING OF THE SAR
!     ASSIMILATION.

!     AUTHOR
!     ------
!     J. BIDLOT  ECMWF APRIL 2001.

!     EXTERNALS
!     ---------

!**********************************************************************

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : DFIM
      USE YOWPCONS , ONLY : EPSMIN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: LSPEC,NSPEC,NANG,NFRE
      REAL(KIND=JWRB), INTENT(INOUT) :: SPEC(LSPEC:NSPEC,NANG,NFRE)

      INTEGER(KIND=JWIM) :: ISPEC,K,M,KM,KP

      REAL(KIND=JWRB) :: C16
      REAL(KIND=JWRB) :: WORK(LSPEC:NSPEC,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(LSPEC:NSPEC) ::  EM, EW, TEMP

!     COMPUTE THE ENERGY OF THE ORIGINAL SPECTRA

      DO ISPEC = LSPEC,NSPEC 
        EM(ISPEC) = EPSMIN 
      ENDDO

      DO M=1,NFRE
        DO ISPEC = LSPEC,NSPEC 
          TEMP(ISPEC) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO ISPEC = LSPEC,NSPEC 
            TEMP(ISPEC) = TEMP(ISPEC)+SPEC(ISPEC,K,M)
          ENDDO
        ENDDO
        DO ISPEC = LSPEC,NSPEC 
          EM(ISPEC) = EM(ISPEC)+DFIM(M)*TEMP(ISPEC)
        ENDDO
      ENDDO


!     SMOOTH SPECTRA

      DO K=1,NANG
        KM  = K - 1
        IF (KM .LT. 1) KM = NANG
        KP = K + 1
        IF (KP.GT. NANG) KP = 1

        DO M=2,NFRE-1
          DO ISPEC = LSPEC,NSPEC 
            WORK(ISPEC,K,M) =       SPEC(ISPEC,KM,M-1) +                &
     &                     2.0_JWRB*SPEC(ISPEC,KM,M)   +                &
     &                              SPEC(ISPEC,KM,M+1) +                &
     &                     2.0_JWRB*SPEC(ISPEC,K,M-1)  +                &
     &                     4.0_JWRB*SPEC(ISPEC,K,M)    +                &
     &                     2.0_JWRB*SPEC(ISPEC,K,M+1)  +                &
     &                              SPEC(ISPEC,KP,M-1) +                &
     &                     2.0_JWRB*SPEC(ISPEC,KP,M)   +                &
     &                              SPEC(ISPEC,KP,M+1)
          ENDDO
        ENDDO

!       We assume zero energy for frequency below M=1
        M=1
          DO ISPEC = LSPEC,NSPEC 
            WORK(ISPEC,K,M) =                                           &
     &                     2.0_JWRB*SPEC(ISPEC,KM,M)   +                &
     &                              SPEC(ISPEC,KM,M+1) +                &
     &                     4.0_JWRB*SPEC(ISPEC,K,M)    +                &
     &                     2.0_JWRB*SPEC(ISPEC,K,M+1)  +                &
     &                     2.0_JWRB*SPEC(ISPEC,KP,M)   +                &
     &                              SPEC(ISPEC,KP,M+1)
          ENDDO
      ENDDO

      C16=1.0_JWRB/16.0_JWRB

      DO K=1,NANG
        KM  = K - 1
        IF (KM .LT. 1) KM = NANG
        KP = K + 1
        IF (KP.GT. NANG) KP = 1

        DO M=1,NFRE-1
          DO ISPEC = LSPEC,NSPEC 
            SPEC(ISPEC,K,M) = C16 * WORK(ISPEC,K,M)
          ENDDO
        ENDDO
      ENDDO

!     COMPUTE THE ENERGY OF THE SMOOTHED SPECTRA

      DO ISPEC = LSPEC,NSPEC 
        EW(ISPEC) = EPSMIN 
      ENDDO

      DO M=1,NFRE
        DO ISPEC = LSPEC,NSPEC 
          TEMP(ISPEC) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO ISPEC = LSPEC,NSPEC 
            TEMP(ISPEC) = TEMP(ISPEC)+SPEC(ISPEC,K,M)
          ENDDO
        ENDDO
        DO ISPEC = LSPEC,NSPEC 
          EW(ISPEC) = EW(ISPEC)+DFIM(M)*TEMP(ISPEC)
        ENDDO
      ENDDO

!     ENSURE THE ENERGY IS CONSERVED BETWEEN ORIGINAL AND SMOOTHED
!     SPECTRA.

      DO M=1,NFRE
        DO K=1,NANG
          DO ISPEC = LSPEC,NSPEC 
            SPEC(ISPEC,K,M) = (EM(ISPEC)/EW(ISPEC))*SPEC(ISPEC,K,M)
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE SMOOTHSARSPEC 
