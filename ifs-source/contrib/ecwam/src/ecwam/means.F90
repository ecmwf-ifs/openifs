! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MEANS(                                                 &
!                   DIMENSIONS
     &                NSPEC , MPART , NPART , NANG, NFRE ,              &
!                   INPUT
     &                SPEC , PART ,                                     &
!                   OUTPUT
     &                MEANE , MEANANG , MEANFRE , PKFRE ,               &
!                   WORK SPACE.
     &                SUM0 , SUMX0 , SUMY0 )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     COMPUTES MEAN DIRECTION AND MEAN FREQUENCY OF 2D SPECTRA.
!     AND THE PEAK FREQUENCY
!     FOR THE CALCULATION OF THE SPECTRAL PARTITIONING.

!-----------------------------------------------------------------------

!     MODULE :
!     --------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,    &
     &           COSTH    ,SINTH     ,FRTAIL

!-----------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER(KIND=JWIM) :: NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!                 ARRAY DIMENSIONS.
      REAL(KIND=JWRB) :: SPEC(0:NSPEC,NANG,NFRE)
!                 2D SPECTRA.
      INTEGER(KIND=JWIM) :: PART(NSPEC,NANG,NFRE)
!          PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
      REAL(KIND=JWRB) :: MEANE(NSPEC,MPART),                            &
     &     MEANANG(NSPEC,MPART), MEANFRE(NSPEC,MPART)
!                 INTEGRATED VALUES OF PARTITIONINGS.
      REAL(KIND=JWRB) ::  PKFRE(NSPEC,MPART)
!                 PEAK FREQUENCY OF PARTITIONINGS.
      REAL(KIND=JWRB) :: SUM0(NSPEC,MPART),                             &
     &     SUMX0(NSPEC,MPART), SUMY0(NSPEC,MPART)
!                 FOR TEMPORARY SUMS

!     LOCAL VARIABLES:
!     ----------------

      INTEGER(KIND=JWIM) :: ISPEC, IPART, IANG, IFRE
!                   LOOP INDEXES.
      REAL(KIND=JWRB) :: FACTOR_TAIL

!-----------------------------------------------------------------------

!     1. COMPUTE MEAN DIRECTIONS.
!     --------------------------

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          SUMX0(ISPEC,IPART) = 0.0_JWRB
          SUMY0(ISPEC,IPART) = 0.0_JWRB
        ENDDO
      ENDDO

      DO IANG = 1,NANG
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,IPART) = 0.0_JWRB
          ENDDO
        ENDDO
        DO IFRE=1,NFRE
          DO ISPEC=1,NSPEC
            SUM0(ISPEC,PART(ISPEC,IANG,IFRE))                           &
     &          = SUM0(ISPEC,PART(ISPEC,IANG,IFRE))                     &
     &          + SPEC(ISPEC,IANG,IFRE) * DFIM(IFRE)
          ENDDO
        ENDDO
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            SUMX0(ISPEC,IPART) = SUMX0(ISPEC,IPART)                     &
     &          + COSTH(IANG) * SUM0(ISPEC,IPART)
            SUMY0(ISPEC,IPART) = SUMY0(ISPEC,IPART)                     &
     &          + SINTH(IANG) * SUM0(ISPEC,IPART)
          ENDDO
        ENDDO
      ENDDO

      DO ISPEC = 1,NSPEC
        DO IPART = 1,NPART(ISPEC)
          MEANANG(ISPEC,IPART)                                          &
     &        = ATAN2(SUMY0(ISPEC,IPART),SUMX0(ISPEC,IPART))
        ENDDO
      ENDDO

!----------------------------------------------------------------------

!    2. COMPUTE MEAN FREQUENCY.
!    --------------------------

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          MEANFRE(ISPEC,IPART) = 0.0_JWRB
        ENDDO
      ENDDO

      DO IFRE = 1,NFRE
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,IPART) = 0.0_JWRB
          ENDDO
        ENDDO
        DO IANG = 1,NANG
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,PART(ISPEC,IANG,IFRE))                           &
     &          = SUM0(ISPEC,PART(ISPEC,IANG,IFRE))                     &
     &          + SPEC(ISPEC,IANG,IFRE)
          ENDDO
        ENDDO
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            MEANFRE(ISPEC,IPART) = MEANFRE(ISPEC,IPART)                 &
     &          + DFIMOFR(IFRE)  * SUM0(ISPEC,IPART)
          ENDDO
        ENDDO
      ENDDO

!     ADD TAIL CONTRIBUTION.

      FACTOR_TAIL = FRTAIL*DELTH

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          MEANFRE(ISPEC,IPART) = MEANFRE(ISPEC,IPART)                   &
     &        + FACTOR_TAIL * SUM0(ISPEC,IPART)
          IF (MEANFRE(ISPEC,IPART).GT.0.0_JWRB) THEN
            MEANFRE(ISPEC,IPART)                                        &
     &          = MEANE(ISPEC,IPART) / MEANFRE(ISPEC,IPART)
          ELSE
            MEANFRE(ISPEC,IPART)=0.0_JWRB
          ENDIF
          MEANFRE(ISPEC,IPART) = MAX(MEANFRE(ISPEC,IPART),FR(1))
        ENDDO
      ENDDO

!----------------------------------------------------------------------

!    3. COMPUTE PEAK FREQUENCY.
!    --------------------------

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          PKFRE(ISPEC,IPART) = FR(NFRE) 
          SUM0(ISPEC,IPART) = 0.0_JWRB
        ENDDO
      ENDDO

      DO IFRE = 1,NFRE
        DO IANG = 1,NANG
          DO ISPEC = 1,NSPEC
            IPART=PART(ISPEC,IANG,IFRE)
            IF (SPEC(ISPEC,IANG,IFRE) .GT. SUM0(ISPEC,IPART)) THEN 
              PKFRE(ISPEC,IPART) = FR(IFRE) 
              SUM0(ISPEC,IPART) = SPEC(ISPEC,IANG,IFRE) 
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE MEANS
