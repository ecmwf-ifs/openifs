! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE MTABS

! ----------------------------------------------------------------------

!**** *MTABS* - ROUTINE TO COMPUTE TABLES USED FOR SHALLOW WATER.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO COMPUTE TABLES USED FOR SHALLOW WATER.

!**   INTERFACE.
!     ----------

!       *CALL* *MTABS*

!     METHOD.
!     -------

!       TABLES FOR GROUP VELOCITY, WAVE NUMBER AND OMEGA/SINH(2KD)
!       ARE YOWPUTED AT ALL FREQUENCIES AND FOR A DEPTH TABLE
!       OF LENGTH NDEPTH, STARTING AT DEPTHA METERS AND INCREMENTED
!        BY DEPTHD METRES.

!     EXTERNALS.
!     ----------

!       *AKI*       - FUNCTION TO COMPUTE WAVE NUMBER.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR
      USE YOWPARAM , ONLY : NFRE     ,NIBLO
      USE YOWPCONS , ONLY : G        ,PI       ,ZPI
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTHA   ,DEPTHD   ,   &
     &            TCGOND   ,TFAK     ,TSIHKD   ,TFAC_ST
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "aki.intfb.h"

      INTEGER(KIND=JWIM) :: M, JD, NAN, NSTP

      REAL(KIND=JWRB) :: GH, OM, AD, AK, AKD, DEPTHE

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      ALLOCATE(TCGOND(NDEPTH,NFRE))
      ALLOCATE(TFAK(NDEPTH,NFRE))
      ALLOCATE(TSIHKD(NDEPTH,NFRE))
      ALLOCATE(TFAC_ST(NDEPTH,NFRE))

!*    1. GROUP VELOCITY AND WAVE NUMBER.
!        -------------------------------

!*    1.1 LOOP OVER FREQUENCIES.
!         ----------------------

      GH = G/(4.0_JWRB*PI)
      DO M=1,NFRE
        OM=ZPI*FR(M)

!*    1.1.1 LOOP OVER DEPTH.
!           -----------------

        DO JD=1,NDEPTH
          AD=DEPTHA*DEPTHD**(JD-1)
          AK=AKI(OM,AD)
          TFAK(JD,M)=AK
          AKD=AK*AD
          IF (AKD.LE.10.0_JWRB) THEN
            TCGOND(JD,M) = 0.5_JWRB*SQRT(G*TANH(AKD)/AK)*               &
     &                     (1.0_JWRB+2.0_JWRB*AKD/SINH(2.0_JWRB*AKD))
            TSIHKD(JD,M) = OM/SINH(2.0_JWRB*AKD)
            TFAC_ST(JD,M) = 2.0_JWRB*G*AK**2/(OM*TANH(2.0_JWRB*AKD))
          ELSE
            TCGOND(JD,M) = GH/FR(M)
            TSIHKD(JD,M) = 0.0_JWRB
            TFAC_ST(JD,M) = 2.0_JWRB/G*OM**3
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. PRINT TABLES.
!        -------------

      NAN  = NDEPTH 
      NSTP = NDEPTH/NAN
      NSTP = MAX(NSTP,1)
      DEPTHE = DEPTHA*DEPTHD**(NDEPTH-1)
      Print *, 'DEPTHA=', DEPTHA
      Print *, 'DEPTHE=', DEPTHE
      Print *, 'DEPTHD=', DEPTHD
      WRITE (IU06,'(1H1, '' SHALLOW WATER TABLES:'',/)')
      WRITE (IU06,'(''  LOGARITHMIC DEPTH FROM: DEPTHA = '',F6.1,      &
     &  '' TO DEPTHE  = '',F6.1, ''IN STEPS OF DEPTHD = '',F6.1)')     &
     &    DEPTHA, DEPTHE, DEPTHD
      WRITE (IU06,'(''  PRINTED IN STEPS OF '',I3,'' ENTRIES'',/)') NSTP
      DO JD=1,NDEPTH,NSTP
        AD=DEPTHA*DEPTHD**(JD-1)
        WRITE (IU06,'(1X,''DEPTH = '',F7.1,'' METRES '')') AD
        WRITE (IU06,'(1X,''GROUP VELOCITY IN METRES/SECOND'')')
        WRITE (IU06,'(1x,13F10.5)') (TCGOND(JD,M),M=1,NFRE)
        WRITE (IU06,'(1X,''WAVE NUMBER IN 1./METRES'')')
        WRITE (IU06,'(1x,13F10.5)') (TFAK(JD,M),M=1,NFRE)
        WRITE (IU06,'(1X,''OMEGA/SINH(2KD) IN 1./SECOND'')')
        WRITE (IU06,'(1x,13F10.5)') (TSIHKD(JD,M),M=1,NFRE)
        WRITE (IU06,'(1X,''2G K**2/(OMEGA*TANH(2KD)) IN 1./(M S)'')')
        WRITE (IU06,'(1x,13F10.5)') (TFAC_ST(JD,M),M=1,NFRE)
      ENDDO

      END SUBROUTINE MTABS
