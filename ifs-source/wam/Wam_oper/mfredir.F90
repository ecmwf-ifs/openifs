! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE MFREDIR

! ----------------------------------------------------------------------

!**** *MFREDIR* - ROUTINE TO COMPUTE FREQUENCY DIRECTION CONSTANTS.

!     H. GUNTHER   ECMWF   2/4/90.
!     B. HANSEN    ECMWF   4/11/97. ACTIVATE THE ROTATION OF SPECTRA.

!*    PURPOSE.
!     --------

!       INITIATES THE FREQUENCY AND DIRECTION CONSTANTS WHICH ARE
!       SAVED IN COMMON FREDIR.

!**   INTERFACE.
!     ----------

!       *CALL* *MFREDIR*

!     METHOD.
!     -------

!       STARTING FROM THE FIRST FREQUENCY THE NEXT  ARE INCREMENTED
!       BY THE FACTOR *FRATIO* (SEE PARAMETER STATEMENT IN YOWFRED).
!       THE DIRECTIONS ARE EQUALLY DISTRIBUTED OVER THE CIRCLE
!       THE CIRCLE STARTING FROM 0.5 DELTH. 

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,    &
     &            DELTH    ,DELTR    ,TH       ,COSTH    ,SINTH    ,    &
     &            FRATIO
      USE YOWPCONS , ONLY : G        ,PI       ,ZPI      ,DEG      ,    &
     &            R
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: M, K

      REAL(KIND=JWRB) :: CO1

!*    1. FREQUENCY DEPENDENT CONSTANTS.
!        ------------------------------

!*    1.1 COMPUTE FREQUENCIES.
!         --------------------

      DO M=2,NFRE
        FR(M) = FRATIO*FR(M-1)
      ENDDO

!*    1.2 COMPUTE DEEP WATER GROUP VELOCITIES.
!        ------------------------------------

      DO M=1,NFRE
        GOM(M) = G/(4.0_JWRB*PI*FR(M))
      ENDDO

!*    1.3 COMPUTE PHASE VELOCITY IN DEEP WATER.
!         -------------------------------------

      DO M = 1,NFRE
        C(M) = G/(ZPI*FR(M))
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF DIRECTIONS, BANDWIDTH, SIN AND COS.
!        --------------------------------------------------

      DELTH = ZPI/REAL(NANG,JWRB)
      DELTR = DELTH*R
      DO K=1,NANG
!CCC        TH(K) = REAL(K-1,JWRB)*DELTH
!CCC the previous line should be used if spectra should not be rotated.
!CCC the next line should be used if rotated spectra are used
        TH(K) = REAL(K-1,JWRB)*DELTH + 0.5_JWRB*DELTH
        COSTH(K) = COS(TH(K))
        SINTH(K) = SIN(TH(K))
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTATION FREQUENCY DIRECTION AREAS
!        -------------------------------------

      CO1 = 0.5_JWRB*(FRATIO-1.0_JWRB)*DELTH
      DFIM(1)= CO1*FR(1)
      DO M=2,NFRE-1
        DFIM(M)=CO1 * (FR(M)+FR(M-1))
      ENDDO
      DFIM(NFRE)=CO1*FR(NFRE-1)

! ----------------------------------------------------------------------

!*    4. PRINTER PROTOCOL
!         ---------------

      WRITE (IU06,'(''1FREQUENCY AND DIRECTION GRID'')')
      WRITE (IU06,'(''0NUMBER OF FREQUENCIES IS  NFRE = '',I3)') NFRE
      WRITE (IU06,'(''0REDUCED NUMBER OF FREQUENCIES IS  NFRE_RED = '',I3)') NFRE_RED
      WRITE (IU06,'('' NUMBER OF DIRECTIONS  IS  NANG = '',I3)') NANG
      WRITE (IU06,'(''0MODEL FREQUENCIES IN HERTZ:'')')
      WRITE (IU06,'(1X,13F10.5)') (FR(M),M=1,NFRE)
      WRITE (IU06,'(''0MODEL FREQUENCY INTERVALLS TIMES DIRECTION'',    &
     &              '' INTERVALL IN HERTZ*RADIENS'')')
      WRITE (IU06,'(1X,13F10.5)') (DFIM(M),M=1,NFRE)
      WRITE (IU06,'(''0MODEL DEEP WATER GROUPVELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') (GOM(M),M=1,NFRE)
      WRITE (IU06,'(''0MODEL DEEP WATER PHASEVELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') (C(M),M=1,NFRE)
      WRITE (IU06,'(''0MODEL DIRECTIONS IN DEGREE'',                    &
     &              '' (CLOCKWISE FROM NORTH):'')')
      WRITE (IU06,'(1X,13F10.5)') (TH(K)*DEG,K=1,NANG)

      END SUBROUTINE MFREDIR
