! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!-----------------------------------------------------------------------
!
      SUBROUTINE SECONDHH_GEN

!----------------------------------------------------------------

!**** *SECONDHH_GEN* - COMPUTATION OF TABLES FOR SECOND ORDER INTERACTION
!                      FOR FINITE DEPTH WATER.

!     P.A.E.M. JANSSEN  JULY 2010

!     PURPOSE.
!     ---------

!          INITIALIZE ALL SECOND HARMONICS AND WAVE SET UP/ SET DOWN
!          INTERACTIONS COEFFICIENTS. PRODUCES PART 3 OF MODULE YOWTABL

!**   INTERFACE.
!     ----------

!          *CALL* *SECONDHH_GEN*

!     METHOD.
!     -------

!          SEE REFERENCE.

!     EXTERNALS.
!     ----------

!         TABLES_2ND

!     REFERENCES.
!     -----------

!          P.A.E.M. JANSSEN (2009), JFM

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DELTH    ,TH       ,FRATIO 
      USE YOWTABL  , ONLY : MR, XMR, MA, XMA, NFREH, NANGH, NMAX,       &
     &                      OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M,      &
     &                      TA, TB, TC_QL, TT_4M, TT_4P, TFAKH
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWPCONS , ONLY : G, PI, ZPI
      USE YOWSHAL  , ONLY : NDEPTH
      USE YOWTEST  , ONLY : IU06
!-----------------------------------------------------------------------

      IMPLICIT NONE
#include "tables_2nd.intfb.h"

      INTEGER(KIND=JWIM) :: M,K,K0,M0,MP

      REAL(KIND=JWRB) :: FRAC,CO1,OMSTART,OM0,OM1,XK0,XK1,XK2,A,B,C_QL
!
!-----------------------------------------------------------------------
! 
!***  1. INITIALISE PARAMETERS.
!     ------------------------
!      
      NFREH = NFRE/2
      NANGH = NANG/2

      FRAC = FRATIO-1.0_JWRB
      OMSTART = ZPI*FR(1)
      MR = NFRE/NFREH
      XMR = 1.0_JWRB/REAL(MR,JWRB)
      MA = NANG/NANGH
      XMA = 1.0_JWRB/REAL(MA,JWRB)
      DELTHH = REAL(MA,JWRB)*DELTH

      IF (.NOT. ALLOCATED(OMEGA)) ALLOCATE(OMEGA(NFREH))
      IF (.NOT. ALLOCATED(THH)) ALLOCATE(THH(NANGH))
      IF (.NOT. ALLOCATED(DFDTH)) ALLOCATE(DFDTH(NFREH))

      DO M=1,NFREH
         OMEGA(M) = ZPI*FR(MR*M)
      ENDDO

      DO K=1,NANGH
         K0 = MA*K+1
         IF (K0 > NANG) K0 = K0-NANG
         THH(K) = TH(K0)
      ENDDO

      CO1   = 1.0_JWRB/2.0_JWRB*DELTHH/ZPI
      DFDTH(1) = CO1*(OMEGA(2)-OMEGA(1))
      DO M=2,NFREH-1
         DFDTH(M)=CO1*(OMEGA(M+1)-OMEGA(M-1))
      ENDDO
      DFDTH(NFREH)=CO1*(OMEGA(NFREH)-OMEGA(NFREH-1))

      NMAX = 1+XMR*(1+NINT(LOG(2._JWRB*OMEGA(NFREH)/OMSTART)/LOG(1._JWRB+FRAC)))
!
!***  2.INITIALISE TABLES
!     -------------------
!
      IF (.NOT.ALLOCATED(TA)) ALLOCATE(TA(NDEPTH,NANGH,NFREH,NFREH))
      IF (.NOT.ALLOCATED(TB)) ALLOCATE(TB(NDEPTH,NANGH,NFREH,NFREH))
      IF (.NOT.ALLOCATED(TC_QL)) ALLOCATE(TC_QL(NDEPTH,NANGH,NFREH,NFREH))
      IF (.NOT.ALLOCATED(TT_4M)) ALLOCATE(TT_4M(NDEPTH,NANGH,NFREH,NFREH))
      IF (.NOT.ALLOCATED(TT_4P)) ALLOCATE(TT_4P(NDEPTH,NANGH,NFREH,NFREH))
      IF (.NOT.ALLOCATED(IM_P)) ALLOCATE(IM_P(NFREH,NFREH))
      IF (.NOT.ALLOCATED(IM_M)) ALLOCATE(IM_M(NFREH,NFREH)) 
      IF (.NOT.ALLOCATED(TFAKH)) ALLOCATE(TFAKH(NFREH,NDEPTH))

      CALL TABLES_2ND(NFREH,NANGH,NDEPTH,OMSTART,FRAC,XMR,              &
     &                DFDTH,OMEGA,THH,TA,TB,TC_QL,TT_4M,TT_4P,          &
     &                IM_P,IM_M,TFAKH)

! ----------------------------------------------------------------------

!***  3. PRINTER PROTOCOL
!         ---------------

      WRITE (IU06,'('' '')')
      WRITE (IU06,'('' CHECK TABLES SECOND-ORDER SPECTRUM'')')
      WRITE (IU06,'(''  NUMBER OF FREQUENCIES NFREH = '',I3)') NFREH
      WRITE (IU06,'(''  NUMBER OF DIRECTIONS NANGH = '',I3)') NANGH
      WRITE (IU06,'(''  MAX FREQUENCY = '',I3)') NMAX
      WRITE (IU06,'(''  THINNING PARAMETER MR = '',I3)') MR
      WRITE (IU06,'(''  THINNING PARAMETER MA = '',I3)') MA

      WRITE (IU06,'('' '')')
      WRITE(IU06,'('' DEEP WATER TABLES FOR SECOND-ORDER '')')
      WRITE(IU06,'(''    M        A NUM           A THEO '')')
      DO M=1,NFREH
         OM0 = OMEGA(M)
         OM1 = OMEGA(2)
         IF (OM1 < OM0/2.0_JWRB) THEN
            XK1 = OM1**2/G
            XK2 = (OM0-OM1)**2/G
            A = (ABS(XK1+XK2)/2.0_JWRB)**2
            WRITE(IU06,'(I5,2F16.9)') M,TA(NDEPTH,NANGH,2,M)/DFDTH(2),A
         ENDIF
      ENDDO 

      WRITE (IU06,'('' '')')
      WRITE(IU06,'(''    M        B NUM           B THEO '')')
      DO M=1,NFREH
         OM0 = OMEGA(M)
         OM1 = OMEGA(M)
         XK1 = OM1**2/G
         XK2 = (OM0+OM1)**2/G
         B = (ABS(XK1-XK2)/2.)**2
         WRITE(IU06,'(I5,2F16.9)') M,TB(NDEPTH,NANGH,M,M)/DFDTH(M),B  
      ENDDO 

      WRITE (IU06,'('' '')')
      WRITE(IU06,'(''    M        C NUM           C THEO '')')
      DO M=1,NFREH
         OM0 = OMEGA(M)
         XK0 = OM0**2/G
         C_QL = -XK0**2
         WRITE(IU06,'(I5,2F16.9)') M,TC_QL(NDEPTH,NANGH,M,M)/DFDTH(M),C_QL 
      ENDDO

      CALL FLUSH (IU06)
!
!----------------------------------------------------------------------
!
      END SUBROUTINE SECONDHH_GEN 
