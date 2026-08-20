! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SECONDHH

!----------------------------------------------------------------

!**** *SECONDHH* - COMPUTATION OF SECOND ORDER HARMONICS AND
!                  RELEVANT TABLES FOR THE ALTIMETER CORRECTIONS.

!     P.A.E.M. JANSSEN

!     PURPOSE.
!     ---------

!          COMPUTE SECOND HARMONICS

!**   INTERFACE.
!     ----------

!          *CALL* *SECONDHH*

!     METHOD.
!     -------

!          SEE REFERENCE.

!     EXTERNALS.
!     ----------

!         VMIN_D
!         VPLUS_D          

!     REFERENCES.
!     -----------

!          V E ZAKHAROV(1967)

!-------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DELTH    ,TH       ,FRATIO   ,    &
     &            COSTH    ,SINTH 
      USE YOWTABL  , ONLY : NFREHF   ,FAC0     ,FAC1     ,FAC2     ,    &
     &            FAC3     ,FAK      ,FRHF     ,DFIMHF
      USE YOWPARAM , ONLY : NANG
      USE YOWPCONS , ONLY : G        ,ZPI

!-------------------------------------------------------------------

      IMPLICIT NONE
#include "vmin_d.intfb.h"
#include "vplus_d.intfb.h"

      INTEGER(KIND=JWIM) :: K, M, K1, M1, K2, M2

      REAL(KIND=JWRB), PARAMETER :: DEL1=1.0E-8_JWRB
      REAL(KIND=JWRB) :: CO1
      REAL(KIND=JWRB) :: XK1, XK1SQ, XK2, XK2SQ, XK3
      REAL(KIND=JWRB) :: COSDIFF
      REAL(KIND=JWRB) :: X12, X13, X32, OM1, OM2, OM3, F1, F2, F3
      REAL(KIND=JWRB) :: VM, VP
      REAL(KIND=JWRB) :: DELOM1, DELOM2
      REAL(KIND=JWRB) :: DELOM321, DELOM312
      REAL(KIND=JWRB) :: C22, S22

      REAL(KIND=JWRB), DIMENSION(NANG,NANG,NFREHF,NFREHF) :: B

!-----------------------------------------------------------------------

!*    1. INITIALISE RELEVANT QUANTITIES.


      IF (.NOT.ALLOCATED(FAC0)) ALLOCATE(FAC0(NANG,NANG,NFREHF,NFREHF))
      IF (.NOT.ALLOCATED(FAC1)) ALLOCATE(FAC1(NANG,NANG,NFREHF,NFREHF))
      IF (.NOT.ALLOCATED(FAC2)) ALLOCATE(FAC2(NANG,NANG,NFREHF,NFREHF))
      IF (.NOT.ALLOCATED(FAC3)) ALLOCATE(FAC3(NANG,NANG,NFREHF,NFREHF))

      IF (.NOT.ALLOCATED(FAK)) ALLOCATE(FAK(NFREHF))
      IF (.NOT.ALLOCATED(FRHF)) ALLOCATE(FRHF(NFREHF))
      IF (.NOT.ALLOCATED(DFIMHF)) ALLOCATE(DFIMHF(NFREHF))


      FRHF(1)  = FR(1)
      DO M=2,NFREHF
        FRHF(M) = FRATIO*FRHF(M-1)
      ENDDO

      DO M=1,NFREHF
         FAK(M) = (ZPI*FRHF(M))**2/G
      ENDDO

      CO1 = 0.5_JWRB*(FRATIO-1.)*DELTH
      DFIMHF(1) = CO1*FRHF(1)
      DO M=2,NFREHF-1
         DFIMHF(M)=CO1*(FRHF(M)+FRHF(M-1))
      ENDDO
      DFIMHF(NFREHF)=CO1*FRHF(NFREHF-1)

      DO M2=1,NFREHF
        XK2 = FAK(M2)
        XK2SQ = FAK(M2)**2
        DO  M1=1,NFREHF
          XK1 = FAK(M1)
          XK1SQ = FAK(M1)**2
          DO K1=1,NANG
            DO K2=1,NANG
              COSDIFF = COS(TH(K1)-TH(K2))
              X12 = XK1*XK2*COSDIFF
              XK3 = XK1SQ + XK2SQ +2.0_JWRB*X12 +DEL1
              XK3 = SQRT(XK3)
              X13 = XK1SQ+X12
              X32 = X12+XK2SQ
              OM1 = SQRT(G*XK1)
              OM2 = SQRT(G*XK2)
              OM3 = SQRT(G*XK3)
              F1 = SQRT(XK1/(2.0_JWRB*OM1))
              F2 = SQRT(XK2/(2.0_JWRB*OM2))
              F3 = SQRT(XK3/(2.0_JWRB*OM3))
              VM = ZPI*VMIN_D(XK3,XK1,XK2,X13,X32,X12,OM3,OM1,OM2)
              VP = ZPI*VPLUS_D(-XK3,XK1,XK2,-X13,-X32,X12,OM3,OM1,OM2)
              DELOM1 = OM3-OM1-OM2+DEL1
              DELOM2 = OM3+OM1+OM2+DEL1
              FAC0(K1,K2,M1,M2) = -F3/(F1*F2)*(VM/(DELOM1)+VP/(DELOM2))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO M2=1,NFREHF
        XK2 = FAK(M2)
        XK2SQ = FAK(M2)**2
        DO  M1=1,NFREHF
          XK1 = FAK(M1)
          XK1SQ = FAK(M1)**2
          DO K1=1,NANG
            DO K2=1,NANG
              COSDIFF = COS(TH(K1)-TH(K2))
              X12 = XK1*XK2*COSDIFF
              XK3 = XK1SQ + XK2SQ - 2.0_JWRB*X12 + DEL1
              XK3 = SQRT(XK3)
              X13 = XK1SQ-X12
              X32 = X12-XK2SQ
              OM1 = SQRT(G*XK1)
              OM2 = SQRT(G*XK2)
              OM3 = SQRT(G*XK3)+DEL1
              F1 = SQRT(XK1/(2.0_JWRB*OM1))
              F2 = SQRT(XK2/(2.0_JWRB*OM2))
              F3 = SQRT(ABS(XK3)/(2.0_JWRB*OM3))
              VM = ZPI*VMIN_D(XK1,XK3,XK2,X13,X12,X32,OM1,OM3,OM2)
              VP = ZPI*VMIN_D(XK2,-XK3,XK1,-X32,X12,-X13,OM2,OM3,OM1)
              DELOM321 = OM3+OM2-OM1+DEL1
              DELOM312 = OM3+OM1-OM2+DEL1
              B(K1,K2,M1,M2) = -F3/(F1*F2)*(VM/(DELOM321)+VP/(DELOM312))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO M2=1,NFREHF
        XK2SQ = FAK(M2)**2
        DO M1=1,NFREHF
          XK1SQ = FAK(M1)**2
          DO K2=1,NANG
            DO K1=1,NANG
              C22 = FAC0(K1,K2,M1,M2)+B(K1,K2,M1,M2)
              S22 = B(K1,K2,M1,M2)-FAC0(K1,K2,M1,M2)
              FAC1(K1,K2,M1,M2) =                                       &
     &             (XK1SQ*COSTH(K1)**2 + XK2SQ*COSTH(K2)**2)*C22        &
     &             -FAK(M1)*FAK(M2)*COSTH(K1)*COSTH(K2)*S22
              FAC2(K1,K2,M1,M2) =                                       &
     &             (XK1SQ*SINTH(K1)**2 + XK2SQ*SINTH(K2)**2)*C22        &
     &             -FAK(M1)*FAK(M2)*SINTH(K1)*SINTH(K2)*S22
              FAC3(K1,K2,M1,M2) =                                       &
     &             (XK1SQ*SINTH(K1)*COSTH(K1) +                         &
     &              XK2SQ*SINTH(K2)*COSTH(K2))*C22                      &
     &             -FAK(M1)*FAK(M2)*COSTH(K1)*SINTH(K2)*S22
              FAC0(K1,K2,M1,M2) = C22
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!     -----------------------------------------------------------------

      END SUBROUTINE SECONDHH
