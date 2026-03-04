! (C) Copyright 1998- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SET99B(TRIGS,IFAX,N,LDUSEFFT992)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM),INTENT(IN) :: N
      REAL(KIND=JPRB),INTENT(OUT) :: TRIGS(N)
      INTEGER(KIND=JPIM),INTENT(OUT) :: IFAX(*)
      LOGICAL,INTENT(OUT) :: LDUSEFFT992

      INTEGER(KIND=JPIM) :: I,IFAC,IL,IXXX,K,BFAX,NHL,NIL,NFAX,NU
      REAL(KIND=JPRB) :: ANGLE,DEL
      INTEGER(KIND=JPIM) :: JFAX(10),NLFAX(7)
!
!     SUBROUTINE 'SET99B' - COMPUTES FACTORS OF N & TRIGONOMETRIC
!     FUNCTIONS REQUIRED BY FFT992.
!     BASED ON SET99, SET99B ALSO RETURNS VIA LUSEFFT992 WHETHER
!     FACTORS HAVE BEEN FOUND THAT CAN PERMIT (OR NOT) FFT992 TO BE USED.
!
      SAVE NLFAX
!
      DATA NLFAX/6,8,5,4,3,2,1/
!
      IXXX=1
!
      DEL=4.0E0_JPRB*ASIN(1.0E0_JPRB)/REAL(N,KIND=JPRB)
      NIL=0
      NHL=(N/2)-1
      DO 10 K=NIL,NHL
      ANGLE=REAL(K,KIND=JPRB)*DEL
      TRIGS(2*K+1)=COS(ANGLE)
      TRIGS(2*K+2)=SIN(ANGLE)
   10 CONTINUE
!
!     FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
!     LOOK FOR SIXES FIRST, STORE FACTORS IN DESCENDING ORDER
      NU=N
      IFAC=6
      K=0
      IL=1
   20 CONTINUE
      IF (MOD(NU,IFAC).NE.0) GO TO 30
      K=K+1
      JFAX(K)=IFAC
      IF (IFAC.NE.8) GO TO 25
      IF (K.EQ.1) GO TO 25
      JFAX(1)=8
      JFAX(K)=6
   25 CONTINUE
      NU=NU/IFAC
      IF (NU.EQ.1) GO TO 50
      IF (IFAC.NE.8) GO TO 20
   30 CONTINUE
      IL=IL+1
      IFAC=NLFAX(IL)
      IF (IFAC.GT.1) GO TO 20
!
      LDUSEFFT992=.FALSE.
      RETURN
!
!     NOW REVERSE ORDER OF FACTORS
   50 CONTINUE
      NFAX=K
      IFAX(1)=NFAX
      DO 60 I=1,NFAX
      IFAX(NFAX+2-I)=JFAX(I)
   60 CONTINUE
      IFAX(10)=N
      LDUSEFFT992=.TRUE.
      END SUBROUTINE SET99B
