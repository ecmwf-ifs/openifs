! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

REAL(KIND=JWRB) FUNCTION AKI_ICE (G,XK,DEPTH,RHOW,CITH)

! ----------------------------------------------------------------------

!**** *AKI_ICE* - FUNCTION TO COMPUTE WAVE NUMBER UNDER THE ICE.
!                 FOR A GIVEN OPEN WATER WAVE NUMBER.

!*    PURPOSE.
!     -------

!       *AKI_ICE* COMPUTES THE REAL WAVE NUMBER UNDER THE ICE AS FUNCTION OF
!                 OPEN OCEAN WAVE NUMBER, THE WATER DEPTH AND THE ICE THICKNESS.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *AKI_ICE (G,XK,DEPTH,RHOW,CITH))*
!          *G*        - ACCELERATION OF GRAVITY (m/s**2).
!          *XK*       - OPEN OCEAN WAVE NUMBER (1/m).
!          *DEPTH*    - WATER DEPTH (m).
!          *RHOW*     - WATER DENSITY (kg/m**3).
!          *CITH*     - SEA ICE THICKNESS (m).

!     METHOD.
!     -------

!       NEWTONS METHOD TO SOLVE THE LINEAR DISPERSION RELATION IN 
!       SHALLOW WATER UNDER AN INFINITELY LONG ELASTIC FLOATING PLATE
!       REPRESENTING THE SEA ICE. THE MECHANICAL PROPERTIES OF THE SEA
!       ICE IS GIVEn BY ITS YOUNG MODULUS, THE POISSON'S RATIO AND ITS
!       DENSITY (these are fixed for now, see below). 

!       IF F(x)=0, then solve iteratively
!       x(n+1) = x(n) - F(x(n))/F'(x(n))

!       WHERE F'(x) IS THE FIRST DERIVATIVE OF F WITH RESPECT TO x.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOX AND SQUIRE, 1991, JGR 96, C3, 4531-4547.

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: G, XK, DEPTH, RHOW, CITH

!     ICE PROPERTIES (assumed fixed for now) 
      REAL(KIND=JWRB), PARAMETER :: YMICE=5.5E+9_JWRB   ! typical value of Young modulus of sea ice
      REAL(KIND=JWRB), PARAMETER :: RMUICE=0.3_JWRB      ! Poisson's ratio of sea ice
      REAL(KIND=JWRB), PARAMETER :: RHOI=922.5_JWRB     ! typical value of the sea ice density

!     RELATIVE ERROR LIMIT OF NEWTONS METHOD.
      REAL(KIND=JWRB), PARAMETER :: EBS = 0.000001_JWRB
!     MAXIMUM WAVE NUMBER
      REAL(KIND=JWRB), PARAMETER :: AKI_MAX = 20.0_JWRB

      REAL(KIND=JWRB) :: FICSTF, RDH
      REAL(KIND=JWRB) :: OM2, AKI, AKIOLD, F, FPRIME, AKID 


      IF (CITH <= 0.0_JWRB) THEN
        AKI=XK
      ELSE
!       BENDING STIFFNESS / WATER DENSITY
        FICSTF=(YMICE*CITH**3/(12*(1-RMUICE**2)))/RHOW

!       DENSITY RATIO * ICE THICKNESS
        RDH=(RHOI/RHOW)*CITH

!       SQUARE OF THE OPEN OCEAN ANGULAR FREQUENCY
        OM2=G*XK*TANH(XK*DEPTH)

!*      2. ITERATION LOOP.
!          ---------------

        AKIOLD=0.0_JWRB
        AKI=MIN(XK ,(OM2/MAX(FICSTF,1.0_JWRB))**0.2_JWRB)

        DO WHILE (ABS(AKI-AKIOLD) > EBS*AKIOLD .AND.                   &
     &            AKI < AKI_MAX                                        &
     &           )
          AKIOLD=AKI 
          AKID=MIN(DEPTH*AKI,50.0_JWRB)
          F=FICSTF*AKI**5+G*AKI-OM2*(RDH*AKI+1./TANH(AKID))
          FPRIME=5._JWRB*FICSTF*AKI**4+G-OM2*(RDH-DEPTH/(SINH(AKID)**2))
          AKI = AKI-F/FPRIME
!         in case of overshoot because it is trying to find a very large wave number
          IF (AKI <= 0.0_JWRB) AKI=AKI_MAX
        ENDDO

      ENDIF

      AKI_ICE = AKI

END FUNCTION AKI_ICE
