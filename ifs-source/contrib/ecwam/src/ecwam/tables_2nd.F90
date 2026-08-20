! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!--------------------------------------------------------------------
!
      SUBROUTINE TABLES_2ND(NFRE,NANG,NDEPTH,OMSTART,FRAC,XMR,          &
     &                      DFDTH,OMEGA,TH,TA,TB,TC_QL,                 &
     &                      TT_4M,TT_4P,IM_P,IM_M,TFAK)
!
!--------------------------------------------------------------------
!
!*****TABLES** COMPUTES TABLES FOR SECOND ORDER SPECTRUM IN FREQUENCY SPACE.
!
!     P.JANSSEN DECEMBER 2008
!
!     PURPOSE
!     -------
!             DETERMINES TABLES, BASED ON JANSSEN (2008)
!             THERE ARE THREE CORRECTIONS:
!                   1) GENERATION OF SECOND-HARMONICS
!                   2) QUASI-LINEAR EFFECT
!                   3) SHIFT OF SPECTRUM BECAUSE OF STOKES FREQUENCY
!                      CORRECTION.
!
!     INTERFACE
!     ---------
!             *CALL* *TABLES(NFRE,NANG,NDEPTH,OMSTART,FRAC,XMR,
!                            OMEGA,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M,
!                            TFAK)*
!
!
!     PARAMETER   TYPE      PURPOSE.
!     ---------   ----      -------
!
!       NFRE      INTEGER   NUMBER OF FREQUENCIES
!       NANG      INTEGER   NUMBER OF DIRECTIONS
!       NDEPTH    INTEGER   NUMBER OF ENTRIES IN THE DEPTH TABLE 
!       OMSTART   REAL      START FREQUENCY
!       FRAC      REAL      FRACTIONAL INCREASE IN FREQUENCY SPACE
!       XMR       REAL      INVERSE OF THINNING FACTOR IN FREQUENCY SPACE
!       DFDTH     REAL      FREQ. DIRECTION INCREMENT
!       OMEGA     REAL      ANGULAR FREQUENCY ARRAY
!       TH        REAL      DIRECTION ARRAY
!       TA        REAL      TABLE FOR MINUS INTERACTIONS
!       TB        REAL      TABLE FOR PLUS INTERACTIONS
!       TC_QL     REAL      TABLE FOR QUASI-LINEAR INTERACTIONS
!       TT_4M     REAL      TABLE FOR STOKES FREQUENCY CORRECTION
!       TT_4P     REAL      TABLE FOR STOKES FREQUENCY CORRECTION
!       IM_P      INTEGER   TABLE FOR WAVENUMBER M2 PLUS
!       IM_M      INTEGER   TABLE FOR WAVENUMBER M2 MIN
!       TFAK      REAL      WAVENUMBER TABLE
!
!
!     METHOD
!     ------
!             
!     EXTERNALS
!     ---------
!             NONE
!
!     REFERENCES
!     ----------
!             V.E. ZAKHAROV, HAMILTONIAN APPROACH (1968) 
!             M.A. SROKOSZ, J.G.R.,91,995-1006 (1986)
!             P.A.E.M. JANSSEN, ECMWF TECH MEMO (2008),JFM PAPER (2009)
!
!
!--------------------------------------------------------------------
 
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWSHAL, ONLY: DEPTHA, DEPTHD
      USE YOWPCONS , ONLY : G, ZPI
      USE YOWCONST_2ND, ONLY: DPTH 

!--------------------------------------------------------------------

      IMPLICIT NONE
include "aki.intfb.h"
   
      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE, NANG, NDEPTH
      INTEGER(KIND=JWIM), DIMENSION(NFRE,NFRE), INTENT(OUT) :: IM_P, IM_M

      REAL(KIND=JWRB), INTENT(IN) :: OMSTART, FRAC, XMR
      REAL(KIND=JWRB), DIMENSION(NFRE), INTENT(IN) :: DFDTH, OMEGA
      REAL(KIND=JWRB), DIMENSION(NANG), INTENT(IN) :: TH
      REAL(KIND=JWRB), DIMENSION(NDEPTH,NANG,NFRE,NFRE), INTENT(OUT) :: TA,TB,TC_QL,TT_4M,TT_4P
      REAL(KIND=JWRB), DIMENSION(NFRE,NDEPTH), INTENT(OUT) :: TFAK


      INTEGER(KIND=JWIM) :: JD,M,K,M1,K1,MP,MM,L

      REAL(KIND=JWRB) :: OM0,TH0,XK0,OM1,TH1,XK1,OM2,XK2,OM0P,XK0P,OM0M
      REAL(KIND=JWRB) :: XK0M,XM2,A,B,C_QL,FAC,W2,V2 

!
!     1. COMPUTATION OF WAVENUMBER ARRAY TFAK
!     ---------------------------------------
!
!
      DO JD=1,NDEPTH
         DPTH = DEPTHA*DEPTHD**(JD-1)
         DO M=1,NFRE
            OM0 = OMEGA(M)
            TFAK(M,JD) = AKI(OM0,DPTH)
         ENDDO
!
!     2. COMPUTATION OF THE 2nd ORDER COEFFICIENTS.
!     ---------------------------------------------
!
!
         K1 = 0
         TH1 = TH(NANG)
         DO M=1,NFRE
            OM0 = OMEGA(M)
            XK0 = TFAK(M,JD)
          
            MP   = MIN(M+1,NFRE)
            OM0P = OMEGA(MP)
            XK0P = TFAK(MP,JD) 

            MM   = MAX(M-1,1) 
            OM0M = OMEGA(MM)
            XK0M = TFAK(MM,JD) 
                  
            DO M1=1,NFRE
                        
               OM1 = OMEGA(M1)
           
               DO L=1,NANG
!
!              XK0-XK1 CASE 
!           
                  K = K1+L
                  TH0 = TH(K)
                  OM2 = OM0-OM1
         
             
                  IF (ABS(OM1) < OM0/2.0_JWRB) THEN
                     XM2  = LOG(OM2/OMSTART)/LOG(1.0_JWRB+FRAC)
                     IM_M(M1,M) = NINT(XMR*(XM2+1.0_JWRB))
                     XK1 = TFAK(M1,JD)
                     XK2 = AKI(OM2,DPTH) 

                     TA(JD,L,M1,M) = DFDTH(M1)*A(XK1,XK2,TH1,TH0)**2
                  ELSE
                     TA(JD,L,M1,M) = 0.0_JWRB
                     IM_M(M1,M) = 1
                  ENDIF 
!
!              XK1+XK0 CASE 
!           
                  OM2 = OM1+OM0 
                  XM2  = LOG(OM2/OMSTART)/LOG(1.0_JWRB+FRAC)
                  IM_P(M1,M) = NINT(XMR*(XM2+1.0_JWRB))
                  XK1 = TFAK(M1,JD)
                  XK2 = AKI(OM2,DPTH)

                  TB(JD,L,M1,M) = DFDTH(M1)*B(XK1,XK2,TH1,TH0)**2
!
!              QUASI-LINEAR EFFECT
!           
!  
                  TC_QL(JD,L,M1,M) = DFDTH(M1)*C_QL(XK0,XK1,TH0,TH1)
!  
!              STOKES-FREQUENCY CORRECTION 
!
!      
                  FAC = 2.0_JWRB*G/OM1*DFDTH(M1)
                  TT_4M(JD,L,M1,M) =                                    &
     &                  FAC*(W2(XK0M,XK1,XK1,XK0M,TH0,TH1,TH1,TH0)+     &
     &                       V2(XK0M,XK1,XK1,XK0M,TH0,TH1,TH1,TH0))
                  TT_4P(JD,L,M1,M) =                                    &
     &                  FAC*(W2(XK0P,XK1,XK1,XK0P,TH0,TH1,TH1,TH0)+     &
     &                       V2(XK0P,XK1,XK1,XK0P,TH0,TH1,TH1,TH0))

               ENDDO
            ENDDO
         ENDDO
       ENDDO
! 
      END SUBROUTINE TABLES_2ND
