! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WSMFEN (FSEA, EW, FM, USTT, DPT)

!---------------------------------------------------------------------- 

!**** *WSMFEN* - COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA      
!****            SPECTRUM.                                              

!     P.LIONELLO     ECMWF       APRIL 1990                             

!     PURPOSE.                                                          
!     --------                                                          

!       COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA SPECTRUM.     

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *WSMFEN (FSEA, EW, FM, USTT, DPT)*                            

!        *FSEA*   REAL     WINDSEA PART OF THE WAVE SPECTRUM.           
!        *EW*     REAL     WINDSEA ENERGY .                             
!        *FM*     REAL     MEAN FREQUENCY OF THE WINDSEA.               
!        *USTT*   REAL     FRICTION VELOCITY.                           
!        *DPT*    REAL     DEPTH.                           

!     METHOD.                                                           
!     -------                                                           

!       ENERGY AND MEAN FREQUENCY ARE FIRST COMPUTED BY INTEGRATION     
!       OF THE WINDSEA PART OF THE SPECTRUM:                            
!       THIS IMPLIES AN UNDERESTIMATION BOTH OF ENERGY AND MEAN         
!       FREQUENCY. THE UNDERESTIMATED VALUES ARE USED IN THE GROWTH     
!       CURVE PRODUCING OVERSTIMATES OF BOTH ENERGY AND MEAN FREQUENCY. 
!       THE AVERAGE OF THE TWO ESTIMATES IS TAKEN TO PROVIDE A BEST     
!       ESTIMATE.                                                       

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : AFCRV   ,BFCRV     ,                        &
     &            ESH      ,ASH     ,BSH
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH   ,DFIMFR    ,    &
     &            WETAIL   ,WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: USTT, DPT
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN) :: FSEA
      REAL(KIND=JWRB), INTENT(OUT) :: EW, FM

      INTEGER(KIND=JWIM) :: M, K
      REAL(KIND=JWRB) :: EN, YNU, EMAX, X
      REAL(KIND=JWRB) :: DELT25, COEF1, XTEMP, TEMP, SPFB
      REAL(KIND=JWRB) :: ESTAR, DSTAR, SPINTDI, FREQDI 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!     INLINE FUNCTIONS.
!     -----------------

!     ENERGY AS FUNCTION OF MEAN FREQUENCY.

      EN(X)=AFCRV*X**BFCRV

!     MEAN FREQUENCY AS FUNCTION OF ENERGY.

      YNU(X)=(X/AFCRV)**(1.0_JWRB/BFCRV)

!     DIMENSIONLESS ENERGY LIMIT AS FUNCTION OF DIMENSIONLESS DEPTH 

      EMAX(X) = ESH*TANH(ASH*X**BSH)

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WSMFEN',0,ZHOOK_HANDLE)

!*    1. INTEGRATING THE WINDSEA PART OF THE SPECTRUM.                  
!*       (UNDERESTIMATE IS IMPLIED)                                     
!        ---------------------------------------------                  


!*    1.1 COMPUTATION OF THE WIND-SEA ENERGY.                           
!*    1.2 COMPUTATION OF THE MEAN FREQUENCY.                            
!         -----------------------------------                           

      EW = EPSMIN
      FM = EPSMIN
      DO M=1,NFRE
        TEMP = 0.0_JWRB
        DO K=1,NANG
          TEMP = TEMP+FSEA(K,M)
        ENDDO
        EW = EW+TEMP*DFIM(M)
        FM = FM+TEMP*DFIMFR(M)
      ENDDO
      DELT25 = WETAIL*FR(NFRE)*DELTH
      EW = EW+DELT25*TEMP                                               
      COEF1 = WP1TAIL*DELTH*FR(NFRE)**2
      FM = FM+COEF1*TEMP
      FM = FM/EW

! ----------------------------------------------------------------------

!*    2. ESTIMATES FROM MODEL RELATIONS.                                
!        -------------------------------                                

!*    2.1 ENERGY IS DERIVED FROM UNDERESTIMATED MEAN FREQUENCY.         
!         -----------------------------------------------------         

      XTEMP = FM*USTT/G

      ESTAR=EN(FM*USTT/G)
      DSTAR = DPT*G/USTT**2
      ESTAR = MIN(EMAX(DSTAR),ESTAR)
      SPINTDI = USTT**4/G**2 * ESTAR

!*    2.2 MEAN FREQUENCY IS DERIVED FROM THE UNDERESTIMATED ENERGY.     
!         ---------------------------------------------------------     

      ESTAR=EW*G**2/USTT**4
      DSTAR = DPT*G/USTT**2
      ESTAR = MIN(EMAX(DSTAR),ESTAR)

      FREQDI = G/USTT * YNU(ESTAR)

! ----------------------------------------------------------------------

!*    3. FINAL ESTIMATE
!        ---------------

!*    3.1 AVERAGING THE ENERGY ESTIMATES
!         -------------------------------

      EW = (EW+SPINTDI)*0.5_JWRB

!*    3.2 AVERAGING THE MEAN FREQUENCY ESTIMATES.                       
!         ---------------------------------------                       

      FM = (FM+FREQDI)*0.5_JWRB
 
      IF (LHOOK) CALL DR_HOOK('WSMFEN',1,ZHOOK_HANDLE)

      END SUBROUTINE WSMFEN
