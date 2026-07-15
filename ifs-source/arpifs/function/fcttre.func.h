! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!*
!     ------------------------------------------------------------------

!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!       Two sets of functions are available. In the first set only the
!       cases water or ice are distinguished by temperature.  This set 
!       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
!       The second set considers, besides the two cases water and ice 
!       also a mix of both for the temperature range YDTHF%RTICE < T < YDTHF%RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
!       FKOOP modifies the ice saturation mixing ratio for homogeneous 
!       nucleation. FOE_DEWM_DT provides an approximate first derivative
!       of FOEEWM.

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.

!     ------------------------------------------------------------------
!     *****************************************************************

!                NO CONSIDERATION OF MIXED PHASES

!     *****************************************************************
REAL(KIND=JPRB) :: FOEDELTA
REAL(KIND=JPRB) :: PTARE
FOEDELTA (PTARE) = MAX (0.0_JPRB,SIGN(1.0_JPRB,PTARE-YDCST%RTT))

!                  FOEDELTA = 1    water
!                  FOEDELTA = 0    ice

!     THERMODYNAMICAL FUNCTIONS .

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEW,FOEDE,FOEDESU,FOELH,FOELDCP
FOEEW ( PTARE ) = YDTHF%R2ES*EXP (&
  &(YDTHF%R3LES*FOEDELTA(PTARE)+YDTHF%R3IES*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE-YDCST%RTT)&
&/ (PTARE-(YDTHF%R4LES*FOEDELTA(PTARE)+YDTHF%R4IES*(1.0_JPRB-FOEDELTA(PTARE)))))

FOEDE ( PTARE ) = &
  &(FOEDELTA(PTARE)*YDTHF%R5ALVCP+(1.0_JPRB-FOEDELTA(PTARE))*YDTHF%R5ALSCP)&
&/ (PTARE-(YDTHF%R4LES*FOEDELTA(PTARE)+YDTHF%R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2

FOEDESU ( PTARE ) = &
  &(FOEDELTA(PTARE)*YDTHF%R5LES+(1.0_JPRB-FOEDELTA(PTARE))*YDTHF%R5IES)&
&/ (PTARE-(YDTHF%R4LES*FOEDELTA(PTARE)+YDTHF%R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2

FOELH ( PTARE ) =&
         &FOEDELTA(PTARE)*YDCST%RLVTT + (1.0_JPRB-FOEDELTA(PTARE))*YDCST%RLSTT

FOELDCP ( PTARE ) = &
         &FOEDELTA(PTARE)*YDTHF%RALVDCP + (1.0_JPRB-FOEDELTA(PTARE))*YDTHF%RALSDCP

!     *****************************************************************

!           CONSIDERATION OF MIXED PHASES

!     *****************************************************************

!     FOEALFA is calculated to distinguish the three cases:

!                       FOEALFA=1            water phase
!                       FOEALFA=0            ice phase
!                       0 < FOEALFA < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEALFA
FOEALFA (PTARE) = MIN(1.0_JPRB,((MAX(YDTHF%RTICE,MIN(YDTHF%RTWAT,PTARE))-YDTHF%RTICE)&
 &*YDTHF%RTWAT_RTICE_R)**2) 


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEWM,FOEDEM,FOELDCPM,FOELHM,FOE_DEWM_DT
FOEEWM ( PTARE ) = YDTHF%R2ES *&
     &(FOEALFA(PTARE)*EXP(YDTHF%R3LES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4LES))+&
  &(1.0_JPRB-FOEALFA(PTARE))*EXP(YDTHF%R3IES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4IES)))

FOE_DEWM_DT( PTARE ) = YDTHF%R2ES * ( &
     & YDTHF%R3LES*FOEALFA(PTARE)*EXP(YDTHF%R3LES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4LES)) &
     &    *(YDCST%RTT-YDTHF%R4LES)/(PTARE-YDTHF%R4LES)**2 + &
     & YDTHF%R3IES*(1.0-FOEALFA(PTARE))*EXP(YDTHF%R3IES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4IES)) &
     &    *(YDCST%RTT-YDTHF%R4IES)/(PTARE-YDTHF%R4IES)**2)

FOEDEM ( PTARE ) = FOEALFA(PTARE)*YDTHF%R5ALVCP*(1.0_JPRB/(PTARE-YDTHF%R4LES)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))*YDTHF%R5ALSCP*(1.0_JPRB/(PTARE-YDTHF%R4IES)**2)

FOELDCPM ( PTARE ) = FOEALFA(PTARE)*YDTHF%RALVDCP+&
            &(1.0_JPRB-FOEALFA(PTARE))*YDTHF%RALSDCP

FOELHM ( PTARE ) =&
         &FOEALFA(PTARE)*YDCST%RLVTT+(1.0_JPRB-FOEALFA(PTARE))*YDCST%RLSTT


!     Temperature normalization for humidity background change of variable
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOETB
FOETB ( PTARE )=FOEALFA(PTARE)*YDTHF%R3LES*(YDCST%RTT-YDTHF%R4LES)*(1.0_JPRB/(PTARE-YDTHF%R4LES)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))*YDTHF%R3IES*(YDCST%RTT-YDTHF%R4IES)*(1.0_JPRB/(PTARE-YDTHF%R4IES)**2)

!     ------------------------------------------------------------------
!     *****************************************************************

!           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV

!     *****************************************************************

!     FOEALFCU is calculated to distinguish the three cases:

!                       FOEALFCU=1            water phase
!                       FOEALFCU=0            ice phase
!                       0 < FOEALFCU < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEALFCU 
FOEALFCU (PTARE) = MIN(1.0_JPRB,((MAX(YDTHF%RTICECU,MIN(YDTHF%RTWAT,PTARE))&
&-YDTHF%RTICECU)*YDTHF%RTWAT_RTICECU_R)**2) 


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEWMCU,FOEDEMCU,FOELDCPMCU,FOELHMCU
FOEEWMCU ( PTARE ) = YDTHF%R2ES *&
     &(FOEALFCU(PTARE)*EXP(YDTHF%R3LES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4LES))+&
  &(1.0_JPRB-FOEALFCU(PTARE))*EXP(YDTHF%R3IES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4IES)))

FOEDEMCU ( PTARE )=FOEALFCU(PTARE)*YDTHF%R5ALVCP*(1.0_JPRB/(PTARE-YDTHF%R4LES)**2)+&
             &(1.0_JPRB-FOEALFCU(PTARE))*YDTHF%R5ALSCP*(1.0_JPRB/(PTARE-YDTHF%R4IES)**2)

FOELDCPMCU ( PTARE ) = FOEALFCU(PTARE)*YDTHF%RALVDCP+&
            &(1.0_JPRB-FOEALFCU(PTARE))*YDTHF%RALSDCP

FOELHMCU ( PTARE ) =&
         &FOEALFCU(PTARE)*YDCST%RLVTT+(1.0_JPRB-FOEALFCU(PTARE))*YDCST%RLSTT
!     ------------------------------------------------------------------

!     Pressure of water vapour at saturation
!     This one is for the WMO definition of saturation, i.e. always
!     with respect to water.
!     
!     Duplicate to FOEELIQ and FOEEICE for separate ice variable
!     FOEELIQ always respect to water 
!     FOEEICE always respect to ice 
!     FOEELIQ2ICE is the ratio of vapour pressures over water and ice 
!       (analytically simplified to avoid a suspected Cray compiler bug and using a single call to exp)
!     (could use FOEEW and FOEEWMO, but naming convention unclear)
!     FOELSON returns e wrt liquid water using D Sonntag (1994, Met. Zeit.)
!      - now recommended for use with radiosonde data (WMO CIMO guide, 2014)
!      unlike the FOEE functions does not include 1/(YDCST%RETV+1.0_JPRB) factor

REAL(KIND=JPRB) :: FOEEWMO, FOEELIQ, FOEEICE, FOEELIQ2ICE, FOELSON 
FOEEWMO( PTARE ) = YDTHF%R2ES*EXP(YDTHF%R3LES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4LES))
FOEELIQ( PTARE ) = YDTHF%R2ES*EXP(YDTHF%R3LES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4LES))
FOEEICE( PTARE ) = YDTHF%R2ES*EXP(YDTHF%R3IES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4IES))
FOEELIQ2ICE( PTARE )  = EXP(      (PTARE-YDCST%RTT)*(YDTHF%R3LES/(PTARE-YDTHF%R4LES) - YDTHF%R3IES/(PTARE-YDTHF%R4IES)))
FOELSON( PTARE ) = EXP( -6096.9385_JPRB/PTARE + 21.2409642_JPRB &
	             - 2.711193E-2_JPRB * PTARE    &
                     + 1.673952E-5_JPRB * PTARE**2 &
		     + 2.433502_JPRB * LOG(PTARE))

REAL(KIND=JPRB) :: FOEEWM_V,FOEEWMCU_V,FOELES_V,FOEIES_V
REAL(KIND=JPRB) :: EXP1,EXP2
      FOELES_V(PTARE)=YDTHF%R3LES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4LES)
      FOEIES_V(PTARE)=YDTHF%R3IES*(PTARE-YDCST%RTT)/(PTARE-YDTHF%R4IES)
      FOEEWM_V( PTARE,EXP1,EXP2 )=YDTHF%R2ES*(FOEALFA(PTARE)*EXP1+ &
          & (1.0_JPRB-FOEALFA(PTARE))*EXP2)
      FOEEWMCU_V ( PTARE,EXP1,EXP2 ) = YDTHF%R2ES*(FOEALFCU(PTARE)*EXP1+&
          &(1.0_JPRB-FOEALFCU(PTARE))*EXP2)

