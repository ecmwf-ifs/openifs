! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
!       also a mix of both for the temperature range RTICE < T < RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.

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
REAL(KIND=JPRB) :: FOEEW, FOEDESU
FOEEW ( PTARE ) = R2ES*EXP (&
  &(R3LES*FOEDELTA(PTARE)+R3IES*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE-YDCST%RTT)&
&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE)))))

FOEDESU ( PTARE ) = &
  &(FOEDELTA(PTARE)*R5LES+(1.0_JPRB-FOEDELTA(PTARE))*R5IES)&
&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2
