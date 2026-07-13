MODULE SPPGUST_MOD
CONTAINS
SUBROUTINE SPPGUST(KIDIA, KFDIA, KLON &
 & , PZ0MM, PBUOM, PUSTAR, PU10M, PV10M &
 & , YDCST, YDEXC &
 ! OUTPUTS
 & , PGUST)  

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_EXCS , ONLY : RCHBCD, RCHBBCD, RCHBB, RCHBD, RCHBA, RCHBHDL, &
 & RCDHALF, RCHETB, RCHB23A, RCHETA, RCDHPI2
USE YOS_CST  , ONLY : TCST
USE YOS_EXC  , ONLY : TEXC

! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *SPPGUST* - COMPUTES THE area averaged 10 m wind and the gust

!     Author.
!     -------
!     A. Beljaars       E.C.M.W.F.    24/02/2000
!     A. Beljaars       E.C.M.W.F.    15/11/2001 Fix orography problem

!     Modifications.
!     --------------
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     P. Viterbo  ECMWF  12/05/2005 Externalize SURF (based on vdfppgust)
!     A. Beljaars ECMWF  18/02/2006 Revised gust to accomodate stochastic physics
!     N.Semane+P.Bechtold 04-10-2012 Add RPARZI

!     PURPOSE
!     -------

!     Compute wind gusts

!     INTERFACE
!     ---------

!     *SPPGUST* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PZ0MM*        AERODYNAMIC ROUGHNESS LENGTH
!     *PBUOM*        BUOYANCY FLUX
!     *PUSTAR*       FRICTION VELOCITY
!      PU10M         U-COMPONENT WIND AT 10 M                         m/s
!      PV10M         V-COMPONENT WIND AT 10 M                         m/s

!     OUTPUT PARAMETERS (REAL):

!     *PGUST*        WIND GUST AT 10 M

!     METHOD
!     ------

!     MO scaling is used to estimate turbulence intensity               

!     MODIFICATIONS
!     -------------

!     18/12/2020 P. Bechtold+A. Beljaars Change gust cefficient ZUGN from 7.72 to 7.2

!     --------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUOM(:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSTAR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU10M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV10M(:) 
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(:) 

!*    LOCAL STORAGE
!     ----- -------

REAL(KIND=JPRB) ::    ZFKLEV(KLON)

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) ::  Z10M, Z10MP, &
 & ZDL, &
 & ZL, ZNLEV, &
 & ZUGN, ZCZI, ZUSTAR,&
 & Z1D3, ZIPBL, ZIDL, ZFF10, ZOROC  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcsvdfs.h"

!     ------------------------------------------------------------------

!*       1.   INITIALIZE CONSTANTS
!             ---------- ----------

!     THIS SPECIFIES THE HEIGHT FOR U,V (10 M) 

IF (LHOOK) CALL DR_HOOK('SPPGUST_MOD:SPPGUST',0,ZHOOK_HANDLE)
ASSOCIATE(RG=>YDCST%RG, RLSTT=>YDCST%RLSTT, RLVTT=>YDCST%RLVTT, RTT=>YDCST%RTT, &
 & REPDU2=>YDEXC%REPDU2, RKAP=>YDEXC%RKAP, RPARZI=>YDEXC%RPARZI)
Z10M=10._JPRB
ZOROC=5._JPRB

!     For the gust model, a dimensionless number is 
!     used ZUG=(UMAX-U)/u*. It is computed from the gust model 
!     with Kaimal spectrum in (see Beljaars 1987, J. Atmos. Oc. Techn., 4, 
!     613-626). The number is reasonably constant for most conditions, 
!     but dependes on Zi/L in the same way as the standard deviation 
!     of horizontal wind (Panosky et al.1977). ZUG depends also on 
!     the assumed probabilty of exceedence P in the following way for 
!     zi/L large: 


!     In view of increased resolution, coefficient ZUGN was adjusted in Cy48r1 to obtain a better
!     match to observations. The new value is ZUGN=7.20.
!     This value is also obtained from the gust model with a time series length of 2200 s.
!     The choice of 2200 s corresponds to a horizontal length scale of 22 km at 10 m/s, following
!     the Taylor frozen turbulence assumption. 
!     for a mean wind of 10 m/s, and an anemometer with 3 seconds averaging:
!
!     P          0.10    0.25    0.50    0.75    0.90
!     -----------------------------------------------
!     ZUGN       8.44    7.81    7.20    6.69    6.29 

!     The stability function is ZUG=ZUGN*( 1+ (0.5/12.)zi/L )^(1./3.)

!     Recomputation with a integration time of 6000 s instead of 600 s
!     The choice of 10 min is for gusts with respect to 10 min averages. 
!     The choice of 100 min is for an area of 60 km at 10 m/s (more model compatible). 

!     P    !  0.10  0.25  0.50  0.75  0.90 
!     ------------------------------------  For an anemometer with 3 s averaging
!     ZUGN !        8.23  7.71  7.26 
!             These numbers apply at 20 m/s

ZUGN=7.20_JPRB
ZCZI=0.5_JPRB/12._JPRB
Z1D3=1._JPRB/3._JPRB
ZIPBL=RPARZI


!     ------------------------------------------------------------------

!        2.   COMPUTE HORIZONTAL WIND AND GUTST
!             ---------------------------------

DO JL=KIDIA,KFDIA

!     AREA AVERAGE OF ABSOLUTE 10 M (TO BE USED FOR GUSTS)

  ZUSTAR=PUSTAR(JL)
  ZIDL=-ZIPBL*RKAP*PBUOM(JL)/ZUSTAR**3
  ZFF10=SQRT(PU10M(JL)**2+PV10M(JL)**2)

  IF (ZIDL >= 0.0_JPRB) THEN
    PGUST(JL)=ZFF10+ZUSTAR*ZUGN
  ELSE
    PGUST(JL)=ZFF10+ZUSTAR*ZUGN*(1.0_JPRB-ZCZI*ZIDL)**Z1D3
  ENDIF
  
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPPGUST_MOD:SPPGUST',1,ZHOOK_HANDLE)
END SUBROUTINE SPPGUST
END MODULE SPPGUST_MOD
