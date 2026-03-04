! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LWV &
 & ( YDELWRAD,YDEPHLI,YDPHNC,KIDIA, KFDIA, KLON , KLEV , KUAER , KTRAER,&
 & PABCU, PB   , PBINT, PBSUR, PBTOP , PDBSL,&
 & PEMIS, PEMIW,&
 & PGA  , PGB  , PGASUR,PGBSUR,PGATOP, PGBTOP,&
 & PCNTRB,PFLUC,&
! for adjoint computation
 & PADJD, PADJU, PDBDT, PDISD, PDISU, PDWFSU, & 
 & PCOND  , PTTPA  , PTTPB , PTTA , PTTB  ,&
 & PZZA   , PZZB   , PXNA  , PXNB , PXDIVA,&
 & PXDIVB , PZZA1  , PZZB1 , PXNA1, PXNB1 ,&
 & PXDIVA1, PXDIVB1, PTTPA1, PTTPB1 &
 & )  

!**** *LWV*   - LONGWAVE RADIATION, VERTICAL INTEGRATION

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
!           FLUXES OR RADIANCES

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
! PB     : (KLON,NSIL,KLEV+1); SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
! PBINT  : (KLON,KLEV+1)     ; HALF-LEVEL PLANCK FUNCTIONS
! PBSUR  : (KLON,NSIL)       ; SURFACE SPECTRAL PLANCK FUNCTION
! PBTOP  : (KLON,NSIL)       ; T.O.A. SPECTRAL PLANCK FUNCTION
! PDBSL  : (KLON,KLEV*2)     ; SUB-LAYER PLANCK FUNCTION GRADIENT
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PGA, PGB                   ; PADE APPROXIMANTS
! PGASUR, PGBSUR             ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP             ; T.O.A. PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PFLUC(KLON,2,KLEV)           ; RADIATIVE FLUXES CLEAR-SKY
!     ==== OUTPUTS FOR ADJOINT COMPUTATION ===
! PADJ.. : (KLON,KLEV+1)       ; CONTRIBUTION OF ADJACENT LAYERS
! PDBDT  : (KLON,NUA,KLEV)     ; LAYER PLANCK FUNCTION GRADIENT
! PDWFSU : (KLON,NSIL)         ; SPECTRAL DOWNWARD FLUX AT SURFACE
! PDIS.. : (KLON,KLEV+1)       ; CONTRIBUTION BY DISTANT LAYERS

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
!     CONTRIBUTIONS BY -  THE NEARBY LAYERS
!                      -  THE DISTANT LAYERS
!                      -  THE BOUNDARY TERMS
!          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.

!     EXTERNALS.
!     ----------

!          *LWVN*, *LWVD*, *LWVB*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        JJ Morcrette 96-06-07 Surface LW window emissivity
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M. Janiskova  22-Nov-2006 call for reduced LW routines 
!                                  (H2O and CO2 only)
!        M. Janiskova   4-Apr-2007 additional arguments passed from
!                                  the routine for AD computation
!        M. Janiskova  19-May-2008 using shorter NPROMALW instead of KLON
!        M. Janiskova  02-Mar-2012 initialization for PCNTRB
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOELW    , ONLY : NSIL     ,NIPD     ,NUA
USE YOEPHLI  , ONLY : TEPHLI
USE YOPHNC   , ONLY : TPHNC
USE YOELWRAD , ONLY : TELWRAD

IMPLICIT NONE

TYPE(TELWRAD)     ,INTENT(IN)    :: YDELWRAD
TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
TYPE(TPHNC)       ,INTENT(IN)    :: YDPHNC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KUAER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KLON,NSIL,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBINT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBSUR(KLON,NSIL) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBTOP(KLON,NSIL) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBSL(KLON,NSIL,KLEV*2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGA(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGB(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGASUR(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGBSUR(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGATOP(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGBTOP(KLON,NIPD,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCNTRB(KLON,KLEV+1,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUC(KLON,2,KLEV+1) 
! for adjoint computation
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PADJD(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PADJU(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDBDT(KLON,NSIL,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISD(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISU(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDWFSU(KLON,NSIL)
! for adjoint computation - from lwvdr
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOND(YDELWRAD%NPROMALW,KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPA(YDELWRAD%NPROMALW,KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPB(YDELWRAD%NPROMALW,KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZA1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZB1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNA1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNB1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVA1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVB1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPA1(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPB1(KLON,KLEV+1)

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZDWFSU(KLON,NSIL)  

INTEGER(KIND=JPIM) :: JA, JK, JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "lwvb.intfb.h"
#include "lwvbr.intfb.h"
#include "lwvd.intfb.h"
#include "lwvdr.intfb.h"
#include "lwvn.intfb.h"
#include "lwvnr.intfb.h"

!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

IF (LHOOK) CALL DR_HOOK('LWV',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMALW=>YDELWRAD%NPROMALW, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, &
 & LH2OCO2=>YDPHNC%LH2OCO2, LWLOPT=>YDPHNC%LWLOPT)
DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    PADJD(JL,JK)=0.0_JPRB
    PADJU(JL,JK)=0.0_JPRB
    PDISD(JL,JK)=0.0_JPRB
    PDISU(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JA=1,NSIL
  DO JL=KIDIA,KFDIA
    ZDWFSU(JL,JA)=0.0_JPRB
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PCNTRB(JL,KLEV+1,KLEV+1)=0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------

!     ------------------------------------------------------------------

!*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
!                  ---------------------------------

IF (LPHYLIN .AND. LH2OCO2) THEN
  CALL  LWVNR &
   & ( KIDIA, KFDIA, KLON  , KLEV , KUAER,&
   & PABCU, PDBSL, PGA   , PGB,&
   & PADJD, PADJU, PCNTRB, PDBDT, ZDWFSU  &
   & )
ELSE
  CALL LWVN &
   & ( KIDIA, KFDIA, KLON  , KLEV , KUAER,&
   & PABCU, PDBSL, PGA   , PGB,&
   & PADJD, PADJU, PCNTRB, PDBDT, ZDWFSU  &
   & )
ENDIF

!     ------------------------------------------------------------------

!*         2.2     CONTRIBUTION FROM DISTANT LAYERS
!                  ---------------------------------

IF (LPHYLIN .AND. LH2OCO2) THEN
  CALL LWVDR &
   & ( YDELWRAD, KIDIA , KFDIA, KLON , KLEV  , KTRAER,&
   & PABCU , PDBDT, PGA  , PGB,&
   & PCNTRB, PDISD, PDISU, ZDWFSU, &
! for adjoint computation
   & PCOND  , PTTPA  , PTTPB , PTTA , PTTB  ,&
   & PZZA   , PZZB   , PXNA  , PXNB , PXDIVA,&
   & PXDIVB , PZZA1  , PZZB1 , PXNA1, PXNB1 ,&
   & PXDIVA1, PXDIVB1, PTTPA1, PTTPB1 &
   & )
ELSE
  CALL LWVD &
   & ( YDPHNC, KIDIA , KFDIA, KLON , KLEV  , KTRAER,&
   & PABCU , PDBDT, PGA  , PGB,&
   & PCNTRB, PDISD, PDISU, ZDWFSU,&
! for adjoint computation
   &  PTTA , PTTB &
   & )
ENDIF

IF (.NOT. LWLOPT) THEN
  DO JA=1,NSIL
    DO JL=KIDIA,KFDIA
      PDWFSU(JL,JA) = ZDWFSU(JL,JA)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*         2.3     EXCHANGE WITH THE BOUNDARIES
!                  ----------------------------

IF (LPHYLIN .AND. LH2OCO2) THEN
  CALL LWVBR &
   & ( KIDIA , KFDIA , KLON  , KLEV  , KUAER,&
   & PABCU , PADJD , PADJU,&
   & PB    , PBINT , PBSUR , PBTOP,&
   & PDISD , PDISU , PEMIS , PEMIW,&
   & PGASUR, PGBSUR, PGATOP, PGBTOP,&
   & ZDWFSU,PFLUC  &
   & )
ELSE
  CALL LWVB &
   & ( KIDIA , KFDIA , KLON  , KLEV  , KUAER,&
   & PABCU , PADJD , PADJU,&
   & PB    , PBINT , PBSUR , PBTOP,&
   & PDISD , PDISU , PEMIS , PEMIW,&
   & PGASUR, PGBSUR, PGATOP, PGBTOP,&
   & ZDWFSU,PFLUC  &
   & )
ENDIF

IF (LWLOPT) THEN
  DO JA=1,NSIL
    DO JL=KIDIA,KFDIA
      PDWFSU(JL,JA) = ZDWFSU(JL,JA)
    ENDDO
  ENDDO
ENDIF

!-----------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LWV',1,ZHOOK_HANDLE)
END SUBROUTINE LWV
