! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE EGGRVS (PRPI, PRA, PDELX, PDELY, KPROF,&
 & KBEG, KEND, KULOUT, PGELAM, PGELAT, PGM, PGNORX, PGNORY)  
!****
!---------------------------------------------------------------------

!     GEOGRAPHY OF GRID-POINTS, INVERSION FROM GRID TO GEOGRAPHICAL SPHERE
!     ARPEGE-ALADIN
!     --------------------------------------------------------------------

!       ---------------------------------------------------
!     PURPOSE
!     -------
!      KNOWING THE PRECISE GEOGRAPHICAL TRANSFORMATION FROM
!      ARGUMENTS AND COMMON /YEMGGCM/, COMPUTES THE LOCATION
!      ON THEE GEOGRAPHICAL POINTS GIVEN IN INPUT BY THEIR LOCATION
!      ON THE ARPEGE-ALADIN GRID

!      MUST BE USED IN CONNECTION WITH SUBROUTINE EGGX
!      EITHER WITHIN IT OR AFTER IT

!     INPUT PARAMETERS
!     ----------------
!      PRPI : PI (3.14ETC)
!      PRA  : A, RADIUS OF PLANET
!      PDELX, PDELY : GRID SIZE IN M IF PROJECTION, OR IN RADIANS
!      KPROF : SIZE OF (1D) ARRAYS
!      KBEG, KEND : BEGINNING AND END POINTS OF CALCULATIONS
!      KULOUT : LOGICAL UNIT OF LISTING
!      PGELAM(KPROF) : X LOCATION OF POINTS, DISTANCE UNDER PROJECTION,
!                            RELATIVE ROTATED LONGITUDE UNDER ROTATION,
!                      DEFINED AS (JLON-KDLUN)*PDELX
!      PGELAT(KPROF) : Y LOCATION OF POINTS, DISTANCE UNDER PROJECTION,
!                            RELATIVE ROTATED LATITUDE UNDER ROTATION,
!                      DEFINED AS (JLAT-KDGUN)*PDELY
!            UNDER ROTATION, THE POSITION OF THE ORIGIN (XLAT1R,XLON1U)
!            IS HANDLED BY THIS SUBROUTINE : ONLY RELATIVE LOCATION
!            NEED TO BE SPECIFIED

!     IMPLICIT INPUT
!     --------------
!      COMMON /YEMGGCM/ MUST HAVE BEEN INITIALIZED

!     OUTPUT PARAMETERS
!     -----------------
!      PGELAM (KPROF): GEOGRAPHICAL LONGITUDE
!      PGELAT (KPROF): GEOGRAPHICAL LATITUDE
!      PGM    (KRPOF): MAP FACTOR
!      PGNORX (KPROF): PROJECTION OF GEOGRAPHICAL NORTH ON X
!      PGNORY (KPROF): PROJECTION OF GEOGRAPHICAL NORTH ON Y

!     WRITTEN BY
!     ---------- ALAIN JOLY

!      ORIGINAL NORTHERN HEMISPHERE VERSION : 27/2/92
!      SOUTHERN HEMISPHERE VERSION : 27/1/93

!     Modified:
!     --------

!            98-05-07: P. Le Moigne :vectorization of eggrvs (LLSTOP)

!-------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YEMGGCM  , ONLY : NYMGGI   ,NYMGGR   ,NYMGGWH  ,XLATR    ,&
 & XLONR    ,XGGPK    ,XRPKSM   ,XLAT0R   ,XLON0R   ,&
 & XLON0U   ,XIPORE   ,XJPORE   ,XGGM0    ,XLON1U   ,&
 & XLAT1R   ,HSUD     ,XBETA  

!-------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBEG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGELAM(KPROF) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGELAT(KPROF) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGM(KPROF) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGNORX(KPROF) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGNORY(KPROF) 

!-------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JJ

LOGICAL :: LLGWH, LLSTOP

REAL(KIND=JPRB) :: Z2PIPK, ZCOBETA, ZCOLA, ZCOSA, ZCOSOG, ZDIST,&
 & ZFUN, ZGAM, ZGM, ZKDL, ZLAT, ZLATG, ZLON, &
 & ZLONG, ZNORX, ZNORXP, ZNORY, ZNORYP, ZPIS2, &
 & ZPIS4, ZRPKSM2, ZSECAN, ZSECUR, ZSIBETA, &
 & ZSINAG, ZSINOG, ZURA2, ZX, ZY  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------

#include "abor1.intfb.h"

!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EGGRVS',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------

LLSTOP=.FALSE.
IF ( NYMGGI /= 10 ) THEN
  WRITE (KULOUT,*) '*** EGGRVS *** UNINITIALISED MODULE '
  CALL ABOR1(' EGGRVS: NYMGGI /= 10 ')
ENDIF

!     INITIALISE ROTATION ANGLE AND OTHER CONSTANTS
!     ---------------------------------------------
ZPIS2 = PRPI*0.5_JPRB
ZPIS4 = PRPI*0.25_JPRB
ZSECUR = 1.E-12_JPRB
ZSECAN = 1.E-05_JPRB
LLGWH = .FALSE.
IF ( NYMGGWH == 1 ) LLGWH = .TRUE.

DO JJ = KBEG, KEND
  PGM(JJ) = 1.0_JPRB
  PGNORX(JJ) = 0.0_JPRB
  PGNORY(JJ) = 1.0_JPRB
ENDDO

!     CORRECTION OF POSITION UNDER ROTATION ONLY
!     -------------------------------------------

IF ( XGGPK < 0.0_JPRB ) THEN
  DO JJ = KBEG, KEND
    PGELAM(JJ) = XLON1U + PGELAM(JJ)
    IF ( PGELAM(JJ) >= 2.0_JPRB*PRPI )PGELAM(JJ) = PGELAM(JJ) - 2.0_JPRB*PRPI
    PGELAT(JJ) = XLAT1R + PGELAT(JJ)
  ENDDO
ENDIF

!*
!--------------------------------------------------------------------
!     1.- REVERSE STEREO-LAMBERT PROJECTION
!     -------------------------------------
IF ( XGGPK > 0.0_JPRB ) THEN
  ZRPKSM2 = XRPKSM*XRPKSM
  ZURA2 = 1.0_JPRB/( PRA*PRA )
  Z2PIPK = 2.0_JPRB*PRPI*XGGPK

  DO JJ = KBEG, KEND
    ZX = PGELAM(JJ) - XIPORE*PDELX
    ZY = PGELAT(JJ) - XJPORE*PDELY
    ZDIST = (ZRPKSM2*( ZX*ZX + ZY*ZY )*ZURA2)**(1.0_JPRB/(2.0_JPRB*XGGPK))

    ZLAT = ZPIS2 - 2.0_JPRB*ATAN( ZDIST )

    IF ( ZDIST  < ZSECUR ) THEN
      ! THE POLE IS TOO NEAR TO DEFINE LONGITUDE
      ZLON = 0.0_JPRB
      IF ( XGGPK /= 1.0_JPRB ) LLSTOP=.TRUE.
      ! MAP FACTOR AT POLE FOR STEREOGRAPHIC PROJECTION
      ZGM = ( 1.0_JPRB + SIN( XLAT0R ) )/( 1.0_JPRB + SIN( ZLAT ) )
    ELSE
      IF ( HSUD > 0.0_JPRB ) THEN
        ZKDL = ATAN2( ZX,-ZY )
      ELSE
        ZKDL = ATAN2( ZX,ZY )
      ENDIF
      ZLON = XLON0R + (ZKDL + XBETA)/XGGPK
      IF ( ZLON < 0.0_JPRB ) ZLON = 2.0_JPRB*PRPI + ZLON
      ZGM = XGGM0*( COS( ZLAT )**(XGGPK-1.0_JPRB) )*&
       & ( ( 1.0_JPRB + SIN( ZLAT ) )**(-XGGPK) )  
    ENDIF

    PGELAM(JJ) = ZLON
    PGELAT(JJ) = HSUD*ZLAT
    PGM(JJ) = ZGM
    ! COMPONENTS OF VECTOR ROTATION MATRIX
    IF ( LLGWH ) ZLON = ZLON + 2.0_JPRB*PRPI
    IF ( XGGPK < 1.0_JPRB .AND. (ZLON-XLON0U)  >  Z2PIPK )&
     & ZLON = ZLON - 2.0_JPRB*PRPI  
    ZGAM = HSUD*(XGGPK*( ZLON - XLON0U ) - XBETA)
    PGNORX(JJ) = -SIN( ZGAM )
    PGNORY(JJ) = COS( ZGAM )
  ENDDO
  IF ( LLSTOP ) THEN
    WRITE (KULOUT,*) ' POLE WITHIN LAMBERT DOMAIN '
    CALL ABOR1(' EGGRVS: POLE WITHIN LAMBERT DOMAIN ')
  ENDIF
ENDIF

!*
!--------------------------------------------------------------------
!     2.- REVERSE MERCATOR PROJECTION
!     -------------------------------
IF ( XGGPK == 0.0_JPRB ) THEN
  ZSIBETA = SIN( XBETA )
  ZCOBETA = COS( XBETA )
  DO JJ = KBEG, KEND
    ZY = PGELAM(JJ)*ZSIBETA + PGELAT(JJ)*ZCOBETA - XJPORE*PDELY
    ZDIST = EXP( -ZY/( PRA*COS(XLAT0R) ) )
    ZLAT = ZPIS2 - 2.0_JPRB*ATAN( ZDIST )

    ZGM = COS( XLAT0R )/COS( ZLAT )

    ZX = PGELAM(JJ)*ZCOBETA - PGELAT(JJ)*ZSIBETA - XIPORE*PDELX
    ZLON = XLON0U + ZX/( PRA*COS( XLAT0R ) )
    IF ( ZLON >= 2.0_JPRB*PRPI ) ZLON = ZLON - 2.0_JPRB*PRPI

    PGELAM(JJ) = ZLON
    PGELAT(JJ) = ZLAT
    PGM(JJ) = ZGM
    ! COMPONENTS OF VECTOR ROTATION MATRIX
    PGNORX(JJ) = ZSIBETA
    PGNORY(JJ) = ZCOBETA
  ENDDO
ENDIF

!*
!--------------------------------------------------------------------
!     3.- REVERSE ROTATION
!     --------------------
IF ( NYMGGR /= 0 ) THEN
  DO JJ = KBEG, KEND
    ZLON = PGELAM(JJ)
    ZLAT = PGELAT(JJ)
    ZSINAG = SIN( XLATR )*COS( ZLAT )*COS( ZLON ) +COS( XLATR )*SIN( ZLAT )
    ZSINAG = MIN(1.0_JPRB,MAX(-1.0_JPRB,ZSINAG))
    ZLATG = ASIN( ZSINAG )
    IF ( ABS( ZLATG ) >= ZPIS2 ) THEN
      ZLONG = 0.0_JPRB
    ELSE
      ZCOSA = COS( ZLATG )
      ZFUN = COS( XLATR )*COS( ZLAT )*COS( ZLON ) -SIN( XLATR )*SIN( ZLAT )
      ZCOSOG = ( COS( XLONR )*ZFUN - SIN( XLONR )*COS( ZLAT )*&
       & SIN( ZLON ) )/ZCOSA  
      ZCOSOG = MIN(1.0_JPRB,MAX(-1.0_JPRB,ZCOSOG))
      ZSINOG = ( SIN( XLONR )*ZFUN + COS( XLONR )*COS( ZLAT )*&
       & SIN( ZLON ) )/ZCOSA  
      ZSINOG = MIN(1.0_JPRB,MAX(-1.0_JPRB,ZSINOG))
      ZLONG = ACOS( ZCOSOG )
      IF ( ASIN(ZSINOG) < 0.0_JPRB ) ZLONG = 2.0_JPRB*PRPI - ZLONG
    ENDIF
    PGELAM(JJ) = ZLONG
    PGELAT(JJ) = ZLATG
    ! COMPONENTS OF ROTATION MATRIX DUE TO ROTATION
    ZCOLA = SQRT( ABS( 1.0_JPRB - (COS( XLATR )*SIN( ZLATG )-&
     & SIN( XLATR )*COS( ZLATG )*COS( ZLONG-XLONR ))**2 ) )  
    IF ( ZCOLA < ZSECUR ) THEN
      WRITE (KULOUT,*) ' *** EGGX QUASI ERROR ***',&
       & ' DOMAIN EXTENDS UP TO NEW POLE : IT IS PROBABLY TOO LARGE'  
      ZNORY = 1.0_JPRB
      ZNORX = 0.0_JPRB
    ELSE
      ZNORY = ( COS( XLATR )*COS( ZLATG ) + SIN( XLATR )*&
       & SIN( ZLATG )*COS( ZLONG-XLONR ) )/ZCOLA  
      ZNORX = - SIN( XLATR )*SIN( ZLONG-XLONR )/ZCOLA
    ENDIF
    ! COMPOSITION OF THIS ROTATION WITH THE ONE RESULTING FROM PROJECTION
    ZNORYP = PGNORY(JJ)
    ZNORXP = PGNORX(JJ)
    PGNORY(JJ) = ZNORY*ZNORYP - ZNORX*ZNORXP
    PGNORX(JJ) = ZNORX*ZNORYP + ZNORY*ZNORXP
  ENDDO
ENDIF

!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EGGRVS',1,ZHOOK_HANDLE)
END SUBROUTINE EGGRVS
