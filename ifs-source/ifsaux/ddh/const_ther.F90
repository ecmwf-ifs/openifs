MODULE CONST_THER
! --------------------------------------------------------------
! **** const_ther physiques.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Auteur/author:   2001-01, J.M. Piriou d'après ARPEGE.
! Modifications:
!      F. Vana  05-Mar-2015  Support for single precision
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
SAVE
!
!-------------------------------------------------
! Nombre d'itérations de la boucle de Newton.
!-------------------------------------------------
!
INTEGER(KIND=JPIM), PARAMETER :: NBITER=2
!      -----------------------------------------------------------------
!
!*       1.    DEFINE FUNDAMENTAL CONSTANTS.
!              -----------------------------
!
REAL(KIND=JPRB), PARAMETER :: RPI=3.14159265358979 ! pi.
REAL(KIND=JPRB), PARAMETER :: RCLUM=299792458. ! célérité de la lumière.
REAL(KIND=JPRB), PARAMETER :: RHPLA=6.6260755E-34 ! cte de Planck.
REAL(KIND=JPRB), PARAMETER :: RKBOL=1.380658E-23 ! cte de Bolzman.
REAL(KIND=JPRB), PARAMETER :: RNAVO=6.0221367E+23 ! nombre d'Avogadro.
!
!     ------------------------------------------------------------------
!
!*       2.    DEFINE ASTRONOMICAL CONSTANTS.
!              ------------------------------
!
REAL(KIND=JPRB), PARAMETER :: RDAY=86400. ! jour solaire.
REAL(KIND=JPRB), PARAMETER :: REA=149597870000. ! demi-grand axe de rév. terrestre.
REAL(KIND=JPRB), PARAMETER :: RSIYEA=365.25*RDAY*2.*RPI/6.283076 ! année sidérale.
REAL(KIND=JPRB), PARAMETER :: RSIDAY=RDAY/(1._JPRB+RDAY/RSIYEA) ! jour sidéral.
REAL(KIND=JPRB), PARAMETER :: ROMEGA=2.*RPI/RSIDAY ! vitesse angulaire terrestre.
!
! ------------------------------------------------------------------
!
! *       3.    DEFINE GEOIDE.
! --------------
!
REAL(KIND=JPRB), PARAMETER :: RG=9.80665 ! accélération de la pesanteur.
REAL(KIND=JPRB), PARAMETER :: RA=6371229. ! rayon terrestre.
!
! ------------------------------------------------------------------
!
! *       4.    DEFINE RADIATION CONSTANTS.
! ---------------------------
!
! Stefan-Bolzman const : 2. * rpi**5 * rkbol**4 /(15.* rclum**2 * rhpla**3) 
REAL(KIND=JPRB), PARAMETER :: RSIGMA=2. * RPI**5 * (RKBOL/RHPLA)**3 * RKBOL/(15.* RCLUM**2)
!real, parameter :: rsigma=5.670509E-08
REAL(KIND=JPRB), PARAMETER :: RI0=1370. ! cte solaire.
!
! ------------------------------------------------------------------
!
! *       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
! ------------------------------------------
!
REAL(KIND=JPRB), PARAMETER :: R=RNAVO*RKBOL ! cte des gaz parfaits.
REAL(KIND=JPRB), PARAMETER :: RMD=28.9644
REAL(KIND=JPRB), PARAMETER :: RMV=18.0153
REAL(KIND=JPRB), PARAMETER :: RMO3=47.9942
REAL(KIND=JPRB), PARAMETER :: RD=1000.*R/RMD ! cte spécifique de l'air sec.
REAL(KIND=JPRB), PARAMETER :: RV=1000.*R/RMV ! cte spécifique de la vapeur d'eau.
REAL(KIND=JPRB), PARAMETER :: RCPD=3.5*RD ! chaleur massique de l'air sec.
REAL(KIND=JPRB), PARAMETER :: RCVD=RCPD-RD
REAL(KIND=JPRB), PARAMETER :: RCPV=4. *RV ! chaleur massique de la vapeur d'eau.
REAL(KIND=JPRB), PARAMETER :: RCVV=RCPV-RV
REAL(KIND=JPRB), PARAMETER :: RKAPPA=RD/RCPD
REAL(KIND=JPRB), PARAMETER :: RETV=RV/RD-1._JPRB
!
REAL(KIND=JPRB), PARAMETER :: RALPW =  .6022274788E+02
REAL(KIND=JPRB), PARAMETER :: RBETW =  .6822400210E+04
REAL(KIND=JPRB), PARAMETER :: RGAMW =  .5139266694E+01
!
REAL(KIND=JPRB), PARAMETER :: RALPS =  .3262117981E+02
REAL(KIND=JPRB), PARAMETER :: RBETS =  .6295421339E+04
REAL(KIND=JPRB), PARAMETER :: RGAMS =  .5631331575E+00
!
REAL(KIND=JPRB), PARAMETER :: RALPD = -.2760156808E+02
REAL(KIND=JPRB), PARAMETER :: RBETD = -.5269788712E+03
REAL(KIND=JPRB), PARAMETER :: RGAMD = -.4576133537E+01
! ------------------------------------------------------------------
!
! *       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
! ---------------------------------------------
!
REAL(KIND=JPRB), PARAMETER :: RCW=4218. ! chaleur massique de l'eau liquide.
!
! ------------------------------------------------------------------
!
! *       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
! --------------------------------------------
!
REAL(KIND=JPRB), PARAMETER :: RCS=2106. ! chaleur massique de la glace.
!
! ------------------------------------------------------------------
!
! *       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
! ----------------------------------------------------
!
REAL(KIND=JPRB), PARAMETER :: RTT=273.16 ! point triple de l'eau.
REAL(KIND=JPRB), PARAMETER :: RDT=11.82
REAL(KIND=JPRB), PARAMETER :: RLVTT=2.5008E+6 ! chaleur latente eau vapeur > eau liquide.
REAL(KIND=JPRB), PARAMETER :: RLSTT=2.8345E+6 ! chaleur latente eau vapeur > eau glace.
REAL(KIND=JPRB), PARAMETER :: RLVZER=RLVTT+RTT*(RCW-RCPV) ! chaleur latente de fusion à 0°K!
REAL(KIND=JPRB), PARAMETER :: RLSZER=RLSTT+RTT*(RCS-RCPV) ! chaleur latente de sublimation à 0°K!
REAL(KIND=JPRB), PARAMETER :: RLMLT=RLSTT-RLVTT ! chaleur latente eau liquide > eau glace.
REAL(KIND=JPRB), PARAMETER :: RATM=100000. ! pression standard.
!
! ------------------------------------------------------------------
!
! *       9.    SATURATED VAPOUR PRESSURE.
! --------------------------
!
REAL(KIND=JPRB), PARAMETER :: RESTT=611.14
!
!-------------------------------------------------
! Constante de Joule.
!-------------------------------------------------
!
REAL(KIND=JPRB), PARAMETER :: RJOULE=4.184

END MODULE CONST_THER
