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

MODULE YOMSCC

USE PARKIND1 , ONLY : JPIM, JPRB
USE PARCMA   , ONLY : JPMXOCT
USE PARDIMO  , ONLY : JPXAMVPROD
USE RTTOV_CONST, ONLY: NINST

IMPLICIT NONE

SAVE

!*     YOMSCC - OBSERVATION SCREENING CONTROL PARAMETERS
!                     SETUP IN SUBROUTINE*SUOBS/NAMELIST*NAMSCC

!        HEIKKI JARVINEN   ECMWF     15/3/96
!        ELISABETH GERARD  ECMWF     15/01/1999 (Modification)
!            SSM/I wind speed and total column water vapour
!            have their own thinning parameter.
!        DAVID TAN           ECMWF     25/11/2003 (Modification)
!            Variables for Analysis Ensembles
!        CHRIS PAYAN       MF        ??/04/2005
!            addition of NICESEAMASK integer key
!        HANS HERSBACH       ECMWF     15/11/2007
!            introduce NSCAWSOLMAX: max number of scat ambiguities
!            introduce LQSCATAZI
!        CHRIS PAYAN       MF        ??/04/2007
!            addition of NSCAT5_RQC integer key
!        PAUL POLI         MF        06/03/2007
!            addition of thinning parameters for GPSRO
!        CHRIS PAYAN       MF          /04/2010
!            36t2 ice-sea QC now handled in blacklist, NICE_SEAMASK no longer used
!        CHRIS PAYAN       MF          /10/2012
!            base 38t1_bf, LOSCAT_OUTERSWATH_REJ
!        P CHAMBON LF MEUNIER MF       24/04/2014
!            BC_QCSAPHIR: add fixed bias correction coefficient for SAPHIR quality control
!        BRUCE INGLEBY     ECMWF     25/04/2014
!            addition of LBTEMDUP, RDUPD_BTEM/SURF, RDUPT_BTEM and RFACV_BTEM
!        CRISTINA LUPU    ECMWF     19/11/2014
!            Updates to sensors and platform lists from rttov_consts
!        GIOVANNA DE CHIARA  ECMWF     01/02/2015
!            Updates for RapiSCAT and HY-2A Scat sensors
!        CHRIS PAYAN       MF         11/12/2014
!            base 40_op2, LHSCAT_OUTERSWATH_REJ, LRSCAT_OUTERSWATH_REJ added (H/RSCAT)
!        CHRIS PAYAN       MF         13/09/2017
!            base 43t2, LSSCAT_OUTERSWATH_REJ added for ScatSAT-1
!        CHRIS PAYAN       MF         11/07/2017
!            base 43t2, updates for ScatSat-1 Scat sensor
!        GIOVANNA DE CHIARA           04/2018
!            base 45r1, LFSCAT_OUTERSWATH_REJ and LWSCAT_OUTERSWATH_REJ added for CFOSAT and WindRAD
!                       HY-2A changed to HY-2B
!        Bruce Ingleby  2019-03-28  Add LAIRVSWITCH
!        Marcin Chrust  2022-06-20  Add REDNMC

!     NAME      TYPE                  MEANING
!     ----      ----                  -------

!     LSCDPR      L    SCREENING DIAGNOSTIC PRINTOUT (DEFAULT = .FALSE.)
!                      - CONTROLS THE AMOUNT OF PRINTOUT IN SCREENING RUN
!     LSCRE4D     L    OBSERVATION SCREENING IN 4D-MODE (DEFAULT .F. = 3D-MODE)
!     LSRERUN     L    RE-RUNNING IFS WITH A SECOND-HAND CMA-FILE (F)
!     LTARGDROP   L    RELAX BACKGROUND CHECK FOR TARGETED DROPSONDES
!     LAIRVSWITCH L    Change sign of v-wind if it seems wrong - B787 aircraft

!     REDNMC     - Global scaling factor for background error standard deviations
!     RMIND_TOV  - Minimum distance ( degrees ) between TOVS obs
!     RFIND_TOV  - Average distance ( degrees ) between TOVS obs
!     RMIND_SSMI - Minimum distance ( degrees ) between SSMI obs
!     RMIND_SATOB - Minimum distance ( degrees ) between SATOB obs
!     RFIND_SATOB - Average distance ( degrees ) between SATOB obs
!     RMIND_RAD1C  - Minimum distance (degrees) between RAD1C obs
!     RFIND_RAD1C  - Average distance (degrees) between RAD1C obs
!     RFIND_SSMI_PWC - Average distance (degrees) between SSMI obs
!                      for total column water vapour
!     RFIND_SSMI_WSP - Average distance (degrees) between SSMI obs
!                      for surface wind speed
!     RMIND_SATAM - Minimum distance (degrees) between SATAM obs
!     RFIND_SATAM - Average distance (degrees) between SATAM obs
!     RFIND_AIREP - Average distance (degrees) between AIREP obs
!     RMIND_SCATT - Minimum distance (degrees) between SCAT obs
!     RFIND_SCATT - Average distance (degrees) between SCAT obs
!     RMIND_GPSRO - Minimum distance (degrees) between GPSRO obs
!     RFIND_GPSRO - Average distance (degrees) between GPSRO obs
!     RMIND_RADAR - Minimum distance (degrees) between RADAR obs
!     RFIND_RADAR - Average distance (degrees) between RADAR obs
!     RDUPD_SURF  - Maximum distance (degrees) in surface duplicate check
!     RDUPD_BTEM  - Maximum distance (degrees) in BTEM duplicate check
!     RDUPT_BTEM  - Maximum time diff (minutes) in BTEM duplicate check
!     RDUPD_DROP  - Maximum distance (degrees) in BTEM DROP duplicate check
!     RDUPT_DROP  - Maximum time diff (minutes) in BTEM DROP duplicate check
!     RFACV_BTEM  - Vertical thinning factor for BTEM
!     NBTEMDUP    - if 0, then no special BTEM duplicate check
!                 - if 1, then BTEM/TEMP duplicate check
!                 - if 2, then BTEM/TEMPorBTEM duplicate check + BTEM thinning
!     NTHINSCA    I    Thinning parameter for ERS scatterometer
!     NSCAWSOLMAX - Maximum allowed number or ambigious scatterometer wind solutions
!     LQSCATAZI   - if T, perform QC on azimuth diversity for QuikSCAT
!     LOSCATAZI   - if T, perform QC on azimuth diversity for OCEANSAT2/3
!     LHSCATAZI   - if T, perform QC on azimuth diversity for HY-2A/B/C/D Scatterometer
!     LRSCATAZI   - if T, perform QC on azimuth diversity for RapidSCAT
!     LSSCATAZI   - if T, perform QC on azimuth diversity for ScatSAT
!     LFSCATAZI   - if T, perform QC on azimuth diversity for RFSCAT
!     LWSCATAZI   - if T, perform QC on azimuth diversity for WindRAD
!     LOSCAT_OUTERSWATH_REJ - if T, OSCAT wvc 1,2,3,4 and 33,34,35,36 rejected
!     LHSCAT_OUTERSWATH_REJ - if T, HSCAT wvc 1,2,3,4 and 35,36,37,38 rejected
!     LRSCAT_OUTERSWATH_REJ - if T, RSCAT wvc 1,2 and 20,21 rejected
!     LSSCAT_OUTERSWATH_REJ - if T, SSCAT wvc 1,2,3,4 and 35,36,37,38 rejected
!     LFSCAT_OUTERSWATH_REJ - if T, FSCAT wvc 1,2,3,4 and 35,36,37,38 rejected !GDC for testing purposes
!     LWSCAT_OUTERSWATH_REJ - if T, WSCAT wvc 1,2,3,4 and 35,36,37,38 rejected !GDC for testing purposes
!     NSCAT5_RQC  - /1,2/ Quikscat RainQC based on /raincontamination,dist2cone/
!     LDFS - if T, perturbation of all observations (active + passive) in order to perturbate also radar reflecivity data for AROME analysis
!     FL_MAXSPEED - Maximum speed expected for aircrafts in km/hour
!                   (used to remove unrealistic AIREP reports)
!     SELAIREPPRE - Pressure levels used for aircraft thinning
!     RAIREPTHIN  - Vertical thinning distance in pressure [Pa] used for aircraft thinning (default 1500.)
!     RAIREPPCENTTHIN  - Minimum difference of pressure ratio between levels used for aircraft thinning (default 5%==0.05)
!     RAIREPTOPPRES - Top pressure level [Pa] for aircraft thinning (default 10000.)
!     NSATAM_THINPRODLST - list of producers which are thinned for SATAM (AMV)
!     NSATAM_CHOSENQIBYPROD - QI choice function of producer
!                /1,2,3/ = /QI1(MPEF FG-QI),QI2(RFF),QI3(MPEF noFG-QI)/
!     RGPSROTHIN - Vertical thinning distance [m] used for GPSRO thinning

!     LPERTURB          L  TO PERTURB IN ANALYSIS ENSEMBLE (DEFAULT = .FALSE.)
!     NAENSEMBLE        I  ANALYSIS ENSEMBLE GROUP NUMBER
!     NAEMEMBER         I  ANALYSIS ENSEMBLE MEMBER NUMBER
!     LPERT_SATOB_CORR  L  SPATIAL CORRELATIONS IN SATOB PERTURBATIONS (DEFAULT=.T. IF.NOT.LELAM)

!     RGPSROTHIN - Vertical thinning distance [m] used for GPSRO thinning
!     BC_QCSAPHIR - fixed bias correction coefficient for SAPHIR/Megha-Tropiques MW sounder quality control

LOGICAL :: LSCDPR
LOGICAL :: LSCRE4D
LOGICAL :: LSRERUN
LOGICAL :: LTARGDROP
LOGICAL :: LPERTURB
LOGICAL :: LPERT_SATOB_CORR
LOGICAL :: LQSCATAZI
LOGICAL :: LOSCATAZI
LOGICAL :: LHSCATAZI
LOGICAL :: LRSCATAZI
LOGICAL :: LSSCATAZI
LOGICAL :: LFSCATAZI
LOGICAL :: LWSCATAZI
LOGICAL :: LAIRVSWITCH
LOGICAL :: LOSCAT_OUTERSWATH_REJ
LOGICAL :: LHSCAT_OUTERSWATH_REJ
LOGICAL :: LRSCAT_OUTERSWATH_REJ
LOGICAL :: LSSCAT_OUTERSWATH_REJ
LOGICAL :: LFSCAT_OUTERSWATH_REJ
LOGICAL :: LWSCAT_OUTERSWATH_REJ
LOGICAL :: LDFS

REAL(KIND=JPRB) :: REDNMC
REAL(KIND=JPRB) :: RMIND_TOV
REAL(KIND=JPRB) :: RFIND_TOV
REAL(KIND=JPRB) :: RMIND_RAD1C(0:NINST-1)
REAL(KIND=JPRB) :: RFIND_RAD1C(0:NINST-1)
REAL(KIND=JPRB) :: RMIND_SSMI
REAL(KIND=JPRB) :: RFIND_SSMI_PWC
REAL(KIND=JPRB) :: RFIND_SSMI_WSP
REAL(KIND=JPRB) :: RMIND_SATOB
REAL(KIND=JPRB) :: RFIND_SATOB
REAL(KIND=JPRB) :: RMIND_SATAM
REAL(KIND=JPRB) :: RFIND_SATAM
REAL(KIND=JPRB) :: RFIND_AIREP
REAL(KIND=JPRB) :: RMIND_SCATT
REAL(KIND=JPRB) :: RFIND_SCATT
REAL(KIND=JPRB) :: RMIND_GPSRO
REAL(KIND=JPRB) :: RFIND_GPSRO
REAL(KIND=JPRB) :: RMIND_RADAR
REAL(KIND=JPRB) :: RFIND_RADAR
REAL(KIND=JPRB) :: RDUPD_SURF
REAL(KIND=JPRB) :: RDUPD_BTEM
REAL(KIND=JPRB) :: RDUPT_BTEM
REAL(KIND=JPRB) :: RDUPD_DROP
REAL(KIND=JPRB) :: RDUPT_DROP
REAL(KIND=JPRB) :: RFACV_BTEM
REAL(KIND=JPRB) :: FL_MAXSPEED
REAL(KIND=JPRB), ALLOCATABLE :: SELAIREPPRE(:)
REAL(KIND=JPRB) :: RAIREPTHIN
REAL(KIND=JPRB) :: RAIREPPCENTTHIN
REAL(KIND=JPRB) :: RAIREPTOPPRES
REAL(KIND=JPRB) :: RGPSROTHIN
REAL(KIND=JPRB) :: BC_QCSAPHIR(1:6)
REAL(KIND=JPRB) :: BC_QCMTVZAGY(1:29)
REAL(KIND=JPRB) :: BC_QCAMSR2(1:14)

INTEGER(KIND=JPIM) :: NTHINSCA
INTEGER(KIND=JPIM) :: NSCAWSOLMAX(JPMXOCT)
INTEGER(KIND=JPIM) :: NSCAT5_RQC
INTEGER(KIND=JPIM) :: NAENSEMBLE, NAEMEMBER
INTEGER(KIND=JPIM) :: NSATAM_THINPRODLST(JPXAMVPROD)
INTEGER(KIND=JPIM) :: NSATAM_CHOSENQIBYPROD(JPXAMVPROD)
INTEGER(KIND=JPIM) :: NBTEMDUP

!-----------------------------------------------------------------------

END MODULE YOMSCC
