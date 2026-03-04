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

MODULE YOMTAG

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

! Modifications :
! P.Marguinaud : 11-Sep-2012 : Add tags for MF IOs
! P.Marguinaud : 10-Oct-2013 : Add/remove tags for MF IOs

!     ------------------------------------------------------------------
!*    Tag identifiers used in message passing communication

! MTAGLM: tag for transpositions done in TRLTOM.
! MTAGMS: tag for transpositions done in TRMTOS.
! MTAGSM: tag for transpositions done in TRSTOM.
! MTAGSPNO: tag for communications done in COMMSPNORM and COMMSPNORM1.
! MTAGSLAG: tag for halo constitution (horizontal interpolations
!  in the semi-Lagrangian scheme, the observation interpolator or FULLPOS).
! MTAGRAD : tag for communications done in SUECRAD (ECMWF physics)
! MTAGPART: tag for communications done in DICOMOUT and GATHFLNM.
! MTAGDISTSP: tag for communications done in DISSPEC, DISSPEC0 and DIWRSPE.
! MTAGDISTGP: tag for communications done in
!  DISGRID, DISGRID_C, DISGRIDFP, DIWRGRFP, DIWRGRID, IRCVGPF, ISNDGPF,
!  ORCVGPF, OSNDGPF (in FULL-POS, DM-distribution in arrival geometry).
! MTAGDISTGP_DEP: tag for communications done in
!  DISGRIDFP, DIWRGRFP (in FULL-POS, DM-distribution in departure geometry).
! MTAGCAIN: tag for communications done in GATHERSPA.
! MTAGCOST: tag for communications done in
!  GATHERCOST1, GATHERCOST2, GATHERCOSTO and GATHERJCVERT.
! MTAGGSUM: tag for communications done in CASND1, CASNDR1 and GATHERSUM.
! MTAGGLOBSI: tag for communications done in CAEXCO and CAUPDO.
! MTAGGOMMP: tag for communications done in GOM_MESSAGE_PASSING.
! MTAGGOMMPAD: tag for communications done in GOM_MESSAGE_PASSING_AD.
! MTAGFCE: tag for communications done in
!  COMMFCE1, COMMFCE2, COMMJBBAL and COMMJBDAT.
! MTAGSIG: tag for communications done in SIGCHECK.
! MTAGBRPR: tag for communications done in BRPTOB and GATHERT.
! MTAGGPNORM: tag for communications done in GPNORM1.
! MTAGSPNORM: tag for communications done in COMMSPNORM.
! MTAGWAMNORM: tag for communications done in IFSTOWAM.
! MTAGDDHRES: tag for communications done in DDHRCV and DDHSND.
! MTAGDDH1: tag for communications done in DISTDDH.
! MTAGDDH2: tag for communications done in DLADDH.
! MTAGDDH3: tag for communications done in DMADDH.
! MTAGDDH4: tag for communications done in DRESDDH.
! MTAGGETV: tag for communications done in SUHESS.
! MT_DISTRIBUTED_VECTOR: tag for communications done in SUMPINI.
! MTAGLCZ: tag for communications done in COMMNSEC1.
! MTAGFREQ: tag for communications done in GATHERFREQ.
! MTAGEIGMD: tag for communications done in GATHEREIGMD.
! MTAGKE: tag for communications done in VMODEENERGY.

!      YOMTAG

!  Tag identifier range 1-1999 are reserved for diagnostics

INTEGER(KIND=JPIM), PARAMETER :: MTAGLM                =   600
INTEGER(KIND=JPIM), PARAMETER :: MTAGMS                =  1000
INTEGER(KIND=JPIM), PARAMETER :: MTAGSM                =  1200
INTEGER(KIND=JPIM), PARAMETER :: MTAGSPNO              =  2200
INTEGER(KIND=JPIM), PARAMETER :: MTAGSLAG              =  2400
INTEGER(KIND=JPIM), PARAMETER :: MTAGRAD               =  2800
INTEGER(KIND=JPIM), PARAMETER :: MTAGPART              =  3600
INTEGER(KIND=JPIM), PARAMETER :: MTAGCAIN              =  4200
INTEGER(KIND=JPIM), PARAMETER :: MTAGCOST              =  4400
INTEGER(KIND=JPIM), PARAMETER :: MTAGGSUM              =  4600
INTEGER(KIND=JPIM), PARAMETER :: MTAGFCE               =  5000
INTEGER(KIND=JPIM), PARAMETER :: MTAGPE                =  5200
INTEGER(KIND=JPIM), PARAMETER :: MTAGGPNORM            =  5400
INTEGER(KIND=JPIM), PARAMETER :: MTAGBRPR              =  5600
INTEGER(KIND=JPIM), PARAMETER :: MTAGLCZ               =  6000
INTEGER(KIND=JPIM), PARAMETER :: MTAGFREQ              =  6400
INTEGER(KIND=JPIM), PARAMETER :: MTAGEIGMD             =  6700
INTEGER(KIND=JPIM), PARAMETER :: MTAGKE                =  6800
INTEGER(KIND=JPIM), PARAMETER :: MTAGSPNORM            =  7000
INTEGER(KIND=JPIM), PARAMETER :: MTAGWAMNORM           =  7200
INTEGER(KIND=JPIM), PARAMETER :: MTAGDISTSP            = 10000
INTEGER(KIND=JPIM), PARAMETER :: MTAGDISTGP            = 11000
INTEGER(KIND=JPIM), PARAMETER :: MTAGDISTGP_DEP        = 12000
INTEGER(KIND=JPIM), PARAMETER :: MTAGGLOBSI            = 22001
INTEGER(KIND=JPIM), PARAMETER :: MTAGGOMMP             = 23000
INTEGER(KIND=JPIM), PARAMETER :: MTAGGOMMPAD           = 23001
INTEGER(KIND=JPIM), PARAMETER :: MTAGDDHRES            = 24000
INTEGER(KIND=JPIM), PARAMETER :: MTAGSIG               = 26000
INTEGER(KIND=JPIM), PARAMETER :: MTAGDDH1              = 27000
INTEGER(KIND=JPIM), PARAMETER :: MTAGDDH2              = 28000
INTEGER(KIND=JPIM), PARAMETER :: MTAGDDH3              = 29000
INTEGER(KIND=JPIM), PARAMETER :: MTAGDDH4              = 30000
INTEGER(KIND=JPIM), PARAMETER :: MTAGGETV              = 30200
INTEGER(KIND=JPIM), PARAMETER :: MT_DISTRIBUTED_VECTOR = 32000


! Tags for MF IO; suffix is subroutine name
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRFU        = 40000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRGRIDA     = 42000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRGRIDALL   = 43000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRGRIDUA    = 44000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRSFX       = 45000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRSPECA_GP  = 46000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRXFU       = 47000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRSPECA     = 48000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_FP2SX1FA    = 49000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_FP2SX2      = 50000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRHFP       = 51000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRSFP       = 52000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_SUSPECA     = 52000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_SUGRIDUA    = 53000
INTEGER(KIND=JPIM), PARAMETER :: MTAG_MFIO_WRSPECSPA   = 54000

END MODULE YOMTAG

