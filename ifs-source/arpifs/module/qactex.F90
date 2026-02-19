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

MODULE QACTEX
!--------------------------------------------------------------------
!  QACTEX : Contient les variables logiques controlant les
!  ------   differentes phases de CANARI
!           Logical variables controling the different CANARI steps.
!--------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!  * LAEINC  : .T. incremental (3 input files)
!  * LAEOMF  : .T. computation and taking of residuals through NCMOMF
!  * LAEOMN  : .T. computation and taking of residuals through NCMOMN
!  * LAECHK  : .T. observation quality control
!  * LAEPDS  : .T. surface pressure analysis
!  * LAEUVT  : .T. wind and temperature analysis
!  * LAEHUM  : .T. humidity analysis
!  * LAET2M  : .T. 2 meter temperature analysis
!  * LAEH2M  : .T. 2 meter humidity analysis
!  * LAEV1M  : .T. 10 meter wind analysis
!  * LAESNM  : .T. snow analysis (Cressman)
!  * LAESST  : .T. SST analysis
!  * LECSST  : .T. use ECMWF SST
!  * LAEICS  : .T. surface fields initialisation
!  * LAEICS_SX : .T. run inline surface assimilation inline
!  * LL_SODA : .T. Use SODA interface for surface assimilation
!  * LAECDS  : .T. surface diagnostic fields calculation
!  * LAESTU  : .T. use assimilated statistics
!  * LAESTA  : .T. statistics assimilation
!  * LAEWIO  : .T. write an ARPEGE file containing analysed fields
!  * LAERFO  : .T. copy file CMAFOC when ending analysis
!  * LAEUPFLG: .T. unplug the update of the observation use flag
!  * LVERAL  : .T. non-update of the observation quality flags
!  * NAEINC  :  0 --> OMF (default), 1 --> OMF + FC1, 2 --> OMN
!  * RCLIMCA : rappel coefficient towards surface fields climatology
!  * RCLISST : rappel coefficient towards SST climatology
!  * NSSTLIS : rappel towards US SST (and number of possible delay days)
!  * NSEAICE : use SSM/I ice limit (and number of possible delay days)
!  * RSNSA   : A-coefficient for snow analysis
!  * RSNSB   : B-coefficient for snow analysis
!  * RWPIA   : A-coefficient for ice analysis
!  * RWPIB   : B-coefficient for ice analysis

LOGICAL :: LAEINC
LOGICAL :: LAEOMF
LOGICAL :: LAEOMN
LOGICAL :: LAECHK
LOGICAL :: LAEPDS
LOGICAL :: LAEUVT
LOGICAL :: LAEHUM
LOGICAL :: LAET2M
LOGICAL :: LAEH2M
LOGICAL :: LAEV1M
LOGICAL :: LAESNM
LOGICAL :: LAESST
LOGICAL :: LECSST
LOGICAL :: LAEICS
LOGICAL :: LAECDS
LOGICAL :: LAESTU
LOGICAL :: LAESTA
LOGICAL :: LAEWIO
LOGICAL :: LAERFO
LOGICAL :: LAEUPFLG
LOGICAL :: LVERAL
LOGICAL :: LAEICS_SX
LOGICAL :: LL_SODA
INTEGER(KIND=JPIM) :: NAEINC
REAL(KIND=JPRB) :: RCLIMCA
REAL(KIND=JPRB) :: RCLISST
INTEGER(KIND=JPIM) :: NSSTLIS
INTEGER(KIND=JPIM) :: NSEAICE
REAL(KIND=JPRB) :: RSNSA
REAL(KIND=JPRB) :: RSNSB
REAL(KIND=JPRB) :: RWPIA
REAL(KIND=JPRB) :: RWPIB

!--------------------------------------------------------------------
END MODULE QACTEX
