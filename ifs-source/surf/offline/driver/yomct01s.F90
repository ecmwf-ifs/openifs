MODULE YOMCT01S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!*    ------------------------------------------------------------------

!**   Control variables for the job - constant within job

!========== TECHNICAL SWITCHES ======================

! LNF     : .T. = start, .F. = restart
! LMPLOT  : .T. = plotting requested
! NFRPLT  : plotting frequency
! NCYCLE  : number of the experiment
! CNMEXP  : name of the experiment
!           An experiment is identified by its name (16 characters)
!           and its cycle (typically same experiment but without a bug)

!========== MODEL     SWITCHES ======================

! NSTART  : first timestep of model
! NSTOP   : last timestep of model
! NFRPOS  : frequency of post-processing events (time-steps)
! NFRRES  ! frequency of writing intermediate restart files (time-steps)
! NFRHIS  : frequency of history write_ups (time-steps)
! NPOSTS  : array containing postprocessing steps
! NHISTS  : array containing history write-up steps

!     EXPLANATION :
!     1) IF NXXXTS(0)=0 ACTION IF MOD(JSTEP,NFRXXX)=0
!     2) IF NXXXTS(0)>0 NXXXTS(0) SIGNIFICANT NUMBERS IN
! NXXXTS ARE THEN CONSIDERED AND :
!       ACTION FOR JSTEP=NXXXTS(.)*NFRXXX
!     3) IF NXXXTS(0)<0
!       ACTION FOR JSTEP= (NXXXTS(.)/DELTAT)*NFRXXX

!========== ECMWF Single Column Model =========================================
! LSCMEC  : .T. = ECMWF Single Column Model
! LSFCFLX : .T. = forcing with surface fluxes (latent and sensible).
! REXTSHF : externally supplied sensible heat flux [W/m^2]
! REXTLHF : externally supplied latent   heat flux [W/m^2]
! LROUGH  : .T. = surface roughness length is externally specified
! REXTZ0M : externally supplied roughness length for momentum [m]
! REXTZ0H : externally supplied roughness length for heat [m]



INTEGER(KIND=JPIM), PARAMETER :: JPNPST=240
INTEGER(KIND=JPIM) :: NPOSTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NHISTS(0:JPNPST)

CHARACTER*16 CNMEXP

INTEGER(KIND=JPIM) :: NFRPLT
INTEGER(KIND=JPIM) :: NCYCLE
INTEGER(KIND=JPIM) :: NSTART
INTEGER(KIND=JPIM) :: NSTOP
INTEGER(KIND=JPIM) :: NFRPOS
INTEGER(KIND=JPIM) :: NFRRES
INTEGER(KIND=JPIM) :: NFRHIS
LOGICAL LNF
LOGICAL LMPLOT
! * ECMWF Single Column Model:
LOGICAL :: LSCMEC
LOGICAL :: LSFCFLX
REAL(KIND=JPRB) :: REXTSHF
REAL(KIND=JPRB) :: REXTLHF
LOGICAL :: LROUGH
REAL(KIND=JPRB) :: REXTZ0M
REAL(KIND=JPRB) :: REXTZ0H

!     ------------------------------------------------------------------
END MODULE YOMCT01S
