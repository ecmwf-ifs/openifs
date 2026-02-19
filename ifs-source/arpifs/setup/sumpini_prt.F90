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

SUBROUTINE SUMPINI_PRT

!**** *SUMPINI_PRT*   - Printings of quantities computed in SUMPINI

!     Purpose.
!     --------
!      Prints quantities computed in SUMPINI.
!      Printings cannot be put in SUMPINI, because SUMPINI is called before SULUN,
!      and these printings must be done after SULUN in order to use the right value of NULOUT.
!      Printed variables are in YOMMP0 and YOMGSTATS.

!**   Interface.
!     ----------
!        CALL SUMPINI_PRT

!        explicit arguments : 
!        --------------------

!        implicit arguments :
!        --------------------

!     method.
!     -------
!        see documentation

!     externals.
!     ----------

!     reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     author.
!     -------
!      K. Yessad (July 2013)

! Modifications
! -------------
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
! End Modifications
!------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : MYPROC, MYSETA, MYSETB, MYSETV, MYSETW, MYSETM, MYSETN, &
 & MP_TYPE, MBX_SIZE, NPROC, NPRGPNS, NPRGPEW, NPRTRW, NPRTRV, LMPDIAG, LMPOFF, &
 & NSPECRESMIN, NPRTRNS, NPRTRN, NOUTPUT, LOUTPUT, LOPT_SCALAR, LOPT_RS6K, LSCMEC, NPRINTLEV
USE YOMGSTATS, ONLY : LSTATS, LSTATSCPU, LSYNCSTATS, LDETAILED_STATS, LXML_STATS, &
 & LSTATS_MEM, LSTATS_ALLOC, LBARRIER_STATS, LBARRIER_STATS2, NSTATS_MEM, &
 & NPRNT_STATS, LSTATS_OMP, LSTATS_COMMS, LSTATS_MPL, LTRACE_STATS, LGSTATS_LABEL

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMPINI_PRT',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!     Printings.
!     ----------

WRITE(NULOUT,*) ''
WRITE(NULOUT,*) ' ----- PRINTINGS IN SUMPINI_PRT: '
WRITE(NULOUT,FMT='(''  NPROC = '',I7)') NPROC
WRITE(NULOUT,FMT='(''  NPRGPNS = '',I7)') NPRGPNS
WRITE(NULOUT,FMT='(''  NPRGPEW = '',I7)') NPRGPEW
WRITE(NULOUT,FMT='(''  NPRTRW = '',I7)') NPRTRW
WRITE(NULOUT,FMT='(''  NPRTRV = '',I7)') NPRTRV
WRITE(NULOUT,FMT='(''  NPRTRN = '',I7)') NPRTRN
WRITE(NULOUT,FMT='(''  NPRTRNS = '',I7)') NPRTRNS
WRITE(NULOUT,FMT='(''  MYPROC = '',I7)') MYPROC
WRITE(NULOUT,FMT='(''  MYSETW = '',I7)') MYSETW
WRITE(NULOUT,FMT='(''  MYSETV = '',I7)') MYSETV
WRITE(NULOUT,FMT='(''  MYSETA = '',I7)') MYSETA
WRITE(NULOUT,FMT='(''  MYSETB = '',I7)') MYSETB
WRITE(NULOUT,FMT='(''  MYSETM = '',I7)') MYSETM
WRITE(NULOUT,FMT='(''  MYSETN = '',I7)') MYSETN
WRITE(NULOUT,FMT='(''  NSPECRESMIN = '',I7)') NSPECRESMIN
WRITE(NULOUT,FMT='(''  MP_TYPE = '',I7)') MP_TYPE
WRITE(NULOUT,FMT='(''  MBX_SIZE = '',I12)') MBX_SIZE
WRITE(NULOUT,FMT='(''  NOUTPUT = '',I7)') NOUTPUT
WRITE(NULOUT,FMT='(''  LOUTPUT = '',L2)') LOUTPUT
WRITE(NULOUT,FMT='(''  LMPDIAG = '',L2)') LMPDIAG
WRITE(NULOUT,FMT='(''  LMPOFF = '',L2)') LMPOFF
WRITE(NULOUT,FMT='(''  LOPT_SCALAR = '',L2)') LOPT_SCALAR
WRITE(NULOUT,FMT='(''  LOPT_RS6K = '',L2)') LOPT_RS6K
WRITE(NULOUT,FMT='(''  LSCMEC = '',L2)') LSCMEC
WRITE(NULOUT,FMT='(''  NPRINTLEV = '',I7)') NPRINTLEV

! * GSTATS printing (YOMGSTATS variables, some of them are in NAMPAR0):
WRITE(NULOUT,FMT='(''  LSTATS = '',L2)') LSTATS
WRITE(NULOUT,FMT='(''  LSTATSCPU = '',L2)') LSTATSCPU
WRITE(NULOUT,FMT='(''  LSYNCSTATS = '',L2)') LSYNCSTATS
WRITE(NULOUT,FMT='(''  LDETAILED_STATS = '',L2)') LDETAILED_STATS
WRITE(NULOUT,FMT='(''  LXML_STATS = '',L2)') LXML_STATS
WRITE(NULOUT,FMT='(''  LSTATS_MEM = '',L2)') LSTATS_MEM
WRITE(NULOUT,FMT='(''  LSTATS_ALLOC = '',L2)') LSTATS_ALLOC
WRITE(NULOUT,FMT='(''  LBARRIER_STATS = '',L2)') LBARRIER_STATS
WRITE(NULOUT,FMT='(''  LBARRIER_STATS2 = '',L2)') LBARRIER_STATS2
WRITE(NULOUT,FMT='(''  NSTATS_MEM = '',I7)') NSTATS_MEM
WRITE(NULOUT,FMT='(''  NPRNT_STATS = '',I7)') NPRNT_STATS
WRITE(NULOUT,FMT='(''  LSTATS_OMP = '',L2)') LSTATS_OMP
WRITE(NULOUT,FMT='(''  LSTATS_COMMS = '',L2)') LSTATS_COMMS
WRITE(NULOUT,FMT='(''  LSTATS_MPL = '',L2)') LSTATS_MPL
WRITE(NULOUT,FMT='(''  LTRACE_STATS = '',L2)') LTRACE_STATS
WRITE(NULOUT,FMT='(''  LGSTATS_LABEL = '',L2)') LGSTATS_LABEL

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMPINI_PRT',1,ZHOOK_HANDLE)
END SUBROUTINE SUMPINI_PRT
