! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE YOMGSTATS

USE EC_PARKIND  ,ONLY : JPRD, JPIM

IMPLICIT NONE

SAVE

PRIVATE :: JPRD, JPIM

!     ------------------------------------------------------------------
! Module for timing statistics. Module is internal to the GSTATS package -
! routines GSTATS, SUSTATS and STATS_OUTPUT. The logical switches are
! re-initialized in SUMPINI

! LSTATS - TRUE for gathering timing statistics
! LSTATSCPU - TRUE for gathering CPU timing  statistics
! LSYNCSTATS - TRUE for syncronization (call to barrier) at the 
!              start of timing event
! LDETAILED_STATS - TRUE for more detail in output
! LXML_STATS - TRUE for stats output in XML
! LSTATS_OMP - TRUE for gathering timing statistics on OpenMP regions
!                 1001-1999
! LSTATS_COMMS - TRUE for gathering detailed timing of Message passing
!                 501-1000
! LSTATS_MPL   - TRUE for gathering detailed info on message passing
! NTRACE_STATS    - max number of entries in trace
! LTRACE_STATS    - True for trace of all calls to gstats
! LGSTATS_LABEL   - True after GSTATS-labels have been set
! JPMAXSTAT - max number of separate  timers in gstats
! JPOBCOUNT_BASE - first counter for obs types
! NCALLS - number of times a timer has been switched on
! TIMESUM - total time spent with timer on
! TIMESQSUM - sum of the squares of times
! TIMEMAX - max time of all calls
! TIMESUMB - sum of times between previous timer was invoked and this
!            timer was switched on ( to be used for finding out which parts
!            of the code that is not being timed)
! TIMELCALL - time when event was switched on or resumed
! TTCPUSUM - total cpu time
! TVCPUSUM - total vector cpu time
! THISTIME - total accumulated time for this call to timing event (necessary
!            to be able to suspend and resume timer and still have it counted
!            as one timing event)
! THISTCPU - as THISTIME but for CPU time
! THISVCPU - as THISTIME but for vector CPU time
! TTCPULCALL - as TIMELCALL but for CPU time
! TVCPULCALL - as TIMELCALL but for vector CPU time
! TIME_LAST_CALL - last time GSTATS was called
! TIME_START - used for recording parallel startup time
!
! NSWITCHVAL - for detecting overlapping counters


LOGICAL :: LSTATS = .TRUE.
LOGICAL :: LSTATS_OMP = .FALSE.
LOGICAL :: LSTATS_COMMS = .FALSE.
LOGICAL :: LSTATS_MPL = .FALSE.
LOGICAL :: LSTATS_MEM = .FALSE.
LOGICAL :: LSTATS_ALLOC = .FALSE.
LOGICAL :: LSTATSCPU = .TRUE.
LOGICAL :: LSYNCSTATS = .FALSE.
LOGICAL :: LXML_STATS = .FALSE.
LOGICAL :: LDETAILED_STATS = .TRUE.
LOGICAL :: LBARRIER_STATS = .FALSE.
LOGICAL :: LBARRIER_STATS2 = .FALSE.
LOGICAL :: LTRACE_STATS = .FALSE.
LOGICAL :: LGSTATS_LABEL = .FALSE.

INTEGER(KIND=JPIM),PARAMETER :: JBMAXBASE=2500
INTEGER(KIND=JPIM),PARAMETER :: JPMAXBARS=500
INTEGER(KIND=JPIM),PARAMETER :: JPMAXSTAT=JBMAXBASE+JPMAXBARS

INTEGER(KIND=JPIM),PARAMETER :: JPOBCOUNT_BASE=201
INTEGER(KIND=JPIM) :: NTRACE_STATS=0
INTEGER(KIND=JPIM) :: NCALLS(0:JPMAXSTAT)
INTEGER(KIND=JPIM) :: NSWITCHVAL(0:JPMAXSTAT)
INTEGER(KIND=JPIM) :: NCALLS_TOTAL=0
INTEGER(KIND=JPIM) :: LAST_KSWITCH=0
INTEGER(KIND=JPIM) :: LAST_KNUM=0
INTEGER(KIND=JPIM) :: NHOOK_MESSAGES=0
INTEGER(KIND=JPIM) :: NBAR_PTR(0:JPMAXSTAT)=0
INTEGER(KIND=JPIM) :: NBAR2=JBMAXBASE+1
INTEGER(KIND=JPIM),ALLOCATABLE :: NCALL_TRACE(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NUMSEND(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NUMRECV(:)
REAL(KIND=JPRD),ALLOCATABLE :: SENDBYTES(:)
REAL(KIND=JPRD),ALLOCATABLE :: RECVBYTES(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: UNKNOWN_NUMSEND(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: UNKNOWN_NUMRECV(:)
REAL(KIND=JPRD),ALLOCATABLE :: UNKNOWN_SENDBYTES(:)
REAL(KIND=JPRD),ALLOCATABLE :: UNKNOWN_RECVBYTES(:)

REAL(KIND=JPRD) :: TIMESUM(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TIMESQSUM(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TIMEMAX(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TIMESUMB(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TIMELCALL(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TTCPUSUM(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TVCPUSUM(0:JPMAXSTAT)
REAL(KIND=JPRD) :: THISTIME(0:JPMAXSTAT)
REAL(KIND=JPRD) :: THISTCPU(0:JPMAXSTAT)
REAL(KIND=JPRD) :: THISVCPU(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TTCPULCALL(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TVCPULCALL(0:JPMAXSTAT)
REAL(KIND=JPRD) :: TIME_LAST_CALL

REAL(KIND=JPRD),ALLOCATABLE :: TIME_START(:)
REAL(KIND=JPRD),ALLOCATABLE :: TIME_TRACE(:)
INTEGER(KIND=JPIM),PARAMETER :: JPERR=0
INTEGER(KIND=JPIM),PARAMETER :: JPTAGSTAT=20555

INTEGER(KIND=JPIM),PARAMETER :: JPMAXDELAYS=1000
INTEGER(KIND=JPIM) :: NDELAY_COUNTER(1:JPMAXDELAYS)
REAL(KIND=JPRD)    :: TDELAY_VALUE(1:JPMAXDELAYS)
CHARACTER*10       :: CDELAY_TIME(1:JPMAXDELAYS)
INTEGER(KIND=JPIM) :: NDELAY_INDEX = 0

CHARACTER*50 :: CCDESC(0:JPMAXSTAT) = ""
CHARACTER*3  :: CCTYPE(0:JPMAXSTAT) = ""

INTEGER(KIND=JPIM) :: NPROC_STATS = 1
INTEGER(KIND=JPIM) :: MYPROC_STATS = 1
INTEGER(KIND=JPIM),ALLOCATABLE :: NPRCIDS_STATS(:)

INTEGER(KIND=JPIM) :: NTMEM(0:JPMAXSTAT,5)
INTEGER(KIND=JPIM) :: NSTATS_MEM=0

INTEGER(KIND=JPIM) :: NPRNT_STATS=3

END MODULE YOMGSTATS




