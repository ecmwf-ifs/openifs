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

      INTEGER FUNCTION IUSRCL()
      USE PARKIND1, ONLY : JPRD, JPIM
      INTEGER(KIND=JPIM) :: TIMENOWHIGH
!RJ       REAL(KIND=JPRB) :: AMICCLK
!
!   MACHINE-DEPENDENT TIMING ROUTINE
!
!----------------------------------------------------------------------
! Author.
! -------
!    Ute Gaertel, Wolfgang Joppich, Anton Schueller   *GMD*
!
! Modifications.
! -------------
!    Original : 93-01-15
!    Modified : 93-03-15
!    Modified : 93-08-25
!      F. Vana  05-Mar-2015  Support for single precision
!----------------------------------------------------------------------
!
!------CM 5
!
!CM5  integer cmmd_node_timer_clear,cmmd_node_timer_start,
!CM5 *        cmmd_node_timer_stop,i,istat
!CM5  double precision time,cmmd_node_timer_elapsed
!CM5  data istat /0/
!CM5  save istat
!CM5
!CM5  if (istat.eq.0) then
!CM5    istat = 1
!CM5    i = cmmd_node_timer_clear(0)
!CM5    i = cmmd_node_timer_start(0)
!CM5  endif
!CM5    i = cmmd_node_timer_stop(0)
!CM5    time = cmmd_node_timer_elapsed(0)
!CM5    iusrcl = nint(time*1.0d3)
!CM5    i = cmmd_node_timer_start(0)
      REAL(KIND=JPRD),EXTERNAL :: UTIL_WALLTIME
      REAL(KIND=JPRD),EXTERNAL :: UTIL_CPUTIME
      IUSRCL = 0
!
!------iPSC/2,  iPSC/860
!
!     IUSRCL = MCLOCK()
!
!------Parsytec GCel
!
!     IUSRCL = TIMENOWHIGH()
!
!------IBM RS 6000, Meiko, SGI
!
!CSB  IUSRCL = 1000 * SECOND()
!cli  IUSRCL = 10 * MCLOCK()
!     IUSRCL =      MCLOCK()
#if defined(SGI) || defined(VPP)
      IUSRCL = 1000 * UTIL_WALLTIME()
!      IUSRCL = 1000 * util_cputime()
#endif
!
!
!------SUN
!
!3LI  IUSRCL = 1000
!
!
!------nCUBE/2
!
!     IUSRCL = INT (AMICCLK()/1000)
#ifdef CRAY
#ifndef T3D
!
!------Cray C90
!
      IUSRCL = INT( RTC()*4.186E-6_JPRD )
#else
!------Cray T3D
!
      IUSRCL = INT( RTC()*6.667E-6_JPRD )
!--      IUSRCL = 1000 * util_walltime()
#endif
#endif
!


      ENDFUNCTION IUSRCL
