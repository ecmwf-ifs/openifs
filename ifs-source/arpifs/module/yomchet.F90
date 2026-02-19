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

MODULE YOMCHET

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE   

! **** *YOMCHET* 

!---------------------------------------------------------------- 
!   Purpose
!   -------
!   Definitions of types and global variables for diagnostics
!   on physical tendencies.

!**   Interface.
!     ----------
!        *CALL* *YOMCHET(...)*

!   References
!   ----------

!   Author
!   ------
!   2004-03-01: T.Kovacic, J.M.Piriou, F.Bouyssel

!   Modifications
!   -------------

!----------------------------------------------------------------

! Global parameters
! -----------------
! MAXPHVAR    : maximal number of diagnosed physical tendencies
! NAMELENTH   : length of names of diagnosed physical tendencies

! Components of TY_PHYS_TEND_PT
! -----------------------------
! %ALAT       : point's latitude 
! %ALONG      : point's longitude
! %HEIGHTSEA  : point's height
! %TEND_VAL   : value of tendency exceeding thresholds
! %NCHETLONG  : position in NPROMA vector

! Components of TYPE_CHETN
! -----------------------------
! %LFREQD     : switch for frequency distrbution calculations
! %LCOORD     : switch for extraction of geographical points
! %LPROFV     : switch for extraction of vertical profiles
! %CPHYST     : determine physical tendencies on which diagnostics are made
! %CPTSEP     : diagnostics on a single or all physical parametrizations
! %NCOORMAX   : max number of tend. exceeding the thresholds in NPROMA vectors
! %NPROFMAX   : max number of vertical profiles extracted in NPROMA vectors
! %NLEV1      : top level on the vertical
! %NLEV2      : bottom level on the vertical
! %APREC      : precision, antilogarithm of class width on log-scale
! %AWIDTH     : max value on the log-scale
! %TENDTHRESH : thresholds on physical tendencies
! %STARTCLASS : starting values for classes

! Components of TYPE_CHET
! -----------------------------
! %CFICFREQ   : output file name for writing frequency distribution
! %CFICCOOR   : output file name for writing geographical coordinates
! %CVARNA     : array with names of variables
! %MAXCLASS   : maximal number of classes for frequency distribution
! %NULFREQ    : logical unit number for writing frequency distribution
! %NULCOOR    : logical unit number for writing geographical coordinates
! %NULPROF    : logical unit number for writing vertical profiles
! %NVPHT      : number of physical tendencies on which diagnostics are made
! %NPHTEXTH   : number of points with phys. tend. exceeding the thresholds
! %NWIDTH     : number of classes on the log-scale
! %NQDIST     : array with frequency distr. of phys. tend. for chosen variables
! %ACLASS     : array with class values
! %APOINT     : array with geographical coordinates exceeding the thresholds

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!      GLOBAL PARAMETER DEFINITIONS
!------------------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: MAXPHVAR  = 8
INTEGER(KIND=JPIM), PARAMETER :: NAMELENTH = 20

!------------------------------------------------------------------------------
!     DEFINITION OF THE NEW TYPES 
!------------------------------------------------------------------------------

TYPE TY_PHYS_TEND_PT
REAL(KIND=JPRB)                                  :: ALAT
REAL(KIND=JPRB)                                  :: ALONG
REAL(KIND=JPRB)                                  :: HEIGHTSEA
REAL(KIND=JPRB)                                  :: TEND_VAL
INTEGER(KIND=JPIM)                               :: NCHETLONG
END TYPE TY_PHYS_TEND_PT

TYPE TYPE_CHETN
LOGICAL                                          :: LFREQD
LOGICAL                                          :: LCOORD
LOGICAL                                          :: LPROFV
CHARACTER(20)                                    :: CPHYST
CHARACTER(4)                                     :: CPTSEP
INTEGER(KIND=JPIM)                               :: NCOORMAX
INTEGER(KIND=JPIM)                               :: NPROFMAX
INTEGER(KIND=JPIM)                               :: NLEV1
INTEGER(KIND=JPIM)                               :: NLEV2
REAL(KIND=JPRB)                                  :: APREC
REAL(KIND=JPRB)                                  :: AWIDTH
REAL(KIND=JPRB),DIMENSION(MAXPHVAR)              :: TENDTHRESH
REAL(KIND=JPRB),DIMENSION(MAXPHVAR)              :: STARTCLASS
END TYPE TYPE_CHETN

TYPE TYPE_CHET
CHARACTER(80)                                    :: CFICFREQ
CHARACTER(80)                                    :: CFICCOOR
CHARACTER(NAMELENTH), DIMENSION(MAXPHVAR)        :: CVARNA
INTEGER(KIND=JPIM)                               :: MAXCLASS
INTEGER(KIND=JPIM)                               :: NULFREQ
INTEGER(KIND=JPIM)                               :: NULCOOR
INTEGER(KIND=JPIM)                               :: NULPROF
INTEGER(KIND=JPIM)                               :: NVPHT
INTEGER(KIND=JPIM)                               :: NPHTEXTH
INTEGER(KIND=JPIM)                               :: NWIDTH
INTEGER(KIND=JPIM),POINTER,DIMENSION(:,:)        :: NQDIST=>NULL()
REAL(KIND=JPRB),POINTER,DIMENSION(:,:)           :: ACLASS=>NULL()
TYPE(TY_PHYS_TEND_PT),POINTER,DIMENSION(:)       :: APOINT=>NULL()
END TYPE TYPE_CHET

!------------------------------------------------------------------------------
!      GLOBAL VARIABLE DEFINITIONS
!------------------------------------------------------------------------------

TYPE(TYPE_CHETN)                                 :: GCHETN
TYPE(TYPE_CHET)                                  :: GCHET

END MODULE YOMCHET
