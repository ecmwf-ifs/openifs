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

MODULE YOEWCOU

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*    ** *YOEWCOU* - VARIABLES FOR COUPLING WITH THE WAVE MODEL

!     P. VITERBO     E.C.M.W.F.       07/10/88
!     P. VITERBO     E.C.M.W.F.       03/02/92
!     J. DOYLE       E.C.M.W.F.       21/11/96 
!     J. BIDLOT      E.C.M.W.F.       13/06/97 
!     J. BIDLOT      E.C.M.W.F.       11/08/98 
!     G. MOZDZYNSKI  E.C.M.W.F.       14/01/05
!     M. Drusch      E.C.M.W.F.       25/09/07
!     J. BIDLOT      E.C.M.W.F.       22/08/08 
!     H. Hersbach    E.C.M.W.F.       13/10/08 add LWCUR
!     J. BIDLOT      E.C.M.W.F.       01/06/09 add LWSTOKES
!     J. BIDLOT      E.C.M.W.F.       11/07/10 add NGRIB_HANDLE_FOR_WAM 
!     P. JANSSEN     E.C.M.W.F.       20/06/11 add LWFLUX 

!      NAME      TYPE       PURPOSE
!      ----      ----       -------

!     *NATMFLX*  INTEGER    TO SPECIFY WHICH COMPONENT OF THE MOMEUMTUM FLUX
!                           IS PASSED TO THE WAVE MODEL:
!                           1: FRICTION VELOCITY VECTOR
!                           2: SURFACE STRESS VECTOR
!                           3: 10m NEUTRAL WIND VECTOR
!     *NGAUSSW*  INTEGER    = 1 : WAVE GRID IS GAUSSIAN
!                           = 0 : WAVE GRID IS LAT-LON
!     *NLONW*    INTEGER    NUMBER OF POINTS IN A LATITUDE LINE IN
!                           THE WAVE MODEL
!     *NLATW*    INTEGER    NUMBER OF LATITUDES IN THE WAVE MODEL.
!     *NLON1W*   INTEGER    *NLONW*
!     *NLAT1W*   INTEGER    TOTAL NUMBER OF LATITUDES WITH THE WAVE
!                           MODEL RESOLUTION.
!     *NNORXW*   INTEGER    NUMBER OF EXTRA POINTS NORTHWARDS OF THE
!                           NORTHERN BOUNDARY OF THE WAVE MODEL.
!     *CBEGDAT*  CHARACTER  INITIAL DATE OF FORECAST (YYYYMMDDHHmm)
!     *NSTPW*    INTEGER    FREQUENCY OF CALL TO THE WAVE MODEL.
!     *LWCOU*    LOGICAL    TRUE IF THE WAVE MODEL IS TO BE RUN.
!     *LWCOU2W*  LOGICAL    TRUE IF TWO-WAY INTERACTION WITH THE WAVE MODEL.
!                           FALSE IF ONE-WAY INTERACTION WITH THE WAVE MODEL
!     *LWCOUNORMS* LOGICAL  TRUE IF NORMS OF COUPLED FIELDS ARE REQUIRED
!     *LWCUR*    LOGICAL    TRUE IF OCEAN CURRENTS ARE PASSED TO THE WAVE MODEL.
!     *LWCOURNW* LOGICAL    TRUE IF 10m NEUTRAL WINDS RELATIVE TO OCEAN CURRENTS ARE PASSED THE THE WAVE MODEL
!!                          OTHERWISE STRESS EQUIVALENT 10m NEUTRAL WINDS WILL BE PASSED.
!     *LWCOUHMF* LOGICAL    TRUE IF OCEAN WAVES WILL ALSO DIRECTLY INFLUENCE THE HEAT AND MOISTURE FLUXES
!                           (ONLY ACTIVE IF LWCOU2W IS TRUE)  
!     *LWSTOKES* LOGICAL    TRUE IF SURFACE STOKES DRIFT IS RETURNED TO IFS.
!     *LWFLUX*   LOGICAL    TRUE IF WAVE FLUXES/PARAMETERS NEEDED FOR WHITECAP PARAMETERISATION ARE RETURNED 
!     *NGRIB_HANDLE_FOR_WAM GRIB HANDLE USED TO EXCHANGE GRIB INFORMATION
!                           AND IFS MODEL GRID INFORMATION. 
!     *RSOUTW*   REAL       SOUTH BOUNDARY OF THE WAVE MODEL.
!     *RNORTW*   REAL       NORTH BOUNDARY OF THE WAVE MODEL.
!     *RDEGREW*  REAL       RESOLUTION OF THE WAVE MODEL (DEGREES).

!     *RMISSW*   REAL       MISSING DATA INDICATOR FOR THE WAVE MODEL.

!     *MASK_WAVE_IN*  INTEGER  COMMS MASK FOR INPUT TO WAVE MODEL
!     *MASK_WAVE_OUT* INTEGER  COMMS MASK FOR OUTPUT FROM WAVE MODEL

!     *LWVIN_MASK_NOT_SET* LOGICAL indicates whether mask_wave_in
!                           has been updated on the first call to the 
!                           wave model
!     *LWVOUT_MASK_NOT_SET* LOGICAL indicates whether mask_wave_out
!                           has been updated following the first call to the 
!                           wave model
!     *LWVIN_UNINITIALISED* LOGICAL indicates whether the mwvin_* data
!                           structures are initialised
!     *MWVIN_SENDCNT* INTEGER nproc sized array describing how many grid 
!                           points need to be sent (or copied) by the 
!                           local task to remote tasks (or this task)
!     *MWVIN_RECVCNT* INTEGER nproc sized array describing how many grid 
!                           points need to be received by the local task 
!                           from a remote task
!     *MWVIN_SENDTOT* INTEGER total number of grid points to be sent to 
!                           remote tasks
!     *MWVIN_RECVTOT* INTEGER total number of grid points to be received 
!                           from remote tasks
!     *MWVIN_SENDOFF* INTEGER nproc sized array containing offsets into 
!                           the MWVIN_SENDBUF and MWVIN_SENDIND arrays
!     *MWVIN_SENDBUF* INTEGER local indexes of data on remote tasks that
!                           the local task needs
!     *MWVIN_SENDIND* INTEGER global indexes of data on remote tasks that
!                           the local task needs

!     *MWVIN_RECVOFF* INTEGER nproc sized array containing offsets into 
!                           the MWVIN_RECVBUF array
!     *MWVIN_RECVBUF* INTEGER local indexes of data on the local task that
!                           remote tasks need
!     *NWV_W2IWGHT*   INTEGER FIRST DIMENSION OF WV_W2IWGHT 
!     *IWV_LON*       INTEGER FIRST INDEX (X) IN INTERPOLATION OF WAVE GRID TO
!                             IFS GRID.
!     *IWV_LAT*       INTEGER SECOND INDEX (Y) IN INTERPOLATION OF WAVE GRID TO
!                             IFS GRID.
!     *WV_W2IWGHT*    REAL  WEIGHT FOR THE INTERPOLATION TO THE IFS GRID
!                           OF WAVE MODEl FIELD(S) RETURNED FROM WAM.

TYPE :: TEWCOU
INTEGER(KIND=JPIM) :: NATMFLX
INTEGER(KIND=JPIM) :: NGAUSSW
INTEGER(KIND=JPIM) :: NLONW
INTEGER(KIND=JPIM) :: NLATW
INTEGER(KIND=JPIM) :: NLON1W
INTEGER(KIND=JPIM) :: NLAT1W
INTEGER(KIND=JPIM) :: NNORXW
INTEGER(KIND=JPIM) :: NSTPW

INTEGER(KIND=JPIM) :: NGRIB_HANDLE_FOR_WAM

REAL(KIND=JPRB) :: RSOUTW
REAL(KIND=JPRB) :: RNORTW
REAL(KIND=JPRB) :: RDEGREW

REAL(KIND=JPRB) :: RMISSW

LOGICAL :: LWCOU
LOGICAL :: LWCOU2W
LOGICAL :: LWCOUNORMS
LOGICAL :: LWCUR
LOGICAL :: LWCOURNW 
LOGICAL :: LWCOUHMF
LOGICAL :: LWSTOKES
LOGICAL :: LWFLUX
CHARACTER :: CBEGDAT*14

INTEGER(KIND=JPIM), ALLOCATABLE :: MASK_WAVE_IN(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MASK_WAVE_OUT(:,:)

LOGICAL :: LWVIN_MASK_NOT_SET
LOGICAL :: LWVOUT_MASK_NOT_SET
LOGICAL :: LWVIN_UNINITIALISED

INTEGER(KIND=JPIM) :: MWVIN_SENDTOT
INTEGER(KIND=JPIM) :: MWVIN_RECVTOT

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDCNT(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_RECVCNT(:)

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_RECVOFF(:)

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDBUF(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_RECVBUF(:)

INTEGER(KIND=JPIM),ALLOCATABLE :: MWVIN_SENDIND(:)

LOGICAL :: LWFRSTIME
LOGICAL :: LWINIT

INTEGER(KIND=JPIM) :: NADV_WAVE

INTEGER(KIND=JPIM) :: NWV_W2IWGHT 
INTEGER(KIND=JPIM),ALLOCATABLE :: IWV_LON(:,:) 
INTEGER(KIND=JPIM),ALLOCATABLE :: IWV_LAT(:,:) 
REAL(KIND=JPRB),ALLOCATABLE :: WV_W2IWGHT(:,:) 

CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TEWCOU

!!TYPE(TEWCOU), POINTER :: YREWCOU => NULL()

  !---------------------------------------------------------------------

  CONTAINS 

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TEWCOU),      INTENT(IN) :: SELF
  INTEGER(KIND=JPIM), INTENT(IN) :: KDEPTH
  INTEGER(KIND=JPIM), INTENT(IN) :: KOUTNO
  END SUBROUTINE PRINT_CONFIGURATION

!     ------------------------------------------------------------------
END MODULE YOEWCOU
