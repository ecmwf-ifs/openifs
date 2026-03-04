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

MODULE YOMJBPAR1DECV
!-----------------------------------------------------------------------
!   YOMJBPAR1DECV - Module for physics parameters optimization.

!   Purpose.
!   --------
!   Set up data structures and controls for physics parameters optimization.

!   Author.
!   -------
!   Philippe Lopez (ECMWF)

!   Modifications.
!   --------------
!   Original    2012/02/24
!   S. Massart  2019/01     from parphy_setup.F90
!------------------------------------------------------------------------

USE PARKIND1      , ONLY : JPIM, JPRB
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN        , ONLY : NULOUT, NULNAM, RESERVE_LUN, FREE_LUN
USE YOMMP0        , ONLY : NPROC
USE YOMJBECV      , ONLY : YRECV0, YRECV5, YRECVCONFIG, CKNOWNECV, ECV_1D_TO_2DSP, ECV_2DSP_TO_1D, ECV_2DSP_TO_1D_AD
USE TYPE_ECV      , ONLY : TECVDIM, ECV_CONTAINER

IMPLICIT NONE

PRIVATE :: MODNAME, CDFNAME, ZMISSING, READNAMPARECV
PUBLIC :: SUPARECVTRAJ, SUPARECVMIN, LOAD_TABLE, PARECV_SUM_AD, PARECV_PRINT, PARECV_SAVE

SAVE

TYPE, PUBLIC :: TPAR1DECV
  !     ------------------------------------------------------------------
  ! DATA     physics parameters - latest estimate (not optimal to have it here )
  ! BKGERR   background error std dev
  ! ACHGVAR  change-of-variable operator
  ! ACVARIN  inverse change-of-variable
  ! AGRAD    gradient w/r to physics parameters
  ! AGRAD0   initial gradient w/r to physics parameters
  !     ------------------------------------------------------------------
  REAL(KIND=JPRB) :: DATA
  REAL(KIND=JPRB) :: DATA0
  REAL(KIND=JPRB) :: BKGERR
  REAL(KIND=JPRB) :: ACHGVAR
  REAL(KIND=JPRB) :: ACVARIN
  REAL(KIND=JPRB) :: AGRAD
  REAL(KIND=JPRB) :: AGRAD0
  LOGICAL         :: LSET = .FALSE.
END TYPE TPAR1DECV

TYPE(TPAR1DECV), ALLOCATABLE :: YRPAR1DECV(:)

!     ------------------------------------------------------------------
! MODNAME       Name of the module
! CDFNAME       Default name for the cycling file
! JPMXNPARAM    Max number of physics parameters
! ZMISSING      Missing value indicator
!     ------------------------------------------------------------------

CHARACTER(LEN=*),   PARAMETER :: MODNAME = 'YOMJBPAR1DECV'
CHARACTER(LEN=*),   PARAMETER :: CDFNAME = 'PAR1DECV.cycle'
INTEGER(KIND=JPIM), PARAMETER :: JPMXNPARAM = 5
REAL(KIND=JPRB),    PARAMETER :: ZMISSING=-99999999._JPRB
CHARACTER(LEN=20),  PARAMETER :: CMISSING='UNKNOWN             '
CHARACTER(LEN =20), PARAMETER :: CKNOWNDESC(JPMXNPARAM) = (/  &
  &   'SOLAR_CONSTANT      ', &
  &   CMISSING, &
  &   CMISSING, &
  &   CMISSING, &
  &   CMISSING/)

!-----------------------------------------------------------------------
! Namelist parameters
!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NPAR1DECV            ! total number of physics parameters
REAL(KIND=JPRB)    :: RBGSTDV(JPMXNPARAM)  ! Bkg error std dev for each physics parameter.
REAL(KIND=JPRB)    :: RCSVALS(JPMXNPARAM)  ! Cold start value of each physics parameter.
REAL(KIND=JPRB)    :: RCSVALS0(JPMXNPARAM) ! Cold start value of each physics parameter (background).
CHARACTER(LEN =20) :: CDESC(JPMXNPARAM)    ! Desciption of the parameter
LOGICAL            :: LCOLDSTART           ! Cold start on all physics parameters.
LOGICAL            :: LPRECOND             ! Preconditionning

#include "namparecv.nam.h"

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

SUBROUTINE READNAMPARECV()
  call abor1("oifs/fc-only - READNAMPARECV should never be called")

END SUBROUTINE READNAMPARECV


SUBROUTINE SUPARECVTRAJ(YDGEOMETRY,YDDIMECV)
!-----------------------------------------------------------------------
!   SUPARECVTRAJ - Set up physics parameters optimization

!   Purpose.
!   --------

!   Interface.
!   ----------
!   CALL SUPARECVTRAJ ()

!   Author.
!   -------
!   Philippe Lopez (ECMWF)

!   Modifications.
!   --------------
!   Original    2012/02/24
!   S. Massart  2019/01 Change of the name (from parphy_setup_traj)

!------------------------------------------------------------------------

USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE YOMCT0        , ONLY : LSCREEN

!     -----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY) ,INTENT(IN)   :: YDGEOMETRY
TYPE(TECVDIM)  ,INTENT(IN)   :: YDDIMECV
  call abor1("oifs/fc-only - SUPARECVTRAJ should never be called")

END SUBROUTINE SUPARECVTRAJ

!-----------------------------------------------------------------------
SUBROUTINE SUPARECVMIN (YDGEOMETRY)
!-----------------------------------------------------------------------

! During minimisation, parameter increments are defined with
! respect to the estimates used in the latest trajectory run.
! The minimisation therefore always starts with zero increments
! for the parameters (see SUVAZX). In order to correctly
! compute the background term (if any), the change-of-variable
! routines for the parameter section of the control variable
! therefore require the background parameters as well as the
! estimates used in the latest trajectory run.
!
USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE YOMCT0        , ONLY : LSCREEN

!     -----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN)  :: YDGEOMETRY
call abor1("oifs/fc-only - SUPARECVMIN should never be called")


END SUBROUTINE SUPARECVMIN

!-----------------------------------------------------------------------
SUBROUTINE PARECV_SAVE(CDFILENAME)
!-----------------------------------------------------------------------
!   PARECV_SAVE - Write file with physics parameters information

!   Purpose.
!   --------
!   Cycling of physics parameter information

!   Optional: cdfilename  - Name of physics parameters file where to write

!   Author.
!   -------
!   Philippe Lopez (ECMWF)

!   Modifications.
!   --------------
!   Original    2012/02/27
!-----------------------------------------------------------------------

USE YOMMP0  , ONLY : NPROC, MYPROC
USE YOMCT0  , ONLY : CNMEXP
USE YOMANA  , ONLY : NANDAT, NANTIM
USE ALGORITHM_STATE_MOD  , ONLY : GET_NUPTRA

IMPLICIT NONE

CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDFILENAME
  call abor1("oifs/fc-only - PARECV_SAVE should never be called")


END SUBROUTINE PARECV_SAVE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE LOAD_TABLE ()
!-----------------------------------------------------------------------
call abor1("oifs/fc-only - LOAD_TABLE should never be called")


END SUBROUTINE LOAD_TABLE

SUBROUTINE PARECV_PRINT (KTABLE, LDVERBOSE)
!-----------------------------------------------------------------------
!   PARECV_PRINT - Print physics parameters information

!   Purpose.
!   --------
!   Print physics parameters information

!   Author.
!   -------
!   Philippe Lopez (ECMWF)

!   Modifications.
!   --------------
!   Original    2012/02/29
!   Sebastien Massart 2019/01
!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KTABLE
LOGICAL,INTENT(IN),OPTIONAL   :: LDVERBOSE
call abor1("oifs/fc-only - PARECV_PRINT should never be called")


END SUBROUTINE PARECV_PRINT
!-----------------------------------------------------------------------------

SUBROUTINE PARECV_SUM_AD(YDDIMECV,YDECV)
!-----------------------------------------------------------------------
!   PARECV_SUM_AD - Sum-up adjoint contributions on physics parameters
!                   across processors.

!   Author.
!   -------
!   Philippe Lopez (ECMWF)

!   Modifications.
!   --------------
!   Original    2012/03/02
!   Sebastien Massart 2019/01
!-----------------------------------------------------------------------

USE YOMVAR  , ONLY : LREPRO4DVAR
USE MPL_MODULE, ONLY: MPL_ALLREDUCE

IMPLICIT NONE

TYPE(TECVDIM)  ,INTENT(IN)   :: YDDIMECV
!INTEGER(KIND=JPIM), INTENT(IN)    :: K1D
TYPE(ECV_CONTAINER),INTENT(INOUT) :: YDECV
call abor1("oifs/fc-only - PARECV_SUM_AD should never be called")


END SUBROUTINE PARECV_SUM_AD
!-----------------------------------------------------------------------------

SUBROUTINE GET_VALUE_FROM_ECV(YDGEOMETRY, YDECV, CDDESC, PP1D, LDFOUND)
!-----------------------------------------------------------------------
!   GET_VALUE_FROM_ECV - Get the parameter value from ECV object

!   Author.
!   -------
!   Sebastien Massart (ECMWF)

!   Modifications.
!   --------------
!   Original    2019/01
!-----------------------------------------------------------------------

USE GEOMETRY_MOD  , ONLY : GEOMETRY

IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN)  :: YDGEOMETRY
TYPE(ECV_CONTAINER),INTENT(IN)  :: YDECV
CHARACTER(LEN =*),  INTENT(IN)  :: CDDESC
REAL(KIND=JPRB),    INTENT(OUT) :: PP1D
LOGICAL,            INTENT(OUT) :: LDFOUND
call abor1("oifs/fc-only - GET_VALUE_FROM_ECV should never be called")


END SUBROUTINE GET_VALUE_FROM_ECV
!-----------------------------------------------------------------------------

SUBROUTINE GET_VALUE_FROM_ECV_AD(YDGEOMETRY, YDECV, CDDESC, PP1D, LDFOUND)
!-----------------------------------------------------------------------
!   GET_VALUE_FROM_ECV_AD - Get the parameter value from ECV object (AD VERSION)

!   Author.
!   -------
!   Sebastien Massart (ECMWF)

!   Modifications.
!   --------------
!   Original    2019/01
!-----------------------------------------------------------------------

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE MPL_MODULE  , ONLY : MPL_ALLREDUCE
USE YOMMP0  , ONLY : MYPROC

IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN)    :: YDGEOMETRY
TYPE(ECV_CONTAINER),INTENT(INOUT) :: YDECV
CHARACTER(LEN =*),  INTENT(IN)    :: CDDESC
REAL(KIND=JPRB),    INTENT(INOUT) :: PP1D
LOGICAL,            INTENT(OUT)   :: LDFOUND
call abor1("oifs/fc-only - GET_VALUE_FROM_ECV_AD should never be called")


END SUBROUTINE GET_VALUE_FROM_ECV_AD
!-----------------------------------------------------------------------------

END MODULE YOMJBPAR1DECV
