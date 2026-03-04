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

MODULE ALGORITHM_STATE_MOD
! Module containing the current state of the algorithm
! For incremental 4Dvar typically outer loop number, iteration count etc.
USE PARKIND1, ONLY : JPIM, JPRB
IMPLICIT NONE

PRIVATE

PUBLIC :: SETUP_ALGORITHM_STATE,GET_ALGOR_TYPE,GET_NUPTRA,SET_NUPTRA, &
 & GET_MUPTRA,SET_MUPTRA,GET_NSIM4D,SET_NSIM4D,GET_NITER4D,SET_NITER4D,&
 & L_IN_TEST,SET_IN_TEST,L_IN_FINAL,SET_IN_FINAL,L_OBS_IN_FC,SET_OBS_IN_FC,SET_IN_WRITE_INCREMENT,&
 & L_IN_WRITE_INCREMENT,SET_IN_WRITE_STATE,L_IOSTREAM_FDB_SET,SET_IOSTREAM_FDB,&
 & L_IN_MINIM,SET_IN_MINIM,&
 & L_IN_WRITE_STATE,L_TRAJ_COMPLETE,SET_L_TRAJ_COMPLETE,L_OUTER_COMPLETE,SET_L_OUTER_COMPLETE,&
 & L_4V_WRITTEN,SET_L_4V_WRITTEN,L_LATE4D_CHECK_COMPLETE,SET_L_LATE4D_CHECK_COMPLETE, &
 & L_IN_LOWRES_STATE_CTOR, SET_L_IN_LOWRES_STATE_CTOR, &
 & L_IN_RESTART_CTLVEC, SET_L_IN_RESTART_CTLVEC, &
 & L_IFSTRAJ_CHKPT_COMPLETE, SET_L_IFSTRAJ_CHKPT_COMPLETE

TYPE :: ALGORITHM_STATE_TYPE
CHARACTER(LEN=:),ALLOCATABLE :: CALGOR_TYPE
LOGICAL :: L_ACTIVE=.FALSE.                       ! TRUE when initialized
LOGICAL :: L_IN_TEST_CONFIG=.FALSE.               ! Currently in some test configuration
LOGICAL :: L_IN_FINAL_CONFIG=.FALSE.              ! Currently in final non-linear integration
LOGICAL :: L_OBS_IN_FC_CONFIG=.FALSE.
LOGICAL :: L_IN_WRITE_INCREMENT=.FALSE.
LOGICAL :: L_IN_WRITE_STATE=.FALSE.
LOGICAL :: L_IOSTREAM_FDB_SET=.FALSE.
LOGICAL :: L_IN_MINIM=.FALSE.                     ! TRUE while in minimization
LOGICAL :: L_TRAJ_COMPLETE=.FALSE.                ! Nonlinear HR trajectory for  current outer loop is complete
LOGICAL :: L_OUTER_COMPLETE=.FALSE.               ! Current outer loop is complete
LOGICAL :: L_4V_WRITTEN=.FALSE.                   ! 4V fields have been written
LOGICAL :: L_IN_LOWRES_STATE_CTOR=.FALSE.         ! Construction of low resolution state
LOGICAL :: L_LATE4D_CHECK_COMPLETE=.FALSE.        ! late4dstart check has finished
LOGICAL :: L_IN_RESTART_CTLVEC=.FALSE.            ! Restart of ctl vecs
LOGICAL :: L_IFSTRAJ_CHKPT_COMPLETE=.FALSE.       ! Checkpoint allowing to start from minimization done
INTEGER(KIND=JPIM) :: NO_OF_UPDATES_ADDED=-99     ! Number of increments added (in incremental VAR) (was NUPTRA)
INTEGER(KIND=JPIM) :: MAX_NO_OF_UPDATES_ADDED=-99 ! Total number of increments added (in incremental VAR) (was MUPTRA)
INTEGER(KIND=JPIM) :: NGRADIENTS_COMPUTED=0       ! Number of gradients computed in minimuzation (was NSIM4D)
INTEGER(KIND=JPIM) :: NCURRENT_INNER_LOOP=0       ! Inner (minimuzation) loop number (was NITER4D)
END TYPE ALGORITHM_STATE_TYPE

TYPE(ALGORITHM_STATE_TYPE),PROTECTED :: ALGORITHM_STATE
CONTAINS
!======================================================================
SUBROUTINE SETUP_ALGORITHM_STATE(CDALGOR_TYPE)
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN) :: CDALGOR_TYPE
ALGORITHM_STATE%CALGOR_TYPE = CDALGOR_TYPE
ALGORITHM_STATE%L_ACTIVE = .TRUE.
END SUBROUTINE SETUP_ALGORITHM_STATE
!======================================================================
PURE FUNCTION GET_ALGOR_TYPE()
IMPLICIT NONE
CHARACTER(LEN=:),ALLOCATABLE :: GET_ALGOR_TYPE
IF (ALGORITHM_STATE%L_ACTIVE) THEN
  GET_ALGOR_TYPE = ALGORITHM_STATE%CALGOR_TYPE
ELSE
  GET_ALGOR_TYPE = 'OTHER'
ENDIF
END FUNCTION GET_ALGOR_TYPE
!======================================================================
PURE FUNCTION GET_NUPTRA()
IMPLICIT NONE
INTEGER(KIND=JPIM) :: GET_NUPTRA
GET_NUPTRA = ALGORITHM_STATE%NO_OF_UPDATES_ADDED
END FUNCTION GET_NUPTRA
!======================================================================
SUBROUTINE SET_NUPTRA(KUPTRA)
INTEGER(KIND=JPIM),INTENT(IN) :: KUPTRA
ALGORITHM_STATE%NO_OF_UPDATES_ADDED = KUPTRA
END SUBROUTINE SET_NUPTRA
!======================================================================
PURE FUNCTION GET_MUPTRA()
IMPLICIT NONE
INTEGER(KIND=JPIM) :: GET_MUPTRA
GET_MUPTRA = ALGORITHM_STATE%MAX_NO_OF_UPDATES_ADDED
END FUNCTION GET_MUPTRA
!======================================================================
SUBROUTINE SET_MUPTRA(KMUPTRA)
INTEGER(KIND=JPIM),INTENT(IN) :: KMUPTRA
ALGORITHM_STATE%MAX_NO_OF_UPDATES_ADDED = KMUPTRA
END SUBROUTINE SET_MUPTRA
!======================================================================
PURE FUNCTION GET_NSIM4D()
IMPLICIT NONE
INTEGER(KIND=JPIM) :: GET_NSIM4D
GET_NSIM4D = ALGORITHM_STATE%NGRADIENTS_COMPUTED
END FUNCTION GET_NSIM4D
!======================================================================
SUBROUTINE SET_NSIM4D(KSIM4D)
INTEGER(KIND=JPIM),INTENT(IN) :: KSIM4D
ALGORITHM_STATE%NGRADIENTS_COMPUTED = KSIM4D
END SUBROUTINE SET_NSIM4D
!======================================================================
PURE FUNCTION GET_NITER4D()
IMPLICIT NONE
INTEGER(KIND=JPIM) :: GET_NITER4D
GET_NITER4D = ALGORITHM_STATE%NCURRENT_INNER_LOOP
END FUNCTION GET_NITER4D
!======================================================================
SUBROUTINE SET_NITER4D(KITER4D)
INTEGER(KIND=JPIM),INTENT(IN) :: KITER4D
ALGORITHM_STATE%NCURRENT_INNER_LOOP = KITER4D
END SUBROUTINE SET_NITER4D
!======================================================================
PURE FUNCTION L_IN_TEST()
LOGICAL :: L_IN_TEST
L_IN_TEST = ALGORITHM_STATE%L_IN_TEST_CONFIG
END FUNCTION L_IN_TEST
!======================================================================
SUBROUTINE SET_IN_TEST(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IN_TEST_CONFIG = LD
END SUBROUTINE SET_IN_TEST
!======================================================================
PURE FUNCTION L_IN_FINAL()
LOGICAL :: L_IN_FINAL
L_IN_FINAL = ALGORITHM_STATE%L_IN_FINAL_CONFIG
END FUNCTION L_IN_FINAL
!======================================================================
SUBROUTINE SET_IN_FINAL(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IN_FINAL_CONFIG = LD
END SUBROUTINE SET_IN_FINAL
!======================================================================
PURE FUNCTION L_OBS_IN_FC()
LOGICAL :: L_OBS_IN_FC
L_OBS_IN_FC = ALGORITHM_STATE%L_OBS_IN_FC_CONFIG
END FUNCTION L_OBS_IN_FC
!======================================================================
SUBROUTINE SET_OBS_IN_FC(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_OBS_IN_FC_CONFIG = LD
END SUBROUTINE SET_OBS_IN_FC
!======================================================================
SUBROUTINE SET_IN_WRITE_INCREMENT(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IN_WRITE_INCREMENT = LD
END SUBROUTINE SET_IN_WRITE_INCREMENT
!======================================================================
PURE FUNCTION L_IN_WRITE_INCREMENT()
LOGICAL :: L_IN_WRITE_INCREMENT
L_IN_WRITE_INCREMENT = ALGORITHM_STATE%L_IN_WRITE_INCREMENT
END FUNCTION L_IN_WRITE_INCREMENT
!======================================================================
SUBROUTINE SET_IN_WRITE_STATE(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IN_WRITE_STATE = LD
END SUBROUTINE SET_IN_WRITE_STATE
!======================================================================
PURE FUNCTION L_IN_WRITE_STATE()
LOGICAL :: L_IN_WRITE_STATE
L_IN_WRITE_STATE = ALGORITHM_STATE%L_IN_WRITE_STATE
END FUNCTION L_IN_WRITE_STATE
!======================================================================
SUBROUTINE SET_IOSTREAM_FDB(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IOSTREAM_FDB_SET = LD
END SUBROUTINE SET_IOSTREAM_FDB
!======================================================================
PURE FUNCTION L_IOSTREAM_FDB_SET()
LOGICAL :: L_IOSTREAM_FDB_SET
L_IOSTREAM_FDB_SET = ALGORITHM_STATE%L_IOSTREAM_FDB_SET
END FUNCTION L_IOSTREAM_FDB_SET
!======================================================================
SUBROUTINE SET_IN_MINIM(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IN_MINIM = LD
END SUBROUTINE SET_IN_MINIM
!======================================================================
PURE FUNCTION L_IN_MINIM()
LOGICAL :: L_IN_MINIM
L_IN_MINIM = ALGORITHM_STATE%L_IN_MINIM
END FUNCTION L_IN_MINIM
!======================================================================
SUBROUTINE SET_L_TRAJ_COMPLETE(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_TRAJ_COMPLETE = LD
END SUBROUTINE SET_L_TRAJ_COMPLETE
!======================================================================
PURE FUNCTION L_TRAJ_COMPLETE()
LOGICAL :: L_TRAJ_COMPLETE
L_TRAJ_COMPLETE = ALGORITHM_STATE%L_TRAJ_COMPLETE
END FUNCTION L_TRAJ_COMPLETE
!======================================================================
SUBROUTINE SET_L_OUTER_COMPLETE(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_OUTER_COMPLETE = LD
END SUBROUTINE SET_L_OUTER_COMPLETE
!======================================================================
PURE FUNCTION L_OUTER_COMPLETE()
LOGICAL :: L_OUTER_COMPLETE
L_OUTER_COMPLETE = ALGORITHM_STATE%L_OUTER_COMPLETE
END FUNCTION L_OUTER_COMPLETE
!======================================================================
SUBROUTINE SET_L_4V_WRITTEN(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_4V_WRITTEN = LD
END SUBROUTINE SET_L_4V_WRITTEN
!======================================================================
PURE FUNCTION L_4V_WRITTEN()
LOGICAL :: L_4V_WRITTEN
L_4V_WRITTEN  = ALGORITHM_STATE%L_4V_WRITTEN
END FUNCTION L_4V_WRITTEN
!======================================================================
SUBROUTINE SET_L_LATE4D_CHECK_COMPLETE(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_LATE4D_CHECK_COMPLETE = LD
END SUBROUTINE SET_L_LATE4D_CHECK_COMPLETE
!======================================================================
PURE FUNCTION L_LATE4D_CHECK_COMPLETE()
LOGICAL :: L_LATE4D_CHECK_COMPLETE
L_LATE4D_CHECK_COMPLETE  = ALGORITHM_STATE%L_LATE4D_CHECK_COMPLETE
END FUNCTION L_LATE4D_CHECK_COMPLETE
!======================================================================
SUBROUTINE SET_L_IN_LOWRES_STATE_CTOR(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IN_LOWRES_STATE_CTOR = LD
END SUBROUTINE SET_L_IN_LOWRES_STATE_CTOR
!======================================================================
PURE FUNCTION L_IN_LOWRES_STATE_CTOR()
LOGICAL :: L_IN_LOWRES_STATE_CTOR
L_IN_LOWRES_STATE_CTOR = ALGORITHM_STATE%L_IN_LOWRES_STATE_CTOR
END FUNCTION L_IN_LOWRES_STATE_CTOR
!======================================================================
SUBROUTINE SET_L_IN_RESTART_CTLVEC(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IN_RESTART_CTLVEC = LD
END SUBROUTINE SET_L_IN_RESTART_CTLVEC
!======================================================================
PURE FUNCTION L_IN_RESTART_CTLVEC()
LOGICAL :: L_IN_RESTART_CTLVEC
L_IN_RESTART_CTLVEC = ALGORITHM_STATE%L_IN_RESTART_CTLVEC
END FUNCTION L_IN_RESTART_CTLVEC
!======================================================================
SUBROUTINE SET_L_IFSTRAJ_CHKPT_COMPLETE(LD)
LOGICAL,INTENT(IN) :: LD
ALGORITHM_STATE%L_IFSTRAJ_CHKPT_COMPLETE = LD
END SUBROUTINE SET_L_IFSTRAJ_CHKPT_COMPLETE
!======================================================================
PURE FUNCTION L_IFSTRAJ_CHKPT_COMPLETE()
LOGICAL :: L_IFSTRAJ_CHKPT_COMPLETE
L_IFSTRAJ_CHKPT_COMPLETE  = ALGORITHM_STATE%L_IFSTRAJ_CHKPT_COMPLETE
END FUNCTION L_IFSTRAJ_CHKPT_COMPLETE
!======================================================================

END MODULE ALGORITHM_STATE_MOD
