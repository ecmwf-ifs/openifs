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

MODULE TRAJECTORY_MOD

!     Purpose.
!     --------
!       Manage trajectory for TL and AD runs in a unified way.

!     Author.
!     -------
!        Y. Tremolet *ECMWF*
!        Original : 2000-08-04

!     Modifications.
!     --------------
!        Modified : 02-05-06 C. Fischer - ltrajsave=false when deallocate
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modified : 03-12-19 C. Temperton - case LSLAG=.F. but LSPRT=.T.
!        Y.Tremolet    18-Mar-2004 Add lprttraj option for diagnostics
!        L.Isaksen: 2004-11-11 Upperair grid point trajectory option added
!        G. Desroziers and K. Yessad (sept 2005):
!         - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!         - adapt option LTRAJHR to METEO-FRANCE configurations.
!        M.Hamrud      01-Jul-2006  Revised surface fields
!        F. Vana  28-Nov-2013 : Redesigned trajectory handling
!        K. Yessad (July 2014): Move some variables.
!        E. Holm  22-Jul-2016 : Hacked BACKGR04-10, do properly!
!     ------------------------------------------------------------------

USE PARKIND1          , ONLY : JPIM, JPRB
USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0            , ONLY : LIFSTRAJ, LIFSMIN, LARPEGEF_TRAJHR, LARPEGEF_TRAJBG, LECMWF
USE YOMDIM            , ONLY : TDIM
USE YOMVAR            , ONLY : MUPTRA
USE ALGORITHM_STATE_MOD  , ONLY : GET_NUPTRA
USE YOMLUN            , ONLY : NULOUT
USE SURFACE_FIELDS_MIX, ONLY : NSTRAJGRIB
USE YOM_GRIB_CODES           , ONLY : NGRBLNSP
USE IOSTREAM_MIX      , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST,&
 &                             CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST, IO_INQUIRE,&
 &                             CLOSE_IOREQUEST
USE YOMTRAJ           , ONLY : LTRAJSAVE, LTRAJSLAG, LTRAJPHYS, LTRAJCST,&
 &                             LTRAJGP, LTRAJHR, LREADGPTRAJ, TSTEP_TRAJ, NSMAX_TRAJ,&
 &                             NSMAX_BACKGR00, NSMAX_BACKGR01, NSMAX_BACKGR02, NSMAX_BACKGR03,&
 &                             NSMAX_BACKGR04, NSMAX_BACKGR05, NSMAX_BACKGR06, NSMAX_BACKGR07,&
 &                             NSMAX_BACKGR08, NSMAX_BACKGR09, NSMAX_BACKGR10,&
 &                             LTRAJHR_ALTI,LTRAJHR_SURF, LSNOWTRAJCONS, &
 &                             NGP5, NTRAJ_CST, LTRAJALLOC, MSTART,&
 &                             NSTEPTRAJ, MSTEPTRAJW, MSTEPTRAJR, LPRTTRAJ, LTRAJRESET,&
 &                             MTYPE_SURF_TRAJ, MTRAJ_GRIB,MTRAJ_GRIB2, MTYPE_MAIN3_TRAJ, MTYPE_MAIN2_TRAJ,&
 &                             MAIN_GRIB, MBACKGR_GRIB, TRAJEC, TRAJ_TYPE
USE TRAJ_MAIN_MOD






IMPLICIT NONE

PRIVATE
PUBLIC LTRAJSAVE, LTRAJSLAG, LTRAJPHYS, LTRAJCST,&
 & LTRAJGP, LTRAJHR, LREADGPTRAJ, TSTEP_TRAJ, NSMAX_TRAJ,&
 & LTRAJHR_ALTI,LTRAJHR_SURF, LSNOWTRAJCONS, &
 & NSMAX_BACKGR00, NSMAX_BACKGR01, NSMAX_BACKGR02, NSMAX_BACKGR03,&
 & NSMAX_BACKGR04, NSMAX_BACKGR05, NSMAX_BACKGR06, NSMAX_BACKGR07,&
 & NSMAX_BACKGR08, NSMAX_BACKGR09, NSMAX_BACKGR10,&
 & ALLOCATE_TRAJECTORY, DEALLOCATE_TRAJECTORY,&
 & READ_TRAJECTORY, READ_BACKGROUND,&
 & STORE_MAIN_TRAJ, GET_TRAJ_SPEC, GET_TRAJ_GRID,&
 & NGP5, NTRAJ_CST, LPRTTRAJ, LTRAJRESET

#include "abor1.intfb.h"





!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE ALLOCATE_TRAJECTORY(YDGEOMETRY,YDGMV,YDGMV5,YDSURF,YDMODEL)
USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGMV             , ONLY : TGMV

TYPE(GEOMETRY),INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)    ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)    ,INTENT(INOUT) :: YDGMV5
TYPE(TSURF)   ,INTENT(INOUT) :: YDSURF
TYPE(MODEL)   ,INTENT(IN)    :: YDMODEL
call abor1("OIFS - allocate_trajectory should never be called")

END SUBROUTINE ALLOCATE_TRAJECTORY

!-----------------------------------------------------------------------

SUBROUTINE DEALLOCATE_TRAJECTORY(YDDIM,YDGMV,YDGMV5,YDSURF,YDMODEL)
USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGMV             , ONLY : TGMV

TYPE(TDIM),  INTENT(IN)    :: YDDIM
TYPE(TGMV),  INTENT(INOUT) :: YDGMV
TYPE(TGMV),  INTENT(INOUT) :: YDGMV5
TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(IN)    :: YDMODEL
call abor1('OIFS - deallocate_trajectory should never be called')

END SUBROUTINE DEALLOCATE_TRAJECTORY

!-----------------------------------------------------------------------
SUBROUTINE READ_TRAJECTORY(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,KSTEP,LD_LASTRAJ,LDREADGPTRAJ)

USE TYPE_MODEL  , ONLY : MODEL
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE FIELDS_MOD  , ONLY : FIELDS
USE MTRAJ_MOD   , ONLY : MTRAJ

TYPE(GEOMETRY)    ,INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS)      ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)       ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
LOGICAL           ,INTENT(IN)    :: LD_LASTRAJ,LDREADGPTRAJ
call abor1('OIFS - read_trajectory should never be called')

END SUBROUTINE READ_TRAJECTORY

!-----------------------------------------------------------------------

SUBROUTINE READ_BACKGROUND(YDGEOMETRY,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,YD_BCKG,YDFIELDS,YDMTRAJ,KSTEP,LD_LASTRAJ,LDREADGPTRAJ)

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE YEMLBC_MODEL   , ONLY : TELBC_MODEL
USE FIELDS_MOD   , ONLY : FIELDS
USE MTRAJ_MOD    , ONLY : MTRAJ

TYPE(GEOMETRY)    , INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN) :: YDML_GCONF
TYPE(TDYN)        , INTENT(IN)    :: YDDYN
TYPE(TDYNA)       , INTENT(IN)    :: YDDYNA
TYPE(TELBC_MODEL)     , INTENT(IN)    :: YDML_LBC
TYPE(TRAJ_TYPE)   , INTENT(INOUT) :: YD_BCKG
TYPE(FIELDS)      , INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)       , INTENT(INOUT) :: YDMTRAJ
INTEGER(KIND=JPIM), INTENT(IN)    :: KSTEP
LOGICAL           , INTENT(IN)    :: LD_LASTRAJ,LDREADGPTRAJ
  call abor1('OIFS - read_background should never be called')

END SUBROUTINE READ_BACKGROUND

!-----------------------------------------------------------------------
END MODULE TRAJECTORY_MOD
