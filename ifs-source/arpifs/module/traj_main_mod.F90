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

MODULE TRAJ_MAIN_MOD

!     Purpose.
!     --------
!       Manage main trajectory

!     Author.
!     -------
!        Y. Tremolet *ECMWF*
!        Original : 09-10-01

!     Modifications.
!     --------------
!        Modified : 02-09-16 G. Desroziers - also treat Aladin mean wind trajectory
!        O.Spaniel    : 03-04-15 cleaning-ordering of declaration
!        Modified : 03-03-31 C. Fischer - readapted Aladin mean wind to control_vectors
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        G. Desroziers & C. Fischer: 04-06-08 Secure gridpoint copy to/from TRAJ_grids
!        Y.Tremolet    18-Mar-2004 Add lprttraj option for diagnostics
!        J.Haseler     11-Oct-2005 Generalise gridpoint trajectory contents
!        M.Hamrud  15-Jan-2006  Revised GPRCP
!        G. Desroziers and K. Yessad (sept 2005):
!         - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!         - adapt option LTRAJHR to METEO-FRANCE configurations.
!        J.Haseler     27-Feb-2007 Simplify logic
!        Y.Tremolet    27-Nov-2008 Long window 4D-Var
!        N Wedi (Nov 2011) : bugfix
!        F. Vana  28-Nov-2013 : Redesigned trajectory handling
!        K. Yessad (July 2014): Move some variables.
!        M. Fisher  9 July 2015: Disentangle background from trajectory
!        A. Trojakova (May 2018): Fix for ALARO 3DVAR minimization using
!                                 spectral Q in the background field input.
!        E. Holm    27-Mar-2019: Write out BG different from FG if L_FGNOTBG
!        Y. Michel (June 2018) : Second fix for LAM when LSPRT=.TRUE.
!        E. Holm    03-Dec-2019: Specialised STORE_BACKGROUND for writing BG only,
!                                for re-centring background (L_BGREC=T in RDFPINC).
!                                STORE_MAIN_TRAJ: L_BGREC=T => don't write BG
!     ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN,   ONLY : NULOUT, NULERR
USE YOMCT0,   ONLY : LIFSTRAJ, NINTERPTRAJ, LARPEGEF_TRAJHR, LARPEGEF_TRAJBG, LELAM
USE YOMVAR,   ONLY : MUPTRA, FILTERFACTOR, FILTEREXPO, FILTERRESOL, L_FGNOTBG, L_BGREC
USE ALGORITHM_STATE_MOD  , ONLY : GET_NUPTRA
USE YOMMODERR,ONLY : N_COUPLED_WINDOWS
USE YOMDIM,   ONLY : TDIM
USE YOMMP0,   ONLY : MYSETV
USE YOM_GRIB_CODES,  ONLY : NGRBVO, NGRBD, NGRBT, NGRBQ, NGRBO3, NGRBLNSP
USE YOMTRAJ,  ONLY : TRAJEC, NSMAX_TRAJ,&
 &                   LPRTTRAJ, LTRAJHR_ALTI,&
 &                   LTRAJHR, LTRAJGP, LREADGPTRAJ, NSTEPTRAJ, MSTART, NSPTRAJ, NGPTRAJ,&
 &                   MSTEPTRAJW, MSTEPTRAJR, MIOTRAJMAIN, TRAJ_SPEC_TMP0, LTRAJRESET,&
 &                   TRAJ_TYPE, MTYPE_MAIN3_TRAJ, MTYPE_MAIN2_TRAJ, MAIN_GRIB, MBACKGR_GRIB
USE YOMGMV,   ONLY : TGMV
USE YOMCST,   ONLY : RD
USE YOMOPH0,  ONLY : CFNTRAJHRGRID, CFNTRAJBGGRID, CFNTRAJHRSPEC, CFNTRAJBGSPEC
USE SPECTRAL_FIELDS_MOD

IMPLICIT NONE

PRIVATE
!OIFS - ifdef included to remove trajectory code from OpenIFS, since code
!       should not be distributed. #else at end of module replaces code
!       with interfaced dummies
PUBLIC STORE_MAIN_TRAJ, GET_TRAJ_SPEC, GET_TRAJ_GRID

CONTAINS

SUBROUTINE STORE_MAIN_TRAJ(YDGEOMETRY,YDGMV,YDGMV5,YDMODEL,YDSPEC,PGMV,PGMVS,PGFL,KSTEP,LDREADGPTRAJ,PTRAJEC,YD_BG)

USE TYPE_MODEL             , ONLY : MODEL
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE FIELDS_MOD             , ONLY : FIELDS, FIELDS_CREATE, FIELDS_DELETE
USE VARIABLES_MOD          , ONLY : VARIABLES, VARIABLES_CREATE, VARIABLES_DELETE

TYPE(GEOMETRY)      ,INTENT(INOUT) :: YDGEOMETRY
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV5
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
TYPE(SPECTRAL_FIELD),INTENT(IN)    :: YDSPEC
REAL(KIND=JPRB)     ,INTENT(IN)    :: PGMV(:,:,:,:)
REAL(KIND=JPRB)     ,INTENT(IN)    :: PGMVS(:,:,:)
REAL(KIND=JPRB)     ,INTENT(IN)    :: PGFL(:,:,:,:)
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KSTEP
LOGICAL             ,INTENT(IN)    :: LDREADGPTRAJ
TYPE(TRAJ_TYPE)     ,INTENT(INOUT) :: PTRAJEC
TYPE(TRAJ_TYPE)     ,INTENT(INOUT) :: YD_BG

call abor1('OIFS - store_main_traj should never be called')

END SUBROUTINE STORE_MAIN_TRAJ

SUBROUTINE GET_TRAJ_SPEC(YDGEOMETRY,YDDIMF,YDTRAJEC,YDGMV,YDGMV5,YDSPEC,KSTEP,LDBACKGR)

USE YOMDIMF     , ONLY : TDIMF
USE GEOMETRY_MOD, ONLY : GEOMETRY

TYPE(GEOMETRY)      , INTENT(IN)    :: YDGEOMETRY
TYPE(TDIMF)         , INTENT(IN)    :: YDDIMF
TYPE(TRAJ_TYPE)     , INTENT(IN)    :: YDTRAJEC
TYPE(TGMV)          , INTENT(INOUT) :: YDGMV
TYPE(TGMV)          , INTENT(INOUT) :: YDGMV5
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSPEC
INTEGER(KIND=JPIM)  , INTENT(IN)    :: KSTEP
LOGICAL, OPTIONAL   , INTENT(IN)    :: LDBACKGR

call abor1('OIFS - get_traj_spec should never be called')

END SUBROUTINE GET_TRAJ_SPEC

SUBROUTINE GET_TRAJ_GRID(YDGEOMETRY,YDML_GCONF,YDTRAJEC,YDGMV,YDGMV5,PGMV,PGMVS,PGFL,KUPTRA)

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY

TYPE(GEOMETRY)  , INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(TRAJ_TYPE) , INTENT(IN)    :: YDTRAJEC
TYPE(TGMV)      , INTENT(INOUT) :: YDGMV
TYPE(TGMV)      , INTENT(INOUT) :: YDGMV5
REAL(KIND=JPRB) , INTENT(OUT)   :: PGMV(:,:,:,:)
REAL(KIND=JPRB) , INTENT(OUT)   :: PGMVS(:,:,:)
REAL(KIND=JPRB) , INTENT(OUT)   :: PGFL(:,:,:,:)
INTEGER(KIND=JPIM), INTENT(IN)    :: KUPTRA

call abor1('OIFS - get_traj_grid should never be called')

END SUBROUTINE GET_TRAJ_GRID



END MODULE TRAJ_MAIN_MOD
