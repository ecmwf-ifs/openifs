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

SUBROUTINE SUTRAJP

!**** *SUTRAJP*   - Set-up parameters for adjoint of physics

!     Purpose.
!     --------

!        Set-up physical parametrizations to be activated in
!        the tangent-linear and adjoints versions of the IFS
!        and options related to the management of the
!        trajectory at t-dt

!**   Interface.
!     ----------
!        *CALL* *SUTRAJP*

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        COMMON YOMIOP, YOPHNC, YOMSIMPHL

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      J.-F. Mahfouf
!      Original : 96-09-27

!     Modifications.
!     --------------
!      1997-1999  M. Janiskova (set-up for Meteo-France adjoint physics)
!      2002-10-30 M. Janiskova (3D fields for cloud parameters)
!      2003-12-10 P. Lopez     (3D fields for convective parameters)
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      2008-05-01 M. Janiskova (updated set-up for lw radiation)
!      2009-09-21 M. Janiskova (3D fields for non-orog.gravity wave)
!      O.Riviere   Nov 09 Increase NG3PR95 from 27 to 28 if LTRAJPST
!                          for LGWDSPNL in simpl physics.
!      2010-05-18 M. Janiskova (updated set-up when lemwave)
!      O.Riviere Feb 11 Additional storage of QL/QI and tend. for MF's setting
!      2012-05-03 M. Janiskova (updated set-up when lesurf2)
!      F. Vana  3-Sep-2020 Pruning
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMIOPNH , ONLY : NG3NH95

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUTRAJP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       2.    Set-default values for the number of 3D/2D fields to store
!              ----------------------------------------------------------

! NOTE: This routine represents the original trajectory handling being used until
!       around 2012. Since then all the useful trajectory code was converted
!       to new style defined in YOMTRAJ.
!       The only code still theoretically relying on this old fashioned style
!       trajectory is the NH dynamics TL code from Meteo-France being for some
!       reason implemented into the trajectory code of physics. Being broken
!       since long ago this NH TL code has its relevant trajectory parts 
!       commented out at the places of potential use (i.e. CPG_GP or CPG5_GP).
!       Logically there is only little reason to keep this totally obsolete setup.


!* ECMWF set-up
!  ------------

!   Not useful for any configuration operated at ECMWF


!* Meteo-France set-up
!  -------------------

NG3NH95 = 9  ! never used anywhere

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUTRAJP',1,ZHOOK_HANDLE)
END SUBROUTINE SUTRAJP
