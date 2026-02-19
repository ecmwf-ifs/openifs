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

SUBROUTINE SU3YOM

!**** *SU3YOM*  - INITIALIZE LEVEL 3 COMMON

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SU3YOM

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

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
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!        K. Yessad (Jan 2010): remove useless variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCT3   , ONLY : NSTEP

!      -----------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU3YOM',0,ZHOOK_HANDLE)
!      -----------------------------------------------------------

!*       1.    Initialize LEVEL 3 COMMONS.
!              ---------------------------

NSTEP=-999

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU3YOM',1,ZHOOK_HANDLE)
END SUBROUTINE SU3YOM
