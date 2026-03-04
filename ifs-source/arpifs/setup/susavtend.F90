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

SUBROUTINE SUSAVTEND(KLEVDIM,LDSLPHY,YDSLPHY)

!**** *SUSAVTEND*   - Initialize pointers for SAVTEND

!     Purpose.
!     --------
!           Initialize pointers for SAVTEND

!**   Interface.
!     ----------
!        *CALL* *SUSAVTEND()

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Called by SUDYN.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!         K.Yessad
!         Original : 10-Feb-2006

!     Modifications
!     -------------
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      F. Vana  11-Nov-2014: Fix for Eulerian advection
!      F. Vana  23-Oct-2018: Cleaning
!      F. Vana  11-Sep-2020: Vertical flexibility in implicitness
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : NUNDEFLD
USE YOMLUN   , ONLY : NULOUT
USE YOMSLPHY , ONLY : TSLPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KLEVDIM
LOGICAL,INTENT(IN) :: LDSLPHY
TYPE(TSLPHY),INTENT(INOUT):: YDSLPHY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSAVTEND',0,ZHOOK_HANDLE)
ASSOCIATE(MT_SAVTEND=>YDSLPHY%MT_SAVTEND, &
 & MU_SAVTEND=>YDSLPHY%MU_SAVTEND, MV_SAVTEND=>YDSLPHY%MV_SAVTEND, &
 & NVTEND=>YDSLPHY%NVTEND)
!     ------------------------------------------------------------------

! Calculation of NVTEND and M.._SAVTEND.

IF( LDSLPHY ) THEN
  NVTEND=3
  MU_SAVTEND=1
  MV_SAVTEND=2
  MT_SAVTEND=3
  ALLOCATE(YDSLPHY%RSLWX(KLEVDIM))
ELSE
  NVTEND=0 ! set to 1 when causing a memory issue 
  MU_SAVTEND=NUNDEFLD
  MV_SAVTEND=NUNDEFLD
  MT_SAVTEND=NUNDEFLD
ENDIF

! Printings.

WRITE(NULOUT,*) ' -- YOMSLPHY CONTENT: '
WRITE(NULOUT,*) ' NVTEND= ',NVTEND

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSAVTEND',1,ZHOOK_HANDLE)
END SUBROUTINE SUSAVTEND
