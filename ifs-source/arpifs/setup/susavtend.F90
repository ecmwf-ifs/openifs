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
!      F. Vana  January 2023: more flexibility to phys-dyn interface
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
!USE YOMLUN   , ONLY : NULOUT
USE YOMSLPHY , ONLY : TSLPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KLEVDIM
LOGICAL,INTENT(IN) :: LDSLPHY
TYPE(TSLPHY),INTENT(INOUT):: YDSLPHY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSAVTEND',0,ZHOOK_HANDLE)
ASSOCIATE(MT_SAVTEND=>YDSLPHY%MT_SAVTEND, MQ_SAVTEND=>YDSLPHY%MQ_SAVTEND, &
 & MU_SAVTEND=>YDSLPHY%MU_SAVTEND, MV_SAVTEND=>YDSLPHY%MV_SAVTEND, &
 & MTS_SAVTEND=>YDSLPHY%MTS_SAVTEND, MQS_SAVTEND=>YDSLPHY%MQS_SAVTEND, &
 & MUS_SAVTEND=>YDSLPHY%MUS_SAVTEND, MVS_SAVTEND=>YDSLPHY%MVS_SAVTEND, &
 & MSPPTCLEAR_SAVTEND=>YDSLPHY%MSPPTCLEAR_SAVTEND,&
 & M_RAD_SW=>YDSLPHY%M_RAD_SW    , M_RAD_LW=>YDSLPHY%M_RAD_LW    , &
 & M_CONV_MFD=>YDSLPHY%M_CONV_MFD, M_CONV_MFU=>YDSLPHY%M_CONV_MFU, &
 & M_CONV_LUDELI=>YDSLPHY%M_CONV_LUDELI)
!     ------------------------------------------------------------------

! Set M.._SAVTEND look up table.
! Pointers are set in any case (physics is not being set at this moment)

!Pointers commonly used in SAVTEND and SAVTEND_VD (for different quantities though)
MU_SAVTEND=1   ! SLAVEPP tendency of u-wind
MV_SAVTEND=2   ! SLAVEPP tendency of v-wind
MT_SAVTEND=3   ! SLAVEPP tendency of temperature

! Pointers exclusive to SAVTEND
! following are corrections for standard SLAVEPP tendencies by convection tendency
MUS_SAVTEND=4  ! SLAVEPP tendency of u-wind for simplified physics steps (LPC_SPHY=T)
MVS_SAVTEND=5  ! SLAVEPP tendency of u-wind for simplified physics steps (LPC_SPHY=T)
MTS_SAVTEND=6  ! SLAVEPP tendency of u-wind for simplified physics steps (LPC_SPHY=T)
MQS_SAVTEND=7  ! SLAVEPP tendency of u-wind for simplified physics steps (LPC_SPHY=T)
! Only usefull with SPPT
MSPPTCLEAR_SAVTEND=8 ! Stores clear sky radiation from previous timstep (provisional value, see SUSC2B)

! Pointers exclusive to SAVTEND_VD
! pointers storing tendencies between timesteps, 
MQ_SAVTEND=4   ! vertical diffusion tendency from previous time step (M[U,V,T]_SAVTEND are also used here)
M_RAD_SW=5     ! SW tendency from radiation
M_RAD_LW=6     ! LW tendency from radiation
M_CONV_MFD=7   ! mass flux tendency from convection
M_CONV_MFU=8   ! mass flux tendency from convection
M_CONV_LUDELI(1:4)=(/9,10,11,12/) ! detrained liquid, ice, vapor & T

! Allocate space to handle implicitness
IF( LDSLPHY ) ALLOCATE(YDSLPHY%RSLWX(KLEVDIM))

! Do we really need to print anything here?

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSAVTEND',1,ZHOOK_HANDLE)
END SUBROUTINE SUSAVTEND
