! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GPCTY1C(YDDIMV,KFLEV,PVVEL,PVVELH)

!**** *GPCTY* - Computes vertical velocities.

!     Purpose.
!     --------
!           COMPUTES VERTICAL VELOCITIES AT HALF LEVELS.

!**   Interface.
!     ----------
!        *CALL* *GPCTY1C(...)

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          KFLEV            : number of layers
!          PVVEL(KFLEV)     : (omega/pressure)

!        * OUTPUT:
!          PVVELH(0:KFLEV)  : (ETA DOT (D P)/(D ETA))

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the single column model

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original    94-05-06
!        M. Ko"hler  6-6-2006 Single Column Model integration within IFS 

!     ------------------------------------------------------------------


USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
INTEGER(KIND=JPIM) :: KFLEV

REAL(KIND=JPRB) :: PVVEL(YDDIMV%NFLEVG) , PVVELH(0:YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: PADV(YDDIMV%NFLEVG)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JLEV


!     ------------------------------------------------------------------

!*       1.    COMPUTES W AT HALF LEVELS.
!              --------------------------


PVVELH(0)=0.0_JPRB
PVVELH(KFLEV)=0.0_JPRB

DO JLEV=1,KFLEV-1

  PVVELH(JLEV)=0.5_JPRB*(PVVEL(JLEV)+PVVEL(JLEV+1))

! omega - Vh grad (p)  IMPORTANT in the presence of orography!
! (this option has now been replaced with supplying etadotdpdeta on half levels)
  
! PVVELH(JLEV)=PVVELH(JLEV)+0.5_JPRB*(PADV(JLEV)+PADV(JLEV+1))

ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE GPCTY1C
