! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWNEMOFLDS

      USE YOWDRVTYPE  , ONLY : WAVE2OCEAN, OCEAN2WAVE

      IMPLICIT NONE

!*     ** *NEMOFLDS* NEMO FIELDS FOR COUPLED RUNS

      TYPE(WAVE2OCEAN) :: WAM2NEMO
      TYPE(OCEAN2WAVE) :: NEMO2WAM

      LOGICAL :: LNEMOCITHICK, LNEMOICEREST

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     LNEMOCITHICK LOGICAL   SET TO TRUE IF SEA ICE THICKNESS IS 
!                            AVAILABLE FROM NEMO (E.G. LIM2 ACTIVE).
!     LNEMOICEREST LOGICAL   SET TO TRUE IF SEA ICE IS NOT RESCALED BY ICE COVER
!---------------------------------------------------------------------
      END MODULE YOWNEMOFLDS
