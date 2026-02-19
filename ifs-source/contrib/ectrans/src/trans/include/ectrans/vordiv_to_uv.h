! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE VORDIV_TO_UV(PSPVOR,PSPDIV,PSPU,PSPV,KSMAX,KVSETUV)

!**** *VORDIV_TO_UV* - Convert spectral vorticity and divergence to spectral U (u*cos(theta)) and V (v*cos(theta).

!     Purpose.
!     --------
!        Interface routine for  Convert spectral vorticity and divergence to spectral U  and V 

!**   Interface.
!     ----------
!     CALL VORDIV_TO_UV(...)

!     Explicit arguments :
!     --------------------
!     PSPVOR(:,:) - spectral vorticity (input)
!     PSPDIV(:,:) - spectral divergence (input)
!     PSPU(:,:)   - spectral U (u*cos(theta) (output)
!     PSPV(:,:)   - spectral V (v*cos(theta) (output)
!     KSMAX       - spectral resolution (input)
!     KVSETUV(:)  - Optionally indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.

!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  VD2UV_CTL   - control vordiv to uv

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 15-06-15


!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB), INTENT(IN) :: PSPVOR(:,:)
REAL(KIND=JPRB), INTENT(IN) :: PSPDIV(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPU(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPV(:,:)
INTEGER(KIND=JPIM) , INTENT(IN) :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)

END SUBROUTINE VORDIV_TO_UV
END INTERFACE