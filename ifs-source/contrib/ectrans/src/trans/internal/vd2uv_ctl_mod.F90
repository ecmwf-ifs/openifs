! (C) Copyright 2015- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VD2UV_CTL_MOD
CONTAINS
SUBROUTINE VD2UV_CTL(KF_UV,PSPVOR,PSPDIV,PU,PV)

!**** *VD2UV_CTL* - Control routine for going from vor/div to spectral U and V.

!     Purpose.
!     --------
!        Control routine for computing spectral U (u*cos(theta)) and V 

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTL(...)
!     KF_UV        - local number of spectral u-v fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PU(:,:)      - U (out)
!     PV(:,:)      - V (out)

!     Method.
!     -------

!     Externals.
!     ----------


!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : July 2015

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR       ,ONLY : D

USE VD2UV_MOD       ,ONLY : VD2UV

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV
REAL(KIND=JPRB),INTENT(IN)    :: PSPVOR(:,:)
REAL(KIND=JPRB),INTENT(IN)    :: PSPDIV(:,:)
REAL(KIND=JPRB),INTENT(OUT)   :: PU(:,:)
REAL(KIND=JPRB),INTENT(OUT)   :: PV(:,:)

INTEGER(KIND=JPIM) :: JM,IM,ILEI2

!     ------------------------------------------------------------------

CALL GSTATS(102,0)
ILEI2 = 8*KF_UV

CALL GSTATS(1647,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)
  CALL VD2UV(IM,JM,KF_UV,ILEI2,PSPVOR,PSPDIV,PU,PV)
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1647,1)
CALL GSTATS(102,1)

!     ------------------------------------------------------------------

END SUBROUTINE VD2UV_CTL
END MODULE VD2UV_CTL_MOD
