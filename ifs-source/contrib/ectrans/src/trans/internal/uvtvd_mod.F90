! (C) Copyright 1991- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UVTVD_MOD
CONTAINS
SUBROUTINE UVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!**** *UVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX.

!**   Interface.
!     ----------
!        CALL UVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!        Explicit arguments :  KM - zonal wave-number
!        --------------------  KFIELD - number of fields (levels)
!                              PEPSNM - REPSNM for wavenumber KM
!                              PU - u wind component for zonal
!                                   wavenumber KM
!                              PV - v wind component for zonal
!                                   wavenumber KM
!                              PVOR - vorticity for zonal
!                                     wavenumber KM
!                              PDIV - divergence for zonal
!                                     wavenumber KM


!     Method.  See ref.
!     -------

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 91-07-01
!        D. Giard : NTMAX instead of NSMAX
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
!USE TPM_DISTR
!

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELD
INTEGER(KIND=JPIM), INTENT(IN)  :: KM

REAL(KIND=JPRB), INTENT(IN)     :: PEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB), INTENT(OUT)    :: PVOR(:,:),PDIV(:,:)
REAL(KIND=JPRB), INTENT(INOUT)  :: PU  (:,:),PV  (:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IN, IR, J, JN, ITMAX

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZN(-1:R%NTMAX+3)


!     ------------------------------------------------------------------


!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ZKM = KM
ITMAX = R%NTMAX
ZN(KM-1:ITMAX+3) = F%RN(KM-1:ITMAX+3)
!*       1.1      SET N=KM-1 COMPONENT TO 0 FOR U AND V

IN = F%NLTN(KM-1)
DO J=1,2*KFIELD
  PU(IN,J) = 0.0_JPRB
  PV(IN,J) = 0.0_JPRB
ENDDO

!*       1.2      COMPUTE VORTICITY AND DIVERGENCE.

IF(KM /= 0) THEN
  DO JN=KM,ITMAX
    IN = ITMAX+2-JN
!DIR$ IVDEP
!OCL NOVREC
    DO J=1,KFIELD
      IR = 2*J-1
      II = IR+1
      PVOR(IN,IR) = -ZKM*PV(IN,II)-&
       &ZN(JN)*PEPSNM(JN+1)*PU(IN-1,IR)+&
       &ZN(JN+1)*PEPSNM(JN)*PU(IN+1,IR)
      PVOR(IN,II) = +ZKM*PV(IN,IR)-&
       &ZN(JN)*PEPSNM(JN+1)*PU(IN-1,II)+&
       &ZN(JN+1)*PEPSNM(JN)*PU(IN+1,II)
      PDIV(IN,IR) = -ZKM*PU(IN,II)+&
       &ZN(JN)*PEPSNM(JN+1)*PV(IN-1,IR)-&
       &ZN(JN+1)*PEPSNM(JN)*PV(IN+1,IR)
      PDIV(IN,II) = +ZKM*PU(IN,IR)+&
       &ZN(JN)*PEPSNM(JN+1)*PV(IN-1,II)-&
       &ZN(JN+1)*PEPSNM(JN)*PV(IN+1,II)
    ENDDO
  ENDDO
ELSE
  DO JN=KM,ITMAX
    IN = ITMAX+2-JN
    DO J=1,KFIELD
      IR = 2*J-1
      PVOR(IN,IR) = -&
       &ZN(JN)*PEPSNM(JN+1)*PU(IN-1,IR)+&
       &ZN(JN+1)*PEPSNM(JN)*PU(IN+1,IR)
      PDIV(IN,IR) = &
       &ZN(JN)*PEPSNM(JN+1)*PV(IN-1,IR)-&
       &ZN(JN+1)*PEPSNM(JN)*PV(IN+1,IR)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UVTVD
END MODULE UVTVD_MOD
