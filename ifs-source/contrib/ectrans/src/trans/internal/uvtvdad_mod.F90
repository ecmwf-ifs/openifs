! (C) Copyright 1991- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UVTVDAD_MOD
CONTAINS
SUBROUTINE UVTVDAD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!**** *UVTVDAD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX.

!**   Interface.
!     ----------
!        CALL UVTVDAD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

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
REAL(KIND=JPRB), INTENT(IN)     :: PVOR(:,:),PDIV(:,:)
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

!*       1.2      COMPUTE VORTICITY AND DIVERGENCE.

IF(KM /= 0) THEN
  DO JN=KM,ITMAX
    IN = ITMAX+2-JN
!DIR$ IVDEP
!OCL NOVREC
    DO J=1,KFIELD
      IR = 2*J-1
      II = IR+1

      PV(IN,II)   = PV(IN,II)-ZKM*PVOR(IN,IR)
      PU(IN-1,IR) = PU(IN-1,IR)-ZN(JN)*PEPSNM(JN+1)*PVOR(IN,IR)
      PU(IN+1,IR) = PU(IN+1,IR)+ZN(JN+1)*PEPSNM(JN)*PVOR(IN,IR)

      PV(IN,IR)   = PV(IN,IR)+ZKM*PVOR(IN,II)
      PU(IN-1,II) = PU(IN-1,II)-ZN(JN)*PEPSNM(JN+1)*PVOR(IN,II)
      PU(IN+1,II) = PU(IN+1,II)+ZN(JN+1)*PEPSNM(JN)*PVOR(IN,II)

      PU(IN,II)   = PU(IN,II)-ZKM*PDIV(IN,IR)
      PV(IN-1,IR) = PV(IN-1,IR)+ZN(JN)*PEPSNM(JN+1)*PDIV(IN,IR)
      PV(IN+1,IR) = PV(IN+1,IR)-ZN(JN+1)*PEPSNM(JN)*PDIV(IN,IR)

      PU(IN,IR)   = PU(IN,IR)+ZKM*PDIV(IN,II)
      PV(IN-1,II) = PV(IN-1,II)+ZN(JN)*PEPSNM(JN+1)*PDIV(IN,II)
      PV(IN+1,II) = PV(IN+1,II)-ZN(JN+1)*PEPSNM(JN)*PDIV(IN,II)
    ENDDO
  ENDDO
ELSE
  DO JN=KM,ITMAX
    IN = ITMAX+2-JN
    DO J=1,KFIELD
      IR = 2*J-1
      PU(IN-1,IR) = PU(IN-1,IR)-ZN(JN  )*PEPSNM(JN+1)*PVOR(IN,IR)
      PU(IN+1,IR) = PU(IN+1,IR)+ZN(JN+1)*PEPSNM(JN  )*PVOR(IN,IR)
      PV(IN-1,IR) = PV(IN-1,IR)+ZN(JN  )*PEPSNM(JN+1)*PDIV(IN,IR)
      PV(IN+1,IR) = PV(IN+1,IR)-ZN(JN+1)*PEPSNM(JN  )*PDIV(IN,IR)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UVTVDAD
END MODULE UVTVDAD_MOD
