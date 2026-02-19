! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VDTUV_MOD
CONTAINS
SUBROUTINE VDTUV(KM,KFIELD,PEPSNM,PVOR,PDIV,PU,PV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F


!**** *VDTUV* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL VDTUV(...)

!        Explicit arguments :  KM -zonal wavenumber (input-c)
!        --------------------  KFIELD - number of fields (input-c)
!                              PEPSNM - REPSNM for wavenumber KM (input-c)
!                              PVOR(NLEI1,2*KFIELD) - vorticity (input)
!                              PDIV(NLEI1,2*KFIELD) - divergence (input)
!                              PU(NLEI1,2*KFIELD)   - u wind (output)
!                              PV(NLEI1,2*KFIELD)   - v wind (output)
!        Organisation within NLEI1:
!        NLEI1 = NSMAX+4+mod(NSMAX+4+1,2)
!                        overdimensioning
!        1        : n=NSMAX+2
!        2        : n=NSMAX+1
!        3        : n=NSMAX
!        .        :
!        .        :
!        NSMAX+3  : n=0
!        NSMAX+4  : n=-1

!        Implicit arguments :  Eigenvalues of inverse Laplace operator
!        --------------------  from YOMLAP

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From VDTUV in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM,KFIELD
REAL(KIND=JPRB), INTENT(IN)    :: PEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB), INTENT(IN)    :: PVOR(:,:),PDIV(:,:)
REAL(KIND=JPRB), INTENT(OUT)   :: PU  (:,:),PV  (:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, ISMAX,JI

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZN(-1:R%NTMAX+4)
REAL(KIND=JPRB) :: ZLAPIN(-1:R%NSMAX+4)
REAL(KIND=JPRB) :: ZEPSNM(-1:R%NSMAX+4)



!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ZKM = KM
ISMAX = R%NSMAX
DO JN=KM-1,ISMAX+2
  IJ = ISMAX+3-JN
  ZN(IJ) = F%RN(JN)
  ZLAPIN(IJ) = F%RLAPIN(JN)
  IF( JN >= 0 ) ZEPSNM(IJ) = PEPSNM(JN)
ENDDO
ZN(0) = F%RN(ISMAX+3)

!*       1.1      U AND V (KM=0) .

IF(KM == 0) THEN
  DO J=1,KFIELD
    IR = 2*J-1
    DO JI=2,ISMAX+3-KM
      PU(JI,IR) = +&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(JI+1,IR)-&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(JI-1,IR)
      PV(JI,IR) = -&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(JI+1,IR)+&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(JI-1,IR)
    ENDDO
  ENDDO

!*       1.2      U AND V (KM!=0) .

ELSE
  DO J=1,KFIELD
    IR = 2*J-1
    II = IR+1
    DO JI=2,ISMAX+3-KM
      PU(JI,IR) = -ZKM*ZLAPIN(JI)*PDIV(JI,II)+&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(JI+1,IR)-&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(JI-1,IR)
      PU(JI,II) = +ZKM*ZLAPIN(JI)*PDIV(JI,IR)+&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(JI+1,II)-&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(JI-1,II)
      PV(JI,IR) = -ZKM*ZLAPIN(JI)*PVOR(JI,II)-&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(JI+1,IR)+&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(JI-1,IR)
      PV(JI,II) = +ZKM*ZLAPIN(JI)*PVOR(JI,IR)-&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(JI+1,II)+&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(JI-1,II)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE VDTUV
END MODULE VDTUV_MOD

