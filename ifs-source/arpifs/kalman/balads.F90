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

SUBROUTINE BALADS(YDLAP,LDADD,PNOR,KFIELD,KFLSUR,KUMP,KSMAX,KASM0,KSPEC2,PVOR,PZ)

!**** *BALADS*  - Selfadjoint linear balance equation

!     Purpose.
!     --------
!           vorticity to equivalent perturbation height conversion

!**   Interface.
!     ----------
!        *CALL* *BALADS(...)

!     Explicit arguments :
!     --------------------
!      Input:
!        LDADD switch : true if operator is additive,
!                       false if PZ is overwritten
!        PNOR normalization factor of PZ (1 for Phi, 1/RG for height)
!        KFIELD number of fields to process
!        KFLSUR array dimension
!        KUMP number of spectral zonal waves handled by this PE
!        KSMAX truncation n of spectral arrays
!        KASM0(0:KSMAX) index of spectral array coeff (n,m) for each m
!        KSPEC2 dimension of spectral arrays=(ksmax+1)*(ksmax+2)
!        PVOR(KFLSUR,KSPEC2) input vorticity fields
!        PZ(KFLSUR,KSPEC2) equivalent height perturbations in geostrophic
!                      framework including beta effect

!     Implicit arguments :
!     --------------------
!        None

!     Method.
!     -------
!       Solve linear balance equation in spectral space
!       to convert vorticity into geopotential.
!       WARNING: the mean height PZ(KASM0(0)) is set to zero regardless
!          of LDADD. One should set it to 5560m for z500 field.
!          (In the adjoint, the mean vorticity gradient is always zero.)
!       Method. See Rochas,la Meteorologie, V-17-1971,p.32
!       solve lapin.div(f.grad.lapin.dzeta)=phi
!       using normalized spectral functions, expressions for
!       lapin,muYmn,dYmn/dmu,f in spectral coordinates,
!       spheric div-> 2romega.mu.dzeta+2romega(1-mu2)dpsi/dmu=lapphi
!       Nota Bene : the discretisation is such that the operator is
!        self adjoint and invertible (though it is ill-conditioned)

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        Rochas,la Meteorologie, V-17-1971,p.32 (in French).
!        METEO-FRANCE/GMAP internal documentation

!     Author.
!     -------
!      Francois Bouttier *Meteo-France*
!      Original : 92-03-27 Francois Bouttier

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE YOMLAP   , ONLY : TLAP
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(TLAP)        ,INTENT(IN)    :: YDLAP
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSUR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2 
LOGICAL           ,INTENT(IN)    :: LDADD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNOR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KUMP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KASM0(0:KSMAX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOR(KFLSUR,KSPEC2) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ(KFLSUR,KSPEC2) 
REAL(KIND=JPRB) :: ZVOR(2,0:KSMAX),ZZ(2,0:KSMAX)

INTEGER(KIND=JPIM) :: IM, JF, JI, JMLOC, JN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "baladsm.intfb.h"

!     -----------------------------------------------------------------

!*    1. SOLVE LINEAR BALANCE EQUATION.
!        ------------------------------

IF (LHOOK) CALL DR_HOOK('BALADS',0,ZHOOK_HANDLE)
ASSOCIATE(MYMS=>YDLAP%MYMS)
DO JF=1,KFIELD
  DO JMLOC=1,KUMP
    IM=MYMS(JMLOC)

!       Interface to input array

    DO JI=1,2
      DO JN=IM,KSMAX
        ZVOR(JI,JN)=PVOR(JF,KASM0(IM)+(JN-IM)*2+JI-1)
      ENDDO
    ENDDO

!         Linear balance equation

    CALL BALADSM(PNOR,IM,KSMAX,ZVOR(1,IM),ZZ(1,IM))

!         Interface to output array

    IF(LDADD)THEN
      DO JI=1,2
        DO JN=IM,KSMAX
          PZ(JF,KASM0(IM)+(JN-IM)*2+JI-1)=&
           & PZ(JF,KASM0(IM)+(JN-IM)*2+JI-1)+ZZ(JI,JN)  
        ENDDO
      ENDDO
    ELSE
      DO JI=1,2
        DO JN=IM,KSMAX
          PZ(JF,KASM0(IM)+(JN-IM)*2+JI-1)=ZZ(JI,JN)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDDO
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('BALADS',1,ZHOOK_HANDLE)
END SUBROUTINE BALADS
