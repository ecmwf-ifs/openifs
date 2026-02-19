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

MODULE DILATATION_MOD 

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE PARKIND2  ,ONLY : JPRH
USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

PRIVATE
PUBLIC DILAT_CALC, DILAT_DEVIATION, DILAT_CONTROL, DILAT_MAPPING

! Dilatation package package.

! To compute dilatation and contraction matrixes from Legendre polynomials and
! gaussian weights and latitudes on a streched sphere ; and control the
! deviation of the product of these matrixes against the identity matrix.

CONTAINS
!================================================================================
SUBROUTINE DILAT_MAPPING(PSTRET,PMU,PMAPP)

!**** *DILAT_MAPPING* - Compute the map factor for each latitudes

!      Purpose.
!      --------
!           To  compute the map factor, given the stretching
!           factor and the sines of the standard gaussian latitudes

!**    Interface.
!      ----------
!         *CALL *DILAT_MAPPING*

!      Explicit arguments.
!      -------------------
!            PSTRET : stretching factor
!            PMU    : sines of the gaussian latitudes
!            PMAPP  : map factor

!      Implicit arguments.
!      -------------------

!      Method.
!      -------
!         See documentation.

!      Reference.
!      ----------
!         Arpege note 19 (in French).

!      Externals.
!      ---------

!      Author.
!      -------
!         R. El Khatib 19-Jun-2013 from Michel Rochas, DMN (Original 91-01-28).

!      Modifications.
!      --------------
!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB), INTENT(IN)  :: PSTRET
REAL(KIND=JPRD), INTENT(IN)  :: PMU(:)
REAL(KIND=JPRB), INTENT(OUT) :: PMAPP(:)

INTEGER(KIND=JPIM) :: IGLS, IDGNH, JGL

REAL(KIND=JPRH) :: Z_DLSINE,Z_DLTAN

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_MAPPING',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF ( (SIZE(PMU,DIM=1) /= SIZE(PMAPP,DIM=1)) ) THEN
  CALL ABOR1('DILAT_MAPPING : SIZE MISMATCH BETWEEN PMU AND PMAPP')
ENDIF

! Remark : this is an inconsistent mixture of double and quadruple precision. REK
Z_DLTAN=(1.0_JPRB-PSTRET**2)/(1.0_JPRB+PSTRET**2)
Z_DLSINE=2.0_JPRB*PSTRET/(1.0_JPRB+PSTRET**2)
IDGNH=(SIZE(PMU,DIM=1)+1)/2
DO JGL=1,IDGNH
  PMAPP(JGL)=REAL((Z_DLSINE/(1.0_JPRB+Z_DLTAN*PMU(JGL)))**2,JPRB)
ENDDO
DO JGL=1,IDGNH
  IGLS=2*IDGNH-JGL+1
  PMAPP(IGLS)=REAL((Z_DLSINE/(1.0_JPRB-Z_DLTAN*PMU(JGL)))**2,JPRB)
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_MAPPING',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE DILAT_MAPPING
!====================================================================
SUBROUTINE DILAT_CALC(KM,KNSMAX,KNSMIN,KDGNH,PNM,PW,PMAPP,PRPNM,PALFA,PBETA)  

!**** *DILAT_CALC* - Computes dilatation matrix.

!      Purpose.
!      --------
!               Computes dilatation matrix for schmidt transform.

!**    Interface.
!      ----------
!         *CALL *DILAT_CALC(...)

!      Explicit Arguments :
!      --------------------
!           INPUT:
!            KM      - Zonal wave number
!            KNSMAX  - Larger truncation
!            KNSMIN  - smaller truncation
!            KDGNH   - number of points on an hemisphere
!            PNM     - Legendre polynomials small truncation
!            PW      - Gaussian weights
!            PMAPP   - Mapping factor
!            PRPNM   - Legendre polynomials large truncation

!           OUTPUT:
!            PALFA  - Contraction (0) matrix
!            PBETA  - Dilatation  (1) matrix

!      Implicit arguments :
!      --------------------

!      Method.
!      -------
!         See documentation.

!      Reference.
!      ----------
!         Arpege note 19 (in French).

!      Externals.
!      ----------
!         Calls MXMAOP.
!         Is called by SUDIL.

!      Author.
!      -------
!         Michel Rochas, DMN.

!      Modifications.
!      --------------
!        Original : 91-01-28.
!        =07-08-91= Philippe Courtier. Changes Leg. polynomials mapping
!        K. YESSAD: 93-05-11 : cleaning, comments put into English.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad: Aug 2005 : A-level distribution of conf 911. 
!        K. Yessad (June 2009): externalisation.
!        R. El Khatib 20-Jun-2013 Optimization
!        R. El Khatib 07-Mar-2016 Simplification/Optimization
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNSMIN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGNH 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNM(2*KDGNH,KNSMIN-KM+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PW(KDGNH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMAPP(2*KDGNH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPNM(KDGNH,KNSMAX-KM+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALFA(KNSMAX-KM+1,KNSMIN-KM+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBETA(KNSMAX-KM+1,KNSMIN-KM+1)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZP(KNSMAX-KM+1)
REAL(KIND=JPRB) :: ZPNM(KDGNH+1,KNSMIN-KM+1)
REAL(KIND=JPRB) :: ZNOR(2*((KNSMAX-KM+1)/2)+1,KNSMIN-KM+1)
REAL(KIND=JPRB) :: ZSUD(2*((KNSMAX-KM+1)/2)+1,KNSMIN-KM+1)
REAL(KIND=JPRB) :: ZRPNM(2*((KNSMAX-KM+1)/2)+1,KDGNH+1)

INTEGER(KIND=JPIM) :: IGLS, IR, IR2, JGL, JN, JS, JI, II, IOFF
INTEGER(KIND=JPIM) :: ISMAXSUR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_CALC',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*        0. Preparations
!            ------------

ISMAXSUR=2*((KNSMAX-KM+1)/2)+1
IOFF=KNSMAX+2-KM

! Initialize a parity array

DO JS=1,KNSMAX-KM+1
  ZP(JS)=REAL(2*MOD(JS+KNSMAX-KM,2)-1,JPRB)
ENDDO

!   multiplication by gaussian weights
!   of the Legendre polynomials at high truncation

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JGL,JS)
DO JGL=1,KDGNH
  DO JS=1,KNSMAX-KM+1
    ZRPNM(JS,JGL)=PRPNM(JGL,JS)*PW(JGL)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

!*    1. Matrix ALPHA


!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JN,JGL)
DO JN=1,KNSMIN-KM+1
  DO JGL=1,KDGNH
    ZPNM(JGL,JN)=PNM(JGL,JN)*PMAPP(JGL)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

CALL MXMAOP(ZRPNM,1,ISMAXSUR,ZPNM(:,:),1,KDGNH+1,ZNOR(:,:),1,ISMAXSUR,&
 & KNSMAX-KM+1,KDGNH,KNSMIN-KM+1)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JN,JGL)
DO JN=1,KNSMIN-KM+1
  DO JGL=1,KDGNH
    ZPNM(JGL,JN)=PNM(2*KDGNH-JGL+1,JN)*PMAPP(2*KDGNH-JGL+1)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

CALL MXMAOP(ZRPNM,1,ISMAXSUR,ZPNM(:,:),1,KDGNH+1,ZSUD(:,:),1,ISMAXSUR,&
 & KNSMAX-KM+1,KDGNH,KNSMIN-KM+1)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JN,JS)
DO JN=1,KNSMIN-KM+1
  DO JS=1,KNSMAX-KM+1
    PALFA(IOFF-JS,KNSMIN-KM+2-JN)=ZNOR(JS,JN)+ZP(JS)*ZSUD(JS,JN)
  ENDDO
ENDDO
!$OMP END PARALLEL DO


!*    2. Matrix BETA


!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JN,JGL)
DO JN=1,KNSMIN-KM+1
  DO JGL=1,KDGNH
    ZPNM(JGL,JN)=PNM(JGL,JN)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

CALL MXMAOP(ZRPNM,1,ISMAXSUR,ZPNM(:,:),1,KDGNH+1,ZNOR(:,:),1,ISMAXSUR,&
 & KNSMAX-KM+1,KDGNH,KNSMIN-KM+1)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JN,JGL)
DO JN=1,KNSMIN-KM+1
  DO JGL=1,KDGNH
    ZPNM(JGL,JN)=PNM(2*KDGNH-JGL+1,JN)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

CALL MXMAOP(ZRPNM,1,ISMAXSUR,ZPNM(:,:),1,KDGNH+1,ZSUD(:,:),1,ISMAXSUR,&
 & KNSMAX-KM+1,KDGNH,KNSMIN-KM+1)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JN,JS)
DO JN=1,KNSMIN-KM+1
  DO JS=1,KNSMAX-KM+1
    PBETA(IOFF-JS,KNSMIN-KM+2-JN)=ZNOR(JS,JN)+ZP(JS)*ZSUD(JS,JN)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_CALC',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DILAT_CALC
!=============================================================================
SUBROUTINE DILAT_DEVIATION(PALFA,PBETA,PMAX)

!**** *DILAT_DEVIATION* - compute the deviation of dilatation/contraction matrixes.

!      Purpose.
!      --------
!           Compute the deviation from identity of the product contraction o dilatation
!           for a given wave number

!**    Interface.
!      ----------
!         *CALL *DILAT_DEVIATION*

!      Explicit arguments.
!      -------------------
!            PALFA      - Matrix Alfa (Contraction)
!            PBETA      - Matrix Beta (Dilatation)
!            PMAX         - deviation from identity 

!      Implicit arguments.
!      -------------------

!      Method.
!      -------
!         See documentation.

!      Reference.
!      ----------
!         Arpege note 19 (in French).

!      Externals.
!      ---------

!      Author.
!      -------
!         R. El Khatib 19-Jun-2013 from Michel Rochas, DMN (Original 91-01-28).

!      Modifications.
!      --------------
!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB), INTENT(IN)  :: PALFA(:,:)
REAL(KIND=JPRB), INTENT(IN)  :: PBETA(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PMAX

INTEGER(KIND=JPIM) :: ID1, ID2, JN1, JN2

REAL(KIND=JPRB), ALLOCATABLE :: ZRESUL(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_DEVIATION',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

ID1=SIZE(PALFA,DIM=1)
ID2=SIZE(PALFA,DIM=2)

IF ( (SIZE(PBETA,DIM=1) /= ID1) .OR. (SIZE(PBETA,DIM=2) /= ID2) ) THEN
  CALL ABOR1('DILAT_DEVIATION : SIZES MISMATCH BETWEEN PALFA AND PBETA')
ENDIF

PMAX=-HUGE(1._JPRB)
ALLOCATE(ZRESUL(ID2,ID2))
CALL MXMAOP(PBETA(:,:),ID1,1,PALFA(:,:),1,ID1,ZRESUL,1,ID2,ID2,ID1,ID2)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN1)
DO JN1=1,ID2
  ZRESUL(JN1,JN1)=ZRESUL(JN1,JN1)-1.0_JPRB
ENDDO
!$OMP END PARALLEL DO
DO JN1=1,ID2
  DO JN2=1,ID2
    PMAX=MAX(PMAX,ABS(ZRESUL(JN1,JN2)))
  ENDDO
ENDDO
DEALLOCATE(ZRESUL)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_DEVIATION',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE DILAT_DEVIATION
!====================================================================
SUBROUTINE DILAT_CONTROL(KMYPROC,KOUTPROC,KWSET,KULOUT,PMAXALL)

!**** *DILAT_CONTROL* - Control the dilatation/contraction matrixes.

!      Purpose.
!      --------
!           Print out the deviation from identity of the product contraction o dilatation
!           for a contiguous set of distributed  wave numbers

!**    Interface.
!      ----------
!         *CALL *DILAT_CONTROL*

!      Explicit arguments.
!      -------------------
!            KMYPROC      - Current mpi task
!            KOUTPROC     - task in charge of writing out
!            KWSET        - wave set for all wave numbers
!            KULOUT       - Output logical unit number
!            PMAXALL      - deviation from identity of all wave numbers

!      Implicit arguments.
!      -------------------

!      Method.
!      -------
!         See documentation.

!      Reference.
!      ----------
!         Arpege note 19 (in French).

!      Externals.
!      ---------

!      Author.
!      -------
!         R. El Khatib 19-Jun-2013 from Michel Rochas, DMN (Original 91-01-28).

!      Modifications.
!      --------------
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KWSET(0:)
INTEGER(KIND=JPIM),INTENT(IN) :: KMYPROC
INTEGER(KIND=JPIM),INTENT(IN) :: KOUTPROC
INTEGER(KIND=JPIM),INTENT(IN) :: KULOUT
REAL(KIND=JPRB),   INTENT(IN) :: PMAXALL(0:)

INTEGER(KIND=JPIM) :: JM
REAL(KIND=JPRB) :: ZMAXRECV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_CONTROL',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF ( (UBOUND(PMAXALL,DIM=1) /= UBOUND(KWSET,DIM=1)) ) THEN
  CALL ABOR1('DILAT_CONTROL : UBOUNDS MISMATCH BETWEEN PMAXALL AND KWSET')
ENDIF

DO JM=0,UBOUND(PMAXALL,DIM=1)
  IF (KWSET(JM) == KMYPROC) THEN
    IF (KMYPROC == KOUTPROC) THEN
      WRITE(KULOUT,'('' ZONAL WAVE NUMBER '',I4, &
       & '' DEVIATION FROM IDENTITY MATRIX '',E10.3)') JM,PMAXALL(JM)
    ELSE
      CALL MPL_SEND(PMAXALL(JM),KDEST=KOUTPROC,KTAG=JM,CDSTRING='DILAT_CONTROL:')
    ENDIF
  ELSEIF(KMYPROC == KOUTPROC) THEN 
    CALL MPL_RECV(ZMAXRECV,KSOURCE=KWSET(JM),KTAG=JM,CDSTRING='DILAT_CONTROL:')
    WRITE(KULOUT,'('' ZONAL WAVE NUMBER '',I4, &
     & '' DEVIATION FROM IDENTITY MATRIX '',E10.3)') JM,ZMAXRECV
  ENDIF
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DILAT_CONTROL',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DILAT_CONTROL
!=============================================================================
END MODULE DILATATION_MOD
