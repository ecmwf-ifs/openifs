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

MODULE SUVFE_HLP

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

REAL(KIND=JPRB) :: RTMIN, RTMAX

CONTAINS

! ======================================================================================================
! ======================================================================================================

!     ------------------------------------------------------------------
!     Additional subroutines and functions: MULPOL, IPOL, EVPOL, DPOL,
!      GLOBAL2LOCAL, SETTM, SAMPLE_VALUE, SETGP2VFE
!     ------------------------------------------------------------------

!----------------------------------
! multiplication of two polynomials
!----------------------------------
! c(t) = a(t)*b(t)
! input:
! a(t) = sum_(i,1,m) a_i t^(i-1)
! b(t) = sum_(i,1,n) b_i t^(i-1)
! output:
! c(t) = sum_(i,1,n+m-1) c_i t^(i-1)
!-----------------------------------
SUBROUTINE MULPOL(KM,PA,KN,PB,KNP,PC)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

! input
INTEGER(KIND=JPIM),INTENT(IN)           :: KM
REAL   (KIND=JPRB),INTENT(IN)           :: PA(0:KM-1)
INTEGER(KIND=JPIM),INTENT(IN)           :: KN
REAL   (KIND=JPRB),INTENT(IN)           :: PB(0:KN-1)
INTEGER(KIND=JPIM),INTENT(IN)           :: KNP
! output
REAL   (KIND=JPRB),INTENT(OUT)          :: PC(0:KNP-1)

! local vars
INTEGER(KIND=JPIM)           :: JI, JJ
INTEGER(KIND=JPIM)           :: IN, IM, INP, IK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MULPOL',0,ZHOOK_HANDLE)

! having indexes from 0 is easier to code
IM  = KM-1
IN  = KN-1
INP = KNP-1

! general convolution
DO JI=0,INP
  PC(JI) = 0.0_JPRB
  DO JJ=0,JI
    IK = JI-JJ
    IF( JJ <= IM .AND. IK <= IN )THEN
      PC(JI) = PC(JI) + PA(JJ)*PB(IK)
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('MULPOL',1,ZHOOK_HANDLE)
END SUBROUTINE MULPOL


!------------------------
! integrate polynomial
! PIP(t) = PC +  int_(0,t) P(t') dt'
! PC - additive constant
!------------------------
SUBROUTINE IPOL(KN,P,PC,KNI,PIP)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMLUN    ,ONLY : NULOUT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

! input:
INTEGER(KIND=JPIM),INTENT(IN)  :: KN
REAL   (KIND=JPRB),INTENT(IN)  :: P(KN)
REAL   (KIND=JPRB),INTENT(IN)  :: PC
INTEGER(KIND=JPIM),INTENT(IN)  :: KNI

! output:
REAL   (KIND=JPRB),INTENT(OUT) :: PIP(KNI)

! aux
INTEGER(KIND=JPIM) :: JI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('IPOL',0,ZHOOK_HANDLE)

IF( KNI /= (KN+1) )THEN
  WRITE(NULOUT,*) "(ERROR) IPOL DIM OF PIP IS:",KNI," BUT SHALL BE:",KN+1
ENDIF

! additive constant 
PIP(1) = PC
DO JI=2,KN+1
  PIP(JI) = P(JI-1)/REAL(JI-1,JPRB)
ENDDO

IF (LHOOK) CALL DR_HOOK('IPOL',1,ZHOOK_HANDLE)
END SUBROUTINE IPOL


! ----------------------------------
! evaluate polynomial at given value
! EVPOL = sum_{i,1,N} POL(i) PX^{i-1}
! ----------------------------------
FUNCTION EVPOL(KORDER,P,PX)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL   (KIND=JPRB) :: EVPOL

! input:
INTEGER(KIND=JPIM),INTENT(IN) :: KORDER
REAL   (KIND=JPRB),INTENT(IN) :: P(KORDER)
REAL   (KIND=JPRB),INTENT(IN) :: PX

! aux
REAL   (KIND=JPRB) :: ZP(KORDER), ZQ(KORDER - 1)
INTEGER(KIND=JPIM) :: I, N
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EVPOL',0,ZHOOK_HANDLE)


! evaluation vis synthetic division
! synthetic division says:
! if you divide polynomial p(eta) by (eta - e0)
! then remainder after division is polynomial
! evaluated at point e0 -> p(e0)

! leading coef
DO I = 1, KORDER
  ZP(I) = P(KORDER - I + 1)
ENDDO

N = KORDER - 1

ZQ(1) = ZP(1)
DO I = 2, N
  ZQ(I) = (PX*ZQ(I-1)) + ZP(I)
ENDDO

! EVPOL = remainder
EVPOL = PX*ZQ(N) + ZP(N+1)

IF (LHOOK) CALL DR_HOOK('EVPOL',1,ZHOOK_HANDLE)
END FUNCTION EVPOL


!------------------------
! derivate polynomial
! PDP(t) = dP/dt
!------------------------
SUBROUTINE DPOL(KN,P,KND,PDP)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMLUN    ,ONLY : NULOUT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

! input:
INTEGER(KIND=JPIM),INTENT(IN)  :: KN
REAL   (KIND=JPRB),INTENT(IN)  :: P(KN)
INTEGER(KIND=JPIM),INTENT(IN)  :: KND

! output:
REAL   (KIND=JPRB),INTENT(OUT) :: PDP(KN-1)

! aux
INTEGER(KIND=JPIM) :: JI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DPOL',0,ZHOOK_HANDLE)

IF( KND /= (KN-1) )THEN
  WRITE(NULOUT,*) "(ERROR) DPOL DIM OF PDP IS:",KND," BUT SHALL BE:",KN-1 
ENDIF

DO JI=1,KN-1
  PDP(JI) = P(JI+1)*REAL(JI,JPRB)
ENDDO

IF (LHOOK) CALL DR_HOOK('DPOL',1,ZHOOK_HANDLE)
END SUBROUTINE DPOL


!-------------------------------------------
! global to local coordinate transformation:
! p(x) in <xmin;xmax) -> p(t) in <0,1)
! 
! t = (x - xmin)/(xmax - xmin)
!
! IN:  pxmin, pxmax, korder, pol_x
! OUT: pol_t
!-------------------------------------------
SUBROUTINE GLOBAL2LOCAL(KORDER,PTM,PXMIN,PXMAX,POLX,POLT)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

! input:
INTEGER(KIND=JPIM),INTENT(IN)  :: KORDER
REAL   (KIND=JPRB),INTENT(IN)  :: PTM(KORDER,KORDER)
REAL   (KIND=JPRB),INTENT(IN)  :: PXMIN
REAL   (KIND=JPRB),INTENT(IN)  :: PXMAX
REAL   (KIND=JPRB),INTENT(IN)  :: POLX(KORDER)
! output:
REAL   (KIND=JPRB),INTENT(OUT) :: POLT(KORDER)

REAL   (KIND=JPRB) :: ZT, ZX
REAL   (KIND=JPRB) :: ZSAMPLE(KORDER)
INTEGER(KIND=JPIM) :: JI

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GLOBAL2LOCAL',0,ZHOOK_HANDLE)

!-------------------------
DO JI = 1, KORDER
  ZT         = SAMPLE_VALUE(KORDER,JI)
  ZX         = FT2X(PXMIN, PXMAX, ZT)
  ZSAMPLE(JI)= EVPOL(KORDER,POLX,ZX)
ENDDO
POLT = MATMUL(PTM,ZSAMPLE)

IF (LHOOK) CALL DR_HOOK('GLOBAL2LOCAL',1,ZHOOK_HANDLE)
END SUBROUTINE GLOBAL2LOCAL


!-------------------------
! setup transform matrix
! (pxmin,pxmax) interval of input values for which polynomial is computed
! in : korder  - dimension of PTM
! out: ptm     - transformation matrix
!-------------------------
SUBROUTINE SETTM(PXMIN, PXMAX, KORDER,PTM)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL(KIND=JPRB),INTENT(IN)     :: PXMIN, PXMAX
INTEGER(KIND=JPIM),INTENT(IN)  :: KORDER
REAL(KIND=JPRB),INTENT(OUT)    :: PTM(KORDER,KORDER)

INTEGER(KIND=JPIM) :: JI,JJ
REAL   (KIND=JPRB) :: ZTM(KORDER,KORDER)
REAL   (KIND=JPRB) :: ZX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "minv.h"

IF (LHOOK) CALL DR_HOOK('SETTM',0,ZHOOK_HANDLE)

DO JI = 1, KORDER
  ZX = FT2X(PXMIN, PXMAX, SAMPLE_VALUE(KORDER,JI))
  PTM(JI,1) = 1.0_JPRB
  DO JJ=2,KORDER
    PTM(JI,JJ) = ZX*PTM(JI,JJ-1)
  ENDDO
ENDDO

CALL MINV_CALLER(.TRUE.,KORDER,PTM,ZTM)
PTM(:,:)=ZTM(:,:)

! CALL PRINT_MATRIX("SETTM : ",KORDER,KORDER,PTM)

IF (LHOOK) CALL DR_HOOK('SETTM',1,ZHOOK_HANDLE)
END SUBROUTINE SETTM


!--------------------------------
! sample valued in interval <0,1)
!--------------------------------
FUNCTION SAMPLE_VALUE(KN,KI) 

USE PARKIND1  ,ONLY : JPIM , JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

REAL    (KIND=JPRB) :: SAMPLE_VALUE
INTEGER (KIND=JPIM),INTENT(IN)  :: KN,KI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SAMPLE_VALUE',0,ZHOOK_HANDLE)

SAMPLE_VALUE = REAL(KI-1,JPRB)/REAL(KN-1,JPRB)
SAMPLE_VALUE = RTMIN + SAMPLE_VALUE * (RTMAX - RTMIN)

IF (LHOOK) CALL DR_HOOK('SAMPLE_VALUE',1,ZHOOK_HANDLE)
END FUNCTION SAMPLE_VALUE

!--------------------------------
! sample valued in interval <-1/2, +1/2)
!--------------------------------
FUNCTION FX2T(PXMIN, PXMAX, PX) 

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL    (KIND=JPRB) :: FX2T
REAL    (KIND=JPRB),INTENT(IN) :: PXMIN, PXMAX, PX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FX2T',0,ZHOOK_HANDLE)

FX2T = (PX - PXMIN) / (PXMAX - PXMIN)
! FX2T = (PX - 0.5_JPRB * (PXMAX + PXMIN)) / (PXMAX - PXMIN)

IF (LHOOK) CALL DR_HOOK('FX2T',1,ZHOOK_HANDLE)
END FUNCTION FX2T

FUNCTION FT2X(PXMIN, PXMAX, PT) 

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL    (KIND=JPRB) :: FT2X
REAL    (KIND=JPRB),INTENT(IN) :: PXMIN, PXMAX, PT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FT2X',0,ZHOOK_HANDLE)

FT2X = PXMIN + PT * (PXMAX - PXMIN)
! FT2X = 0.5_JPRB * (PXMAX + PXMIN) + PT * (PXMAX - PXMIN)

IF (LHOOK) CALL DR_HOOK('FT2X',1,ZHOOK_HANDLE)
END FUNCTION FT2X


!----------------------------------
! set GP to VFE (or back)
!----------------------------------
SUBROUTINE SETGP2VFE(LDAPPROX,PRMINDETA,KIN,KOUT,PETA,KTYPE,KORDER,KBASIS,PSPLINE,KOFF,KKNOT,PKNOT,PMTF)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULOUT

IMPLICIT NONE

! PMTF(M,N) - transformation matrix 
!     - transforms from KIN dimensional space into KOUT dimensional space 

LOGICAL,INTENT(IN)              :: LDAPPROX         ! if function is interpolated or approximated
REAL(KIND=JPRB),INTENT(IN)      :: PRMINDETA
INTEGER(KIND=JPIM),INTENT(IN)   :: KIN              ! dimension of input
INTEGER(KIND=JPIM),INTENT(IN)   :: KOUT             ! dimension of output
REAL(KIND=JPRB),INTENT(IN)      :: PETA (KOUT)      ! eta locations of output points
INTEGER(KIND=JPIM),INTENT(IN)   :: KTYPE(KOUT)      ! KTYPE=0 value, KTYPE=1 derivative 
INTEGER(KIND=JPIM),INTENT(IN)   :: KORDER           ! order of basis being evaluated in order
                                                    ! to define PMTF
INTEGER(KIND=JPIM),INTENT(IN)   :: KBASIS           ! number of basis functions of PSPLINE array
REAL(KIND=JPRB),INTENT(IN)      :: PSPLINE(KBASIS,KORDER,KORDER)    ! spline coefficients on local coordinates
INTEGER(KIND=JPIM),INTENT(IN)   :: KOFF             ! index of first non-zero function in PSPLINE
INTEGER(KIND=JPIM),INTENT(IN)   :: KKNOT            ! dimension of knot vector
REAL(KIND=JPRB),INTENT(IN)      :: PKNOT(KKNOT)     ! knot vector
REAL(KIND=JPRB),INTENT(OUT)     :: PMTF(KOUT,KIN)   ! transformation matrix

! local variables
INTEGER(KIND=JPIM)  :: JR,JC,JI
INTEGER(KIND=JPIM)  :: IK,IMIN,IMAX,ISEG,IFRST,ILAST
REAL(KIND=JPRB)     :: ZTT, ZDTDETA
REAL(KIND=JPRB)     :: ZPOL_DER(KORDER-1)
REAL(KIND=JPRB)     :: ZPOL_DDER(KORDER-2)
REAL(KIND=JPHOOK)     :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('SETGP2VFE',0,ZHOOK_HANDLE)

! index offset
IFRST =   1+KOFF
ILAST = KIN+KOFF

! initialize matrix
PMTF   = 0.0_JPRB

DO JR=1,KOUT

  ! find in which interval of knots is eta_full(i)
  ! interval is: t(k) < eta_full(ii) < t(k+1)
  IK = -1
  DO JI=1,KKNOT-1
    IF( ABS(PKNOT(JI+1)-PKNOT(JI)) > PRMINDETA )THEN
      IF( (PETA(JR) >= PKNOT(JI)) .AND. (PETA(JR) <= PKNOT(JI+1)) )THEN
        IK = JI
      ENDIF
    ENDIF
  ENDDO
  IF( IK == -1 )THEN
    WRITE(NULOUT,*) "(E) ERROR IN SETGP2VFE. ETA NOT FOUND WITHIN KNOTS"
    CALL  EXIT(-1)
  ENDIF

  ! minimum and maximum basis functions above interval
  ! take care to index only basis functions i=IFRST,ILAST
  IMIN = MAX(IFRST     ,IK + 1 - KORDER)
  IMAX = MIN(ILAST     ,IK + 1 - 1     )

  ! d p(eta) / deta = d p(deta)/dt * dt/deta
  ! this is factor dtdeta
  ZDTDETA = 1.0_JPRB / (PKNOT(IK+1) - PKNOT(IK))
  ZTT     = FX2T(PKNOT(IK), PKNOT(IK+1), PETA(JR))

  IF( KTYPE(JR) == 0 )THEN

    ! if VALUE is prescibed
    IF( LDAPPROX )THEN
      DO JC=IMIN,IMAX
        IF( JR == (JC-KOFF) )THEN
          PMTF(JR,JC-KOFF) = 1.0_JPRB
        ENDIF
      ENDDO
    ELSE
      ! evaluate all functions above this interval 
      ! and put them at the same row JR
      DO JC=IMIN,IMAX
        ISEG = IK + 1 - JC
        PMTF(JR,JC-KOFF) = EVPOL(KORDER,PSPLINE(JC,ISEG,:),ZTT)
      ENDDO
    ENDIF

  ELSEIF( KTYPE(JR) == 1 )THEN

    ! if DERIVATIVE is prescibed
    DO JC=IMIN,IMAX
      ISEG = IK + 1 - JC 
      ZPOL_DER = 0.0_JPRB
      CALL DPOL(KORDER,PSPLINE(JC,ISEG,:),KORDER-1,ZPOL_DER)      
      PMTF(JR,JC-KOFF) = EVPOL(KORDER-1,ZPOL_DER,ZTT) * ZDTDETA
    ENDDO

  ELSEIF( KTYPE(JR) == 2 )THEN

    ! if SECOND DERIVATIVE is prescibed
    DO JC=IMIN,IMAX
      ISEG = IK + 1 - JC
      ZPOL_DER = 0.0_JPRB
      CALL DPOL(KORDER,PSPLINE(JC,ISEG,:),KORDER-1,ZPOL_DER)
      CALL DPOL(KORDER-1,ZPOL_DER,KORDER-2,ZPOL_DDER)
      PMTF(JR,JC-KOFF) = EVPOL(KORDER-2,ZPOL_DDER,ZTT) * (ZDTDETA**2)
    ENDDO

  ELSE

    CALL ABOR1("ERROR: UNDEFINED KTYPE IN SETGP2VFE")

  ENDIF

ENDDO ! JR

IF (LHOOK) CALL DR_HOOK('SETGP2VFE',1,ZHOOK_HANDLE)
END SUBROUTINE SETGP2VFE

END MODULE SUVFE_HLP

