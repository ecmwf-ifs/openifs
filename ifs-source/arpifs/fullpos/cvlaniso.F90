! (C) Copyright 1989- Meteo-France.

SUBROUTINE CVLANISO(KST,KEND,KPROMA,LDGAMA,LDALFA,PGAMA,PALFA,PU,PV,PSDO)  

!**** *CVLANISO*  - CONVOLUTION of the vector 
!                 (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))
!                 where :
!                 sigma = standard deviation of orography    
!                 gamma = Anisotropy coefficient of topography    
!                 alpha = Direction of the principal axis of the topography

!     PURPOSE.
!     --------
!        Convolution of the vector 
!        (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))
!        which components can be interpolated/extrapolated in place of 
!        Anisotropy coefficient or Direction of the principal axis 
!        of the topography.

!**   INTERFACE.
!     ----------
!       *CALL* *CVLANISO*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KST    : first point in row.
!         KEND   : last point in row.
!         KPROMA : length of the row.
!         LDGAMA : .TRUE. <=> U component
!         LDALFA : .TRUE. <=> V component
!         PGAMA  : Anisotropy coefficient of topography
!         PALFA  : Direction of the principal axis of the topography
!         PSDO   : standart deviation of orography times g

!        INPUT/OUTPUT:
!         PU,PV  : components of vector 
!                 (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))

!        IMPLICIT ARGUMENTS
!        --------------------
!          See module above

!     METHOD.
!     -------
!          See documentation.

!     EXTERNALS.
!     ----------
!      None

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-09-10 (from HPOS and FPINTPHY)
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDGAMA 
LOGICAL           ,INTENT(IN)    :: LDALFA 
REAL(KIND=JPRB)   ,INTENT(IN) :: PGAMA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PALFA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PU(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PV(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDO(KPROMA) 
!     ZSDO2: PSDO**2

REAL(KIND=JPRB) :: ZSDO2(KPROMA)

INTEGER(KIND=JPIM) :: JI

LOGICAL :: LLU, LLV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CVLANISO',0,ZHOOK_HANDLE)


!*       1. CONVOLUTION
!           -----------

  DO JI=KST,KEND
    ZSDO2(JI)=PSDO(JI)**2
  ENDDO
  LLU=LDGAMA
  LLV=LDALFA

!       first component of vector 
  IF (LLU) THEN
    DO JI=KST,KEND
      PU(JI)=ZSDO2(JI)*((1.0_JPRB-PGAMA(JI))/(1.0_JPRB+PGAMA(JI)))*&
       & ( COS(PALFA(JI))**2 - SIN(PALFA(JI))**2 )  
    ENDDO
  ENDIF

!       Second component of vector 
  IF (LLV) THEN
    DO JI=KST,KEND
      PV(JI)=ZSDO2(JI)*((1.0_JPRB-PGAMA(JI))/(1.0_JPRB+PGAMA(JI)))*&
       & ( 2.0_JPRB*SIN(PALFA(JI))*COS(PALFA(JI)) )  
    ENDDO
  ENDIF


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CVLANISO',1,ZHOOK_HANDLE)
END SUBROUTINE CVLANISO
