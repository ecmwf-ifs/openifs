! (C) Copyright 1989- Meteo-France.

SUBROUTINE DECVLANISO(KST,KEND,KPROMA,LDGAMA,LDALFA,PGAMA,PALFA,PU,PV,PSDO)  

!**** *DECVLANISO*  - DECONVOLUTION of the vector 
!                 (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))
!                 where :
!                 sigma = standard deviation of orography    
!                 gamma = Anisotropy coefficient of topography    
!                 alpha = Direction of the principal axis of the topography

!     PURPOSE.
!     --------
!        Deconvolution of the vector 
!        (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))
!        which components can be interpolated/extrapolated in place of 
!        Anisotropy coefficient or Direction of the principal axis 
!        of the topography.

!**   INTERFACE.
!     ----------
!       *CALL* *DECVLANISO*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KST    : first point in row.
!         KEND   : last point in row.
!         KPROMA : length of the row.
!         LDGAMA : .TRUE. <=> anisotropy
!         LDALFA : .TRUE. <=> direction
!         PU,PV  : components of vector 
!                 (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))
!         PSDO   : standart deviation of orography times g

!        OUTPUT:
!         PGAMA  : Anisotropy coefficient of topography
!         PALFA  : Direction of the principal axis of the topography

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

USE YOMCST, ONLY : RPI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDGAMA 
LOGICAL           ,INTENT(IN)    :: LDALFA 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PGAMA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PALFA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PU(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PV(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDO(KPROMA) 
!     ZMOD : (sigma**2)*(1-gamma)/(1+gamma)
!     ZSDO2: PSDO**2

REAL(KIND=JPRB) :: ZMOD(KPROMA), ZSDO2(KPROMA)

INTEGER(KIND=JPIM) :: JI

REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DECVLANISO',0,ZHOOK_HANDLE)


!*       2. DECONVOLUTION
!           -------------

  ZEPS=EPSILON(1.0_JPRB)*10000._JPRB
  DO JI=KST,KEND
    ZMOD(JI)=SQRT( PU(JI)**2 + PV(JI)**2 )
  ENDDO

!       Anisotropy of orography :
  IF (LDGAMA) THEN
!         When the tensor is isotropic, the anisotropy is 1.
    PGAMA(KST:KEND)=1.0_JPRB
    ZSDO2(KST:KEND)=PSDO(KST:KEND)*PSDO(KST:KEND)
    DO JI=KST,KEND
      IF (ZMOD(JI) >= ZEPS) THEN
        PGAMA(JI)=MAX(0.0_JPRB, (ZSDO2(JI)-ZMOD(JI))&
         & /MAX(ZEPS,ZSDO2(JI)+ZMOD(JI)))  
      ENDIF
    ENDDO
  ENDIF

!       Direction of orography :
  IF (LDALFA) THEN
!         When the tensor is isotropic, the direction is 0.
    PALFA(KST:KEND)=0.0_JPRB
    DO JI=KST,KEND
      IF (ZMOD(JI) >= ZEPS) THEN
        IF (ABS(PU(JI)) <= ZEPS) THEN
          IF (PV(JI) > 0) THEN
            PALFA(JI)=RPI*0.25_JPRB
          ELSE
            PALFA(JI)=-RPI*0.25_JPRB
          ENDIF
        ELSE
          IF (PU(JI) > 0) THEN
            PALFA(JI)=0.5_JPRB*ATAN(PV(JI)/PU(JI))
          ELSEIF (PV(JI) > 0) THEN
            PALFA(JI)=0.5_JPRB*(RPI+ATAN(PV(JI)/PU(JI)))
          ELSE
            PALFA(JI)=0.5_JPRB*(ATAN(PV(JI)/PU(JI))-RPI)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
!         Set final direction between -pi/2 and pi/2 :
    DO JI=KST,KEND
      PALFA(JI)=ATAN(TAN(PALFA(JI)))
    ENDDO
  ENDIF


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DECVLANISO',1,ZHOOK_HANDLE)
END SUBROUTINE DECVLANISO
