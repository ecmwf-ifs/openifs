! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GRADI (KIJS, KIJL, NINF, NSUP, IREFRA,   &
 &                  BLK2GLO,                        &
 &                  DPTHEXT, UEXT, VEXT,            &
 &                  DDPHI, DDLAM, DUPHI,            &
 &                  DULAM, DVPHI, DVLAM)

! ----------------------------------------------------------------------

!**** *GRADI* - CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS.

!     K.P. HUBBERT              AUGUST   1988
!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS OF A BLOCK.

!**   INTERFACE.
!     ----------

!       *CALL* *GRADI (KIJS, KIJL, NINF, NSUP, IREFRA,
!                      BLK2GLO,
!                      DPTHEXT, UEXT, VEXT,
!                      DDPHI, DDLAM, DUPHI,
!                      DULAM, DVPHI, DVLAM)*
!          *KIJS*   - STARTING INDEX
!          *KIJL*   - ENDING INDEX
!          *NINF:NSUP+1* : DIMENSION OF DPTHEXT, UEXT, VEXT 
!          *IREFRA* - REFRACTION OPTION.
!          *BLK2GLO*- BLOCK TO GRID TRANSFORMATION
!          *DDPHI*  - LATITUDE DEPTH GRADIENT.
!          *DDLAM*  - LONGITUDE DEPTH GRADIENT.
!          *DUPHI*  - LATITUDE  U-COMPONENT GRADIENT.
!          *DULAM*  - LONGITUDE U-COMPONENT GRADIENT.
!          *DVPHI*  - LATITUDE  V-COMPONENT GRADIENT.
!          *DVLAM*  - LONGITUDE V-COMPONENT GRADIENT.

!     METHOD.
!     ------

!       CENTRAL DIFFERENCING FOR DEPTH AND CURRENT VELOCITY GRADIENTS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWCURR  , ONLY : CURRENT_GRADIENT_MAX
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,COSPH
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,WLAT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: NINF 
      INTEGER(KIND=JWIM), INTENT(IN) :: NSUP 
      INTEGER(KIND=JWIM), INTENT(IN) :: IREFRA
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: DPTHEXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: UEXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: VEXT 

      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: DDPHI
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: DDLAM
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: DUPHI
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: DULAM
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: DVPHI
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: DVLAM


      INTEGER(KIND=JWIM) :: NLAND, IJ, IPP, IPM, IPP2, IPM2, ILP, ILM, KX
      REAL(KIND=JWRB) :: DPTP, DPTM, UP, UM, VP, VM 
      REAL(KIND=JWRB) :: CGMAX
      REAL(KIND=JWRB) :: ONEO2DELPHI
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRADI',0,ZHOOK_HANDLE)


!*    1. INITIALISE.
!        -----------

      NLAND=NSUP+1
      ONEO2DELPHI = 0.5_JWRB/DELPHI

! ----------------------------------------------------------------------

!*    2. CALCULATE DEPTH GRADIENTS.
!        --------------------------

      IF (IREFRA == 1 .OR. IREFRA == 3) THEN
        DO IJ=KIJS,KIJL
          IPP = KLAT(IJ,2,1)
          IPM = KLAT(IJ,1,1)
          IPP2 = KLAT(IJ,2,2)
          IPM2 = KLAT(IJ,1,2)
          IF (IPP /= NLAND  .AND. IPM /= NLAND .AND.                    &
     &        IPP2 /= NLAND .AND. IPM2 /= NLAND ) THEN
            DPTP=WLAT(IJ,2)*DPTHEXT(IPP)+(1.0_JWRB-WLAT(IJ,2))*DPTHEXT(IPP2)
            DPTM=WLAT(IJ,1)*DPTHEXT(IPM)+(1.0_JWRB-WLAT(IJ,1))*DPTHEXT(IPM2)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSEIF (IPP /= NLAND .AND. IPM /= NLAND) THEN
            DPTP=DPTHEXT(IPP)
            DPTM=DPTHEXT(IPM)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSEIF (IPP2 /= NLAND .AND. IPM2 /= NLAND) THEN
            DPTP=DPTHEXT(IPP2)
            DPTM=DPTHEXT(IPM2)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSE
            DDPHI(IJ) = 0.0_JWRB
          ENDIF
          ILP = KLON(IJ,2)
          ILM = KLON(IJ,1)
          KX  = BLK2GLO%KXLT(IJ)
          IF (ILP /= NLAND .AND. ILM /= NLAND) THEN
            DDLAM(IJ)=(DPTHEXT(ILP)-DPTHEXT(ILM))/(2._JWRB*DELLAM(KX))
          ELSE
            DDLAM(IJ) = 0.0_JWRB 
          ENDIF
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
          DDPHI(IJ) = 0.0_JWRB
          DDLAM(IJ) = 0.0_JWRB
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    3. CALCULATE CURRENT VELOCITY GRADIENTS.
!        -------------------------------------

      IF (IREFRA == 2 .OR. IREFRA == 3) THEN
        DO IJ=KIJS,KIJL
          IPP = KLAT(IJ,2,1)
!         exact 0 means that the current field was not defined, hence
!         no gradient should be extrapolated
          IF (UEXT(IPP) == 0.0_JWRB .AND. VEXT(IPP) == 0.0_JWRB) IPP = NLAND
          IPM = KLAT(IJ,1,1)
          IF (UEXT(IPM) == 0.0_JWRB .AND. VEXT(IPM) == 0.0_JWRB) IPM = NLAND
          IPP2 = KLAT(IJ,2,2)
          IF (UEXT(IPP2) == 0.0_JWRB .AND. VEXT(IPP2) == 0.0_JWRB) IPP2 = NLAND
          IPM2 = KLAT(IJ,1,2)
          IF (UEXT(IPM2) == 0.0_JWRB .AND. VEXT(IPM2) == 0.0_JWRB) IPM2 = NLAND

          IF (IPP /= NLAND .AND. IPM /= NLAND .AND.                     &
     &        IPP2 /= NLAND .AND. IPM2 /= NLAND) THEN
            UP = WLAT(IJ,2)*UEXT(IPP)+(1.0_JWRB-WLAT(IJ,2))*UEXT(IPP2)
            VP = WLAT(IJ,2)*VEXT(IPP)+(1.0_JWRB-WLAT(IJ,2))*VEXT(IPP2)
            UM = WLAT(IJ,1)*UEXT(IPM)+(1.0_JWRB-WLAT(IJ,1))*UEXT(IPM2)
            VM = WLAT(IJ,1)*VEXT(IPM)+(1.0_JWRB-WLAT(IJ,1))*VEXT(IPM2)
            DUPHI(IJ) = (UP-UM)*ONEO2DELPHI
            DVPHI(IJ) = (VP-VM)*ONEO2DELPHI
          ELSEIF (IPP /= NLAND .AND. IPM /= NLAND) THEN
            UP = UEXT(IPP)
            VP = VEXT(IPP)
            UM = UEXT(IPM)
            VM = VEXT(IPM)
            DUPHI(IJ) = (UP-UM)*ONEO2DELPHI
            DVPHI(IJ) = (VP-VM)*ONEO2DELPHI
          ELSE
            DUPHI(IJ) = 0.0_JWRB
            DVPHI(IJ) = 0.0_JWRB
          ENDIF

          ILP = KLON(IJ,2)
!         exact 0 means that the current field was not defined, hence
!         no gradient should be extrapolated
          IF (UEXT(ILP) == 0.0_JWRB .AND. VEXT(ILP) == 0.0_JWRB) ILP=NLAND
          ILM = KLON(IJ,1)
          IF (UEXT(ILM) == 0.0_JWRB .AND. VEXT(ILM) == 0.0_JWRB) ILM=NLAND
          KX  = BLK2GLO%KXLT(IJ)
          IF (ILP /= NLAND .AND. ILM /= NLAND) THEN
            DULAM(IJ) = (UEXT(ILP)-UEXT(ILM))/(2.0_JWRB*DELLAM(KX))
            DVLAM(IJ) = (VEXT(ILP)-VEXT(ILM))/(2.0_JWRB*DELLAM(KX))
          ELSE
            DULAM(IJ) = 0.0_JWRB
            DVLAM(IJ) = 0.0_JWRB
          ENDIF
        ENDDO

        DO IJ=KIJS,KIJL
          KX  = BLK2GLO%KXLT(IJ)
          CGMAX = CURRENT_GRADIENT_MAX*COSPH(KX)
          DUPHI(IJ) = SIGN(MIN(ABS(DUPHI(IJ)),CGMAX),DUPHI(IJ))
          DVPHI(IJ) = SIGN(MIN(ABS(DVPHI(IJ)),CGMAX),DVPHI(IJ))
          DULAM(IJ) = SIGN(MIN(ABS(DULAM(IJ)),CGMAX),DULAM(IJ))
          DVLAM(IJ) = SIGN(MIN(ABS(DVLAM(IJ)),CGMAX),DVLAM(IJ))
        ENDDO

      ELSE
        DO IJ=KIJS,KIJL
          DUPHI(IJ) = 0.0_JWRB
          DVPHI(IJ) = 0.0_JWRB
          DULAM(IJ) = 0.0_JWRB
          DVLAM(IJ) = 0.0_JWRB
        ENDDO
      ENDIF

IF (LHOOK) CALL DR_HOOK('GRADI',1,ZHOOK_HANDLE)

END SUBROUTINE GRADI
