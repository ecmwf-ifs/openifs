!     path:      $Source:rrtmg_sw_vrtqdr /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_vrtqdr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 22:22:22 $
!
      MODULE RRTMG_SW_VRTQDR

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!      use parrrsw, only: ngptsw

      IMPLICIT NONE

      CONTAINS

! --------------------------------------------------------------------------
      SUBROUTINE VRTQDR_SW(KLEV, KW, &
                           & PREF, PREFD, PTRA, PTRAD, &
                           & PDBT, PRDND, PRUP, PRUPD, PTDBT, &
                           & PFD, PFU)
! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006

!!!!!!!!MODIFICATIONS FOR UV PROCESSOR!!!!!!!
! ABozzo Dec 2014. Small changes to work as monochromatic UV code
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      INTEGER(KIND=JPIM), INTENT (IN) :: KLEV                   ! number of model layers
      INTEGER(KIND=JPIM), INTENT (IN) :: KW                     ! g-point index

      REAL(KIND=JPRB), INTENT(IN) :: PREF(:)                    ! direct beam reflectivity
                                                              !   Dimensions: (nlayers+1)
      REAL(KIND=JPRB), INTENT(IN) :: PREFD(:)                   ! diffuse beam reflectivity
                                                              !   Dimensions: (nlayers+1)
      REAL(KIND=JPRB), INTENT(IN) :: PTRA(:)                    ! direct beam transmissivity
                                                              !   Dimensions: (nlayers+1)
      REAL(KIND=JPRB), INTENT(IN) :: PTRAD(:)                   ! diffuse beam transmissivity
                                                              !   Dimensions: (nlayers+1)

      REAL(KIND=JPRB), INTENT(IN) :: PDBT(:)
                                                              !   Dimensions: (nlayers+1)
      REAL(KIND=JPRB), INTENT(IN) :: PTDBT(:)
                                                              !   Dimensions: (nlayers+1)

      REAL(KIND=JPRB), INTENT(INOUT) :: PRDND(:)
                                                              !   Dimensions: (nlayers+1)
      REAL(KIND=JPRB), INTENT(INOUT) :: PRUP(:)
                                                              !   Dimensions: (nlayers+1)
      REAL(KIND=JPRB), INTENT(INOUT) :: PRUPD(:)
                                                              !   Dimensions: (nlayers+1)

! Output
      REAL(KIND=JPRB), INTENT(OUT) :: PFD(:,:)                  ! downwelling flux (W/m2)
                                                              !   Dimensions: (nlayers+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
      REAL(KIND=JPRB), INTENT(OUT) :: PFU(:,:)                  ! upwelling flux (W/m2)
                                                              !   Dimensions: (nlayers+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle

! Local

      INTEGER(KIND=JPIM) :: IKP, IKX, JK

      REAL(KIND=JPRB) :: ZREFLECT
      REAL(KIND=JPRB) :: ZTDN(KLEV+1)  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
IF (LHOOK) CALL DR_HOOK('rrtmg_sw_vrtqdr',0,ZHOOK_HANDLE)

      ZREFLECT = 1._JPRB / (1._JPRB - PREFD(KLEV+1) * PREFD(KLEV))
      PRUP(KLEV) = PREF(KLEV) + (PTRAD(KLEV) * &
                 & ((PTRA(KLEV) - PDBT(KLEV)) * PREFD(KLEV+1) + &
                 &  PDBT(KLEV) * PREF(KLEV+1))) * ZREFLECT
      PRUPD(KLEV) = PREFD(KLEV) + PTRAD(KLEV) * PTRAD(KLEV) * &
                    & PREFD(KLEV+1) * ZREFLECT

! Pass from bottom to top 

      DO JK = 1,KLEV-1
         IKP = KLEV+1-JK                       
         IKX = IKP-1
         ZREFLECT = 1._JPRB / (1._JPRB -PRUPD(IKP) * PREFD(IKX))
         PRUP(IKX) = PREF(IKX) + (PTRAD(IKX) * &
                   & ((PTRA(IKX) - PDBT(IKX)) * PRUPD(IKP) + &
                   &  PDBT(IKX) * PRUP(IKP))) * ZREFLECT
         PRUPD(IKX) = PREFD(IKX) + PTRAD(IKX) * PTRAD(IKX) * &
                      & PRUPD(IKP) * ZREFLECT
      ENDDO
    
! Upper boundary conditions

      ZTDN(1) = 1._JPRB
      PRDND(1) = 0._JPRB
      ZTDN(2) = PTRA(1)
      PRDND(2) = PREFD(1)

! Pass from top to bottom

      DO JK = 2,KLEV
         IKP = JK+1
         ZREFLECT = 1._JPRB / (1._JPRB - PREFD(JK) * PRDND(JK))
         ZTDN(IKP) = PTDBT(JK) * PTRA(JK) + &
                    & (PTRAD(JK) * ((ZTDN(JK) - PTDBT(JK)) + &
                    & PTDBT(JK) * PREF(JK) * PRDND(JK))) * ZREFLECT
         PRDND(IKP) = PREFD(JK) + PTRAD(JK) * PTRAD(JK) * &
                      & PRDND(JK) * ZREFLECT
      ENDDO
    
! Up and down-welling fluxes at levels

      DO JK = 1,KLEV+1
         ZREFLECT = 1._JPRB / (1._JPRB - PRDND(JK) * PRUPD(JK))
         PFU(JK,KW) = (PTDBT(JK) * PRUP(JK) + &
                      & (ZTDN(JK) - PTDBT(JK)) * PRUPD(JK)) * ZREFLECT
         PFD(JK,KW) = PTDBT(JK) + (ZTDN(JK) - PTDBT(JK)+ &
                      & PTDBT(JK) * PRUP(JK) * PRDND(JK)) * ZREFLECT
      ENDDO
IF (LHOOK) CALL DR_HOOK('rrtmg_sw_vrtqdr',1,ZHOOK_HANDLE)
      END SUBROUTINE VRTQDR_SW

      END MODULE RRTMG_SW_VRTQDR
