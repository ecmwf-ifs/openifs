MODULE SUSOCEAN_ML_MOD

CONTAINS
SUBROUTINE SUSOCEAN_ML(LD_LEOCML,YDOCEAN_ML)

! (C) Copyright 2008- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! Purpose :
! -------
!   Setup for the ocean mixed layer model (KPP).

! Interface :
! ---------

! Method :
! ------

! Externals :

! ---------

! Reference :
! ---------
!   Large, W. G., J. C. McWilliams, and S. C. Doney (1994), Rev. Geophys.
!   Bernie, D. J., S. J. Woolnough, and J. M. Slingo (2005), J. Climate

! Modifications :
! -------------
!   07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
!---------------------------------------------------------------------


USE PARKIND1,     ONLY : JPIM, JPRB
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML

USE KPP_WSCALE_MOD

IMPLICIT NONE

TYPE(TOCEAN_ML), INTENT(INOUT) :: YDOCEAN_ML

LOGICAL :: LD_LEOCML

INTEGER(KIND=JPIM) :: ITMP(3)
REAL(KIND=JPRB) :: Z0          ! 0.0
REAL(KIND=JPRB) :: Z1          ! 1.0
REAL(KIND=JPRB) :: ZTMP(6) 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE 

LOGICAL :: LLTMP(2)

!#include "namocean_ml.h"

!---------------------------------------------------------------------

Z0= 0.0_JPRB
Z1= 1.0_JPRB

IF (LHOOK) CALL DR_HOOK('SUSOCEAN_ML_MOD:SUSOCEAN_ML',0,ZHOOK_HANDLE)
ASSOCIATE(LCDIAG_KPP=>YDOCEAN_ML%LCDIAG_KPP, LDD_KPP=>YDOCEAN_ML%LDD_KPP, &
 & LEOCML=>YDOCEAN_ML%LEOCML, LGEO_KPP=>YDOCEAN_ML%LGEO_KPP, &
 & LINIT_WSCALE=>YDOCEAN_ML%LINIT_WSCALE, LKPP_KPP=>YDOCEAN_ML%LKPP_KPP, &
 & LNBFLX_KPP=>YDOCEAN_ML%LNBFLX_KPP, LOCDEPTH_KPP=>YDOCEAN_ML%LOCDEPTH_KPP, &
 & LRI_KPP=>YDOCEAN_ML%LRI_KPP, LSF_NEW_KPP=>YDOCEAN_ML%LSF_NEW_KPP, &
 & LVL_NEW_KPP=>YDOCEAN_ML%LVL_NEW_KPP, NDIFFO=>YDOCEAN_ML%NDIFFO, &
 & NITERMAX_KPP=>YDOCEAN_ML%NITERMAX_KPP, NJERLOV_KPP=>YDOCEAN_ML%NJERLOV_KPP, &
 & NMRI_KPP=>YDOCEAN_ML%NMRI_KPP, NOCNSTEP_KPP=>YDOCEAN_ML%NOCNSTEP_KPP, &
 & NPD_KPP=>YDOCEAN_ML%NPD_KPP, NUTBLO=>YDOCEAN_ML%NUTBLO, &
 & NZTBLO=>YDOCEAN_ML%NZTBLO, RA1_JERLOV=>YDOCEAN_ML%RA1_JERLOV, &
 & RA2_JERLOV=>YDOCEAN_ML%RA2_JERLOV, RAM_KPP=>YDOCEAN_ML%RAM_KPP, &
 & RAS_KPP=>YDOCEAN_ML%RAS_KPP, RC1_KPP=>YDOCEAN_ML%RC1_KPP, &
 & RC2_KPP=>YDOCEAN_ML%RC2_KPP, RC3_KPP=>YDOCEAN_ML%RC3_KPP, &
 & RCEKMAN_KPP=>YDOCEAN_ML%RCEKMAN_KPP, RCMONOB_KPP=>YDOCEAN_ML%RCMONOB_KPP, &
 & RCM_KPP=>YDOCEAN_ML%RCM_KPP, RCSTAR_KPP=>YDOCEAN_ML%RCSTAR_KPP, &
 & RCS_KPP=>YDOCEAN_ML%RCS_KPP, RCV_KPP=>YDOCEAN_ML%RCV_KPP, &
 & RDEPTH_KPP=>YDOCEAN_ML%RDEPTH_KPP, RDIFM_BOT=>YDOCEAN_ML%RDIFM_BOT, &
 & RDIFM_IW=>YDOCEAN_ML%RDIFM_IW, RDIFM_MAX=>YDOCEAN_ML%RDIFM_MAX, &
 & RDIFS_BOT=>YDOCEAN_ML%RDIFS_BOT, RDIFS_IW=>YDOCEAN_ML%RDIFS_IW, &
 & RDIFS_MAX=>YDOCEAN_ML%RDIFS_MAX, RDMOL=>YDOCEAN_ML%RDMOL, &
 & RDSCALE_KPP=>YDOCEAN_ML%RDSCALE_KPP, RDSFMAX=>YDOCEAN_ML%RDSFMAX, &
 & REPSILON_KPP=>YDOCEAN_ML%REPSILON_KPP, REPS_KPP=>YDOCEAN_ML%REPS_KPP, &
 & RFAC_JERLOV=>YDOCEAN_ML%RFAC_JERLOV, RICR_KPP=>YDOCEAN_ML%RICR_KPP, &
 & RLAMBDA_KPP=>YDOCEAN_ML%RLAMBDA_KPP, RMIN_JERLOV=>YDOCEAN_ML%RMIN_JERLOV, &
 & RRAMSICE=>YDOCEAN_ML%RRAMSICE, RRIINFTY=>YDOCEAN_ML%RRIINFTY, &
 & RRRHO0=>YDOCEAN_ML%RRRHO0, RRRHO0_D06=>YDOCEAN_ML%RRRHO0_D06, &
 & RSICE=>YDOCEAN_ML%RSICE, RTICE=>YDOCEAN_ML%RTICE, &
 & RTOLFAC_KPP=>YDOCEAN_ML%RTOLFAC_KPP, RUMAX_KPP=>YDOCEAN_ML%RUMAX_KPP, &
 & RUMIN_KPP=>YDOCEAN_ML%RUMIN_KPP, RVC_KPP=>YDOCEAN_ML%RVC_KPP, &
 & RVDEL_KPP=>YDOCEAN_ML%RVDEL_KPP, RVD_KPP=>YDOCEAN_ML%RVD_KPP, &
 & RVE_KPP=>YDOCEAN_ML%RVE_KPP, RVONK=>YDOCEAN_ML%RVONK, RVTC=>YDOCEAN_ML%RVTC, &
 & RZETAM=>YDOCEAN_ML%RZETAM, RZETAS=>YDOCEAN_ML%RZETAS, &
 & RZETMAX_KPP=>YDOCEAN_ML%RZETMAX_KPP, RZETMIN_KPP=>YDOCEAN_ML%RZETMIN_KPP)


! 1. Define ocean mixed layer constants
!    ----------------------------------

NDIFFO = 3    ! 3 for UV, T, S 
NZTBLO = 892
NUTBLO = 482

REPS_KPP     = 1.0E-12_JPRB ! small number

! bldepth
!RICR_KPP     = 0.3_JPRB
RICR_KPP     = 0.4_JPRB ! ROMS,0.45 
REPSILON_KPP = 0.1_JPRB
RCEKMAN_KPP  = 0.7_JPRB
RCMONOB_KPP  = 1.0_JPRB
!RCV_KPP      = 1.6_JPRB
RCV_KPP      = 1.8_JPRB  ! 1 < RCV_KPP < 2, see HYCOM manual.
RVONK        = 0.4_JPRB
RCS_KPP      = 98.96_JPRB
RVTC =  RCV_KPP * SQRT( 0.2_JPRB / RCS_KPP / REPSILON_KPP ) &
     &  / RVONK**2 / RICR_KPP

! blmix
RCSTAR_KPP   = 5.0_JPRB   ! See Bernie et al. (2005)  10.0 -> 5.0
!RCSTAR_KPP   = 10.0_JPRB  ! See Bernie et al. (2005)  10.0 -> 5.0
!RCSTAR_KPP   = 15.8_JPRB  ! Large et al. 1997

! ddmix
RRRHO0      = 1.9_JPRB        
RRRHO0_D06  = 2.55_JPRB   ! Danabasoglu et al. (2006)        
RDSFMAX     = 1.0E-4_JPRB ! Changed from 1.0E-3 in LMD94
!RDSFMAX     = 1.0E-3_JPRB   
RDMOL       = 1.5E-6_JPRB    

! ri_iwmix
RRIINFTY   = 0.8_JPRB     ! Changed from 0.7 in LMD94
RDIFM_MAX  = 5.0E-3_JPRB
RDIFS_MAX  = 5.0E-3_JPRB
!RDIFM_IW   = 1.0E-4_JPRB ! Changed from 1.0E-5 in LMD94  
!RDIFM_IW   = 1.0E-5_JPRB
RDIFM_IW   = 1.5E-4_JPRB 
!RDIFS_IW   = 1.0E-5_JPRB 
RDIFS_IW   = 1.5E-5_JPRB 

!Jerlov water type :      I         Ia         Ib          II        III
!      NJERLOV_KPP :      1          2          3           4          5
RFAC_JERLOV =  (/ 0.58_JPRB, 0.62_JPRB, 0.67_JPRB,  0.77_JPRB, 0.78_JPRB /) 
RA1_JERLOV  =  (/ 0.35_JPRB,  0.6_JPRB,  1.0_JPRB,   1.5_JPRB, 1.4_JPRB  /)  
RA2_JERLOV  =  (/ 23.0_JPRB, 20.0_JPRB, 17.0_JPRB,  14.0_JPRB, 7.9_JPRB  /)  
RMIN_JERLOV = -100.0_JPRB

! vmix
RDIFM_BOT = 1.0E-4_JPRB
RDIFS_BOT = 1.0E-5_JPRB
RSICE     = 4.0_JPRB 
RTICE     = -1.8_JPRB
RRAMSICE  = 4.32E5_JPRB !5 days, parameter for sea-ice flux [s]

! wscale
RZETMIN_KPP = -4.0E-7_JPRB ! m3/s3
RZETMAX_KPP = 0.0_JPRB          
RUMIN_KPP   = 0.0_JPRB     ! m/s
RUMAX_KPP   = 0.04_JPRB

RC1_KPP =  5.0_JPRB
RC2_KPP = 16.0_JPRB
RC3_KPP = 16.0_JPRB
RZETAM  = -0.2_JPRB
RZETAS  = -1.0_JPRB
RAM_KPP = 1.257_JPRB
RAS_KPP = -28.86_JPRB
RCM_KPP = 8.380_JPRB

! implicit time integration
NITERMAX_KPP = 20          
RLAMBDA_KPP = 0.5_JPRB     ! parameter for Euler-backward scheme) ( 0 =< RLAMBDA_KPP <= 1
RTOLFAC_KPP = 0.5_JPRB     ! factor of the criterior for convergence 

! 2. Initialize ocean mixed layer settings
!    -------------------------------------

LEOCML = LD_LEOCML         ! Ocean mixed layer on

LCDIAG_KPP = .FALSE.       ! .TRUE.  : compute diagnostic variables
LOCDEPTH_KPP  = .TRUE.     ! .TRUE.  : ocean topography
                           ! .FALSE. : fixed depth

RDEPTH_KPP    = 250.0_JPRB ! used for fixed ocean depth (>0) [m]
RDSCALE_KPP   = 10.0E0_JPRB ! scale parameter, = -lambda in LMD94 Appendix D. 
                           ! > 0.0 for stretch grid
                           ! < 0.0 for constant grid 
!RDSCALE_KPP    = -1.0E0_JPRB !constant grid interval

LVL_NEW_KPP = .TRUE.       ! new vertical level definition
!RVC_KPP     = 15.0_JPRB    ! PDO = ZA*ZX1 + ZB*(ZX2**RVC_KPP) + RVD_KPP*(ZX2**RVE_KPP)
!RVD_KPP     = 35.0_JPRB
!RVE_KPP     = 2.0_JPRB                          
!RVDEL_KPP   = 0.6_JPRB     ! thickness of the top layer

LGEO_KPP     = .TRUE.     ! geostrophic flow decomposition

RVC_KPP     = 8.0_JPRB    ! PDO = ZA*ZX1 + ZB*(ZX2**RVC_KPP) + RVD_KPP*(ZX2**RVE_KPP)
RVD_KPP     = 10.0_JPRB
RVE_KPP     = 2.0_JPRB                          
RVDEL_KPP   = 0.25_JPRB     ! thickness of the top layer

NOCNSTEP_KPP = 1           ! The number of the time step. A time step of AGCM 
                           ! is smaller than KPP needs, so usually set to 1.
NJERLOV_KPP  = 1           ! Jerlov water type
                           ! 1: I,  2: Ia,  3: Ib,  4: II,  5:III
NMRI_KPP = 1               ! number of vertical smoothing in bulk Ri number
NPD_KPP  = 1               ! no. of layers distributed surface fluxes

!vmix
!LDD_KPP    = .TRUE.        ! double diffusion
LDD_KPP    = .FALSE.       ! double diffusion
LSF_NEW_KPP= .TRUE.        ! new salt fingering, Danabasoglu et al. 2006 JC
LKPP_KPP   = .TRUE.        ! KPP non-local mixing 
LNBFLX_KPP = .TRUE.        ! .TRUE. -> diff.=0.0 at the bottom 
LRI_KPP    = .TRUE.        ! Richardson instability mixing 

! tables
LINIT_WSCALE  = .TRUE.     ! flag for universal function 


! 3. Allocate variables
!    ------------------
IF(.NOT.ALLOCATED(YDOCEAN_ML%RWMT)) ALLOCATE(YDOCEAN_ML%RWMT(0:NZTBLO+1,0:NUTBLO+1))
IF(.NOT.ALLOCATED(YDOCEAN_ML%RWST)) ALLOCATE(YDOCEAN_ML%RWST(0:NZTBLO+1,0:NUTBLO+1))

! 4. Read namelist
!    -------------

!READ(NULNAM,NAMOML)
!WRITE(NULOUT,NAMOML)

! 5. Initialize the table for universal function
!    -------------------------------------------
ITMP = 1
ZTMP = Z1
LLTMP = .TRUE.
CALL KPP_WSCALE &
 & ( ITMP(1)  ,ITMP(2)  ,ITMP(3)  ,LLTMP    ,ZTMP(1)  ,&
 &   ZTMP(2)  ,ZTMP(3)  ,ZTMP(4)  ,ZTMP(5)  ,ZTMP(6)  ,& ! All dummy arguments.
 &   YDOCEAN_ML )
LINIT_WSCALE = .FALSE.

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSOCEAN_ML_MOD:SUSOCEAN_ML',1,ZHOOK_HANDLE)

!---------------------------------------------------------------------

END SUBROUTINE SUSOCEAN_ML
END MODULE SUSOCEAN_ML_MOD
