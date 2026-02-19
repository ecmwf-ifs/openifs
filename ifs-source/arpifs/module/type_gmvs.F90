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

MODULE TYPE_GMVS

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE
SAVE

!--------------------------------------------------------------------
! Type declarations used for GMV (see yomgmv.F90)

! Modifications:
! N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
! K. Yessad (Sep 2008): remove lpc_xidt
! K. Yessad (Jan 2011): add attributes MCSPDPT and MCSVDPT.
! K. Yessad (June 2017): Introduce NHQE model.
! K. Yessad (Feb 2018): remove deep-layer formulations.
! End Modifications
!--------------------------------------------------------------------

TYPE TYPE_T0
! Pointers to arrays at time t and some dimensions
INTEGER(KIND=JPIM) :: NDIM   = 0
INTEGER(KIND=JPIM) :: NDIMS  = 0
INTEGER(KIND=JPIM) :: NDIMUV = 0
! Field pointers in GMV
INTEGER(KIND=JPIM) :: MU     = - HUGE(1_JPIM) 
INTEGER(KIND=JPIM) :: MV     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MT     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MTL    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MTM    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MDIV   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVOR   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MUL    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVL    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPD   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPDL  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPDM  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVD   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVDL  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVDM  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MNHX   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MNHXL  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MNHXM  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MEDOT  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSGRTL = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSGRTM = - HUGE(1_JPIM)
! Field pointers in GMVS
INTEGER(KIND=JPIM) :: MSP    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPL   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPM   = - HUGE(1_JPIM)
END TYPE TYPE_T0
!-------------------------------------------------------------------------
TYPE TYPE_T9
! Pointers to arrays at time t-dt and some dimensions
INTEGER(KIND=JPIM) :: NDIM       = 0
INTEGER(KIND=JPIM) :: NDIMS      = 0
! Field pointers in GMV
INTEGER(KIND=JPIM) :: MU         = - HUGE(1_JPIM) 
INTEGER(KIND=JPIM) :: MV         = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MT         = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MTL        = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MTM        = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MEDOT      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MUNL       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MUNL_SI    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVNL       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVNL_SI    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MTNL       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MTNL_SI    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPNL      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPNL_SI   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVWVNL     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVDNL_SI  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPDNL     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPDNL_SI  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MDIV       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVOR       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MUL        = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVL        = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPD       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPDL      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPDM      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVD       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVDL      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVDM      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MNHX       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MGW        = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCURHS     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCVRHS     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCTRHS     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSPRHS    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSPDRHS   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSVDRHS   = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCUNL      = - HUGE(1_JPIM) 
INTEGER(KIND=JPIM) :: MCVNL      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCTNL      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSPNL     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCVWVNL    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSPDNL    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCUPT      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCVPT      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCTPT      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSPDPT    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSVDPT    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MDPHI      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MNHXNL     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCNHXNL    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSGRTL     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSGRTM     = - HUGE(1_JPIM)
! Field pointers in GMVS
INTEGER(KIND=JPIM) :: MSP        = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPL       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPM       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSPPT     = - HUGE(1_JPIM) 
INTEGER(KIND=JPIM) :: MDBBC      = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MGWS       = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPNL2     = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MCSPNL2    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MPREHYDS   = - HUGE(1_JPIM)

END TYPE TYPE_T9
!-------------------------------------------------------------------------
TYPE TYPE_T1
! Pointers to arrays at time t+dt and some dimensions
INTEGER(KIND=JPIM) :: NDIM  = 0
INTEGER(KIND=JPIM) :: NDIMS = 0
! Field pointers in GMVT1
INTEGER(KIND=JPIM) :: MU    = - HUGE(1_JPIM) 
INTEGER(KIND=JPIM) :: MV    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MT    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSPD  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MSVD  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MVOR  = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MDIV  = - HUGE(1_JPIM)
! MNHX    - pointer to the term "X" at time "t+dt" 
INTEGER(KIND=JPIM) :: MNHX  = - HUGE(1_JPIM)
! Field pointers in GMVT1S
INTEGER(KIND=JPIM) :: MSP   = - HUGE(1_JPIM)
END TYPE TYPE_T1
!-------------------------------------------------------------------------
TYPE TYPE_PH9
! Pointers to arrays at time "t-dt" for ECMWF physics
! Field pointers in GMV
INTEGER(KIND=JPIM) :: MU    = - HUGE(1_JPIM) 
INTEGER(KIND=JPIM) :: MV    = - HUGE(1_JPIM)
INTEGER(KIND=JPIM) :: MT    = - HUGE(1_JPIM)
! Field pointers in GMVS
INTEGER(KIND=JPIM) :: MSP   = - HUGE(1_JPIM)
END TYPE TYPE_PH9
!-------------------------------------------------------------------------
TYPE TYPE_GP
INTEGER(KIND=JPIM) :: NDIM  = 0              ! Number of grid-point GMV fields
! Pointers (offsets) for individual fields in grid-point arrays
INTEGER(KIND=JPIM) :: MU    = -HUGE(1_JPIM)  ! U-wind
INTEGER(KIND=JPIM) :: MV    = -HUGE(1_JPIM)  ! V-wind
INTEGER(KIND=JPIM) :: MT    = -HUGE(1_JPIM)  ! Temperature
INTEGER(KIND=JPIM) :: MSPD  = -HUGE(1_JPIM)  ! Pressure departure variable
INTEGER(KIND=JPIM) :: MSVD  = -HUGE(1_JPIM)  ! Vertical div or velocity variable
INTEGER(KIND=JPIM) :: MNHX  = -HUGE(1_JPIM)  ! NHX term
END TYPE TYPE_GP
!-------------------------------------------------------------------------
END MODULE TYPE_GMVS
