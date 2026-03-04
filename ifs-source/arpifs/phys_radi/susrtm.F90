! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUSRTM

!     Adapted from E.J. Mlawer, J. Delamere, Atmospheric & Environmental Research.
!     by JJMorcrette, ECMWF
!     Modified to add arrays relevant to mapping for g-point reduction,
!     M.J. Iacono, Atmospheric & Environmental Research, Inc. 
!     JJMorcrette 20010610 Flexible configuration for number of g-points
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB ,   JPIM
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOESRTM  , ONLY : JPGPT, NGBSW, NGN
USE YOESRTWN , ONLY : NG      , NSPA, NSPB   , NMPSRTM, &
 & PREF    , PREFLOG , TREF   , &
 & NGM     , WT      , NGC    , NGS
! & NGM     , WT      , NGC    , NGS , NGN    , NGBSW
! & WAVENUM1, WAVENUM2, DELWAVE, PREF, PREFLOG, TREF   , &

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC56(14), IGC112(14) , IGC224(14)
INTEGER(KIND=JPIM) :: IGS56(14), IGS112(14) , IGS224(14)

INTEGER(KIND=JPIM) :: IGM56(224),IGM112(224), IGM224(224)

INTEGER(KIND=JPIM) :: IGN56(56), IGN112(112), IGN224(224)
INTEGER(KIND=JPIM) :: IGB56(56), IGB112(112), IGB224(224)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSRTM',0,ZHOOK_HANDLE)

NG(:)     =(/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16 /)
NSPA(:)   =(/  9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1 /)
NSPB(:)   =(/  1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1 /)
NMPSRTM(:)=(/  6, 6, 5, 5, 5, 5, 5, 4, 4, 3, 2, 2, 1, 6 /)

!WAVENUM1( :) = (/&
! & 2600._JPRB, 3250._JPRB, 4000._JPRB, 4650._JPRB, 5150._JPRB, 6150._JPRB, 7700._JPRB &
! & , 8050._JPRB,12850._JPRB,16000._JPRB,22650._JPRB,29000._JPRB,38000._JPRB,  820._JPRB /)  
!WAVENUM2( :) = (/&
! & 3250._JPRB, 4000._JPRB, 4650._JPRB, 5150._JPRB, 6150._JPRB, 7700._JPRB, 8050._JPRB &
! & ,12850._JPRB,16000._JPRB,22650._JPRB,29000._JPRB,38000._JPRB,50000._JPRB, 2600._JPRB /)  
!DELWAVE( :) = (/&
! & 650._JPRB,  750._JPRB,  650._JPRB,  500._JPRB, 1000._JPRB, 1550._JPRB,  350._JPRB &
! & , 4800._JPRB, 3150._JPRB, 6650._JPRB, 6350._JPRB, 9000._JPRB,12000._JPRB, 1780._JPRB /)  

!=====================================================================
! Set arrays needed for the g-point reduction from 224 to 
! - either 112 for the high-resolution forecast model configuration
! - or 56 for the EPS-type configuration  
! in the 14 SW bands:

! NB: This mapping from 224 to 112 points has been carefully selected to
! minimize the effect on the resulting fluxes and cooling rates, and
! caution should be used if the mapping is modified.
!     The further reduction to 56 for EPS configuration is considered 
! acceptable, only because of the random perturbations introduced on 
! the total heating rates produced by the physical parametrization package.
! While a reduction to 56 obviously speeds up the model, it as obviously 
! reduces the accuracy that could be expected from the radiation scheme.

! JPGPT   The total number of new g-points (NGPT)
! NGC     The number of new g-points in each band (14)
! NGS     The cumulative sum of new g-points for each band (14)
! NGM     The index of each new g-point relative to the original
!         16 g-points for each band.
! NGN     The number of original g-points that are combined to make
!         each new g-point in each band.
! NGB     The band index for each new g-point.
! WT      RRTM weights for 16 g-points. (16)

!-- ECMWF EPS model RRTM_SW configuration with 56 g-points
IGC56(:) = (/ 3, 6, 4, 4, 5, 5, 1, 5, 4, 3, 3, 4, 3, 6 /)
IGS56(:) = (/ 3, 9,13,17,22,27,28,33,37,40,43,47, 50, 56 /)

IGM56(:) = (/ 1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3, &            ! Band 16
            & 1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6, &            ! Band 17
            & 1,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4, &            ! Band 18
            & 1,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4, &            ! Band 19
            & 1,1,2,2,3,3,4,4,5,5,5,5,5,5,5,5, &            ! Band 20
            & 1,1,2,2,3,3,4,4,5,5,5,5,5,5,5,5, &            ! Band 21
            & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &            ! Band 22
            & 1,1,1,1,2,2,3,3,4,4,5,5,5,5,5,5, &            ! Band 23
            & 1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4, &            ! Band 24
            & 1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3, &            ! Band 25
            & 1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3, &            ! Band 26
            & 1,1,2,2,3,3,4,4,4,4,4,4,4,4,4,4, &            ! Band 27
            & 1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3, &            ! Band 28
            & 1,1,2,2,3,3,3,3,4,4,4,4,5,5,6,6 /)            ! Band 29

IGN56(:) = (/ 4,4,8, &                                     ! Band 16
            & 2,2,3,3,3,3, &                               ! Band 17
            & 2,2,4,8, &                                   ! Band 18
            & 2,2,4,8, &                                   ! Band 19
            & 2,2,2,2,8, &                                 ! Band 20
            & 2,2,2,2,8, &                                 ! Band 21
            & 16, &                                        ! Band 22
            & 4,2,2,2,6, &                                 ! Band 23
            & 4,4,4,4, &                                   ! Band 24
            & 2,4,10, &                                    ! Band 25
            & 2,4,10, &                                    ! Band 26
            & 2,2,2,10, &                                  ! Band 27
            & 2,4,10, &                                    ! Band 28
            & 2,2,4,4,2,2 /)                               ! Band 29

IGB56(:) = (/ 16,16,16, &                                  ! Band 16
            & 17,17,17,17,17,17, &                         ! Band 17
            & 18,18,18,18, &                               ! Band 18
            & 19,19,19,19, &                               ! Band 19
            & 20,20,20,20,20, &                            ! Band 20
            & 21,21,21,21,21, &                            ! Band 21
            & 22, &                                        ! Band 22
            & 23,23,23,23,23, &                            ! Band 23
            & 24,24,24,24, &                               ! Band 24
            & 25,25,25, &                                  ! Band 25
            & 26,26,26, &                                  ! Band 26
            & 27,27,27,27, &                               ! Band 27
            & 28,28,28, &                                  ! Band 28
            & 29,29,29,29,29,29 /)                         ! Band 29

!-------------------------------------------------------------------------------
!-- ECMWF high-resolution model RRTM_SW configuration with 112 g-points
! Use this NGC, NGS, NGM, and NGN for reduced (112) g-point set
! (A related code change is required in modules parsrtm.F90 and yoesrtwn.F90)

IGC112(:) = (/ 6,12, 8, 8,10,10, 2,10, 8, 6, 6, 8, 6,12 /)
IGS112(:) = (/ 6,18,26,34,44,54,56,66,74,80,86,94,100,112 /)

!NGM(:)
IGM112(:) = (/ 1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! Band 16
             & 1,2,3,4,5,6,6,7,8,8,9,10,10,11,12,12, &      ! Band 17
             & 1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! Band 18
             & 1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! Band 19
             & 1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! Band 20
             & 1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! Band 21
             & 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! Band 22
             & 1,1,2,2,3,4,5,6,7,8,9,9,10,10,10,10, &       ! Band 23
             & 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! Band 24
             & 1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! Band 25
             & 1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! Band 26
             & 1,2,3,4,5,6,7,7,7,7,8,8,8,8,8,8, &           ! Band 27
             & 1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! Band 28
             & 1,2,3,4,5,5,6,6,7,7,8,8,9,10,11,12 /)        ! Band 29
!NGN(:)
IGN112(:) = (/ 2,2,2,2,4,4, &                               ! Band 16
             & 1,1,1,1,1,2,1,2,1,2,1,2, &                   ! Band 17
             & 1,1,1,1,2,2,4,4, &                           ! Band 18
             & 1,1,1,1,2,2,4,4, &                           ! Band 19
             & 1,1,1,1,1,1,1,1,2,6, &                       ! Band 20
             & 1,1,1,1,1,1,1,1,2,6, &                       ! Band 21
             & 8,8, &                                       ! Band 22
             & 2,2,1,1,1,1,1,1,2,4, &                       ! Band 23
             & 2,2,2,2,2,2,2,2, &                           ! Band 24
             & 1,1,2,2,4,6, &                               ! Band 25
             & 1,1,2,2,4,6, &                               ! Band 26
             & 1,1,1,1,1,1,4,6, &                           ! Band 27
             & 1,1,2,2,4,6, &                               ! Band 28
             & 1,1,1,1,2,2,2,2,1,1,1,1 /)                   ! Band 29
!NGBSW(:)
IGB112(:) = (/ 16,16,16,16,16,16, &                         ! Band 16
             & 17,17,17,17,17,17,17,17,17,17,17,17, &       ! Band 17
             & 18,18,18,18,18,18,18,18, &                   ! Band 18
             & 19,19,19,19,19,19,19,19, &                   ! Band 19
             & 20,20,20,20,20,20,20,20,20,20, &             ! Band 20
             & 21,21,21,21,21,21,21,21,21,21, &             ! Band 21
             & 22,22, &                                     ! Band 22
             & 23,23,23,23,23,23,23,23,23,23, &             ! Band 23
             & 24,24,24,24,24,24,24,24, &                   ! Band 24
             & 25,25,25,25,25,25, &                         ! Band 25
             & 26,26,26,26,26,26, &                         ! Band 26
             & 27,27,27,27,27,27,27,27, &                   ! Band 27
             & 28,28,28,28,28,28, &                         ! Band 28
             & 29,29,29,29,29,29,29,29,29,29,29,29 /)       ! Band 29

!-------------------------------------------------------------------------------
!-- original RRTM_SW configuration with 224 (14*16 g-points)
! Use this NGC, NGS, NGM, and NGN for full (224) g-point set
! (A related code change is required in modules parsrtm.F90 and yoesrtwn.F90)
IGC224(:) = (/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16 /)
IGS224(:) = (/ 16,32,48,64,80,96,112,128,144,160,176,192,208,224 /)

IGM224(:) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 16
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 17
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 18
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 19
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 20
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 21
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 22
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 23
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 24
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 25
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 26
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 27
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! Band 28
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 /)    ! Band 29

IGN224(:) = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 16
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 17
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 18
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 19
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 20
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 21
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 22
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 23
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 24
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 25
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 26
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 27
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! Band 28
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /)           ! Band 29

IGB224(:) = (/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16, &   ! Band 16
             & 17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17, &   ! Band 17
             & 18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18, &   ! Band 18
             & 19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19, &   ! Band 19
             & 20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20, &   ! Band 20
             & 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21, &   ! Band 21
             & 22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22, &   ! Band 22
             & 23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23, &   ! Band 23
             & 24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24, &   ! Band 24
             & 25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25, &   ! Band 25
             & 26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26, &   ! Band 26
             & 27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27, &   ! Band 27
             & 28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28, &   ! Band 28
             & 29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,28 /)   ! Band 29

!=============================================================================

WT(:) =  (/ 0.1527534276_JPRB, 0.1491729617_JPRB, 0.1420961469_JPRB, &
          & 0.1316886544_JPRB, 0.1181945205_JPRB, 0.1019300893_JPRB, &
          & 0.0832767040_JPRB, 0.0626720116_JPRB, 0.0424925000_JPRB, &
          & 0.0046269894_JPRB, 0.0038279891_JPRB, 0.0030260086_JPRB, &
          & 0.0022199750_JPRB, 0.0014140010_JPRB, 0.0005330000_JPRB, &
          & 0.0000750000_JPRB /)

!=============================================================================

! These pressures are chosen such that the ln of the first pressure
! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
!  each subsequent ln(pressure) differs from the previous one by 0.2.
PREF = (/ &
 & 1.05363E+03_JPRB,8.62642E+02_JPRB,7.06272E+02_JPRB,5.78246E+02_JPRB,4.73428E+02_JPRB, &
 & 3.87610E+02_JPRB,3.17348E+02_JPRB,2.59823E+02_JPRB,2.12725E+02_JPRB,1.74164E+02_JPRB, &
 & 1.42594E+02_JPRB,1.16746E+02_JPRB,9.55835E+01_JPRB,7.82571E+01_JPRB,6.40715E+01_JPRB, &
 & 5.24573E+01_JPRB,4.29484E+01_JPRB,3.51632E+01_JPRB,2.87892E+01_JPRB,2.35706E+01_JPRB, &
 & 1.92980E+01_JPRB,1.57998E+01_JPRB,1.29358E+01_JPRB,1.05910E+01_JPRB,8.67114E+00_JPRB, &
 & 7.09933E+00_JPRB,5.81244E+00_JPRB,4.75882E+00_JPRB,3.89619E+00_JPRB,3.18993E+00_JPRB, &
 & 2.61170E+00_JPRB,2.13828E+00_JPRB,1.75067E+00_JPRB,1.43333E+00_JPRB,1.17351E+00_JPRB, &
 & 9.60789E-01_JPRB,7.86628E-01_JPRB,6.44036E-01_JPRB,5.27292E-01_JPRB,4.31710E-01_JPRB, &
 & 3.53455E-01_JPRB,2.89384E-01_JPRB,2.36928E-01_JPRB,1.93980E-01_JPRB,1.58817E-01_JPRB, &
 & 1.30029E-01_JPRB,1.06458E-01_JPRB,8.71608E-02_JPRB,7.13612E-02_JPRB,5.84256E-02_JPRB, &
 & 4.78349E-02_JPRB,3.91639E-02_JPRB,3.20647E-02_JPRB,2.62523E-02_JPRB,2.14936E-02_JPRB, &
 & 1.75975E-02_JPRB,1.44076E-02_JPRB,1.17959E-02_JPRB,9.65769E-03_JPRB /)  
PREFLOG = (/ &
 & 6.9600E+00_JPRB, 6.7600E+00_JPRB, 6.5600E+00_JPRB, 6.3600E+00_JPRB, 6.1600E+00_JPRB, &
 & 5.9600E+00_JPRB, 5.7600E+00_JPRB, 5.5600E+00_JPRB, 5.3600E+00_JPRB, 5.1600E+00_JPRB, &
 & 4.9600E+00_JPRB, 4.7600E+00_JPRB, 4.5600E+00_JPRB, 4.3600E+00_JPRB, 4.1600E+00_JPRB, &
 & 3.9600E+00_JPRB, 3.7600E+00_JPRB, 3.5600E+00_JPRB, 3.3600E+00_JPRB, 3.1600E+00_JPRB, &
 & 2.9600E+00_JPRB, 2.7600E+00_JPRB, 2.5600E+00_JPRB, 2.3600E+00_JPRB, 2.1600E+00_JPRB, &
 & 1.9600E+00_JPRB, 1.7600E+00_JPRB, 1.5600E+00_JPRB, 1.3600E+00_JPRB, 1.1600E+00_JPRB, &
 & 9.6000E-01_JPRB, 7.6000E-01_JPRB, 5.6000E-01_JPRB, 3.6000E-01_JPRB, 1.6000E-01_JPRB, &
 & -4.0000E-02_JPRB,-2.4000E-01_JPRB,-4.4000E-01_JPRB,-6.4000E-01_JPRB,-8.4000E-01_JPRB, &
 & -1.0400E+00_JPRB,-1.2400E+00_JPRB,-1.4400E+00_JPRB,-1.6400E+00_JPRB,-1.8400E+00_JPRB, &
 & -2.0400E+00_JPRB,-2.2400E+00_JPRB,-2.4400E+00_JPRB,-2.6400E+00_JPRB,-2.8400E+00_JPRB, &
 & -3.0400E+00_JPRB,-3.2400E+00_JPRB,-3.4400E+00_JPRB,-3.6400E+00_JPRB,-3.8400E+00_JPRB, &
 & -4.0400E+00_JPRB,-4.2400E+00_JPRB,-4.4400E+00_JPRB,-4.6400E+00_JPRB /)  
! These are the temperatures associated with the respective 
! pressures for the MLS standard atmosphere. 
TREF = (/ &
 & 2.9420E+02_JPRB, 2.8799E+02_JPRB, 2.7894E+02_JPRB, 2.6925E+02_JPRB, 2.5983E+02_JPRB, &
 & 2.5017E+02_JPRB, 2.4077E+02_JPRB, 2.3179E+02_JPRB, 2.2306E+02_JPRB, 2.1578E+02_JPRB, &
 & 2.1570E+02_JPRB, 2.1570E+02_JPRB, 2.1570E+02_JPRB, 2.1706E+02_JPRB, 2.1858E+02_JPRB, &
 & 2.2018E+02_JPRB, 2.2174E+02_JPRB, 2.2328E+02_JPRB, 2.2479E+02_JPRB, 2.2655E+02_JPRB, &
 & 2.2834E+02_JPRB, 2.3113E+02_JPRB, 2.3401E+02_JPRB, 2.3703E+02_JPRB, 2.4022E+02_JPRB, &
 & 2.4371E+02_JPRB, 2.4726E+02_JPRB, 2.5085E+02_JPRB, 2.5457E+02_JPRB, 2.5832E+02_JPRB, &
 & 2.6216E+02_JPRB, 2.6606E+02_JPRB, 2.6999E+02_JPRB, 2.7340E+02_JPRB, 2.7536E+02_JPRB, &
 & 2.7568E+02_JPRB, 2.7372E+02_JPRB, 2.7163E+02_JPRB, 2.6955E+02_JPRB, 2.6593E+02_JPRB, &
 & 2.6211E+02_JPRB, 2.5828E+02_JPRB, 2.5360E+02_JPRB, 2.4854E+02_JPRB, 2.4348E+02_JPRB, &
 & 2.3809E+02_JPRB, 2.3206E+02_JPRB, 2.2603E+02_JPRB, 2.2000E+02_JPRB, 2.1435E+02_JPRB, &
 & 2.0887E+02_JPRB, 2.0340E+02_JPRB, 1.9792E+02_JPRB, 1.9290E+02_JPRB, 1.8809E+02_JPRB, &
 & 1.8329E+02_JPRB, 1.7849E+02_JPRB, 1.7394E+02_JPRB, 1.7212E+02_JPRB /)  
!     -----------------------------------------------------------------

IF (JPGPT == 56) THEN

!- 14
  NGC(:)=IGC56(:)
  NGS(:)=IGS56(:)
!- 14*16=224
  NGM(:)=IGM56(:)

  NGN(1:56)=IGN56(1:56)
  NGBSW(1:56)=IGB56(1:56)

ELSEIF (JPGPT == 112) THEN
!- 14
  NGC(:)=IGC112(:)
  NGS(:)=IGS112(:)
!- 14*16=224
  NGM(:)=IGM112(:)

  NGN(1:112)=IGN112(1:112)
  NGBSW(1:112)=IGB112(1:112)

ELSEIF (JPGPT == 224) THEN
!- 14
  NGC(:)=IGC224(:)
  NGS(:)=IGS224(:)
!- 14*16=224
  NGM(:)=IGM224(:)

  NGN(1:224)=IGN224(1:224)
  NGBSW(1:224)=IGB224(1:224)

ENDIF

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSRTM',1,ZHOOK_HANDLE)
END SUBROUTINE SUSRTM

