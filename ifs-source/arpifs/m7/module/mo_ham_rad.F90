!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_rad.f90
!!
!! \brief
!! mo_ham_rad holds the parameters and routines for
!! the calculation of the optical parameters
!! for the aerosol distribution simulated
!! by the ECHAM/HAM aerosol module.
!!
!! \author Olivier Boucher (Univ. Lille)
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! [ John Doe, john.doe@blabla.com -Compulsory- ]
!!
!! \revision_history
!!   -# O. Boucher (Univ. Lille) - original csource of some radiation routines (2003)
!!   -# P. Stier (MPI-Met) - (2003-05-08)
!!   -# P. Stier (Caltech) - adaption to ECHAM5/HAM, additional routines, refractive indices (2007)
!!   -# S. Rast (MPI-Met) - adaptation to echam5.3.2 (2007-03)
!!   -# D. O'Donnell (MPI-Met) - SOA (XXXX)
!!   -# K. Zhang (MPI-Met) - submodel interface (2009-07)
!!   -# M.G. Schultz (FZ Juelich) - cleanup (XXXX)
!!   -# P. Stier (Uni Oxford) - adaptation to RRTM-SW (2010)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!   -#  H. Kokkola (FMI) - modified to include SALSA (2013-06)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!!   - Stier et al., ACP,  2005 (SW, volume weighted refractive indices)
!!   - Stier et al., ACPD, 2007 (LW, mixing rules,...)
!!
!! \belongs_to
!!  HAMMOZ
!!
!!
!! \belongs_to
!!  HAMMOZ
!!
!!  SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_rad

  USE mo_ham_rad_data, ONLY: Nwv_sw, Nwv_sw_opt, Nwv_lw,        &
                              Nwv_tot, Nwv_sw_tot,              &
                              cnr, cni,                         &
                              log_x0_min, log_x0_max,           &
                              log_ni_min,                       &
                              x0_min, x0_max,                   &
                              nnrmax, nnimax, ndismax,          &
                              nr_min, nr_max, ni_min, ni_max,   &
                              inc_nr, inc_ni,                   &
                              lambda, lambda_sw_opt,            &
                              ham_rad_data_initialize
  USE mo_ham,           ONLY: naerocomp,           &
                              aerocomp,          &
                              aerowater,         &
                              nraddiag, nrad,    &
                              sizeclass,            &
                              nclass
  USE mo_ham,           ONLY: subm_aerospec 
  USE mo_kind,          ONLY: dp, THRESHOLD
  USE mo_species,       ONLY: speclist, naerospec, nmaxspec
  !>>dod soa
  USE mo_ham_species,   ONLY: id_oc, id_wat    !!mgs!!   , naerospec=>ham_naerospec, aerospec=>ham_aerospec
  !<<dod
!#ifdef _OPENMP
!    use omp_lib
!#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC ham_rad,                  & !
#ifdef HAMMOZ
         ham_rad_diag,             & !
#endif
         ham_rad_cache,            & !
         ham_rad_cache_cleanup,    & ! 
         ham_rad_mem,              & !
         ham_rad_mem_cleanup,      & !
         ham_rad_initialize

  !--- Optical properties at GCM wavebands:

  !REAL(dp), PUBLIC, ALLOCATABLE :: sigma(:,:,:,:)
  !REAL(dp), PUBLIC, ALLOCATABLE :: omega(:,:,:,:)
  !REAL(dp), PUBLIC, ALLOCATABLE :: asym(:,:,:,:)
  !REAL(dp), PUBLIC, ALLOCATABLE :: nr(:,:,:,:)
  !REAL(dp), PUBLIC, ALLOCATABLE :: ni(:,:,:,:)


  INTEGER, PUBLIC    :: nhamaer=0            !--Aerosol modes potentially 
                                             !  considered in radiation

  INTEGER :: aero_ridx(nmaxspec)             ! reverse index from speclist to naerospec
 
  !--- Look-up tables:

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: lut1_sigma, lut1_g, lut1_omega, lut1_pp180
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: lut2_sigma, lut2_g, lut2_omega, lut2_pp180
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: lut3_sigma
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: lut4_sigma

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: znum ! aerosol number per mode per unit area for each layer [m-2]

  !@@@ Currently lut_pp180 are always allocated and read - need for switch?
  !!$OMP THREADPRIVATE (lambda, sigma, omega, asym, nr, ni, lut1_sigma,lut1_g, lut1_omega, lut1_pp180,lut2_sigma, lut2_g, lut2_omega, lut2_pp180, lut3_sigma, lut4_sigma, znum)

CONTAINS

!----------------------------------------------------------------------------------------------------------------

  SUBROUTINE ham_rad_refrac(kproma, kbdim, klev, krow, ktrac, kmod, kwv, &
                             pxtm1, pnr,  pni )

    ! *ham_rad_refrac* calculates average refractive indices
    !                   for the internally mixed aerosol  
    !                   modes of ECHAM5-HAM   
    ! 
    ! Author:
    ! -------
    ! Philip Stier, MPI-Met, Hamburg,  09/04/2003
    !
    ! Modifications:
    ! --------------
    ! Philip Stier, Caltech, Pasadena, 11/2005 - modularisation
    !                                          - effective medium approximations
    !
    ! Declan O'Donnell, MPI-Met, Hamburg, 2007 - openMP bugfix (removed need for
    !                                            allocatable arrays in subordinate routines)
    ! Method:
    ! -------
    ! The real and imaginary parts of the refractive index 
    ! are obtained by taking a volume-weighted average over
    ! all compounds in the respecive mode.
    !
    ! Interface:
    ! ----------
    ! *ham_rad_refrac* is called by *ham_rad*

    USE mo_kind,         ONLY: dp
    USE mo_ham,          ONLY: nradmix

    IMPLICIT NONE

    !---subroutine interface:
    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow      ! grid parameters
    INTEGER,  INTENT(in)  :: ktrac                          ! number of tracers
    INTEGER,  INTENT(in)  :: kmod                           ! current m7 mode
    INTEGER,  INTENT(in)  :: kwv                            ! current wavelength
    REAL(dp), INTENT(in)  :: pxtm1(kbdim,klev,ktrac)        ! tracer concentrations
    REAL(dp), INTENT(out) :: pnr(kbdim,klev), pni(kbdim,klev)


    !--- 1) Calculate effective refractive index with mixing rules:

    SELECT CASE (nradmix(kmod)) 

    CASE (1)

       CALL  ham_rad_refrac_volume(kproma, kbdim, klev, krow, ktrac, kmod, kwv, &
                                    pxtm1, pnr,  pni )

    CASE (2)
       CALL  ham_rad_refrac_maxgar(kproma, kbdim, klev, krow, ktrac, kmod, kwv, &
                                    pxtm1, pnr,  pni )

    CASE (3)

       CALL  ham_rad_refrac_brugge(kproma, kbdim, klev, krow, ktrac, kmod, kwv, &
                                    pxtm1, pnr,  pni )
       
    END SELECT
    !<<dod

  END SUBROUTINE ham_rad_refrac

!----------------------------------------------------------------------------------------------------------------
  !>>dod omp bugfix
  SUBROUTINE ham_rad_refrac_volume(kproma, kbdim, klev, krow, ktrac, kmod, kwv, &
                                    pxtm1, pnr,  pni )
  !<<dod  
    ! *ham_rad_refrac_volume* calculates volume averaged 
    !                          refractive indices for the
    !                          internally mixed aerosol  
    !                          modes of ECHAM5-HAM   
    ! 
    ! Author:
    ! -------
    ! Philip Stier, Caltech, Pasadena, 11/2005 - modularisation
    !
    ! Modified
    ! --------
    ! Declan O'Donnell MPI-M 2008, SOA and OpenMP bugfix...major changes
    !
    ! Method:
    ! -------
    ! The real and imaginary parts of the refractive index 
    ! are obtained by taking a volume-weighted average over
    ! all compounds in the respecive mode.
    !
    ! Interface:
    ! ----------
    ! *ham_rad_refrac_volume* is called by *ham_rad_refrac*

    !---inherited types, data and functions
    USE mo_kind,         ONLY: dp

    IMPLICIT NONE

    !---subroutine interface:
    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow      ! grid parameters
    INTEGER,  INTENT(in)  :: ktrac                          ! number of tracers
    INTEGER,  INTENT(in)  :: kmod                           ! current m7 mode
    INTEGER,  INTENT(in)  :: kwv                            ! current wavelength
    REAL(dp), INTENT(in)  :: pxtm1(kbdim,klev,ktrac)        ! tracer concentrations
    REAL(dp), INTENT(out) :: pnr(kbdim,klev), pni(kbdim,klev)

    !--- Local variables:

    INTEGER :: jl, jk, jt, jn
    REAL(dp)    :: zdensity, zeps, zv
    INTEGER :: ikey
    !>>dod openmp bugfix
    REAL(dp) :: zvsum(kbdim,klev),  &
                znrsum(kbdim,klev), &
                znisum(kbdim,klev)
    !<<dod
    
    !---executable procedure

    zeps = THRESHOLD ! Used for both physical cutoff and safe division

    !>>dod openmp bugfix removed allocation of arrays


    !---sum over aerosol compounds:

    zvsum(1:kproma,:) =0._dp 
    znrsum(1:kproma,:)=0._dp
    znisum(1:kproma,:)=0._dp
       
    !>>dod soa
    DO jn = 1,naerocomp
       IF (aerocomp(jn)%iclass == kmod .AND. nrad(kmod) > 0) THEN
          jt = aerocomp(jn)%idt
          zdensity=aerocomp(jn)%species%density
          ikey = aerocomp(jn)%species%iaerorad

          DO jk=1,klev
             DO jl=1,kproma
                IF(pxtm1(jl,jk,jt)>zeps) THEN

                   zv=pxtm1(jl,jk,jt) / zdensity
                   
                   znrsum(jl,jk)=znrsum(jl,jk)+cnr(kwv,ikey)*zv
                   znisum(jl,jk)=znisum(jl,jk)+cni(kwv,ikey)*zv

                   zvsum(jl,jk) =zvsum(jl,jk)+zv

                END IF
             END DO
          END DO

       END IF

    END DO

    ! Add aerosol water
    IF (sizeclass(kmod)%lsoluble .AND. nrad(kmod) > 0) THEN
       jt = aerowater(kmod)%idt
       zdensity = aerowater(kmod)%species%density
       ikey = aerowater(kmod)%species%iaerorad
             
       DO jk=1, klev
          DO jl=1, kproma
             IF(pxtm1(jl,jk,jt)>zeps) THEN
                      
                zv=pxtm1(jl,jk,jt)/zdensity

                znrsum(jl,jk)=znrsum(jl,jk)+cnr(kwv,ikey)*zv
                znisum(jl,jk)=znisum(jl,jk)+cni(kwv,ikey)*zv

                zvsum(jl,jk) =zvsum(jl,jk)+zv

             END IF
          END DO
       END DO

    END IF


    !---Weighted averaging:
    DO jk=1,klev
       DO jl=1,kproma
          IF(zvsum(jl,jk)>zeps) THEN ! safe division

             pnr(jl,jk)=znrsum(jl,jk)/zvsum(jl,jk)
             pni(jl,jk)=znisum(jl,jk)/zvsum(jl,jk)

          ELSE

             pnr(jl,jk)=0._dp
             pni(jl,jk)=0._dp

          END IF
       END DO
    END DO

    !>>dod openmp bugfix removed deallocation of arrays
    !<<dod

  END SUBROUTINE ham_rad_refrac_volume

 !----------------------------------------------------------------------------------------------------------------

  !>>dod omp bugfix
  SUBROUTINE ham_rad_refrac_maxgar(kproma, kbdim, klev, krow, ktrac, kmod, kwv, &
                                    pxtm1, pnr,  pni )
  !<<dod  

    ! *ham_rad_refrac_maxgar* calculates effective refractive 
    !                          indices for the internally mixed   
    !                          aerosol modes of ECHAM5-HAM using
    !                          using the Maxwell-Garnett effective
    !                          medium approach.
    ! 
    ! Author:
    ! -------
    ! Philip Stier, Caltech, Pasadena, 11/2005
    !
    ! Modified
    ! --------
    ! Declan O'Donnell MPI-M 2008, SOA and OpenMP bugfix...major changes
    !
    ! Method:
    ! -------
    ! The component refractive indices are converted in the
    ! respective dielectric constants for which the 
    ! Maxwell-Garnett mixing rule (Garnett, 1904, 1906) is 
    ! applied. The Maxwell-Garnett mixing rule requires the
    ! choice of a host medium in which the other components are 
    ! embedded. It is assumed that 
    !
    ! Interface:
    ! ----------
    ! *ham_rad_refrac_maxgar* is called by *ham_rad_refrac*

    USE mo_kind,         ONLY: dp

    IMPLICIT NONE

    !---subroutine interface:
    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow      ! grid parameters
    INTEGER,  INTENT(in)  :: ktrac                          ! number of tracers
    INTEGER,  INTENT(in)  :: kmod                           ! current m7 mode
    INTEGER,  INTENT(in)  :: kwv                            ! current wavelength
    REAL(dp), INTENT(in)  :: pxtm1(kbdim,klev,ktrac)        ! tracer concentrations
    REAL(dp), INTENT(out) :: pnr(kbdim,klev), pni(kbdim,klev)

    !--- Local variables:

    INTEGER  :: jl, jk, jt, jn
    REAL(dp) :: zdensity, zeps, zv, zvfrac, zepsrel

    INTEGER :: ikey

    !>>dod openmp bugfix removed allocatable property of arrays
    LOGICAL  :: lcore(kbdim,klev)

    REAL(dp) :: zvsum(kbdim,klev),  zvcore(kbdim,klev), &
                znrsum(kbdim,klev), znisum(kbdim,klev)

    COMPLEX  :: ce             ! component dielectric constant

    COMPLEX  :: cn_eff(kbdim,klev), & ! mode effective refractive index
                ce_eff(kbdim,klev), & ! mode effective dielectric constant
                cn_0(kbdim,klev),   & ! host medium effective refractive index
                ce_0(kbdim,klev),   & ! host medium effective dielectric constant
                csum(kbdim,klev)      ! local summation term
    !<<dod
 
    !---executable procedure

    zeps = THRESHOLD   ! Used for both physical cutoff and safe division
    zepsrel = EPSILON(1._dp)*100. ! Used for bound checking

    !>>dod openmp bugfix removed allocation of arrays
    !<<dod

    zvsum(1:kproma,:)  = 0.0_dp
    zvcore(1:kproma,:) = 0.0_dp
    znrsum(1:kproma,:) = 0.0_dp
    znisum(1:kproma,:) = 0.0_dp

    !>>dod soa
    DO jn = 1,naerocomp
       IF (aerocomp(jn)%iclass == kmod .AND. nrad(kmod) > 0) THEN
          jt = aerocomp(jn)%idt
          zdensity=aerocomp(jn)%species%density
          ikey = aerocomp(jn)%species%iaerorad

          !---volume of insoluble core:

          IF (.NOT. aerocomp(jn)%species%lwatsol) THEN
             zvcore(1:kproma,:) = zvcore(1:kproma,:)+pxtm1(1:kproma,:,jt)/zdensity
          END IF
       
          !--- Total mode volume and summation of refractive indices:

          DO jk=1,klev
             DO jl=1,kproma
                IF(pxtm1(jl,jk,jt)>zeps) THEN

                   zv=pxtm1(jl,jk,jt)/zdensity

                   znrsum(jl,jk)=znrsum(jl,jk)+cnr(kwv,ikey)*zv
                   znisum(jl,jk)=znisum(jl,jk)+cni(kwv,ikey)*zv

                   zvsum(jl,jk) =zvsum(jl,jk)+zv

                END IF
             END DO
          END DO

       END IF
    END DO

    ! Add aerosol water
    IF (sizeclass(kmod)%lsoluble .AND. nrad(kmod) > 0) THEN
       jt = aerowater(kmod)%idt
       zdensity = aerowater(kmod)%species%density
       ikey = aerowater(kmod)%species%iaerorad
             
       DO jk=1, klev
          DO jl=1, kproma
             IF(pxtm1(jl,jk,jt)>zeps) THEN
                      
                zv=pxtm1(jl,jk,jt)/zdensity

                znrsum(jl,jk)=znrsum(jl,jk)+cnr(kwv,ikey)*zv
                znisum(jl,jk)=znisum(jl,jk)+cni(kwv,ikey)*zv

                zvsum(jl,jk) =zvsum(jl,jk)+zv

             END IF
          END DO
       END DO

    END IF
    !<<dod soa

    !--- Volume weighted averaging for all regions:

    !>>dod openmp bugfix
    DO jk=1,klev
       DO jl=1,kproma
    !<<dod
          IF(zvsum(jl,jk)>zeps) THEN ! Safe division
             cn_eff(jl,jk)=CMPLX( znrsum(jl,jk)/zvsum(jl,jk) , znisum(jl,jk)/zvsum(jl,jk), kind=dp )
          ELSE
             cn_eff(jl,jk)=CMPLX( 0.0_dp , 0.0_dp, kind=dp )
          END IF
          
       END DO
    END DO

    !--- 2) Apply Maxwell-Garnett for regions with insoluble core:

    !--- Find regions with insoluble core:

    lcore(1:kproma,:)=.FALSE.

    !>>dod openmp bugfix
    DO jk=1,klev
       DO jl=1,kproma
    !<<dod
          IF (zvsum(jl,jk)>0._dp) THEN
             IF (zvcore(jl,jk) > zvsum(jl,jk)*zepsrel .AND. zvcore(jl,jk) < zvsum(jl,jk)*(1.0_dp-zepsrel)) THEN
                lcore(jl,jk)=.TRUE.
             END IF
          END IF
       END DO
    END DO

    !--- 2.1) Calculate volume weighted refractive index of host medium:  

    zvsum(1:kproma,:)  = 0.0_dp
    znrsum(1:kproma,:) = 0.0_dp
    znisum(1:kproma,:) = 0.0_dp

    DO jn = 1,naerocomp
       IF (aerocomp(jn)%iclass == kmod .AND. nrad(kmod) > 0) THEN
          IF (aerocomp(jn)%species%lwatsol) THEN
             jt = aerocomp(jn)%idt
             zdensity=aerocomp(jn)%species%density
             ikey = aerocomp(jn)%species%iaerorad

             DO jk=1,klev
                DO jl=1,kproma
                   IF(pxtm1(jl,jk,jt)>zeps .AND. lcore(jl,jk)) THEN

                      zv=pxtm1(jl,jk,jt)/zdensity

                      znrsum(jl,jk)=znrsum(jl,jk)+cnr(kwv,ikey)*zv
                      znisum(jl,jk)=znisum(jl,jk)+cni(kwv,ikey)*zv

                      zvsum(jl,jk) =zvsum(jl,jk)+zv

                   END IF
                END DO
             END DO

          END IF
       END IF

    END DO

    ! Add aerosol water
    IF (sizeclass(kmod)%lsoluble .AND. nrad(kmod) > 0) THEN
       jt = aerowater(kmod)%idt
       zdensity = aerowater(kmod)%species%density
       ikey = aerowater(kmod)%species%iaerorad
             
       DO jk=1, klev
          DO jl=1, kproma
             IF(pxtm1(jl,jk,jt)>zeps) THEN !Physical cutoff
                      
                zv=pxtm1(jl,jk,jt)/zdensity

                znrsum(jl,jk)=znrsum(jl,jk)+cnr(kwv,ikey)*zv
                znisum(jl,jk)=znisum(jl,jk)+cni(kwv,ikey)*zv

                zvsum(jl,jk) =zvsum(jl,jk)+zv

             END IF
          END DO
       END DO

    END IF
    !<<dod soa

    !--- Volume weighted averaging of host medium refractive index:

    DO jk=1,klev
       DO jl=1,kproma
          IF(zvsum(jl,jk)>zeps .AND. lcore(jl,jk)) THEN ! Safe division
             cn_0(jl,jk)=CMPLX( znrsum(jl,jk)/zvsum(jl,jk) , znisum(jl,jk)/zvsum(jl,jk), kind=dp )
          ELSE
             cn_0(jl,jk)=CMPLX( 0.0_dp , 0.0_dp, kind=dp )
          END IF

          !--- Dielectric constant:

          ce_0(jl,jk)=cn_0(jl,jk)**2

       END DO
    END DO

    !--- 2.2) Apply Maxwell-Garnett mixing rule for insoluble core:
    
    !--- Calculate total mode volume

    zvsum(1:kproma,:) = 0.0_dp

    DO jn = 1,naerocomp
       IF (aerocomp(jn)%iclass == kmod .AND. nrad(kmod) > 0) THEN
          IF (aerocomp(jn)%species%lwatsol) THEN
             jt = aerocomp(jn)%idt
             zdensity=aerocomp(jn)%species%density
             ikey = aerocomp(jn)%species%iaerorad

             DO jk=1,klev
                DO jl=1,kproma
                   IF(pxtm1(jl,jk,jt)>zeps .AND. lcore(jl,jk)) THEN
                      zvsum(jl,jk) =zvsum(jl,jk)+pxtm1(jl,jk,jt)/zdensity
                   END IF
                END DO
             END DO

          END IF
       END IF
    END DO

    ! Add aerosol water
    IF (sizeclass(kmod)%lsoluble .AND. nrad(kmod) > 0) THEN
       jt = aerowater(kmod)%idt
       zdensity = aerowater(kmod)%species%density
       ikey = aerowater(kmod)%species%iaerorad
             
       DO jk=1, klev
          DO jl=1, kproma
             IF(pxtm1(jl,jk,jt)>zeps .AND. lcore(jl,jk)) THEN
                zvsum(jl,jk) =zvsum(jl,jk)+pxtm1(jl,jk,jt)/zdensity
             END IF
          END DO
       END DO
    END IF

    !--- Apply M&G for the insoluble core components embedded in soluble host medium:

    csum(1:kproma,:)=CMPLX(0.0_dp,0.0_dp, kind=dp)

    !>>dod soa
    DO jn = 1,naerocomp
       IF (aerocomp(jn)%iclass == kmod .AND. nrad(kmod) > 0) THEN
          
          IF (.NOT. aerocomp(jn)%species%lwatsol) THEN
             jt = aerocomp(jn)%idt

             zdensity=aerocomp(jn)%species%density
             ikey = aerocomp(jn)%species%iaerorad


             ce=CMPLX(cnr(kwv,ikey),cni(kwv,ikey), kind=dp)**2

             DO jk=1,klev
                DO jl=1,kproma
                   IF(pxtm1(jl,jk,jt)>zeps .AND. lcore(jl,jk)) THEN !Physical cutoff/safe division implied

                      zvfrac=(pxtm1(jl,jk,jt)/zdensity) / zvsum(jl,jk)

                      csum(jl,jk)=csum(jl,jk)+zvfrac*(ce-ce_0(jl,jk))/(ce+2.0_dp*ce_0(jl,jk))

                   END IF
                END DO
             END DO

          END IF
       END IF

    END DO
    !<<dod soa

    !--- Calculation of effective dielectric constant:

    !>>dod openmp bugfix
    DO jk=1,klev
       DO jl=1,kproma
    !<<dod      
          IF(lcore(jl,jk)) THEN
             ce_eff(jl,jk)=ce_0(jl,jk)*(1.0_dp+2.0_dp*csum(jl,jk))/(1.0_dp-csum(jl,jk))
          END IF
       END DO
    END DO

    !--- Calculation of effective refractive index:

    DO jk=1,klev
       DO jl=1,kproma

          !--- Replace volume averaged first guess for regions with insoluble core:

          IF (lcore(jl,jk)) THEN 
             cn_eff(jl,jk)=SQRT(ce_eff(jl,jk))
          END IF

          !--- Store in real variables: 

          pnr(jl,jk)=REAL(cn_eff(jl,jk))
          pni(jl,jk)=AIMAG(cn_eff(jl,jk))

       END DO
    END DO
    
  END SUBROUTINE ham_rad_refrac_maxgar

 !----------------------------------------------------------------------------------------------------------------

  !>>dod omp bugfix
  SUBROUTINE ham_rad_refrac_brugge(kproma, kbdim, klev, krow, ktrac, kmod, kwv, &
                                    pxtm1, pnr,  pni )
  !<<dod  

    ! *ham_rad_refrac_brugge* calculates effective refractive 
    !                          indices for the internally mixed   
    !                          aerosol modes of ECHAM5-HAM using
    !                          using the Bruggeman effective
    !                          medium approach.
    ! 
    ! Author:
    ! -------
    ! Philip Stier, Caltech, Pasadena, 11/2005
    !
    ! Modified
    ! --------
    ! Declan O'Donnell MPI-M 2008, SOA and OpenMP bugfix...major changes
    !
    ! Method:
    ! -------
    ! The component refractive indices are converted in the
    ! respective dielectric constants for which the 
    ! Bruggeman mixing rule (Bruggeman, 1935) is 
    ! applied. The Bruggeman multicomponent mixing rule is an 
    ! implicit equation for the complex dielectric constants    
    ! that is solved with a Newton-Raphson iteration procedure.
    !
    ! Interface:
    ! ----------
    ! *ham_rad_refrac_brugge* is called by *ham_rad_refrac*

    USE mo_kind,         ONLY: dp

    IMPLICIT NONE

    !---subroutine interface:
    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow      ! grid parameters
    INTEGER,  INTENT(in)  :: ktrac                          ! number of tracers
    INTEGER,  INTENT(in)  :: kmod                           ! current m7 mode
    INTEGER,  INTENT(in)  :: kwv                            ! current wavelength
    REAL(dp), INTENT(in)  :: pxtm1(kbdim,klev,ktrac)        ! tracer concentrations
    REAL(dp), INTENT(out) :: pnr(kbdim,klev), pni(kbdim,klev)

    !--- Local variables:
    INTEGER  :: jl, jk, jt, jiter, jn
    REAL(dp) :: zdensity, zeps, zvfrac
    !>>dod soa
    INTEGER :: ikey
    !<<dod

    !>>dod openmp bugfix removed allocatable property of arrays

    REAL(dp) :: zvsum(kbdim,klev)

    COMPLEX :: ce, cn_eff_old ! component dielectric constant

    COMPLEX :: cn_eff(kbdim,klev), & ! mode effective refractive index
               ce_eff(kbdim,klev), & ! mode effective dielectric constant
               cfe(kbdim,klev),    & ! f(e)
               cfep(kbdim,klev)      ! f'(e)

    INTEGER, PARAMETER   :: niter=7        ! maximum number of iterations
                                           ! (1 year test run showed convergence
                                           ! after 6 iterations everywhere      )


    !---executable procedure

    zeps=THRESHOLD   ! Used for both physical cutoff and safe division 

    !>>dod deleted allocation of arrays
    !<<dod

    zvsum(1:kproma,:) =0._dp
    cn_eff(1:kproma,:)=CMPLX(0.0_dp,0.0_dp, kind=dp)

    !--- 1) Calculate total volume of the mode:

    !>>dod soa
    DO jn = 1,naerocomp
       IF (aerocomp(jn)%iclass == kmod .AND. nrad(kmod) > 0) THEN
          jt = aerocomp(jn)%idt
          zdensity=aerocomp(jn)%species%density
          ikey = aerocomp(jn)%species%iaerorad

          zvsum(1:kproma,:) =zvsum(1:kproma,:)+pxtm1(1:kproma,:,jt)/zdensity

       END IF

    END DO

    ! Add aerosol water
    IF (sizeclass(kmod)%lsoluble .AND. nrad(kmod) > 0) THEN
       jt = aerowater(kmod)%idt
       zdensity = aerowater(kmod)%species%density
       ikey = aerowater(kmod)%species%iaerorad
             
       zvsum(1:kproma,:) =zvsum(1:kproma,:)+pxtm1(1:kproma,:,jt)/zdensity

    END IF


    !--- 2) Calculate f(e0) and f'(e0)
    !       and start Newtonian iteration:

    !--- First guess for effective dielectric constant:


    !--- Newtonian iteration:

    !--- Choose POM values as initial guess (approximate zero) for Newton-Raphson method:
    !    (Radiative properties of POM lie within range of other compounds)

    ce_eff(1:kproma,:)=CMPLX(cnr(kwv,speclist(id_oc)%iaerorad),cni(kwv,speclist(id_oc)%iaerorad), kind=dp)**2

    DO jiter=1, niter

       cfe(1:kproma,:)  = CMPLX(0.0_dp,0.0_dp, kind=dp)
       cfep(1:kproma,:) = CMPLX(0.0_dp,0.0_dp, kind=dp)

       !>>dod soa
       DO jn = 1,naerocomp
          IF (aerocomp(jn)%iclass == kmod .AND. nrad(kmod) > 0) THEN
             jt = aerocomp(jn)%idt
             zdensity=aerocomp(jn)%species%density
             ikey = aerocomp(jn)%species%iaerorad

             !--- Component dielectric constant:
             
             ce=CMPLX(cnr(kwv,ikey),cni(kwv,ikey), kind=dp)**2

             DO jk=1,klev
                DO jl=1,kproma
                   IF(pxtm1(jl,jk,jt)>zeps) THEN

                      !--- Component volume fraction:

                      zvfrac = (pxtm1(jl,jk,jt)/zdensity) / zvsum(jl,jk)

                      !--- Calculate f(e):

                      cfe(jl,jk)=cfe(jl,jk) + &
                                 zvfrac*(ce-ce_eff(jl,jk))/(ce+2.0_dp*ce_eff(jl,jk))

                      !--- Calculate f'(e):

                      cfep(jl,jk)=cfep(jl,jk) + &
                                  zvfrac*( (-3.0_dp*ce) / ( (ce+2.0_dp*ce_eff(jl,jk))**2 ) )

                   END IF
                END DO
             END DO

             !--- Solve for new approximation of effective dielectric constant:

             DO jk=1,klev
                DO jl=1,kproma
                   IF(zvsum(jl,jk)>zeps) THEN

                      !                  IF (CABS(cfep(jl,jk))<zeps ) cfep(jl,jk)=CMPLX(zeps,zeps, kind=dp)

                      !                   cn_eff_old=SQRT(ce_eff(jl,jk))
                      
                      ce_eff(jl,jk)=SQRT(ce_eff(jl,jk)-cfe(jl,jk)/cfep(jl,jk))

                   END IF
                END DO
             END DO

          END IF

       END DO     ! naerocomp

    END DO        ! niter

    !--- Store refractive index in real variables:

    DO jk=1,klev
       DO jl=1,kproma

          pnr(jl,jk)=REAL(cn_eff(jl,jk))
          pni(jl,jk)=AIMAG(cn_eff(jl,jk))

       END DO
    END DO

       
  END SUBROUTINE ham_rad_refrac_brugge

!----------------------------------------------------------------------------------------------------------------

  SUBROUTINE ham_rad(kproma, kbdim, klev, krow, kpband, kb_sw,                 &
       pxtm1,         ppd_hl,                                    &
       aer_tau_sw_vr, aer_piz_sw_vr, aer_cg_sw_vr, aer_tau_lw_vr, rwet_m7, &
       & ldiag_aeropt, kb_diag, ntype_diaf, &
       & lambda_diag, zaer_tau_diag, zaer_ssa_diag, zaer_asym_diag, zaer_lidarR_diag)
    ! *ham_rad* calculates optical properties for
    !            aerosol distributions from look-up
    !            tables.
    !
    ! Author:
    ! -------
    ! Oliver Boucher, Univ. Lille,        2003
    !
    ! Modifications:
    ! --------------
    ! Philip Stier, MPI-MET, Hamburg,     2003
    !
    ! Declan O'Donnell MPI-MET, Hamburg, 2008
    ! restuctured to loop around the wavelengths and removed allocatable 
    ! arrays to enable running under OpenMP. 
    !
    ! Method:
    ! -------
    ! To be done!
    !
    ! Interface:
    ! ----------
    ! *ham_rad* is called by *radiation*

    USE mo_ham,           ONLY: naerorad, nrad, nham_subm, HAM_SALSA, HAM_M7, &
                                sigma_fine, sigma_coarse !SF #320
    USE mo_ham_m7ctl,     ONLY: modesigma=>sigma
    USE mo_math_constants, ONLY: pi
    USE mo_physical_constants, ONLY: grav
    USE mo_kind,          ONLY: dp
    USE mo_exception,     ONLY: finish
    USE mo_tracdef,       ONLY: ntrac
#ifdef HAMMOZ
    USE mo_ham_streams,   ONLY: rwet
#endif
#ifdef SALSA
    USE mo_ham_salsa,     ONLY: rwet_salsa
    USE mo_ham_salsactl,  ONLY: fn2a, fn2b, nbin3
#endif
    USE mo_control,       ONLY: ltimer
    !>>dod split of mo_timer (#51)
#ifdef HAMMOZ
    USE mo_hammoz_timer,  ONLY: timer_start, timer_stop, &
         timer_ham_rad_fitplus,  &
         timer_ham_rad_refrac
#endif
    USE mo_ham_rad_data, ONLY: nraddiagwv
#ifdef HAMMOZ
    USE mo_ham_streams,  ONLY: tau_mode
#endif


    IMPLICIT NONE

    !--- Arguments:

    INTEGER,INTENT(IN)      :: kproma , kbdim, klev, krow, kpband, kb_sw,kb_diag,ntype_diaf

    REAL(dp), INTENT(in)    :: ppd_hl(kbdim,klev)                 ! pressure diff between half levels [Pa]

    REAL(dp),INTENT(IN)     :: pxtm1(kbdim,klev,ntrac)            ! tracer mass/number mixing ratio (t-dt)   [kg/kg]/[#/kg]

    REAL(dp), INTENT(inout) :: aer_tau_lw_vr(kbdim,klev,kpband),& !< LW optical thickness of aerosols
                               aer_tau_sw_vr(kbdim,klev,kb_sw), & !< aerosol optical thickness
                               aer_cg_sw_vr(kbdim,klev,kb_sw),  & !< aerosol asymmetry factor
                               aer_piz_sw_vr(kbdim,klev,kb_sw)    !< aerosol single scattering albedo

    ! diagnostic aerosol optical properties
    logical,intent(in)            :: ldiag_aeropt ! logical for aerosol optics
    real(dp),intent(in)    :: lambda_diag(kb_diag)
    real(dp),intent(inout) :: zaer_tau_diag(kbdim,klev,kb_diag)
    real(dp),intent(inout) :: zaer_ssa_diag(kbdim,klev,kb_diag)
    real(dp),intent(inout) :: zaer_asym_diag(kbdim,klev,kb_diag)
    real(dp),intent(inout) :: zaer_lidarR_diag(kbdim,klev,kb_diag)
    real(dp) :: zaer_tau180_diag(kbdim,klev,kb_diag)

    REAL(dp) :: sigma_diag(kbdim,klev,kb_diag,nclass),    &
                omega_diag(kbdim,klev,kb_diag,nclass), &
                asym_diag (kbdim,klev,kb_diag,nclass), &
                pp180_diag(kbdim,klev,kb_diag,nclass), &
                nr_diag(kbdim,klev,kb_diag,nclass),       &
                ni_diag(kbdim,klev,kb_diag,nclass)
    !--- Local Variables:
#ifdef HAMMOZ       

#else
    REAL(dp), INTENT(in)    ::    rwet_m7(kbdim,klev,nclass) 
#endif  
    INTEGER  :: jclass, jl, jk, jwv, jwv_diag, itable, itrac, ikl

    REAL(dp) :: zeps

    REAL(dp) :: zxx(kbdim,klev),                           & ! size parameter
                zdpg(kbdim,klev)                             ! auxiliary parameter dp/grav

    REAL(dp) :: zaer_tau_sw_vr(kbdim,klev,Nwv_sw_tot,nclass),& ! SW optical depth for each band and mode
                zaer_tau_lw_vr(kbdim,klev,Nwv_lw,nclass),&      ! LW optical depth for each band and mode
                zaer_tau_diag_vr(kbdim,klev,kb_diag,nclass)    ! diag optical depth for each band and mode

!>>gf: needed to avoid architecture-dependent problems (Cray XT5)
    REAL(dp) :: znr2d(kbdim,klev),                         & ! 2D subset of 4D array nr
                zni2d(kbdim,klev),                         & ! 2D subset of 4D array ni
                zsigma2d(kbdim,klev),                      & ! 2D subste of 4D array sigma
                zomega2d(kbdim,klev),                      & ! 2D subste of 4D array omega
                zasym2d(kbdim,klev),                       & ! 2D subste of 4D array asym
                zpp180deg2d(kbdim,klev)                      ! Phase function at 180 degrees
!<<gf


    REAL(dp) :: sigma(kbdim,klev,Nwv_tot,nclass),    &
                omega(kbdim,klev,Nwv_sw_tot,nclass), &
                asym (kbdim,klev,Nwv_sw_tot,nclass), &
                nr(kbdim,klev,Nwv_tot,nclass),       &
                ni(kbdim,klev,Nwv_tot,nclass),       &
                znum(kbdim,klev,nclass)

    INTEGER :: jlwv                                          ! RRTM-LW band number
    
    !--- Stream:
#ifdef HAMMOZ
    REAL(dp), POINTER     :: rwet_p(:,:,:)
#endif
    !---executable procedure

    !--- 0) Initialization:

    zeps=EPSILON(1._dp)*10._dp   ! Used for comparisons with sigma_coarse and sigma_fine 

    sigma(1:kproma,:,:,:)=0._dp
    omega(1:kproma,:,:,:)=0._dp
    asym(1:kproma,:,:,:) =0._dp

    sigma_diag(1:kproma,:,:,:)=0._dp
    omega_diag(1:kproma,:,:,:)=0._dp
    asym_diag (1:kproma,:,:,:)=0._dp
    !nr(1:kproma,:,:,:) =0._dp
    !ni(1:kproma,:,:,:) =0._dp
    !znum(1:kproma,:,:) =0._dp

    zdpg(1:kproma,:)=ppd_hl(1:kproma,:)/grav

    DO jclass=1, nclass
       itrac=sizeclass(jclass)%idt_no
       ! [#/kg] * [Pa] / [m/s2] -> #/m2 - Layer integrated
       znum(1:kproma,:,jclass)=pxtm1(1:kproma,:,itrac)*zdpg(1:kproma,:)
       ! IDEA: minimum density [#/m3] = 1e-15 (as used in
       ! tm5m7_optics_aop_get.F90). Applying to layer makes it a
       ! smaller limit.
       WHERE(znum(1:kproma,:,jclass) .lt. 1.e-15_dp ) znum(1:kproma,:,jclass) = 1.e-15_dp
    END DO

    !--- 1) Calculate optical properties for GCM SW bands:

    IF (ANY(nrad(:)==1) .OR. ANY(nrad(:)==3)) THEN

       DO jclass=1, nclass
          IF (nrad(jclass)==1 .OR. nrad(jclass)==3) THEN

#ifdef HAMMOZ
             rwet_p => rwet(jclass)%ptr
#endif

             DO jwv=1,Nwv_sw+Nwv_sw_opt

                !--- 1.1) Calculate volume averaged refractive index nr and ni:
#ifdef HAMMOZ
                IF (ltimer) CALL timer_start(timer_ham_rad_refrac)
#endif  
                !gf: the former usage of nr(1:kproma,:,jwv,jclass) and ni(1:kproma,:,jwv,jclass)
                !    directly in the call to ham_rad_refrac is causing architecture-dependent problems (Cray XT5)
                !    Therefore intermediate variables znr2d and zni2d are introduced
                !pls: no need to intialize them, they are intent(out) in ham_rad_refrac. Comment these 2 lines:
                ! znr2d(1:kproma,:) = nr(1:kproma,:,jwv,jclass)
                ! zni2d(1:kproma,:) = ni(1:kproma,:,jwv,jclass)

                CALL ham_rad_refrac(kproma, kbdim, klev, krow, &
                                     ntrac,  jclass,  jwv,     &
                                     pxtm1,  znr2d, zni2d )

#ifdef HAMMOZ
                IF (ltimer) CALL timer_stop(timer_ham_rad_refrac)
#endif  
                !--- 1.1) Calculate size parameter:
#ifdef HAMMOZ          
                zxx(1:kproma,:) = 2._dp*pi*rwet_p(1:kproma,:,krow)/lambda(jwv)
#else
                zxx(1:kproma,:) = 2._dp*pi*rwet_m7(1:kproma,:,jclass)/lambda(jwv) 
#endif  

                !--- 1.2) Table-lookup for optical properties:
#ifdef HAMMOZ 
                IF (ltimer) CALL timer_start(timer_ham_rad_fitplus)
#endif 
                !gf: same as in the call to ham_rad_refrac, for the call to ham_rad_fitplus

                zsigma2d(1:kproma,:) = sigma(1:kproma,:,jwv,jclass)
                zomega2d(1:kproma,:) = omega(1:kproma,:,jwv,jclass)
                zasym2d(1:kproma,:)  = asym(1:kproma,:,jwv,jclass)

                SELECT CASE(nham_subm)

                CASE(HAM_M7)

                   IF (ABS(modesigma(jclass)-sigma_fine) < zeps*sigma_fine) THEN
                      itable=1
                   ELSE IF (ABS(modesigma(jclass)-sigma_coarse) < zeps*sigma_coarse) THEN
                      itable=2
                   ELSE 
                      CALL finish('ham_rad','incompatible standard deviation in modal setup')
                   END IF

                CASE(HAM_SALSA)

#ifdef SALSA
                   zxx(1:kproma,:) = 2._dp*pi*rwet_salsa(1:kproma,:,jclass)/lambda(jwv)

                   IF(jclass < 6 .OR. (jclass > fn2a .AND. jclass < fn2b-(nbin3-1))) THEN
                      itable=1                      
                   ELSE                     
                      itable=2                      
                   END IF
#endif
                END SELECT

                IF (itable == 1) THEN
                   CALL ham_rad_fitplus(kproma, kbdim,      klev,                &
                                         zxx,   znr2d,     zni2d,                &
                                         itable, lut1_sigma, zsigma2d,           & 
                                                 lut1_omega, zomega2d,           &
                                                 lut1_g,     zasym2d             )

                ELSE 

                   CALL ham_rad_fitplus(kproma, kbdim,      klev,                & 
                                         zxx,   znr2d,     zni2d,                &
                                         itable, lut2_sigma, zsigma2d,           &
                                                 lut2_omega, zomega2d,           &
                                                 lut2_g,     zasym2d             )

                END IF
#ifdef HAMMOZ
                IF (ltimer) CALL timer_stop(timer_ham_rad_fitplus)
#endif
                !>>gf: update the original 4d arrays
                sigma(1:kproma,:,jwv,jclass) = zsigma2d(1:kproma,:)*lambda(jwv)*lambda(jwv)
                omega(1:kproma,:,jwv,jclass) = zomega2d(1:kproma,:)
                asym(1:kproma,:,jwv,jclass)  = zasym2d(1:kproma,:)
                nr(1:kproma,:,jwv,jclass)    = znr2d(1:kproma,:)
                ni(1:kproma,:,jwv,jclass)    = zni2d(1:kproma,:)
                !<<gf

          !ALLOCATE(sigma(kbdim,klev,Nwv_tot,nclass))
          !ALLOCATE(omega(kbdim,klev,Nwv_sw_tot,nclass))
          !ALLOCATE(asym (kbdim,klev,Nwv_sw_tot,nclass))
          !ALLOCATE(nr(kbdim,klev,Nwv_tot,nclass))
          !ALLOCATE(ni(kbdim,klev,Nwv_tot,nclass))
                !sigma(1:kproma,:,jwv,jclass) = sigma(1:kproma,:,jwv,jclass)*lambda(jwv)*lambda(jwv)
         
             END DO ! jwv

          END IF ! nrad
       END DO ! jclass

       !--- Summation over modes, conversion into total extinction and calculation 
       !    of weighted single scattering albedo and assymetry factor. 

       zaer_tau_sw_vr(1:kproma,:,:,:)=0.0_dp
 
       DO jclass=1, nclass
          IF(nrad(jclass)>0) THEN
             DO jwv=1, Nwv_sw+Nwv_sw_opt
                zaer_tau_sw_vr(1:kproma,:,jwv,jclass)=znum(1:kproma,:,jclass)*sigma(1:kproma,:,jwv,jclass)
             END DO
          END IF
       END DO

       !--- Diagnose AOD for requested each mode (nrad) and wavelength (nraddiagwv):
#ifdef HAMMOZ
       DO jclass=1, nclass
          IF(nrad(jclass)>0) THEN
             DO jwv=1, Nwv_sw+Nwv_sw_opt
                IF (nraddiagwv(jwv)>0) THEN
                   tau_mode(jclass,jwv)%ptr(1:kproma,:,krow)=zaer_tau_sw_vr(1:kproma,:,jwv,jclass)
                END IF
             END DO
          END IF
       END DO
#endif
       !--- Calculation of weighted properties and vertical reordering to RRTM structure:

       DO jwv=1, Nwv_sw !ham_ps +Nwv_sw_opt

          DO jclass=1, nclass
             DO jk=1, klev
#ifdef HAMMOZ       
               ikl=klev+1-jk
#else       
               !No vertical reordering here
               ikl=jk
#endif               
                DO jl=1, kproma
                   aer_tau_sw_vr(jl,jk,jwv)=aer_tau_sw_vr(jl,jk,jwv) + &
                                            zaer_tau_sw_vr(jl,ikl,jwv,jclass)
                   aer_piz_sw_vr(jl,jk,jwv)=aer_piz_sw_vr(jl,jk,jwv) + &
                                            zaer_tau_sw_vr(jl,ikl,jwv,jclass)*omega(jl,ikl,jwv,jclass)
                   aer_cg_sw_vr(jl,jk,jwv) =aer_cg_sw_vr(jl,jk,jwv) + &
                                            zaer_tau_sw_vr(jl,ikl,jwv,jclass)*omega(jl,ikl,jwv,jclass)*asym(jl,ikl,jwv,jclass)
                END DO
             END DO
          END DO

          DO jk=1, klev
             DO jl=1, kproma
                ! PLS - TODO: SAFE DIVISION
                IF (aer_piz_sw_vr(jl,jk,jwv) > THRESHOLD ) THEN 
                   aer_cg_sw_vr(jl,jk,jwv) =aer_cg_sw_vr(jl,jk,jwv)/aer_piz_sw_vr(jl,jk,jwv)
                END IF
                IF (aer_tau_sw_vr(jl,jk,jwv) > THRESHOLD ) THEN 
                   aer_piz_sw_vr(jl,jk,jwv)=aer_piz_sw_vr(jl,jk,jwv)/aer_tau_sw_vr(jl,jk,jwv)
                END IF
             END DO
          END DO

       END DO

    END IF

    !--- 2) Calculate optical properties for GCM LW bands:

    IF (ANY(nrad(:)==2) .OR. ANY(nrad(:)==3)) THEN

       DO jclass=1, nclass
          IF (nrad(jclass)==2 .OR. nrad(jclass)==3) THEN

#ifdef HAMMOZ
             rwet_p => rwet(jclass)%ptr
#endif  

             DO jwv=1, Nwv_lw

                jlwv = Nwv_sw+Nwv_sw_opt+jwv ! Total SW wavelengths + LW

                !--- 1.1) Calculate volume averaged refractive index nr and ni:
#ifdef HAMMOZ
                IF (ltimer) CALL timer_start(timer_ham_rad_refrac)
#endif
                !gf: the former usage of nr(1:kproma,:,jwv,jclass) and ni(1:kproma,:,jwv,jclass)
                !    directly in the call to ham_rad_refrac is causing architecture-dependent problems (Cray XT5)
                !    Therefore intermediate variables znr2d and zni2d are introduced
                !pls: no need to intialize them, they are intent(out) in ham_rad_refrac. Comment these 2 lines:
                ! znr2d(1:kproma,:) = nr(1:kproma,:,jlwv,jclass)
                ! zni2d(1:kproma,:) = ni(1:kproma,:,jlwv,jclass)

                CALL ham_rad_refrac(kproma, kbdim, klev, krow,                                &
                     ntrac, jclass,  jlwv,                                       &
                     pxtm1, znr2d, zni2d )
#ifdef HAMMOZ
                IF (ltimer) CALL timer_stop(timer_ham_rad_refrac)
#endif
                !--- 1.1) Calculate size parameter:
#ifdef HAMMOZ         
                zxx(1:kproma,:)=2._dp*pi*rwet_p(1:kproma,:,krow)/lambda(jlwv)
#else
                zxx(1:kproma,:) = 2._dp*pi*rwet_m7(1:kproma,:,jclass)/lambda(jlwv) 
#endif  

                !--- 1.2) Table-lookup for optical properties:
#ifdef HAMMOZ
                IF (ltimer) CALL timer_start(timer_ham_rad_fitplus)
#endif
                !gf: same as in the call to ham_rad_refrac, for the call to ham_rad_fitplus

                zsigma2d(1:kproma,:) = sigma(1:kproma,:,jlwv,jclass)

                SELECT CASE(nham_subm)

                CASE(HAM_M7)


                   IF (ABS(modesigma(jclass)-sigma_fine) < zeps*sigma_fine) THEN
                      itable=3
                   ELSE IF (ABS(modesigma(jclass)-sigma_coarse) < zeps*sigma_coarse) THEN
                      itable=4
                   ELSE 
                      CALL finish('ham_rad','incompatible standard deviation in modal setup')
                   END IF

                CASE(HAM_SALSA)
#ifdef SALSA
                   zxx(1:kproma,:)=2._dp*pi*rwet_salsa(1:kproma,:,krow)/lambda(jlwv)

                   IF(jclass < 6 .OR. (jclass > fn2a .AND. jclass < fn2b-(nbin3-1))) THEN
                      itable=3
                   ELSE                     
                      itable=4                      
                   END IF
#endif
                END SELECT

                IF (itable == 3) THEN

                   CALL ham_rad_fitplus(kproma,  kbdim,      klev,        & 
                                         zxx,    znr2d,     zni2d,        &
                                         itable, lut3_sigma, zsigma2d     )

                ELSE

                   CALL ham_rad_fitplus(kproma, kbdim,      klev,         & 
                                         zxx,    znr2d,    zni2d,         &
                                         itable, lut4_sigma, zsigma2d     )

                END IF
#ifdef HAMMOZ
                IF (ltimer) CALL timer_stop(timer_ham_rad_fitplus)
#endif
                !>>gf: update the original 4d arrays
                sigma(1:kproma,:,jlwv,jclass) = zsigma2d(1:kproma,:)
                nr(1:kproma,:,jlwv,jclass)    = znr2d(1:kproma,:)
                ni(1:kproma,:,jlwv,jclass)    = zni2d(1:kproma,:)
                !<<gf

                sigma(1:kproma,:,jlwv,jclass)=sigma(1:kproma,:,jlwv,jclass)*lambda(jlwv)*lambda(jlwv)

             END DO !jwv

          END IF !nrad
       END DO !nclass
          
       DO jclass=1, nclass
          IF(nrad(jclass)>0) THEN
             DO jwv=1, Nwv_lw
                jlwv = Nwv_sw+Nwv_sw_opt+jwv ! Total SW wavelengths + LW

                DO jk=1, klev
#ifdef HAMMOZ       
                  ikl=klev+1-jk
#else       
                  !No vertical reordering in openifs (here)
                  ikl=jk
#endif       
                   DO jl=1, kproma
                      zaer_tau_lw_vr(jl,jk,jwv,jclass)=znum(jl,jk,jclass)*sigma(jl,jk,jlwv,jclass)
                      aer_tau_lw_vr(jl,ikl,jwv)=aer_tau_lw_vr(jl,ikl,jwv) + &
                           zaer_tau_lw_vr(jl,jk,jwv,jclass)
                   END DO
                END DO

                !--- Diagnose AOD for requested each mode (nrad) and wavelength (nraddiagwv):
#ifdef HAMMOZ
                IF (nraddiagwv(jwv)>0) THEN
                   tau_mode(jclass,jwv)%ptr(1:kproma,:,krow)=zaer_tau_lw_vr(1:kproma,:,jwv,jclass)
                END IF
#endif
             END DO
          END IF
       END DO

    END IF





 !--- 3) Calculate optical properties for diagnostic bands:

      IF (ldiag_aeropt) THEN
       zaer_tau180_diag(1:kproma,:,:) = 0.0_dp
       zaer_lidarR_diag(1:kproma,:,:) = 0.0_dp
       pp180_diag(1:kproma,:,:,:) = 0.0_dp
       IF (ANY(nrad(:)==1) .OR. ANY(nrad(:)==3)) THEN

       DO jclass=1, nclass
          IF (nrad(jclass)==1 .OR. nrad(jclass)==3) THEN

#ifdef HAMMOZ
             rwet_p => rwet(jclass)%ptr
#endif

             DO jwv=1,kb_diag

                !--- 1.1) Calculate volume averaged refractive index nr and ni:
#ifdef HAMMOZ
                IF (ltimer) CALL timer_start(timer_ham_rad_refrac)
#endif
                !  Fill znr2d(1:kproma,:) with interpolated values; to becopied into nr_diag(1:kproma,:,jwv,jclass)
                !
                !    3.46, 2.79, 2.33, 2.05, 1.78, 1.46, 1.27, 1.01, 0.70, 0.53, 0.39, 0.30, 0.23, 8.02 [um]
                call interp_refr_index(kproma, kbdim, klev, lambda, nr(:,:,:,jclass), lambda_diag(jwv), znr2d)

                !  Fill zni2d(1:kproma,:) with interpolated values; to becopied into ni_diag(1:kproma,:,jwv,jclass)
                !
                ! should we do logarithmic interpolation for imaginary part?
                call interp_refr_index(kproma, kbdim, klev, lambda, ni(:,:,:,jclass), lambda_diag(jwv), zni2d)

#ifdef HAMMOZ
                IF (ltimer) CALL timer_stop(timer_ham_rad_refrac)
#endif  
                !--- 1.1) Calculate size parameter:
#ifdef HAMMOZ          
                zxx(1:kproma,:) = 2._dp*pi*rwet_p(1:kproma,:,krow)/lambda_diag(jwv)
#else
                zxx(1:kproma,:) = 2._dp*pi*rwet_m7(1:kproma,:,jclass)/lambda_diag(jwv) 
#endif  

                !--- 1.2) Table-lookup for optical properties:
#ifdef HAMMOZ 
                IF (ltimer) CALL timer_start(timer_ham_rad_fitplus)
#endif 
                !gf: same as in the call to ham_rad_refrac, for the call to ham_rad_fitplus

                zsigma2d(1:kproma,:) = sigma_diag(1:kproma,:,jwv,jclass)
                zomega2d(1:kproma,:) = omega_diag(1:kproma,:,jwv,jclass)
                zasym2d(1:kproma,:)  = asym_diag(1:kproma,:,jwv,jclass)

                SELECT CASE(nham_subm)

                CASE(HAM_M7)

                   IF (ABS(modesigma(jclass)-sigma_fine) < zeps*sigma_fine) THEN
                      itable=1
                   ELSE IF (ABS(modesigma(jclass)-sigma_coarse) < zeps*sigma_coarse) THEN
                      itable=2
                   ELSE 
                      CALL finish('ham_rad','incompatible standard deviation in modal setup')
                   END IF

                CASE(HAM_SALSA)

#ifdef SALSA
                   zxx(1:kproma,:) = 2._dp*pi*rwet_salsa(1:kproma,:,jclass)/lambda(jwv)

                   IF(jclass < 6 .OR. (jclass > fn2a .AND. jclass < fn2b-(nbin3-1))) THEN
                      itable=1                      
                   ELSE                     
                      itable=2                      
                   END IF
#endif
                END SELECT

                IF (itable == 1) THEN
                   CALL ham_rad_fitplus(kproma, kbdim,      klev,                &
                                         zxx,   znr2d,     zni2d,                &
                                         itable, lut1_sigma, zsigma2d,           & 
                                                 lut1_omega, zomega2d,           &
                                                 lut1_g,     zasym2d,            &
                                                 lut1_pp180, zpp180deg2d)

                ELSE 

                   CALL ham_rad_fitplus(kproma, kbdim,      klev,                & 
                                         zxx,   znr2d,     zni2d,                &
                                         itable, lut2_sigma, zsigma2d,           &
                                                 lut2_omega, zomega2d,           &
                                                 lut2_g,     zasym2d,            &
                                                 lut2_pp180, zpp180deg2d)

                END IF
#ifdef HAMMOZ
                IF (ltimer) CALL timer_stop(timer_ham_rad_fitplus)
#endif
                !>>gf: update the original 4d arrays
                sigma_diag(1:kproma,:,jwv,jclass) = zsigma2d(1:kproma,:)*lambda_diag(jwv)*lambda_diag(jwv)
                omega_diag(1:kproma,:,jwv,jclass) = zomega2d(1:kproma,:)
                asym_diag(1:kproma,:,jwv,jclass)  = zasym2d(1:kproma,:)
                pp180_diag(1:kproma,:,jwv,jclass) = zpp180deg2d(1:kproma,:)
                nr_diag(1:kproma,:,jwv,jclass)    = znr2d(1:kproma,:)
                ni_diag(1:kproma,:,jwv,jclass)    = zni2d(1:kproma,:)
                !<<gf

          !ALLOCATE(sigma(kbdim,klev,Nwv_tot,nclass))
          !ALLOCATE(omega(kbdim,klev,Nwv_sw_tot,nclass))
          !ALLOCATE(asym (kbdim,klev,Nwv_sw_tot,nclass))
          !ALLOCATE(nr(kbdim,klev,Nwv_tot,nclass))
          !ALLOCATE(ni(kbdim,klev,Nwv_tot,nclass))
                !sigma(1:kproma,:,jwv,jclass) = sigma(1:kproma,:,jwv,jclass)*lambda(jwv)*lambda(jwv)
         
             END DO ! jwv diag

          END IF ! nrad
       END DO ! jclass

       !--- Summation over modes, conversion into total extinction and calculation 
       !    of weighted single scattering albedo and assymetry factor. 

       zaer_tau_diag_vr(1:kproma,:,:,:)=0.0_dp
 
       DO jclass=1, nclass
          IF(nrad(jclass)>0) THEN
             DO jwv=1, kb_diag
                zaer_tau_diag_vr(1:kproma,:,jwv,jclass)=znum(1:kproma,:,jclass)*sigma_diag(1:kproma,:,jwv,jclass)
             END DO
          END IF
       END DO

       !--- Diagnose AOD for requested each mode (nrad) and wavelength (nraddiagwv):
!#ifdef HAMMOZ
!       DO jclass=1, nclass
!          IF(nrad(jclass)>0) THEN
!             DO jwv=1, kb_diag
!                IF (nraddiagwv(jwv)>0) THEN
!                   tau_mode(jclass,jwv)%ptr(1:kproma,:,krow)=zaer_tau_sw_vr(1:kproma,:,jwv,jclass)
!                END IF
!             END DO
!          END IF
!       END DO
!#endif
       !--- Calculation of weighted properties and vertical reordering to RRTM structure:

       DO jwv=1, kb_diag !ham_ps +Nwv_sw_opt

          DO jclass=1, nclass
             DO jk=1, klev
#ifdef HAMMOZ       
               ikl=klev+1-jk
#else       
               !No vertical reordering here
               ikl=jk
#endif               
                DO jl=1, kproma


    !real(kind=jprb),intent(inout) :: zaer_tau_diag(kbdim,klev,kb_diag,ntype_diaf)
    !real(kind=jprb),intent(inout) :: zaer_ssa_diag(kbdim,klev,kb_diag,ntype_diaf)
    !real(kind=jprb),intent(inout) :: zaer_asym_diag(kbdim,klev,kb_diag,ntype_diaf)

                   zaer_tau_diag(jl,jk,jwv)  =zaer_tau_diag(jl,jk,jwv) + &
                                            zaer_tau_diag_vr(jl,ikl,jwv,jclass)
                   zaer_ssa_diag(jl,jk,jwv)  =zaer_ssa_diag(jl,jk,jwv) + &
                                            zaer_tau_diag_vr(jl,ikl,jwv,jclass)*omega_diag(jl,ikl,jwv,jclass)
                   zaer_asym_diag(jl,jk,jwv) =zaer_asym_diag(jl,jk,jwv) + &
                                            zaer_tau_diag_vr(jl,ikl,jwv,jclass)*omega_diag(jl,ikl,jwv,jclass)*asym_diag(jl,ikl,jwv,jclass)
                   zaer_tau180_diag(jl,jk,jwv) = zaer_tau180_diag(jl,jk,jwv) + &
                                            zaer_tau_diag_vr(jl,ikl,jwv,jclass)*omega_diag(jl,ikl,jwv,jclass)*pp180_diag(jl,ikl,jwv,jclass)
                END DO
             END DO
          END DO

          DO jk=1, klev
             DO jl=1, kproma
               ! PLS - TODO: SAFE DIVISION
               IF (zaer_ssa_diag(jl,jk,jwv) > THRESHOLD) THEN 
                 zaer_asym_diag(jl,jk,jwv) = zaer_asym_diag(jl,jk,jwv)/zaer_ssa_diag(jl,jk,jwv)
                 zaer_ssa_diag(jl,jk,jwv)  = zaer_ssa_diag(jl,jk,jwv)/ zaer_tau_diag(jl,jk,jwv)
                END IF
               ! The phase function is normalized to 4*pi; guard the backscatter denominator.
               IF (zaer_tau180_diag(jl,jk,jwv) > THRESHOLD) THEN
                 zaer_lidarR_diag(jl,jk,jwv) = 4.0_dp*pi*zaer_tau_diag(jl,jk,jwv)/zaer_tau180_diag(jl,jk,jwv)
               END IF
             END DO
          END DO

       END DO

    END IF
    END IF ! ldiag_aeropt


    !--- 4) Set echam fields to zero for diagnostic aerosol radiative properties only:

    IF (naerorad==2) THEN 
       aer_tau_lw_vr(1:kproma,:,:) = 0.0_dp
       aer_tau_sw_vr(1:kproma,:,:) = 0.0_dp
       aer_piz_sw_vr(1:kproma,:,:) = 0.0_dp
       aer_cg_sw_vr(1:kproma,:,:)  = 0.0_dp
    END IF

  END SUBROUTINE ham_rad

  !----------------------------------------------------------------------------------------------------------------
  
  SUBROUTINE interp_refr_index(kproma, kbdim, klev, wl_in, nr_in, wl_out, nr_out)
    !! Interpolates refractive index nr_in(kbdim,klev,Nwv_tot)
    !! at a single wavelength wl_out, result in nr_out(kbdim,klev)

    implicit none

    ! Input
    integer,  intent(in) :: kproma, kbdim, klev
    real(dp), intent(in) :: wl_in(Nwv_tot)           ! (Nwv_tot), unsorted
    real(dp), intent(in) :: nr_in(kbdim,klev,Nwv_tot)       ! (kbdim,klev,Nwv_tot)
    real(dp), intent(in) :: wl_out             ! single query wavelength

    ! Output
    real(dp), intent(out) :: nr_out(kbdim, klev)

    ! Locals
    integer :: idx_low, idx_high
    real(dp) :: val_low, val_high, weight

    ! ---- Sort wl_in (insertion sort) ----

    ! nearest smaller/equal
    val_low = maxval(pack(wl_in, wl_in <= wl_out))
    idx_low = maxloc(wl_in, mask = wl_in == val_low, dim=1)

    ! nearest larger/equal
    val_high = minval(pack(wl_in, wl_in > wl_out))
    idx_high = maxloc(wl_in, mask = wl_in == val_high, dim=1)

    weight = (wl_out-val_low)/(val_high-val_low)

    nr_out(1:kproma,1:klev) = (1.0_dp - weight) * nr_in(1:kproma,:,idx_low)+ weight*nr_in(1:kproma,:,idx_high)

  END SUBROUTINE interp_refr_index

  !----------------------------------------------------------------------------------------------------------------
  !>>dod removed wavelength from subroutine interface
  SUBROUTINE ham_rad_fitplus(kproma, kbdim, klev,    &
                              pxx,    pnr,   pni,     &
                              ktable, plut1,  pfit1,  &
                                      plut2,  pfit2,  &
                                      plut3,  pfit3,  &
                                      plut4,  pfit4   )
  !<<dod
  
    ! *ham_rad_fitplus* returns linear interpolated fit of the 
    !                    look-up tables for the aerosol optical 
    !                    properties.
    !
    ! Authors:
    ! --------
    ! Olivier Boucher, Univ. Lille,                                  2003
    !    (original source)
    ! Philip Stier, MPI-Met, Hamburg,                          08/05/2003
    !    (adaption to ECHAM5/HAM, reduced memory expense)
    ! Luis Kornblueh, MPI-Met, Hamburg                         2006-07-07
    !    (bugfix for uninitialized variables by adding kproma)
    ! 
    ! Interface:
    ! ----------
    !

    USE mo_kind,      ONLY: dp

    IMPLICIT NONE

    !--- Arguments:

    INTEGER,  INTENT(in)  :: ktable, kproma, kbdim, klev 

    !>>dod changed arrays to 2-dimensional (removed wavelength dimension since the
    !      calling subroutine now contains the loop over the wavelengths)
    REAL(dp), INTENT(out) :: pfit1(kbdim,klev)
    REAL(dp), INTENT(out), OPTIONAL :: pfit2(kbdim,klev)
    REAL(dp), INTENT(out), OPTIONAL :: pfit3(kbdim,klev)
    REAL(dp), INTENT(out), OPTIONAL :: pfit4(kbdim,klev)

    REAL(dp), INTENT(in)  :: pnr(kbdim,klev),  pni(kbdim,klev), pxx(kbdim,klev)

    REAL(dp), INTENT(in)  :: plut1(0:Nnrmax(ktable), 0:Nnimax(ktable), 0:Ndismax(ktable))
    REAL(dp), INTENT(in), OPTIONAL  :: plut2(0:Nnrmax(ktable), 0:Nnimax(ktable), 0:Ndismax(ktable))
    REAL(dp), INTENT(in), OPTIONAL  :: plut3(0:Nnrmax(ktable), 0:Nnimax(ktable), 0:Ndismax(ktable))
    REAL(dp), INTENT(in), OPTIONAL  :: plut4(0:Nnrmax(ktable), 0:Nnimax(ktable), 0:Ndismax(ktable))

    !--- Local variables:

    LOGICAL, PARAMETER :: laerocom_diag=.FALSE. !-- Diagnostic for values out of range
    LOGICAL, PARAMETER :: loint =.FALSE. !-- Linear interpolation in look-up table

    INTEGER            :: jl,jk 

    REAL(dp)           :: zeps

    INTEGER  :: Ndis, Nnr, Nni
    !>>dod deleted security check
    REAL(dp) :: xx1, xx2, nr1, nr2, ni1, ni2
    REAL(dp) :: fitndisnr, fitndisp1nr, fitndisnrp1, fitndisp1nrp1, fitndis, fitndisp1

    !--- 0) 
    !>>dod moved calculation of inc_nr and inc_ni to mo_ham_rad_data
    !<<dod
    !--- Security check:
    !>>dod deleted: clearly untested in parallel environment, so deleted rather than modified

    !--- 1) Quick and dirty linear interpolation: 
    !
    IF(loint) THEN
       !>>dod WARNING : NEEDS UPDATING! DOES NOT WORK! 
       !>>dod removed loop over wavelengths
       DO jk=1, klev
          DO jl=1, kproma
             IF ((pxx(jl,jk)>=x0_min(ktable)) .AND. (pxx(jl,jk)<=x0_max(ktable)) .AND. &
                 (pnr(jl,jk)>=nr_min(ktable)) .AND. (pnr(jl,jk)<=nr_max(ktable)) .AND. &
                 (pni(jl,jk)>=ni_min(ktable)) .AND. (pni(jl,jk)<=ni_max(ktable))       ) THEN

                Ndis=INT( (LOG(pxx(jl,jk))-LOG(x0_min(ktable))) / &
                          (LOG(x0_max(ktable)) -LOG(x0_min(ktable)))*REAL(Ndismax(ktable),dp) )
                Ndis=MIN(Ndismax(ktable)-1,MAX(0,Ndis))

                xx1=EXP(LOG(x0_min(ktable))+(LOG(x0_max(ktable))-LOG(x0_min(ktable)))*REAL(Ndis,dp)/REAL(Ndismax(ktable),dp))
                xx2=EXP(LOG(x0_min(ktable))+(LOG(x0_max(ktable))-LOG(x0_min(ktable)))*REAL(Ndis+1,dp)/REAL(Ndismax(ktable),dp))

                Nnr=INT((pnr(jl,jk)-nr_min(ktable))/inc_nr(ktable))
                Nnr=MIN(Nnrmax(ktable)-1,MAX(0,Nnr))

                nr1=nr_min(ktable)+REAL(Nnr,dp)*inc_nr(ktable)
                nr2=nr_min(ktable)+REAL(Nnr+1,dp)*inc_nr(ktable)

                Nni=INT((LOG(pni(jl,jk))-LOG(ni_min(ktable)))/inc_ni(ktable))
                Nni=MIN(Nnimax(ktable)-1,MAX(0,Nni))

                ni1=EXP(LOG(ni_min(ktable))+REAL(Nni,dp)*inc_ni(ktable))
                ni2=EXP(LOG(ni_min(ktable))+REAL(Nni+1,dp)*inc_ni(ktable))

                fitndisnr      =plut1(Nnr,Nni,Ndis)+(pni(jl,jk)-ni1)/(ni2-ni1)*  & 
                                (plut1(Nnr,Nni+1,Ndis)-plut1(Nnr,Nni,Ndis))
                fitndisnrp1    =plut1(Nnr+1,Nni,Ndis)+(pni(jl,jk)-ni1)/(ni2-ni1)*  &
                                (plut1(Nnr+1,Nni+1,Ndis)-plut1(Nnr+1,Nni,Ndis))
                fitndisp1nr    =plut1(Nnr,Nni,Ndis+1)+(pni(jl,jk)-ni1)/(ni2-ni1)*  &
                                (plut1(Nnr,Nni+1,Ndis+1)-plut1(Nnr,Nni,Ndis+1))
                fitndisp1nrp1  =plut1(Nnr+1,Nni,Ndis+1)+(pni(jl,jk)-ni1)/(ni2-ni1)* & 
                                (plut1(Nnr+1,Nni+1,Ndis+1)-plut1(Nnr+1,Nni,Ndis+1))
                fitndis        =fitndisnr+(pnr(jl,jk)-nr1)/(nr2-nr1)*(fitndisnrp1-fitndisnr)
                fitndisp1      =fitndisp1nr+(pnr(jl,jk)-nr1)/(nr2-nr1)*(fitndisp1nrp1-fitndisp1nr)

                pfit1(jl,jk)=fitndis+(pxx(jl,jk)-xx1)/(xx2-xx1)*(fitndisp1-fitndis)

             ELSE
                
                pfit1(jl,jk)=0._dp

             END IF
          END DO
       END DO
       !<<dod

    ELSE ! .NOT. loint

       !--- 2) Quicker and dirtier table look-up:
       !>>dod removed loop over wavelengths
       DO jk=1, klev
          DO jl=1, kproma
             IF ((pxx(jl,jk)>=x0_min(ktable)) .AND. (pxx(jl,jk)<=x0_max(ktable)) .AND. &
                 (pnr(jl,jk)>=nr_min(ktable)) .AND. (pnr(jl,jk)<=nr_max(ktable)) .AND. &
                 (pni(jl,jk)>=ni_min(ktable)) .AND. (pni(jl,jk)<=ni_max(ktable))       ) THEN
                Ndis=NINT( (LOG(pxx(jl,jk))-log_x0_min(ktable)) / &
                           (log_x0_max(ktable) - log_x0_min(ktable))*REAL(Ndismax(ktable),dp) )
                Ndis=MIN(Ndismax(ktable)-1,MAX(0,Ndis))

                Nnr=NINT((pnr(jl,jk)-nr_min(ktable))/inc_nr(ktable))
                Nnr=MIN(Nnrmax(ktable)-1,MAX(0,Nnr))

                Nni=NINT((LOG(pni(jl,jk))-log_ni_min(ktable))/inc_ni(ktable))
                Nni=MIN(Nnimax(ktable)-1,MAX(0,Nni))

                pfit1(jl,jk)=plut1(Nnr,Nni,Ndis)
                
                IF (PRESENT(plut2)) pfit2(jl,jk) = plut2(Nnr, Nni,Ndis)
                IF (PRESENT(plut3)) pfit3(jl,jk) = plut3(Nnr, Nni,Ndis)
                IF (PRESENT(plut4)) pfit4(jl,jk) = plut4(Nnr, Nni,Ndis)

             ELSE
                
                pfit1(jl,jk)=0._dp
                IF (PRESENT(plut2)) pfit2(jl,jk) = 0._dp
                IF (PRESENT(plut3)) pfit3(jl,jk) = 0._dp
                IF (PRESENT(plut4)) pfit4(jl,jk) = 0._dp

             END IF
          END DO
       END DO
       !<<dod

    END IF

  END SUBROUTINE ham_rad_fitplus

!----------------------------------------------------------------------------------------------------------------

  SUBROUTINE ham_rad_initialize(nclass)

    ! *ham_rad_initialize* initializes the HAM aerosol radiation module
    !                       mo_ham_rad. It sets module variables and 
    !                       reads in the look-up tables for the optical
    !                       properties.
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-Met, Hamburg          05/2003
    ! Anton Laakso, FMI, Modifications for OIFS  02/2023
    !
    ! Interface:
    ! ----------
    ! *ham_rad_initialize* is called from *ham_initialize*
    !
#ifdef HAMMOZ
    USE mo_mpi,           ONLY: p_io, p_parallel_io, p_bcast
    USE mo_srtm_setup,    ONLY: wavenum1_sw=>wavenum1, wavenum2_sw=>wavenum2
    USE mo_rrtm_params,   ONLY: nbndlw
    USE mo_lrtm_setup,    ONLY: wavenum1, wavenum2
    USE mo_ham,           ONLY: nrad, nham_subm, HAM_SALSA, HAM_M7
    USE mo_read_netcdf77, ONLY: read_var_nf77_3d
    USE mo_util_string,   ONLY: separator
#else !OIFS
    !USE YOERAD    ,ONLY : nbndlw=> STRATO_CMIP6_NTB !this is really stupid solution
    ! commented by Lwu
    USE YOMMP0      ,ONLY : MYPROC
    USE MPL_MODULE ,ONLY : MPL_BROADCAST 
    USE TM5M7_DATA, ONLY : TM5M7_DATADIR
    USE TM5M7_OPTICS_DATA, ONLY : NASWBAND,ASWBAND,wavenum1=>ALWWN1, wavenum2=>ALWWN2
#endif
    USE mo_read_netcdf77, ONLY: read_var_nf77_3d
    USE mo_exception,     ONLY: finish, message, message_text, em_param, em_error, &
         em_info
    USE mo_species,       ONLY: aero_idx





    IMPLICIT NONE

    INTEGER, INTENT(in)       :: nclass    ! number of aerosol bins/modes

    INTEGER, PARAMETER :: jpsw = 14 !SF new rad scheme to check

    INTEGER :: jwv, iwv, jt, ierr

    LOGICAL :: lex

    CHARACTER(len=256) :: cfile

    INTEGER, PARAMETER :: RPRC = 1  !for OIFS
    INTEGER, PARAMETER :: ITAG = 98784 !for OIFS
#ifdef HAMMOZ
    CALL message('','')
    CALL message('',separator)
    CALL message('ham_rad_initialize','Parameter settings for the aerosol radiation interaction', &
                 level=em_info)
    CALL message('','',level=em_param)
#endif
   
    nhamaer = nclass

    ! construct species indices from speclist into aerosol species list
    aero_ridx(:) = 0
    DO jt = 1, naerospec
      aero_ridx(aero_idx(jt)) = jt
    END DO

    CALL ham_rad_data_initialize

    !--- 1) Set total number of wavelengths (GCM + optional):

    !--- 2) Consistency checks:

    !IF ( Nwv_sw /= jpsw .OR. Nwv_lw /= nbndlw) THEN! stratosphere aerosol should be updated
    IF ( Nwv_sw /= jpsw ) THEN

       CALL finish('ham_rad_initialize','inconsistent number of wavelengths')

    END IF

    IF (ANY(nrad(:)==1)) THEN
       IF ( SIZE(cnr,1)<Nwv_sw+Nwv_sw_opt .OR. SIZE(cni,1)<Nwv_sw+Nwv_sw_opt) THEN
          CALL finish('ham_rad_initialize', &
                      'insufficient number of refractive indices defined for nrad=1')
       END IF
    END IF
    IF (ANY(nrad(:)==2)) THEN
       IF ( SIZE(cnr,1)<Nwv_lw .OR. SIZE(cni,1)<Nwv_lw) THEN
          CALL finish('ham_rad_initialize', &
                      'insufficient number of refractive indices defined for nrad=2')
       END IF
    END IF
    IF (ANY(nrad(:)==3)) THEN
       IF ( SIZE(cnr,1)<Nwv_sw+Nwv_sw_opt+Nwv_lw .OR. SIZE(cni,1)<Nwv_sw+Nwv_sw_opt+Nwv_lw) THEN
          CALL finish('ham_rad_initialize', &
                      'insufficient number of refractive indices defined for nrad=3')
       END IF
    END IF

    !--- Allocate memory for look-up tables:
    !    (needs to be done before read-in, ham_init_memory is called after)

    CALL ham_rad_mem

    !--- 2) Read in look-up tables:

    IF (ANY(nrad(:)==1) .OR. ANY(nrad(:)==3)) THEN
#ifdef HAMMOZ
       IF (p_parallel_io) THEN
          cfile='lut_optical_properties.nc'
#else
          IF (MYPROC==RPRC) THEN !alaakso MUUTA TAMA
             cfile=TRIM(TM5M7_DATADIR)//'lut_optical_properties_M7.nc'
#endif

             CALL message('ham_rad_initialize', 'Reading lookup table from '//TRIM(ADJUSTL(cfile)), level=em_info)
             INQUIRE (file=cfile,exist=lex)

             IF (lex) THEN
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "sigma_1", lut1_sigma, ierr                 )
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "sigma_2", lut2_sigma, ierr                 )
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "omega_1", lut1_omega, ierr                 )
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "omega_2", lut2_omega, ierr                 )
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "asym_1",  lut1_g,     ierr                 )
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "asym_2",  lut2_g,     ierr                 )
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "pp180_1", lut1_pp180, ierr                 )
                CALL read_var_nf77_3d (cfile,     "nr",       "ni",     "dis",     &
                     "pp180_2", lut2_pp180, ierr                 )

             ELSE
                CALL message('ham_rad_initialize','file '//TRIM(ADJUSTL(cfile))//' not available', level=em_error)
                CALL finish('ham_rad_initialize','file '//TRIM(ADJUSTL(cfile))//' missing!',1)
             END IF

          END IF
#ifdef HAMMOZ
          CALL p_bcast(lut1_sigma,p_io)
          CALL p_bcast(lut1_g    ,p_io)
          CALL p_bcast(lut1_omega,p_io)
          CALL p_bcast(lut1_pp180,p_io)
          CALL p_bcast(lut2_sigma,p_io)
          CALL p_bcast(lut2_g    ,p_io)
          CALL p_bcast(lut2_omega,p_io)
          CALL p_bcast(lut2_pp180,p_io)
#else
          CALL MPL_BROADCAST(lut1_sigma,KTAG=ITAG+1,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-1: ')
          CALL MPL_BROADCAST(lut1_g,KTAG=ITAG+2,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-2: ')
          CALL MPL_BROADCAST(lut1_omega,KTAG=ITAG+3,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-3: ')
          CALL MPL_BROADCAST(lut1_pp180,KTAG=ITAG+4,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-4: ')
          CALL MPL_BROADCAST(lut2_sigma,KTAG=ITAG+5,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-5: ')
          CALL MPL_BROADCAST(lut2_g,KTAG=ITAG+6,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-6: ')
          CALL MPL_BROADCAST(lut2_omega,KTAG=ITAG+7,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-7: ')
          CALL MPL_BROADCAST(lut2_pp180,KTAG=ITAG+8,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-8: ')
#endif
       END IF

       IF (ANY(nrad(:)==2) .OR. ANY(nrad(:)==3)) THEN
#ifdef HAMMOZ
          IF (p_parallel_io) THEN
             cfile='lut_optical_properties_lw.nc'
#else
             IF (MYPROC==RPRC) THEN !alaakso MUUTA TAMA
                cfile=TRIM(TM5M7_DATADIR)//'lut_optical_properties_lw_M7.nc'
#endif
                CALL message('ham_rad_initialize', 'Reading lookup table from '//TRIM(ADJUSTL(cfile)), level=em_info)
                INQUIRE (file=cfile,exist=lex)
                IF (lex) THEN
                   CALL read_var_nf77_3d (cfile,        "nr",       "ni",     "dis",     &
                        "sigma_1_lw", lut3_sigma, ierr                 )
                   CALL read_var_nf77_3d (cfile,        "nr",       "ni",     "dis",     &
                        "sigma_2_lw", lut4_sigma, ierr                 )
                ELSE
                   CALL message('ham_rad_initialize','file '//TRIM(ADJUSTL(cfile))//' not available', level=em_error)
                   CALL finish('ham_rad_initialize','file '//TRIM(ADJUSTL(cfile))//' missing!',1)
                END IF

             END IF
#ifdef HAMMOZ
             CALL p_bcast(lut3_sigma,p_io)
             CALL p_bcast(lut4_sigma,p_io)
#else
             CALL MPL_BROADCAST(lut3_sigma,KTAG=ITAG+9,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-9: ')
             CALL MPL_BROADCAST(lut4_sigma,KTAG=ITAG+10,KROOT=RPRC,CDSTRING='AER_M7orS_PROP-10: ')
#endif
          END IF



          !--- 3) SW initializations:

          IF (ANY(nrad==1) .OR. ANY(nrad==3)) THEN 

             !---  Interpolate wavelengths [m] from RRTM-SW wavenumbers [cm-1]:
             !     (Hardcoded to RRTM-SW wavenumber indices)

             CALL message('','',level=em_param)
             CALL message('', ' SW wavelengths [um] : ', level=em_param)
#ifdef HAMMOZ
             DO jwv=1, Nwv_sw
                !--- Mid-wavelength (mid-wavenumber wavelength):
                !lambda(jwv)= 2._dp/(wavenum1_sw(15+jwv)+wavenum2_sw(15+jwv)) * 0.01_dp
                !--- Mid-wavelength (mid-wavelength wavelength):
                lambda(jwv)=(1._dp/wavenum1_sw(15+jwv)+1._dp/wavenum2_sw(15+jwv))/2._dp * 0.01_dp
                WRITE(message_text,fmt='(a,i3,a,f8.2)') '      lambda(', jwv, ') = ', lambda(jwv)*1.E6_dp
                CALL message('', message_text, level=em_param)
             ENDDO

             !--- Add optional SW wavelengths
             CALL message('','',level=em_param)
             CALL message('', ' SW wavelengths (optional) [um] : ', level=em_param)

             DO jwv=1, Nwv_sw_opt
                iwv=Nwv_sw+jwv
                lambda(iwv)=lambda_sw_opt(jwv)
                WRITE(message_text,fmt='(a,i3,a,f8.2)') '      lambda(', iwv, ') = ', lambda(iwv)*1.E6_dp
                CALL message('', message_text, level=em_param)
             END DO
#else

             DO jwv=1, Nwv_sw  !Laakso: note different order than HAM (here same as RRTM)        
                lambda(jwv)=ASWBAND(jwv)%wl*1.E-6_dp
             END DO
             !!! include diagnostic wavelengths, lhw
             DO jwv=1, Nwv_sw_opt
                iwv=Nwv_sw+jwv
                lambda(iwv)=lambda_sw_opt(jwv)
                WRITE(message_text,fmt='(a,i3,a,f8.2)') '      lambda(', iwv, ') = ', lambda(iwv)*1.E6_dp
                CALL message('', message_text, level=em_param)
             END DO
#endif
          END IF

          !--- 4) LW initializations:

          !--- Set up RRTM bands (convert from wavenumber [cm-1] to wavelength [m]:
          IF (ANY(nrad==2) .OR. ANY(nrad==3)) THEN

             CALL message('','',level=em_param)
             CALL message('', ' LW wavelengths [um] : ', level=em_param)

             DO jwv=1, Nwv_lw
                iwv=Nwv_sw+Nwv_sw_opt+jwv
                !--- Mid-wavelength (mid-wavenumber wavelength):
                lambda(iwv)=2.0_dp/(wavenum1(jwv)+wavenum2(jwv)) * 0.01_dp
                WRITE(message_text,fmt='(a,i3,a,f8.2)') '      lambda(', iwv, ') = ', lambda(iwv)*1.E6_dp
                CALL message('', message_text, level=em_param)
             END DO

          END IF
#ifdef HAMMOZ
          CALL message('','',level=em_param)
          CALL message('',separator)
#endif


        END SUBROUTINE ham_rad_initialize

        !----------------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------------
#ifdef HAMMOZ
        SUBROUTINE ham_rad_diag(kproma, kbdim, klev, krow, pxtm1)

          ! *ham_rad_diag* calculated type specific diagnostics
          !                 of aerosol optical properties
          !
          ! Authors:
          ! --------
          ! Philip Stier, MPI-Met, Hamburg,       05/2003
          ! 
          ! Modified:
          ! ---------
          ! Declan O'Donnell MPI-Met Hamburg : modifications for SOA
          ! Interface:
          ! ----------
          ! *ham_rad_diag* is called from *radiation*


          USE mo_kind,         ONLY: dp
          USE mo_tracdef,      ONLY: ntrac 
          USE mo_ham,          ONLY: nrad 
          USE mo_ham_rad_data, ONLY: nraddiagwv, nradang
          USE mo_ham_streams,  ONLY: tau_mode, abs_mode, omega_mode, sigma_mode, asym_mode,    &
               tau_comp, abs_comp, tau_2d,     abs_2d,     ang,          &
               nr_mode,    ni_mode,                  &
               omega_2d_mode, sigma_2d_mode, asym_2d_mode,               &
               nr_2d_mode,    ni_2d_mode
          USE mo_species,      ONLY: aero_idx, speclist

          IMPLICIT NONE

          INTEGER   :: kproma, kbdim, klev, krow

          REAL(dp)  :: pxtm1(kbdim,klev,ntrac)

          INTEGER   :: jl, jk, jclass, ikey, jt

          REAL(dp)  :: zv, zdensity, zeps

          REAL(dp)  :: ztaucomp(kbdim), zabscomp(kbdim)

          REAL(dp)  :: zomega(kbdim), zsigma(kbdim), zasym(kbdim), &
               znr(kbdim),    zni(kbdim),    ztau(kbdim),  &
               znr_2d(kbdim), zni_2d(kbdim)

          REAL(dp)  :: zvsum(kbdim,klev,nclass), znivsum(kbdim,klev,nclass)

          REAL(dp)  :: zvcomp(kbdim,klev,naerospec,nclass)

          REAL(dp), POINTER :: tau_2d_p(:,:)
          REAL(dp), POINTER :: tau_p(:,:,:),    abs_p(:,:,:)

          REAL(dp) :: ztmp1(kbdim), ztmp2(kbdim) !SF #458 temporary vars.

          LOGICAL  :: ll1(kbdim)   !SF #458 temporary var.

          INTEGER :: jspec, jwv

          !--- 0)

          zeps = EPSILON(1.0_dp)  ! Usage must be checked. Not done yet since we do not use HAMMOZ

          !--- Optical thickness for optional wavelengths:

          DO jwv=1, Nwv_tot

             IF ( nraddiagwv(jwv) > 0 ) THEN

                tau_2d_p => tau_2d(jwv)%ptr
                tau_2d_p(1:kproma,krow) = 0._dp

                DO jclass=1, nclass
                   IF( nrad(jclass) > 0 )THEN

                      tau_p    => tau_mode(jclass,jwv)%ptr

                      !--- Optical thickness per mode at optional wavelengths:

                      tau_p(1:kproma,:,krow) = znum(1:kproma,:,jclass)*sigma(1:kproma,:,jwv,jclass)

                      !--- 2) Vertical integral summed over all modes:


                      DO jk=1, klev
                         DO jl=1, kproma
                            tau_2d_p(jl,krow)=tau_2d_p(jl,krow)+tau_p(jl,jk,krow)
                         END DO
                      END DO

                   END IF

                END DO

             END IF ! nraddiagwv(jwv)>0

             !--- 2) 2D extended diatnostics of mode radiative parameters:

             IF (nraddiagwv(jwv)>1) THEN

                IF (nraddiag>0) THEN

                   !--- Integrate mode radiative properties to 2D and weight with aerosol optical depth:

                   znr_2d(1:kproma)=0.0_dp
                   zni_2d(1:kproma)=0.0_dp

                   DO jclass=1, nclass
                      IF( nrad(jclass) > 0 )THEN

                         tau_p => tau_mode(jclass,jwv)%ptr

                         zomega(1:kproma)=0.0_dp
                         zsigma(1:kproma)=0.0_dp
                         zasym(1:kproma) =0.0_dp
                         znr(1:kproma)   =0.0_dp
                         zni(1:kproma)   =0.0_dp
                         ztau(1:kproma)  =0.0_dp

                         DO jk=1, klev
                            DO jl=1, kproma
                               zomega(jl)=zomega(jl)+omega(jl,jk,jwv,jclass)*tau_p(jl,jk,krow)
                               zsigma(jl)=zsigma(jl)+sigma(jl,jk,jwv,jclass)*tau_p(jl,jk,krow)
                               zasym(jl) =zasym(jl) +asym(jl,jk,jwv,jclass) *tau_p(jl,jk,krow)
                               znr(jl)   =znr(jl)   +nr(jl,jk,jwv,jclass)   *tau_p(jl,jk,krow)
                               zni(jl)   =zni(jl)   +ni(jl,jk,jwv,jclass)   *tau_p(jl,jk,krow)
                               znr_2d(jl)=znr_2d(jl)+nr(jl,jk,jwv,jclass)   *tau_p(jl,jk,krow)
                               zni_2d(jl)=zni_2d(jl)+ni(jl,jk,jwv,jclass)   *tau_p(jl,jk,krow)
                               ztau(jl)  =ztau(jl)  +tau_p(jl,jk,krow)
                            END DO
                         END DO

                         !>>SF #458 (replacing WHERE statements)
                         ll1(1:kproma) = (ztau(1:kproma) > zeps)
                         ztmp1(1:kproma) = MERGE(ztau(1:kproma), 1._dp, ll1(1:kproma)) !SF 1._dp is a dummy val.

                         ztmp1(1:kproma) = 1._dp / ztmp1(1:kproma)

                         omega_2d_mode(jclass,jwv)%ptr(1:kproma,krow) = &
                              MERGE(zomega(1:kproma)*ztmp1(1:kproma), 0._dp, ll1(1:kproma))

                         sigma_2d_mode(jclass,jwv)%ptr(1:kproma,krow) = &
                              MERGE(zsigma(1:kproma)*ztmp1(1:kproma), 0._dp, ll1(1:kproma))

                         asym_2d_mode(jclass,jwv)%ptr(1:kproma,krow) = &
                              MERGE(zasym(1:kproma)*ztmp1(1:kproma), 0._dp, ll1(1:kproma))

                         nr_2d_mode(jclass,jwv)%ptr(1:kproma,krow) = &
                              MERGE(znr(1:kproma)*ztmp1(1:kproma), 0._dp, ll1(1:kproma))

                         ni_2d_mode(jclass,jwv)%ptr(1:kproma,krow) = &
                              MERGE(zni(1:kproma)*ztmp1(1:kproma), 0._dp, ll1(1:kproma))

                         !<<SF #458 (replacing WHERE statements)
                      END IF
                   END DO

                END IF ! nraddiag

                !--- 3) 3D extended diatnostics of mode radiative parameters:

                IF (nraddiag==2) THEN
                   DO jclass=1, nclass
                      IF( nrad(jclass) > 0 )THEN
                         omega_mode(jclass,jwv)%ptr(1:kproma,:,krow)=omega(1:kproma,:,jwv,jclass)
                         sigma_mode(jclass,jwv)%ptr(1:kproma,:,krow)=sigma(1:kproma,:,jwv,jclass)
                         asym_mode(jclass,jwv)%ptr(1:kproma,:,krow) =asym(1:kproma,:,jwv,jclass)
                         nr_mode(jclass,jwv)%ptr(1:kproma,:,krow)   =nr(1:kproma,:,jwv,jclass)
                         ni_mode(jclass,jwv)%ptr(1:kproma,:,krow)   =ni(1:kproma,:,jwv,jclass)
                      END IF
                   END DO
                END IF

                !--- 4) Calculate absorption optical depth:

                abs_2d(jwv)%ptr(1:kproma,krow)=0.0_dp

                DO jclass=1, nclass
                   IF(nrad(jclass)>0)THEN

                      abs_p   => abs_mode(jclass,jwv)%ptr
                      tau_p   => tau_mode(jclass,jwv)%ptr

                      !--- For each mode:

                      abs_p(1:kproma,:,krow)  =(1.0_dp-omega(1:kproma,:,jwv,jclass))*tau_p(1:kproma,:,krow)

                      !--- Total vertical integral:

                      DO jk=1, klev
                         abs_2d(jwv)%ptr(1:kproma,krow)=abs_2d(jwv)%ptr(1:kproma,krow)+abs_p(1:kproma,jk,krow)
                      END DO

                   END IF
                END DO

                !--- 5) Split up according to compounds 
                !       (based on volume average for optical thickness,
                !        additionally weighted with ni for absorption ):

                zvcomp(1:kproma,:,:,:)  = 0._dp
                zvsum(1:kproma,:,:)     = 0._dp
                znivsum(1:kproma,:,:)   = 0._dp

                DO jclass=1,nclass
                   DO jspec=1,naerospec

                      !--- Check if species jspec exists in mode jclass:

                      !ham_ps: this has been included until consistent definition of aerosol water as species within HAM / M7:
                      !        iaerocomp of water is currently set to -1 in mo_ham_init
                      IF (speclist(aero_idx(jspec))%iaerocomp(jclass)>0) THEN

                         jt=aerocomp( speclist(aero_idx(jspec))%iaerocomp(jclass) )%idt

                         zdensity=speclist(aero_idx(jspec))%density

                         IF (nrad(jclass)>0) THEN

                            !                jspec = aerocomp(jn)%spid
                            !                ikey = aerocomp(jn)%species%iaerorad
                            ikey=speclist(aero_idx(jspec))%iaerorad

                            !--- Sum volume of compound weighted by volume:

                            DO jk=1, klev
                               DO jl=1, kproma
                                  IF(pxtm1(jl,jk,jt)>zeps) THEN

                                     !ham_ps:redundant
                                     !                         zv=pxtm1(jl,jk,jt)*zmassfac/zdensity
                                     zv=pxtm1(jl,jk,jt)/zdensity

                                     zvcomp(jl,jk,jspec,jclass)=zvcomp(jl,jk,jspec,jclass) + zv
                                     zvsum(jl,jk,jclass)       =zvsum(jl,jk,jclass)        + zv
                                     znivsum(jl,jk,jclass)     =znivsum(jl,jk,jclass)      + zv * cni(jwv,ikey)

                                  END IF
                               END DO
                            END DO
                         END IF

                      END IF

                   END DO
                END DO

                ! add aerosol water

                DO jclass=1,nclass
                   IF (nrad(jclass) > 0 .AND. sizeclass(jclass)%lsoluble) THEN
                      jt = aerowater(jclass)%idt
                      ikey = speclist(id_wat)%iaerorad
                      zdensity = speclist(id_wat)%density

                      DO jk=1, klev
                         DO jl=1, kproma
                            zv=pxtm1(jl,jk,jt)/zdensity

                            zvcomp(jl,jk,aero_ridx(id_wat),jclass)=zvcomp(jl,jk,aero_ridx(id_wat),jclass) + zv
                            zvsum(jl,jk,jclass)       =zvsum(jl,jk,jclass)        + zv
                            znivsum(jl,jk,jclass)     =znivsum(jl,jk,jclass)      + zv * cni(jwv,ikey)
                         END DO
                      END DO
                   END IF
                END DO

                DO jspec=1,naerospec

                   ztaucomp(1:kproma)    =0._dp
                   zabscomp(1:kproma)    =0._dp
                   !ham_ps:why m7 specific?
                   !ikey = speclist(subm_aerospec(jspec))%iaerorad
                   ikey=speclist(aero_idx(jspec))%iaerorad

                   !--- Weighted averaging and vertical integration:

                   DO jclass=1, nclass
                      IF (nrad(jclass) > 0) THEN

                         tau_p     => tau_mode(jclass,jwv)%ptr
                         abs_p     => abs_mode(jclass,jwv)%ptr

                         DO jk=1, klev
                            DO jl=1, kproma
                              IF (zvsum(jl,jk,jclass)>zeps) THEN
                                  ztaucomp(jl)=ztaucomp(jl) + &
                                       tau_p(jl,jk,krow)*zvcomp(jl,jk,jspec,jclass)/zvsum(jl,jk,jclass)
                                  zabscomp(jl)=zabscomp(jl) + &
                                       abs_p(jl,jk,krow)*zvcomp(jl,jk,jspec,jclass)*cni(jwv,ikey) / &
                                       znivsum(jl,jk,jclass)
                               END IF
                            END DO
                         END DO

                      END IF
                   END DO     !jclass

                   !--- Store in output streams:

                   tau_comp(jspec,jwv)%ptr(1:kproma,krow)=ztaucomp(1:kproma)
                   abs_comp(jspec,jwv)%ptr(1:kproma,krow)=zabscomp(1:kproma)

                END DO     !jspec   

             END IF !nraddiagwv(jwv)>1

          END DO !jwv

          !--- 5) Calculate Angstroem parameter between two wavelengths:

          IF (nradang(1)/=0 .AND. nradang(2)/=0) THEN 

             !>>SF #458 (replacing WHERE statements)
             ll1(1:kproma) = (tau_2d(nradang(1))%ptr(1:kproma,krow)>zeps) &
                  .AND. (tau_2d(nradang(2))%ptr(1:kproma,krow)>zeps)

             ztmp1(1:kproma) = MERGE(tau_2d(nradang(1))%ptr(1:kproma,krow), 1._dp, ll1(1:kproma)) !SF 1. is dummy
             ztmp2(1:kproma) = MERGE(tau_2d(nradang(2))%ptr(1:kproma,krow), 1._dp, ll1(1:kproma)) !SF 1. is dummy

             ang(1:kproma,krow) = MERGE( &
                  LOG(ztmp2(1:kproma)/ztmp1(1:kproma)) / LOG(lambda(nradang(1))/lambda(nradang(2))), &
                  ang(1:kproma,krow), &                           
                  ll1(1:kproma))

             !<<SF #458 (replacing WHERE statements)

          END IF

        END SUBROUTINE ham_rad_diag
#endif
        !----------------------------------------------------------------------------------------------------------------

        !ham_ps:radiation This should be part of the echam-6 (optional) standard diagnostics
        !                 Talk to Bjorn

!!$  SUBROUTINE ham_rad_diag_clearsky(kbdim,  klev,   kidia, kfdia, knu, krow, &
!!$                                    ptauaz, ppizaz, pcgaz                    )
!!$
!!$    ! *ham_rad_diag_clearsky* diagnoses applied clear-sky aerosol 
!!$    !                          radiative properties
!!$    !
!!$    ! Authors:
!!$    ! --------
!!$    ! Philip Stier, MPI-Met, Hamburg,       05/2003
!!$    !
!!$    ! Interface:
!!$    ! ----------
!!$    ! *ham_rad_diag_clearsky* is called from *swclr*
!!$
!!$
!!$    USE mo_kind,         ONLY: dp
!!$    USE mo_ham_rad_data,ONLY: ivis
!!$    USE mo_ham_rad_mem, ONLY: tau_sw_2d, abs_sw_2d
!!$
!!$    INTEGER,  INTENT(in) :: kbdim, klev, kidia, kfdia, knu
!!$
!!$    REAL(dp), INTENT(in) :: ptauaz(kbdim,klev), &
!!$                            ppizaz(kbdim,klev), &
!!$                            pcgaz(kbdim,klev)
!!$
!!$    INTEGER :: jk,jl,krow
!!$
!!$    IF (knu==ivis) THEN 
!!$
!!$      !--- 0) Initialization:
!!$
!!$      tau_sw_2d(kidia:kfdia,krow)=0.0_dp
!!$      abs_sw_2d(kidia:kfdia,krow)=0.0_dp
!!$
!!$      !--- SW visible (1st band):
!!$
!!$      DO jk=1, klev
!!$        DO jl=kidia, kfdia
!!$          tau_sw_2d(jl,krow)=tau_sw_2d(jl,krow)+ptauaz(jl,jk)
!!$          abs_sw_2d(jl,krow)=abs_sw_2d(jl,krow)+ptauaz(jl,jk)*(1._dp-ppizaz(jl,jk))
!!$        END DO
!!$      END DO
!!$
!!$    END IF
!!$
!!$  END SUBROUTINE ham_rad_diag_clearsky

        !----------------------------------------------------------------------------------------------------------------

        !----------------------------------------------------------------------------------------------------------------

        SUBROUTINE ham_rad_cache(kbdim,klev)

          ! *ham_rad_mem* allocates local (blocked) memory for the 
          !                HAM aerosol optical properties
          !
          ! Author:
          ! -------
          ! Philip Stier, MPI-Met, Hamburg          05/2003
          !
          ! Interface:
          ! ----------
          ! *ham_rad_cache* is called from *radiation*

          IMPLICIT NONE

          INTEGER, INTENT(in) :: kbdim,klev

          !ALLOCATE(sigma(kbdim,klev,Nwv_tot,nclass))
          !ALLOCATE(omega(kbdim,klev,Nwv_sw_tot,nclass))
          !ALLOCATE(asym (kbdim,klev,Nwv_sw_tot,nclass))
          !ALLOCATE(nr(kbdim,klev,Nwv_tot,nclass))
          !ALLOCATE(ni(kbdim,klev,Nwv_tot,nclass))

          !ALLOCATE(znum(kbdim,klev,nclass))
           RETURN

        END SUBROUTINE ham_rad_cache

        !----------------------------------------------------------------------------------------------------------------

        SUBROUTINE ham_rad_cache_cleanup

          ! *ham_rad_mem_cleanup* de-allocates local (blocked) memory for the 
          !                        HAM aerosol optical properties
          !
          ! Author:
          ! -------
          ! Philip Stier, MPI-Met, Hamburg          05/2003
          !
          ! Interface:
          ! ----------
          ! *ham_rad_cache_cleanup* is called from *radiation*

          IMPLICIT NONE

          !DEALLOCATE(sigma,omega,asym, nr, ni, znum)

          RETURN
        END SUBROUTINE ham_rad_cache_cleanup

        !----------------------------------------------------------------------------------------------------------------
        SUBROUTINE ham_rad_mem

          ! *ham_rad_mem* allocates memory for the 
          !                HAM aerosol optical properties
          !
          ! Author:
          ! -------
          ! Philip Stier, MPI-Met, Hamburg          05/2003
          !
          ! Interface:
          ! ----------
          ! *ham_rad_mem* is called from *init_subm_memory*
          ! in *mo_submodel_interface*

          USE mo_ham,  ONLY: nrad

          IMPLICIT NONE

          IF (ANY(nrad(:)==1) .OR. ANY(nrad(:)==3))THEN

             ALLOCATE(lut1_sigma(0:Nnrmax(1), 0:Nnimax(1), 0:Ndismax(1)), &
                  lut1_g    (0:Nnrmax(1), 0:Nnimax(1), 0:Ndismax(1)), &
                  lut1_omega(0:Nnrmax(1), 0:Nnimax(1), 0:Ndismax(1)), &
                  lut1_pp180(0:Nnrmax(1), 0:Nnimax(1), 0:Ndismax(1))  )

             ALLOCATE(lut2_sigma(0:Nnrmax(2), 0:Nnimax(2), 0:Ndismax(2)), &
                  lut2_g    (0:Nnrmax(2), 0:Nnimax(2), 0:Ndismax(2)), &
                  lut2_omega(0:Nnrmax(2), 0:Nnimax(2), 0:Ndismax(2)), &
                  lut2_pp180(0:Nnrmax(2), 0:Nnimax(2), 0:Ndismax(2))  )

          END IF

          IF (ANY(nrad(:)==2) .OR. ANY(nrad(:)==3))THEN

             ALLOCATE(lut3_sigma(0:Nnrmax(3), 0:Nnimax(3), 0:Ndismax(3)), &
                  lut4_sigma(0:Nnrmax(4), 0:Nnimax(4), 0:Ndismax(4))  )

          END IF

        END SUBROUTINE ham_rad_mem

        !----------------------------------------------------------------------------------------------------------------

        SUBROUTINE ham_rad_mem_cleanup

          ! *ham_rad_cleanup* de-allocates memory of the 
          !                    module mo_ham_rad
          !
          ! Author:
          ! -------
          ! Philip Stier, MPI-Met, Hamburg          05/2003
          !
          ! Interface:
          ! ----------
          ! *ham_rad_cleanup* is called from *free_subm_memory*
          ! in *mo_submodel_interface*

          USE mo_ham,  ONLY: nrad

          IMPLICIT NONE

          IF (ANY(nrad(:)==1) .OR. ANY(nrad(:)==3))THEN

             DEALLOCATE(lut1_sigma, lut1_g, lut1_omega, lut1_pp180, &
                  lut2_sigma, lut2_g, lut2_omega, lut2_pp180  )

          END IF

          IF (ANY(nrad(:)==2) .OR. ANY(nrad(:)==3))THEN

             DEALLOCATE(lut3_sigma, lut4_sigma)

          END IF

        END SUBROUTINE ham_rad_mem_cleanup

      END MODULE mo_ham_rad
