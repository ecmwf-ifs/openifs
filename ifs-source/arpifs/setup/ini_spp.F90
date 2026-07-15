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

SUBROUTINE INI_SPP(YDGEOMETRY,YDRIP,YDCONF,YDSPP)

! Purpose :
! -------
!    Initialises the SPP scheme, which represents model uncertainties
!    with stochastically perturbed parameterisations.
! 
!    Pointers in YDSPP are allocated and the variables are initialised.

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------

! Author :
! ------
!    M. Leutbecher (ECMWF)
!    Original : December 2014

! Modifications :
! -------------
!  
!    SJ Lock  :  Jan-2016 Added perturbations for clear-skies heating rates
!    S. Lang  :  Feb-2017 Added modification to cycle pattern in data assimilation
!    SJ Lock  :  Oct-2017 Option to reduce frequency of SPP pattern updates
!    SJ Lock  :  Oct-2017 Enabled options for new SPP microphysics perturbations
!    M Leutbecher: Oct-2020 SPP abstraction
!    Ulf Andrae : Sep-2022 Populate lam settings
!    M Leutbecher: Aug-2024 abstract level type (AL) for random fields  
!  
!-----------------------------------------------------------------------------
USE GEOMETRY_MOD        , ONLY : GEOMETRY
USE PARKIND1            , ONLY : JPIM, JPRB
USE YOMGRIB             , ONLY : NENSFNB 
USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST              , ONLY : RA
USE YOMLUN              , ONLY : NULOUT
USE YOMRIP              , ONLY : TRIP
USE SPECTRAL_ARP_MOD    , ONLY : ALLOCATE_ARP, SET_ARP2D, EVOLVE_ARP  ! , EVOLVE_ARP_SPG removed by someone
USE SPECTRAL_FIELDS_MOD , ONLY : SPECTRAL_NORM
USE GRIDPOINT_FIELDS_MIX, ONLY : ALLOCATE_GRID
USE SPP_GEN_MOD         , ONLY : JPMXSCALES
USE SPP_MOD             , ONLY : JPSPPLABLEN,TSPP_CONFIG, TSPP_DATA, KGET_SEED_SPP

!     ------------------------------------------------------------------
IMPLICIT NONE

TYPE(GEOMETRY)     , INTENT(INOUT) :: YDGEOMETRY
TYPE(TRIP)         , INTENT(INOUT) :: YDRIP
TYPE(TSPP_CONFIG)  , INTENT(INOUT) :: YDCONF
TYPE(TSPP_DATA)    , INTENT(INOUT) :: YDSPP


INTEGER(KIND=JPIM) :: I2D
INTEGER(KIND=JPIM) :: JPERT, JRF, JSCALE
INTEGER(KIND=JPIM) :: IRF
INTEGER(KIND=JPIM) :: ISEED
INTEGER(KIND=JPIM), DIMENSION(JPMXSCALES) :: IGRIB
INTEGER(KIND=JPIM) :: JN
INTEGER(KIND=JPIM) :: IN0

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB), DIMENSION(JPMXSCALES) :: ZTAU, ZLCOR, ZSDEV
REAL(KIND=JPRB), DIMENSION(:),     ALLOCATABLE:: ZCHI
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE:: ZSDEV2D, ZPHI2D
REAL(KIND=JPRB) :: ZPHI, ZWC_KAPPA_T, ZGAMMAN, ZCONSTF0

CHARACTER(LEN=JPSPPLABLEN) :: CLLABEL
CHARACTER(LEN=3)           :: CLJRF
CHARACTER(LEN=1024) :: CLTEMP
!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "read_spec_grib.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INI_SPP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDMP=>YDGEOMETRY%YRMP)

IF (YDCONF%LSPP) THEN
  I2D=YDCONF%SM%NRFTOTAL ! for convencience
  WRITE(NULOUT,'(''INI_SPP: number of 2D-patterns = '',I3)') YDCONF%SM%NRFTOTAL

  IF (ASSOCIATED(YDSPP%SP_ARP)) DEALLOCATE(YDSPP%SP_ARP)
  ALLOCATE(YDSPP%SP_ARP(I2D))

  IF(YDCONF%LSPG_SPP) THEN
    !IF (ASSOCIATED(YDSPP%SP_ARP_SPG)) DEALLOCATE(YDSPP%SP_ARP_SPG)
    ! ALLOCATE(YDSPP%SP_ARP_SPG(I2D,3)) Removed by someone
  ENDIF


  IF (ASSOCIATED(YDSPP%GP_ARP )) DEALLOCATE(YDSPP%GP_ARP )
  ALLOCATE(YDSPP%GP_ARP (I2D))
  IF (ASSOCIATED(YDSPP%GP_ARP0)) DEALLOCATE(YDSPP%GP_ARP0)
  ALLOCATE(YDSPP%GP_ARP0(I2D))
  IF (ASSOCIATED(YDSPP%GP_ARP1)) DEALLOCATE(YDSPP%GP_ARP1)
  ALLOCATE(YDSPP%GP_ARP1(I2D))

  IF (ASSOCIATED(YDSPP%LAB_ARP)) DEALLOCATE(YDSPP%LAB_ARP)
  ALLOCATE(YDSPP%LAB_ARP(I2D))  

  WRITE(NULOUT,'(''INI_SPP: allocated YDSPP%SP_ARP, YDSPP%GP_ARP, YDSPP%LAB_ARP'',I3)')

  ! set labels
  YDSPP%LAB_ARP(:)="_unset_"
  DO JPERT=1, YDCONF%SM%NACT
    IRF = YDCONF%SM%PN(JPERT)%MP - 1 ! points before random field for this perturbations
    DO JRF=1, YDCONF%SM%PN(JPERT)%NRF
      IRF = IRF + 1
      YDSPP%LAB_ARP( IRF ) = YDCONF%SM%PN(JPERT)%LABEL
    ENDDO
  ENDDO

  IF  (ANY(YDSPP%LAB_ARP(:)=="_unset_")) THEN
    CALL ABOR1("All elements of YDSPP%LAB_ARP(:) have to be initialised.")
  ENDIF
  ! allocate spectral ARPs
  ALLOCATE(ZCHI(0:YDDIM%NSMAX))
  ALLOCATE(ZSDEV2D(0:YDDIM%NSMAX,JPMXSCALES))
  ALLOCATE(ZPHI2D( 0:YDDIM%NSMAX,JPMXSCALES))

  ALLOCATE(YDSPP%IGRIBCODE(I2D))

  DO JPERT=1, YDCONF%SM%NACT
    IRF  = YDCONF%SM%PN(JPERT)%MP - 1 ! points before random field for this perturbation
    ZTAU(:) = YDCONF%SM%PN(JPERT)%TAU(:)
    ZLCOR(:)= YDCONF%SM%PN(JPERT)%XLCOR(:)
    ZSDEV(:)= YDCONF%SM%PN(JPERT)%SDEV(:)
    DO JRF=1, YDCONF%SM%PN(JPERT)%NRF
      IRF = IRF + 1
      WRITE( CLJRF, '(I3)') JRF
      DO JSCALE=1, JPMXSCALES
        !
        !   determine ARP characteristics for this random field
        !
        ZCHI(:)    =  0._JPRB
        ZWC_KAPPA_T= 0.5_JPRB * (ZLCOR(JSCALE)/RA)**2
        DO JN=0,YDDIM%NSMAX
          !ML-note: 0.5 has no _JPRB, fix this!
          ZCHI(JN)=EXP(-0.5_JPRB * ZWC_KAPPA_T * REAL(JN*(JN+1),KIND=JPRB))
        ENDDO

        ZGAMMAN  =  0._JPRB
        !ML-fix-begin: a quick-and-dirty fix that should be cleaned up
        IF ( ZLCOR(JSCALE) > RA) THEN
          IN0=0
        ELSE
          IN0=1
        ENDIF
        !ML-fix-end
        DO JN=IN0,YDDIM%NSMAX
          ZGAMMAN= ZGAMMAN + (2*JN+1)*ZCHI(JN)**2
        ENDDO
        ZPHI         = EXP(-YDRIP%TDT*YDCONF%NPATFR/ZTAU(JSCALE))
        ZCONSTF0     = SQRT(0.5_JPRB*(1.0_JPRB - ZPHI*ZPHI))*ZSDEV(JSCALE)/SQRT(ZGAMMAN) 
        ZSDEV2D(:,JSCALE) = ZCONSTF0*ZCHI(:)
        ZPHI2D( :,JSCALE) = ZPHI
      ENDDO
      CLLABEL=YDCONF%SM%PN(JPERT)%LABEL
      ISEED= KGET_SEED_SPP(YDCONF, JPERT, JRF, KMEMBER=NENSFNB, KSHIFT=YDCONF%SHIFTSEED,&
        & LABSTIME=YDCONF%LABSTIMSEED, LDVERBOSE=.TRUE.,PTSTEP=YDRIP%TSTEP)
      WRITE(NULOUT,'(''   pattern '',I3,'' for '',A16,'' using seed '',I12)')  IRF,CLLABEL,ISEED 

      YDSPP%IGRIBCODE(IRF)=213100+IRF  ! YDCONF%SM%NACT
      IGRIB(1)=YDSPP%IGRIBCODE(IRF) ! used for grid point
      DO JSCALE=2, JPMXSCALES
        IGRIB(JSCALE)=IGRIB(1)+YDCONF%SM%NRFTOTAL*(JSCALE-1)
        WRITE(NULOUT, *) 'SPP GRIB CODES SCALES ', IGRIB(1), JSCALE, IGRIB(JSCALE)
      ENDDO
      CALL ALLOCATE_ARP(YDGEOMETRY,YDSPP%SP_ARP(IRF),0,JPMXSCALES,0,IGRIB,ISEED,LDCLIP=.TRUE.,LDSUM=.TRUE.)
      IF(YDCONF%LSPG_SPP) THEN
#ifdef FIXME
        CALL ALLOCATE_ARP(YDGEOMETRY,YDSPP%SP_ARP_SPG(IRF,1),0,1,0,IGRIB,ISEED,LDCLIP=.TRUE.,LDSUM=.FALSE.)
        CALL ALLOCATE_ARP(YDGEOMETRY,YDSPP%SP_ARP_SPG(IRF,2),0,1,0,IGRIB,ISEED,LDCLIP=.TRUE.,LDSUM=.FALSE.)
        CALL ALLOCATE_ARP(YDGEOMETRY,YDSPP%SP_ARP_SPG(IRF,3),0,1,0,IGRIB,ISEED,LDCLIP=.TRUE.,LDSUM=.FALSE.)
        CALL SET_ARP2D(YDSPP%SP_ARP(JARP), SPGMU=YDCONF%SPGMU, SPGLAMBDA=YDCONF%SPGLAMBDA, SPGSIGMA=YDCONF%SPGSIGMA,&
       &               SDEVTOT=YDCONF%SDEV, SPGQ=YDCONF%SPGQ, SPGA_DT_MAX=YDCONF%SPGADTMAX, SPGA_DT_MIN=YDCONF%SPGADTMIN,&
       &               SPGTDT=YDCONF%SPGTDT)
        CALL EVOLVE_ARP_SPG(YDSPP%SP_ARP(IRF), &
       &                    YDSPP%SP_ARP_SPG(IRF,1), &
       &                    YDSPP%SP_ARP_SPG(IRF,2), &
       &                    YDSPP%SP_ARP_SPG(IRF,3), &
       &                    LDINIT=.TRUE.)
#endif
    ELSE
      CALL SET_ARP2D(YDSPP%SP_ARP(IRF), ZSDEV2D, ZPHI2D )

      IF (YDCONF%LRDPATINIT) THEN
        WRITE(CLTEMP, "(I10)") YDSPP%IGRIBCODE(IRF)
        WRITE(NULOUT,*) TRIM(YDCONF%SPP_RDPATINIT)//'_'//TRIM(ADJUSTL(CLTEMP))
        CALL READ_SPEC_GRIB(YDMP, TRIM(YDCONF%SPP_RDPATINIT)//'_'//TRIM(ADJUSTL(CLTEMP)), YDSPP%SP_ARP(IRF)%SF, CDLEVTYPE="AL")
      ELSE
        CALL EVOLVE_ARP(YDSPP%SP_ARP(IRF), LDINIT=.TRUE.)
      ENDIF
      ENDIF
      CALL SPECTRAL_NORM(YDSPP%SP_ARP(IRF)%SF, "YDSPP%SP_ARP ["//TRIM(CLLABEL)//"]"//CLJRF)
      CALL ALLOCATE_GRID(YDGEOMETRY,YDSPP%GP_ARP (IRF), 0, 1, IGRIB )
      IF (YDCONF%NPATFR/=1) THEN
        CALL ALLOCATE_GRID(YDGEOMETRY,YDSPP%GP_ARP0(IRF), 0, 1, IGRIB )
        CALL ALLOCATE_GRID(YDGEOMETRY,YDSPP%GP_ARP1(IRF), 0, 1, IGRIB )
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(ZSDEV2D, ZPHI2D, ZCHI)
  !
  !   set up pointers
  !
  DO JPERT=1, YDCONF%SM%NACT
    CALL YDCONF%PPTR%SET( YDCONF%SM%PN(JPERT)%LABEL, JPERT )
  ENDDO
  WRITE( NULOUT, * ) 'YDCONF%PPTR = ', YDCONF%PPTR  
ELSE
  YDSPP%SP_ARP  => NULL()
  !YDSPP%SP_ARP_SPG  => NULL()
  YDSPP%GP_ARP  => NULL()
  YDSPP%GP_ARP0 => NULL()
  YDSPP%GP_ARP1 => NULL()
  YDSPP%LAB_ARP => NULL()
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INI_SPP',1,ZHOOK_HANDLE)

END SUBROUTINE INI_SPP
