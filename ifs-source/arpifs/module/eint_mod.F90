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

MODULE EINT_MOD

!--------------------------------------------------------------------------
! EINT - Externalisable part of horizonal interpolators and halo management
!        used for example in the semi-Lagrangian scheme
!
!   Modifications.
!   --------------
!     M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!--------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC SL_STRUCT, DEALLO_SL_STRUCT, UNUSED_HALO_STATS, SLHALO_DEBUG, &
       & JPDUP

TYPE SL_STRUCT

  ! CVER: two letter code describing the SL_STRUCT version
  !             'SL': SL advection
  !             'FP': Fullpos
  !             'OB': Observation
  !             'OA': Observation (adjoint)
  !             'RI': Radiation input
  !             'RO': Radiation output
  CHARACTER(LEN=2) :: CVER='XX'

  ! NSLGROUP: flag to allow different coding
  INTEGER(KIND=JPIM) :: NSLGROUP=0

  ! NSLSTA: interpolation buffer start position of grid columns.
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLSTA(:)   

  ! NSLONL: interpolation buffer number of grid column on latitudes.
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLONL(:)   

  ! NSLOFF: interpolation buffer offset to start of each row in.
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLOFF(:)   

  ! NSLPTSWEST: number of grid points on western side of core points
  !             only for LCOMPLAT=F
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLPTSWEST(:)   

  ! NSLPTSEAST: number of grid points on esatern side of core points
  !             only for LCOMPLAT=F
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLPTSEAST(:)   

  ! NSLEXT: pointer that makes sure addressing of points in the east-west
  !         extension zone is correct. It also handles the half latitude 
  !         shift of extension latitudes at the poles.
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLEXT(:,:) 

  ! LSLWHOLELAT: records whether complete latitude is present
  LOGICAL,ALLOCATABLE :: LCOMPLAT(:) 

  ! NLATGLO: global latitude
  INTEGER(KIND=JPIM),ALLOCATABLE :: NLATGLO(:) 

  ! DIST1GP: distance in kilometres of one grid point
  REAL(KIND=JPRB),ALLOCATABLE :: DIST1GP(:) 

  ! NSLSENDPOS: the addresses within the interpolation buffer of point
  !             sent from this PE (computed only if NPROC>1).
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLSENDPOS(:)

  ! NSLRECVPOS: the addresses within the interpolation buffer of point
  !             received on this PE (computed only if NPROC>1).
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLRECVPOS(:)

  ! NSLSENDPTR: allocatable to the first point for each of the PE's that has
  !             to receive interpolation halo-data from this.
  !             Used for addressing NSLSENDPOS(). 
  !             Computed only if NPROC > 1.
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLSENDPTR(:)

  ! NSLRECVPTR: pointer to the first point for each of the PE's that are
  !             sending interpolation halo-data to this PE.
  !             Used for addressing NSLRECVPOS(). 
  !             Computed only if NPROC > 1.
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLRECVPTR(:)


  ! NSLCOMM: list of the processors this processor has to communicate
  !          with (computed only if NPROC>1).
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLCOMM(:)  

  ! LSLCOMM: flag indicating whether this processor needs to communicate 
  !          with a specific processor
  LOGICAL           ,ALLOCATABLE :: LSLCOMM(:,:)

  ! NASLB1: local inner dimension of interpolation buffer.
  INTEGER(KIND=JPIM) :: NASLB1=0

  ! NASLB1_TRUE: local inner dimension of interpolation buffer,
  !              without padding for reducing memory conflicts
  INTEGER(KIND=JPIM) :: NASLB1_TRUE=0

  ! NSLPAD: number of pad words initialised to a huge number at either
  !         of side of the sl halo, used to trap halo problems.
  INTEGER(KIND=JPIM) :: NSLPAD=0

  ! LSLT_ARRAYS_INIT : true if on-demand per timestep masks are allocated
  LOGICAL            :: LSLT_ARRAYS_INIT=.FALSE.

  ! LSLONDEM: on-demand SL communications switch
  LOGICAL            :: LSLONDEM=.FALSE.

  ! LSLONDEM_ACTIVE: on-demand SL communications active switch
  LOGICAL            :: LSLONDEM_ACTIVE=.FALSE.

  ! NSAFEHALO: minimum of unused halo grid points on a latitude (West or East)
  ! over all latitudes of the halo
  INTEGER(KIND=JPIM) :: NUNUSEDHALO = 999999

  ! DISTSAFEHALO: minimum distance in km of unused halo grid points on a latitude
  ! (West or East) over all latitudes of the halo
  REAL(KIND=JPRB) :: DISTUNUSEDHALO = 999999.0_JPRB

  ! MASKS for on-demand SL comms
  INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_SL1(:)   
  INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_SL2(:)   
  INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_SL2T(:,:)

  ! MASK for SL debugging (LSLDEBUG=T)
  INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_SLD(:)   

  ! NSLPROCS: number of processors to communicate with.
  INTEGER(KIND=JPIM) :: NSLPROCS=0

  ! NSLRPT: dimension for NSLRECVPOS (the number of columns received
  !         from other PE's when computing the halo for interpolations).
  INTEGER(KIND=JPIM) :: NSLRPT=0

  ! NSLSPT: dimension for NSLSENDPOS (the number of columns sent
  !         to other PE's when computing the halo for interpolations).
  INTEGER(KIND=JPIM) :: NSLSPT=0

  ! NSLWIDEN,NSLWIDES,NSLWIDEE,NSLWIDEW: number of grid points required for halo
  !                                      (resp. north, south, east, west).
  INTEGER(KIND=JPIM) :: NSLWIDEN=0
  INTEGER(KIND=JPIM) :: NSLWIDES=0
  INTEGER(KIND=JPIM) :: NSLWIDEE=0
  INTEGER(KIND=JPIM) :: NSLWIDEW=0

  ! NSLWIDE number of grid points required for halo
  INTEGER(KIND=JPIM) :: NSLWIDE=0

  ! NMAP: temporary data structure only for SL adjoint reproducibility
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLMAP(:,:)     

  ! NSLCORE: pointer to this processors core region points within the
  !          interpolation buffer.
  INTEGER(KIND=JPIM),ALLOCATABLE :: NSLCORE(:)  

  ! LSLCORE: T if NSLCORE point, F if not
  LOGICAL,ALLOCATABLE            :: LSLCORE(:)  

  ! MASK_SLTOT: sum of sl buffer points used in trajectory calc.
  !             used to optimise ngptotad and thereby improve performance
  INTEGER(KIND=JPIM),ALLOCATABLE :: MASK_SLTOT(:)

  ! Copy of some variables to simplify the call tree for the SL interface
  INTEGER(KIND=JPIM) :: NDGLG=0
  INTEGER(KIND=JPIM) :: NDLON=0
  INTEGER(KIND=JPIM) :: NDGSAG=0
  INTEGER(KIND=JPIM) :: NDGENG=0
  INTEGER(KIND=JPIM) :: NDGSAL=0
  INTEGER(KIND=JPIM) :: NDGENL=0
  INTEGER(KIND=JPIM) :: NDGSAH=0
  INTEGER(KIND=JPIM) :: NDGENH=0
  INTEGER(KIND=JPIM) :: NGPTOT=0
  INTEGER(KIND=JPIM) :: NDGUXL=0
  INTEGER(KIND=JPIM) :: NDLUNG=0
  INTEGER(KIND=JPIM) :: NDLUXG=0
  INTEGER(KIND=JPIM) :: NDGUNG=0
  INTEGER(KIND=JPIM) :: NDGUXG=0
  INTEGER(KIND=JPIM) :: NDSUR1=0
  INTEGER(KIND=JPIM) :: NDLSUR=0
  INTEGER(KIND=JPIM) :: NDGSUR=0
  INTEGER(KIND=JPIM) :: NPTRFLOFF=0
  INTEGER(KIND=JPIM) :: NFRSTLOFF=0
  INTEGER(KIND=JPIM) :: MYFRSTACTLAT=0
  INTEGER(KIND=JPIM) :: MYLSTACTLAT=0
  INTEGER(KIND=JPIM),ALLOCATABLE :: NLOENG(:)

CONTAINS

  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION

END TYPE SL_STRUCT


! structures used to compute and fill halo for interpolations
! YRSL : semi-lagrangian advection
!! TYPE(SL_STRUCT),POINTER :: YRSL !! moved to MODEL_DYNAMICS_TYPE
! YRAD : semi-lagrangian advection (adjoint)
!! TYPE(SL_STRUCT),POINTER :: YRAD !! moved to MODEL_DYNAMICS_TYPE
! YRRI : model grid to radiation grid
!! TYPE(SL_STRUCT),POINTER :: YRRI !! moved to MODEL_PHYSICS_RADIATION_TYPE
! YRRO : radiation grid to model grid
!! TYPE(SL_STRUCT),POINTER :: YRRO !! moved to MODEL_PHYSICS_RADIATION_TYPE

! -----------------------------------------------------------------------------

! * Duplication of some dummy or local arrays for optimisation on NEC platform.
#ifdef NECSX
!Note: To make this way of preventing memory bank conflict on vector machines to work again
!      some coding effort is required for LASCAW/TL/AD and routines bellow.
!      Possibly this problem is no longer an issue for current computers.
INTEGER(KIND=JPIM), PARAMETER :: JPDUP=33
#else
INTEGER(KIND=JPIM), PARAMETER :: JPDUP=1
#endif
 
! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO, CNAME)
  IMPLICIT NONE
  CLASS(SL_STRUCT), INTENT(IN) :: SELF
  INTEGER         , INTENT(IN) :: KDEPTH
  INTEGER         , INTENT(IN) :: KOUTNO
  CHARACTER(LEN=*), INTENT(IN) :: CNAME

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH+2

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_gconf%' // CNAME // ' : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CVER = ', SELF%CVER
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLGROUP = ', SELF%NSLGROUP
  IF (ALLOCATED(SELF%NSLSTA)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLSTA allocated of shape ', SHAPE(SELF%NSLSTA), &
 &        ' and sum ',SUM(SELF%NSLSTA)
  IF (ALLOCATED(SELF%NSLONL)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLONL allocated of shape ', SHAPE(SELF%NSLONL), &
 &        ' and sum ',SUM(SELF%NSLONL)
  IF (ALLOCATED(SELF%NSLOFF)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLOFF allocated of shape ', SHAPE(SELF%NSLOFF), &
 &        ' and sum ',SUM(SELF%NSLOFF)
  IF (ALLOCATED(SELF%NSLPTSWEST)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLPTSWEST allocated of shape ', SHAPE(SELF%NSLPTSWEST), &
 &        ' and sum ',SUM(SELF%NSLPTSWEST)
  IF (ALLOCATED(SELF%NSLPTSEAST)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLPTSEAST allocated of shape ', SHAPE(SELF%NSLPTSEAST), &
 &        ' and sum ',SUM(SELF%NSLPTSEAST)
  IF (ALLOCATED(SELF%NSLEXT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLEXT allocated of shape ', SHAPE(SELF%NSLEXT), &
 &        ' and sum ',SUM(SELF%NSLEXT)
  IF (ALLOCATED(SELF%LCOMPLAT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCOMPLAT allocated of shape ', SHAPE(SELF%LCOMPLAT), &
 &        ' values ',SELF%LCOMPLAT
  IF (ALLOCATED(SELF%NLATGLO)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLATGLO allocated of shape ', SHAPE(SELF%NLATGLO), &
 &        ' and sum ',SUM(SELF%NLATGLO)
  IF (ALLOCATED(SELF%DIST1GP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'DIST1GP allocated of shape ', SHAPE(SELF%DIST1GP), &
 &        ' and sum ',SUM(SELF%DIST1GP)
  IF (ALLOCATED(SELF%NSLSENDPOS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLSENDPOS allocated of shape ', SHAPE(SELF%NSLSENDPOS), &
 &        ' and sum ',SUM(SELF%NSLSENDPOS)
  IF (ALLOCATED(SELF%NSLRECVPOS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLRECVPOS allocated of shape ', SHAPE(SELF%NSLRECVPOS), &
 &        ' and sum ',SUM(SELF%NSLRECVPOS)
  IF (ALLOCATED(SELF%NSLSENDPTR)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLSENDPTR allocated of shape ', SHAPE(SELF%NSLSENDPTR), &
 &        ' and sum ',SUM(SELF%NSLSENDPTR)
  IF (ALLOCATED(SELF%NSLRECVPTR)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLRECVPTR allocated of shape ', SHAPE(SELF%NSLRECVPTR), &
 &        ' and sum ',SUM(SELF%NSLRECVPTR)
  IF (ALLOCATED(SELF%NSLCOMM)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLCOMM allocated of shape ', SHAPE(SELF%NSLCOMM), &
 &        ' and sum ',SUM(SELF%NSLCOMM)
  IF (ALLOCATED(SELF%LSLCOMM)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSLCOMM allocated of shape ', SHAPE(SELF%LSLCOMM), &
 &        ' values ',SELF%LSLCOMM

  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NASLB1 = ', SELF%NASLB1
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NASLB1_TRUE = ', SELF%NASLB1_TRUE
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLPAD = ', SELF%NSLPAD
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSLT_ARRAYS_INIT = ', SELF%LSLT_ARRAYS_INIT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSLONDEM = ', SELF%LSLONDEM
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSLONDEM_ACTIVE = ', SELF%LSLONDEM_ACTIVE
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NUNUSEDHALO = ', SELF%NUNUSEDHALO
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'DISTUNUSEDHALO = ', SELF%DISTUNUSEDHALO
  IF (ALLOCATED(SELF%MASK_SL1)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MASK_SL1 allocated of shape ', SHAPE(SELF%MASK_SL1), &
 &        ' and sum ',SUM(SELF%MASK_SL1)
  IF (ALLOCATED(SELF%MASK_SL2)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MASK_SL2 allocated of shape ', SHAPE(SELF%MASK_SL2), &
 &        ' and sum ',SUM(SELF%MASK_SL2)
  IF (ALLOCATED(SELF%MASK_SL2T)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MASK_SL2T allocated of shape ', SHAPE(SELF%MASK_SL2T), &
 &        ' and sum ',SUM(SELF%MASK_SL2T)
  IF (ALLOCATED(SELF%MASK_SLD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MASK_SLD allocated of shape ', SHAPE(SELF%MASK_SLD), &
 &        ' and sum ',SUM(SELF%MASK_SLD)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLPROCS = ', SELF%NSLPROCS
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLRPT = ', SELF%NSLRPT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLSPT = ', SELF%NSLSPT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLWIDEN = ', SELF%NSLWIDEN
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLWIDES = ', SELF%NSLWIDES
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLWIDEE = ', SELF%NSLWIDEE
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLWIDEW = ', SELF%NSLWIDEW
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLWIDE = ', SELF%NSLWIDE
  IF (ALLOCATED(SELF%NSLMAP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLMAP allocated of shape ', SHAPE(SELF%NSLMAP), &
 &        ' and sum ',SUM(SELF%NSLMAP)
  IF (ALLOCATED(SELF%NSLCORE)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLCORE allocated of shape ', SHAPE(SELF%NSLCORE), &
 &        ' and sum ',SUM(SELF%NSLCORE)
  IF (ALLOCATED(SELF%LSLCORE)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSLCORE allocated of shape ', SHAPE(SELF%LSLCORE), &
 &        ' values ',SELF%LSLCORE
  IF (ALLOCATED(SELF%MASK_SLTOT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MASK_SLTOT allocated of shape ', SHAPE(SELF%MASK_SLTOT), &
 &        ' and sum ',SUM(SELF%MASK_SLTOT)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGLG = ', SELF%NDGLG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDLON = ', SELF%NDLON
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSAG = ', SELF%NDGSAG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGENG = ', SELF%NDGENG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSAL = ', SELF%NDGSAL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGENL = ', SELF%NDGENL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSAH = ', SELF%NDGSAH
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGENH = ', SELF%NDGENH
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGPTOT = ', SELF%NGPTOT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGUXL = ', SELF%NDGUXL
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDLUNG = ', SELF%NDLUNG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDLUXG = ', SELF%NDLUXG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGUNG = ', SELF%NDGUNG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGUXG = ', SELF%NDGUXG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDSUR1 = ', SELF%NDSUR1
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDLSUR = ', SELF%NDLSUR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDGSUR = ', SELF%NDGSUR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPTRFLOFF = ', SELF%NPTRFLOFF
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFRSTLOFF = ', SELF%NFRSTLOFF
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MYFRSTACTLAT = ', SELF%MYFRSTACTLAT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MYLSTACTLAT = ', SELF%MYLSTACTLAT
  IF (ALLOCATED(SELF%NLOENG)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLOENG allocated of shape ', SHAPE(SELF%NLOENG), &
 &        ' and sum ',SUM(SELF%NLOENG)

 
END SUBROUTINE PRINT_CONFIGURATION

! ------------------------------------------------------------------

SUBROUTINE DEALLO_SL_STRUCT(YDSL)

  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

  TYPE(SL_STRUCT), INTENT(INOUT) :: YDSL
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('EINT_MOD:DEALLO_SL_STRUCT',0,ZHOOK_HANDLE)
  IF(ALLOCATED(YDSL%NSLSTA))     DEALLOCATE(YDSL%NSLSTA    )
  IF(ALLOCATED(YDSL%NSLONL))     DEALLOCATE(YDSL%NSLONL    )
  IF(ALLOCATED(YDSL%NSLOFF))     DEALLOCATE(YDSL%NSLOFF    )
  IF(ALLOCATED(YDSL%NSLPTSWEST)) DEALLOCATE(YDSL%NSLPTSWEST)
  IF(ALLOCATED(YDSL%NSLPTSEAST)) DEALLOCATE(YDSL%NSLPTSEAST)
  IF(ALLOCATED(YDSL%NSLEXT))     DEALLOCATE(YDSL%NSLEXT    )
  IF(ALLOCATED(YDSL%LCOMPLAT))   DEALLOCATE(YDSL%LCOMPLAT  )
  IF(ALLOCATED(YDSL%NLATGLO))    DEALLOCATE(YDSL%NLATGLO   )
  IF(ALLOCATED(YDSL%DIST1GP))    DEALLOCATE(YDSL%DIST1GP   )
  IF(ALLOCATED(YDSL%NSLSENDPOS)) DEALLOCATE(YDSL%NSLSENDPOS)
  IF(ALLOCATED(YDSL%NSLRECVPOS)) DEALLOCATE(YDSL%NSLRECVPOS)
  IF(ALLOCATED(YDSL%NSLSENDPTR)) DEALLOCATE(YDSL%NSLSENDPTR)
  IF(ALLOCATED(YDSL%NSLRECVPTR)) DEALLOCATE(YDSL%NSLRECVPTR)
  IF(ALLOCATED(YDSL%NSLCORE))    DEALLOCATE(YDSL%NSLCORE   )
  IF(ALLOCATED(YDSL%LSLCORE))    DEALLOCATE(YDSL%LSLCORE   )
  IF(ALLOCATED(YDSL%NSLCOMM))    DEALLOCATE(YDSL%NSLCOMM   )
  IF(ALLOCATED(YDSL%LSLCOMM))    DEALLOCATE(YDSL%LSLCOMM   )
  IF(ALLOCATED(YDSL%MASK_SL1))   DEALLOCATE(YDSL%MASK_SL1  )
  IF(ALLOCATED(YDSL%MASK_SL2))   DEALLOCATE(YDSL%MASK_SL2  )
  IF(ALLOCATED(YDSL%MASK_SL2T))  DEALLOCATE(YDSL%MASK_SL2T )
  IF(ALLOCATED(YDSL%MASK_SLD))   DEALLOCATE(YDSL%MASK_SLD  )
  IF(ALLOCATED(YDSL%NSLMAP))    DEALLOCATE(YDSL%NSLMAP    )
  IF(ALLOCATED(YDSL%NLOENG))     DEALLOCATE(YDSL%NLOENG    )
  IF (LHOOK) CALL DR_HOOK('EINT_MOD:DEALLO_SL_STRUCT',1,ZHOOK_HANDLE)
END SUBROUTINE DEALLO_SL_STRUCT

! ------------------------------------------------------------------

SUBROUTINE UNUSED_HALO_STATS(YDSL,YDAD,YDRI,YDRO)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULOUT

TYPE(SL_STRUCT), INTENT(INOUT) :: YDSL
TYPE(SL_STRUCT), INTENT(INOUT) :: YDAD
TYPE(SL_STRUCT), INTENT(INOUT) :: YDRI
TYPE(SL_STRUCT), INTENT(INOUT) :: YDRO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EINT_MOD:UNUSED_HALO_STATS',0,ZHOOK_HANDLE)
WRITE(NULOUT,'("MINIMUM UNUSED HALO GRID POINTS AND MINIMUM UNUSED HALO DISTANCE ",&
  &"PER LATITUDE OVER ALL HALO LATITUDES")')
CALL UNUSED_HALO_STATS1(YDSL)
CALL UNUSED_HALO_STATS1(YDAD)
CALL UNUSED_HALO_STATS1(YDRI)
CALL UNUSED_HALO_STATS1(YDRO)
!CALL UNUSED_HALO_STATS1(YROB)
!CALL UNUSED_HALO_STATS1(YROA)
IF (LHOOK) CALL DR_HOOK('EINT_MOD:UNUSED_HALO_STATS',1,ZHOOK_HANDLE)
END SUBROUTINE UNUSED_HALO_STATS

! ------------------------------------------------------------------

SUBROUTINE UNUSED_HALO_STATS1(YDSL)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULOUT

TYPE(SL_STRUCT), INTENT(IN) :: YDSL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EINT_MOD:UNUSED_HALO_STATS1',0,ZHOOK_HANDLE)
IF( YDSL%NUNUSEDHALO /= 999999 )THEN
  WRITE(NULOUT,'("EINT_MOD: ",A," GRID POINTS=",I3," DISTANCE IN KM=",F8.1)')&
   & YDSL%CVER,YDSL%NUNUSEDHALO,YDSL%DISTUNUSEDHALO
ENDIF
IF (LHOOK) CALL DR_HOOK('EINT_MOD:UNUSED_HALO_STATS1',1,ZHOOK_HANDLE)
END SUBROUTINE UNUSED_HALO_STATS1

SUBROUTINE SLHALO_DEBUG(YDSL,KIND)
USE YOMLUN   , ONLY : NULERR
USE YOMMP0   , ONLY : MYPROC
IMPLICIT NONE
TYPE(SL_STRUCT),   INTENT(IN) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN) ::KIND
#include "abor1.intfb.h"
IF( YDSL%MASK_SLD(KIND) /= 888888 .AND. YDSL%MASK_SLD(KIND) /= 999999 )THEN
  WRITE(NULERR,'("SLHALO_DEBUG: HALO INDEX INVALID, MYPROC=",I6," KIND=",I10)')MYPROC,KIND
  CALL ABOR1('SLHALO_DEBUG: HALO INDEX INVALID')
ENDIF
END SUBROUTINE SLHALO_DEBUG

! ------------------------------------------------------------------

END MODULE EINT_MOD
