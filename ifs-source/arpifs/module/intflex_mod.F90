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

MODULE INTFLEX_MOD

! INTFLEX_MOD
!
!   derived data types and routines for flexible interface
!
!   author: 2013-11, D. Degrauwe
!

USE PARKIND1, ONLY : JPIM, JPRB
USE YOMLUN,   ONLY : NULOUT
USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK, JPHOOK

SAVE

LOGICAL :: LINTFLEX                ! use flexible interface cptend_flex/cputqy instead of cptend_new/cputqy or cputqy_arome
LOGICAL :: LENTHPREC               ! account for heat transport by precipitation in AROME

LOGICAL :: LRADFLEX                ! use flexible interface for AROME radiation separately

! 
! 1. DERIVED DATA TYPES
! 

! 1.1. TYPE_INTFIELD is one field (flux/tendency) that is passed from the physics to the interface
TYPE TYPE_INTFIELD
  CHARACTER(LEN=20)  :: CNAME                     ! name of the field
  CHARACTER(LEN=1)   :: CVAR                      ! 'U', 'V', 'H' (enthalpy) or 'G' (GFL variable)
  INTEGER(KIND=JPIM) :: JGFL                      ! index of GFL variable in YGFL%YCOMP array
  INTEGER(KIND=JPIM) :: JGFLTARGET                ! index of target GFL variable (only for pseudofluxes)
  CHARACTER(LEN=1)   :: CTYPE                     ! for all fields: 'F' (flux), 'T' (tendency)
                                                  ! for water species: 'P' (precipitation flux), 'D' (diffusive flux)
                                                  !                    'S' (pseudoflux), 'C' (corrective flux)
  LOGICAL            :: LDDH                      ! to be included in DDH.
  REAL(KIND=JPRB), POINTER  :: RVAL(:,:)=>NULL()  ! pointer to 2D array
END TYPE TYPE_INTFIELD

! 1.2. TYPE_INTPROC is a set of INTFIELDS
TYPE TYPE_INTPROC
  CHARACTER(LEN=20)            :: CNAME
  INTEGER(KIND=JPIM)           :: NFIELD
  TYPE(TYPE_INTFIELD), POINTER :: RFIELD(:)=>NULL()
END TYPE TYPE_INTPROC

! 1.3. TYPE_INTPROCSET is a set of INTPROCS
TYPE TYPE_INTPROCSET
  INTEGER(KIND=JPIM)          :: NPROCESS
  TYPE(TYPE_INTPROC), POINTER :: RPROCESS(:)=>NULL()
END TYPE TYPE_INTPROCSET

CONTAINS

! 
! 2. CONSTRUCTOR FUNCTIONS
! 

! 2.1. NEWINTFIELD: define a new field
!      NOTE: output is pointer to numerical values!

FUNCTION NEWINTFIELD(YDPROCESS, KPROMA, KFLEV, CDNAME, CDVAR, CDTYPE, KGFL, KGFLTARGET, LDDH)
  
  IMPLICIT NONE
  ! output: pointer to the (numerical values of the) field
  REAL(KIND=JPRB), POINTER :: NEWINTFIELD(:,:)
  ! arguments
  TYPE(TYPE_INTPROC), INTENT(INOUT) :: YDPROCESS    ! process that will contain the field
  INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA          ! horizontal dimension
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV           ! vertical dimension
  CHARACTER(LEN=*), INTENT(IN) :: CDNAME            ! name of the field
  CHARACTER(LEN=1), INTENT(IN) :: CDVAR             ! variable: 'U', 'V', 'H' (enthalpy) or 'G' (GFL variable)
  CHARACTER(LEN=1), INTENT(IN) :: CDTYPE            ! type of field: 'F' (flux), 'T' (tendency), 'S' (pseudoflux), 
                                                    ! 'C' (corrective flux), 'D' (diffusive flux) or 'P' (precipitation flux)
  LOGICAL, INTENT(IN), OPTIONAL :: LDDH             ! field should be passed to ddh
  INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KGFL  ! GFL variable under consideration
  INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KGFLTARGET    ! target GFL variable of pseudoflux

  TYPE(TYPE_INTFIELD), DIMENSION(:), ALLOCATABLE :: YLFIELD   ! auxiliary variable
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:NEWINTFIELD',0,ZHOOK_HANDLE)

  ! 1. Increase process's fields array
  !
  ! NOTE: it is not known in advance how many fields the process will contain, hence this dirty method of reallocation.
  !       Since only metadata are copied and pointers are moved, it should be fast though.
  
  IF (YDPROCESS%NFIELD > 0) THEN
    ! transfer to temporary array
    ALLOCATE(YLFIELD(YDPROCESS%NFIELD))
    YLFIELD(:)=YDPROCESS%RFIELD(:)
    ! reallocate
    DEALLOCATE(YDPROCESS%RFIELD)
    YDPROCESS%NFIELD=YDPROCESS%NFIELD+1
    ALLOCATE(YDPROCESS%RFIELD(YDPROCESS%NFIELD))
    ! transfer back
    YDPROCESS%RFIELD(1:YDPROCESS%NFIELD-1)=YLFIELD(:)
    ! clean up
    DEALLOCATE(YLFIELD)
  ELSE
    YDPROCESS%NFIELD=1
    ALLOCATE(YDPROCESS%RFIELD(YDPROCESS%NFIELD))
  ENDIF

  ! 2. Add field metadata
  YDPROCESS%RFIELD(YDPROCESS%NFIELD)%CNAME = CDNAME
  YDPROCESS%RFIELD(YDPROCESS%NFIELD)%CVAR  = CDVAR
  YDPROCESS%RFIELD(YDPROCESS%NFIELD)%CTYPE = CDTYPE
  YDPROCESS%RFIELD(YDPROCESS%NFIELD)%LDDH = .TRUE. ! default; modified by the next line
  IF (PRESENT(LDDH)) YDPROCESS%RFIELD(YDPROCESS%NFIELD)%LDDH=LDDH
  
  ! 2.1. extra data for GFL variable
  IF (CDVAR=='G') THEN
  
    ! GFL index must be present
    IF (.NOT.PRESENT(KGFL)) THEN
      WRITE (NULOUT,*) 'NEWINTFIELD: trying to add GFL field, but index KGFL not specified'
      CALL ABOR1('NEWINTFIELD: missing KGFL argument')
    ENDIF
    YDPROCESS%RFIELD(YDPROCESS%NFIELD)%JGFL=KGFL
  
    ! for pseudoflux, target GFL index must be present
    IF (CDTYPE=='S') THEN
      IF (.NOT. PRESENT(KGFLTARGET)) THEN
        WRITE (NULOUT,*) 'NEWINTFIELD: trying to add pseudoflux, but index KGFLTARGET not specified'
        CALL ABOR1('NEWINTFIELD: missing KGFLTARGET argument')
      ENDIF
      YDPROCESS%RFIELD(YDPROCESS%NFIELD)%JGFLTARGET=KGFLTARGET
    ENDIF

  ENDIF

  ! 3. allocate and set pointer to data
  NULLIFY(YDPROCESS%RFIELD(YDPROCESS%NFIELD)%RVAL)
  IF (CDTYPE=='T') THEN
    ALLOCATE(YDPROCESS%RFIELD(YDPROCESS%NFIELD)%RVAL(KPROMA,KFLEV))
  ELSE
    ALLOCATE(YDPROCESS%RFIELD(YDPROCESS%NFIELD)%RVAL(KPROMA,0:KFLEV))
  ENDIF
  NEWINTFIELD=>YDPROCESS%RFIELD(YDPROCESS%NFIELD)%RVAL

  ! initialize to zeros ???
  !NEWINTFIELD=0.0_JPRB
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:NEWINTFIELD',1,ZHOOK_HANDLE)

END FUNCTION NEWINTFIELD

! 2.2. NEWINTPROC : define a new process
!      NOTE: output is pointer to process structure inside procset

FUNCTION NEWINTPROC(YDPROCSET,CDNAME)

  IMPLICIT NONE

  ! result
  TYPE(TYPE_INTPROC),    POINTER       :: NEWINTPROC    ! output is pointer to the defined process

  ! arguments
  TYPE(TYPE_INTPROCSET), INTENT(INOUT) :: YDPROCSET     ! set of processes to which this process is added
  CHARACTER(LEN=*),      INTENT(IN)    :: CDNAME        ! name of the process
  
  ! local variables
  TYPE(TYPE_INTPROC), ALLOCATABLE :: YLPROC(:)
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:NEWINTPROC',0,ZHOOK_HANDLE)
     
  ! 1. Increase procset's process array
  !
  ! NOTE: it is not known in advance how many processes the procset will contain, hence this dirty method of reallocation.
  !       Since only metadata are copied and pointers are moved, it should be fast though.

  IF (YDPROCSET%NPROCESS > 0) THEN
    ! transfer to temporary array
    ALLOCATE(YLPROC(YDPROCSET%NPROCESS))
    YLPROC(:)=YDPROCSET%RPROCESS(:)
    ! reallocate
    DEALLOCATE(YDPROCSET%RPROCESS)
    YDPROCSET%NPROCESS=YDPROCSET%NPROCESS+1
    ALLOCATE(YDPROCSET%RPROCESS(YDPROCSET%NPROCESS))
    ! transfer back
    YDPROCSET%RPROCESS(1:YDPROCSET%NPROCESS-1)=YLPROC(:)
    ! clean up
    DEALLOCATE(YLPROC)
  ELSE
    YDPROCSET%NPROCESS=1
    NULLIFY(YDPROCSET%RPROCESS)
    ALLOCATE(YDPROCSET%RPROCESS(YDPROCSET%NPROCESS))
  ENDIF

  ! 2. set metadata for new process
  YDPROCSET%RPROCESS(YDPROCSET%NPROCESS)%CNAME=CDNAME
  YDPROCSET%RPROCESS(YDPROCSET%NPROCESS)%NFIELD=0
  NULLIFY(YDPROCSET%RPROCESS(YDPROCSET%NPROCESS)%RFIELD)
  ALLOCATE(YDPROCSET%RPROCESS(YDPROCSET%NPROCESS)%RFIELD(0))
  
  NEWINTPROC=>YDPROCSET%RPROCESS(YDPROCSET%NPROCESS)
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:NEWINTPROC',1,ZHOOK_HANDLE)

END FUNCTION NEWINTPROC

! 2.3. NEWINTPROCSET : setup a set of processes
!
FUNCTION NEWINTPROCSET()
  
  IMPLICIT NONE
  
  ! result
  TYPE(TYPE_INTPROCSET) :: NEWINTPROCSET
  
  ! auxiliary variables
  REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:NEWINTPROCSET',0,ZHOOK_HANDLE)
  
  ! set number of processes to zero
  NEWINTPROCSET%NPROCESS=0
  ! set pointer to zero-sized array
  NULLIFY(NEWINTPROCSET%RPROCESS)
  ALLOCATE(NEWINTPROCSET%RPROCESS(0))
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:NEWINTPROCSET',1,ZHOOK_HANDLE)

END FUNCTION NEWINTPROCSET

! 
! 3. DESTRUCTOR SUBROUTINES
! 

! 3.1. CLEANINTFIELD: clean a field structure
! 
SUBROUTINE CLEANINTFIELD(YDFIELD)

  IMPLICIT NONE
  
  ! arguments
  TYPE(TYPE_INTFIELD), INTENT(INOUT) :: YDFIELD
  
  ! auxiliary variables
  REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:CLEANINTFIELD',0,ZHOOK_HANDLE)

  ! deallocate the numerical array
  DEALLOCATE(YDFIELD%RVAL)
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:CLEANINTFIELD',1,ZHOOK_HANDLE)

END SUBROUTINE CLEANINTFIELD

! 3.1. CLEANINTPROC: clean a process structure
! 
SUBROUTINE CLEANINTPROC(YDPROC)

  IMPLICIT NONE
  
  ! arguments
  TYPE(TYPE_INTPROC), INTENT(INOUT) :: YDPROC
  
  ! auxiliary variables
  INTEGER(KIND=JPIM) :: JFIELD
  REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:CLEANINTPROC',0,ZHOOK_HANDLE)
  
  ! clean up contained fields
  DO JFIELD=1,YDPROC%NFIELD
    CALL CLEANINTFIELD(YDPROC%RFIELD(JFIELD))
  ENDDO
  
  ! deallocate the array of pointers
  DEALLOCATE(YDPROC%RFIELD)
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:CLEANINTPROC',1,ZHOOK_HANDLE)
  
END SUBROUTINE CLEANINTPROC


! 3.3. CLEANINTPROCSET: clean up a procset
!
SUBROUTINE CLEANINTPROCSET(YDPROCSET)

  IMPLICIT NONE

  ! arguments
  TYPE(TYPE_INTPROCSET), INTENT(INOUT) :: YDPROCSET  ! set of processes to clean
  
  ! auxiliary variables
  INTEGER(KIND=JPIM) :: JPROCESS
  REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
  
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:CLEANINTPROCSET',0,ZHOOK_HANDLE)

  ! clean up contained processes
  DO JPROCESS=1,YDPROCSET%NPROCESS
    CALL CLEANINTPROC(YDPROCSET%RPROCESS(JPROCESS))
  ENDDO
  
  ! deallocate the array of pointers
  DEALLOCATE(YDPROCSET%RPROCESS)
 
  IF (LHOOK) CALL DR_HOOK('INTFLEX_MOD:CLEANINTPROCSET',1,ZHOOK_HANDLE)
  
END SUBROUTINE CLEANINTPROCSET

END MODULE INTFLEX_MOD
