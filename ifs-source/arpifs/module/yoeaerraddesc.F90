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

MODULE YOEAERRADDESC

! YOEAERRADDESC - Flexible aerosol-radiation description
!
! PURPOSE
! -------
!   Store information about which aerosol species are to be used in
!   the radiation scheme, which are prognostic and which diagnostic
!   (climatological), and information to map from the relevant GFL
!   fields or variables in the climatology file to the optical
!   properties in ecRad.
!
! INTERFACE
! ---------
!   All information is stored in a TAER_RAD_DESC object. You can use
!   the READ_NAMELIST member routine to read its contents from
!   namelist, which contains multiple groups of the same name, e.g.:
!
! &NAM_AEROSOL_RADIATION
!    LPROGNOSTIC = true, ! This aerosol type will be prognostic
!    NAME = "Dust1",     ! GFL field to search for; abort if not found
!    LHYDROPHILIC = false,
!    SCALING = 1.0,      ! Optional scaling of global aerosol load
!    AER_TYPE = "DU",    ! Two character code for type
!    OPTICS_MODEL = "Dubovic", ! ...and the specific optical model
!    SIZE_BIN = 1,       ! Size bin for optical properties
! /
! &NAM_AEROSOL_RADIATION
!    LPROGNOSTIC = false, ! This aerosol type will be climatological
!    NAME = "Sulphates",  ! Variable name in climatology file
!    LHYDROPHILIC = true,
!    SCALING = 1.0,
!    AER_TYPE = "SS",
!    OPTICS_MODEL = "OPAC",
! /
!
!   Alternatively, call the ALLOCATE member routine followed by
!   repeated calls of ADD_PROGNOSTIC, ADD_CLIMATOLOGICAL and/or
!   ADD_PARAMETRIC to build up the list of aerosol species required,
!   and their properties.
!
! AUTHOR
! ------
!    Robin Hogan, ECMWF
!    Original: 2022-01-24
!
! MODIFICATIONS
! -------------
!
! ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM), PARAMETER, PRIVATE :: JPMAX_STRING_LEN = 512

! Aerosols can come from a climatology (read from file), be a
! prognostic variable, or be parameterized (e.g. for volcanoes or
! background aerosols)
ENUM, BIND(C)
  ENUMERATOR ICLIMATOLOGICAL, IPROGNOSTIC, IPARAMETRIC
END ENUM

! ------------------------------------------------------------------
! Structure for storing a description of how a single aerosol species
! will be prescribed and how it interacts with radiation
TYPE :: TAER_RAD_DESC1
  CHARACTER(LEN=JPMAX_STRING_LEN) :: NAME
  INTEGER(KIND=JPIM)              :: NSOURCE      = ICLIMATOLOGICAL
  LOGICAL                         :: LHYDROPHILIC = .FALSE.
  REAL(KIND=JPRB)                 :: SCALING      = 1.0_JPRB
  CHARACTER(2)                    :: AER_TYPE     = '  '
  CHARACTER(LEN=JPMAX_STRING_LEN) :: OPTICS_MODEL = ' '
  INTEGER(KIND=JPIM)              :: SIZE_BIN     = 0 ! 0 = not binned
  ! PARAMETRIC_SCHEME is only used if NSOURCE==IPARAMETRIC: -1 (not
  ! set, error), 0 (tropospheric background), 1 (stratospheric
  ! background), 2 (volcanic)
  INTEGER(KIND=JPIM)              :: PARAMETRIC_SCHEME = -1
END TYPE TAER_RAD_DESC1

! ------------------------------------------------------------------
! Structure for describing how all aerosol species are prescribed and
! how they interact with radiation
TYPE :: TAER_RAD_DESC
  TYPE(TAER_RAD_DESC1), DIMENSION(:), ALLOCATABLE :: DATA
  INTEGER(KIND=JPIM) :: NAEROSOL    ! Number of species described
  INTEGER(KIND=JPIM) :: NDIAGNOSTIC ! Number of diagnostic species (=climatological+parametric)
  INTEGER(KIND=JPIM) :: NPROGNOSTIC ! Number of prognostic species
  INTEGER(KIND=JPIM) :: NPARAMETRIC ! Number of parametric species
  INTEGER(KIND=JPIM) :: NCLIMATOLOGICAL ! Number of climatological species
  ! Index to GFL aerosol field, of length NPROGNOSTIC
  INTEGER(KIND=JPIM), ALLOCATABLE :: PROG_INDEX(:)
CONTAINS

  PROCEDURE :: READ_NAMELIST
  PROCEDURE :: RESERVE
  PROCEDURE :: CLEAR
  PROCEDURE :: PRINT
  PROCEDURE :: PRINT_NAMELIST
  PROCEDURE :: ADD_CLIMATOLOGICAL
  PROCEDURE :: ADD_PARAMETRIC
  PROCEDURE :: ADD_PROGNOSTIC
  PROCEDURE :: ADD_ALL_PROGNOSTICS
  PROCEDURE :: FIND_PROGNOSTICS
  PROCEDURE :: SET_HYDROPHOBIC
  PROCEDURE :: PARAMETRIC_INDEX
  PROCEDURE :: GET_INDEX

END TYPE TAER_RAD_DESC

CONTAINS

! ------------------------------------------------------------------
! Read an aerosol-radiation description from a namelist with unit
! KUNIT, looking for group-name NAM_AEROSOL_RADIATION.  Multiple
! entries are possible, each interpreted as a new aerosol species
SUBROUTINE READ_NAMELIST(SELF, KUNIT)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,  ONLY : NULERR
  
  CLASS(TAER_RAD_DESC), INTENT(INOUT) :: SELF
  INTEGER(KIND=JPIM),   INTENT(IN)    :: KUNIT

  CHARACTER(LEN=JPMAX_STRING_LEN) :: NAME
  LOGICAL                         :: LPROGNOSTIC
  LOGICAL                         :: LPARAMETRIC
  LOGICAL                         :: LHYDROPHILIC
  REAL(KIND=JPRB)                 :: SCALING
  CHARACTER(LEN=JPMAX_STRING_LEN) :: OPTICS_MODEL
  CHARACTER(2)                    :: AER_TYPE
  INTEGER(KIND=JPIM)              :: SIZE_BIN
  INTEGER(KIND=JPIM)              :: PARAMETRIC_SCHEME

  INTEGER(KIND=JPIM)              :: ISTAT
  CHARACTER(LEN=JPMAX_STRING_LEN) :: DUMMY_LINE
  INTEGER(KIND=JPIM)              :: NAER ! Number of aerosol species
  INTEGER(KIND=JPIM)              :: JAER ! Loop counter

  REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE

  NAMELIST /NAM_AEROSOL_RADIATION/ LPROGNOSTIC, LPARAMETRIC, NAME, &
       &   SCALING, AER_TYPE, OPTICS_MODEL, SIZE_BIN, LHYDROPHILIC, &
       &   PARAMETRIC_SCHEME

#include "posname.intfb.h"
#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:READ_NAMELIST',0,ZHOOK_HANDLE)

  CALL SELF%CLEAR()
  NAER = 0

  ! Count the number of aerosol species
  REWIND(KUNIT)
  DO 
    CALL POSNAME(KUNIT, 'NAM_AEROSOL_RADIATION', ISTAT, LDNOREWIND=.TRUE.)
    IF (ISTAT == 1) EXIT ! No more entries found
    NAER = NAER + 1
    READ(KUNIT, '(A)', IOSTAT=ISTAT) DUMMY_LINE
  ENDDO

  IF (NAER > 0) THEN
    ! At least one aerosol species found
    ALLOCATE(SELF%DATA(NAER))
    REWIND(KUNIT)
    DO JAER = 1,NAER
      CALL POSNAME(KUNIT, 'NAM_AEROSOL_RADIATION', ISTAT, LDNOREWIND=.TRUE.)
      LPROGNOSTIC  = .FALSE.
      LPARAMETRIC  = .FALSE.
      LHYDROPHILIC = .FALSE.
      NAME         = ' '
      AER_TYPE     = 'XX'
      OPTICS_MODEL = ' '
      SIZE_BIN     = 0
      SCALING      = 1.0_JPRB
      PARAMETRIC_SCHEME = -1
      ! Note that there is no IOSTAT argument in the following, in
      ! order that if the read fails the message reports which
      ! namelist entry is incorrect.  Because we have already counted
      ! how many NAM_AEROSOL_RADIATION groups are present, it should
      ! always find the next one.
      READ(UNIT=KUNIT, NML=NAM_AEROSOL_RADIATION)
      SELF%DATA(JAER)%NAME         = TRIM(NAME)
      SELF%DATA(JAER)%LHYDROPHILIC = LHYDROPHILIC
      SELF%DATA(JAER)%SCALING      = SCALING
      SELF%DATA(JAER)%AER_TYPE     = AER_TYPE
      SELF%DATA(JAER)%OPTICS_MODEL = TRIM(OPTICS_MODEL)
      SELF%DATA(JAER)%SIZE_BIN     = SIZE_BIN
      SELF%DATA(JAER)%PARAMETRIC_SCHEME = PARAMETRIC_SCHEME
      IF (LPROGNOSTIC) THEN
        SELF%DATA(JAER)%NSOURCE = IPROGNOSTIC
        SELF%NPROGNOSTIC = SELF%NPROGNOSTIC + 1
      ELSEIF (LPARAMETRIC) THEN
        SELF%DATA(JAER)%NSOURCE = IPARAMETRIC
        SELF%NPARAMETRIC = SELF%NPARAMETRIC + 1
        SELF%NDIAGNOSTIC = SELF%NDIAGNOSTIC + 1
        IF (PARAMETRIC_SCHEME < 0) THEN
          ! Check if name specifies what parametric scheme to use
          IF (LEN_TRIM(NAME) >= 12) THEN
            IF (NAME(1:12) == 'Tropospheric') THEN
              SELF%DATA(JAER)%PARAMETRIC_SCHEME = 0
            ELSEIF (NAME(1:12) == 'Stratospheri') THEN
              SELF%DATA(JAER)%PARAMETRIC_SCHEME = 1
            ELSE
              WRITE(NULERR,'(A,A,A)') 'YEAERRADDESC: NAM_AEROSOL_RADIATION entry "', &
                   &  TRIM(NAME), '" must specify PARAMETRIC_SCHEME'
              CALL ABOR1('YEAERRADDESC: Error reading namelist')
            ENDIF
          ELSE
            WRITE(NULERR,'(A,A,A)') 'YEAERRADDESC: NAM_AEROSOL_RADIATION entry "', &
                 &  TRIM(NAME), '" must specify PARAMETRIC_SCHEME'
            CALL ABOR1('YEAERRADDESC: Error reading namelist')
          ENDIF
        ENDIF
      ELSE
        SELF%DATA(JAER)%NSOURCE = ICLIMATOLOGICAL
        SELF%NCLIMATOLOGICAL    = SELF%NCLIMATOLOGICAL + 1
        SELF%NDIAGNOSTIC        = SELF%NDIAGNOSTIC + 1
      ENDIF
    ENDDO

    SELF%NAEROSOL = NAER
    IF (SELF%NPROGNOSTIC > 0) THEN
      ALLOCATE(SELF%PROG_INDEX(SELF%NPROGNOSTIC))
    ENDIF

  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:READ_NAMELIST',1,ZHOOK_HANDLE)

END SUBROUTINE READ_NAMELIST


! ------------------------------------------------------------------
! Pre-allocate a certain number of aerosol species in a TAER_RAD_DESC
! object
SUBROUTINE RESERVE(SELF, KALLOC)
  CLASS(TAER_RAD_DESC), INTENT(INOUT) :: SELF
  INTEGER(KIND=JPIM),   INTENT(IN)    :: KALLOC
  CALL SELF%CLEAR()
  ALLOCATE(SELF%DATA(KALLOC))
  ALLOCATE(SELF%PROG_INDEX(KALLOC))
  SELF%PROG_INDEX = 0
END SUBROUTINE RESERVE


! ------------------------------------------------------------------
! Deallocate the contents of a TAER_RAD_DESC object
SUBROUTINE CLEAR(SELF)
  CLASS(TAER_RAD_DESC), INTENT(INOUT) :: SELF
  IF (ALLOCATED(SELF%DATA)) THEN
    DEALLOCATE(SELF%DATA)
  ENDIF
  SELF%NAEROSOL    = 0
  SELF%NDIAGNOSTIC = 0
  SELF%NPROGNOSTIC = 0
  IF (ALLOCATED(SELF%PROG_INDEX)) THEN
    DEALLOCATE(SELF%PROG_INDEX)
  ENDIF
END SUBROUTINE CLEAR


! ------------------------------------------------------------------
! Print the contents of a TAER_RAD_DESC object as text to the
! specified unit
SUBROUTINE PRINT(SELF, KUNIT)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

  CLASS(TAER_RAD_DESC), INTENT(IN) :: SELF
  INTEGER(KIND=JPIM),   INTENT(IN) :: KUNIT

  INTEGER(KIND=JPIM)               :: JAER ! Loop counter

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:PRINT',0,ZHOOK_HANDLE)

  IF (SELF%NAEROSOL > 0) THEN
    WRITE(KUNIT,'(A)') 'Aerosol species seen by radiation'
    !WRITE(KUNIT,'(A,I0,A)') '  ', SELF%NAEROSOL, ' total'
    !WRITE(KUNIT,'(A,I0,A,I0,A,I0,A)') '  ', SELF%NDIAGNOSTIC, ' diagnostic (', &
    !     &  SELF%NCLIMATOLOGICAL, ' climatological + ', SELF%NPARAMETRIC, ' parametric)'
    !WRITE(KUNIT,'(A,I0,A)') '  ', SELF%NPROGNOSTIC, ' prognostic'
    DO JAER = 1,SELF%NAEROSOL
      ASSOCIATE(THIS => SELF%DATA(JAER))
        WRITE(KUNIT,'(I3,A,A,A)',advance='no') JAER, '. ', TRIM(THIS%NAME), ': '
        IF (THIS%NSOURCE == IPROGNOSTIC) THEN
          WRITE(KUNIT,'(A)',advance='no') 'prognostic, '
        ELSEIF (THIS%NSOURCE == ICLIMATOLOGICAL) THEN
          WRITE(KUNIT,'(A)',advance='no') 'climatological, '
        ELSE 
          WRITE(KUNIT,'(A)',advance='no') 'parametric, '
        ENDIF
        IF (THIS%LHYDROPHILIC) THEN
          WRITE(KUNIT,'(A)',advance='no') 'hydrophilic'
        ELSE
          WRITE(KUNIT,'(A)',advance='no') 'hydrophobic'
        ENDIF
        IF (THIS%SCALING /= 1.0) THEN
          WRITE(KUNIT,'(A,F7.3)',advance='no') ', scaling=', THIS%SCALING
        END IF
        WRITE(KUNIT,'(A,A)',advance='no') ', type=', THIS%AER_TYPE
        IF (THIS%OPTICS_MODEL /= ' ') THEN
          WRITE(KUNIT,'(A,A)',advance='no') ', optics_model=', TRIM(THIS%OPTICS_MODEL)
        END IF
        IF (THIS%SIZE_BIN > 0) THEN
          WRITE(KUNIT,'(A,I0)',advance='no') ', bin=', THIS%SIZE_BIN
        ENDIF
        IF (THIS%NSOURCE == IPARAMETRIC) THEN
          WRITE(KUNIT,'(A,I0)',advance='no') ', parametric_scheme=', THIS%PARAMETRIC_SCHEME
        ENDIF
        WRITE(KUNIT,'(A)') ' '
      END ASSOCIATE
    ENDDO
  ELSE
    WRITE(KUNIT,'(A)') 'No aerosol species seen by radiation'
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:PRINT',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT


! ------------------------------------------------------------------
! Print the contents of a TAER_RAD_DESC object as namelist-formatted
! text to the specified unit
SUBROUTINE PRINT_NAMELIST(SELF, KUNIT)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

  CLASS(TAER_RAD_DESC), INTENT(IN) :: SELF
  INTEGER(KIND=JPIM),   INTENT(IN) :: KUNIT
  INTEGER(KIND=JPIM)               :: JAER ! Loop counter

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:PRINT_NAMELIST',0,ZHOOK_HANDLE)

  IF (SELF%NAEROSOL > 0) THEN
    WRITE(KUNIT,'(A)') '! Aerosol species seen by radiation'
    DO JAER = 1,SELF%NAEROSOL
      ASSOCIATE(THIS => SELF%DATA(JAER))
        WRITE(KUNIT,'(A)') '&NAM_AEROSOL_RADIATION'
        WRITE(KUNIT,'(A,A,A)') '  NAME = "', TRIM(THIS%NAME), '",'
        IF (THIS%NSOURCE == IPROGNOSTIC) THEN
          WRITE(KUNIT,'(A)') '  LPROGNOSTIC = true,'
        ELSEIF (THIS%NSOURCE == ICLIMATOLOGICAL) THEN
          WRITE(KUNIT,'(A)') '  LPROGNOSTIC = false,'
        ELSE 
          WRITE(KUNIT,'(A)') '  LPARAMETRIC = true,'
          WRITE(KUNIT,'(A,I0,A)') '  PARAMETRIC_SCHEME = ', THIS%PARAMETRIC_SCHEME, ','
        ENDIF
        IF (THIS%LHYDROPHILIC) THEN
          WRITE(KUNIT,'(A)') '  LHYDROPHILIC = true,'
        ELSE
          WRITE(KUNIT,'(A)') '  LHYDROPHILIC = false,'
        ENDIF
        IF (THIS%SCALING /= 1.0) THEN
          WRITE(KUNIT,'(A,F7.3,A)') '  SCALING = ', THIS%SCALING, ','
        END IF
        WRITE(KUNIT,'(A,A,A)') '  AER_TYPE = "', THIS%AER_TYPE, '",'
        IF (THIS%OPTICS_MODEL /= ' ') THEN
          WRITE(KUNIT,'(A,A,A)') '  OPTICS_MODEL = "', TRIM(THIS%OPTICS_MODEL), '",'
        END IF
        IF (THIS%SIZE_BIN > 0) THEN
          WRITE(KUNIT,'(A,I0,A)') '  SIZE_BIN = ', THIS%SIZE_BIN, ','
        ENDIF
        WRITE(KUNIT,'(A)') '/'
      END ASSOCIATE
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:PRINT_NAMELIST',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_NAMELIST


! ------------------------------------------------------------------
! Add a climatological aerosol to a TAER_RAD_DESC object
SUBROUTINE ADD_CLIMATOLOGICAL(SELF, NAME, AER_TYPE, LHYDROPHILIC, &
                           &  BIN, OPTICS_MODEL, SCALING, KINDEX, LWRITEOVER)
  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

  CLASS(TAER_RAD_DESC), INTENT(INOUT)        :: SELF
  CHARACTER(LEN=*),     INTENT(IN)           :: NAME
  CHARACTER(LEN=2),     INTENT(IN)           :: AER_TYPE
  LOGICAL,              INTENT(IN)           :: LHYDROPHILIC
  INTEGER(KIND=JPIM),   INTENT(IN), OPTIONAL :: BIN
  CHARACTER(LEN=*),     INTENT(IN), OPTIONAL :: OPTICS_MODEL
  REAL(KIND=JPRB),      INTENT(IN), OPTIONAL :: SCALING
  INTEGER(KIND=JPIM),   INTENT(OUT),OPTIONAL :: KINDEX
  LOGICAL           ,   INTENT(IN), OPTIONAL :: LWRITEOVER
  INTEGER(KIND=JPIM)                         :: AINDEX, CBIN
  LOGICAL                                    :: LWRITE
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_CLIMATOLOGICAL',0,ZHOOK_HANDLE)

  IF (.NOT. ALLOCATED(SELF%DATA)) THEN
    CALL ABOR1('Attempt to ADD_CLIMATOLOGICAL aerosol before ALLOCATE')
  ELSEIF (SIZE(SELF%DATA) <= SELF%NAEROSOL) THEN
    CALL ABOR1('TAER_RAD_DESC has insufficient allocation to ADD_CLIMATOLOGICAL aerosol again') 
  ENDIF

  IF (.NOT. PRESENT(LWRITEOVER)) THEN
    LWRITE=.FALSE.
  ELSE
    LWRITE=LWRITEOVER
  ENDIF

  IF (.NOT. PRESENT(BIN)) THEN
    CBIN=0
  ELSE
    CBIN=BIN
  ENDIF
  AINDEX = SELF%GET_INDEX(AER_TYPE, LHYDROPHILIC, CBIN)

  ! If aerosol not added, return zero
  IF (PRESENT(KINDEX)) THEN
    KINDEX = 0
  ENDIF
    
  IF (AINDEX == 0 .OR. LWRITE) THEN
    SELF%NAEROSOL        = SELF%NAEROSOL + 1
    SELF%NCLIMATOLOGICAL = SELF%NCLIMATOLOGICAL + 1
    SELF%NDIAGNOSTIC     = SELF%NDIAGNOSTIC + 1
    ASSOCIATE(THIS => SELF%DATA(SELF%NAEROSOL))
      THIS%NAME = NAME
      THIS%NSOURCE = ICLIMATOLOGICAL
      THIS%AER_TYPE = AER_TYPE
      THIS%LHYDROPHILIC = LHYDROPHILIC
      IF (PRESENT(OPTICS_MODEL)) THEN
        THIS%OPTICS_MODEL = OPTICS_MODEL
      ENDIF
      IF (PRESENT(BIN)) THEN
        THIS%SIZE_BIN = BIN
      ENDIF
      IF (PRESENT(SCALING)) THEN
        THIS%SCALING = SCALING
      END IF
    END ASSOCIATE

    IF (PRESENT(KINDEX)) THEN
      KINDEX = SELF%NAEROSOL
    ENDIF
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_CLIMATOLOGICAL',1,ZHOOK_HANDLE)

END SUBROUTINE ADD_CLIMATOLOGICAL


! ------------------------------------------------------------------
! Add a parametric aerosol to a TAER_RAD_DESC object
SUBROUTINE ADD_PARAMETRIC(SELF, NAME, PARAMETRIC_SCHEME, AER_TYPE, LHYDROPHILIC, &
     &                    BIN, OPTICS_MODEL, SCALING, KINDEX)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

  CLASS(TAER_RAD_DESC), INTENT(INOUT)        :: SELF
  CHARACTER(LEN=*),     INTENT(IN)           :: NAME
  INTEGER(KIND=JPIM),   INTENT(IN)           :: PARAMETRIC_SCHEME
  CHARACTER(LEN=2),     INTENT(IN)           :: AER_TYPE
  LOGICAL,              INTENT(IN)           :: LHYDROPHILIC
  INTEGER(KIND=JPIM),   INTENT(IN), OPTIONAL :: BIN
  CHARACTER(LEN=*),     INTENT(IN), OPTIONAL :: OPTICS_MODEL
  REAL(KIND=JPRB),      INTENT(IN), OPTIONAL :: SCALING
  INTEGER(KIND=JPIM),   INTENT(OUT),OPTIONAL :: KINDEX

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_PARAMETRIC',0,ZHOOK_HANDLE)

  IF (.NOT. ALLOCATED(SELF%DATA)) THEN
    CALL ABOR1('Attempt to ADD_PARAMETRIC aerosol before ALLOCATE')
  ELSEIF (SIZE(SELF%DATA) <= SELF%NAEROSOL) THEN
    CALL ABOR1('TAER_RAD_DESC has insufficient allocation to ADD_PARAMETRIC aerosol again') 
  ENDIF

  SELF%NAEROSOL    = SELF%NAEROSOL + 1
  SELF%NPARAMETRIC = SELF%NPARAMETRIC + 1
  SELF%NDIAGNOSTIC = SELF%NDIAGNOSTIC + 1
  ASSOCIATE(THIS => SELF%DATA(SELF%NAEROSOL))
    THIS%NAME = NAME
    THIS%NSOURCE = IPARAMETRIC
    THIS%PARAMETRIC_SCHEME = PARAMETRIC_SCHEME
    THIS%AER_TYPE = AER_TYPE
    THIS%LHYDROPHILIC = LHYDROPHILIC
    IF (PRESENT(OPTICS_MODEL)) THEN
      THIS%OPTICS_MODEL = OPTICS_MODEL
    ENDIF
    IF (PRESENT(BIN)) THEN
      THIS%SIZE_BIN = BIN
    ENDIF
    IF (PRESENT(SCALING)) THEN
      THIS%SCALING = SCALING
    END IF
  END ASSOCIATE

  IF (PRESENT(KINDEX)) THEN
    KINDEX = SELF%NAEROSOL
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_PARAMETRIC',1,ZHOOK_HANDLE)

END SUBROUTINE ADD_PARAMETRIC


! ------------------------------------------------------------------
! Add a prognostic aerosol species to a TAER_RAD_DESC object
SUBROUTINE ADD_PROGNOSTIC(SELF, NAME, AER_TYPE, LHYDROPHILIC, &
     &                    BIN, OPTICS_MODEL, SCALING, KINDEX)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

  CLASS(TAER_RAD_DESC), INTENT(INOUT)        :: SELF
  CHARACTER(LEN=*),     INTENT(IN)           :: NAME
  CHARACTER(LEN=2),     INTENT(IN)           :: AER_TYPE
  LOGICAL,              INTENT(IN)           :: LHYDROPHILIC
  INTEGER(KIND=JPIM),   INTENT(IN), OPTIONAL :: BIN
  CHARACTER(LEN=*),     INTENT(IN), OPTIONAL :: OPTICS_MODEL
  REAL(KIND=JPRB),      INTENT(IN), OPTIONAL :: SCALING
  INTEGER(KIND=JPIM),   INTENT(OUT),OPTIONAL :: KINDEX

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_PROGNOSTIC',0,ZHOOK_HANDLE)

  IF (.NOT. ALLOCATED(SELF%DATA)) THEN
    CALL ABOR1('Attempt to ADD_PROGNOSTIC aerosol before ALLOCATE')
  ELSEIF (SIZE(SELF%DATA) <= SELF%NAEROSOL) THEN
    CALL ABOR1('TAER_RAD_DESC has insufficient allocation to ADD_PROGNOSTIC aerosol again') 
  ENDIF

  SELF%NAEROSOL    = SELF%NAEROSOL + 1
  SELF%NPROGNOSTIC = SELF%NPROGNOSTIC + 1
  ASSOCIATE(THIS => SELF%DATA(SELF%NAEROSOL))
    THIS%NAME = NAME
    THIS%NSOURCE = IPROGNOSTIC
    THIS%AER_TYPE = AER_TYPE
    THIS%LHYDROPHILIC = LHYDROPHILIC
    IF (PRESENT(OPTICS_MODEL)) THEN
      THIS%OPTICS_MODEL = OPTICS_MODEL
    ENDIF
    IF (PRESENT(BIN)) THEN
      THIS%SIZE_BIN = BIN
    ENDIF
    IF (PRESENT(SCALING)) THEN
      THIS%SCALING = SCALING
    END IF
  END ASSOCIATE

  IF (PRESENT(KINDEX)) THEN
    KINDEX = SELF%NAEROSOL
  ENDIF

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_PROGNOSTIC',1,ZHOOK_HANDLE)

END SUBROUTINE ADD_PROGNOSTIC


! ------------------------------------------------------------------
! Add all prognostic aerosols in the YAERO_DESC array to
! SELF
SUBROUTINE ADD_ALL_PROGNOSTICS(SELF, YAERO_DESC)

  USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOEAERATM, ONLY : TYPE_AERO_DESC
  USE YOMLUN,    ONLY : NULERR

  CLASS(TAER_RAD_DESC), INTENT(INOUT) :: SELF
  TYPE(TYPE_AERO_DESC), INTENT(IN)    :: YAERO_DESC(:)

  INTEGER(KIND=JPIM) :: JAER

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_ALL_PROGNOSTICS',0,ZHOOK_HANDLE)

  IF (ALLOCATED(SELF%DATA)) THEN
    ! Check that enough slots have been added
    IF (SIZE(SELF%DATA)-SELF%NAEROSOL < SIZE(YAERO_DESC)) THEN
      WRITE(NULERR,'(A)') 'YOEAERRADDESC: Not enough slots in SELF to hold all prognostic aerosols; call RESERVE with more'
      CALL ABOR1('YOEAERRADDESC: abort')
    ENDIF
  ELSE
    CALL SELF%RESERVE(SIZE(YAERO_DESC))
  ENDIF
  
  DO JAER = 1,SIZE(YAERO_DESC)
    SELF%NAEROSOL    = SELF%NAEROSOL+1
    SELF%NPROGNOSTIC = SELF%NPROGNOSTIC+1
    ASSOCIATE(THIS => SELF%DATA(SELF%NAEROSOL))
      SELF%PROG_INDEX(JAER) = JAER
      THIS%NAME         = TRIM(YAERO_DESC(JAER)%CNAME)
      THIS%NSOURCE      = IPROGNOSTIC
      ! Default settings for size and hydrophilic/hydrophobic; might
      ! be overwritten
      THIS%SIZE_BIN     = YAERO_DESC(JAER)%NBIN
      IF (YAERO_DESC(JAER)%CHYGCLASS == "Hydrophobic") THEN
        THIS%LHYDROPHILIC = .FALSE.
      ELSE
        THIS%LHYDROPHILIC = .TRUE.
      ENDIF

      SELECT CASE (YAERO_DESC(JAER)%NTYP)
        CASE (1)
          THIS%AER_TYPE = 'SS' ! Sea salt
        CASE (2)
          THIS%AER_TYPE = 'DD' ! Desert dust
        CASE (3)
          THIS%AER_TYPE = 'OM' ! Organic matter
        CASE (4)
          THIS%AER_TYPE = 'BC' ! Black carbon
          ! We can't do hydrophilic BC in radiation yet, but need to
          ! keep the possibility of the flag being TRUE in order that
          ! the BC hydrophilic/hydrophobic types are
          ! distinguishable. A subsequent call to
          ! SET_HYDROPHOBIC('BC') ensures the radiation scheme treats
          ! it correctly. Therefore the following is commented.
          !THIS%LHYDROPHILIC = .FALSE.
        CASE (5)
          THIS%AER_TYPE = 'SU' ! Sulfate
        CASE (6)
          THIS%AER_TYPE = 'NI' ! Nitrate
        CASE (7)
          THIS%AER_TYPE = 'AM' ! Ammonium
        CASE (8) ! Secondary organic aerosol (SOA)...
          IF (YAERO_DESC(JAER)%NBIN == 1) THEN
            THIS%AER_TYPE = 'OB' ! ...of biogenic origin
          ELSE
            THIS%AER_TYPE = 'OA' ! ...of anthropogenic origin
          ENDIF
        CASE (9) ! Volcanic fly ash, treat as coarse desert dust
          THIS%AER_TYPE = 'DD'
          THIS%SIZE_BIN = 3
          THIS%LHYDROPHILIC = .FALSE.
        CASE (10) ! Volcanic SO4/SO2
          IF (YAERO_DESC(JAER)%NBIN == 1) THEN
            THIS%AER_TYPE = 'SU' ! SO4: Treat as ammonium sulphate
          ELSE
            ! SO2 precursor: radiatively inactive, so decrement the
            ! counter
            SELF%NAEROSOL    = SELF%NAEROSOL-1
            SELF%NPROGNOSTIC = SELF%NPROGNOSTIC-1
          ENDIF
        CASE DEFAULT
          WRITE(NULERR,'(A,A,A,I0,A)') 'YOEAERRADDESC: ', TRIM(YAERO_DESC(JAER)%CNAME), &
               &  ' has NTYP=', YAERO_DESC(JAER)%NTYP, &
               &  ': not understood'
          CALL ABOR1('YOEAERRADDESC: abort')
      END SELECT
    END ASSOCIATE

  ENDDO

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:ADD_ALL_PROGNOSTICS',1,ZHOOK_HANDLE)

END SUBROUTINE ADD_ALL_PROGNOSTICS


! ------------------------------------------------------------------
! Find the prognostic aerosols requested by a TAER_RAD_DESC object in
! the list of aerosol GFL fields in YAERO_DESC and store their indices
SUBROUTINE FIND_PROGNOSTICS(SELF, YAERO_DESC)

  USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOEAERATM, ONLY : TYPE_AERO_DESC
  USE YOMLUN,    ONLY : NULERR

  CLASS(TAER_RAD_DESC), INTENT(INOUT) :: SELF
  TYPE(TYPE_AERO_DESC), INTENT(IN)    :: YAERO_DESC(:)

  INTEGER(KIND=JPIM) :: JAER, JGFL, IPROG
  LOGICAL :: LFOUND

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:FIND_PROGNOSTICS',0,ZHOOK_HANDLE)

  IPROG = 1
  DO JAER = 1,SELF%NAEROSOL
    IF (SELF%DATA(JAER)%NSOURCE == IPROGNOSTIC) THEN
      LFOUND = .FALSE.
      DO JGFL = 1,SIZE(YAERO_DESC)
        IF (TRIM(SELF%DATA(JAER)%NAME) == TRIM(YAERO_DESC(JGFL)%CNAME)) THEN
          SELF%PROG_INDEX(IPROG) = JGFL
          LFOUND = .TRUE.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT. LFOUND) THEN
        WRITE(NULERR,'(A,A,A)') 'YOEAERRADDESC: requested prognostic aerosol ', &
             &    TRIM(SELF%DATA(JAER)%NAME), ' not found in GFL fields'
        CALL ABOR1('YOEAERRADDESC: abort')
      ENDIF
      IPROG = IPROG+1
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:FIND_PROGNOSTICS',1,ZHOOK_HANDLE)

END SUBROUTINE FIND_PROGNOSTICS


! Some species, specifically black carbon, cannot be represented as
! hydrophilic in radiation even though the prognostic aerosol package
! classes them as such
SUBROUTINE SET_HYDROPHOBIC(SELF, AER_TYPE)

  USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN,  ONLY : NULOUT
  
  CLASS(TAER_RAD_DESC), INTENT(INOUT)        :: SELF
  CHARACTER(LEN=2),     INTENT(IN)           :: AER_TYPE

  INTEGER(KIND=JPIM) :: JAER
  
  LOGICAL :: LFOUND
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:SET_HYDROPHOBIC',0,ZHOOK_HANDLE)

  IF (.NOT. ALLOCATED(SELF%DATA)) THEN
    CALL ABOR1('Attempt to SET_HYDROPHOBIC before ALLOCATE')
  ENDIF

  LFOUND = .FALSE.
  DO JAER = 1,SELF%NAEROSOL
    IF (SELF%DATA(JAER)%AER_TYPE == AER_TYPE &
         &  .AND. SELF%DATA(JAER)%LHYDROPHILIC) THEN
      LFOUND = .TRUE.
      SELF%DATA(JAER)%LHYDROPHILIC = .FALSE.
    ENDIF
  ENDDO

  IF (LFOUND) THEN
    WRITE(NULOUT,'(A)') 'Setting all aerosols of type "', AER_TYPE, '" to hydrophobic in radiation'
  ENDIF
    
  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:SET_HYDROPHOBIC',1,ZHOOK_HANDLE)

END SUBROUTINE SET_HYDROPHOBIC


! ------------------------------------------------------------------
! Return the index to the first parametric aerosol type with a
! particular PARAMETRIC_SCHEME, or zero if not found
FUNCTION PARAMETRIC_INDEX(SELF, PARAMETRIC_SCHEME)

  USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK

  CLASS(TAER_RAD_DESC), INTENT(IN) :: SELF
  INTEGER(KIND=JPIM),   INTENT(IN) :: PARAMETRIC_SCHEME
  INTEGER(KIND=JPIM)               :: PARAMETRIC_INDEX

  INTEGER(KIND=JPIM) :: JAER
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:PARAMETRIC_INDEX',0,ZHOOK_HANDLE)

  PARAMETRIC_INDEX = 0
  DO JAER = 1,SELF%NAEROSOL
    IF (SELF%DATA(JAER)%NSOURCE == IPARAMETRIC &
         &  .AND. SELF%DATA(JAER)%PARAMETRIC_SCHEME == PARAMETRIC_SCHEME) THEN
      PARAMETRIC_INDEX = JAER
      EXIT
    ENDIF
  ENDDO
  
  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:PARAMETRIC_INDEX',1,ZHOOK_HANDLE)
  
END FUNCTION PARAMETRIC_INDEX

! ------------------------------------------------------------------
! Return the index to the first aerosol type with the specified
! type/bin/hydrophilicity, noting that if KBIN is zero then any bin
! will be considered a match
FUNCTION GET_INDEX(SELF, AER_TYPE, LHYDROPHILIC, KBIN)

  USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK

  CLASS(TAER_RAD_DESC), INTENT(IN) :: SELF
  CHARACTER(2),         INTENT(IN) :: AER_TYPE
  LOGICAL,              INTENT(IN) :: LHYDROPHILIC
  INTEGER(KIND=JPIM),   INTENT(IN) :: KBIN
  INTEGER(KIND=JPIM)               :: GET_INDEX

  LOGICAL            :: LCONDTYP, LCONDPHI
  INTEGER(KIND=JPIM) :: JAER
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:GET_INDEX',0,ZHOOK_HANDLE)
 
  GET_INDEX = 0
  DO JAER = 1,SELF%NAEROSOL
    LCONDTYP = SELF%DATA(JAER)%AER_TYPE == AER_TYPE 
    LCONDPHI = SELF%DATA(JAER)%LHYDROPHILIC .EQV. LHYDROPHILIC
    IF (LCONDTYP .AND. LCONDPHI) THEN
      IF (KBIN > 0) THEN
        IF (SELF%DATA(JAER)%SIZE_BIN == KBIN) THEN
          GET_INDEX = JAER
        ENDIF
      ELSE
        GET_INDEX = JAER
      ENDIF
    ENDIF
  ENDDO
  
  IF (LHOOK) CALL DR_HOOK('YOEAERRADDESC:GET_INDEX',1,ZHOOK_HANDLE)
  
END FUNCTION GET_INDEX


END MODULE YOEAERRADDESC
