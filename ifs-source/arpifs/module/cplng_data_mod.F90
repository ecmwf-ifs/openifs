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

MODULE CPLNG_DATA_MOD
   
   USE PARKIND1, ONLY: JPIM
   USE PARKIND_OCEAN, ONLY: JPRO

   USE CPLNG_LOG_MOD
   USE CPLNG_TYPES_MOD
   USE YOMMCC, ONLY : TMCC

#ifdef WITH_OASIS
   USE MOD_OASIS
#endif

   IMPLICIT NONE

   PRIVATE

   PUBLIC CPLNG_IS_ACTIVE

   PUBLIC CPLNG_FLD_TYPE
   
   PUBLIC CPLNG_ADD_FLD
   PUBLIC CPLNG_ADD_FLD_COMPLETED

   PUBLIC CPLNG_IDX

   PUBLIC CPLNG_FLD_TYPE_GRIDPOINT
   PUBLIC CPLNG_FLD_TYPE_SPECTRAL

   PUBLIC :: CPL_IN
   PUBLIC :: CPL_OUT
   PUBLIC :: CPL_OUTINST
   
   INTEGER(KIND=JPIM), PARAMETER :: CPLNG_FLD_TYPE_GRIDPOINT = 0
   INTEGER(KIND=JPIM), PARAMETER :: CPLNG_FLD_TYPE_SPECTRAL  = 1

#ifdef WITH_OASIS
   INTEGER(KIND=JPIM), SAVE :: CPL_IN  = OASIS_IN
   INTEGER(KIND=JPIM), SAVE :: CPL_OUT = OASIS_OUT
   INTEGER(KIND=JPIM), SAVE :: CPL_OUTINST = OASIS_OUT
#else
   INTEGER(KIND=JPIM), SAVE :: CPL_IN  = 1
   INTEGER(KIND=JPIM), SAVE :: CPL_OUT = 2
   INTEGER(KIND=JPIM), SAVE :: CPL_OUTINST = 3
#endif

CONTAINS

   FUNCTION CPLNG_IS_ACTIVE(YRMCC)

      LOGICAL :: CPLNG_IS_ACTIVE
      TYPE(TMCC), INTENT(IN) :: YRMCC

      CPLNG_IS_ACTIVE = YRMCC%CPLNG_ACTIVE

   END FUNCTION CPLNG_IS_ACTIVE

   SUBROUTINE CPLNG_ADD_FLD(YRMCC,CDNAME,KTYPE,KINOUT,KSTAGE,KLVL,KCAT,KIDX)
      
      ! PARAMETERS
      INTEGER(KIND=JPIM),PARAMETER :: IALLOCATE_CHUNK = 20
      
      ! ARGUMENTS
      TYPE(TMCC), INTENT(INOUT)              :: YRMCC
      CHARACTER(LEN=*),  INTENT(IN)          :: CDNAME
      INTEGER(KIND=JPIM),INTENT(IN)          :: KTYPE
      INTEGER(KIND=JPIM),INTENT(IN)          :: KINOUT
      INTEGER(KIND=JPIM),INTENT(IN)          :: KSTAGE
      INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KLVL
      INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KCAT
      INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL:: KIDX

      ! LOCALS
      INTEGER(KIND=JPIM)               :: INUM_LVL,INUM_CAT
      TYPE(CPLNG_FLD_TYPE), ALLOCATABLE :: CPLNG_FLD_TMP(:)
      
      ! Early return if CPLNG not active
      IF (.NOT.YRMCC%CPLNG_ACTIVE) THEN
         CALL CPLNG_WARNING('CPLNG_ADD_FLD called but YRMCC%CPLNG_ACTIVE is false.')
         RETURN
      ENDIF
      
      ! Set number of levels/categories to default or optional arguments
      INUM_LVL = 1
      INUM_CAT = 1
      
      IF (PRESENT(KLVL)) INUM_LVL = KLVL
      IF (PRESENT(KCAT)) INUM_CAT = KCAT
      
      ! Make initial allocation, if necessary
      IF (.NOT.ALLOCATED(YRMCC%CPLNG_FLD)) THEN
         ALLOCATE(YRMCC%CPLNG_FLD(IALLOCATE_CHUNK))
      ENDIF

      ! Increase (by reallocation) size of YRMCC%CPLNG_FLD, if necessary
      IF (YRMCC%CPLNG_NUM_FIELDS==SIZE(YRMCC%CPLNG_FLD)) THEN
         ALLOCATE(CPLNG_FLD_TMP(SIZE(YRMCC%CPLNG_FLD)+IALLOCATE_CHUNK))
         CPLNG_FLD_TMP(1:YRMCC%CPLNG_NUM_FIELDS) = YRMCC%CPLNG_FLD
         CALL MOVE_ALLOC(CPLNG_FLD_TMP,YRMCC%CPLNG_FLD) ! Fortran 2003
      ENDIF

      YRMCC%CPLNG_NUM_FIELDS=YRMCC%CPLNG_NUM_FIELDS+1

      YRMCC%CPLNG_FLD(YRMCC%CPLNG_NUM_FIELDS)%NAME    = CDNAME
      YRMCC%CPLNG_FLD(YRMCC%CPLNG_NUM_FIELDS)%TYPE    = KTYPE
      YRMCC%CPLNG_FLD(YRMCC%CPLNG_NUM_FIELDS)%INOUT   = KINOUT
      YRMCC%CPLNG_FLD(YRMCC%CPLNG_NUM_FIELDS)%STAGE   = KSTAGE
      YRMCC%CPLNG_FLD(YRMCC%CPLNG_NUM_FIELDS)%NUM_LVL = INUM_LVL
      YRMCC%CPLNG_FLD(YRMCC%CPLNG_NUM_FIELDS)%NUM_CAT = INUM_CAT

      IF (PRESENT(KIDX)) THEN
         KIDX = YRMCC%CPLNG_NUM_FIELDS
      ENDIF

   END SUBROUTINE CPLNG_ADD_FLD

   FUNCTION CPLNG_IDX(YRMCC,CDFLD_NAME)
      ! Argument
      TYPE(TMCC), INTENT(IN) :: YRMCC
      CHARACTER(LEN=*), INTENT(IN) :: CDFLD_NAME
      ! Result
      INTEGER(KIND=JPIM) :: CPLNG_IDX
      ! Locals
      INTEGER(KIND=JPIM) :: II

      ! Early return if CPLNG not active
      IF (.NOT.YRMCC%CPLNG_ACTIVE) THEN
         CALL CPLNG_WARNING('CPLNG_IDX called but YRMCC%CPLNG_ACTIVE is false.')
         RETURN
      ENDIF

      IF (.NOT.ALLOCATED(YRMCC%CPLNG_FLD)) THEN
         CALL ABOR1("CPLNG_IDX: YRMCC%CPLNG_FLD not allocated upon call")
      ENDIF

      ! Simple loop searching for the right field name
      DO II=1,YRMCC%CPLNG_NUM_FIELDS
         IF (TRIM(CDFLD_NAME) == TRIM(YRMCC%CPLNG_FLD(II)%NAME)) EXIT
      ENDDO

      IF (II>YRMCC%CPLNG_NUM_FIELDS) THEN
         CALL ABOR1("CPLNG_IDX: Field name not found: "//TRIM(CDFLD_NAME))
      ENDIF

      CPLNG_IDX = II
   END FUNCTION CPLNG_IDX

   SUBROUTINE CPLNG_ADD_FLD_COMPLETED(YRMCC,YDGEOMETRY)

      USE GEOMETRY_MOD , ONLY : GEOMETRY
      USE YOMLUN, ONLY: NULOUT

      ! Arguments
      TYPE(TMCC), INTENT(INOUT)    :: YRMCC
      TYPE(GEOMETRY),INTENT(IN)    :: YDGEOMETRY
      ! Variables needed in OASIS calls
      INTEGER(KIND=JPIM)              :: IERROR
      INTEGER(KIND=JPIM), ALLOCATABLE :: IG_PARAL(:)
      INTEGER(KIND=JPIM)              :: IL_PART_ID, IL_PART_ID_GP, IL_PART_ID_SP
      INTEGER(KIND=JPIM)              :: IVAR_NODIMS(2)
      INTEGER(KIND=JPIM)              :: IVAR_ACTUAL_SHAPE(2)

      LOGICAL            :: LLSEGSTART(YDGEOMETRY%YRGEM%NGPTOT)
      INTEGER(KIND=JPIM) :: II,IP
      INTEGER(KIND=JPIM) :: ILVL,ICAT
      INTEGER(KIND=JPIM) :: IOFFSET
      INTEGER(KIND=JPIM) :: IEXTENT
      INTEGER(KIND=JPIM) :: INUMSEGMENTS
      CHARACTER(LEN=128) :: CLFLDNAME
      CHARACTER(LEN=3)   :: CLERRSTR

      INTEGER(KIND=JPIM) :: ISPM

      IL_PART_ID_GP = 0
      IL_PART_ID_SP = 0

      ! Early return if CPLNG not active
      IF (.NOT.YRMCC%CPLNG_ACTIVE) THEN
         CALL CPLNG_WARNING('CPLNG_ADD_FLD_COMPLETED called but YRMCC%CPLNG_ACTIVE is false.')
         RETURN
      ENDIF

      ! -------------------------------------------------------------------------
      ! * (3) SET UP OASIS PARTITION
      ! -------------------------------------------------------------------------

      ! * (3A) Grid point fields
#ifdef WITH_OASIS
      IF( ANY(CPLNG_FLD(1:YRMCC%CPLNG_NUM_FIELDS)%TYPE == CPLNG_FLD_TYPE_GRIDPOINT) )THEN

         ! Set up logical array LLSEGSTART that indicats where OASIS segments start
         ! Method: Use the NGLOBALINDEX array from YOMMP that contains a
         !         local-to-global mapping of the grid points. If two consecutive
         !         indices in this array differ by any other value than one, it's the
         !         start of a new segment.
         LLSEGSTART(1)        = .TRUE.
         LLSEGSTART(2:YDGEOMETRY%YRGEM%NGPTOT) = &
            & YDGEOMETRY%YRMP%NGLOBALINDEX(2:YDGEOMETRY%YRGEM%NGPTOT) -  &
            & YDGEOMETRY%YRMP%NGLOBALINDEX(1:YDGEOMETRY%YRGEM%NGPTOT-1) &
            & /= 1

         INUMSEGMENTS = COUNT(LLSEGSTART) ! Counts number of segments

         ! For shape and meaning of IG_PARAL, see OASIS documentation
         ALLOCATE(IG_PARAL(2 + 2*INUMSEGMENTS))

         IG_PARAL(1) = 3 ! Value 3 indicates an orange partition, which is
         ! an ensemble of segments of the global domain
         IG_PARAL(2) = INUMSEGMENTS

         ! Compute segment offsets and extents and store them in IG_PARAL for later
         ! use in the OASIS_DEF_PARTITION call
         IOFFSET = YDGEOMETRY%YRMP%NGLOBALINDEX(1)
         IEXTENT = 1
         IP      = 3 ! Pointer into IG_PARAL array
         DO II=2,YDGEOMETRY%YRGEM%NGPTOT
            IF (LLSEGSTART(II)) THEN
               ! OASIS counts offsets starting with zero, hence minus one
               IG_PARAL(IP)   = IOFFSET-1
               IG_PARAL(IP+1) = IEXTENT
               IP = IP+2
               IOFFSET = YDGEOMETRY%YRMP%NGLOBALINDEX(II)
               IEXTENT = 1
            ELSE
               IEXTENT = IEXTENT+1
            ENDIF
         ENDDO
         IG_PARAL(IP)   = IOFFSET-1
         IG_PARAL(IP+1) = IEXTENT

         ! Define partition for OASIS
         CALL OASIS_DEF_PARTITION(IL_PART_ID_GP,IG_PARAL,IERROR)
         IF (IERROR/=OASIS_OK) THEN
            WRITE (CLERRSTR,'(I3)') IERROR
            CALL ABOR1("CPLNG_OCE_CONFIG: Error on OASIS_DEF_PARTITION (gridpoint): "//CLERRSTR)
         ENDIF

         DEALLOCATE(IG_PARAL)

      ENDIF ! Field type == gridpoint
#endif

      ! * (3B) Spectral fields
      IF( ANY(YRMCC%CPLNG_FLD(1:YRMCC%CPLNG_NUM_FIELDS)%TYPE == CPLNG_FLD_TYPE_SPECTRAL) )THEN
         !
         ! Orange partition: each processor owns a number of (not following) wave numbers.
         !
         !  yomdim/NSMAX        : Spectral truncation T
         !  yomdim/NUMP         : Number of wavenumbers locally
         !  yomlap/MYMS(1:NUMP) : Actual wave numbers   (NOTE: not 0:NUMP as in comment!)

         ! Allocate definition array
         ALLOCATE( IG_PARAL(2+2*YDGEOMETRY%YRDIM%NUMP) )

         ! Indicator for orange partition
         IG_PARAL(1) = 3

         ! Number of segments
         IG_PARAL(2) = YDGEOMETRY%YRDIM%NUMP

         ! Global offset and size
         DO II = 1, YDGEOMETRY%YRDIM%NUMP
            ! Global wave number
            ISPM = YDGEOMETRY%YRLAP%MYMS(II)
            ! Offset
            IG_PARAL(1+II*2) = ((YDGEOMETRY%YRDIM%NSMAX+1)+(YDGEOMETRY%YRDIM%NSMAX+2-ISPM)) * ISPM
            ! Size
            IG_PARAL(2+II*2) = (YDGEOMETRY%YRDIM%NSMAX+1-ISPM)*2  ! M:T, Re/Im
         ENDDO

#ifdef WITH_OASIS
         ! Define partition for OASIS
         CALL OASIS_DEF_PARTITION(IL_PART_ID_SP,IG_PARAL,IERROR)
         IF (IERROR/=OASIS_OK) THEN
            WRITE (CLERRSTR,'(I3)') IERROR
            CALL ABOR1("CPLNG_OCE_CONFIG: Error on OASIS_DEF_PARTITION (spectral): "//CLERRSTR)
         ENDIF
#else
         CALL ABOR1('CPLNG_DATA_MOD: SGLEXE Do not support spectral fields')
#endif
         DEALLOCATE(IG_PARAL)

      ENDIF ! Field type == spectral

      ! -------------------------------------------------------------------------
      ! * (4) DEFINE COUPLING FIELDS FOR OASIS
      ! -------------------------------------------------------------------------
      IVAR_NODIMS = (/ 1, 1 /)

      WRITE (NULOUT,'(/,1X,A,/,1X,A,I4)')        &
         &     'CPLNG: Coupling field definition.', &
         &     'CPLNG: ... Number of fields is: ',YRMCC%CPLNG_NUM_FIELDS

      DO II=1,YRMCC%CPLNG_NUM_FIELDS

         WRITE (NULOUT,'(1X,A,I4,A,I2,A,I2,A,A20)') &
            &     'CPLNG: .... Field ',II,             &
            &     '(',YRMCC%CPLNG_FLD(II)%NUM_LVL,' lvl,',   &
            &     YRMCC%CPLNG_FLD(II)%NUM_CAT,' cat) Name: ',TRIM(YRMCC%CPLNG_FLD(II)%NAME)

         ! Allocate for ID's of all levels/categories of the field
         ALLOCATE(YRMCC%CPLNG_FLD(II)%ID(YRMCC%CPLNG_FLD(II)%NUM_LVL,YRMCC%CPLNG_FLD(II)%NUM_CAT))

         ! Set OASIS partition id and data size according to field type
         ! (gp/spectral). Allocate data member.
         SELECT CASE (YRMCC%CPLNG_FLD(II)%TYPE)

         CASE (CPLNG_FLD_TYPE_GRIDPOINT)
            IL_PART_ID        = IL_PART_ID_GP
            IVAR_ACTUAL_SHAPE = (/ 1, YDGEOMETRY%YRGEM%NGPTOT /)
            ALLOCATE(YRMCC%CPLNG_FLD(II)%D(YDGEOMETRY%YRGEM%NGPTOT, &
               &                        YRMCC%CPLNG_FLD(II)%NUM_LVL, &
               &                        YRMCC%CPLNG_FLD(II)%NUM_CAT  ))

         CASE (CPLNG_FLD_TYPE_SPECTRAL )
            IL_PART_ID        = IL_PART_ID_SP
            IVAR_ACTUAL_SHAPE = (/ 1, YDGEOMETRY%YRDIM%NSPEC2 /)
            ALLOCATE(YRMCC%CPLNG_FLD(II)%D(YDGEOMETRY%YRDIM%NSPEC2, &
               &                        YRMCC%CPLNG_FLD(II)%NUM_LVL, &
               &                        YRMCC%CPLNG_FLD(II)%NUM_CAT  ))

         CASE DEFAULT
            CALL ABOR1("CPLNG_OCE_CONFIG: Wrong field type for " // TRIM(YRMCC%CPLNG_FLD(II)%NAME))

         END SELECT

         ! Initialise to a strange value
         YRMCC%CPLNG_FLD(II)%D = -HUGE(0.0_JPRO)

#ifdef WITH_OASIS
         DO ICAT=1,YRMCC%CPLNG_FLD(II)%NUM_CAT
            DO ILVL=1,YRMCC%CPLNG_FLD(II)%NUM_LVL

               CLFLDNAME = TRIM(YRMCC%CPLNG_FLD(II)%NAME)
               IF (YRMCC%CPLNG_FLD(II)%NUM_CAT>1) WRITE (CLFLDNAME,'(A,".C",I3.3)') TRIM(CLFLDNAME),ICAT
               IF (YRMCC%CPLNG_FLD(II)%NUM_LVL>1) WRITE (CLFLDNAME,'(A,".L",I3.3)') TRIM(CLFLDNAME),ILVL

               CALL OASIS_DEF_VAR(YRMCC%CPLNG_FLD(II)%ID(ILVL,ICAT), &
                  &                  CLFLDNAME,                    &
                  &                  IL_PART_ID,                  &
                  &                  IVAR_NODIMS,                 &
                  &                  YRMCC%CPLNG_FLD(II)%INOUT,         &
                  &                  IVAR_ACTUAL_SHAPE,           &
                  &                  OASIS_REAL,                  &
                  &                  IERROR                       )
               WRITE (NULOUT,'(1X,A,I4,A,I3)') &
                  &     'CPLNG: .... Field ',II,  &
                  &     ' OASIS_DEF_VAR returns ',IERROR

               IF (IERROR/=OASIS_OK) THEN
                  WRITE (CLERRSTR,'(I3)') IERROR
                  CALL ABOR1("CPLNG_OCE_CONFIG: Error in OASIS_DEF_VAR: " // CLERRSTR // &
                     &          " (" // TRIM(CLFLDNAME) // ")")
               ENDIF

            ENDDO
         ENDDO
#endif

      ENDDO

      ! -------------------------------------------------------------------------
      ! * (5) FINALISE OASIS DEFINITION PHASE
      ! -------------------------------------------------------------------------
#ifdef WITH_OASIS
      CALL OASIS_ENDDEF(IERROR)
      IF (IERROR/=OASIS_OK) THEN
         WRITE (CLERRSTR,'(I3)') IERROR
         CALL ABOR1("CPLNG_OCE_CONFIG: Error in OASIS_ENDDEF: "//CLERRSTR)
      ENDIF
#endif

   END SUBROUTINE CPLNG_ADD_FLD_COMPLETED

END MODULE CPLNG_DATA_MOD
