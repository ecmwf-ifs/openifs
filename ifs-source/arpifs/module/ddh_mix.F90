! (C) Copyright 1989- Meteo-France.

!     ##################
      MODULE DDH_MIX
!     ##################


!!    PURPOSE
!!    -------
!       module for type structure in  DDH

!!**  IMPLICIT ARGUMENTS
!!    ------------------

!!    AUTHOR
!!    ------
!!    O.Riviere   *Meteo France*
!!
!!    Modifications :
!!    --------------
!!    F.Voitus     18/09/17 : New DDH flexible structure for OpenMP
  
!-------------------------------------------------------------------------------


USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK,DR_HOOK, JPHOOK
USE YOMLUN_IFSAUX, ONLY : NULOUT

IMPLICIT NONE

SAVE


!*       1.   TYPE DEFINITION


TYPE DDHFLEX
  CHARACTER(LEN=11)::CNAME !name of field
  CHARACTER(LEN=1)::CFLUX !'F' if flux 'V' if variable 'T' if tendency
  CHARACTER(LEN=3)::CMOD  ! 'ARP','ARO': name of model
  LOGICAL:: LKDDH !TRUE if to be stored into DDH

! rfield has to be a pointer because allocatable not allowed in structure type
  REAL(KIND=JPRB),DIMENSION(:,:),POINTER:: RFIELD=>NULL()  ! value of retrieved field
  INTEGER(KIND=JPIM):: NFIELDIND! position of flux in ddh array
END TYPE DDHFLEX


! structure for vertical profile
TYPE TYP_FIELD3D
  REAL(KIND=JPRB), POINTER  :: RVAL(:,:)=>NULL() ! values (domains x vertical)
  CHARACTER(LEN=16)         :: CNAME             ! name
  CHARACTER(LEN=1)          :: CTYPE             ! indicate wether it's a flux, variable or tendency
END TYPE TYP_FIELD3D

TYPE TYP_FIELD2D
  REAL(KIND=JPRB), POINTER :: RVAL(:)=>NULL() ! values (domains)
  CHARACTER(LEN=16)        :: CNAME           ! name
  CHARACTER(LEN=1)         :: CTYPE           ! indicate wether it's a flux, variable or tendency
END TYPE TYP_FIELD2D

TYPE TYP_BUVAR
  REAL(KIND=JPRB), POINTER  :: RVAL(:,:)=>NULL() ! values (domains)
  CHARACTER(LEN=2)          :: CNAME             ! name
END TYPE TYP_BUVAR


! structure containing fields, weights for horizontal averaging, etc.
TYPE TYP_DDH
  INTEGER(KIND=JPIM)          :: NLEV                  ! vertical dimension
  INTEGER(KIND=JPIM)          :: NPROMA                ! horizontal dimension                                                    
  INTEGER(KIND=JPIM)          :: KST                   ! first point
  INTEGER(KIND=JPIM)          :: KEND                  ! last point

  ! 3D field info
  INTEGER(KIND=JPIM)          :: NFIELDS3D             ! number of fields
  INTEGER(KIND=JPIM)          :: NFIELDS3D_OFFSET      ! offset value for the fields when 
                                                       ! putting in the global arrays
  INTEGER(KIND=JPIM)          :: NFIELDS3D_AUTO        ! number of fields that should be 
                                                       ! automatically allocated 
                                                       ! (maximum of cpg and cpglag)
  TYPE(TYP_FIELD3D), POINTER  :: YFIELD3D(:)=>NULL()   ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVAL3D(:,:,:)=>NULL() ! auxiliary array for fast 
                                                       ! (automatic) allocation
                                                       ! (vertical x fields x domains)
  ! 2D field info
  INTEGER(KIND=JPIM)          :: NFIELDS2D             ! number of fields
  INTEGER(KIND=JPIM)          :: NFIELDS2D_OFFSET      ! offset value for the fields when 
                                                       ! putting in the global arrays
  INTEGER(KIND=JPIM)          :: NFIELDS2D_AUTO        ! number of fields that should be 
                                                       ! automatically allocated 
                                                       ! (maximum of cpg and cpglag)
  TYPE(TYP_FIELD2D), POINTER  :: YFIELD2D(:)=>NULL()   ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVAL2D(:,:)=>NULL()   ! auxiliary array for fast 
                                                       ! (automatic) allocation
                                                       ! (vertical x fields x domains)
  ! horizontal info
  REAL(KIND=JPRB), POINTER    :: WEIGHT(:)=>NULL()     ! weights inside one NPROMA block
  INTEGER(KIND=JPIM), POINTER :: NDDHI(:)=>NULL()      ! cfr. KDDHI in cpg_dia


  TYPE(TYP_BUVAR), POINTER    :: YVARMULT(:)=>NULL()    ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVARSM(:,:,:,:)=>NULL()! auxiliary array for fast 
                                                        ! (automatic) allocation

END TYPE TYP_DDH



!*       2.   DECLARATIONS

REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE,TARGET::RDDH_FIELD
! ddh_field is the target of rfield

TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: RDDH_DESCR ! array of several structure type
INTEGER(KIND=JPIM):: NFLEVGDDH=0   ! number of processed field
INTEGER(KIND=JPIM):: NPROMADDH=0   ! number of processed field
INTEGER(KIND=JPIM):: NFIELDDDH=0   ! number of processed field
INTEGER(KIND=JPIM):: NTOTFIELD=0 !nb of fields in structure
INTEGER(KIND=JPIM):: NTOTVAR=0 ! nb of var. in structure
INTEGER(KIND=JPIM):: NTOTFSL=0 ! nb of SL fields. in structure
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDH ! mask for operations in ddh operators
CHARACTER(LEN=13)::CNAMEINI ! name of first field
LOGICAL:: LALLOC=.FALSE.  !if TRUE allocation of arrrays completed


REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE,TARGET::RDDHSURF_FIELD
! surface ddh_field is the target of rfield
TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: RDDHSURF_DESCR ! array of several structure type
INTEGER(KIND=JPIM):: NSURFDDH=0   ! number of processed surface field
INTEGER(KIND=JPIM):: NTOTSURF=0 !nb of surface fields in structure
INTEGER(KIND=JPIM):: NTOTSVAR=0 ! nb of surface var. in structure
INTEGER(KIND=JPIM):: NTOTSVFS=0 ! nb of surface free style var. in structure
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDHVAS ! mask for operations in ddh operators
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDHVFS ! mask for operations in ddh operators
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDHSURF ! mask for operations in ddh operators
CHARACTER(LEN=13)::CNAMESURFINI ! name of first surface field
LOGICAL:: LALLOCSURF=.FALSE.  !if TRUE allocation of arrrays completed


CONTAINS

!-----------------------------------------------
!*       3.   SUBROUTINES
!-----------------------------------------------

!-----------------------------------------------
! SUBROUTINE ADD_FIELD_3D
! allocates on the fly the flexible DDH structure
! during first time step and then stores the
! velues of retrieved fields into it
! PMAT= Value of field to retrieve
! CDNAME= Name of field
! CDFLUX='F','V','T' if flux,variable or tendency
! LDINST=TRUE if instantaneous field
! LDDH=TRUE if field to be stored into ddh files


SUBROUTINE ADD_FIELD_3D(YDLDDH,PMAT,CDNAME,CDFLUX,CDMOD,LDINST,LDDH,LDSLFD)
USE YOMLDDH , ONLY : TLDDH
TYPE(TLDDH)                   ,INTENT(IN)   :: YDLDDH
REAL(KIND=JPRB),DIMENSION(:,:),INTENT(IN)   :: PMAT
CHARACTER(LEN=*)              ,INTENT(IN)   :: CDNAME
CHARACTER(LEN=1)              ,INTENT(IN)   :: CDFLUX
CHARACTER(LEN=3)              ,INTENT(IN)   :: CDMOD
LOGICAL                       ,INTENT(IN)   :: LDINST,LDDH
LOGICAL              ,OPTIONAL,INTENT(IN)   :: LDSLFD

REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE::ZDDH_FIELD_BKP
INTEGER(KIND=JPIM)::IDIM1,IDIM2,IDIM3,ITEST
TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: YLDDH_DESCR_BKP
LOGICAL::LLSLFD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('DDH_MIX:ADD_FIELD_3D',0,ZHOOK_HANDLE)

IF (YDLDDH%LDDH_OMP) THEN
  !WRITE (NULOUT,*) 'WARNING: ADD_FIELD_3D disabled for LDDH_OMP=.TRUE.'
  IF (LHOOK) CALL DR_HOOK('DDH_MIX:ADD_FIELD_3D',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IF (PRESENT(LDSLFD)) THEN
  LLSLFD=LDSLFD
ELSE
  LLSLFD=.FALSE.
ENDIF


IF(LALLOC) THEN
  ITEST=NFIELDDDH+1
  IF(NFIELDDDH==NTOTFIELD) ITEST=NFIELDDDH+1-NTOTFIELD
  IF(CDNAME /= RDDH_DESCR(ITEST)%CNAME) THEN
    WRITE(NULOUT,*) 'PB NOMS ADD_FIELD3D ',CDNAME,' vs ',RDDH_DESCR(ITEST)%CNAME
    CALL ABOR1('ADD_FIELD_3D')
  ENDIF
ENDIF

IF(.NOT.LALLOC) THEN
  IF (CNAMEINI == CDNAME) THEN
    LALLOC=.TRUE.
  ENDIF
ENDIF

IF (LALLOC .AND. NFIELDDDH == NTOTFIELD) THEN
  NFIELDDDH=NFIELDDDH-NTOTFIELD
  RDDH_FIELD(:,:,:)=0.0
ENDIF



IDIM1=NPROMADDH

IDIM2=NFLEVGDDH+1! arrays are stored on (nflevg+1) levels
                 ! with eventually last level equal to zero
                 ! for tendencies and variables


IF (.NOT.LALLOC) THEN ! (re-)allocation of arrays
  NTOTFIELD=NTOTFIELD+1 !updates total number of fields
  IF(CDFLUX=='V') NTOTVAR=NTOTVAR+1
  IF(LLSLFD) NTOTFSL=NTOTFSL+1

  IF (NFIELDDDH == 0)  THEN
    ALLOCATE(RDDH_DESCR(1))


    ALLOCATE(RDDH_FIELD(IDIM1,IDIM2,1))
    RDDH_FIELD(:,:,:)=0.
    CNAMEINI=CDNAME ! stores name of first field
    IDIM3=0

  ELSE !NFIELDDDH=0
        ! add a field to rddh_descr
    ALLOCATE(YLDDH_DESCR_BKP(NFIELDDDH))
    YLDDH_DESCR_BKP=RDDH_DESCR
    DEALLOCATE(RDDH_DESCR)
    ALLOCATE(RDDH_DESCR(NFIELDDDH+1))
    RDDH_DESCR(1:NFIELDDDH)=YLDDH_DESCR_BKP(1:NFIELDDDH)
    DEALLOCATE(YLDDH_DESCR_BKP)

    IDIM3=SIZE(RDDH_FIELD,3)
    ! reallocation of rddh_field in order to increase dimension
    ALLOCATE(ZDDH_FIELD_BKP(IDIM1,IDIM2,IDIM3))
    ZDDH_FIELD_BKP=RDDH_FIELD ! temporary bkp of rddh_field
    DEALLOCATE(RDDH_FIELD)
    ALLOCATE(RDDH_FIELD(IDIM1,IDIM2,IDIM3+1)) !one field added
    RDDH_FIELD(:,:,:)=0.0
    RDDH_FIELD(1:IDIM1,1:IDIM2,1:IDIM3)=ZDDH_FIELD_BKP(1:IDIM1,1:IDIM2,1:IDIM3)
    DEALLOCATE(ZDDH_FIELD_BKP)

  ENDIF !NFIELDDDH=0

  IF (CDFLUX =='F') THEN
    RDDH_FIELD(1:SIZE(PMAT,1),1:IDIM2,IDIM3+1)=PMAT(1:SIZE(PMAT,1),1:IDIM2) !stores new flux
  ELSE
    RDDH_FIELD(1:SIZE(PMAT,1),1:IDIM2-1,IDIM3+1)=PMAT(1:SIZE(PMAT,1),1:IDIM2-1) !stores new tend/var
    RDDH_FIELD(1:SIZE(PMAT,1),IDIM2,IDIM3+1)=0.0
  ENDIF

  ! stores new attributes
  RDDH_DESCR(NFIELDDDH+1)%CNAME=CDNAME
  RDDH_DESCR(NFIELDDDH+1)%CMOD=CDMOD
  RDDH_DESCR(NFIELDDDH+1)%CFLUX=CDFLUX
  RDDH_DESCR(NFIELDDDH+1)%LKDDH=LDDH

  ALLOCATE(RDDH_DESCR(NFIELDDDH+1)%RFIELD(IDIM1,IDIM2))
  RDDH_DESCR(NFIELDDDH+1)%RFIELD=> RDDH_FIELD(1:IDIM1,1:IDIM2,IDIM3+1)
  RDDH_DESCR(NFIELDDDH+1)%NFIELDIND=NFIELDDDH+1

ELSE

!-------------------------------


  IF (CDFLUX =='F') THEN
    RDDH_FIELD(1:SIZE(PMAT,1),1:SIZE(RDDH_FIELD,2),NFIELDDDH+1)=PMAT(1:SIZE(PMAT,1),1:SIZE(RDDH_FIELD,2))
  ELSE
    RDDH_FIELD(1:SIZE(PMAT,1),1:SIZE(RDDH_FIELD,2)-1,NFIELDDDH+1)=PMAT(1:SIZE(PMAT,1),1:SIZE(RDDH_FIELD,2)-1)
    RDDH_FIELD(1:SIZE(PMAT,1),SIZE(RDDH_FIELD,2),NFIELDDDH+1)=0.0
  ENDIF

  RDDH_DESCR(NFIELDDDH+1)%RFIELD => RDDH_FIELD(:,:,NFIELDDDH+1)

ENDIF


NFIELDDDH=NFIELDDDH+1

!IF (LHOOK) CALL DR_HOOK('ADD_FIELD_3D',1,ZHOOK_HANDLE)
IF (LHOOK) CALL DR_HOOK('DDH_MIX:ADD_FIELD_3D',1,ZHOOK_HANDLE)

END SUBROUTINE ADD_FIELD_3D




!-*-------------------------------------------------------------------------------

!-----------------------------------------------
! SUBROUTINE ADD_FIELD_2D
! allocates on the fly the flexible DDH structure
! during first time step and then stores the
! velues of retrieved fields into it
! PMAT= Value of field to retrieve
! CDNAME= Name of field
! CDFLUX='F','V','S' if flux,variable or tendency
! LDINST=TRUE if instantaneous field
! LDDH=TRUE if field to be stored into ddh files


SUBROUTINE ADD_FIELD_2D(YDLDDH,PMAT,CDNAME,CDFLUX,CDMOD,LDINST,LDDH)
USE YOMLDDH , ONLY : TLDDH
TYPE(TLDDH)                 ,INTENT(IN)   :: YDLDDH
REAL(KIND=JPRB),DIMENSION(:),INTENT(IN)   :: PMAT
CHARACTER(LEN=*)            ,INTENT(IN)   :: CDNAME
CHARACTER(LEN=1)            ,INTENT(IN)   :: CDFLUX
CHARACTER(LEN=3)            ,INTENT(IN)   :: CDMOD
LOGICAL                     ,INTENT(IN)   :: LDINST,LDDH
REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE::ZSDDH_FIELD_BKP
INTEGER(KIND=JPIM)::IDIM1,IDIM2,IDIM3,ITEST
TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: YSLDDH_DESCR_BKP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('DDH_MIX:ADD_FIELD_2D',0,ZHOOK_HANDLE)

IF (YDLDDH%LDDH_OMP) THEN
  !WRITE (NULOUT,*) 'WARNING: ADD_FIELD_2D disabled for LDDH_OMP=.TRUE.'
  IF (LHOOK) CALL DR_HOOK('DDH_MIX:ADD_FIELD_2D',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IF(LALLOCSURF) THEN
  ITEST=NSURFDDH+1
  IF(NSURFDDH==NTOTSURF) ITEST=NSURFDDH+1-NTOTSURF
  IF(CDNAME /= RDDHSURF_DESCR(ITEST)%CNAME) THEN
    WRITE(NULOUT,*) 'PB NOMS ADD_FIELD2D ',CDNAME,' vs ',RDDHSURF_DESCR(ITEST)%CNAME
    CALL ABOR1('ADD_FIELD_2D')
  ENDIF
ENDIF

IF(.NOT.LALLOCSURF) THEN
  IF (CNAMESURFINI == CDNAME) THEN
    LALLOCSURF=.TRUE.
  ENDIF
ENDIF

IF (LALLOCSURF .AND. NSURFDDH == NTOTSURF) THEN
  NSURFDDH=NSURFDDH-NTOTSURF
  RDDHSURF_FIELD(:,:,:)=0.0
ENDIF

IDIM1=NPROMADDH

IDIM2=1          ! arrays are stored on surface level


IF (.NOT.LALLOCSURF) THEN ! (re-)allocation of arrays
  NTOTSURF=NTOTSURF+1 !updates total number of fields
  IF(CDFLUX=='V') NTOTSVAR=NTOTSVAR+1
  IF(CDFLUX=='S') NTOTSVFS=NTOTSVFS+1

  IF (NSURFDDH == 0)  THEN

    ALLOCATE(RDDHSURF_DESCR(1))
    ALLOCATE(RDDHSURF_FIELD(IDIM1,IDIM2,1))
    RDDHSURF_FIELD(:,:,:)=0.
    CNAMESURFINI=CDNAME ! stores name of first field
    IDIM3=0

  ELSE !NSURFDDH=0

    ! add a field to rddh_descr
    ALLOCATE(YSLDDH_DESCR_BKP(NSURFDDH))
    YSLDDH_DESCR_BKP=RDDHSURF_DESCR
    DEALLOCATE(RDDHSURF_DESCR)
    ALLOCATE(RDDHSURF_DESCR(NSURFDDH+1))
    RDDHSURF_DESCR(1:NSURFDDH)=YSLDDH_DESCR_BKP(1:NSURFDDH)
    DEALLOCATE(YSLDDH_DESCR_BKP)

    IDIM3=SIZE(RDDHSURF_FIELD,3)
    ! reallocation of rddh_field in order to increase dimension
    ALLOCATE(ZSDDH_FIELD_BKP(IDIM1,IDIM2,IDIM3))
    ZSDDH_FIELD_BKP=RDDHSURF_FIELD ! temporary bkp of rddh_field
    DEALLOCATE(RDDHSURF_FIELD)
    ALLOCATE(RDDHSURF_FIELD(IDIM1,IDIM2,IDIM3+1)) !one field added
    RDDHSURF_FIELD(:,:,:)=0.0
    RDDHSURF_FIELD(1:IDIM1,1:IDIM2,1:IDIM3)=ZSDDH_FIELD_BKP(1:IDIM1,1:IDIM2,1:IDIM3)
    DEALLOCATE(ZSDDH_FIELD_BKP)

  ENDIF !NSURFDDH=0

  RDDHSURF_FIELD(1:SIZE(PMAT,1),1,IDIM3+1)=PMAT(1:SIZE(PMAT,1)) !stores new flux

  ! stores new attributes
  RDDHSURF_DESCR(NSURFDDH+1)%CNAME=CDNAME
  RDDHSURF_DESCR(NSURFDDH+1)%CMOD=CDMOD
  RDDHSURF_DESCR(NSURFDDH+1)%CFLUX=CDFLUX
  RDDHSURF_DESCR(NSURFDDH+1)%LKDDH=LDDH

  ALLOCATE(RDDHSURF_DESCR(NSURFDDH+1)%RFIELD(IDIM1,IDIM2))
  RDDHSURF_DESCR(NSURFDDH+1)%RFIELD=> RDDHSURF_FIELD(1:IDIM1,1:IDIM2,IDIM3+1)
  RDDHSURF_DESCR(NSURFDDH+1)%NFIELDIND=NSURFDDH+1

ELSE

!-------------------------------

  RDDHSURF_FIELD(1:SIZE(PMAT,1),1,NSURFDDH+1)=PMAT(1:SIZE(PMAT,1))
  RDDHSURF_DESCR(NSURFDDH+1)%RFIELD => RDDHSURF_FIELD(:,:,NSURFDDH+1)

ENDIF

NSURFDDH=NSURFDDH+1

IF (LHOOK) CALL DR_HOOK('DDH_MIX:ADD_FIELD_2D',1,ZHOOK_HANDLE)

END SUBROUTINE ADD_FIELD_2D



!-*------------------------------------------
! INITIALISE ET DESALLOUE LES STRUCTURES DDH

SUBROUTINE RESET_DDHFLEX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDH_MIX:RESET_DDHFLEX',0,ZHOOK_HANDLE)

!IF (LHOOK) CALL DR_HOOK('RESET_DDHFLEX',0,ZHOOK_HANDLE)
IF(ALLOCATED(RDDH_DESCR)) DEALLOCATE(RDDH_DESCR)
IF(ALLOCATED(RDDH_FIELD)) DEALLOCATE(RDDH_FIELD)
NFIELDDDH=0
NTOTFIELD=0
NTOTVAR=0
NTOTFSL=0
CNAMEINI='XXXXXXXXX'
LALLOC=.FALSE.

IF(ALLOCATED(RDDHSURF_DESCR)) DEALLOCATE(RDDHSURF_DESCR)
IF(ALLOCATED(RDDHSURF_FIELD)) DEALLOCATE(RDDHSURF_FIELD)
NSURFDDH=0
NTOTSURF=0
NTOTSVAR=0
NTOTSVFS=0
CNAMESURFINI='XXXXXXXXX'
LALLOCSURF=.FALSE.

!IF (LHOOK) CALL DR_HOOK('RESET_DDH_FLEX',1,ZHOOK_HANDLE)
IF (LHOOK) CALL DR_HOOK('DDH_MIX:RESET_DDHFLEX',1,ZHOOK_HANDLE)

END SUBROUTINE RESET_DDHFLEX


!------------------------------------------------
! creates the mask needed for computation on DDH ARRAYS
! mask has same dimension than pddhcv in ppsydh
! namely (nflevg+1)*(ndhcvsu+ndhvv)
! or (nflevg+1)*(ntotfield+ntotvar) within new struct.

SUBROUTINE MK_MASK
INTEGER(KIND=JPIM)::JFIELD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDH_MIX:MK_MASK',0,ZHOOK_HANDLE)

ALLOCATE(LMASKDDH((NTOTFIELD+NTOTVAR)*(NFLEVGDDH+1)))
LMASKDDH(:)=.FALSE.
DO JFIELD=1,NTOTFIELD
  IF(RDDH_DESCR(JFIELD)%CFLUX=='V') THEN
    LMASKDDH((JFIELD-1)*(NFLEVGDDH+1)+1:JFIELD*(NFLEVGDDH+1))=.TRUE.
  ENDIF
ENDDO

! variables at initial time
LMASKDDH(NTOTFIELD*(NFLEVGDDH+1)+1:(NTOTFIELD+NTOTVAR)*(NFLEVGDDH+1))=.TRUE.

IF (LHOOK) CALL DR_HOOK('DDH_MIX:MK_MASK',1,ZHOOK_HANDLE)

END SUBROUTINE MK_MASK



!------------------------------------------------
! creates the mask needed for computation on DDH ARRAYS
! mask has same dimension than pddhcv in ppsydh namely
! (ndhcssu+ndhvs)or(ntotsurf+ntotsvar) within new struct.

SUBROUTINE MK_MASKSURF
INTEGER(KIND=JPIM)::JFIELD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDH_MIX:MK_MASKSURF',0,ZHOOK_HANDLE)

! 2D surface fields
ALLOCATE(LMASKDDHVAS(NTOTSURF))
ALLOCATE(LMASKDDHVFS(NTOTSURF))
ALLOCATE(LMASKDDHSURF(NTOTSURF+NTOTSVAR))
LMASKDDHVAS(:)=.FALSE.
LMASKDDHVFS(:)=.FALSE.
LMASKDDHSURF(:)=.FALSE.
DO JFIELD=1,NTOTSURF
  IF(RDDHSURF_DESCR(JFIELD)%CFLUX=='V') THEN
    LMASKDDHVAS(JFIELD)=.TRUE.
  ENDIF
  IF(RDDHSURF_DESCR(JFIELD)%CFLUX=='S') THEN
    LMASKDDHVFS(JFIELD)=.TRUE.
  ENDIF
  LMASKDDHSURF(JFIELD)=LMASKDDHVAS(JFIELD).OR.LMASKDDHVFS(JFIELD)
ENDDO

! surface variables at initial time
LMASKDDHSURF(NTOTSURF+1:NTOTSURF+NTOTSVAR)=.TRUE.

IF (LHOOK) CALL DR_HOOK('DDH_MIX:MK_MASKSURF',1,ZHOOK_HANDLE)

END SUBROUTINE MK_MASKSURF



!*----------------------------
!stores nlon,nflevg into module

SUBROUTINE SUDDHFLEX(KPROMA,KLEV)

INTEGER(KIND=JPIM), INTENT(IN) ::KPROMA,KLEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDH_MIX:SUDDHFLEX',0,ZHOOK_HANDLE)

NFLEVGDDH=KLEV
NPROMADDH=KPROMA

IF (LHOOK) CALL DR_HOOK('DDH_MIX:SUDDHFLEX',1,ZHOOK_HANDLE)

END SUBROUTINE SUDDHFLEX
!------------------------------------------------


!====================================================================
!====================================================================
!====================================================================

SUBROUTINE NEW_ADD_FIELD_3D(YDMDDH,PVAL,CDNAME,YDDDH,CDTYPE)

  ! routine to add a single field to the DDH superstructure

USE YOMMDDH , ONLY : TMDDH
  IMPLICIT NONE

 TYPE(TMDDH)      , INTENT(IN)    :: YDMDDH
  REAL(KIND=JPRB) , INTENT(IN)    :: PVAL(:,:)          ! numerical values
  CHARACTER(LEN=*), INTENT(IN)    :: CDNAME             ! name of the field
  TYPE(TYP_DDH)   , INTENT(INOUT) :: YDDDH              ! DDH superstructure
  CHARACTER(LEN=1), OPTIONAL, INTENT(IN)   :: CDTYPE   ! type of the field (flux, tendency or variable)
  TYPE(TYP_FIELD3D), ALLOCATABLE   :: YLFIELDS_BKP(:)  ! array of fields, necessary for reallocation
  TYPE(TYP_FIELD3D), POINTER       :: YLFIELD          ! pointer to the added field

  INTEGER(KIND=JPIM) :: JROF, JFLD

  CHARACTER(LEN=1) :: CLTYPE
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:NEW_ADD_FIELD_3D',0,ZHOOK_HANDLE)

  YDDDH%NFIELDS3D=YDDDH%NFIELDS3D+1
  JFLD=YDDDH%NFIELDS3D

  ! reallocate metadata if necessary
  IF (JFLD>SIZE(YDDDH%RVAL3D,3)) THEN
    !WRITE (NULOUT,*) 'Warning: reallocating metadata'
    ALLOCATE(YLFIELDS_BKP(JFLD-1)); YLFIELDS_BKP=YDDDH%YFIELD3D; DEALLOCATE(YDDDH%YFIELD3D)
    ALLOCATE(YDDDH%YFIELD3D(JFLD)); YDDDH%YFIELD3D(1:JFLD-1)=YLFIELDS_BKP;DEALLOCATE(YLFIELDS_BKP)
  ENDIF
  YLFIELD=>YDDDH%YFIELD3D(JFLD)

  ! allocate or point numerical data
  IF (JFLD>SIZE(YDDDH%RVAL3D,3)) THEN
    ! starting at 1 instead of 0 because pointers always start at 1 if
    ! pointer bound mapping is not supported by the compiler
    ALLOCATE(YLFIELD%RVAL(1:YDDDH%NLEV+1,YDMDDH%NDHIDH))
  ELSE
    YLFIELD%RVAL=>YDDDH%RVAL3D(:,:,JFLD)
  ENDIF

  ! set numerical data (average)
  YLFIELD%RVAL(:,:)=0.
  DO JROF=YDDDH%KST,YDDDH%KEND
    YLFIELD%RVAL(1:SIZE(PVAL,2),YDDDH%NDDHI(JROF))=YLFIELD%RVAL(1:SIZE(PVAL,2),YDDDH%NDDHI(JROF))&
     & +PVAL(JROF,:)*YDDDH%WEIGHT(JROF)
  ENDDO

  IF (PRESENT(CDTYPE)) THEN
    CLTYPE=CDTYPE
  ELSE
    CLTYPE=CDNAME(1:1)
  ENDIF

  ! set metadata
  YLFIELD%CNAME=CDNAME
  YLFIELD%CTYPE=CLTYPE

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:NEW_ADD_FIELD_3D',1,ZHOOK_HANDLE)

END SUBROUTINE NEW_ADD_FIELD_3D


!====================================================================

SUBROUTINE NEW_ADD_FIELD_2D(YDMDDH,PVAL,CDNAME,YDDDH,CDTYPE)

  ! routine to add a single field to the DDH superstructure

USE YOMMDDH , ONLY : TMDDH
  IMPLICIT NONE

 TYPE(TMDDH)      , INTENT(IN)    :: YDMDDH
  REAL(KIND=JPRB) , INTENT(IN)    :: PVAL(:)            ! numerical values
  CHARACTER(LEN=*), INTENT(IN)    :: CDNAME             ! name of the field
  TYPE(TYP_DDH)   , INTENT(INOUT) :: YDDDH              ! DDH superstructure
  CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: CDTYPE    ! type of the field (flux, tendency or variable)
  TYPE(TYP_FIELD2D), ALLOCATABLE   :: YLFIELDS_BKP(:)  ! array of fields, necessary for reallocation
  TYPE(TYP_FIELD2D), POINTER       :: YLFIELD          ! pointer to the added field

  INTEGER(KIND=JPIM) :: JROF, JFLD

  CHARACTER(LEN=1) :: CLTYPE
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:NEW_ADD_FIELD_2D',0,ZHOOK_HANDLE)

  YDDDH%NFIELDS2D=YDDDH%NFIELDS2D+1
  JFLD=YDDDH%NFIELDS2D

  ! reallocate metadata if necessary
  IF (JFLD>SIZE(YDDDH%RVAL2D,2)) THEN
    !WRITE (NULOUT,*) 'Warning: reallocating metadata'
    ALLOCATE(YLFIELDS_BKP(JFLD-1)); YLFIELDS_BKP=YDDDH%YFIELD2D; DEALLOCATE(YDDDH%YFIELD2D)
    ALLOCATE(YDDDH%YFIELD2D(JFLD)); YDDDH%YFIELD2D(1:JFLD-1)=YLFIELDS_BKP;DEALLOCATE(YLFIELDS_BKP)
  ENDIF
  YLFIELD=>YDDDH%YFIELD2D(JFLD)

  ! allocate or point numerical data
  IF (JFLD>SIZE(YDDDH%RVAL2D,2)) THEN
    ALLOCATE(YLFIELD%RVAL(YDMDDH%NDHIDH))
  ELSE
    YLFIELD%RVAL=>YDDDH%RVAL2D(:,JFLD)
  ENDIF


  ! set numerical data (average)
  YLFIELD%RVAL(:)=0.
  DO JROF=YDDDH%KST,YDDDH%KEND
    YLFIELD%RVAL(YDDDH%NDDHI(JROF))=YLFIELD%RVAL(YDDDH%NDDHI(JROF))&
     & +PVAL(JROF)*YDDDH%WEIGHT(JROF)
  ENDDO

  IF (PRESENT(CDTYPE)) THEN
    CLTYPE=CDTYPE
  ELSE
    CLTYPE=CDNAME(1:1)
  ENDIF

  ! set metadata
  YLFIELD%CNAME=CDNAME
  YLFIELD%CTYPE=CLTYPE

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:NEW_ADD_FIELD_2D',1,ZHOOK_HANDLE)

END SUBROUTINE NEW_ADD_FIELD_2D

!====================================================================

SUBROUTINE SETDDH(YDMDDH,YDDDH,KST,KEND,KPROMA,KLEV,KITER,PXW,KDDHI,PAUX3D,PAUX2D,PAUXSM)

  ! initialize the DDH superstructure

USE YOMMDDH , ONLY : TMDDH
  IMPLICIT NONE

  TYPE(TMDDH)       , INTENT(IN)         :: YDMDDH
  TYPE(TYP_DDH)     , INTENT(INOUT)      :: YDDDH           ! DDH superstructure
  INTEGER(KIND=JPIM), INTENT(IN)         :: KST             ! index of first point in this NPROMA block
  INTEGER(KIND=JPIM), INTENT(IN)         :: KEND            ! index of first point in this NPROMA block
  INTEGER(KIND=JPIM), INTENT(IN)         :: KPROMA          ! horizontal dimension
  INTEGER(KIND=JPIM), INTENT(IN)         :: KLEV            ! vertical dimension
  INTEGER(KIND=JPIM), INTENT(IN)         :: KITER           ! Nb of Corrector iterations
  REAL(KIND=JPRB),    INTENT(IN), TARGET :: PXW(KPROMA)     ! averaging weights
  INTEGER(KIND=JPIM), INTENT(IN), TARGET :: KDDHI(:)        ! DDH domain index
  REAL(KIND=JPRB),    INTENT(IN), TARGET :: PAUX3D(:,:,:)   ! automatically allocated array for 3d fields
  REAL(KIND=JPRB),    INTENT(IN), TARGET :: PAUX2D(:,:)     ! automatically allocated array for 2d fields
  REAL(KIND=JPRB),    INTENT(IN), TARGET :: PAUXSM(:,:,:,:) ! automatically allocated array for variables

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:SETDDH',0,ZHOOK_HANDLE)
 
  IF (KITER == 0) THEN
  
   ! automatically allocated auxiliary array
   YDDDH%RVAL3D=>PAUX3D
   YDDDH%RVAL2D=>PAUX2D
   YDDDH%RVARSM=>PAUXSM

   ! start counting at zero ...
   YDDDH%NFIELDS3D=0
   YDDDH%NFIELDS2D=0

   ! ... but with an offset
   YDDDH%NFIELDS3D_OFFSET=YDMDDH%NFIELDS3D_OFFSET
   YDDDH%NFIELDS2D_OFFSET=YDMDDH%NFIELDS2D_OFFSET

   ! set weights
   YDDDH%WEIGHT=>PXW

   ! set domains
   YDDDH%NDDHI=>KDDHI
 
   ! set dimensions
   YDDDH%NLEV=KLEV
   YDDDH%NPROMA=KPROMA
   YDDDH%KST=KST
   YDDDH%KEND=KEND

   ! allocate space for metadata
   ALLOCATE(YDDDH%YFIELD3D(SIZE(PAUX3D,3)))
   ALLOCATE(YDDDH%YFIELD2D(SIZE(PAUX2D,2)))
   !ALLOCATE(YDDDH%YVARMULT(SIZE(PAUXSM,4)))
   ALLOCATE(YDDDH%YVARMULT(13))
 
 
   YDDDH%YVARMULT(1)%CNAME  = "UU"
   YDDDH%YVARMULT(2)%CNAME  = "VV"
   YDDDH%YVARMULT(3)%CNAME  = "WW"
   YDDDH%YVARMULT(4)%CNAME  = "CT"
   YDDDH%YVARMULT(5)%CNAME  = "KK"
   YDDDH%YVARMULT(6)%CNAME  = "QV"
   YDDDH%YVARMULT(7)%CNAME  = "QL"
   YDDDH%YVARMULT(8)%CNAME  = "QR"
   YDDDH%YVARMULT(9)%CNAME  = "QI"
   YDDDH%YVARMULT(10)%CNAME = "QS"
   YDDDH%YVARMULT(11)%CNAME = "QG"
   YDDDH%YVARMULT(12)%CNAME = "QH"
   YDDDH%YVARMULT(13)%CNAME = "SV"

   ! check dimensions of PAUX
   IF ( SIZE(PAUX3D,1) /= KLEV+1 ) CALL ABOR1("dimension of PAUX should be 0:KLEV")

  ELSE
    
    ! set dimensions
    YDDDH%NLEV=KLEV
    YDDDH%NPROMA=KPROMA
    YDDDH%KST=KST
    YDDDH%KEND=KEND      
    YDDDH%NFIELDS3D= NTOTFIELD  
    YDDDH%NFIELDS2D= NTOTSURF
    YDDDH%NFIELDS3D_OFFSET= 0
    YDDDH%NFIELDS2D_OFFSET= 0

  END IF   

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:SETDDH',1,ZHOOK_HANDLE)

END SUBROUTINE SETDDH

!====================================================================

SUBROUTINE CLEANDDH(YDTDDH,YDDDH,KSTG,KITER)

  ! clean the DDH superstructure

USE YOMMDDH , ONLY : TMDDH
USE YOMTDDH , ONLY : TTDDH
  IMPLICIT NONE

  TYPE(TTDDH)   , INTENT(INOUT) :: YDTDDH
  TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH      ! DDH superstructure
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTG     ! global index of first point
  INTEGER(KIND=JPIM), INTENT(IN) :: KITER    ! Nbe of corrector iterations

  INTEGER(KIND=JPIM)           :: JFLD

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:CLEANDDH',0,ZHOOK_HANDLE)

  IF (KITER == 0) THEN  
  ! store constant data (names and types)
   IF (KSTG==1) THEN
     DO JFLD=1,YDDDH%NFIELDS3D
       YDTDDH%CFLDNAMES3D(YDDDH%NFIELDS3D_OFFSET+JFLD)=YDDDH%YFIELD3D(JFLD)%CNAME
       YDTDDH%CFLDTYPES3D(YDDDH%NFIELDS3D_OFFSET+JFLD)=YDDDH%YFIELD3D(JFLD)%CTYPE
     ENDDO
     DO JFLD=1,YDDDH%NFIELDS2D
       YDTDDH%CFLDNAMES2D(YDDDH%NFIELDS2D_OFFSET+JFLD)=YDDDH%YFIELD2D(JFLD)%CNAME
       YDTDDH%CFLDTYPES2D(YDDDH%NFIELDS2D_OFFSET+JFLD)=YDDDH%YFIELD2D(JFLD)%CTYPE
     ENDDO
   ENDIF
 
   ! deallocate manually allocated fields
   DO JFLD=SIZE(YDDDH%RVAL3D,3)+1,YDDDH%NFIELDS3D
     DEALLOCATE(YDDDH%YFIELD3D(JFLD)%RVAL)
   ENDDO
   DO JFLD=SIZE(YDDDH%RVAL2D,2)+1,YDDDH%NFIELDS2D
     DEALLOCATE(YDDDH%YFIELD2D(JFLD)%RVAL)
   ENDDO

   ! deallocate array of fields
   DEALLOCATE(YDDDH%YFIELD3D)
   DEALLOCATE(YDDDH%YFIELD2D)
   DEALLOCATE(YDDDH%YVARMULT)

   ! nullify pointers for safety
   YDDDH%RVAL3D=>NULL()
   YDDDH%YFIELD3D=>NULL()
   YDDDH%RVAL2D=>NULL()
   YDDDH%YFIELD2D=>NULL()
   YDDDH%RVARSM=>NULL()
   YDDDH%YVARMULT=>NULL()

   ! update number of fields from this routine
   YDDDH%NFIELDS3D_AUTO=MAX(YDDDH%NFIELDS3D,YDDDH%NFIELDS3D_AUTO)
   YDDDH%NFIELDS2D_AUTO=MAX(YDDDH%NFIELDS2D,YDDDH%NFIELDS2D_AUTO)

  ELSE

   YDDDH%NFIELDS3D= NTOTFIELD  
   YDDDH%NFIELDS2D= NTOTSURF
   YDDDH%NFIELDS3D_OFFSET= 0
   YDDDH%NFIELDS2D_OFFSET= 0

  END IF 

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:CLEANDDH',1,ZHOOK_HANDLE)

END SUBROUTINE CLEANDDH

!====================================================================

SUBROUTINE STOREDDH(YDMDDH,YDTDDH,YDDDH)

  ! store the constant features of the DDH superstructure in global variables
  ! WARNING: this routine should not be called inside an OpenMP parallel region!

USE YOMMDDH , ONLY : TMDDH
USE YOMTDDH , ONLY : TTDDH
  IMPLICIT NONE

 TYPE(TMDDH)   , INTENT(INOUT) :: YDMDDH
 TYPE(TTDDH)   , INTENT(IN)    :: YDTDDH
  TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH      ! DDH superstructure
  INTEGER(KIND=JPIM)           :: JFLD, IPLST, IPLEND

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:STOREDDH',0,ZHOOK_HANDLE)

  ! store total number of fields and variables
  NTOTFIELD=YDDDH%NFIELDS3D_OFFSET+YDDDH%NFIELDS3D
  NTOTVAR=COUNT(YDTDDH%CFLDTYPES3D(1:NTOTFIELD)=='V')
  NTOTSURF=YDDDH%NFIELDS2D_OFFSET+YDDDH%NFIELDS2D
  NTOTSVAR=COUNT(YDTDDH%CFLDTYPES2D(1:NTOTSURF)=='V')
  NTOTSVFS=COUNT(YDTDDH%CFLDTYPES2D(1:NTOTSURF)=='S')

  ! create variable mask
  IF (ALLOCATED(LMASKDDH)) DEALLOCATE(LMASKDDH)
  ALLOCATE(LMASKDDH((NTOTFIELD+NTOTVAR)*(YDDDH%NLEV+1)))
  LMASKDDH=.FALSE.
  DO JFLD=1,NTOTFIELD
    IPLST =(1+YDDDH%NLEV)*(JFLD-1)+1
    IPLEND=(1+YDDDH%NLEV)*(JFLD)
    IF (YDTDDH%CFLDTYPES3D(JFLD)=='V') THEN
      LMASKDDH(IPLST:IPLEND)=.TRUE.
    ENDIF
  ENDDO
  !LMASKDDH(IPLEND+1:(NTOTFIELD+NTOTVAR)*(YDDDH%NLEV+1))=.TRUE.
  IPLST=IPLEND+1
  IPLEND=(NTOTFIELD+NTOTVAR)*(YDDDH%NLEV+1)
  DO JFLD=IPLST,IPLEND
    LMASKDDH(JFLD)=.TRUE.
  ENDDO
 
  ! create surface variable mask
  IF (ALLOCATED(LMASKDDHSURF)) DEALLOCATE(LMASKDDHSURF)
  IF (ALLOCATED(LMASKDDHVAS)) DEALLOCATE(LMASKDDHVAS)
  IF (ALLOCATED(LMASKDDHVFS)) DEALLOCATE(LMASKDDHVFS)
  ALLOCATE(LMASKDDHSURF((NTOTSURF+NTOTSVAR)))
  ALLOCATE(LMASKDDHVAS((NTOTSURF)))
  ALLOCATE(LMASKDDHVFS((NTOTSURF)))
  LMASKDDHSURF=.FALSE.
  LMASKDDHVAS=.FALSE.
  LMASKDDHVFS=.FALSE.
  DO JFLD=1,NTOTSURF
    IF (YDTDDH%CFLDTYPES2D(JFLD)=='V') THEN
      LMASKDDHVAS(JFLD)=.TRUE.
    ELSEIF (YDTDDH%CFLDTYPES2D(JFLD)=='S') THEN
      LMASKDDHVFS(JFLD)=.TRUE.
    ENDIF
    LMASKDDHSURF(JFLD)= LMASKDDHVAS(JFLD).OR.LMASKDDHVFS(JFLD)
    !IF (YDTDDH%CFLDTYPES2D(JFLD)=='V') THEN
    !  LMASKDDHSURF(JFLD)=.TRUE.
    !ENDIF
  ENDDO
  LMASKDDHSURF(NTOTSURF+1:NTOTSURF+NTOTSVAR)=.TRUE.

  ! set offset value before the second loop
  YDMDDH%NFIELDS3D_OFFSET=NTOTFIELD
  YDMDDH%NFIELDS2D_OFFSET=NTOTSURF

  ! store number of fields that should be allocated automatically
  YDMDDH%NFIELDS3D_AUTO=YDDDH%NFIELDS3D_AUTO
  YDMDDH%NFIELDS2D_AUTO=YDDDH%NFIELDS2D_AUTO

  IF (LHOOK) CALL DR_HOOK('DDH_MIX:STOREDDH',1,ZHOOK_HANDLE)

END SUBROUTINE STOREDDH

!------------------------------------------------

END MODULE DDH_MIX

