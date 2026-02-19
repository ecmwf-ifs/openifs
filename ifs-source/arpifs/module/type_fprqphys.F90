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

MODULE TYPE_FPRQPHYS

! Purpose :
! -------
!    To define the internal type "TRQFP" for the control of 
!    post-processed physical surface fields :
!      - %ICOD : Internal code number of the field
!      - %IVEC : Information about tensor aspect of the field in spectral space :
!        -j : vector 2nd momentum corresponding to the field nr j
!         0 : "scalar" 
!         1 : vector 1st momentum
!      - %IDOM : the number of subdomains for each field
!      - %IDMP : the subdomains pointers for each field

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
!    Fullpos technical & users guide.

! Author :
! ------
!    Ryad El Khatib *METEO-FRANCE*

! Modifications :
! -------------
! Original : 2003-03-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARFPOS, ONLY : JPOSPHY
USE TYPE_FPDSPHYS, ONLY : FPDSPHY
USE YOMFP4L, ONLY : TRQFP

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC ALLOCATE_FPRQPHY, DEALLOCATE_FPRQPHY, INQUIRE_FPRQPHY

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE ALLOCATE_FPRQPHY(YDGFP_PHYDS,YDTYPE,KDOM,KFIELDS,KFPCOD,KULOUT,LDALLOPR,KFPDOM,KDMPTR)

TYPE(FPDSPHY),      INTENT(IN)    :: YDGFP_PHYDS(JPOSPHY)
TYPE(TRQFP), INTENT(INOUT) :: YDTYPE
INTEGER(KIND=JPIM), INTENT(IN)    :: KDOM
INTEGER(KIND=JPIM), INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN)    :: KFPCOD(KFIELDS)
INTEGER(KIND=JPIM), INTENT(IN)    :: KULOUT
LOGICAL,            INTENT(IN)    :: LDALLOPR
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KFPDOM(KFIELDS)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KDMPTR(KDOM,KFIELDS)

INTEGER(KIND=JPIM) :: JF, JV, JU, JD, IUCOD, IVCOD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! KDOM : number of subdomains
! KFIELDS : number of fields
! KFPCOD : internal fields codes
! KULOUT : output unit number
! LDALLOPR : .TRUE. to report about allocation
! KFPDOM : number of domains for each field
! KDMPTR : domains pointers for each field

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQPHYS:ALLOCATE_FPRQPHY',0,ZHOOK_HANDLE)

YDTYPE%NFIELDG=KFIELDS
ALLOCATE(YDTYPE%ICOD(KFIELDS))
ALLOCATE(YDTYPE%IGRIB(KFIELDS))
ALLOCATE(YDTYPE%LLSURF(KFIELDS))
ALLOCATE(YDTYPE%CLNAME(KFIELDS))
ALLOCATE(YDTYPE%CLPREF(KFIELDS))
ALLOCATE(YDTYPE%IVEC(KFIELDS))
ALLOCATE(YDTYPE%IDOM(KFIELDS))
ALLOCATE(YDTYPE%IDMP(KDOM,KFIELDS))

IF (LDALLOPR) THEN
  CALL INQUIRE_FPRQPHY(YDTYPE,'YDTYPE',KULOUT)
ENDIF

DO JF=1,KFIELDS
  YDTYPE%ICOD(JF)=KFPCOD(JF)
  YDTYPE%CLPREF(JF)=YDGFP_PHYDS(KFPCOD(JF))%CLNAME(1:4)
  YDTYPE%CLNAME(JF)=YDGFP_PHYDS(KFPCOD(JF))%CLNAME(5:)      
  YDTYPE%IGRIB(JF)=YDGFP_PHYDS(KFPCOD(JF))%IGRIB
  YDTYPE%LLSURF(JF)=YDGFP_PHYDS(KFPCOD(JF))%LLSRF
  YDTYPE%IVEC(JF)=0
ENDDO
DO JV=1,KFIELDS
  IVCOD=YDTYPE%ICOD(JV)
  IF (YDGFP_PHYDS(IVCOD)%IORDR < 0) THEN
    ! V => find U :
    DO JU=1,KFIELDS
      IUCOD=YDTYPE%ICOD(JU)
      IF (YDGFP_PHYDS(IUCOD)%IORDR > 0) THEN
        IF (YDGFP_PHYDS(IUCOD)%CLPAIR == YDGFP_PHYDS(IVCOD)%CLPAIR) THEN
          YDTYPE%IVEC(JU)=1
          YDTYPE%IVEC(JV)=-JU
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDDO
IF (PRESENT(KFPDOM)) THEN
  ! Select domains
  DO JF=1,KFIELDS
    YDTYPE%IDOM(JF)=KFPDOM(JF)
  ENDDO
  IF (PRESENT(KDMPTR)) THEN
    DO JF=1,KFIELDS
      DO JD=1,KFPDOM(JF)
        YDTYPE%IDMP(JD,JF)=KDMPTR(JD,JF)
      ENDDO
    ENDDO
  ELSE
    ! First domains
    DO JF=1,KFIELDS
      DO JD=1,KFPDOM(JF)
        YDTYPE%IDMP(JD,JF)=JD
      ENDDO
    ENDDO
  ENDIF
ELSE
  ! All domains
  DO JF=1,KFIELDS
    YDTYPE%IDOM(JF)=KDOM
  ENDDO
  DO JF=1,KFIELDS
    DO JD=1,KDOM
      YDTYPE%IDMP(JD,JF)=JD
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQPHYS:ALLOCATE_FPRQPHY',1,ZHOOK_HANDLE)

END SUBROUTINE ALLOCATE_FPRQPHY

!-----------------------------------------------------------------------------

SUBROUTINE DEALLOCATE_FPRQPHY(YDTYPE)

TYPE(TRQFP), INTENT(INOUT) :: YDTYPE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQPHYS:DEALLOCATE_FPRQPHY',0,ZHOOK_HANDLE)

YDTYPE%NFIELDG=0
DEALLOCATE(YDTYPE%ICOD)
DEALLOCATE(YDTYPE%IGRIB)
DEALLOCATE(YDTYPE%LLSURF)
DEALLOCATE(YDTYPE%CLNAME)
DEALLOCATE(YDTYPE%CLPREF)
DEALLOCATE(YDTYPE%IVEC)
DEALLOCATE(YDTYPE%IDOM)
DEALLOCATE(YDTYPE%IDMP)

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQPHYS:DEALLOCATE_FPRQPHY',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLOCATE_FPRQPHY

!-----------------------------------------------------------------------------

SUBROUTINE INQUIRE_FPRQPHY(YDTYPE,CDNAME,KULOUT)

! YDTYPE  : structure variable
! CDNAME : structure name
! KULOUT : output logical unit number

TYPE(TRQFP), INTENT(IN) :: YDTYPE
CHARACTER(LEN=*),   INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),          INTENT(IN) :: KULOUT

! CLFMT  : format for output

CHARACTER(LEN=36) :: CLFMT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQPHYS:INQUIRE_FPRQPHY',0,ZHOOK_HANDLE)

CLFMT='(1X,''ARRAY '',A,A7,'' ALLOCATED '',8I8)'

WRITE(KULOUT,CLFMT) CDNAME,'%ICOD  ',SIZE(YDTYPE%ICOD),SHAPE(YDTYPE%ICOD)
WRITE(KULOUT,CLFMT) CDNAME,'%IVEC  ',SIZE(YDTYPE%IVEC),SHAPE(YDTYPE%IVEC)
WRITE(KULOUT,CLFMT) CDNAME,'%IDOM  ',SIZE(YDTYPE%IDOM),SHAPE(YDTYPE%IDOM)
WRITE(KULOUT,CLFMT) CDNAME,'%IDMP  ',SIZE(YDTYPE%IDMP),SHAPE(YDTYPE%IDMP)

IF (LHOOK) CALL DR_HOOK('TYPE_FPRQPHYS:INQUIRE_FPRQPHY',1,ZHOOK_HANDLE)

END SUBROUTINE INQUIRE_FPRQPHY

!-----------------------------------------------------------------------------

END MODULE TYPE_FPRQPHYS


