! (C) Copyright 1989- Meteo-France.

SUBROUTINE CPVPOSPR(KFPCONF,YDAFN,YDQTYPE)

!**** *CPVPOSPR*  - COMPUTE FIELD POINTERS FOR VERTICAL POST-PROCESSING

!     PURPOSE.
!     --------
!        Compute the dimensions and the pointers of the post-processing fields so that
!        the fields are stored one after the others in the arrays
!        - GT1 (spectrally fitted fields) 
!        - GAUX (gridpoint fields).
!        Fields not requested are pointing at adress 0
!         To fit the spectral transforms ordering, fields are ordered in GT1 as follows : 
!          all "U", then all "V", then all scalars

!**   INTERFACE.
!     ----------
!       *CALL* *CPVPOSPR*

!        EXPLICIT ARGUMENTS
!        --------------------
!           INPUT  :
!           KFPCONF : configuration of the post-processing :
!                     0 : vertical interpolation only (<CFPFMT='MODEL'>)
!                     1 : gridpoint post-processing, possibly with spectral filters (<NFPOS=1>)
!                     2 : gridpoint/spectral post-processing (spectral outputs possible) (<NFPOS=2>)

!        IMPLICIT ARGUMENTS
!        -----------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!    R. El Khatib : 01-03-28 compute YDQTYPE%IGT1 in the proper order to fit the spectral transforms
!    R. El Khatib : 01-04-04 YDQTYPE%IGT1(0) after inverse spectral transform
!    R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!    R. El Khatib : 03-04-17 Fullpos improvemnts
!    M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARFPOS  , ONLY : JPOSDYN
USE YOMMP0   , ONLY : MYSETV, NPRTRV 
USE YOMAFN   , ONLY : TAFN
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KFPCONF
TYPE(TAFN),        INTENT(IN)    :: YDAFN
TYPE(TYPE_FPRQDYN),INTENT(INOUT) :: YDQTYPE

INTEGER(KIND=JPIM) :: IVPTR(NPRTRV)   ! Pointer of fields for each V-set in output array
INTEGER(KIND=JPIM) :: IVCOUNT(NPRTRV) ! Number of fields for each V-set in spectral transform
INTEGER(KIND=JPIM) :: ILOC(1)
INTEGER(KIND=JPIM) :: J, INC, IJ, JL, ISPREAD, IVSET, ISDIR, ISINV, IMOMDIR, IMOMDIRG, IMOMINV, IMOMINVG, I, JJ
INTEGER(KIND=JPIM) :: IPTRGPX  ! last pointer for GAUX
INTEGER(KIND=JPIM) :: IVOR     ! last pointer of vorticity fields in GT1
INTEGER(KIND=JPIM) :: IDIV     ! last pointer of divergence fields in GT1
INTEGER(KIND=JPIM) :: ISCA1    ! last pointer of scalar fields in GT1
INTEGER(KIND=JPIM) :: IUMOM    ! last pointer of U momenta in GT0
INTEGER(KIND=JPIM) :: IVMOM    ! last pointer of V momenta in GT0
INTEGER(KIND=JPIM) :: ISCA0    ! last pointer of scalar fields in GT0

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPVPOSPR',0,ZHOOK_HANDLE)
ASSOCIATE(TFP=>YDAFN%TFP, TFP_DYNDS=>YDAFN%TFP_DYNDS)

!*       3.    CONTROL REQUESTS
!              ----------------

!        3.1 Tensor aspects

YDQTYPE%ILED(:)=0
! Vorticities, divergence :
YDQTYPE%ILED(TFP%ABS%ICOD)=2
YDQTYPE%ILED(TFP%VOR%ICOD)=2
YDQTYPE%ILED(TFP%DIV%ICOD)=3
! Psi and Khi are computed from U and V in spectral space :
YDQTYPE%ILED(TFP%KHI%ICOD)=1
YDQTYPE%ILED(TFP%PSI%ICOD)=-TFP%KHI%ICOD
! Surface Moisture Convergence is computed as DIV(q.Ucls,q.Vcls) ;
YDQTYPE%ILED(TFP%SMC%ICOD)=3
YDQTYPE%ILED(TFP%ASMC%ICOD)=3
YDQTYPE%ILED(TFP%VSMC%ICOD)=3

! Vectors :
DO J=1,JPOSDYN
  IF (TFP_DYNDS(J)%IORDR < 0) THEN
    ! V => find U :
    DO JJ=1,JPOSDYN
      IF (TFP_DYNDS(JJ)%CLPAIR == TFP_DYNDS(J)%CLPAIR .AND. JJ /= J) THEN
        YDQTYPE%ILED(JJ)=1
        YDQTYPE%ILED(J)=-JJ
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO

IF (KFPCONF /= 1) THEN
  ! Vectors will be considered as individual "scalars" if spectrally fitted because nothing else follows 
  DO J=1,JPOSDYN
    IF (ABS(TFP_DYNDS(J)%IORDR) == 1) THEN
      YDQTYPE%ILED(J)=YDQTYPE%ILED(J)*(1-MIN(YDQTYPE%ISF(J),1))
    ENDIF
  ENDDO
ENDIF


!        3.3   Number of fields needed in gridpoint space :

YDQTYPE%ISKP(:)=1

YDQTYPE%ISKP(TFP%ABS%ICOD)=1+MIN(YDQTYPE%ISF(TFP%ABS%ICOD),1)
YDQTYPE%ISKP(TFP%VOR%ICOD)=1+MIN(YDQTYPE%ISF(TFP%VOR%ICOD),1)
YDQTYPE%ISKP(TFP%DIV%ICOD)=1+MIN(YDQTYPE%ISF(TFP%DIV%ICOD),1)

! Surface Moisture Convergence is computed as DIV(q.Ucls,q.Vcls) ;
YDQTYPE%ISKP(TFP%SMC%ICOD)=2
YDQTYPE%ISKP(TFP%ASMC%ICOD)=2
YDQTYPE%ISKP(TFP%VSMC%ICOD)=2



!*       2. COMPUTE DIMENSIONS

DO J=1,JPOSDYN
  IF (.NOT.YDQTYPE%LL(J)) THEN
    ! Nullify dimensions for non-requested fields and count requested ones
    YDQTYPE%ISF(J)=0
    YDQTYPE%ISKP(J)=0
    YDQTYPE%ILED(J)=0
    YDQTYPE%ILEV(J)=0
  ENDIF
ENDDO

!        2.1 Global dimensions for fitted/non-fitted fields in gp space

YDQTYPE%NFPDYNB = 0
YDQTYPE%NFPAUXB = 0
DO J=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  YDQTYPE%NFPDYNB = YDQTYPE%NFPDYNB + YDQTYPE%ILEV(IJ)
  YDQTYPE%NFPAUXB = YDQTYPE%NFPAUXB + (1-MIN(YDQTYPE%ISF(IJ),1))*YDQTYPE%ILEV(IJ)
ENDDO

!        3.7 "B-level" distribution (balanced here for direct transforms)

YDQTYPE%ISET(:,:)=1
IF (NPRTRV > 1) THEN
  IVCOUNT(:)=0
! Make 2 scans to improve balance : first those which have more that 1 field : 
  DO J=1,YDQTYPE%NFPOSDYN
    IJ=YDQTYPE%NFPTRDYN(J)
    IF (YDQTYPE%ISF(IJ) >= 1 .AND. YDQTYPE%ISKP(IJ) > 1) THEN
      DO JL=1,YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ILED(IJ) /= 1) THEN
!         Any field but U-vector momenta in direct transforms :
          ILOC=MINLOC(IVCOUNT)
          IVSET=ILOC(1)
          YDQTYPE%ISET(JL,IJ)=IVSET
          IVCOUNT(IVSET)=IVCOUNT(IVSET)+YDQTYPE%ISKP(IJ)
          IF (YDQTYPE%ILED(IJ) < 0) THEN
!           V-vector momenta => put counterpart U on the same set
            I=-YDQTYPE%ILED(IJ)
            YDQTYPE%ISET(JL,I)=IVSET
            IVCOUNT(IVSET)=IVCOUNT(IVSET)+YDQTYPE%ISKP(I)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  DO J=1,YDQTYPE%NFPOSDYN
    IJ=YDQTYPE%NFPTRDYN(J)
    IF (YDQTYPE%ISF(IJ) >= 1 .AND. YDQTYPE%ISKP(IJ) == 1) THEN
      DO JL=1,YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ILED(IJ) /= 1) THEN
!         Any field but U-vector momenta in direct transforms :
          ILOC=MINLOC(IVCOUNT)
          IVSET=ILOC(1)
          YDQTYPE%ISET(JL,IJ)=IVSET
          IVCOUNT(IVSET)=IVCOUNT(IVSET)+YDQTYPE%ISKP(IJ)
          IF (YDQTYPE%ILED(IJ) < 0) THEN
!           V-vector momenta => put counterpart U on the same set
            I=-YDQTYPE%ILED(IJ)
            YDQTYPE%ISET(JL,I)=IVSET
            IVCOUNT(IVSET)=IVCOUNT(IVSET)+YDQTYPE%ISKP(I)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDIF

!        2.2 Global/Local dimensions

YDQTYPE%NFPGT1  = 0
YDQTYPE%NFPGT0B = 0
YDQTYPE%NFPSCAG = 0
IMOMDIRG        = 0
YDQTYPE%NFPISCAG= 0
IMOMINVG        = 0
YDQTYPE%NFPSPD  = 0
YDQTYPE%NFPSPB  = 0
YDQTYPE%NFPSCA  = 0
IMOMDIR         = 0
YDQTYPE%NFPISCA = 0
IMOMINV         = 0
YDQTYPE%NFPUVMN = 0
YDQTYPE%ILOC(:,:)=0
IVPTR(:)=0
YDQTYPE%ISPD(:,:)=0
DO J=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  ! YDQTYPE%ISPD : Count how many arrays fields are spread over
  ISPREAD=MAX(0,YDQTYPE%ISF(IJ)-2)
  DO JL=1,YDQTYPE%ILEV(IJ)
    YDQTYPE%ISPD(JL,IJ)=ISPREAD*(YDQTYPE%IDOM(JL,IJ))+(1-ISPREAD)
  ENDDO
  IF (YDQTYPE%ISF(IJ) >= 1) THEN
    ISDIR=MAX(0,(1-ABS(YDQTYPE%ILED(IJ))))
    IF (YDQTYPE%ILED(IJ)==1.OR.YDQTYPE%ILED(IJ)<0) THEN
      ISINV=0
    ELSE
      ISINV=1
    ENDIF
    DO JL=1,YDQTYPE%ILEV(IJ)
      INC      = YDQTYPE%ISPD(JL,IJ)
      YDQTYPE%NFPGT1     = YDQTYPE%NFPGT1     + YDQTYPE%ISKP(IJ)
      YDQTYPE%NFPSCAG  = YDQTYPE%NFPSCAG  + YDQTYPE%ISKP(IJ)*ISDIR
      IMOMDIRG = IMOMDIRG + YDQTYPE%ISKP(IJ)*(1-ISDIR)
      IF (YDQTYPE%IFIT /= 2) YDQTYPE%NFPGT0B     = YDQTYPE%NFPGT0B     + INC
      YDQTYPE%NFPISCAG  = YDQTYPE%NFPISCAG  + INC*ISINV
      IMOMINVG = IMOMINVG + INC*(1-ISINV)
      IVSET=YDQTYPE%ISET(JL,IJ)
      IVPTR(IVSET)=IVPTR(IVSET)+1
      YDQTYPE%ILOC(JL,IJ)=IVPTR(IVSET)
      IF (IVSET==MYSETV) THEN
        ! V-set local dimensions
        YDQTYPE%NFPSPD    = YDQTYPE%NFPSPD    + MAX(0,YDQTYPE%ISF(IJ)-2)
        YDQTYPE%NFPSCA  = YDQTYPE%NFPSCA  + YDQTYPE%ISKP(IJ)*ISDIR
        IMOMDIR = IMOMDIR + YDQTYPE%ISKP(IJ)*(1-ISDIR)
        YDQTYPE%NFPSPB    = YDQTYPE%NFPSPB    + INC
        YDQTYPE%NFPISCA  = YDQTYPE%NFPISCA  + INC*ISINV
        IMOMINV = IMOMINV + INC*(1-ISINV)
        YDQTYPE%NFPUVMN   = YDQTYPE%NFPUVMN   + INC*MAX(0,-SIGN(1,YDQTYPE%ILED(IJ)))
      ENDIF
    ENDDO
  ENDIF
ENDDO
IF (MOD(IMOMDIRG,2) /= 0) THEN
  CALL ABOR1('CPVPOSPR : INTERNAL ERROR ON DIR VECTOR FIELDS => SEE SUVPOS')
ELSEIF (MOD(IMOMINVG,2) /= 0) THEN
  CALL ABOR1('CPVPOSPR : INTERNAL ERROR ON INV VECTOR FIELDS => SEE SUVPOS')
ELSE
  YDQTYPE%NFPVECG = IMOMDIRG/2
  YDQTYPE%NFPIVECG = IMOMINVG/2
ENDIF
IF (MOD(IMOMDIR,2) /= 0) THEN
  CALL ABOR1('CPVPOSPR : INTERNAL ERROR ON LOCAL DIR VECTOR FIELDS => SEE SUVPOS')
ELSEIF (MOD(IMOMINV,2) /= 0) THEN
  CALL ABOR1('CPVPOSPR : INTERNAL ERROR ON LOCAL INV VECTOR FIELDS => SEE SUVPOS')
ELSE
  YDQTYPE%NFPVEC = IMOMDIR/2
  YDQTYPE%NFPIVEC = IMOMINV/2
ENDIF


!*       3.    COMPUTE POINTERS
!              ----------------

! Compute pointers for GT1/GAUX (splitting the flow gp data against spectrally fitted data)

IVOR=1
IDIV=YDQTYPE%NFPVECG+1
ISCA1=2*YDQTYPE%NFPVECG+1
IUMOM=1
IVMOM=YDQTYPE%NFPIVECG+1
ISCA0=2*YDQTYPE%NFPIVECG+1
YDQTYPE%IGT1(:,:)=0
IPTRGPX=1
YDQTYPE%IGPX(:)=0
DO J=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  IF (YDQTYPE%ISF(IJ) >= 1) THEN
    ! Pointers % transforms : 
    INC=YDQTYPE%ILEV(IJ)
    ISPREAD=0
    IF (YDQTYPE%IFIT /= 2) THEN
      ! No "stretching effect" on the target grids
      DO  JL=1, INC
        ISPREAD=ISPREAD+YDQTYPE%ISPD(JL,IJ)
      ENDDO
    ENDIF
    SELECT CASE (YDQTYPE%ILED(IJ))
    CASE (0)
      ! Scalar
      YDQTYPE%IGT1(IJ,1)=ISCA1
      ISCA1=ISCA1+INC
      YDQTYPE%IGT1(IJ,0)=ISCA0
      ISCA0=ISCA0+ISPREAD
    CASE (:-1)
      ! V from (U,V)<=>(Vor,Div) : set Vor then Div for direct transforms
      YDQTYPE%IGT1(-YDQTYPE%ILED(IJ),1)=IVOR
      IVOR=IVOR+INC
      YDQTYPE%IGT1(IJ,1)=IDIV
      IDIV=IDIV+INC
      ! set U then V for inverse transforms
      YDQTYPE%IGT1(-YDQTYPE%ILED(IJ),0)=IUMOM
      IUMOM=IUMOM+ISPREAD
      YDQTYPE%IGT1(IJ,0)=IVMOM
      IVMOM=IVMOM+ISPREAD
    CASE (2,3)
      !  (U,V) for individual Vor or Div : set Vor then Div
      YDQTYPE%IGT1(IJ,1)=IVOR
      IVOR=IVOR+INC
      YDQTYPE%IGT1(IJ,2)=IDIV
      IDIV=IDIV+INC
      ! set as scalar in inverse transforms because only 1 field selected
      YDQTYPE%IGT1(IJ,0)=ISCA0
      ISCA0=ISCA0+ISPREAD
    END SELECT
  ELSE
    YDQTYPE%IGPX(IJ)= IPTRGPX
    IPTRGPX         = IPTRGPX + YDQTYPE%ILEV(IJ)
  ENDIF
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPVPOSPR',1,ZHOOK_HANDLE)
END SUBROUTINE CPVPOSPR
