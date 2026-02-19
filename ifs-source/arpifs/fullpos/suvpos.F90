! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUVPOS(YDQTYPE,YD3D,YD2D,YDAFN,KFPCONF,CDVER,CDHOR,LDSURF,LDLOWPASS,LDSP,LDSORT)

!**** *SUVPOS*  - INITIALIZE COMMON YOMVPOS

!     PURPOSE.
!     --------
!        Initialize common YOMVPOS which contains control variables
!        for vertical/horizontal post-processing
!        "B-level" parallelization is manage as follows :
!        the attributed V-set to a given level of a given field is "the one that is the less used" ;
!        Exceptions for vector fields : U/V corresponding to the same level are attributed the same
!        V-set

!**   INTERFACE.
!     ----------
!       *CALL* *SUVPOS*

!        EXPLICIT ARGUMENTS
!        --------------------
!           KFPCONF : configuration of the post-processing :
!                     0 : vertical interpolation only (<CFPFMT='MODEL'>)
!                     1 : gridpoint post-processing, possibly with spectral filters (<NFPOS=1>)
!                     2 : gridpoint/spectral post-processing (spectral outputs possible) (<NFPOS=2>)
!        CDVER : configuration of vertical interpolations
!        CDHOR : configuration of horizontal interpolations
!        LDSURF : .TRUE. to compute surface fields in this scan
!        LDLOWPASS  : .T. if any low-pass filter active
!        LDSP : enable/disable spectral output representation
!        LDSORT : sort output fields following the user request (T) or use the internal post-processing order (F)

!          OUTPUT
!     ==== 3D VARIABLES ===

!        IMPLICIT ARGUMENTS
!        --------------------
!          YOMAFN (input)

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
!      R. El Khatib : 01-05-07 control HUn,HUx
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvments
!      R. El Khatib : 03-11-13 Fix on YDQTYPE%ISPD
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib: 19-Jul-2012 LFIT* => NFIT* + NFITS, NFITH, NSPFIL*
!      R. El Khatib: 22-Aug-2012 Replace tests on cfpfmt
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!      P. Marquet   23-May-2014 Allow filtering of (u,v) on PV level
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib : 07-Jul-2014 remove exception on filtering U,V,Psi,Khi for PV surfaces
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El Khatib 26-Nov-2014 Bugfix for isothermic levels when Vor or Div are in the request list
!      R. El Khatib 07-Mar-2018 Optional arguments LDSP and LDSORT
!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARFPOS  , ONLY : JPOSDYN
USE YOMAFN   , ONLY : TAFN
USE YOMLUN   , ONLY : NULOUT
USE YOMFPC   , ONLY : LTRACEFP
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN
USE YOM4FPOS , ONLY : TRQ3FP, TRQ2FP

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(INOUT) :: YDQTYPE
TYPE(TRQ3FP),      INTENT(IN) :: YD3D
TYPE(TRQ2FP),      INTENT(IN) :: YD2D
TYPE (TAFN),       INTENT(IN)    :: YDAFN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPCONF
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDVER
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDHOR
LOGICAL           ,INTENT(IN)    :: LDSURF
LOGICAL           ,INTENT(IN)    :: LDLOWPASS
LOGICAL           ,INTENT(IN), OPTIONAL :: LDSP
LOGICAL           ,INTENT(IN), OPTIONAL :: LDSORT

INTEGER(KIND=JPIM) :: JDOM, J, JLEV, ICOD, I, IB, JJ, IU, IV, ISFPPLEV

! ISFPPLEV : Filter indicator of the current pp. level

LOGICAL :: LLSUPL  ! a supplementary condition
LOGICAL :: LLSP, LLSORT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1.    PRESET REQUEST
!              --------------

IF (LHOOK) CALL DR_HOOK('SUVPOS',0,ZHOOK_HANDLE)
ASSOCIATE(TFP=>YDAFN%TFP, TFP_DYNDS=>YDAFN%TFP_DYNDS)

LLSP=.TRUE.
IF (PRESENT(LDSP)) THEN
  LLSP=LDSP
ENDIF

LLSORT=.FALSE.
IF (PRESENT(LDSORT)) THEN
  LLSORT=LDSORT
ENDIF

YDQTYPE%LL(:)=.FALSE.

!*       2.    SETUP REQUEST
!              -------------

!*       2.2  Request on 2D fields

  IF (LDSURF) THEN
    DO J=1,YD2D%NPPFIELDG
      ICOD=YD2D%ICOD(J)
!     The next line enable to remove automatically the uncomputable fields
!     We give a chance to surface fields if .not.lhpos because any pp level
!     should be suitable
      YDQTYPE%LL(ICOD)=(INDEX(TFP_DYNDS(ICOD)%CLNIL,CDVER) == 0).OR.(CDHOR == '0')
!       The next block is for debugging purpose
        IF (LTRACEFP .AND. .FALSE.) THEN
          IF (YDQTYPE%LL(ICOD)) THEN
            WRITE(NULOUT,'(''SUVPOS : POST-PROCESSING ON FIELD '',A16, &
             & '' FOR CONFIGURATION '',A1)') TFP_DYNDS(ICOD)%CLNAME(1:16),CDVER
          ELSE 
            WRITE(NULOUT,'(''SUVPOS : MISSED FIELD '',A16, '' FOR CONFIGURATION '',A1,L2,L2)') &
             & TFP_DYNDS(ICOD)%CLNAME(1:16),CDVER,(INDEX(TFP_DYNDS(ICOD)%CLNIL,CDVER) == 0)
          ENDIF
        ENDIF
      IF (TFP_DYNDS(ICOD)%LLGP .OR. .NOT.LLSP) THEN
        IF (ICOD == TFP%ASMC%ICOD .OR. ICOD == TFP%SMC%ICOD .OR. ICOD == TFP%VSMC%ICOD) THEN
          IF (LLSP) THEN
            ! Force spectral fit
            YDQTYPE%ISF(ICOD)=MAX(TFP_DYNDS(ICOD)%ISF,1)
          ELSE
            ! there is a potential confusion between spectral fit and spectral output here :
            ! spectral fit is needed to build these fields 
            ! better remove them for now :
            YDQTYPE%LL(ICOD)=.FALSE.
            WRITE(NULOUT,'(''SUVPOS : MISSED FIELD '',A16, '' FOR CONFIGURATION '',A1,L2,L2, '' LLSP = '',L2)') &
             & TFP_DYNDS(ICOD)%CLNAME(1:16),CDVER,(INDEX(TFP_DYNDS(ICOD)%CLNIL,CDVER) == 0),LLSP
          ENDIF
        ELSE
!         No spectral fit
          YDQTYPE%ISF(ICOD)=0
        ENDIF
      ELSE
        YDQTYPE%ISF(ICOD)=TFP_DYNDS(ICOD)%ISF
        IF (YDQTYPE%ISF(ICOD) == 3 .AND. .NOT.LDLOWPASS) THEN
          YDQTYPE%ISF(ICOD)=1
        ENDIF
!       Monitoring of the spectral fit on 2D fields :
        IF (LTRACEFP) THEN
          IF (YDQTYPE%LL(ICOD) .AND. YDQTYPE%ISF(ICOD) > 1) THEN
            WRITE(NULOUT,'(''SUVPOS : SPECTRAL FILTER ACTIVE (ISF='',I1,'') ON FIELD '',A16, &
             & '' FOR CONFIGURATION '',A1)') &
             & YDQTYPE%ISF(ICOD),TFP_DYNDS(ICOD)%CLNAME(1:16),CDVER
          ENDIF
        ENDIF
      ENDIF
      YDQTYPE%ILEV(ICOD)=1
      YDQTYPE%ILVP(1,ICOD)=1
      YDQTYPE%IDOM(1,ICOD)=YD2D%IDOM(J)
    ENDDO
    DO J=1,YD2D%NPPFIELDG
      ICOD=YD2D%ICOD(J)
      DO JDOM=1,YD2D%IDOM(J)
        YDQTYPE%IDMP(JDOM,1,ICOD)=YD2D%IDMP(JDOM,J)
      ENDDO
    ENDDO
  ENDIF

!*       2.3  Request on 3D fields 

  YDQTYPE%IFIT=YD3D%NFIT

! Low-pass and gaussian filters are exclusive
  IF (YD3D%NSPFIL== 3 .AND. .NOT.LDLOWPASS) THEN
    ISFPPLEV=1
  ELSE
    ISFPPLEV=YD3D%NSPFIL
  ENDIF

  DO J=1,YD3D%NPPFIELDG
    ICOD=YD3D%ICOD(J)
!   The next line enables to remove automatically the uncomputable fields
    YDQTYPE%LL(ICOD)=(INDEX(TFP_DYNDS(ICOD)%CLNIL,CDVER) == 0)
    IF (TFP_DYNDS(ICOD)%LLGP .OR. .NOT.LLSP) THEN
!     No spectral fit
      YDQTYPE%ISF(ICOD)=0
    ELSE
      SELECT CASE (YDQTYPE%IFIT)
      CASE (1,2)
        YDQTYPE%ISF(ICOD)=MAX(TFP_DYNDS(ICOD)%ISF,ISFPPLEV)
        IF (YDQTYPE%ISF(ICOD)==3 .AND. .NOT.LDLOWPASS) THEN
          YDQTYPE%ISF(ICOD)=1
        ENDIF
      CASE DEFAULT
!       The filter is switched off
        YDQTYPE%ISF(ICOD)=0
      END SELECT
      IF (LTRACEFP) THEN
        IF (YDQTYPE%LL(ICOD) .AND. YDQTYPE%ISF(ICOD) > 1) THEN
          WRITE(NULOUT,'(''SUVPOS : SPECTRAL FILTER ACTIVE (ISF='',I1,'') ON FIELD '',A16, &
           & '' FOR CONFIGURATION '',A1)') &
           & YDQTYPE%ISF(ICOD),TFP_DYNDS(ICOD)%CLNAME(1:16),CDVER
        ENDIF
      ENDIF
    ENDIF
    YDQTYPE%ILEV(ICOD)=YD3D%NPPLEVG(J)
  ENDDO
  DO J=1,YD3D%NPPFIELDG
    ICOD=YD3D%ICOD(J)
    DO JLEV=1,YDQTYPE%ILEV(ICOD)
      YDQTYPE%ILVP(JLEV,ICOD)=YD3D%ILEV(JLEV,J)
      YDQTYPE%IDOM(JLEV,ICOD)=YD3D%IDOM(JLEV,J)
    ENDDO
  ENDDO
  DO J=1,YD3D%NPPFIELDG
    ICOD=YD3D%ICOD(J)
    DO JLEV=1,YD3D%NPPLEVG(J)
      DO JDOM=1,YDQTYPE%IDOM(JLEV,ICOD)
        YDQTYPE%IDMP(JDOM,JLEV,ICOD)=YD3D%IDMP(JDOM,JLEV,J)
      ENDDO
    ENDDO
  ENDDO

!        3.4   Special difficulties related to spectral space :

  ! Khi, Psi can be computed only in spectral space :
  I=TFP%KHI%ICOD
  YDQTYPE%LL(I)=YDQTYPE%LL(I).AND.(YDQTYPE%IFIT >= 1)
  I=TFP%PSI%ICOD
  YDQTYPE%LL(I)=YDQTYPE%LL(I).AND.(YDQTYPE%IFIT >= 1)
  ! Vor and Div can be computed from U and V if spectral space is possible
  I=TFP%VOR%ICOD
  YDQTYPE%LL(I)=YDQTYPE%LL(I).AND. ((YDQTYPE%IFIT == 2).OR.(CDHOR /= 'I'))
  I=TFP%DIV%ICOD
  YDQTYPE%LL(I)=YDQTYPE%LL(I).AND. ((YDQTYPE%IFIT == 2).OR.(CDHOR /= 'I'))
  ! Omega vertical velocity may be computed if no horizontal interpolations on surface-dependent levels
  I=TFP%VV%ICOD
  YDQTYPE%LL(I)=YDQTYPE%LL(I).AND. (CDHOR /= 'I')

  IF (KFPCONF /= 1) THEN
  ! Absolute Vorticity as spectral field would need to add the Coriolis parameter
  ! just before writing out to file and thus extra spectral transforms.
  ! Better remove it for the time being :
  ! ?? Not very convincing ...
    I=TFP%ABS%ICOD
    YDQTYPE%LL(I)=YDQTYPE%LL(I).AND.(YDQTYPE%ISF(I) == 0)
  ! Computation of Extrema are out of object since cls fields are available (??)
    YDQTYPE%LL(TFP%TX%ICOD)=.FALSE.
    YDQTYPE%LL(TFP%TN%ICOD)=.FALSE.
    YDQTYPE%LL(TFP%HUX%ICOD)=.FALSE.
    YDQTYPE%LL(TFP%HUN%ICOD)=.FALSE.
  ! Moisture convergence as a gridpoint field would need
  ! an extra inverse spectral transform. Better remove it for the time being :
    I=TFP%SMC%ICOD
    YDQTYPE%LL(I)=YDQTYPE%LL(I).AND.(YDQTYPE%ISF(I) >= 1)
    I=TFP%ASMC%ICOD
    YDQTYPE%LL(I)=YDQTYPE%LL(I).AND.(YDQTYPE%ISF(I) >= 1)
    I=TFP%VSMC%ICOD
    YDQTYPE%LL(I)=YDQTYPE%LL(I).AND.(YDQTYPE%ISF(I) >= 1)
  ELSE
    ! Vectors : add missing components (for spectral fit or compass)
    DO J=1,JPOSDYN
      IF (TFP_DYNDS(J)%IORDR == 1) THEN
        ! U => find V :
        DO JJ=1,JPOSDYN
          IF (TFP_DYNDS(JJ)%CLPAIR == TFP_DYNDS(J)%CLPAIR .AND. JJ /= J) THEN
            IU=J
            IV=JJ
            LLSUPL=.FALSE.
            IF (YDQTYPE%LL(IU).AND..NOT.YDQTYPE%LL(IV)) THEN
              I=IU
              IB=IV
              LLSUPL=.TRUE.
            ELSEIF(YDQTYPE%LL(IV).AND..NOT.YDQTYPE%LL(IU)) THEN
              IB=IU
              I=IV
              LLSUPL=.TRUE.
            ENDIF
            IF (LLSUPL) THEN
              YDQTYPE%LL(IB)=.TRUE.
              YDQTYPE%ISF(IB)=YDQTYPE%ISF(I)
              YDQTYPE%ILEV(IB)=YDQTYPE%ILEV(I)
              DO JLEV=1,YDQTYPE%ILEV(I)
                YDQTYPE%ILVP(JLEV,IB)=YDQTYPE%ILVP(JLEV,I)
                YDQTYPE%IDOM(JLEV,IB)=0
              ENDDO
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (LLSORT) THEN
    I=0
    DO J=1,YD3D%NPPFIELDG
      ICOD=YD3D%ICOD(J)
      IF (YDQTYPE%LL(ICOD)) THEN
        I=I+1
        YDQTYPE%NFPTRDYN(I)=ICOD
      ENDIF
    ENDDO
    IF (LDSURF) THEN
      DO J=1,YD2D%NPPFIELDG
        ICOD=YD2D%ICOD(J)
        IF (YDQTYPE%LL(ICOD)) THEN
          I=I+1
          YDQTYPE%NFPTRDYN(I)=ICOD
        ENDIF
      ENDDO
    ENDIF
    YDQTYPE%NFPOSDYN=COUNT(YDQTYPE%LL(:).EQV..TRUE.)
    IF (YDQTYPE%NFPOSDYN /= I) THEN
      CALL ABOR1('SUVPOS: SORTING FIELDS HAS FAILED BECAUSE CERTAIN FIELDS WERE REJECTED')
    ENDIF
  ELSE
    DO J=1,JPOSDYN
      IF (YDQTYPE%LL(J)) THEN
        YDQTYPE%NFPOSDYN=YDQTYPE%NFPOSDYN+1
        YDQTYPE%NFPTRDYN(YDQTYPE%NFPOSDYN)=J
      ENDIF
    ENDDO
  ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVPOS',1,ZHOOK_HANDLE)
END SUBROUTINE SUVPOS
