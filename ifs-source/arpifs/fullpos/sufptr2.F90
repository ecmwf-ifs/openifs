! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPTR2(YDRQPHY,YDRQAUX,KFPINPHY,YDNAMFPSCI,YDAFN,LDPHYS,LDHPOS,LDFPOSHOR,YDRQDYN)

!**** *SUFPTR2*  - POST-PROCESSING SETUP POINTERS OF BUFFER #2 -

!     PURPOSE.
!     --------
!        To compute the number of fields in auxilary fullpos buffer and 
!        their descriptors

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPTR2*

!        EXPLICIT ARGUMENTS
!        ------------------
!        KFPINPHY : = 0 if no actual interpolations on physical fields
!        LDPHYS : physics is active in the model
!        LDHPOS : any horizontal post-processing
!        LDFPOSHOR : actual horizontal interpolations

!        IMPLICIT ARGUMENTS
!        ------------------
!        See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION ABOUT FULL-POS.
!        NB : the auxilary buffer will help for the post-processing of physical
!        fields and dynamic fields on surface-dependent levels (height, eta). 
!        The auxilary buffer must contain all the fields needed for itself !

!     EXTERNALS.
!     ----------
!        NONE.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      Original : 98-10-08 from SUFPSC2B

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-01-05 plug-in FPCICA
!      R. El Khatib : 01-04-23 MAQ2M -> relative moisture ; HUn,HUx
!      R. El Khatib : 01-06-14 Interpolated Ts+alt.corr for iso-0 & NFPLAKE
!      R. El Khatib : 02-09-04 Iso -10 Celsius
!      O.Spaniel    : 03-04-15 cleaning-a same named entity from modules
!      R. El Khatib : 03-04-17 Fullpos improvments
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      S. Ivatek-Sahdan : 04-11-30 pre29 dbg lsm not needed for O3 and aeros in e923
!      R. El Khatib : 09-May-2003 preparatory fields for SURFEX
!      F. Bouyssel : 17-Jul-2009 LFPCAPEX
!      Y. Seity    : 4-April-2010 bf for running prepsurfex in the same geometry
!      Y. Seity    : 4-Feb-2011 bf for activating Ts if proftemperature required
!      R. El Khatib 22-Aug-2012 cfpfmt=model replaced by .not.lhpos
!      Y. Seity    : 4-Sep-2012 bf for prepsurfex (PHI was over-written)
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!      R. El Khatib 01-Sep-2014 Generalize option LFPCLSTOGMV
!      T. Aspelien : 09-Sep-2013 Make LSM visible for SURFEX 
!      R. El Khatib 07-Apr-2015 Bugfix to allow the computation of CAPE under
!      certain specific condtions (no climatology, no interpolation).
!      R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!      R. El Khatib 17-Aug-2016 Better interoperability GRIB2 vs FA
!      R. El Khatib 09-Sep-2016 Allow NFPCLI=2 again for interoperability
!      R. Brozkova Jul-2018 Added convective temperature
!      R. El Khatib 16-Apr-2020 Bugfix for KFPINPHY=0
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARFPOS   , ONLY : JPOSPHY
USE YOMLUN    , ONLY : NULOUT
USE YOMCT0    , ONLY : LARPEGEF
USE YOMAFN    , ONLY : TAFN
USE YOMFPC       , ONLY : TNAMFPSCI, LTRACEFP, LALLOFP
USE YOMFP4L , ONLY : TRQFP, IFPSEARCH
USE YOM4FPOS, ONLY : TRQFPDYN
USE TYPE_FPRQPHYS , ONLY :  ALLOCATE_FPRQPHY

IMPLICIT NONE

TYPE(TRQFP), INTENT(IN) :: YDRQPHY
TYPE(TRQFP), INTENT(OUT) :: YDRQAUX
INTEGER(KIND=JPIM), INTENT(IN) :: KFPINPHY
TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
TYPE (TAFN),  INTENT(IN) :: YDAFN
LOGICAL      , INTENT(IN) :: LDPHYS
LOGICAL      , INTENT(IN) :: LDHPOS
LOGICAL      , INTENT(IN) :: LDFPOSHOR
TYPE(TRQFPDYN), INTENT(IN), OPTIONAL :: YDRQDYN

INTEGER(KIND=JPIM) :: JFLD, IPTR, IFPDOM

INTEGER(KIND=JPIM) :: IFPVT0(JPOSPHY) ! temporary list of auxilary fields codes

LOGICAL :: LLUPPERAIR ! upper air vertical post-processing on surface-dependent levels is expected
LOGICAL :: LLDYPP ! Surface-dependent vertical post-processing is expected
LOGICAL :: LLTS   ! Physical surface temperature requested
LOGICAL :: LLTCLS ! Physical PBL temperature requested
LOGICAL :: LLRHCLS! Physical PBL relative moisture requested
LOGICAL :: LLRHSRF! Surface relative moisture requested
LOGICAL :: LLORO  ! Interpolated orography requested for physical fields
LOGICAL :: LLPBL  ! Recomputation of PBL fields
LOGICAL :: LLCAPE ! Computation of cape
LOGICAL :: LLDY2D ! any dynamical surface-dependant 2D-fields
LOGICAL :: LLFPSX1! Fullpos to Surfex fields

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!     ------------------------------------------------------------------

!*       1.    PRELIMINAR CONTROLS
!              -------------------

IF (LHOOK) CALL DR_HOOK('SUFPTR2',0,ZHOOK_HANDLE)
ASSOCIATE(NFPLAKE=>YDNAMFPSCI%NFPLAKE, NFPCLI=>YDNAMFPSCI%NFPCLI, NFPSURFEX=>YDNAMFPSCI%NFPSURFEX, LFPCAPEX=>YDNAMFPSCI%LFPCAPEX, &
 & LFPCLSTOGMV=>YDNAMFPSCI%LFPCLSTOGMV, NFPMASK=>YDNAMFPSCI%NFPMASK, TFP=>YDAFN%TFP, GFP=>YDAFN%GFP, GFP_PHYDS=>YDAFN%GFP_PHYDS)


!*       2.  COMPUTE MAIN LOGICAL KEYS
!            -------------------------

LLCAPE=.FALSE.
IF (PRESENT(YDRQDYN)) THEN
  DO JFLD=1,YDRQDYN%Y2D%NPPFIELDG
    IF (YDRQDYN%Y2D%ICOD(JFLD)==TFP%CAPE%ICOD.OR. &
     &  YDRQDYN%Y2D%ICOD(JFLD)==TFP%CIEN%ICOD.OR. &
     &  YDRQDYN%Y2D%ICOD(JFLD)==TFP%TCVS%ICOD) THEN
      LLCAPE=.TRUE.
      EXIT
    ENDIF
  ENDDO
ENDIF

LLPBL=LLCAPE.AND.(.NOT.LFPCAPEX)
IF (PRESENT(YDRQDYN)) THEN
  DO JFLD=1,YDRQDYN%Y2D%NPPFIELDG
    IF (YDRQDYN%Y2D%ICOD(JFLD)==TFP%UCLS%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%VCLS%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%TCLS%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%QCLS%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%RCLS%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%TX%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%TN%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%HUX%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%HUN%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%UGST%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%VGST%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%FGST%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%FCLS%ICOD ) THEN  
      IF (LDPHYS) THEN
        IF (.NOT.LDHPOS .OR. NFPCLI >= 3) THEN
          LLPBL=.TRUE.
          EXIT
        ELSE
          CALL ABOR1('SUFPTR2 : PBL POST-PROCESSING NEEDS NFPCLI=3')
        ENDIF
      ELSE
        CALL ABOR1('SUFPTR2 : PBL POST-PROCESSING NEEDS LDPHYS=.TRUE.')
      ENDIF
    ENDIF
  ENDDO
ENDIF

LLDY2D=.FALSE.
IF (PRESENT(YDRQDYN)) THEN
  DO JFLD=1,YDRQDYN%Y2D%NPPFIELDG
    IF (YDRQDYN%Y2D%ICOD(JFLD)==TFP%SP%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%LNSP%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%WWS%ICOD  .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%UJET%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%VJET%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%PJET%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%TCAO%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%PCAO%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%HTPW%ICOD .OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%HTPW1%ICOD.OR.&
       & YDRQDYN%Y2D%ICOD(JFLD)==TFP%HTPW2%ICOD ) THEN  
      LLDY2D=.TRUE.
      EXIT
    ENDIF
  ENDDO
ENDIF

LLUPPERAIR=.FALSE.
IF (PRESENT(YDRQDYN)) THEN
  IF (YDRQDYN%YHO%NPPFIELDG > 0) THEN
    LLUPPERAIR=LLUPPERAIR .OR. (MAXVAL(YDRQDYN%YHO%NPPLEVG(:)) > 0)
  ENDIF
  IF (YDRQDYN%YIT%NPPFIELDG > 0) THEN
    LLUPPERAIR=LLUPPERAIR .OR. (MAXVAL(YDRQDYN%YIT%NPPLEVG(:)) > 0)
  ENDIF
  IF (YDRQDYN%YFL%NPPFIELDG > 0) THEN
    LLUPPERAIR=LLUPPERAIR .OR. (MAXVAL(YDRQDYN%YFL%NPPLEVG(:)) > 0)
  ENDIF
  IF (YDRQDYN%YES%NPPFIELDG > 0) THEN
    LLUPPERAIR=LLUPPERAIR .OR. (MAXVAL(YDRQDYN%YES%NPPLEVG(:)) > 0)
  ENDIF
ENDIF
LLDYPP = LDHPOS .AND. (LLUPPERAIR .OR. LLDY2D .OR. LLPBL .OR. LLCAPE)

LLFPSX1=(NFPSURFEX==1)

!*       2.  COMPUTE POINTERS, DESCRIPTORS & AND RESULTING DIMENSION
!            -------------------------------------------------------

IPTR=0

!        2.1 Orography

! Interpolated orography needed for surface-dependent vertical post-processing
IF (LLDYPP) THEN
  CALL SET_NEXT(GFP%SFIS%ICOD)
ELSEIF (LLFPSX1) THEN
  IF (NFPCLI == 0) THEN
    ! Or for Surfex if no climatology dataset present:
    CALL SET_NEXT(GFP%SFIS%ICOD)
  ELSE
    ! Output orography put as well in auxilary buffer for Surfex : makes an easy redirection
     CALL SET_NEXT(GFP%GFIS%ICOD)
  ENDIF
ENDIF

!        2.2 Surface physical fields

IF (LDPHYS) THEN

!        2.2.1 Climatology fields

! Interpolated land-sea mask needed for any correction on surface fields
! except surface geopotential, O3 profiles and aerosols
  IF (NFPCLI == 0 .AND. LDFPOSHOR) THEN
    DO JFLD=1,YDRQPHY%NFIELDG
      IF (YDRQPHY%ICOD(JFLD) /= GFP%SFIS%ICOD .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%GFIS%ICOD .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%O3A%ICOD  .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%O3B%ICOD  .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%O3C%ICOD  .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%ASEA%ICOD .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%ALAN%ICOD .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%ASOO%ICOD .AND. &
       & YDRQPHY%ICOD(JFLD) /= GFP%ADES%ICOD )  THEN
        CALL SET_NEXT(GFP%LSM%ICOD)
        EXIT
      ENDIF
    ENDDO
  ENDIF

! if no climatology, anisotropy and direction of topography are computed
! from the vector "anisotropy"
  IF (IFPSEARCH(YDRQPHY,GFP%DPAT%ICOD) /= 0 .OR. IFPSEARCH(YDRQPHY,GFP%ACOT%ICOD) /= 0) THEN
    IF (NFPCLI <= 1 .AND. LDFPOSHOR) THEN
      CALL SET_NEXT(GFP%SDOG%ICOD)
      CALL SET_NEXT(GFP%PADOU%ICOD)
      CALL SET_NEXT(GFP%PADOV%ICOD)
    ENDIF
  ENDIF

! if no climatology and actual interpolations, Leaf area index is needed for Interception content
  IF (NFPCLI <= 1 .AND. LDFPOSHOR) THEN
    IF (IFPSEARCH(YDRQPHY,GFP%IC%ICOD) /= 0) THEN
      CALL SET_NEXT(GFP%LAI%ICOD)
    ENDIF
  ENDIF

!        2.2.1 Pronostic fields

! Interpolated surface temperature is needed for surface-dependent 
! vertical interpolations or a lake option
  IF (LLDYPP.OR.(NFPLAKE == -1)) THEN
    CALL SET_NEXT(GFP%RDST%ICOD)
!   Supplementary fields if 2 masks:
    IF (NFPMASK >=2 ) THEN
      CALL SET_NEXT(GFP%LAN%ICOD)
    ENDIF
  ENDIF

! SST :
  IF (((LLDYPP .OR. NFPLAKE==-1) .AND. NFPMASK >=2) .OR. LLFPSX1) THEN
    CALL SET_NEXT(GFP%SST%ICOD)
  ENDIF

! Cls temperature needed for XFU extreme temperatures, gmv-to-cls option and sometimes CAPE
  LLTCLS=(LLCAPE.AND.LFPCAPEX).OR.LFPCLSTOGMV
  IF (NFPCLI > 0 .OR. LDFPOSHOR) THEN
    IF (IFPSEARCH(YDRQPHY,GFP%XX2T%ICOD) /= 0 .OR. IFPSEARCH(YDRQPHY,GFP%XN2T%ICOD) /= 0) THEN
      LLTCLS=.TRUE.
    ENDIF
  ENDIF
  IF (LLTCLS) THEN
    CALL SET_NEXT(GFP%X2T%ICOD)
  ENDIF

! Cls relative moisture needed for XFU extreme moistures, gmv-to-cls option, and CAPE sometimes
  LLRHCLS=(LLCAPE.AND.LFPCAPEX).OR.LFPCLSTOGMV
  IF (NFPCLI > 0 .OR. LDFPOSHOR) THEN
    IF (IFPSEARCH(YDRQPHY,GFP%XX2HU%ICOD) /= 0 .OR. IFPSEARCH(YDRQPHY,GFP%XN2HU%ICOD) /= 0) THEN
      LLRHCLS=.TRUE.
    ENDIF
  ENDIF
  IF (LLRHCLS) THEN
    CALL SET_NEXT(GFP%X2RH%ICOD)
  ENDIF

! Output surface temperature is needed for surface-dependent vertical interpolations, Tcls, snow depth, Tp and Surfex :
  LLTS=LLDYPP.OR.LLTCLS.OR.LLFPSX1
  IF (NFPCLI > 0 .OR. LDFPOSHOR) THEN
    IF (IFPSEARCH(YDRQPHY,GFP%SD%ICOD) /= 0 .OR. IFPSEARCH(YDRQPHY,GFP%DST%ICOD) /= 0 &
     & .OR. IFPSEARCH(YDRQPHY,GFP%X2T%ICOD) /= 0) THEN
      LLTS=.TRUE.
    ENDIF
  ENDIF
  IF (LLTS) THEN
    CALL SET_NEXT(GFP%ST%ICOD)
  ENDIF

! Interpolated orography may be needed for surface temperature :
  LLORO=.FALSE.
  IF (.NOT.LLDYPP .AND. (NFPCLI > 0 .OR. NFPLAKE == -1)) THEN  
    IF (LLTS .OR. IFPSEARCH(YDRQPHY,GFP%ST%ICOD) /= 0) THEN
      LLORO=.TRUE.
    ENDIF
  ENDIF
  IF (LLORO) THEN
    CALL SET_NEXT(GFP%SFIS%ICOD) 
  ENDIF

! Cls wind intensity is computed from the Cls wind components
  IF (IFPSEARCH(YDRQPHY,GFP%X10FF%ICOD) /= 0) THEN
    CALL SET_NEXT(GFP%X10U%ICOD)
    CALL SET_NEXT(GFP%X10V%ICOD)
  ENDIF

! Gusts intensity is computed from the Gust components
  IF (IFPSEARCH(YDRQPHY,GFP%XGUST%ICOD) /= 0) THEN
    CALL SET_NEXT(GFP%XUGST%ICOD)
    CALL SET_NEXT(GFP%XVGST%ICOD)
  ENDIF

! Surface relative moisture needed for cls relative moisture
  LLRHSRF=LLRHCLS
  IF (NFPCLI > 0 .OR. LDFPOSHOR) THEN
    IF (IFPSEARCH(YDRQPHY,GFP%X2RH%ICOD) /= 0) THEN
      LLRHSRF=(GFP%X2RH%IANO /= 0)
    ENDIF
  ENDIF
  IF (LLRHSRF) THEN
    CALL SET_NEXT(GFP%PSRHU%ICOD)
  ENDIF

! Advanced PBL post-processing or Surfex :
  IF (LLPBL .OR. LLFPSX1) THEN
!   Model snow depth :
    CALL SET_NEXT(GFP%SD%ICOD)
!   Model surface water content :
    CALL SET_NEXT(GFP%SSW%ICOD)
!   surface frost :
    CALL SET_NEXT(GFP%FSSW%ICOD)
  ENDIF

! Advanced PBL post-processing for U,V,T,Q & HU :
  IF (LLPBL) THEN
!   Interpolated dynamic surface roughness length :
    CALL SET_NEXT(GFP%IDZ0%ICOD)
!   Resistance to evapotranspiration : 
    CALL SET_NEXT(GFP%HV%ICOD)
!   Interpolated thermal surface roughness length :
    CALL SET_NEXT(GFP%ITZ0%ICOD)
  ENDIF
  IF (PRESENT(YDRQDYN)) THEN
    DO JFLD=1,YDRQDYN%Y2D%NPPFIELDG
      IF (YDRQDYN%Y2D%ICOD(JFLD)==TFP%TX%ICOD) THEN
        ! Increment to maxi temperature :
        CALL SET_NEXT(GFP%INCTX%ICOD)
        EXIT
      ENDIF
    ENDDO
  ENDIF
  IF (PRESENT(YDRQDYN)) THEN
    DO JFLD=1,YDRQDYN%Y2D%NPPFIELDG
      IF (YDRQDYN%Y2D%ICOD(JFLD)==TFP%TN%ICOD) THEN
        ! Increment to mini temperature :
        CALL SET_NEXT(GFP%INCTN%ICOD)
        EXIT
      ENDIF
    ENDDO
  ENDIF
  IF (PRESENT(YDRQDYN)) THEN
    DO JFLD=1,YDRQDYN%Y2D%NPPFIELDG
      IF (YDRQDYN%Y2D%ICOD(JFLD)==TFP%HUX%ICOD) THEN
        ! Increment to maxi relative moisture :
        CALL SET_NEXT(GFP%INCHX%ICOD)
        EXIT
      ENDIF
    ENDDO
  ENDIF
  IF (PRESENT(YDRQDYN)) THEN
    DO JFLD=1,YDRQDYN%Y2D%NPPFIELDG
      IF (YDRQDYN%Y2D%ICOD(JFLD)==TFP%HUN%ICOD) THEN
        ! Increment to mini temperature :
        CALL SET_NEXT(GFP%INCHN%ICOD)
        EXIT
      ENDIF
    ENDDO
  ENDIF

! Preparatory fields CAPE under certain conditions (no interpolations) :
  IF (LLCAPE.AND.KFPINPHY==0.AND.NFPCLI==0) THEN
    CALL SET_NEXT(GFP%IVEG%ICOD)
  ENDIF
! if no climatology, vegetation is needed for Interception content
  IF (NFPCLI <= 1 .AND. LDFPOSHOR) THEN
    IF (IFPSEARCH(YDRQPHY,GFP%IC%ICOD) /= 0 .OR. (LLCAPE.AND.KFPINPHY==0.AND.NFPCLI==0)) THEN
      CALL SET_NEXT(GFP%VEG%ICOD)
    ENDIF
  ENDIF

! Other preparatory fields for Surfex or CAPE under certain conditions (no interpolations) :
  IF (LLFPSX1.OR.(LLCAPE.AND.KFPINPHY==0.AND.NFPCLI==0)) THEN
    CALL SET_NEXT(GFP%LSM%ICOD)
    CALL SET_NEXT(GFP%ARG%ICOD)
    CALL SET_NEXT(GFP%SAB%ICOD)
    CALL SET_NEXT(GFP%D2%ICOD)
  ENDIF

! Other preparatory fields for Surfex :
  IF (LLFPSX1) THEN
    CALL SET_NEXT(GFP%DST%ICOD)
    CALL SET_NEXT(GFP%DSW%ICOD)
    CALL SET_NEXT(GFP%FDSW%ICOD)
  ENDIF

ENDIF

IF (IPTR > 0) THEN
  IFPDOM=1
  CALL ALLOCATE_FPRQPHY(GFP_PHYDS,YDRQAUX,IFPDOM,IPTR,IFPVT0(1:IPTR),NULOUT,LALLOFP)
ENDIF


!*       3.  PRINT OUT FINAL VALUES
!            ----------------------

IF (LTRACEFP) THEN
  WRITE(NULOUT,'(/,'' MODULE PTRFPB2 '')')
  WRITE(UNIT=NULOUT,FMT='('' NFPVT0 = '',I2)') YDRQAUX%NFIELDG
  IF (YDRQAUX%NFIELDG > 0) THEN
    IF (LARPEGEF) THEN
      WRITE(NULOUT,FMT='(''(JFLD GFP_PHYDS%CLNAME)'')')
      WRITE(UNIT=NULOUT,FMT='(6('' ('',I2,1X,A16,'')''))')&
       & (JFLD,GFP_PHYDS(YDRQAUX%ICOD(JFLD))%CLNAME, JFLD=1,YDRQAUX%NFIELDG)  
    ELSE
      WRITE(NULOUT,FMT='(''(JFLD GFP_PHYDS%IGRIB)'')')
      WRITE(UNIT=NULOUT,FMT='(6('' ('',I2,1X,I5,'')''))')&
       & (JFLD,GFP_PHYDS(YDRQAUX%ICOD(JFLD))%IGRIB, JFLD=1,YDRQAUX%NFIELDG)   
    ENDIF
  ENDIF
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPTR2',1,ZHOOK_HANDLE)

!------------------------------------------------------------------------------

CONTAINS

SUBROUTINE SET_NEXT(KCOD,CD)

INTEGER(KIND=JPIM), INTENT(IN) :: KCOD
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: CD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUFPTR2:SET_NEXT',0,ZHOOK_HANDLE)

IPTR=IPTR+1
IFPVT0(IPTR)=KCOD

IF (LHOOK) CALL DR_HOOK('SUFPTR2:SET_NEXT',1,ZHOOK_HANDLE)

END SUBROUTINE SET_NEXT

END SUBROUTINE SUFPTR2


