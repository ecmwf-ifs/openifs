! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CHEM_MASSDIA(YDGEOMETRY,YDGFL,YDSURF,YDML_GCONF,YDCOMPO,YDCHEM,YDDPHY,YDSPEC)

!***** CHEM_MASSDIA
!
!**   PURPOSE.
!     --------
!     Calculate total integral of CHEM 3D tracers, tendencies and fluxes for
!     monitoring mass conversation, use GFLT
!
!**   INTERFACE.
!     ----------
!              *CHEM_MASSDIA
!
!
!        Explicit arguments :
!        --------------------
!        none
!
!
!     AUTHOR.
!     -------
!     2009-12-02  J. Flemming
!
!     MODIFICATIONS.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMCOMPO               , ONLY : TCOMPO
USE YOMCHEM                , ONLY : TCHEM, IEXTR_PH, IEXTR_WD, IEXTR_EM, IEXTR_NG, IEXTR_ROH, IEXTR_CH, IEXTR_CHTR, IEXTR_FE,&
 &                                  IEXTR_ROH_TROP, IEXTR_DD, IEXTR_CHEMX,IEXTR_SDM, IEXTR_COND, IEXTR_NUC, IEXTR_RMOD
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX     , ONLY : TSURF
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK                , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN                 , ONLY : NULOUT
USE YOMMP0                 , ONLY : MYPROC
USE YOMCST                 , ONLY : RPI, RA
USE YOMCT3                 , ONLY : NSTEP
USE YOMDPHY                , ONLY : TDPHY
USE TM5_CHEM_MODULE        , ONLY : NRJ, NBUD_EXTRA_CHEM, NBUD_EXTRA, RATNAM, CHEM_BUDG, NSOA_BUDG
USE TM5_PHOTOLYSIS         , ONLY : NPHOTO




USE SPECTRAL_FIELDS_MOD    , ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL
TYPE(TSURF)         ,INTENT(INOUT) :: YDSURF
TYPE(TCOMPO)        ,INTENT(IN)    :: YDCOMPO
TYPE(TCHEM)         ,INTENT(INOUT) :: YDCHEM
TYPE(TDPHY)         ,INTENT(INOUT) :: YDDPHY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(SPECTRAL_FIELD),INTENT(IN)    :: YDSPEC
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM):: IERROR, ICHEM,  I, IGFL, JGFL, JSTGLO, ICHEMDV, ICHEMFLX, ICHEMSCAV, IPHOTO
INTEGER(KIND=JPIM):: ICEND, ISTC, IBL, JLEV, J, JR, ICASE, ICCASE, IGRIBCODE, IEXTR_MAX
INTEGER(KIND=JPIM):: INPHOTO, NNBUD_EXTRA, INSOA_BUDG

INTEGER(KIND=JPIM), PARAMETER :: ICCASEMAX=29

REAL (KIND=JPRB) ::  ZPRH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)  !
REAL (KIND=JPRB), ALLOCATABLE ::  ZDELP(:,:,:)  ! (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG+1,IBLK)
REAL (KIND=JPRB), ALLOCATABLE ::  ZPRESH(:,:,:)  ! (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,IBLK)
REAL (KIND=JPRB), ALLOCATABLE ::  ZTC(:,:,:)  ! (YDGEOMETRY%YRDIM%NPROMA,IBLK)
REAL(KIND=JPRB),  ALLOCATABLE ::  ZLOCGP(:,:,:)  ! (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, NGPBLKS2)
REAL(KIND=JPRB),  ALLOCATABLE ::  ZWORK(:)
REAL (KIND=JPRB), ALLOCATABLE ::  ZTROPO(:,:)  ! ((YDGEOMETRY%YRDIM%NPROMA,IBLK)

!REAL(KIND=JPRB) :: ZPTROP
REAL(KIND=JPRB) :: ZAREA, ZDELTA, ZRESIDUAL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZMASSADJ=1.0043_JPRB
                           ! adjustment factor when computing total mass
                           ! to conmpensate perfect sphere assumption.

CHARACTER(LEN=8), DIMENSION(ICCASEMAX), PARAMETER :: &
             & CLCASEN =(/'TOT_MASS','TRP_MASS','NEG_MASS',&
! Source terms
             &'EMIS_FLX','DDEP_FLX','WDEP_FLX','CHEM_TND','CHEM_TRO','NEGA_FIX','FLUX_ERR','EMIS_ATM', &
             & 'SEDM_FLX','COND_FLX','RMOD_FLX','NUCL_FLX',&
! Detailed chemistry related budgets (LCHEM_DIAC)
             & 'ROH_TEND','ROH_TROP','PHO_TEND','PHO_TROP',&
             & 'TRPMS_NH','TRPMS_SH',&
! Reaction budgets specifically for NH- and SH-extratropics...
             & 'ROH_TRNH','ROH_TRSH','PHO_TRNH','PHO_TRSH',&
! Additional chemistry budgets (O3 budget, heterogeneous reactions, ...)
             & 'CHX_TEND','CHX_TROP','CHX_TRNH','CHX_TRSH' /)
REAL(KIND=JPRB),  DIMENSION(3,ICCASEMAX)       :: ZTM
REAL(KIND=JPRB),  DIMENSION(220), SAVE  :: ZINIMASS, ZSCALMASS
REAL(KIND=JPRB),  DIMENSION(220,4), SAVE  :: ZCORMASS
CHARACTER(LEN=17) :: CLNAME

LOGICAL,  DIMENSION(ICCASEMAX)       :: LLBDG

!-----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gphpre.intfb.h"
#include "gptco3.intfb.h"
#include "gpnorm1.intfb.h"
#include "speree.intfb.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHEM_MASSDIA',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDVAB=>YDGEOMETRY%YRVAB, &
 & YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,YDCVER=>YDGEOMETRY%YRCVER, &
 & YDSTA=>YDGEOMETRY%YRSTA, YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
  & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDRIP=>YDML_GCONF%YRRIP,  &
  & YGFL=>YDML_GCONF%YGFL)

ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NCHEM=>YGFL%NCHEM, NGHG=>YGFL%NGHG,NUMFLDS=>YGFL%NUMFLDS, &
 & YAERO=>YGFL%YAERO, YCHEM=>YGFL%YCHEM, YGHG=>YGFL%YGHG, YCOMP=>YGFL%YCOMP, &
 & IEXTR_FREE=>YDCHEM%IEXTR_FREE, LCHEM_DIAC=>YDCHEM%LCHEM_DIAC, &
 & RCHEM_DIA_PERIOD=>YDCOMPO%RCHEM_DIA_PERIOD, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME, &
 & CHEM_SCHEME=>YDCHEM%CHEM_SCHEME, SMASSCOR=>YDCHEM%SMASSCOR, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCEXTR=>YDDPHY%NCEXTR, NVEXTR=>YDDPHY%NVEXTR, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GFL=>YDGFL%GFL, &
 & TSTEP=>YDRIP%TSTEP, &
 & SD_XA=>YDSURF%SD_XA)
!-----------------------------------------------------------------------

! test if output interval time
IF ( MOD((NSTEP*YDRIP%TSTEP)/3600.0_JPRB,RCHEM_DIA_PERIOD) /= 0.0_JPRB ) THEN
  IF (LHOOK) CALL DR_HOOK('CHEM_MASSDIA',1,ZHOOK_HANDLE)
  RETURN
ENDIF


! test if extra fields are defined
IF (NACTAERO > 0 .AND. TRIM(AERO_SCHEME)=='glomap') THEN
  IEXTR_MAX=IEXTR_RMOD
ELSE
  IEXTR_MAX=IEXTR_FE
ENDIF
IF ( NVEXTR < IEXTR_MAX ) THEN
  WRITE(NULOUT,*) ' NOT ENOUGH EXTRA FIELDS ',IEXTR_MAX, NVEXTR
  CALL ABOR1(" NOT ENOUGH EXTRA FIELDS FOR LCHEM_DIA ")
ENDIF
! The 6 here is in case they ugly-hack IEXTR_FREE fields that are used as
! standard rather than just for testing couldn't find any gaps to fill.
IF (NCHEM + NACTAERO + NGHG + 6  > NCEXTR) THEN
  WRITE(NULOUT,*) ' NOT ENOUGH LEVELS IN EXTRA field ',NCHEM+NACTAERO+NGHG+6, NCEXTR
  CALL ABOR1(" NOT ENOUGH LEVELS IN EXTRA field ")
ENDIF
IF (LCHEM_DIAC) THEN
  ! check dimensions of extra fields
  IF ( NVEXTR < IEXTR_CHEMX ) THEN
    WRITE(NULOUT,*) ' NOT ENOUGH EXTRA FIELDS ',IEXTR_CHEMX, NVEXTR
    CALL ABOR1(" NOT ENOUGH EXTRA FIELDS FOR LCHEM_DIAC ")
  ENDIF
ENDIF

! set ICCASE - number of diagnostics output for each species

IF (LCHEM_DIAC) THEN
  ! Extra diagnostics enabled
  ! additional diagnostics output below ('CHX_TEND','CHX_TROP','CHX_TRNH','CHX_TRSH' ) ! to be coded, not in species loop
  ICCASE=ICCASEMAX-4
ELSE
  ! standard output, budget contributions only
  IF (NACTAERO > 0 .AND. AERO_SCHEME=='GLOMAP') THEN
    ! GLOMAP has extra contributions
    ICCASE = 15
  ELSE
    ICCASE = 11
  ENDIF
ENDIF

LLBDG(:)=.FALSE.
LLBDG(4:7) =.TRUE.
LLBDG(9) =.TRUE.
LLBDG(11) =.TRUE.

! allocate work array
IF(.NOT.(ALLOCATED(ZLOCGP)))ALLOCATE(ZLOCGP(NPROMA,NFLEVG,NGPBLKS),STAT=IERROR)

! allocate half level pressures for conversions
IF(.NOT.(ALLOCATED(ZPRESH)))ALLOCATE(ZPRESH(NPROMA,0:NFLEVG,NGPBLKS),STAT=IERROR)
IF(.NOT.(ALLOCATED(ZDELP)))ALLOCATE(ZDELP(NPROMA,NFLEVG,NGPBLKS),STAT=IERROR)
IF(.NOT.(ALLOCATED(ZWORK)))ALLOCATE(ZWORK(NGPTOT),STAT=IERROR)

! total column
IF(.NOT.(ALLOCATED(ZTC))) ALLOCATE( ZTC(NPROMA,ICCASE,NGPBLKS), STAT=IERROR )

! tropopause height, calculated in chem_main
IF(.NOT.(ALLOCATED(ZTROPO))) ALLOCATE( ZTROPO(NPROMA,NGPBLKS), STAT=IERROR )

! convert surface pressure to gridpoint, store in zwork
CALL SPEREE(YDGEOMETRY,1,1,YDSPEC%SP,ZWORK)

ZPRESH(:,:,:) = 0.0_JPRB
ZDELP(:,:,:) = 0.0_JPRB
ZPRH(:,:) = 0.0_JPRB

!$OMP PARALLEL DO PRIVATE (JSTGLO, ICEND,ISTC, IBL, J, ZPRH, I)
! calculate p and delta p
DO JSTGLO=1,NGPTOT,NPROMA
  ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
  ISTC=1
  IBL=(JSTGLO-1)/NPROMA+1
  DO J=ISTC,ICEND
    ZPRH(J,NFLEVG) = EXP(ZWORK(JSTGLO-1+J))
! copy tropopause pressure
    ZTROPO(J,IBL) =  SD_XA(J,IEXTR_FREE(1,6),IEXTR_EM,IBL)
  ENDDO

  CALL GPHPRE(NPROMA,NFLEVG,ISTC,ICEND,YDVAB,YDCVER,ZPRH)
  DO I=1, NFLEVG
    ZDELP(ISTC:ICEND,I,IBL)=  ZPRH(ISTC:ICEND,I) - ZPRH(ISTC:ICEND,I-1)
    ZPRESH(ISTC:ICEND,I,IBL)=  ZPRH(ISTC:ICEND,I)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!loop over species
DO ICHEM=1,NCHEM + NACTAERO + NGHG

! find array positions
  IGFL=0
  ICHEMFLX=0
  ICHEMDV=0
  ICHEMSCAV=0
  IPHOTO=0
  CLNAME='XXXXXXXXXX'

  ! chemistry array
  !IF ( TRIM(YCHEM(ICHEM)%CNAME) == "PSC" ) CYCLE

  ! chemistry
  IF ( ICHEM <= NCHEM ) THEN
    CLNAME=YCHEM(ICHEM)%CNAME
    IGRIBCODE=YCHEM(ICHEM)%IGRBCODE
  ! Find emissions flux within surface fields
    !IF( YCHEM(ICHEM)%IGRIBSFC > 0 )
    ICHEMFLX = ICHEM

  ! Find deposition velocity within surface fields
    IF( YCHEM(ICHEM)%IGRIBDV > 0 .OR. TRIM(YCHEM(ICHEM)%CNAME) == 'CH4' ) THEN
      ICHEMDV = ICHEM
    ENDIF

  ! Find wet deposition species in levels in extra field 3
    IF( YCHEM(ICHEM)%HENRYA > 0) ICHEMSCAV = ICHEM

    SELECT CASE (TRIM(CHEM_SCHEME))

    CASE ("tm5","bascoetm5")

    DO JR = 1, NPHOTO
      IF (NRJ(JR) == ICHEM ) IPHOTO = JR
    ENDDO

      INPHOTO=NPHOTO
      INSOA_BUDG=NSOA_BUDG
    CASE DEFAULT

      ! Use dummy-type as default - there is no information content though
      INPHOTO=1
      INSOA_BUDG=2

    END SELECT

   ENDIF


! aerosols
   IF  ( ICHEM > NCHEM .AND. ICHEM <= NCHEM + NACTAERO  ) THEN
     CLNAME=YAERO(ICHEM-NCHEM)%CNAME
     IGRIBCODE=YAERO(ICHEM-NCHEM)%IGRBCODE
!    no trigger available - so always compute
     ICHEMSCAV = ICHEM
     ICHEMFLX = ICHEM
     ICHEMDV = ICHEM

!    SOA photolysis
     IF (CLNAME == "SecOrg_Bio") THEN
       IPHOTO=INPHOTO+1
     ELSEIF (CLNAME == "SecOrg_Anth" ) THEN
       IPHOTO=INPHOTO+2
     ENDIF
   ENDIF

! GHG
   IF  ( ICHEM > NCHEM + NACTAERO .AND. ICHEM <= NCHEM + NACTAERO + NGHG  ) THEN
     CLNAME=TRIM(YGHG(ICHEM-NCHEM-NACTAERO)%CNAME)
     ! Append _GHG if it's not already there
     IF (INDEX(CLNAME,'_GHG',BACK=.TRUE.) /= MAX(1,LEN_TRIM(CLNAME) - 3)) THEN
       CLNAME=TRIM(CLNAME) // '_GHG'
     ENDIF
     IGRIBCODE=YGHG(ICHEM-NCHEM-NACTAERO)%IGRBCODE
!    no trigger available - so always compute
     ICHEMFLX = ICHEM
   ENDIF

 !Find species within GFL range
   DO JGFL=1,NUMFLDS
      IF( YCOMP(JGFL)%IGRBCODE == IGRIBCODE ) IGFL=JGFL
   ENDDO

   IF (IGFL == 0_JPIM) THEN
      WRITE(NULOUT,*) ' massdia not found ', CLNAME, IGRIBCODE, NUMFLDS, ICHEM
      CALL ABOR1(" massdia chem not found " // CLNAME )
   ENDIF

! loop over differen contributions to mass budget
  DO ICASE= 1,ICCASE

    !$OMP PARALLEL DO PRIVATE(IBL)
    DO IBL=1,NGPBLKS
      ZTC(:,ICASE,IBL)=0.0_JPRB
    ENDDO
    !$OMP END PARALLEL DO

    SELECT CASE (TRIM( CLCASEN(ICASE)))

    CASE ('TOT_MASS' )

!   copy GFL field to local array
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,JLEV,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
          DO JLEV=1,NFLEVG
           DO J=ISTC,ICEND
              ZLOCGP(J,JLEV,IBL) = GFL(J,JLEV,YCOMP(IGFL)%MP, IBL)
           ENDDO
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
 ! calculate total col
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        CALL GPTCO3(NPROMA,ISTC,ICEND, NFLEVG,ZDELP(:,:,IBL),&
        & ZLOCGP(:,:,IBL),ZTC(:,ICASE,IBL))

       ENDDO
       !$OMP END PARALLEL DO

  CASE ('TRP_MASS' )
!   copy tropospheric GFL field to local array

       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,JLEV,J)
       DO JSTGLO=1,NGPTOT,NPROMA
         ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
         ISTC=1
         IBL=(JSTGLO-1)/NPROMA+1

         DO JLEV=1,NFLEVG
           DO J=ISTC,ICEND

              IF (ZPRESH(J,JLEV,IBL) < ZTROPO(J,IBL)) THEN
        ! in stratosphere set concentrations to zero
                ZLOCGP(J,JLEV,IBL) = 0._JPRB
              ELSE
                ! in troposphere fill with actual conc.
                ZLOCGP(J,JLEV,IBL) = GFL(J,JLEV,YCOMP(IGFL)%MP, IBL)
              ENDIF
           ENDDO
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
 ! calculate partial tropospheric col
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        CALL GPTCO3(NPROMA,ISTC,ICEND, NFLEVG,ZDELP(:,:,IBL),&
        & ZLOCGP(:,:,IBL),ZTC(:,ICASE,IBL))

       ENDDO
       !$OMP END PARALLEL DO

  CASE ('TRPMS_NH' )
!   copy tropospheric GFL field over extra-tropical  NH to local array
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,JLEV,J)
       DO JSTGLO=1,NGPTOT,NPROMA
         ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
         ISTC=1
         IBL=(JSTGLO-1)/NPROMA+1

         DO JLEV=1,NFLEVG
           DO J=ISTC,ICEND

              IF (ZPRESH(J,JLEV,IBL) < ZTROPO(J,IBL) .OR. YDGSGEOM_NB%GELAT(J+JSTGLO-1) < (0.1667_JPRB * RPI) ) THEN
        ! in stratosphere and outside extra-tropical NH set concentrations to zero
                ZLOCGP(J,JLEV,IBL) = 0._JPRB
              ELSE
                ! in extra-trop NH troposphere fill with actual conc.
                ZLOCGP(J,JLEV,IBL) = GFL(J,JLEV,YCOMP(IGFL)%MP, IBL)
              ENDIF
           ENDDO
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO

 ! calculate partial tropospheric col over NH
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        CALL GPTCO3(NPROMA,ISTC,ICEND, NFLEVG,ZDELP(:,:,IBL),&
        & ZLOCGP(:,:,IBL),ZTC(:,ICASE,IBL))

       ENDDO
       !$OMP END PARALLEL DO

  CASE ('TRPMS_SH' )
!   copy tropospheric GFL field over extra-tropical SH to local array
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,JLEV,J)
       DO JSTGLO=1,NGPTOT,NPROMA
         ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
         ISTC=1
         IBL=(JSTGLO-1)/NPROMA+1

         DO JLEV=1,NFLEVG
           DO J=ISTC,ICEND
             IF (ZPRESH(J,JLEV,IBL) < ZTROPO(J,IBL) .OR. YDGSGEOM_NB%GELAT(J+JSTGLO-1) > (-0.1667_JPRB * RPI) ) THEN
                ! in stratosphere and outside extra-tropical SH set concentrations to zero
                ZLOCGP(J,JLEV,IBL) = 0._JPRB
              ELSE
                ! in extra-trop SH troposphere fill with actual conc.
                ZLOCGP(J,JLEV,IBL) = GFL(J,JLEV,YCOMP(IGFL)%MP, IBL)
              ENDIF
           ENDDO
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
 ! calculate partial tropospheric col over SH
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        CALL GPTCO3(NPROMA,ISTC,ICEND, NFLEVG,ZDELP(:,:,IBL),&
        & ZLOCGP(:,:,IBL),ZTC(:,ICASE,IBL))

       ENDDO
       !$OMP END PARALLEL DO

   CASE ('NEG_MASS' )
!   copy GFL field to local array
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,JLEV,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
          DO JLEV=1,NFLEVG
           DO J=ISTC,ICEND
              ZLOCGP(J,JLEV,IBL) =0.0_JPRB
              ZLOCGP(J,JLEV,IBL) = MIN(ZLOCGP(J,JLEV,IBL) , GFL(J,JLEV,YCOMP(IGFL)%MP, IBL))
           ENDDO
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
 ! caluclate total col
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
         CALL GPTCO3(NPROMA,ISTC,ICEND, NFLEVG,ZDELP(:,:,IBL),&
        & ZLOCGP(:,:,IBL),ZTC(:,ICASE,IBL))
       ENDDO
       !$OMP END PARALLEL DO
   CASE ('EMIS_FLX' )

      IF (ICHEMFLX > 0 ) THEN
! map emissions flux into local array
        !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL)
        DO JSTGLO=1,NGPTOT,NPROMA
          ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
          ISTC=1
          IBL=(JSTGLO-1)/NPROMA+1
          DO J=ISTC,ICEND
              ZTC(J,ICASE,IBL) = (SD_XA(J,ICHEMFLX,IEXTR_EM,IBL))
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF

  CASE ('EMIS_ATM' )
!     add air craft emissions
     IF( TRIM(CLNAME) == "NO")  THEN
      !   3D NO emissions are accumulated in PEXTRA(:,IEXTR_FREE(1,1:2),IEXTR_EM)
      !   map emissions flux into local array
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) =  SD_XA(J,IEXTR_FREE(1,2),IEXTR_EM,IBL)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF
    IF (NCHEM > 0) THEN ! CHEM_SCHEME may not be defined otherwise
!     add AC emissions as NO2
      IF( TRIM(CLNAME) == "NO2")  THEN
        !   3D NO emissions are accumulated in PEXTRA(:,IEXTR_FREE(1,1:2),IEXTR_EM)
        !   map emissions flux into local array
        !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
        DO JSTGLO=1,NGPTOT,NPROMA
          ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
          ISTC=1
          IBL=(JSTGLO-1)/NPROMA+1
          DO J=ISTC,ICEND
            ZTC(J,ICASE,IBL) = SD_XA(J,IEXTR_FREE(1,1),IEXTR_EM,IBL)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
      IF( TRIM(CHEM_SCHEME) == "tm5" .AND.  TRIM(CLNAME) == "CH4")  THEN
      !      see chem_main , UBC CH4 are stored in IEXTR_FREE(1,3)
        !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
        DO JSTGLO=1,NGPTOT,NPROMA
          ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
          ISTC=1
          IBL=(JSTGLO-1)/NPROMA+1
          DO J=ISTC,ICEND
            ZTC(J,ICASE,IBL) = SD_XA(J,IEXTR_FREE(1,3),IEXTR_EM,IBL)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
      IF(TRIM(CHEM_SCHEME) == "tm5" .AND. TRIM(CLNAME) == "HNO3")  THEN
!       see chem_main , UBC HNO3 are strored in IEXTR_FREE(1,4)
        !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
        DO JSTGLO=1,NGPTOT,NPROMA
          ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
          ISTC=1
          IBL=(JSTGLO-1)/NPROMA+1
          DO J=ISTC,ICEND
            ZTC(J,ICASE,IBL) = SD_XA(J,IEXTR_FREE(1,4),IEXTR_EM,IBL)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
    ENDIF
    IF( TRIM(CLNAME) == "CO2_GHG")  THEN
!     see chem_main , CO2 aircraft emissions are stored  IEXTR_FREE(1,5)
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
       ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
       ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) = SD_XA(J,IEXTR_FREE(1,5),IEXTR_EM,IBL)
        ENDDO
      ENDDO
     !$OMP END PARALLEL DO
    ENDIF



  CASE ('DDEP_FLX' )

    IF (ICHEMDV > 0  ) THEN
!   map emissions flux into local array
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
              ZTC(J,ICASE,IBL) = (SD_XA(J,ICHEMDV,IEXTR_DD,IBL))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF


  CASE ('WDEP_FLX' )

    IF (ICHEMSCAV > 0 ) THEN
! map wet depsition flux from extr2d field into local array
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL)= SD_XA(J,ICHEMSCAV,IEXTR_WD,IBL)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF

  CASE ('SEDM_FLX' )

      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL)= SD_XA(J,ICHEMSCAV,IEXTR_SDM,IBL)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

  CASE ('COND_FLX' )

      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL)= SD_XA(J,ICHEMSCAV,IEXTR_COND,IBL)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

  CASE ('RMOD_FLX' )

      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL)= SD_XA(J,ICHEMSCAV,IEXTR_RMOD,IBL)
        ENDDO
      ENDDO

  CASE ('NUCL_FLX' )

      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL)= SD_XA(J,ICHEMSCAV,IEXTR_NUC,IBL)
        ENDDO
      ENDDO

   CASE ('CHEM_TND' )
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
         ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
         ISTC=1
         IBL=(JSTGLO-1)/NPROMA+1
         DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_CH, IBL)
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
   CASE ('CHEM_TRO' )
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
         ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
         ISTC=1
         IBL=(JSTGLO-1)/NPROMA+1
         DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_CHTR, IBL)
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
  CASE ('ROH_TEND' )
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_ROH, IBL)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
 CASE ('ROH_TROP' )
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_ROH_TROP, IBL)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

  CASE ('ROH_TRNH' )
! reaction budget with OH, tropospheric NH-extra-tropics only
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          IF (YDGSGEOM_NB%GELAT(J+JSTGLO-1) > (0.1667_JPRB * RPI)) THEN
            ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_ROH_TROP, IBL)
          ELSE
            ZTC(J,ICASE,IBL) = 0._JPRB
          ENDIF
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
  CASE ('ROH_TRSH' )
! reaction budget with OH, tropospheric SH-extra-tropics only
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          IF (YDGSGEOM_NB%GELAT(J+JSTGLO-1) < (-0.1667_JPRB * RPI)) THEN
            ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_ROH_TROP, IBL)
          ELSE
            ZTC(J,ICASE,IBL) = 0._JPRB
          ENDIF
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

  CASE ('PHO_TEND' )
    IF (IPHOTO > 0 ) THEN
      !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
      DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) = SD_XA(J,IPHOTO,IEXTR_PH, IBL)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF

  CASE ('PHO_TROP' )
     IF (IPHOTO > 0 ) THEN
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) = SD_XA(J,IPHOTO+INPHOTO+INSOA_BUDG,IEXTR_PH, IBL)
        ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

  CASE ('PHO_TRNH' )
 ! reaction budget of photolysis, tropospheric NH-extra-tropics only
     IF (IPHOTO > 0 ) THEN
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          IF (YDGSGEOM_NB%GELAT(J+JSTGLO-1) > (0.1667_JPRB * RPI)) THEN
            ZTC(J,ICASE,IBL) = SD_XA(J,IPHOTO+INPHOTO+INSOA_BUDG,IEXTR_PH, IBL)
          ELSE
            ZTC(J,ICASE,IBL) = 0._JPRB
          ENDIF
        ENDDO
       ENDDO
       !$OMP END PARALLEL DO
     ENDIF

 CASE ('PHO_TRSH' )
! reaction budget of photolysis, tropospheric SH-extra-tropics only
     IF (IPHOTO > 0 ) THEN
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          IF (YDGSGEOM_NB%GELAT(J+JSTGLO-1) < (-0.1667_JPRB * RPI)) THEN
            ZTC(J,ICASE,IBL) = SD_XA(J,IPHOTO+INPHOTO+INSOA_BUDG,IEXTR_PH, IBL)
          ELSE
            ZTC(J,ICASE,IBL) = 0._JPRB
          ENDIF
        ENDDO
       ENDDO
       !$OMP END PARALLEL DO
     ENDIF

  CASE ('NEGA_FIX' )

! negative fix addition is stored in extra field iextr_ng , level indicates species

! copy GFL field to local array
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
         ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
         ISTC=1
         IBL=(JSTGLO-1)/NPROMA+1
         DO J=ISTC,ICEND
            ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_NG, IBL)
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO

  CASE ('FLUX_ERR' )

! negative fix addition is stored in extra field iextr_ng , level indicates species

! copy GFL field to local array
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
         ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
         ISTC=1
         IBL=(JSTGLO-1)/NPROMA+1
         DO J=ISTC,ICEND
           ZTC(J,ICASE,IBL) = SD_XA(J,ICHEM,IEXTR_FE, IBL)
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO

  CASE DEFAULT

    CALL ABOR1(" budget option not found: " // TRIM( CLCASEN(ICASE)))

  END SELECT

  ENDDO ! ICASE

  CALL GPNORM1(YDGEOMETRY,ZTC,ICCASE,.FALSE.,PNORMS=ZTM)


  ! output
  IF (MYPROC == 1) THEN

    ! multiply total ZTC by grid box area (estimate)
    ! does gaw also include reduced gaussian gris size ????
    ! earth surface in m**2
    ! convert in Tg

    ZAREA=4.0_JPRB*RPI*RA*RA*ZMASSADJ *  0.000000001_JPRB

    ZTM(1,1:ICCASE) = ZTM(1,1:ICCASE) * ZAREA


    ! set initial values
    IF (NSTEP == 0 ) THEN
      IF (NCHEM+NGHG+NACTAERO > 220 )  CALL ABOR1(" increase dimension of ZACC in CHEM_MASSDIA " )
      !  ZINIMASS(ICHEM) = PMASS(ICHEM) * ZAREA
      ZINIMASS(ICHEM) =  ZTM(1,1)
      ZCORMASS(ICHEM,1:4) = 0.0_JPRB
      ZSCALMASS(ICHEM) =  ZINIMASS(ICHEM)
      ! time step greater zero
    ELSE

      ! write header
      IF(ICHEM == 1 .AND. NSTEP*TSTEP/3600.0_JPRB == RCHEM_DIA_PERIOD) &
     &  WRITE(NULOUT,'(3a13,33a16)') ' MASSDIA ', ' NAME ', '  SIM_HOUR ',CLCASEN(1:ICCASE), &
     &  'MASS_CHG', 'RESIDUAL','MASS_CHG_r', 'RESIDUAL_r', ' SLMGLO_ER1 ' ,' SLMGLO_ER2 ', ' SLMTRO_ERR1 ' ,' SLMTRO_ERR2 '

      ! calculate total mass change in time step
      ZDELTA=ZTM(1,1)-ZINIMASS(ICHEM)
      ZRESIDUAL=ZDELTA

      ! accumulated mass SL deficit  - should be close to zero with MF
      ! Fails in single-precision mode.. Switch off when running in SP
      IF (JPRB == JPRD) THEN
        ZCORMASS(ICHEM,1)=ZCORMASS(ICHEM,1) + SMASSCOR(IGFL,1) * ZAREA
        ZCORMASS(ICHEM,2)=ZCORMASS(ICHEM,2) + SMASSCOR(IGFL,2) * ZAREA
        ZCORMASS(ICHEM,3)=ZCORMASS(ICHEM,3) + SMASSCOR(IGFL,3) * ZAREA
        ZCORMASS(ICHEM,4)=ZCORMASS(ICHEM,4) + SMASSCOR(IGFL,4) * ZAREA
      ENDIF

      ! calculate residual  (mass increase minus (accumulated) fluxes / tendencies / corrections)
      DO ICASE=1,ICCASE
        IF ( LLBDG(ICASE) )  ZRESIDUAL=ZRESIDUAL-ZTM(1,ICASE)
      ENDDO

      ! IF ( .NOT. ( YCHEM(ICHEM)%LMASSFIX .AND. YCHEM(ICHEM)%LADV ) )  ZRESIDUAL = ZRESIDUAL +  ZCORMASS(ICHEM)

      ! to avoid zero mass or strong dependance on diurnal cycle take average mass between start and nstep
      ZSCALMASS(ICHEM)  =  ( ZSCALMASS(ICHEM)*FLOAT(NSTEP) + ZTM(1,1) ) / FLOAT(NSTEP + 1 )
      IF ( ZSCALMASS(ICHEM)  == 0.0_JPRB)  ZSCALMASS(ICHEM) = -99999.9
      ! write to nulout
      WRITE(NULOUT,'(a10,a17,f10.2,330e16.8)') 'MASSDIA1 ', CLNAME, (NSTEP*TSTEP)/3600.0_JPRB,&
     &  ZTM(1,1:ICCASE), ZDELTA, ZRESIDUAL,&
     & 100.0_JPRB*ZDELTA / ZSCALMASS(ICHEM),100.0_JPRB*ZRESIDUAL / ZSCALMASS(ICHEM) ,ZCORMASS(ICHEM,1:4)
    ENDIF
! do we need a barrier here
  ENDIF

ENDDO ! ICHEM


! separate output for ozone diagnostics - loop over NBUD_EXTRA, all values in one row
IF (LCHEM_DIAC .AND. NSTEP > 0) THEN
! all budget
  NNBUD_EXTRA=NBUD_EXTRA





  DO ICASE = 1, NNBUD_EXTRA
!  ! Extra chemistry reaction budgets
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) = SD_XA(J,ICASE,IEXTR_CHEMX, IBL)
        ENDDO
       ENDDO
       !$OMP END PARALLEL DO
   ENDDO

   CALL GPNORM1(YDGEOMETRY,ZTC,ICCASE,.FALSE.,PNORMS=ZTM)
   IF (MYPROC == 1) THEN
      ! write header
      IF (CHEM_SCHEME == "tm5" .or. CHEM_SCHEME =="bascoetm5") THEN
        IF( NSTEP*TSTEP/3600.0_JPRB == RCHEM_DIA_PERIOD ) WRITE(NULOUT,'(3a10,30a16)')&
         & ' OZONBUD ', ' REGION ', '  SIM_HOUR ',RATNAM(CHEM_BUDG(1:NBUD_EXTRA_CHEM)),"JO2","JA_CH2O"





      ELSE
        IF( NSTEP*TSTEP/3600.0_JPRB == RCHEM_DIA_PERIOD ) WRITE(NULOUT,'(3a10,10I16)')&
         & ' OZONBUD ', ' REGION ', '  SIM_HOUR ' !,ICHEM_BUDG(1:NNBUD_EXTRA)
      ENDIF

      ZTM(1,1:NNBUD_EXTRA) = ZTM(1,1:NNBUD_EXTRA) * ZAREA
      WRITE(NULOUT,'(2a10,f10.2,40e16.8)') 'OZONBUD ', 'ALL', (NSTEP*TSTEP)/3600.0_JPRB,&
        &  ZTM(1,1:NNBUD_EXTRA)
   ENDIF

! trop budget
  DO ICASE = 1, NNBUD_EXTRA
!  ! Extra chemistry reaction budgets
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          ZTC(J,ICASE,IBL) = SD_XA(J,ICASE+NNBUD_EXTRA,IEXTR_CHEMX, IBL)
        ENDDO
       ENDDO
       !$OMP END PARALLEL DO
   ENDDO


   CALL GPNORM1(YDGEOMETRY,ZTC,ICCASE,.FALSE.,PNORMS=ZTM)
   IF (MYPROC == 1) THEN

      ZTM(1,1:NNBUD_EXTRA) = ZTM(1,1:NNBUD_EXTRA) * ZAREA

      WRITE(NULOUT,'(2a10,f10.2,40e16.8)') 'OZONBUD ', 'TROP', (NSTEP*TSTEP)/3600.0_JPRB,&
     &  ZTM(1,1:NNBUD_EXTRA)
   ENDIF

! extra chemistry budget, tropospheric NH-extra-tropics only
  DO ICASE = 1, NNBUD_EXTRA
!  ! Extra chemistry reaction budgets
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
          IF (YDGSGEOM_NB%GELAT(J+JSTGLO-1) > (0.1667_JPRB * RPI)) THEN
            ZTC(J,ICASE,IBL) = SD_XA(J,ICASE+NNBUD_EXTRA,IEXTR_CHEMX, IBL)
          ELSE
            ZTC(J,ICASE,IBL) = 0._JPRB
          ENDIF
        ENDDO
       ENDDO
       !$OMP END PARALLEL DO
   ENDDO


   CALL GPNORM1(YDGEOMETRY,ZTC,ICCASE,.FALSE.,PNORMS=ZTM)
   IF (MYPROC == 1) THEN

      ZTM(1,1:NNBUD_EXTRA) = ZTM(1,1:NNBUD_EXTRA) * ZAREA

      WRITE(NULOUT,'(2a10,f10.2,40e16.8)') 'OZONBUD ', 'TRNH', (NSTEP*TSTEP)/3600.0_JPRB,&
     &  ZTM(1,1:NNBUD_EXTRA)
   ENDIF


! trop budget SH extratropics
  DO ICASE = 1, NNBUD_EXTRA
!  ! Extra chemistry reaction budgets
       !$OMP PARALLEL DO PRIVATE(JSTGLO,ICEND,ISTC,IBL,J)
       DO JSTGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        ISTC=1
        IBL=(JSTGLO-1)/NPROMA+1
        DO J=ISTC,ICEND
           IF (YDGSGEOM_NB%GELAT(J+JSTGLO-1) < (-0.1667_JPRB * RPI)) THEN
             ZTC(J,ICASE,IBL) = SD_XA(J,ICASE+NNBUD_EXTRA,IEXTR_CHEMX, IBL)
           ELSE
             ZTC(J,ICASE,IBL) = 0._JPRB
           ENDIF
        ENDDO
       ENDDO
       !$OMP END PARALLEL DO
   ENDDO

   CALL GPNORM1(YDGEOMETRY,ZTC,ICCASE,.FALSE.,PNORMS=ZTM)
   IF (MYPROC == 1) THEN

      ZTM(1,1:NNBUD_EXTRA) = ZTM(1,1:NNBUD_EXTRA) * ZAREA

      WRITE(NULOUT,'(2a10,f10.2,40e16.8)') 'OZONBUD ', 'TRSH', (NSTEP*TSTEP)/3600.0_JPRB,&
     &  ZTM(1,1:NNBUD_EXTRA)
   ENDIF
ENDIF

IF (ALLOCATED(ZDELP)  ) DEALLOCATE (ZDELP)
IF (ALLOCATED(ZPRESH) ) DEALLOCATE (ZPRESH)
IF (ALLOCATED(ZLOCGP) ) DEALLOCATE (ZLOCGP)
IF (ALLOCATED(ZTC) )    DEALLOCATE (ZTC)
IF (ALLOCATED(ZTROPO) )    DEALLOCATE (ZTROPO)

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_MASSDIA',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_MASSDIA
