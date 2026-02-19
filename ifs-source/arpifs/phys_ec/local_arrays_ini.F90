! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LOCAL_ARRAYS_INI(YDGEOMETRY,YDSURF,YDMODEL,KDIM, LLKEYS, PAUX, AUXL, SURFL, PERTL, GEMSL,&
  & PRAD, PSURF, PGFL, PTENGFL)

!**** *Initialization of derived variables local to CALLPAR* - creator

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------

! Derived arrays           Reserved space
! -----------------------------------------------
! KDIM     : Derived variable for dimensions

! Global derived variables for initialization
! PRAD   - quantities for radiation
! PAUX   - auxiliary quantities
! PSURF  - surface quantities
! PGFL   - GFL
! PTENGFL - tendency of GFL (GFL time T1)


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-12-05  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     JJMorcrette 20131001 DMS-related varaibles
!     S.Remy      20150128 Injection heights for reactive and greenhouse gases 
!     E. Dutra/G.Arduini Jan 2018: Snow mass and density in SURFL are 2D, snow multi-layer 
!     R. Hogan    20190123 Longwave surface emissivity in NLWEMISS intervals
!-----------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK,  DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE, KEYS_LOCAL_TYPE, AUX_DIAG_LOCAL_TYPE,&
  & AUX_RAD_TYPE, SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE, PERTURB_LOCAL_TYPE,&
  & GEMS_LOCAL_TYPE, AUX_TYPE
USE YOMCTESSELDIM, ONLY: IKVTYPES, IKDHVVEGS, IKDHFVEGS

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)                  ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)                     ,INTENT(INOUT) :: YDSURF
TYPE(MODEL)                     ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)           ,INTENT (IN)   :: KDIM
TYPE (KEYS_LOCAL_TYPE)          ,INTENT(INOUT) :: LLKEYS
TYPE (AUX_TYPE)                 ,INTENT(IN)    :: PAUX
TYPE (AUX_DIAG_LOCAL_TYPE)      ,INTENT(INOUT) :: AUXL
TYPE (SURF_AND_MORE_LOCAL_TYPE) ,INTENT(INOUT) :: SURFL
TYPE (PERTURB_LOCAL_TYPE)       ,INTENT(INOUT) :: PERTL
TYPE (GEMS_LOCAL_TYPE)          ,INTENT(INOUT) :: GEMSL
TYPE (AUX_RAD_TYPE)             ,INTENT(IN)    :: PRAD
TYPE (SURF_AND_MORE_TYPE)       ,INTENT(INOUT) :: PSURF
REAL(KIND=JPRB)                 ,INTENT(IN)    :: PGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)                 ,INTENT(IN)    :: PTENGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM1)
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "gems_init.intfb.h"
#include "aero_init.intfb.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LOCAL_ARRAYS_INI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 &  YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO,YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH, &
  & YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
  & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
  & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, &
 & NAERO=>YGFL%NAERO, NCHEM=>YGFL%NCHEM, NCHEM_DV=>YGFL%NCHEM_DV,NGEMS=>YGFL%NGEMS, &
 & LEFLAKE=>YDEPHY%LEFLAKE, LEOCML=>YDEPHY%LEOCML, &
 & NTSW=>YDERAD%NTSW, NLWEMISS=>YDERAD%NLWEMISS, &
 & LSTOPH_CASBS=>YDSTOPH%LSTOPH_CASBS, LSTOPH_SPBS=>YDSTOPH%LSTOPH_SPBS, &
 & YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSD_VN=>YDSURF%YSD_VN, & 
 & NGHG=>YGFL%NGHG,TSPHY=>YDPHY2%TSPHY,LCHEM_DIA=>YDCOMPO%LCHEM_DIA,YGHG=>YGFL%YGHG)
!-----------------------------------------------------------------------

! -----------------------------------------
!*  LLKEYS
! -----------------------------------------

! initialization...
LLKEYS%LLSHCV(KDIM%KIDIA:KDIM%KFDIA)=.TRUE.
LLKEYS%LLCUM(KDIM%KIDIA:KDIM%KFDIA)=.FALSE.
LLKEYS%LLSC(KDIM%KIDIA:KDIM%KFDIA) =.FALSE.
IF (.NOT. LEFLAKE) LLKEYS%LLLAKE(KDIM%KIDIA:KDIM%KFDIA)    = .FALSE.
IF (.NOT. LEOCML)  LLKEYS%LLOCN_KPP(KDIM%KIDIA:KDIM%KFDIA) = .FALSE.

! -----------------------------------------
!*  AUXL  
! -----------------------------------------

! Initialization
DO JL=KDIM%KIDIA,KDIM%KFDIA
  AUXL%ZSUDU (JL)  =PRAD%PSRSWDCS(JL)
  AUXL%ZCEMTR(JL,0)=PRAD%PEMTC(JL,1)
  AUXL%ZCEMTR(JL,1)=PRAD%PEMTC(JL,KDIM%KLEV+1)
  AUXL%ZCTRSO(JL,0)=PRAD%PTRSC(JL,1)
  AUXL%ZCTRSO(JL,1)=PRAD%PTRSC(JL,KDIM%KLEV+1)
  AUXL%ZTRSOD(JL)  =PRAD%PSRSWD(JL) 
  AUXL%ITOPC(JL)=INT(PSURF%PSD_VN(JL,YSD_VN%YTOP%MP)+0.01_JPRB)
  AUXL%IBASC(JL)=INT(PSURF%PSD_VN(JL,YSD_VN%YBAS%MP)+0.01_JPRB)
  AUXL%ZRAINFRAC_TOPRFZ(1:KDIM%KLON)=0.0_JPRB
ENDDO


! -----------------------------------------
!*    SURFL
! -----------------------------------------

! To reduce code changes, point the "traditional" two values of
! emissivity to the first two columns of the spectral emissivity
! matrix
SURFL%ZEMIR => SURFL%ZSPECTRALEMISS(:,1)
SURFL%ZEMIW => SURFL%ZSPECTRALEMISS(:,2)

! initialization
SURFL%ZLAILI    = 0._JPRB
SURFL%ZLAIHI    = 0._JPRB
! Sea-ice
SURFL%ZTHKICE(:) = 0.0_JPRB
SURFL%ZSNTICE(:) = 0.0_JPRB
!CTESSEL dynamics initialized to zero
SURFL%ZANDAYVT(:,:) = 0._JPRB
SURFL%ZANFMVT(:,:)  = 0._JPRB
SURFL%ZAHFTRTI(:,:) = 0.0_JPRB

! -----------------------------------------
!*    PERTL
! -----------------------------------------

!  initialization
IF (LSTOPH_SPBS.OR.LSTOPH_CASBS) THEN
  PERTL%ZDISSGW(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)   = 0.0_JPRB
ENDIF



! -----------------------------------------
!*     GEMSL
! -----------------------------------------

! Integer scalars

! Association of pointers
ALLOCATE(GEMSL%ITRAC)
ALLOCATE(GEMSL%IAERO(MAX(1,NAERO)))
ALLOCATE(GEMSL%ICHEM(MAX(1,NCHEM)))

! Real arrays

! Association of pointers

! These need to be full arrays for the output from the turbulence scheme, even
! when GEMS is turned off, since there is no option to disable turbulence_layer
! from writing them out.
ALLOCATE(GEMSL%ZAZ0M(1:KDIM%KLON))
ALLOCATE(GEMSL%ZAZ0H(1:KDIM%KLON))

! Likewise, for their use in the visibility scheme / with aerosol climatology
! Initialize these to zero since they are not calculated on every time step, but
! find their way into AEROUT fields which are processed by GPNORMs which produces
! SIGFPE if they are NaN.
ALLOCATE(GEMSL%ZCLAERS(1:KDIM%KLON))
ALLOCATE(GEMSL%ZPRAERS(1:KDIM%KLON))
ALLOCATE(GEMSL%ZVISICL(1:KDIM%KLON))
GEMSL%ZCLAERS(:) = 0.0_JPRB
GEMSL%ZPRAERS(:) = 0.0_JPRB
GEMSL%ZVISICL(:) = 0.0_JPRB

! these field could perhaps be allocated (1,1,..), if not needed (NACTAERO=0,  YLRCH4%LGP) 
ALLOCATE(GEMSL%ZLRCH4(1:KDIM%KLON,1:KDIM%KLEV))
ALLOCATE(GEMSL%ZKOZO(1:KDIM%KLON,1:KDIM%KLEV,1:KDIM%KVCLIS))

! Initialization + Allocation
! This code is being inherited from the past.  Some unifications with the other GEMS pointed would be desirable.
GEMSL%ITRAC=0
! Allocate general variables for tracer transport: 
! IMPORTANT: Tracer order is : CO2 - other tracers - react Gases - Aerosol - extra GFL
IF (NGHG > 0 .OR. NCHEM > 0 .OR. NAERO > 0 ) THEN
  !CALL GEMS_INIT( KDIM, GEMSL, PSURF, PGFL, PTENGFL)
  CALL GEMS_INIT(YDSURF, YDMODEL%YRML_CHEM,YDEPHY,YGFL,YDPHY2,YDMODEL%YRML_GCONF%YRRIP,&
      & KDIM%KIDIA, KDIM%KFDIA, KDIM%KLEV, KDIM%KLON, GEMSL%ITRAC,&
      & PAUX%PGELAM, PAUX%PGELAT, PSURF%PSD_VF, GEMSL%IAERO, GEMSL%ICHEM, GEMSL%ZLRCH4, PGFL,&
      & PTENGFL, GEMSL%ZCEN, GEMSL%ZTENC, GEMSL%ZTENC_SKF, GEMSL%ZCFLX, GEMSL%ZSCAV) 

  ALLOCATE(GEMSL%ZCFLXO(1:KDIM%KLON,GEMSL%ITRAC))
  GEMSL%ZCFLXO(:,:)= 0.0_JPRB 
  ALLOCATE(GEMSL%ZDDVLC(1:KDIM%KLON,GEMSL%ITRAC))
  GEMSL%ZDDVLC(:,:)= 0.0_JPRB
 ! allocate dry deposition field 
  IF (NCHEM_DV > 0 ) THEN
    ALLOCATE(  GEMSL%ZCHEMDV(1:KDIM%KLON,NCHEM_DV) )
    GEMSL%ZCHEMDV(:,:)=0.0_JPRB
  ELSE
    ALLOCATE( GEMSL%ZCHEMDV(1,1) )
    GEMSL%ZCHEMDV(1,1)=0.0_JPRB
  ENDIF

ELSE
! ALLOCATE arrays so they can be passed to turbulence scheme etc. in the non-GEMS case
  ALLOCATE(GEMSL%ZCEN(1,1,1))
  ALLOCATE(GEMSL%ZTENC(1,1,1))
  ALLOCATE(GEMSL%ZTENC_SKF(1,1,1))
  ALLOCATE(GEMSL%ZCFLX(1,1))
  ALLOCATE(GEMSL%ZCFLXO(1,1))
  ALLOCATE(GEMSL%ZDDVLC(1,1))
  ALLOCATE(GEMSL%ZSCAV(1))
  ALLOCATE( GEMSL%ZCHEMDV(1,1) )
  GEMSL%ZCEN(1,1,1) = 0.0_JPRB
  GEMSL%ZTENC(1,1,1) = 0.0_JPRB
  GEMSL%ZTENC_SKF(1,1,1)= 0.0_JPRB
  GEMSL%ZCFLX(1,1) = 0.0_JPRB
  GEMSL%ZCFLXO(1,1) = 0.0_JPRB
  GEMSL%ZDDVLC(1,1) = 0.0_JPRB
  GEMSL%ZSCAV(1) = 0.0_JPRB
  GEMSL%ZCHEMDV(1,1)=0.0_JPRB
ENDIF

! Allocate local arrays for aerosol transport
IF (NACTAERO > 0) THEN
  ALLOCATE(GEMSL%ZAERWS(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZAERGUST(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZAERUST(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZDIST(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZAERMAP(1:KDIM%KLON,1:5))

  ! additional ones due to DMSO
  ALLOCATE(GEMSL%ZDMSO(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZLDAY(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZLISS(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZSO2(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZTDMS(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZDMSI(1:KDIM%KLON))
  ALLOCATE(GEMSL%ZODMS(1:KDIM%KLON))

  CALL AERO_INIT(YGFL,YDCOMPO,YDEAERATM,KDIM,GEMSL)
ELSE
  ALLOCATE(GEMSL%ZAERWS(1))
  ALLOCATE(GEMSL%ZAERGUST(1))
  ALLOCATE(GEMSL%ZAERUST(1))
  ALLOCATE(GEMSL%ZDIST(1))
  ALLOCATE(GEMSL%ZAERMAP(1,1))

  ALLOCATE(GEMSL%ZDMSO(1))
  ALLOCATE(GEMSL%ZLDAY(1))
  ALLOCATE(GEMSL%ZLISS(1))
  ALLOCATE(GEMSL%ZSO2(1))
  ALLOCATE(GEMSL%ZTDMS(1))
  ALLOCATE(GEMSL%ZDMSI(1))
  ALLOCATE(GEMSL%ZODMS(1))

  ALLOCATE(GEMSL%ZCAERO(1,1,1))
  ALLOCATE(GEMSL%ZAEROP(1,1,1))

  ALLOCATE(GEMSL%ZAERSRC(1,1))
  ALLOCATE(GEMSL%ZAERDDP(1,1))
  ALLOCATE(GEMSL%ZAERSDM(1,1))

  ALLOCATE(GEMSL%ZAERFLX(1,1))
  ALLOCATE(GEMSL%ZAERLIF(1,1))
  ALLOCATE(GEMSL%ZTAUAER(1,1,1))
  ALLOCATE(GEMSL%ZAERAOT(1,1,1))
  ALLOCATE(GEMSL%ZAERAOTLEV(1,1,1))
  ALLOCATE(GEMSL%ZAERAAOTLEV(1,1,1))
  ALLOCATE(GEMSL%ZAERASYLEV(1,1,1))
  ALLOCATE(GEMSL%ZAERLISI(1,1,1,1))

  ALLOCATE( GEMSL%ZDRYDIAM(1,1,1))
  ALLOCATE( GEMSL%ZWETDIAM(1,1,1))
  ALLOCATE( GEMSL%ZDRYVOL(1,1,1))
  ALLOCATE( GEMSL%ZWETVOL(1,1,1))
  ALLOCATE( GEMSL%ZVOL_WAT(1,1,1))
  ALLOCATE( GEMSL%ZRHOPAR(1,1,1))
  ALLOCATE( GEMSL%ZMD(1,1,1,1) )
  ALLOCATE( GEMSL%ZMDT(1,1,1) )
  ALLOCATE( GEMSL%ZND(1,1,1) )

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LOCAL_ARRAYS_INI',1,ZHOOK_HANDLE)
END SUBROUTINE LOCAL_ARRAYS_INI
