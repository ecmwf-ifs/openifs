! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLDPP ( YDRADIATION,YDERAD,YDERDI,YDECLD,YDRIP,KIDIA  , KFDIA , KLON   , KLEV,&
 & PAPHM1 , PAPM1 , PA     , PL     , PI     , PT , &
 & PGELAT , PGELAM , &
 & LDCLDCOVER     , PCC    , PCH   , PCL    , PCM    , PCT,&
 & PCLC   , PQIWP , PQLWP )  

!**** *CLDPP*  - COMPUTES CLOUD COVER DIAGNOSTICS

!     PURPOSE.
!     --------
!         This routine computes cloud amounts for the postprocessing
!         (Total/High/Medium/Low) from the 3D array of cloud fraction
!         using specified cloud overlap assumptions.
!         It also applies numerical security checks on the 3D arrays of
!         cloud fraction, cloud water and cloud ice contents
!         required by the radiation scheme for the computation 
!         of the fluxes. The inputs are cloud fraction and water contents
!         calculated in CLOUDSC in the previous time step.

!**   INTERFACE.
!     ----------
!              *CLDPP* IS CALLED FROM *CALLPAR*, *RADPAR* and *CALLPARAD*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!     INPUT PARAMETERS (REAL)

!    *PAPHM1*       PRESSURE ON HALF-LEVELS (T-1)                  PA
!    *PAPM1*        PRESSURE ON FULL-LEVELS (T-1)                  PA
!    *PA*           FRACTIONAL CLOUD COVER (T-1)
!    *PL*           LIQUID WATER CONTENT (T-1) (IN GRID VOLUME)   KG/KG
!    *PI*           ICE WATER CONTENT (T-1) (IN GRID VOLUME)      KG/KG
!    *PT*           TEMPERATURE (T-1)                             K
!    *PGELAT*       LATITUDE                                      RADIANS
!    *PGELAM*       LONGITUDE                                     RADIANS

!     INPUT PARAMETERS (LOGICAL)

!    *LDCLDCOVER*   IF TRUE COMPUTE CLOUD COVER DIAGNOSTICS

!     OUTPUT PARAMETERS (REAL)

!     FOR POSTPROCESSING:

!    *PCC*          CONVECTIVE CLOUD COVER (0. AT THE MOMENT)
!    *PCH*          CLOUD COVER HIGH CLOUDS
!    *PCL*          CLOUD COVER LOW CLOUDS
!    *PCM*          CLOUD COVER MID-LEVEL CLOUDS
!    *PCT*          TOTAL CLOUD COVER

!     FOR RADIATION:

!    *PCLC*         FRACTIONAL CLOUD COVER
!    *PQIWP*        CLOUD ICE CONTENT                             KG/KG
!    *PQLWP*        CLOUD WATER CONTENT                           KG/KG

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     CH. JAKOB         *ECMWF*       94-01-27

!     MODIFICATIONS.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R.Forbes      15-Oct-2007 Add MCICA generalised overlap option
!        G.Mozdzynski  30-Jun-2008 Compute cloud cover diagnostics only when needed
!        R.Forbes      01-Apr-2010 Add latitude dependence of decorr scale
!        P.Bechtold    24-Apr-2012 Use correct general parameters RG, RD
!        K. Yessad (July 2014): Move some variables.
!        R.Hogan       09-Jun-2016 Call routines from new radiation scheme
!        R.Hogan       24-Nov-2016 Compute sine of lat just for KIDIA-KFDIA
!        R.Hogan       28-Mar-2017 Use YRADIATION object
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM    ,JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK
USE YOMCT3   , ONLY : NSTEP
USE YOMRIP0  , ONLY : NINDAT
USE YOMRIP   , ONLY : TRIP
USE YOMCST   , ONLY : RPI, RG, RD
USE YOECLD   , ONLY : TECLD
USE YOERDI   , ONLY : TERDI
USE YOERAD   , ONLY : TERAD
USE RANDOM_NUMBERS_MIX  , ONLY : RANDOMNUMBERSTREAM, UNIFORM_DISTRIBUTION

! Modules from new radiation scheme
USE RADIATION_SETUP,       ONLY : TRADIATION
USE RADIATION_CLOUD_COVER, ONLY : CLOUD_COVER


!-----------------------------------------------------------------------

IMPLICIT NONE

! PARAMETER

TYPE(TRADIATION)  ,INTENT(INOUT):: YDRADIATION
TYPE(TECLD)       ,INTENT(INOUT):: YDECLD
TYPE(TERAD)       ,INTENT(INOUT):: YDERAD
TYPE(TERDI)       ,INTENT(INOUT):: YDERDI
TYPE(TRIP)        ,INTENT(INOUT):: YDRIP
REAL(KIND=JPRB), PARAMETER ::&
 & ZM2KM  = 1.0_JPRB/1000.0_JPRB,&! Convert meters to kilometers
 & ZCUT   = 0.001_JPRB              ! Cutoff for cloud amount    
 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCC(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCH(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCM(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQIWP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQLWP(KLON,KLEV) 
LOGICAL           ,INTENT(IN)    :: LDCLDCOVER


!     -----------------------------------------------------------------

!*       0.1   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JK, JL, I      ! Loop indices
INTEGER(KIND=JPIM) :: I_LOCT(KLON)   ! Indicator for total cloud cover
INTEGER(KIND=JPIM) :: I_LOCH(KLON)   ! Indicator for high cloud cover
INTEGER(KIND=JPIM) :: I_LOCM(KLON)   ! Indicator for medium cloud cover
INTEGER(KIND=JPIM) :: I_LOCL(KLON)   ! Indicator for low cloud cover
INTEGER(KIND=JPIM) :: I_TOP(KLON)    ! Index of top most cloud layer
INTEGER(KIND=JPIM) :: I_BASE(KLON)   ! Index of lowest cloud layer
INTEGER(KIND=JPIM) :: ISEED, ITIM, IDAY  ! Random number seed variables

LOGICAL :: LLMAXRAN        ! .T. = Used Maximum-Random overlap assumption 
LOGICAL :: LLH, LLL, LLM
LOGICAL :: LL_CLDT(KLON)
LOGICAL :: LL_CLDH(KLON)   
LOGICAL :: LL_CLDM(KLON)   
LOGICAL :: LL_CLDL(KLON)  

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZZETA                 ! Eta level
REAL(KIND=JPRB) :: ZZETAR                ! Reciprocal of eta level
REAL(KIND=JPRB) :: ZROG                  ! R/G
REAL(KIND=JPRB) :: ZM(KLON,KLEV-1)       ! Depth of layer (km)
REAL(KIND=JPRB) :: ZLC_CF(KLON,KLEV)     ! 3D array of decorrelation scale
REAL(KIND=JPRB) :: ZALPHA(KLON,KLEV)     ! Calculated overlap randomness
REAL(KIND=JPRB) :: ZALPHA_TMP(KLEV)      ! Tmp. array of ZALPHA(*,:) for stride-1 access
REAL(KIND=JPRB) :: ZCLOUD_FRAC(KLEV)     ! Bounded cloud fraction (0-1)
REAL(KIND=JPRB) :: ZRPI                  ! Reciprocal of Pi
REAL(KIND=JPRB) :: ZGLAT(KLON)           ! Latitude in degrees
!REAL(KIND=JPRB) :: ZGLON(KLON)           ! Longitude in degrees
REAL(KIND=JPRB) :: ZCLAT(KLON)           ! Cosine(latitude)
REAL(KIND=JPRB) :: ZSLAT(KLON)           ! Sine(latitude)
REAL(KIND=JPRB) :: ZDECORR_CF(KLON)      ! Decorrelation length for cloud frac (km)

! Random number vectors
REAL(KIND=JPRB) :: ZX(KLON)
REAL(KIND=JPRB) :: ZX1(KLON)
REAL(KIND=JPRB) :: ZX2(KLON)

! Random number vectors
REAL(KIND=JPRB) :: ZZX (KLEV)
REAL(KIND=JPRB) :: ZZX1(KLEV)
REAL(KIND=JPRB) :: ZZX2(KLEV)

TYPE(RANDOMNUMBERSTREAM) :: YL_RANDOM_STREAM

! Parameters
INTEGER(KIND=JPIM), PARAMETER :: JPX_LOC=100      ! Number of sub-columns
INTEGER(KIND=JPIM) :: JKM1

! Index of first model level to be classed as "mid" and first to be
! classed as "low", when counting down from the top
INTEGER(KIND=JPIM) :: I_FIRST_MID, I_FIRST_LOW

! Include variables for random number seed/generator
#include "fcttim.func.h"
#include "setran.intfb.h"
#include "cloud_overlap_decorr_len.intfb.h"

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CLDPP',0,ZHOOK_HANDLE)
ASSOCIATE(RDECORR_CF=>YDECLD%RDECORR_CF, REPSEC=>YDECLD%REPSEC, &
 & RETAHB=>YDECLD%RETAHB, RETAMB=>YDECLD%RETAMB, &
 & NDECOLAT=>YDERAD%NDECOLAT, NMCICA=>YDERAD%NMCICA, &
 & REPCLC=>YDERDI%REPCLC, TSTEP=>YDRIP%TSTEP, &
 & I_OVERLAP_SCHEME=>YDRADIATION%RAD_CONFIG%I_OVERLAP_SCHEME)
! ----------------------------------------------------------------------

! Cloud overlap algorithm options:
! NMcICA == 0   Old version of maximum-random cloud overlap 
! NMcICA == 1   MCICA cloud generator with maximum-random overlap
! NMcICA == 2   MCICA cloud generator with generalised overlap 


IF (LDCLDCOVER) THEN

  IF (.NOT. YDERAD%LUSEPRE2017RAD) THEN
!-----------------------------------------------------------------------
!
!      0.    Use new radiation scheme
!
!-----------------------------------------------------------------------
    ! Compute sine of latitude
    ZSLAT = 0.0_JPRB
    ZSLAT(KIDIA:KFDIA) = SIN(PGELAT(KIDIA:KFDIA))

    ! Compute overlap decorrelation length in km using the scheme
    ! described by NDECOLAT
    CALL CLOUD_OVERLAP_DECORR_LEN(YDECLD,KIDIA,KFDIA,KLON,ZSLAT, &
         &    YDERAD%NDECOLAT,PDECORR_LEN_EDGES_KM=ZDECORR_CF)
    
    ! Compute the separation of adjacent layers in km
    ZROG=RD/RG    
    DO JK = 1, KLEV-1
      DO JL = KIDIA, KFDIA
        ZZETAR = PAPM1(JL,JK+1)/PAPM1(JL,JK)
        ZM(JL,JK) = ZROG*0.5*(PT(JL,JK)+PT(JL,JK+1))*LOG(ZZETAR)*ZM2KM
      ENDDO
    ENDDO

    ! Convert to "alpha" overlap parameter 
    DO JK = 1, KLEV-1
      ZALPHA(KIDIA:KFDIA,JK) &
           &  = EXP(-ZM(KIDIA:KFDIA,JK) / ZDECORR_CF(KIDIA:KFDIA))
    ENDDO

    ! Each column analyzed in turn - this is suboptimal in terms of
    ! array access, but the cloud cover algorithms are much faster
    ! than using random number generators
    DO JL = KIDIA, KFDIA
      ! Default values
      I_FIRST_MID = 2
      I_FIRST_MID = 3
      JK = 2
      ! Locate first model layer classed as "mid-level"
      DO WHILE (JK < KLEV .AND. PAPHM1(JL,JK)/PAPHM1(JL,KLEV+1) < RETAHB)
        JK = JK+1
      ENDDO
      I_FIRST_MID = JK
      JK = JK+1
      ! Locate first model layer classed as "low-level"
      DO WHILE (JK < KLEV .AND. PAPHM1(JL,JK)/PAPHM1(JL,KLEV+1) < RETAMB)
        JK = JK+1
      ENDDO
      I_FIRST_LOW = JK

      ! Cloud fraction for just this column, bounded between 0 and 1
      ! in case of numerical issues beforehand
      ZCLOUD_FRAC(:)=MIN(MAX(PA(JL,:),0.0_JPRB),1.0_JPRB)

      ! Stride-1 access
      ZALPHA_TMP(:) = ZALPHA(JL,:)

      ! Compute total cloud cover
      PCT(JL) = CLOUD_COVER(KLEV, I_OVERLAP_SCHEME,&
           &                ZCLOUD_FRAC,&
           &                ZALPHA_TMP(1:KLEV-1))
      ! Compute the three sub cloud covers
      PCH(JL) = CLOUD_COVER(I_FIRST_MID-1, I_OVERLAP_SCHEME,&
           &                ZCLOUD_FRAC(1:I_FIRST_MID-1),&
           &                ZALPHA_TMP(1:I_FIRST_MID-2))
      PCM(JL) = CLOUD_COVER(I_FIRST_LOW-I_FIRST_MID,&
           &                I_OVERLAP_SCHEME,&
           &                ZCLOUD_FRAC(I_FIRST_MID:I_FIRST_LOW-1),&
           &                ZALPHA_TMP(I_FIRST_MID:I_FIRST_LOW-2))
      PCL(JL) = CLOUD_COVER(KLEV-I_FIRST_LOW+1, I_OVERLAP_SCHEME,&
           &                ZCLOUD_FRAC(I_FIRST_LOW:KLEV),&
           &                ZALPHA_TMP(I_FIRST_LOW:KLEV-1))
    ENDDO

  ELSEIF (NMCICA == 0) THEN

!-----------------------------------------------------------------------

!                  Computationally efficient
!      1.    Maximum-Random cloud overlap assumption  
!                (Geleyn and Hollingsworth, 1979)
!-----------------------------------------------------------------------
! This algorithm assumes MAXIMUM overlap between layers within a vertically 
! continuous cloud when there when there is no cloud fraction minimum. 
! If there is a minimum, cloud overlap has a random component.
! The algorithm assumes RANDOM overlap between vertically separated clouds.
! NOTE: This particular formulation will completely miss low and medium 
! level cloud if a cloud crosses the high/medium/low boundaries and is 
! decreasing with decreasing height !

    ! Initialise arrays
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
      PCLC(JL,1) = MIN(MAX(PA(JL,1),REPCLC),1.0_JPRB-REPCLC)
      PCT(JL)    = 1.0_JPRB-PCLC(JL,1)
      PCH(JL)    = 1.0_JPRB-PCLC(JL,1)
      PCM(JL)    = 1.0_JPRB
      PCL(JL)    = 1.0_JPRB
      PCC(JL)    = 0.0_JPRB
    ENDDO

    DO JK=2,KLEV
!DEC$ IVDEP
      DO JL=KIDIA,KFDIA
        ZZETA = PAPHM1(JL,JK)/PAPHM1(JL,KLEV+1)
        LLL = ZZETA > RETAMB
        LLM = ZZETA <= RETAMB.AND.ZZETA >= RETAHB
        LLH = ZZETA < RETAHB
        IF(PA(JL,JK-1) < 1.0_JPRB-REPSEC) THEN
          PCT(JL)=PCT(JL)*(1.0_JPRB-MAX(PA(JL,JK),PA(JL,JK-1)))    &
           & /(1.0_JPRB-MIN(PA(JL,JK-1),1.0_JPRB-REPSEC))  
  
          IF(LLH)&
           & PCH(JL)=PCH(JL)*(1.0_JPRB-MAX(PA(JL,JK),PA(JL,JK-1))) &
           & /(1.0_JPRB-MIN(PA(JL,JK-1),1.0_JPRB-REPSEC))  
  
          IF(LLM)&
           & PCM(JL)=PCM(JL)*(1.0_JPRB-MAX(PA(JL,JK),PA(JL,JK-1))) &
           & /(1.0_JPRB-MIN(PA(JL,JK-1),1.0_JPRB-REPSEC))  
  
          IF(LLL)&
           & PCL(JL)=PCL(JL)*(1.0_JPRB-MAX(PA(JL,JK),PA(JL,JK-1))) &
           & /(1.0_JPRB-MIN(PA(JL,JK-1),1.0_JPRB-REPSEC))  
        ENDIF
      ENDDO
    ENDDO 
  
    DO JL=KIDIA,KFDIA
      PCT(JL) = 1.0_JPRB-PCT(JL)
      PCH(JL) = 1.0_JPRB-PCH(JL)
      PCM(JL) = 1.0_JPRB-PCM(JL)
      PCL(JL) = 1.0_JPRB-PCL(JL)
    ENDDO
  
  ELSEIF (NMCICA == 1 .OR. NMCICA == 2) THEN
  
!-----------------------------------------------------------------------

!        2.     Generalised Cloud Overlap using McICA Generator  

!-----------------------------------------------------------------------
! Function of maximum overlap and random overlap
! Maximum-random is a special case of the generalised overlap

    IF (NMCICA == 1) LLMAXRAN = .TRUE.   ! Maximum-random overlap
    IF (NMCICA == 2) LLMAXRAN = .FALSE.  ! Generalised overlap
  
    ! ---------------------------------------------------------------
    ! Initialize cloud indicator arrays
    DO JL = KIDIA, KFDIA
      I_LOCT(JL) = 0
      I_LOCH(JL) = 0
      I_LOCM(JL) = 0
      I_LOCL(JL) = 0
      LL_CLDT(JL) = .FALSE.
      LL_CLDH(JL) = .FALSE.
      LL_CLDM(JL) = .FALSE.
      LL_CLDL(JL) = .FALSE.
    ENDDO

    ! ---------------------------------------------------------------
    ! Create a random number seed from date, time, latitude and longitude
  
    ! Get the starting day and number of minutes since start
    IDAY = NDD(NINDAT)
    ITIM = NINT(NSTEP*TSTEP/60._JPRB)

!    ! Calculate latitude and longitude of each point
!    ZRPI=1.0_JPRB/RPI
!    DO JL=KIDIA,KFDIA
!      ZGLON(JL) = PGELAM(JL)*180._JPRB*ZRPI
!      ZGLAT(JL) = PGELAT(JL)*180._JPRB*ZRPI
!    ENDDO

    ! Calculate random number seed
    ISEED = ITIM+IDAY
    CALL SETRAN ( ISEED, YL_RANDOM_STREAM )
  
!    DO JL=KIDIA,KFDIA
!      ISEED = NINT( ZGLON(JL)+ZGLAT(JL)+ITIM+IDAY )
!      CALL SETRAN ( ISEED, YL_RANDOM_STREAM(JL) )
!    ENDDO

    ! ---------------------------------------------------------------
    ! Set decorrelation depth for cloud fraction overlap
    ! e.g. Hogan and Illingworth (2000)
  
    ! The following code sets the cloud overlap decorrelation scale to 
    ! a constant for all points. The same code appears in the radiation
    ! code in routine For consistency with the radiation scheme, should define ZLC_CF
    ! (PRLC_CF in MCICA_CLD_GEN) in only one place. Then calculate 
    ! ZALPHA just once and pass in to this routine.....
    ZRPI=1.0_JPRB/RPI
    DO JL=KIDIA,KFDIA
      ! Latitude in degrees
      ZGLAT(JL)= PGELAT(JL) * 180._JPRB*ZRPI
      ! cosine(latitude in radians)
      ZCLAT(JL)= COS(PGELAT(JL))
      IF (NDECOLAT == 0) THEN
        ZDECORR_CF(JL)=RDECORR_CF
      ELSEIF (NDECOLAT == 1) THEN
        ZDECORR_CF(JL)=2.899_JPRB-0.02759_JPRB*ABS(ZGLAT(JL))
      ELSEIF (NDECOLAT == 2) THEN
        ZDECORR_CF(JL)=0.75_JPRB + 2.149_JPRB*ZCLAT(JL)*ZCLAT(JL)
      ENDIF 
    ENDDO
    
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZLC_CF(JL,JK) = ZDECORR_CF(JL)
      ENDDO
    ENDDO

    ! Compute the separation of adjacent layers
    ZROG=RD/RG    
    DO JK = 1, KLEV-1
      DO JL = KIDIA, KFDIA
        ZZETAR = PAPM1(JL,JK+1)/PAPM1(JL,JK)
        ZM(JL,JK) = ZROG*0.5*(PT(JL,JK)+PT(JL,JK+1))*LOG(ZZETAR)*ZM2KM
      ENDDO
    ENDDO

    ! Calculate overlap factors ZALPHA for cloud fraction 
    ! based on layer midpoint distances and decorrelation depths 
    ! ZM = Full-level (layer midpoint) separation  (km)
    ! ZLC_CF = cloud overlap decorrelation depth scale (km)
    DO JK=1,KLEV-1
      DO JL = KIDIA, KFDIA ! loop over all points
        ZALPHA(JL,JK) = EXP(-ZM(JL,JK) / ZLC_CF(JL,JK))
      ENDDO ! JL
    ENDDO ! JK

    ! ---------------------------------------------------------------
    ! Calculate top and bottom cloudy layers for each grid column 

    ! Find uppermost cloudy layer
    DO JL = KIDIA, KFDIA
      DO JK = 1,KLEV
        I_TOP(JL) = JK
        IF (PA(JL,JK) > ZCUT) EXIT
      ENDDO ! JK
    ENDDO ! JL


    ! Find lowermost cloudy layer
    DO JL = KIDIA, KFDIA
      DO JK=KLEV,1,-1
        I_BASE(JL) = JK
        IF (PA(JL,JK) > ZCUT) THEN
          EXIT
        ENDIF
      ENDDO ! JK
    ENDDO ! JL
  
    ! ---------------------------------------------------------------
    ! Loop over sub-columns to sample a distribution of cloud/nocloud 
    ! profiles based on probabilities.

    DO I = 1,JPX_LOC
  
      ! Set up three different arrays of uncorrelated random numbers
      CALL UNIFORM_DISTRIBUTION (ZZX(1:KLEV) , YL_RANDOM_STREAM )
      CALL UNIFORM_DISTRIBUTION (ZZX1(1:KLEV), YL_RANDOM_STREAM )
      CALL UNIFORM_DISTRIBUTION (ZZX2(1:KLEV), YL_RANDOM_STREAM )
  
      ! Generate all subcolumns for latitude chain
      DO JL = KIDIA, KFDIA
        DO JK = I_TOP(JL), I_BASE(JL)

          ! If at the top of the cloud initialise with a random number
          IF (JK == I_TOP(JL)) THEN
            ZX(JL) = ZZX(JK)
          ENDIF
  
          ZX1(JL) = ZZX1(JK)
          ZX2(JL) = ZZX2(JK)

          IF (LLMAXRAN) THEN
            ! Maximum-random overlap
            IF (ZX(JL) <= 1.0_JPRB - PA(JL,JK-1)) THEN ! It is clear above
              ZX(JL) = ZX1(JL) * (1.0_JPRB - PA(JL,JK-1))
            ENDIF
          ELSE
            ! Generalized overlap based on a decorrelation scale, whether
            ! vertically continuous or separated clouds
            ! (For random overlap between vertically separated clouds
            !  include .OR. PA(JL,JK) < REPSEC)
            IF(JK /= 1)THEN
              JKM1=JK-1
            ELSE
              JKM1=1
            ENDIF
            IF (ZX1(JL) > ZALPHA(JL,JKM1)) ZX(JL) = ZX2(JL)
          ENDIF

          ! If ZX < cloud fraction in this cell then generate cloud
          IF (ZX(JL) > 1.0_JPRB - PA(JL,JK)) THEN 
            ! Determine whether this grid point is in the High, Medium 
            ! or Low cloud category
            ZZETA = PAPHM1(JL,JK)/PAPHM1(JL,KLEV+1)
            LLL = ZZETA > RETAMB
            LLM = ZZETA <= RETAMB.AND.ZZETA >= RETAHB
            LLH = ZZETA < RETAHB
            ! Set cloudy indicators
            LL_CLDT(JL)         = .TRUE.
            IF(LLH) LL_CLDH(JL) = .TRUE.
            IF(LLM) LL_CLDM(JL) = .TRUE.
            IF(LLL) LL_CLDL(JL) = .TRUE.
          ENDIF
  
        ENDDO                ! JK
      ENDDO              ! JL
  
      ! Need to check if a cloudy subcolumn was generated and set
      ! Total/High/Medium/Low cloud indicators as appropriate
      DO JL = KIDIA, KFDIA
        IF (LL_CLDT(JL)) THEN
          I_LOCT(JL)  = I_LOCT(JL) + 1
          LL_CLDT(JL) = .FALSE.
        ENDIF
        IF (LL_CLDH(JL)) THEN
          I_LOCH(JL)  = I_LOCH(JL) + 1
          LL_CLDH(JL) = .FALSE.
        ENDIF
        IF (LL_CLDM(JL)) THEN
          I_LOCM(JL)  = I_LOCM(JL) + 1
          LL_CLDM(JL) = .FALSE.
        ENDIF
        IF (LL_CLDL(JL)) THEN
          I_LOCL(JL)  = I_LOCL(JL) + 1
          LL_CLDL(JL) = .FALSE.
        ENDIF
      ENDDO
  
    ENDDO                  ! I, number of sub-columns

    ! The number of cloudy subcolumns generated divided by 
    ! the total number of subcolumns gives the cloud fraction for
    ! each cloud category
  
    DO JL = KIDIA, KFDIA
      PCT(JL) = FLOAT(I_LOCT(JL)) / FLOAT(JPX_LOC)
      PCH(JL) = FLOAT(I_LOCH(JL)) / FLOAT(JPX_LOC)
      PCM(JL) = FLOAT(I_LOCM(JL)) / FLOAT(JPX_LOC)
      PCL(JL) = FLOAT(I_LOCL(JL)) / FLOAT(JPX_LOC)
      PCC(JL) = 0.0_JPRB
    ENDDO ! JL
  
  ENDIF ! on NMcICA

ENDIF ! LDCLDCOVER


!-----------------------------------------------------------------------

!        3.     Array security for radiation scheme

!-----------------------------------------------------------------------
! Set arrays for radiation scheme removing inconsistencies due to 
! numerical truncation, i.e.: 
!  Cloud Fraction must be between 0.0 and 1.0 and 
!  LWP and IWP are 0.0 if there is no cloud
! Note: These checks could be separated from the above cloud cover 
! diagnostic code in the future. The adjoint/TL calls to this routine 
! or the equivalent only require the following security check code.

DO JK=1,KLEV
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA
    PCLC(JL,JK)=MIN(MAX(PA(JL,JK),REPCLC),1.0_JPRB-REPCLC)
    IF(PCLC(JL,JK) > REPCLC) THEN
      PQLWP(JL,JK)=PL(JL,JK)
      PQIWP(JL,JK)=PI(JL,JK)
    ELSE
      PQLWP(JL,JK)=0.0_JPRB
      PQIWP(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO


!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLDPP',1,ZHOOK_HANDLE)
END SUBROUTINE CLDPP
