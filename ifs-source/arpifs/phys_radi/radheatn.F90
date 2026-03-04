! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADHEATN &
 & (  YDERAD,YDERDI,YDPHY2,YDPHY3,YDSPP_CONFIG,KIDIA  , KFDIA  , KLON   , KLEV,&
 & PAPHM1 ,&
 & PEMIS  , PSPECTRALEMISS, PEMTD, PMU0, PMU0M, &
 & PQM1,&
 & PTE    , PTRSOL , PTRSOD , PTRTHD, PTRSODC, PTRTHDC,PTSM1M ,& 
 & PHRSW  , PHRLW  , PHRSC  , PHRLC,&
 & PFRSO  , PFRTH  , PFRSOD , PFRSODC, PFRTHD, PFRTHDC,& 
 & PEMTEC , PTRSOC , PFRSOC , PFRTHC , PINCSR,&
 & PSUDU  , PSDUR  , PDSRP,&
 & PUVDFI , PPARFI , PPARCFI, PTINCFI, PFDIRI, PCDIRI,& 
 & PUVDF  , PPARF  , PPARCF , PTINCF , PFDIR , PCDIR,  &
 & PLWDERIVATIVE, &
 & PNEB   , PAP    , PGP2DSPP )  

!**** *RADHEATN* - COMPUTES TEMPERATURE CHANGES DUE TO RADIATION.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TENDENCIES OF THE ATMOSPHERE'S
!     TEMPERATURE DUE TO THE EFFECTS OF LONG WAVE AND SHORT WAVE
!     RADIATION. THE COMPUTATION IS DONE ON THE T-1 TIME LEVEL USING
!     VALUES OF ATMOSPHERIC TRANSMISIVITIES AND EMISSIVITIES THAT HAVE
!     BEEN STORED AT THE LAST FULL RADIATION TIME STEP. THE SURFACE
!     SOLAR FLUX LATER TO BE USED IN THE SOIL PROCESS CALCULATIONS IS
!     ALSO STORED.
!     Correction June 2014: longwave fluxes rather than emissivities
!     are used

!**   INTERFACE.
!     ----------

!          *RADHEATN* IS CALLED FROM *CALLPAR* via *RADFLUX_LAYER*.
!          THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE: TS,
!     T AND Q AT T-1 AND P ON LEVEL BOUNDARIES AND IN THE MIDDLE OF THE
!     LAYERS AT THE SAME TIME. ALSO USED ARE THE FOUR PARAMETERS
!     DESCRIBING THE SUN'S POSITION. THE ROUTINE RETURNS ITS OUTPUT TO
!     THE LONG TERM STORAGE: TENDENCIES OF T AND SURFACE SOLAR FLUX.

!     METHOD.
!     -------

!          A CALL TO SUBROUTINE *SOLANG* GIVES FIELDS OF SOLAR ZENITH
!     ANGLES AND RELATIVE DAY LENGTH FROM WHICH AN EFFECTIVE SOLAR
!     INFLUX IS COMPUTED. THE RESULTS ARE OF COURSE DIFFERENT DEPENDING
!     ON THE SWITCH ON OR OFF OF THE DIURNAL CYCLE. PRODUCT OF SOLAR
!     INFLUX BY TRANSMISSIVITIES LEADS TO SOLAR FLUXES. THEN THE
!     TEMPERATURES ARE INTERPOLATED/EXTRAPOLATED TO THE LAYER BOUNDARIES
!     (AT THE BOTTOM ONE TAKES THE SURFACE TEMPERATURE) AND A PRODUCT BY
!     EMISSIVITIES OF SIGMA*T**4 GIVES THERMAL FLUXES. THE TWO FLUXES
!     ARE ADDED AND DIVERGENCES COMPUTED TO GIVE HEATING RATES.
!     Correction June 2014: No temperature interpolation is carried
!     out; the fluxes from the LW radiation scheme are taken.  These
!     are updated based on the surface temperature using the stored
!     partial derivatives.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!          SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION FOR DETAILS
!     ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     AUTHOR.
!     -------
!      J.-F. GELEYN     E.C.M.W.F.    82/06/03.

!     MODIFICATIONS.
!     --------------
!      JJMorcrette      ECMWF      02-09-01    UV & PAR
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      JJMorcrette      20060721   PP clear-sky PAR and TOA incident solar radiation
!      JJMorcrette      20090408   bug-fix to solar zenith angle
!      JJMorcrette      20091201   Total and clear-sky direct SW radiation flux at surface 
!      M Ahlgrimm       31 Oct 2011 Surface downward clear-sky LW and SW fluxes 
!      N.Semane+P.Bechtold     04-10-2012 Add RPLDARE factor for small planet
!      R J Hogan        May/June 2014  Approximate update to LW net fluxes; added comments
!      R J Hogan        24 Oct 2014 Correction for solar zenith angle using Manners et al. (2009)
!      R J Hogan        15 Apr 2015 Use LMANNERSSWUPDATE rather than LAPPROXSWUPDATE
!      A Bozzo          Nov 2015 Computation of sunshine duration using the direct beam SW flux
!                       witch corrected path-length
!      F. Vana          17-Dec-2015  Support for single precision
!      M.Leutbecher&S.-J.Lock Jan 2016  Introduced SPP scheme (LSPP)
!      R J Hogan         5 Feb 2019 Spectral approximate longwave updates
!      M Leutbecher      Oct 2020 SPP abstraction
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHY2  , ONLY : TPHY2
USE YOMPHY3  , ONLY : TPHY3
USE YOMCST   , ONLY : RG, RSIGMA, RCPD
USE YOETHF   , ONLY : RVTMP2
USE YOERDI   , ONLY : TERDI
USE YOMDYNCORE,ONLY : LAPE, RPLDARE
USE YOERAD   , ONLY : TERAD
USE SPP_MOD     , ONLY : TSPP_CONFIG
USE SPP_GEN_MOD , ONLY : SPP_PERT

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TERDI)       ,INTENT(INOUT) :: YDERDI
TYPE(TPHY2)       ,INTENT(INOUT) :: YDPHY2
TYPE(TPHY3)       ,INTENT(INOUT) :: YDPHY3
TYPE(TSPP_CONFIG) ,INTENT(IN)    :: YDSPP_CONFIG
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 

! Inputs based on current conditions
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,KLEV+1) ! Half-level pressure (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV)     ! Specific humidity (kg kg-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSM1M(KLON)        ! Skin temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON)          ! Instantaneous cosine of solar zenith angle
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0M(KLON)         ! Cosine of solar angle used by radiation scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON)         ! Surface longwave emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPECTRALEMISS(KLON,YDERAD%NLWOUT) ! Surface longwave emissivity

! Inputs from previous call of the radiation scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMTD(KLON,KLEV+1)  ! Net longwave flux (W m-2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRSOL(KLON,KLEV+1) ! Net shortwave transmissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMTEC(KLON,KLEV+1) ! Clear-sky net LW flux (W m-2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRSOC(KLON,KLEV+1) ! Clear-sky net SW transmissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRSOD(KLON)  ! Surface downwelling shortwave
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRTHD(KLON,YDERAD%NLWOUT) ! Surface downwelling longwave in emissivity bands (W m-2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRSODC(KLON) ! Surface downwelling clear-sky shortwave
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRTHDC(KLON) ! Surface downwelling clear-sky longwave (W m-2)

! Outputs, some updated based on current conditions
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTE(KLON,KLEV)   ! Temperature tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHRSW(KLON,KLEV) ! Shortwave heating rate
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHRLW(KLON,KLEV) ! Longwave heating rate
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHRSC(KLON,KLEV) ! Shortwave clear-sky heating rate
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHRLC(KLON,KLEV) ! Longwave clear-sky heating rate
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSO(KLON,KLEV+1) ! Net SW accounting for sun angle
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTH(KLON,KLEV+1) ! Net LW updated for surf temperature
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOD(KLON)   ! Surface shortwave downwelling flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSODC(KLON)  ! Clear-sky surface SW downwelling
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTHD(KLON)   ! Surface longwave downwelling flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTHDC(KLON)  ! Clear-sky surface LW downwelling
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOC(KLON,2) ! Clear-sky net SW flux 1=TOA 2=surf
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTHC(KLON,2) ! Clear-sky net LW flux 1=TOA 2=surf
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PINCSR(KLON)   ! SW TOA flux accounting for zenith angle

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSUDU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDUR(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDSRP(KLON)    ! Output surface direct beam shortwave flux (orthogonal to Sun's direction) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUVDFI(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPARFI(KLON)   ! Input photosynthetically active radiation
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPARCFI(KLON), PTINCFI(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVDF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARCF(KLON), PTINCF(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFDIRI(KLON), PCDIRI(KLON) ! Input surface direct shortwave transmissivity (total-/clear-sky)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDIR (KLON), PCDIR (KLON) ! Output surface direct shortwave flux (total-/clear-sky)

! Partial derivatives of upwelling longwave to the surface value; note
! that original dimensioning was (KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLWDERIVATIVE(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PNEB(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PGP2DSPP(KLON, YDSPP_CONFIG%SM%NRFTOTAL)
!-----------------------------------------------------------------------

REAL(KIND=JPRB) ::&
 &    ZFLB(KLON) , ZFLT(KLON) , ZFSO(KLON) , ZFTE(KLON)&
 & ,  ZFSWT(KLON), ZFSWB(KLON), ZFLWT(KLON), ZFLWB(KLON)&
 & ,  ZFSWTC(KLON),ZFSWBC(KLON),ZFLWTC(KLON),ZFLWBC(KLON)&
 & ,  ZI0(KLON)

INTEGER(KIND=JPIM) :: ILONS, JK, JL
INTEGER(KIND=JPIM) :: IPPHRLW, IPPHRSW

REAL(KIND=JPRB) :: ZBUD, ZCONS1, ZCONS3, ZRII0, ZTMST

! Factor to convert from the flux difference across a layer to the
! rate of change of temperature
REAL(KIND=JPRB) :: ZFLUXTOHEATINGRATE

! Difference in surface upwelling and downwelling longwave fluxes
! between the current timestep/gridpoint and the one from the most
! recent call to the radiation scheme; these can be used to update the
! longwave net flux profile
REAL(KIND=JPRB) :: ZLWUPSURFDIFF(KLON)
REAL(KIND=JPRB) :: ZLWDNSURFDIFF(KLON)

! Planck function at surface skin temperature integrated across
! surface emissivity intervals
REAL(KIND=JPRB) :: ZPLANCK(KLON,YDERAD%NLWOUT)

! We assume that due to the coupling of the near-surface air
! temperature to the surface on a rapid timescale, and the opaqueness
! of the lower atmosphere, that in addition to updating the upwelling
! longwave flux according to skin temperature, we need to change the
! downwelling flux in proportion but by a smaller amount: this
! constant is the proportion.
REAL(KIND=JPRB), PARAMETER :: ZLWDNDIFFRATIO = 0.2

! The above implies that there will be a dependence of the change in
! downwelling flux at each height on the change in the surface
! downwelling flux
REAL(KIND=JPRB) :: ZLWDNDERIVATIVE

! Cosine of solar zenith angles accounting for earth curvature
REAL(KIND=JPRB) :: ZMU0CURV, ZMU0MCURV
! Ratio of the old cosine of the solar zenith angle to the new one
REAL(KIND=JPRB) :: ZMU0RATIO
! New downwelling direct shortwave transmissivity at the surface,
! accounting for the correction for path length associated with
! changed solar zenith angle
REAL(KIND=JPRB) :: ZSWDIRECTTRANS

! Change to net shortwave flux (W m-2) at all model levels due to
! correction for path length associated with changed solar zenith
! angle (and clear-sky value)
REAL(KIND=JPRB) :: ZSWNETDIFF(KLON), ZSWNETDIFFC(KLON)
! Scaling to be applied to surface total (net or downward) and direct
! downwelling fluxes associated with correcting the path length due to
! the changed solar zenith angle (and clear-sky values)
REAL(KIND=JPRB) :: ZSWSURFSCALING(KLON),  ZSWDIRECTSCALING(KLON)
REAL(KIND=JPRB) :: ZSWSURFSCALINGC(KLON), ZSWDIRECTSCALINGC(KLON)

! Factor associated with Earth curvature correction
REAL(KIND=JPRB) :: ZCRAE

! Aux variables for SPP
REAL(KIND=JPRB)    :: ZPRESS
LOGICAL            :: LLPERT_PHR
INTEGER(KIND=JPIM) :: IPN
INTEGER(KIND=JPIM) :: IDXSTTR
TYPE(SPP_PERT)     :: PN1

REAL(KIND=JPRB) :: ZEPSILON, ZEPSILON1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------

IF (LHOOK) CALL DR_HOOK('RADHEATN',0,ZHOOK_HANDLE)
ASSOCIATE(LAPPROXLWUPDATE=>YDERAD%LAPPROXLWUPDATE, &
 & LMANNERSSWUPDATE=>YDERAD%LMANNERSSWUPDATE, &
 & RRAE=>YDERDI%RRAE, RSUNDUR=>YDERDI%RSUNDUR, &
 & TSPHY=>YDPHY2%TSPHY, RII0=>YDPHY3%RII0)
ZTMST=TSPHY
ZCONS1=RRAE*(RRAE+2.0_JPRB)
ZCONS3=(RG*RPLDARE)/RCPD
ILONS=KFDIA-KIDIA+1
ZBUD = 1.0_JPRB
ZRII0=RII0

!-- Prepare SPP heating rate perturbation settings
IF (YDSPP_CONFIG%LSPP) THEN
  IPN = YDSPP_CONFIG%PPTR%PHR
  LLPERT_PHR = IPN > 0
  IF (LLPERT_PHR) THEN
    PN1 = YDSPP_CONFIG%SM%PN(IPN)
    IPPHRLW=PN1%MP 
    IPPHRSW=PN1%MP + 1 
  ENDIF
ELSE
  LLPERT_PHR=.FALSE.
ENDIF
ZEPSILON=100._JPRB*EPSILON(ZEPSILON)
ZEPSILON1=10000._JPRB*EPSILON(ZEPSILON)  ! Bit stronger criterion
!     ------------------------------------------------------------------

!*         1.     SECURITIES.
!                 -----------

!     ------------------------------------------------------------------

!*         2.     INCIDENT SOLAR RADIATION
!                 ------------------------

! To convert shortwave transmissivities to fluxes we multiply them by
! ZI0, which contains the total solar irradiance (accounting for
! sun-earth distance) multiplied by the cosine of the solar zenith
! angle.  Here ZI0 is calculated.
DO JL=KIDIA,KFDIA
  IF (PMU0(JL) >= ZEPSILON) THEN
    IF( LAPE ) THEN
      ZI0(JL)=ZRII0*PMU0(JL)
    ELSE
!!      ZI0(JL)=RRAE/(SQRT(PMU0(JL)**2+ZCONS1)-PMU0(JL)) * ZRII0
      ZI0(JL)=ZRII0*PMU0(JL)
    ENDIF
  PDSRP(JL)=1.0_JPRB
  ELSE
    ZI0(JL)=0.0_JPRB
    PDSRP(JL)=0.0_JPRB
  ENDIF
  PINCSR(JL)=ZI0(JL)
ENDDO

!     ------------------------------------------------------------------

!*         3.     TEMPERATURES AT LAYER BOUNDARIES.
!                 -----------------------------------

!     ------------------------------------------------------------------

!*         4.     FLUXES AND THEIR DIVERGENCES.
!                 ------ --- ----- ------------

!*         4.0.1     OPTIONAL UPDATE TO LONGWAVE NET FLUXES

! Longwave net flux PEMTD entering this subroutine is directly from
! the last call of the radiation scheme. Longwave net flux PFRTH
! exiting this subroutine is used by the surface scheme, and in this
! subroutine to compute heating rate profiles.  We may either copy
! PEMTD to PFRTH, or make an approximate update to the longwave net
! flux that accounts for the different surface temperature compared to
! that seen by the radiation scheme, but does not account for the
! different atmospheric temperature.
IF (LAPPROXLWUPDATE) THEN
  ! We want the change to the upwelling longwave flux, defined as
  ! LwUpSurfDiff = LwUp_now - LwUp_rad (where "rad" denotes last
  !                                     call to radiation scheme)
  ! LwUp_now = emissivity*sigma*T^4 + (1-emissivity)*LwDn_rad
  ! LwUp_rad = - LwNet_rad + LwDn_rad
  IF (YDERAD%NLWOUT == 1) THEN
    ! Use the broadband emissivity - this is not as accurate if there
    ! is significant spectral variation in emissivity
    DO JL=KIDIA,KFDIA
      ! Substituting into the above formula for LwUpSurfDiff we get:
      ZLWUPSURFDIFF(JL) = PEMTD(JL,KLEV+1) &
         & + PEMIS(JL)*(RSIGMA*PTSM1M(JL)**4 - PTRTHD(JL,1))
    ENDDO
  ELSE
    ! Use the spectral emissivity.  First compute the Planck function
    ! in each emissivity interval.
    CALL YDERAD%YSPECTPLANCK%CALC(KIDIA,KFDIA,KLON,PTSM1M,ZPLANCK)
    ! Now apply the same formula above but using the spectral
    ! emissivity and downwelling flux
    ZLWUPSURFDIFF(KIDIA:KFDIA) = PEMTD(KIDIA:KFDIA,KLEV+1)
    DO JK=1,YDERAD%NLWOUT
      DO JL=KIDIA,KFDIA
        ZLWUPSURFDIFF(JL) = ZLWUPSURFDIFF(JL) &
             & + PSPECTRALEMISS(JL,JK) * (ZPLANCK(JL,JK)-PTRTHD(JL,JK))
      ENDDO
    ENDDO
  ENDIF

  ! Assume longwave downwelling responds in a fixed ratio to the
  ! change in longwave upwelling
  ZLWDNSURFDIFF(KIDIA:KFDIA) = ZLWDNDIFFRATIO * ZLWUPSURFDIFF(KIDIA:KFDIA)

  ! Update net longwave flux at every half-level
  DO JK=1,KLEV+1
    DO JL=KIDIA,KFDIA
       ! The partial derivative of the change in downwelling longwave
       ! at half-level JL with respect to the change in surface
       ! downwelling longwave has the same shape as the equivalent
       ! upwelling partial derivative, but must go to zero at the top
       ! of atmosphere.
       ZLWDNDERIVATIVE = (PLWDERIVATIVE(JL,JK) - PLWDERIVATIVE(JL,1))&
            &          / (1.0 - PLWDERIVATIVE(JL,1))
       ! Update net flux
       PFRTH(JL,JK) = PEMTD(JL,JK)&
            & + ZLWDNSURFDIFF(JL) * ZLWDNDERIVATIVE&
            & - ZLWUPSURFDIFF(JL) * PLWDERIVATIVE(JL,JK)

       ! Alternatively we could assume net flux change is constant
       ! with height, so longwave atmospheric heating rate is
       ! unchanged, but this has been found to lead to a worse surface
       ! temperature forecast.  It would be done with:
       ! PFRTH(JL,JK) = PEMTD(JL,JK) + LwDnSurfDiff(JL) - LwUpSurfDiff(JL) 
    ENDDO
  ENDDO
  ! Update surface longwave downwelling flux
  DO JL=KIDIA,KFDIA
     PFRTHD(JL) = SUM(PTRTHD(JL,:),1) + ZLWDNSURFDIFF(JL)
  ENDDO
ELSE
   ! No approximate longwave update: use array-wise assignment to copy
   ! input to output for both net and downwelling
   PFRTH (KIDIA:KFDIA,:) = PEMTD(KIDIA:KFDIA,:)
   PFRTHD(KIDIA:KFDIA)   = SUM(PTRTHD(KIDIA:KFDIA,:),2)
ENDIF

!*         4.0.2     OPTIONAL UPDATE TO SHORTWAVE FLUXES

! Updating shortwave fluxes to account for local surface albedo has
! already been done in RADINTG.  The first-order treatment of varying
! solar zenith angle is accounted for by carrying shortwave
! "transmissivities" through to this subroutine and then multiplying
! them by the incoming TOA flux to get shortwave fluxes.  However,
! this does not take into account the change in path-length for the
! direct incoming radiation.  Following Manners et al. (2009, QJRMS
! 135, 457-468), we compute the change to direct-beam solar radiation
! and assume that any excess is scattered half to the surface and half
! to space.
IF (LMANNERSSWUPDATE) THEN
  ! Factor to use in accounting for Earth curvature
  ZCRAE=RRAE*(RRAE+2.0_JPRB)

  DO JL=KIDIA,KFDIA
    ! Account for Earth curvature in same way as in RADINTG
    IF (PMU0(JL) >= ZEPSILON1) THEN
      ZMU0CURV = RRAE / (SQRT(PMU0(JL)**2 + ZCRAE) - PMU0(JL))
    ELSE
      ZMU0CURV = RRAE / SQRT(ZCRAE)
    ENDIF
    IF (PMU0M(JL) >= ZEPSILON1) THEN
      ZMU0MCURV = RRAE / (SQRT(PMU0M(JL)**2 + ZCRAE) - PMU0M(JL))
    ELSE
      ZMU0MCURV = RRAE / SQRT(ZCRAE)
    ENDIF
    ! MU0 is never zero so we don't need to test for this
    ZMU0RATIO = ZMU0MCURV / ZMU0CURV

    ! PFDIRI is the input direct-beam transmissivity,
    ! i.e. exp(-opt_depth/MU0).
    IF (PMU0(JL) >= ZEPSILON1 .AND. PFDIRI(JL) >= ZEPSILON1&
         & .AND. PTRSOD(JL) >= ZEPSILON1) THEN
      ZSWDIRECTTRANS = PFDIRI(JL)**ZMU0RATIO
      ! The first line of the following statement provides the change
      ! to downwelling flux, while the second line then scales this
      ! with the ratio of net to downwelling (equal to 1 - albedo) to
      ! get the change to net flux
      ZSWNETDIFF(JL) = 0.5_JPRB * ZI0(JL) * (ZSWDIRECTTRANS - PFDIRI(JL))&
           & * PTRSOL(JL,KLEV+1) / PTRSOD(JL)

      ZSWDIRECTSCALING(JL) = ZSWDIRECTTRANS / PFDIRI(JL)
      ZSWSURFSCALING(JL) = 1.0_JPRB&
           & + 0.5_JPRB*(ZSWDIRECTTRANS - PFDIRI(JL)) / PTRSOD(JL)
    ELSE
      ZSWNETDIFF(JL)        = 0.0_JPRB
      ZSWDIRECTSCALING(JL)  = 1.0_JPRB
      ZSWSURFSCALING(JL)    = 1.0_JPRB
    ENDIF

    ! Now clear-sky equivalent
    IF (PMU0(JL) >= ZEPSILON1 .AND. PCDIRI(JL) >= ZEPSILON1&
         & .AND. PTRSODC(JL) >= ZEPSILON1) THEN
      ZSWDIRECTTRANS = PCDIRI(JL)**ZMU0RATIO
      ZSWNETDIFFC(JL) = 0.5_JPRB * ZI0(JL) * (ZSWDIRECTTRANS - PCDIRI(JL))&
           & * PTRSOC(JL,KLEV+1) / PTRSODC(JL)
      ZSWDIRECTSCALINGC(JL) = ZSWDIRECTTRANS / PCDIRI(JL)
      ZSWSURFSCALINGC(JL) = 1.0_JPRB&
           & + 0.5_JPRB*(ZSWDIRECTTRANS - PCDIRI(JL)) / PTRSODC(JL)
    ELSE
      ZSWNETDIFFC(JL)       = 0.0_JPRB
      ZSWDIRECTSCALINGC(JL) = 1.0_JPRB
      ZSWSURFSCALINGC(JL)   = 1.0_JPRB
    ENDIF
  ENDDO
ELSE
  ! No correction for solar zenith angle is requested
  ZSWNETDIFF        = 0.0_JPRB
  ZSWSURFSCALING    = 1.0_JPRB
  ZSWDIRECTSCALING  = 1.0_JPRB
  ZSWNETDIFFC       = 0.0_JPRB
  ZSWSURFSCALINGC   = 1.0_JPRB
  ZSWDIRECTSCALINGC = 1.0_JPRB
ENDIF


!*         SCALE SHORTWAVE TRANSMISSIVITIES TO GET NET FLUXES

! Copy inputs to outputs but multiplying by TOA incoming flux
DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    PFRSO(JL,JK) = ZI0(JL)*PTRSOL(JL,JK)   + ZSWNETDIFF(JL)  ! Net flux at all levels
    PFRSOC(JL,1) = ZI0(JL)*PTRSOC(JL,1)    + ZSWNETDIFFC(JL) ! ...clear-sky TOA
    PFRSOC(JL,2) = ZI0(JL)*PTRSOC(JL,KLEV+1)+ZSWNETDIFFC(JL) ! ...clear-sky surface
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PFRSOD(JL)  = ZI0(JL)*PTRSOD(JL) *ZSWSURFSCALING(JL)  ! Surface downwelling
  PFRSODC(JL) = ZI0(JL)*PTRSODC(JL)*ZSWSURFSCALINGC(JL) ! ...clear-sky value
ENDDO


! Set surface solar flux to minimum value for stability in vertical
! diffusion
DO JL=KIDIA,KFDIA
  IF (PFRSO(JL,KLEV+1) < ZEPSILON) THEN
    PFRSO(JL,KLEV+1)   = ZEPSILON
  ENDIF
ENDDO

!*         4.1     TOP FLUX.

DO JL=KIDIA,KFDIA
  ZFSO(JL)  = PFRSO(JL,1)        ! The updated net shortwave flux
  ZFTE(JL)  = PFRTH(JL,1)        ! The updated net longwave flux
  ZFSWT(JL) = ZFSO(JL)           ! Net shortwave at top of layer
  ZFLWT(JL) = ZFTE(JL)           ! Net longwave at top of layer
  ZFLT(JL)  = ZFSO(JL)+ZFTE(JL)  ! Net radiation (SW+LW) at top of layer

  ! Clear-sky equivalent quantities           

  ! The following is unmodified by ZSwNetDiffC but since net flux is
  ! shifted by a constant, the clear-sky heating rate is unchanged
  ZFSWTC(JL)   = ZI0(JL) * PTRSOC(JL,1)

  ZFLWTC(JL)   = PEMTEC(JL,1)
  PFRTHC(JL,1) = PEMTEC(JL,1)
ENDDO

!*         4.2     VERTICAL LOOP, BOTTOM FLUX AND HEATING RATE.

!***
DO JK=1,KLEV
!***
  DO JL=KIDIA,KFDIA
    ZFSO(JL)  = PFRSO(JL,JK+1)   ! Updated net shortwave flux
    ZFTE(JL)  = PFRTH(JL,JK+1)   ! Updated net longwave flux
    ZFSWB(JL) = ZFSO(JL)         ! Net shortwave flux at base of layer
    ZFLWB(JL) = ZFTE(JL)         ! Net longwave flux at base of layer
    ZFLB(JL)  = ZFSO(JL)+ZFTE(JL)! Net radiation (SW+LW) at base of layer
    ZFLUXTOHEATINGRATE = -ZCONS3/((PAPHM1(JL,JK+1)&
         & -PAPHM1(JL,JK))*(1.0_JPRB+RVTMP2*PQM1(JL,JK)))
    ! Increment temperature tendency:
    PTE(JL,JK)   = PTE(JL,JK) + ZFLUXTOHEATINGRATE*(ZFLB(JL)-ZFLT(JL))
    PHRSW(JL,JK) = ZFLUXTOHEATINGRATE*(ZFSWB(JL)-ZFSWT(JL)) ! SW heating rate
    PHRLW(JL,JK) = ZFLUXTOHEATINGRATE*(ZFLWB(JL)-ZFLWT(JL)) ! LW heating rate

    !-- SPP perturbations applied to heating rates:
    IF (LLPERT_PHR) THEN
      IF (PNEB(JL,JK) < 0.01_JPRB) THEN   !applied in clear-skies only
        ZPRESS=PAP(JL,JK)
        IF (ZPRESS < YDSPP_CONFIG%XPRESS_PHR_ST) THEN
          IDXSTTR=0
        ELSE
          IDXSTTR=1
        ENDIF
        PHRLW(JL,JK)=PHRLW(JL,JK)*EXP(PN1%MU(1+IDXSTTR)+PN1%XMAG(1+IDXSTTR)*PGP2DSPP(JL,IPPHRLW))
        PHRSW(JL,JK)=PHRSW(JL,JK)*EXP(PN1%MU(3+IDXSTTR)+PN1%XMAG(3+IDXSTTR)*PGP2DSPP(JL,IPPHRSW))
      ENDIF
    ENDIF

    ! Clear-sky quantities
    ZFSWBC(JL)   = ZI0(JL) * PTRSOC(JL,JK+1)
    ZFLWBC(JL)   = PEMTEC(JL,JK+1)
    PHRSC(JL,JK) = ZFLUXTOHEATINGRATE * (ZFSWBC(JL)-ZFSWTC(JL))
    PHRLC(JL,JK) = ZFLUXTOHEATINGRATE * (ZFLWBC(JL)-ZFLWTC(JL)) 

  ENDDO

!*         4.4     FLUX SWAP, END OF THE LOOP AND SURFACE SOLAR FLUX.

  DO JL=KIDIA,KFDIA
    ZFLT(JL)   = ZFLB(JL)
    ZFSWT(JL)  = ZFSWB(JL)
    ZFLWT(JL)  = ZFLWB(JL)
    ! Clear-sky quantities           
    ZFSWTC(JL) = ZFSWBC(JL)
    ZFLWTC(JL) = ZFLWBC(JL)
  ENDDO
!***
ENDDO
!***

DO JL=KIDIA,KFDIA
  PFRTHC(JL,2)=PEMTEC(JL,KLEV+1)
ENDDO

!*         4.5     SOLAR AND TERRESTRIAL DOWNWARD FLUXES AT SURFACE

!DEC$ IVDEP
DO JL=KIDIA,KFDIA
  PFRTHDC(JL)= PTRTHDC(JL)        
  PUVDF (JL) = ZI0(JL) * PUVDFI(JL) * ZSWSURFSCALING(JL)
  PPARF (JL) = ZI0(JL) * PPARFI(JL) * ZSWSURFSCALING(JL)
  PPARCF(JL) = ZI0(JL) * PPARCFI(JL)* ZSWSURFSCALINGC(JL)
  PTINCF(JL) = ZI0(JL) * PTINCFI(JL)
  PFDIR (JL) = ZI0(JL) * PFDIRI(JL) * ZSWDIRECTSCALING(JL)
  PCDIR (JL) = ZI0(JL) * PCDIRI(JL) * ZSWDIRECTSCALINGC(JL)
  PDSRP (JL)=PDSRP(JL)*ZRII0*PFDIRI(JL)*ZSWDIRECTSCALING(JL)
  
  !computation of the sunshine duration uses now the SW direct beam 
  !with the corrected path-length
  IF (PDSRP(JL) >= RSUNDUR) THEN
    PSDUR(JL)=1.0_JPRB
  ELSE
    PSDUR(JL)=0.0_JPRB
  ENDIF

! In previous versions of the code, surface downwelling longwave
! PFRTHD was estimated as a diagnostic here, via 
!   ZTH4(JL)=RSIGMA*PTSM1M(JL)**4
!   PFRTHD(JL)=(PEMTD(JL,KLEV+1)+PEMIS(JL)*ZTH4(JL))/PEMIS(JL)
! However, the result was unsatisfactory because surface temperature
! PTSM1M changed between radiation calls but PEMTD did not.  In the
! current version of the code, PFRTHD is assigned earlier in this
! subroutine using the value actually calculated in the radiation
! scheme.  This is done even if LApproxLwUpdate=FALSE, but since
! PFRTHD is only a diagnostic, it does not affect the evolution of the
! model.

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADHEATN',1,ZHOOK_HANDLE)
END SUBROUTINE RADHEATN
