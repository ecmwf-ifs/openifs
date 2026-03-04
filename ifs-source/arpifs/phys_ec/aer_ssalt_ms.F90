! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SSALT_MS &
  &( KIDIA, KFDIA, KLON, &
  &  PCI  , PLSM , PCLK, PWIND, &
  &  PFLXSS & 
  &)

!*** * AER_SSALT_MS* - SOURCE TERMS FOR SEA SALT AEROSOLS

!**   INTERFACE.
!     ----------
!          *AER_SSALT* IS CALLED FROM *AER_SRC*.

!     AUTHOR.
!     -------
!        original version
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        modified version:
!        Michael Schulz, LSCE/Gif-sur-Yvette

!     SOURCE.
!     -------

!    Modification to original ss routine based on O. Boucher: 
!     Tabulated wet sea salt mass fluxes at 80% RH merging Monahan 86 and
!     Smith&Harrison 98; the fluxes are integrated for three size bins according
!     to work published by Guelle et al. (2001) JGR 106, pp. 27509-27524.
!     Similar source functions have been proposed by Vignati et al 2001, JGR

!     MODIFICATIONS.
!     --------------
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEAERATM ,ONLY : YREAERATM
!USE YOEAERSRC ,ONLY : RSSFLX  

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLON), PLSM(KLON), PCLK(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWIND(KLON) 

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFLXSS(KLON,3)


!--- precalculated mass flux per wind class of 0.5 m/s
!--------unit kg m**-2 s-1 

! fine bin 0.03 - 0.5 um wet radius
REAL(KIND=JPRB), PARAMETER :: ZMASS1FLUX(0:80) = (/ &
 & 0.000E+00, 5.386E-16, 5.725E-15, 2.282E-14, 6.085E-14, 1.302E-13, 2.425E-13, &
 & 4.102E-13, 6.468E-13, 9.665E-13, 1.384E-12, 1.916E-12, 2.578E-12, 3.387E-12, &
 & 4.361E-12, 5.517E-12, 6.875E-12, 8.454E-12, 1.027E-11, 1.235E-11, 1.471E-11, &
 & 1.738E-11, 2.037E-11, 2.370E-11, 2.740E-11, 3.149E-11, 3.600E-11, 4.094E-11, &
 & 4.635E-11, 5.224E-11, 5.864E-11, 6.558E-11, 7.308E-11, 8.117E-11, 8.986E-11, &
 & 9.920E-11, 1.092E-10, 1.199E-10, 1.313E-10, 1.435E-10, 1.564E-10, 1.702E-10, &
 & 1.847E-10, 2.002E-10, 2.165E-10, 2.337E-10, 2.519E-10, 2.711E-10, 2.913E-10, &
 & 3.125E-10, 3.348E-10, 3.581E-10, 3.827E-10, 4.083E-10, 4.352E-10, 4.633E-10, &
 & 4.927E-10, 5.233E-10, 5.553E-10, 5.886E-10, 6.234E-10, 6.595E-10, 6.971E-10, &
 & 7.362E-10, 7.768E-10, 8.190E-10, 8.628E-10, 9.081E-10, 9.552E-10, 1.004E-09, &
 & 1.054E-09, 1.107E-09, 1.161E-09, 1.217E-09, 1.274E-09, 1.334E-09, 1.396E-09, &
 & 1.459E-09, 1.525E-09, 1.593E-09, 1.663E-09 /)

! middle bin 0.5 - 5 um wet radius
REAL(KIND=JPRB), PARAMETER :: ZMASS2FLUX(0:80) = (/ &
 & 0.000E+00, 4.428E-14, 4.727E-13, 1.888E-12, 5.046E-12, 1.082E-11, 2.016E-11, &
 & 3.415E-11, 5.389E-11, 8.059E-11, 1.155E-10, 1.600E-10, 2.154E-10, 2.831E-10, &
 & 3.647E-10, 4.617E-10, 5.756E-10, 7.081E-10, 8.609E-10, 1.036E-09, 1.234E-09, &
 & 1.458E-09, 1.709E-09, 1.989E-09, 2.301E-09, 2.645E-09, 3.025E-09, 3.441E-09, &
 & 3.897E-09, 4.393E-09, 4.933E-09, 5.518E-09, 6.150E-09, 6.833E-09, 7.566E-09, &
 & 8.354E-09, 9.199E-09, 1.010E-08, 1.107E-08, 1.209E-08, 1.319E-08, 1.435E-08, &
 & 1.558E-08, 1.688E-08, 1.826E-08, 1.972E-08, 2.126E-08, 2.288E-08, 2.459E-08, &
 & 2.639E-08, 2.827E-08, 3.025E-08, 3.233E-08, 3.450E-08, 3.678E-08, 3.916E-08, &
 & 4.165E-08, 4.425E-08, 4.696E-08, 4.978E-08, 5.272E-08, 5.579E-08, 5.898E-08, &
 & 6.229E-08, 6.574E-08, 6.932E-08, 7.303E-08, 7.688E-08, 8.088E-08, 8.501E-08, &
 & 8.930E-08, 9.374E-08, 9.833E-08, 1.031E-07, 1.080E-07, 1.131E-07, 1.183E-07, &
 & 1.237E-07, 1.293E-07, 1.350E-07, 1.410E-07 /)

! coarse bin 5 - 20 um wet radius
REAL(KIND=JPRB), PARAMETER :: ZMASS3FLUX(0:80) = (/ &
 & 0.000E+00, 2.842E-13, 2.688E-12, 1.014E-11, 2.619E-11, 5.483E-11, 1.005E-10, &
 & 1.680E-10, 2.624E-10, 3.892E-10, 5.541E-10, 7.632E-10, 1.023E-09, 1.339E-09, &
 & 1.719E-09, 2.169E-09, 2.697E-09, 3.311E-09, 4.017E-09, 4.823E-09, 5.739E-09, &
 & 6.771E-09, 7.928E-09, 9.218E-09, 1.065E-08, 1.224E-08, 1.398E-08, 1.589E-08, &
 & 1.799E-08, 2.027E-08, 2.275E-08, 2.543E-08, 2.834E-08, 3.147E-08, 3.484E-08, &
 & 3.846E-08, 4.234E-08, 4.649E-08, 5.092E-08, 5.564E-08, 6.066E-08, 6.600E-08, &
 & 7.167E-08, 7.767E-08, 8.401E-08, 9.072E-08, 9.780E-08, 1.053E-07, 1.131E-07, &
 & 1.214E-07, 1.301E-07, 1.392E-07, 1.488E-07, 1.588E-07, 1.693E-07, 1.803E-07, &
 & 1.918E-07, 2.038E-07, 2.163E-07, 2.293E-07, 2.429E-07, 2.571E-07, 2.718E-07, &
 & 2.872E-07, 3.031E-07, 3.197E-07, 3.368E-07, 3.547E-07, 3.732E-07, 3.924E-07, &
 & 4.122E-07, 4.328E-07, 4.541E-07, 4.761E-07, 4.989E-07, 5.224E-07, 5.467E-07, &
 & 5.718E-07, 5.978E-07, 6.245E-07, 6.521E-07 /)

REAL(KIND=JPRB), SAVE :: ZDMASS1(0:80),ZDMASS2(0:80),ZDMASS3(0:80) = 0.
LOGICAL, SAVE :: LLENTERED = .FALSE.
INTEGER(KIND=JPIM) :: IWCL(KLON)
REAL(KIND=JPRB) :: ZDZSPEED(KLON)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZFROC, ZOCEA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SSALT_MS',0,ZHOOK_HANDLE)
ASSOCIATE(RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC)

IF (.NOT. LLENTERED) THEN
  ZDMASS1(0:79) = ZMASS1FLUX(1:80)-ZMASS1FLUX(0:79) 
  ZDMASS2(0:79) = ZMASS2FLUX(1:80)-ZMASS2FLUX(0:79) 
  ZDMASS3(0:79) = ZMASS3FLUX(1:80)-ZMASS3FLUX(0:79) 
  LLENTERED = .TRUE.
ENDIF

! windclass computation, tabulation of mass fluxes is linear, equidistant f(u)
DO JL=KIDIA,KFDIA   
  IWCL(JL)=INT(PWIND(JL)*2)
  IWCL(JL)=MAX(0,MIN(IWCL(JL),80))
  ZDZSPEED(JL)=PWIND(JL)*2-IWCL(JL)
ENDDO

PFLXSS = 0._JPRB
DO JL=KIDIA,KFDIA
  ! LSM treats lakes as ocean, but (with a few exceptions e.g. the Dead Sea) they don't produce sea salt.
  ! This logic should work for either a fractional LSM or a binary one
  IF (PLSM(JL) == 0._JPRB) THEN
    ! Even with a fractional LSM, this must be all ocean or all lake, because there would be no land to divide them.
    ! Thresholding ensures sensible behaviour if LSM is binary but lake cover fractional.
    IF (PCLK(JL) >= 0.5_JPRB) THEN
      ZOCEA= 0._JPRB
    ELSE
      ZOCEA= 1._JPRB
    ENDIF
  ELSE
    ZOCEA= MAX(1._JPRB-PLSM(JL)-PCLK(JL), 0._JPRB)
  ENDIF
  ! Also exclude sea ice
  ZFROC= ZOCEA*(1._JPRB-PCI(JL))
!-- flux is considered only over full ocean grids, and for their open ocean fraction
  IF (ZOCEA == 1._JPRB) THEN

! code is for three sea salt bins!!!
! assuming bin 1 is the finest bin, 3 the coarsest
    PFLXSS(JL,1)=  ZFROC * (ZMASS1FLUX(IWCL(JL))+ZDMASS1(IWCL(JL))*ZDZSPEED(JL)) / RSS_RH80_MASSFAC
    PFLXSS(JL,2)=  ZFROC * (ZMASS2FLUX(IWCL(JL))+ZDMASS2(IWCL(JL))*ZDZSPEED(JL)) / RSS_RH80_MASSFAC
    PFLXSS(JL,3)=  ZFROC * (ZMASS3FLUX(IWCL(JL))+ZDMASS3(IWCL(JL))*ZDZSPEED(JL)) / RSS_RH80_MASSFAC
!   original IFS SS source line: Did not understand 3.41_JPRB factor, omitted
!   PFLXSS(JL,JBIN) = ZFROC * RSSFLX() * (PWIND(JL)**3.41_JPRB)
  ENDIF
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SSALT_MS',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SSALT_MS
