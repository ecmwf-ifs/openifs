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

MODULE YOEAERC

USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

!     ------------------------------------------------------------------
!*    ** *YOEAERC* - TEGEN/GISS AEROSOL CLIMATOLOGY
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NFDECSO4DIM=185 , NLDECSO4DIM=210
INTEGER(KIND=JPIM), PARAMETER :: ISHIST=1850, IEHIST=2100, INLAT=24


TYPE :: TEAERC_TEGEN
REAL(KIND=JPRB) :: REPAER = -HUGE(1._JPRB)
REAL(KIND=JPRB), POINTER, DIMENSION(:,:,:) :: RAERBC => NULL(), RAEROR => NULL(), RAERSD => NULL()
REAL(KIND=JPRB), POINTER, DIMENSION(:,:,:) :: RAERSS => NULL(), RAERSU => NULL()
REAL(KIND=JPRB), POINTER, DIMENSION(:,:,:) :: RAERV  => NULL()
REAL(KIND=JPRB), POINTER, DIMENSION(:)     :: RLATV  => NULL()
REAL(KIND=JPRB) :: RTAEBC(72,46), RTAEOR(72,46), RTAESD(72,46)
REAL(KIND=JPRB) :: RTAESS(72,46), RTAESU(72,46), RTAEVO(46)
INTEGER(KIND=JPIM) :: ISYEAR=IEHIST, IEYEAR=ISHIST

LOGICAL :: IS_INITIALISED = .FALSE.

CONTAINS 

PROCEDURE :: SETUP => SETUP_TEGEN
FINAL     :: FINALIZE_TEGEN

END TYPE TEAERC_TEGEN

!! storage arrays declared as global, to save on memory / initialisation costs
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:), TARGET, PRIVATE :: RAERBC_STORAGE, RAEROR_STORAGE, RAERSD_STORAGE
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:), TARGET, PRIVATE :: RAERSS_STORAGE, RAERSU_STORAGE
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:), TARGET, PRIVATE :: RAERV_STORAGE
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:),     TARGET, PRIVATE :: RLATV_STORAGE
INTEGER(KIND=JPIM) :: SMART_TEGEN_POINTER_COPIES = 0

!     ------------------------------------------------------------------
!*                   MACC-BASED AEROSOL CLIMATOLOGY
!     ------------------------------------------------------------------

TYPE :: TEAERC_MACC
!- based on optical depth

REAL(KIND=JPRB),ALLOCATABLE, DIMENSION(:,:,:) :: RMACBC, RMACOR, RMACSD, RMACSS, RMACSU !! appear not to be used
REAL(KIND=JPRB),ALLOCATABLE, DIMENSION(:,:)   :: RMAEBC, RMAEOR, RMAESD, RMAESS, RMAESU

!- based on vertically integrated mass

REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: RMACBC1, RMACBC2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:)   :: RMAEBC1, RMAEBC2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: RMACOR1, RMACOR2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:)   :: RMAEOR1, RMAEOR2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: RMACSD1, RMACSD2, RMACSD3
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:)   :: RMAESD1, RMAESD2, RMAESD3
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: RMACSS1, RMACSS2, RMACSS3
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:)   :: RMAESS1, RMAESS2, RMAESS3
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:) :: RMACSU1
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:)   :: RMAESU1

!arrays for input 3D climatology on a 3x3 deg and 60 level grid
!lon,lat,lev,month
REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: RMACBC13D => NULL(), RMACBC23D => NULL()
REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: RMACOR13D => NULL(), RMACOR23D => NULL()
REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: RMACSD13D => NULL(), RMACSD23D => NULL(), RMACSD33D => NULL()
REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: RMACSS13D => NULL(), RMACSS23D => NULL(), RMACSS33D => NULL()
REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: RMACSU13D => NULL()
REAL(KIND=JPRB), POINTER, CONTIGUOUS, DIMENSION(:,:,:,:) :: REF_MON_AER_PRS => NULL(),REF_MON_AER_DPRS => NULL()
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:)   :: RMAEBC13D, RMAEBC23D
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:)   :: RMAEOR13D, RMAEOR23D
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:)   :: RMAESD13D, RMAESD23D, RMAESD33D
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:)   :: RMAESS13D, RMAESS23D, RMAESS33D
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:)   :: RMAESU13D
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:)   :: REF_AER_PRS,REF_AER_DPRS

LOGICAL :: IS_INITIALISED = .FALSE.

CONTAINS 

PROCEDURE :: SETUP => SETUP_MACC
FINAL     :: FINALIZE_MACC

END TYPE TEAERC_MACC

!! storage array declared as global, to save on memory / initialisation costs
INTEGER(KIND=JPIM), PARAMETER, PRIVATE :: NUM_MAC3D_ARRAYS=13
INTEGER(KIND=JPIM), PARAMETER, PRIVATE :: NBC1=1, NBC2=2, NOR1=3, NOR2=4, NSD1=5, NSD2=6, NSD3=7, NSS1=8, NSS2=9, NSS3=10, &
 &                                        NSU1=11, NPRS=12, NDPRS=13
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_BC1
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_BC2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_OR1
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_OR2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_SD1
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_SD2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_SD3
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_SS1
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_SS2
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_SS3
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_SU1
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_PRS
REAL(KIND=JPRB), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET, PRIVATE :: RMAC3D_DPRS
INTEGER(KIND=JPIM) :: SMART_MACC3D_POINTER_COPIES = 0
INTEGER(KIND=JPIM) :: SMART_MACC_POINTER_COPIES = 0

!     ----------------------------------------------------------------------------
!*    ** *YOEAERC* - HISTORICAL TROPOSPHERIC SO4 DATA (TOTAL BURDEN IN G/M^2)
!     ----------------------------------------------------------------------------

!TYPE(REGLATLON_FIELD) ::  RAERSO4  !! moved to YOMRIP, because has time-varying values
CHARACTER(LEN=260),SAVE              :: CLISTSO4,FILESO4(NFDECSO4DIM:NLDECSO4DIM)

LOGICAL              , SAVE ::  LVOLCDATA
CHARACTER (LEN = 260), SAVE ::  CVOLCDATA

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/12/21

!  NAME     TYPE    PURPOSE
!  ----  :  ----  : ---------------------------------------------------
! -- for Tegen et al. 1997 climatology
! RAERBC : REAL   : OPTICAL THICKNESS BLACK CARBON TYPE AEROSOL
! RAEROR : REAL   : OPTICAL THICKNESS ORGANIC TYPE AEROSOL
! RAERSD : REAL   : OPTICAL THICKNESS SOIL-DUST TYPE AEROSOL
! RAERSS : REAL   : OPTICAL THICKNESS SEA-SALT TYPE AEROSOL
! RAERSU : REAL   : OPTICAL THICKNESS SULFATE-TYPE AEROSOL
! RAERV  : REAL   : VOLCANIC OPTICAL THICKNESS 
! RTAEBC : REAL   : TIME-INTERPOLATED BLACK CARBON OPTICAL THICKNESS
! RTAEOR : REAL   : TIME-INTERPOLATED ORGANIC OPTICAL THICKNESS
! RTAESD : REAL   : TIME-INTERPOLATED SOIL-DUST OPTICAL THICKNESS
! RTAESS : REAL   : TIME-INTERPOLATED SEA-SALT OPTICAL THICKNESS
! RTAESU : REAL   : TIME-INTERPOLATED SULFATE CARBON OPTICAL THICKNESS
! RLATV  : REAL   : REFERENCE LATITUDES, stored in type for convenience

! -- for MACC-derived aerosol climatology
! RMACBC : REAL   : OPTICAL THICKNESS BLACK CARBON TYPE AEROSOL
! RMACOR : REAL   : OPTICAL THICKNESS ORGANIC TYPE AEROSOL
! RMACSD : REAL   : OPTICAL THICKNESS SOIL-DUST TYPE AEROSOL
! RMACSS : REAL   : OPTICAL THICKNESS SEA-SALT TYPE AEROSOL
! RMACSU : REAL   : OPTICAL THICKNESS SULFATE-TYPE AEROSOL
! RMAEBC : REAL   : TIME-INTERPOLATED BLACK CARBON OPTICAL THICKNESS
! RMAEOR : REAL   : TIME-INTERPOLATED ORGANIC OPTICAL THICKNESS
! RMAESD : REAL   : TIME-INTERPOLATED SOIL-DUST OPTICAL THICKNESS
! RMAESS : REAL   : TIME-INTERPOLATED SEA-SALT OPTICAL THICKNESS
! RMAESU : REAL   : TIME-INTERPOLATED SULFATE CARBON OPTICAL THICKNESS

! RTAEVO : REAL   : TIME-INTERPOLATED STRATOS.VOLCANIC OPTICAL THICKNESS
! RTAESO4: REAL   : TIME-INTERPOLATED SO4 LOAD (G/M^2)

! REPAER : REAL   : SECURITY PARAMETER FOR AEROSOLS

! LVOLCDATA : LOGICAL   : .TRUE.  READ STRATOS.VOLCANIC OPTICAL THICKNESS FROM FILE
!                         .FALSE. USE HARD-CODED VALUES
! CVOLCDATA : CHARACTER : FILENAME FOR STRATOS.VOLCANIC OPTICAL THICKNESS DATA
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE SETUP_TEGEN(SELF, LCALL_SUECAEC)
  USE YOMLUN, ONLY : NULOUT
  !------------
  CLASS(TEAERC_TEGEN), INTENT(INOUT) :: SELF
  LOGICAL, INTENT(IN)                :: LCALL_SUECAEC
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  !------------
#include "abor1.intfb.h"
#include "suecaebc.intfb.h"
#include "suecaeor.intfb.h"
#include "suecaesd.intfb.h"
#include "suecaess.intfb.h"
#include "suecaesu.intfb.h"
  !------------

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SETUP_TEGEN',0,ZHOOK_HANDLE)

!! check whether the main storage arrays have already been initialised
  IF (SMART_TEGEN_POINTER_COPIES > 0) THEN
    SMART_TEGEN_POINTER_COPIES = SMART_TEGEN_POINTER_COPIES + 1 !! this is used to be able to deallocate cleanly only if number reaches 0
  
    !! point derived_type pointers to main storage arrays
    SELF%RAERBC => RAERBC_STORAGE 
    SELF%RAEROR => RAEROR_STORAGE
    SELF%RAERSD => RAERSD_STORAGE
    SELF%RAERSS => RAERSS_STORAGE
    SELF%RAERSU => RAERSU_STORAGE
    SELF%RAERV  => RAERV_STORAGE
    SELF%RLATV  => RLATV_STORAGE

  ELSE
  
    IF (ANY( (/ ALLOCATED(RAERBC_STORAGE), ALLOCATED(RAEROR_STORAGE), ALLOCATED(RAERSD_STORAGE), &
 &              ALLOCATED(RAERSS_STORAGE), ALLOCATED(RAERSU_STORAGE), ALLOCATED(RAERV_STORAGE ) , ALLOCATED(RLATV_STORAGE) /) ) ) THEN
      CALL ABOR1('Global storage arrays for TEGEN type have already been allocated!')
    ENDIF

    ALLOCATE(RAERBC_STORAGE(72,46,12))
    ALLOCATE(RAEROR_STORAGE(72,46,12))
    ALLOCATE(RAERSD_STORAGE(72,46,12))
    ALLOCATE(RAERSS_STORAGE(72,46,12))
    ALLOCATE(RAERSU_STORAGE(72,46,12))
    ALLOCATE(RAERV_STORAGE(INLAT,ISHIST:IEHIST,12))
    ALLOCATE(RLATV_STORAGE(INLAT))
   
    !! point derived_type pointers to main storage arrays
    SELF%RAERBC => RAERBC_STORAGE 
    SELF%RAEROR => RAEROR_STORAGE
    SELF%RAERSD => RAERSD_STORAGE
    SELF%RAERSS => RAERSS_STORAGE
    SELF%RAERSU => RAERSU_STORAGE
    SELF%RAERV  => RAERV_STORAGE
    SELF%RLATV  => RLATV_STORAGE

    CALL SUECAEBC(RAERBC_STORAGE)
    CALL SUECAEOR(RAEROR_STORAGE)
    CALL SUECAESD(RAERSD_STORAGE)
    CALL SUECAESS(RAERSS_STORAGE)
    CALL SUECAESU(RAERSU_STORAGE)
    IF (LCALL_SUECAEC) CALL SUECAEC(SELF)
    
    !! this is the first time the storage arrays are being set up, so set to 1
    SMART_TEGEN_POINTER_COPIES = 1 !! this is used to be able to deallocate cleanly only if number reaches 0

  ENDIF

  SELF%IS_INITIALISED = .TRUE.

  IF (LHOOK) CALL DR_HOOK('YOEAERC:SETUP_TEGEN',1,ZHOOK_HANDLE)

END SUBROUTINE SETUP_TEGEN

SUBROUTINE FINALIZE_TEGEN(SELF)
  USE YOMLUN, ONLY : NULOUT
  !------------
  TYPE(TEAERC_TEGEN), INTENT(INOUT) :: SELF
  !------------
  INTEGER(KIND=JPIM) :: I_PTR_CNT

#include "abor1.intfb.h" 

!! if self has not been initialized, it should not count in the smart pointer references
IF (SELF%IS_INITIALISED) THEN 

  I_PTR_CNT = SMART_TEGEN_POINTER_COPIES

  !! check whether the main storage arrays are referenced by more than one derived type
  IF (I_PTR_CNT > 1) THEN
    SMART_TEGEN_POINTER_COPIES = SMART_TEGEN_POINTER_COPIES - 1 !! this is used to be able to deallocate cleanly only if number reaches 0
  
  ELSEIF (I_PTR_CNT == 1) THEN
  
    WRITE(NULOUT,*) 'Last reference to tegen storage arrays being destroyed, so deallocating main storage'
  
    DEALLOCATE(RAERBC_STORAGE)
    DEALLOCATE(RAEROR_STORAGE)
    DEALLOCATE(RAERSD_STORAGE)
    DEALLOCATE(RAERSS_STORAGE)
    DEALLOCATE(RAERSU_STORAGE)
    DEALLOCATE(RAERV_STORAGE)
    DEALLOCATE(RLATV_STORAGE)

    SMART_TEGEN_POINTER_COPIES = 0

  ELSEIF (I_PTR_CNT <= 0) THEN
  
    CALL ABOR1('Trying to deallocate tegen storage arrays that have already been deallocated')

  ENDIF

  SELF%RAERBC => NULL()
  SELF%RAEROR => NULL()
  SELF%RAERSD => NULL()
  SELF%RAERSS => NULL()
  SELF%RAERSU => NULL()
  SELF%RAERV  => NULL()
  SELF%RLATV  => NULL()

  SELF%IS_INITIALISED = .FALSE.

ENDIF

END SUBROUTINE FINALIZE_TEGEN

!==========================================================================================

SUBROUTINE SETUP_MACC(SELF,LLAER3D)
  USE YOMLUN, ONLY : NULOUT
  !------------
  CLASS(TEAERC_MACC), INTENT(INOUT) :: SELF
  LOGICAL, INTENT(IN)               :: LLAER3D
  !------------

  CALL SU_AER_CLIMATOLOGY(SELF)
  WRITE(NULOUT,*) '2D CAMS AEROSOL CLIMATOLOGY IS LOADED'
  IF (LLAER3D) THEN
    CALL SU_AER_CLIMATOLOGY3D(SELF)
    WRITE(NULOUT,*) "3D CAMS AEROSOL CLIMATOLOGY IS LOADED"
  ENDIF
  SELF%IS_INITIALISED = .TRUE.

END SUBROUTINE SETUP_MACC

SUBROUTINE FINALIZE_MACC(SELF)
  USE YOMLUN, ONLY : NULOUT
  !------------
  TYPE(TEAERC_MACC), INTENT(INOUT) :: SELF
  !------------
  INTEGER(KIND=JPIM) :: I_PTR_CNT
#include "abor1.intfb.h" 

!! if self has not been initialized, it should not count in the smart pointer references
IF (SELF%IS_INITIALISED) THEN 

  I_PTR_CNT = SMART_MACC3D_POINTER_COPIES

  !! check whether the main storage arrays are referenced by more than one derived type
  IF (I_PTR_CNT > 1) THEN
    SMART_MACC3D_POINTER_COPIES = SMART_MACC3D_POINTER_COPIES - 1 !! this is used to be able to deallocate cleanly only if number reaches 0
  
  ELSEIF (I_PTR_CNT == 1) THEN
  
    WRITE(NULOUT,*) 'Last reference to MACC3D storage arrays being destroyed, so deallocating main storage'
  
   DEALLOCATE(RMAC3D_BC1, RMAC3D_BC2, RMAC3D_OR1, RMAC3D_OR2, RMAC3D_SD1, RMAC3D_SD2, RMAC3D_SD3, &
 &            RMAC3D_SS1, RMAC3D_SS2, RMAC3D_SS3, RMAC3D_SU1, RMAC3D_PRS, RMAC3D_DPRS)

    SMART_MACC3D_POINTER_COPIES = 0

  ELSEIF (I_PTR_CNT <= 0) THEN

    IF (ALLOCATED(RMAC3D_BC1)) THEN !! ideally we would check LAER3D, but it belongs to the parent derived type
      CALL ABOR1('Trying to deallocate MACC3D storage arrays that have already been deallocated')
    ENDIF

  ENDIF

  SELF%RMACBC13D => NULL()
  SELF%RMACBC23D => NULL()
  SELF%RMACOR13D => NULL()
  SELF%RMACOR23D => NULL()
  SELF%RMACSD13D => NULL()
  SELF%RMACSD23D => NULL()
  SELF%RMACSD33D => NULL()
  SELF%RMACSS13D => NULL()
  SELF%RMACSS23D => NULL()
  SELF%RMACSS33D => NULL()
  SELF%RMACSU13D => NULL()
  SELF%REF_MON_AER_PRS => NULL()
  SELF%REF_MON_AER_DPRS => NULL()

  SELF%IS_INITIALISED = .FALSE.

ENDIF

!! storage arrays declared as global, to save on memory / initialisation costs

END SUBROUTINE FINALIZE_MACC


#ifdef RS6K
@PROCESS NOOPTIMIZE
#endif
!pgi$r opt=0 
SUBROUTINE SUECAEC(YDEAERC_TEGEN)
!DEC$ NOOPTIMIZE

!**** *SUECAEC* - DEFINES NEW CLIMATOLOGICAL DISTRIBUTION OF AEROSOLS

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECAEC* 

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===
!     AEROSOL TOTAL OPTICAL THICKNESS (0.55 MICRON)

!-- TEGEN ET AL. MONTHLY CLIMATOLOGY OF TROPOSPHERIC AEROSOLS

! RTAEBC : BLACK CARBON AEROSOL
! RTAEOR : ORGANIC AEROSOL
! RTAESD : SAND DUST AEROSOL
! RTAESS : SEA SALT AEROSOL
! RTAESU : SULFATE AEROSOL

!-- HISTORICAL STRATOSPHERIC VOLCANIC AEROSOLS

! RTAEVO : VOLCANIC AEROSOL

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!      J.-J. MORCRETTE  E.C.M.W.F.    98/12/21

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!      JJ.Morcrette 20080206 revision handling GISS volcanic aerosols
!      H.Hersbach  08-11-2010 Option for reading GISS volcanic aerosols from file
!      F. Vana     05-Mar-2015  Support for single precision
!      JJMorcrette  20120813  new 2D (lat-lon) MACC-based aerosol climatology
!      C. Roberts 2017-02-30, added support for perpetual forcing (NCMIPFIXYR)
!      A Bozzo 2018-01-22 new 3D (lat-lon-lev) CAMS-based aerosol climatology
!      O. Marsden 2018-01-30 Split the update out from the setup, new UPDECAEC routine
!-----------------------------------------------------------------------

USE PARKIND1     , ONLY : JPRD, JPIM, JPRB
USE YOMLUN       , ONLY : NULOUT, FOPEN
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TEAERC_TEGEN), INTENT(INOUT) :: YDEAERC_TEGEN

!     -----------------------------------------------------------------

INTEGER(KIND=JPIM) :: JI, JL, JK
INTEGER(KIND=JPIM) :: ILMO(12)

INTEGER(KIND=JPIM) :: IU, INLEV, IYEAR, ILAT, IDUM
REAL(KIND=JPRB)    :: Z_GRID, ZINT, ZLAT, ZRES, ZRESM, ZRESP, ZTIMI, ZXTIME
REAL(KIND=JPRB)    :: ZW1, ZW2
REAL(KIND=JPRD)    :: ZAERV_D(12), ZLATV_D(INLAT)
CHARACTER(LEN=1000):: CLINE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM),SAVE :: ISYEAR=IEHIST, IEYEAR=ISHIST


#include "abor1.intfb.h"
#include "fcttim.func.h"      

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOEAERC:SUECAEC',0,ZHOOK_HANDLE)

INLEV=0
YDEAERC_TEGEN%RAERV(:,:,:)=0.0_JPRB


IF(TRIM(CVOLCDATA)=="")LVOLCDATA=.FALSE.
IF (LVOLCDATA) THEN
!- Read stratospheric AOD from file in the following GISS ascii format:
!
!
!    Optical Thickness in <xx1>to<yy1> km Level
!    Year/Lat Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec 
!    <year1> <lat1>  <optical thickness, Jan> ... <optical thickness, Dec> 
!            <lat2>  <optical thickness, Jan> ... <optical thickness, Dec>
!             .....
!            <lat24> <optical thickness, Jan> ... <optical thickness, Dec>
!    <year2> <lat1>  <optical thickness, Jan> ... <optical thickness, Dec> 
!             .....
!    Optical Thickness in <xx2>to<yy2> km Level
!             .....

   CALL FOPEN(KUNIT=IU,CDFILE=TRIM(CVOLCDATA),CDFORM='FORMATTED',CDTEXT="SUECAEC VOLCANIC HISTORY FILE")

   ILAT=1
   DO IDUM=1,100*INLAT*(IEHIST-ISHIST+1)  ! Dummy counter up to 100 levels to prevent endless loop
     READ(IU,'(A1000)',END=999)CLINE
     IF (INDEX(CLINE,"Optical")/=0) THEN
        INLEV=INLEV+1
        WRITE(NULOUT,'(A,I2,":",X,A)')"SUECAEC: INLEV=",INLEV,TRIM(CLINE)
        CYCLE
     ENDIF
     IF (INDEX(CLINE,"Year")   /=0) CYCLE

     IF (ILAT==1) THEN
        READ(CLINE,*)IYEAR,ZLATV_D(ILAT),ZAERV_D(1:12)
        ISYEAR=MIN(ISYEAR,IYEAR)
        IEYEAR=MAX(IEYEAR,IYEAR)
     ELSE
        READ(CLINE,*)      ZLATV_D(ILAT),ZAERV_D(1:12)
     ENDIF
     YDEAERC_TEGEN%RLATV(ILAT) = REAL(ZLATV_D(ILAT),JPRB)

!    - Sum all vertical levels to get total optical thickness
!    - In future one might wish to retain level-dependent info
     YDEAERC_TEGEN%RAERV(ILAT,IYEAR,1:12)=YDEAERC_TEGEN%RAERV(ILAT,IYEAR,1:12)+REAL(ZAERV_D(1:12),JPRB)

     ILAT=MOD(ILAT,INLAT)+1
   ENDDO
   CALL ABOR1('SUECAEC: AEROSOL OPTICAL DEPTHS FILE TOO LONG')
999  CONTINUE
   CLOSE (IU)

!  - Write out some statistics
   DO IYEAR=ISYEAR,IEYEAR
     WRITE(NULOUT,*)"SUECAEC: AVERAGE VOLCANIC OPTICAL DEPTH FOR YEAR ",&
     &               IYEAR,":, ",SUM(YDEAERC_TEGEN%RAERV(1:INLAT,IYEAR,1:12))/(INLAT*12)
   ENDDO


!- Use hard-coded Volcanic GISS data; runs from 1957 only; not latest available set
ELSE
  ISYEAR=1957
  IEYEAR =2010
  CALL PART1
  CALL PART2
  CALL PART3
  CALL PART4
  CALL PART5
  CALL PART6
  CALL PART7
!   - Rescale to optical depth
  YDEAERC_TEGEN%RAERV(:,:,:)=1.E-04_JPRB*YDEAERC_TEGEN%RAERV(:,:,:)
ENDIF

WRITE(NULOUT,'(A,I3,2I6)')"SUECAEC: INLEV,ISYEAR,IEYEAR=",INLEV,ISYEAR,IEYEAR
 


!-- security for aerosol optical thickness
YDEAERC_TEGEN%REPAER=1.E-12_JPRB

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOEAERC:SUECAEC',1,ZHOOK_HANDLE)

!-----------------------------------------------------------------------------

CONTAINS


SUBROUTINE PART1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOEAERC:PART1',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!     VOLCANIC AEROSOL OPTICAL THICKNESS (GISS)     

!- reference latitudes

YDEAERC_TEGEN%RLATV(:) = (/&
 & 90.0_JPRB, 82.2_JPRB, 74.3_JPRB, 66.5_JPRB, 58.7_JPRB, 50.9_JPRB, 43.0_JPRB, 35.2_JPRB, 27.4_JPRB, 19.6_JPRB, 11.7_JPRB,&
 & 3.9_JPRB, -3.9_JPRB,-11.7_JPRB,-19.6_JPRB,-27.4_JPRB,-35.2_JPRB,-43.0_JPRB,-50.9_JPRB,-58.7_JPRB,-66.5_JPRB,-74.3_JPRB,&
 & -82.2_JPRB,-90.0_JPRB /)   

!- total optical thickness * 10000

YDEAERC_TEGEN%RAERV(: ,1957, 1) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 2) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 3) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 4) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 5) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 6) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1957,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 1) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 2) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 3) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 4) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 5) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 6) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1958,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 1) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 2) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 3) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 4) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 5) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 6) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1959,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1960, 1) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1960, 2) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1960, 3) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1960, 4) = (/&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17/)  
YDEAERC_TEGEN%RAERV(: ,1960, 5) = (/&
 & 33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,&
 & 33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33/)  
YDEAERC_TEGEN%RAERV(: ,1960, 6) = (/&
 & 50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,&
 & 50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50/)  
YDEAERC_TEGEN%RAERV(: ,1960, 7) = (/&
 & 67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67,&
 & 67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67/)  
YDEAERC_TEGEN%RAERV(: ,1960, 8) = (/&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83/)  
YDEAERC_TEGEN%RAERV(: ,1960, 9) = (/&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100/)  
YDEAERC_TEGEN%RAERV(: ,1960,10) = (/&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100/)  
YDEAERC_TEGEN%RAERV(: ,1960,11) = (/&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100/)  
YDEAERC_TEGEN%RAERV(: ,1960,12) = (/&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,&
 & 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100/)  
YDEAERC_TEGEN%RAERV(: ,1961, 1) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 2) = (/&
 & 1,   1,   1,  40,   1,   1,  20,   1,  72, 190, 190, 101,&
 & 30,  48,  69,  90, 110,  94,  78,  63,  47,  31,  16,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 3) = (/&
 & 1,   1,   1,  50,   1,   1,  10,   1, 145, 380, 380, 230,&
 & 110,  85,  57,  28,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 4) = (/&
 & 60,  60,  60,  50,   1,   1,   1,   1, 175, 460, 460, 272,&
 & 120,  93,  62,  30,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 5) = (/&
 & 80,  80,  80,  40,   1,   1,   1,   1, 110, 290, 290, 129,&
 & 1,  18,  39,  60,  80,  68,  57,  46,  34,  23,  11,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 6) = (/&
 & 90,  90,  90,  10,  50,  50,   1,   1,  61, 160, 160,  71,&
 & 1,  25,  53,  82, 110,  94,  78,  63,  47,  31,  16,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 7) = (/&
 & 70,  70,  70,   1,  70,  70,   1,  10,  52, 120, 120,  70,&
 & 30,  23,  15,   8,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 8) = (/&
 & 40,  40,  40,   1,  60,  60,   1,  80,  95, 120, 120, 153,&
 & 180, 142,  97,  53,  10,   9,   7,   6,   4,   3,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961, 9) = (/&
 & 1,   1,   1,   1,  40,  40,   1,  60,  83, 120, 120, 214,&
 & 290, 236, 173, 111,  50,  43,  36,  28,  21,  14,   7,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961,10) = (/&
 & 1,   1,   1,   1,  20,  20,   1,   1,  38, 100, 100, 228,&
 & 330, 258, 175,  91,  10,   9,   7,   6,   4,   3,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961,11) = (/&
 & 1,   1,   1,  10,  10,  10,  80,   1,   4,  10,  10, 193,&
 & 340, 264, 175,  86,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1961,12) = (/&
 & 1,   1,   1,  10,  10,  10,  80,   1,   8,  20,  20, 164,&
 & 280, 217, 144,  71,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 1) = (/&
 & 1,   1,   1,  10,  20,  20,   1,   1,  53, 140, 140, 112,&
 & 90,  70,  46,  23,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 2) = (/&
 & 1,   1,   1,   1,  50,  50,   1,   1,  61, 160, 160, 116,&
 & 80,  78,  75,  73,  70,  60,  50,  40,  30,  20,  10,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 3) = (/&
 & 1,   1,   1,   1,  70,  70,   1,   1,  11,  30,  30, 152,&
 & 250, 259, 269, 280, 290, 248, 206, 165, 124,  83,  41,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 4) = (/&
 & 1,   1,   1,   1,  60,  60,  50,   1,   1,   1,   1, 216,&
 & 390, 363, 332, 300, 270, 231, 192, 154, 115,  77,  38,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 5) = (/&
 & 1,   1,   1,   1,  40,  40, 100,   1,   1,   1,   1,  28,&
 & 50,  46,  40,  35,  30,  26,  21,  17,  13,   9,   4,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 6) = (/&
 & 1,   1,   1,   1,   1,   1, 150,   1,  19,  50,  50,  33,&
 & 20,  20,  20,  20,  20,  17,  14,  11,   9,   6,   3,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 7) = (/&
 & 1,   1,   1,  20,   1,   1, 160,   1,  46, 120, 120, 270,&
 & 390, 311, 220, 129,  40,  34,  28,  23,  17,  11,   6,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 8) = (/&
 & 1,   1,   1,  20,   1,   1, 140,   1,  57, 150, 150, 405,&
 & 610, 473, 314, 155,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962, 9) = (/&
 & 30,  30,  30,  10,  10,  10,  80,  10,  44, 100, 100, 189,&
 & 260, 202, 134,  66,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962,10) = (/&
 & 80,  80,  80,  10,  20,  20,  10, 120,  78,  10,  10, 110,&
 & 190, 147,  98,  48,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962,11) = (/&
 & 130, 130, 130,  50,  20,  20,   1,  50,  39,  20,  20, 142,&
 & 240, 186, 123,  61,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1962,12) = (/&
 & 180, 180, 180,  60,   1,   1,   1,   1,  19,  50,  50, 122,&
 & 180, 140,  93,  46,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 1) = (/&
 & 230, 230, 230,  10,   1,   1, 100,   1,  53, 140, 140,  62,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 2) = (/&
 & 280, 280, 280,  70,   1,   1, 130,   1,  61, 160, 160,  71,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 3) = (/&
 & 290, 290, 290, 150,   1,   1,  50, 150, 154, 160, 160,  71,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 4) = (/&
 & 230, 230, 230, 150,  20,  20,  70,  90, 113, 150, 150, 294,&
 & 410, 379, 342, 306, 270, 231, 192, 154, 115,  77,  38,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 5) = (/&
 & 80,  80,  80, 100,  40,  40, 130,   1,  53, 140, 140, 578,&
 & 930, 874, 809, 743, 680, 581, 484, 387, 290, 194,  97,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 6) = (/&
 & 1,   1,   1,  20,   1,   1, 120,   1,  76, 200, 200, 621,&
 & 960,1025,1101,1176,1250,1067, 889, 711, 534, 356, 178,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 7) = (/&
 & 1,   1,   1,   1,   1,   1,  50, 140, 208, 320, 320, 608,&
 & 840, 975,1131,1288,1440,1229,1025, 820, 615, 410, 205,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 8) = (/&
 & 1,   1,   1,   1,   1,   1,  50, 280, 333, 420, 420, 608,&
 & 760, 946,1163,1379,1590,1357,1131, 905, 679, 452, 226,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963, 9) = (/&
 & 30,  30,  30,   1,  20,  20,  60, 200, 299, 460, 460, 698,&
 & 890,1072,1283,1495,1700,1451,1209, 968, 726, 484, 242,   1/)  
YDEAERC_TEGEN%RAERV(: ,1963,10) = (/&
 & 90,  90,  90,  10,  90,  90, 120, 180, 271, 420, 420,1135,&
 & 1710,1714,1720,1725,1730,1524,1323,1123, 922, 721, 521, 320/)  
YDEAERC_TEGEN%RAERV(: ,1963,11) = (/&
 & 160, 160, 160, 110, 150, 150, 210, 420, 352, 240, 240,1360,&
 & 2260,2116,1949,1782,1620,1560,1502,1443,1385,1327,1268,1210/)  
YDEAERC_TEGEN%RAERV(: ,1963,12) = (/&
 & 220, 220, 220, 180, 210, 210, 330, 490, 364, 160, 160,1180,&
 & 2000,1823,1616,1410,1210,1264,1317,1369,1422,1475,1527,1580/)  
YDEAERC_TEGEN%RAERV(: ,1964, 1) = (/&
 & 280, 280, 280, 230, 320, 320, 370, 460, 346, 160, 160, 753,&
 & 1230,1154,1065, 976, 890, 985,1078,1170,1263,1355,1448,1540/)  
YDEAERC_TEGEN%RAERV(: ,1964, 2) = (/&
 & 340, 340, 340, 300, 390, 390, 350, 370, 317, 230, 230, 707,&
 & 1090,1029, 959, 888, 820, 911, 999,1087,1175,1264,1352,1440/)  
YDEAERC_TEGEN%RAERV(: ,1964, 3) = (/&
 & 390, 390, 390, 360, 380, 380, 290, 240, 270, 320, 320, 808,&
 & 1200,1106, 996, 887, 780, 859, 936,1013,1089,1166,1243,1320/)  
YDEAERC_TEGEN%RAERV(: ,1964, 4) = (/&
 & 430, 430, 430, 390, 310, 310, 190, 150, 245, 400, 400,1049,&
 & 1570,1379,1157, 936, 720, 792, 861, 931,1001,1071,1140,1210/)  
YDEAERC_TEGEN%RAERV(: ,1964, 5) = (/&
 & 450, 450, 450, 270,   1,   1, 100,  90, 223, 440, 440,1188,&
 & 1790,1534,1237, 939, 650, 716, 780, 844, 908, 972,1036,1100/)  

IF (LHOOK) CALL DR_HOOK('YOEAERC:PART1',1,ZHOOK_HANDLE)
END SUBROUTINE PART1


SUBROUTINE PART2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOEAERC:PART2',0,ZHOOK_HANDLE)

YDEAERC_TEGEN%RAERV(: ,1964, 6) = (/&
 & 430, 430, 430, 220,   1,   1,  30,  90, 231, 460, 460,1397,&
 & 2150,1827,1451,1075, 710, 752, 794, 835, 876, 917, 959,1000/)  
YDEAERC_TEGEN%RAERV(: ,1964, 7) = (/&
 & 400, 400, 400, 250,  30,  30,  10, 170, 277, 450, 450,1697,&
 & 2700,2296,1826,1357, 900, 899, 897, 896, 894, 893, 891, 890/)  
YDEAERC_TEGEN%RAERV(: ,1964, 8) = (/&
 & 320, 320, 320, 290,  90,  90,  10, 240, 308, 420, 420,1734,&
 & 2790,2363,1868,1372, 890, 875, 861, 847, 833, 818, 804, 790/)  
YDEAERC_TEGEN%RAERV(: ,1964, 9) = (/&
 & 220, 220, 220, 360, 160, 160,  20, 320, 343, 380, 380,1455,&
 & 2320,1938,1495,1051, 620, 629, 637, 646, 654, 663, 671, 680/)  
YDEAERC_TEGEN%RAERV(: ,1964,10) = (/&
 & 120, 120, 120, 330, 230, 230, 100, 430, 384, 310, 310,1053,&
 & 1650,1356,1014, 672, 340, 396, 450, 504, 558, 612, 666, 720/)  
YDEAERC_TEGEN%RAERV(: ,1964,11) = (/&
 & 80,  80,  80, 260, 320, 320, 230, 520, 417, 250, 250, 727,&
 & 1110, 903, 663, 423, 190, 292, 392, 492, 591, 691, 790, 890/)  
YDEAERC_TEGEN%RAERV(: ,1964,12) = (/&
 & 130, 130, 130, 250, 410, 410, 430, 530, 397, 180, 180, 457,&
 & 680, 561, 423, 284, 150, 273, 392, 512, 631, 751, 870, 990/)  
YDEAERC_TEGEN%RAERV(: ,1965, 1) = (/&
 & 320, 320, 320, 290, 490, 490, 490, 490, 361, 150, 150, 167,&
 & 180, 178, 175, 173, 170, 288, 404, 519, 634, 749, 865, 980/)  
YDEAERC_TEGEN%RAERV(: ,1965, 2) = (/&
 & 570, 570, 570, 430, 500, 500, 450, 420, 317, 150, 150,  67,&
 & 1,  49, 107, 164, 220, 314, 405, 496, 587, 678, 769, 860/)  
YDEAERC_TEGEN%RAERV(: ,1965, 3) = (/&
 & 660, 660, 660, 640, 460, 460, 290, 390, 299, 150, 150,  67,&
 & 1,  56, 121, 187, 250, 316, 380, 444, 508, 572, 636, 700/)  
YDEAERC_TEGEN%RAERV(: ,1965, 4) = (/&
 & 600, 600, 600, 660, 400, 400, 200, 330, 296, 240, 240, 107,&
 & 1,  72, 155, 239, 320, 351, 381, 410, 440, 470, 500, 530/)  
YDEAERC_TEGEN%RAERV(: ,1965, 5) = (/&
 & 520, 520, 520, 630, 270, 270, 170, 190, 239, 320, 320, 143,&
 & 1,  90, 194, 299, 400, 400, 400, 400, 400, 400, 400, 400/)  
YDEAERC_TEGEN%RAERV(: ,1965, 6) = (/&
 & 470, 470, 470, 550, 250, 250,  60,  20, 146, 350, 350, 156,&
 & 1,  90, 194, 299, 400, 380, 360, 340, 320, 300, 280, 260/)  
YDEAERC_TEGEN%RAERV(: ,1965, 7) = (/&
 & 400, 400, 400, 490, 260, 260,   1, 170, 242, 360, 360, 388,&
 & 410, 379, 342, 306, 270, 255, 241, 227, 213, 198, 184, 170/)  
YDEAERC_TEGEN%RAERV(: ,1965, 8) = (/&
 & 300, 300, 300, 420, 280, 280,   1, 210, 240, 290, 290, 440,&
 & 560, 459, 342, 224, 110, 113, 116, 119, 121, 124, 127, 130/)  
YDEAERC_TEGEN%RAERV(: ,1965, 9) = (/&
 & 200, 200, 200, 340, 290, 290,  30, 240, 248, 260, 260, 399,&
 & 510, 411, 296, 182,  70,  80,  90, 100, 110, 120, 130, 140/)  
YDEAERC_TEGEN%RAERV(: ,1965,10) = (/&
 & 100, 100, 100, 250, 290, 290, 220, 270, 262, 250, 250, 355,&
 & 440, 368, 285, 201, 120, 130, 140, 150, 160, 170, 180, 190/)  
YDEAERC_TEGEN%RAERV(: ,1965,11) = (/&
 & 1,   1,   1, 180, 280, 280, 300, 290, 275, 250, 250, 305,&
 & 350, 339, 326, 313, 300, 287, 274, 261, 248, 236, 223, 210/)  
YDEAERC_TEGEN%RAERV(: ,1965,12) = (/&
 & 1,   1,   1, 110, 270, 270, 350, 300, 296, 290, 290, 296,&
 & 300, 338, 383, 427, 470, 417, 366, 315, 264, 212, 161, 110/)  
YDEAERC_TEGEN%RAERV(: ,1966, 1) = (/&
 & 335, 335, 335, 335, 335, 335, 335, 335, 335, 335, 335, 335,&
 & 335, 335, 335, 335, 335, 335, 335, 335, 335, 335, 335, 335/)  
YDEAERC_TEGEN%RAERV(: ,1966, 2) = (/&
 & 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320,&
 & 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320/)  
YDEAERC_TEGEN%RAERV(: ,1966, 3) = (/&
 & 305, 305, 305, 305, 305, 305, 305, 305, 305, 305, 305, 305,&
 & 305, 305, 305, 305, 305, 305, 305, 305, 305, 305, 305, 305/)  
YDEAERC_TEGEN%RAERV(: ,1966, 4) = (/&
 & 290, 290, 290, 290, 290, 290, 290, 290, 290, 290, 290, 290,&
 & 290, 290, 290, 290, 290, 290, 290, 290, 290, 290, 290, 290/)  
YDEAERC_TEGEN%RAERV(: ,1966, 5) = (/&
 & 275, 275, 275, 275, 275, 275, 275, 275, 275, 275, 275, 275,&
 & 275, 275, 275, 275, 275, 275, 275, 275, 275, 275, 275, 275/)  
YDEAERC_TEGEN%RAERV(: ,1966, 6) = (/&
 & 260, 260, 260, 260, 260, 260, 260, 260, 260, 260, 260, 260,&
 & 260, 260, 260, 260, 260, 260, 260, 260, 260, 260, 260, 260/)  
YDEAERC_TEGEN%RAERV(: ,1966, 7) = (/&
 & 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245,&
 & 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245/)  
YDEAERC_TEGEN%RAERV(: ,1966, 8) = (/&
 & 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230,&
 & 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230/)  
YDEAERC_TEGEN%RAERV(: ,1966, 9) = (/&
 & 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215,&
 & 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215/)  
YDEAERC_TEGEN%RAERV(: ,1966,10) = (/&
 & 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,&
 & 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200/)  
YDEAERC_TEGEN%RAERV(: ,1966,11) = (/&
 & 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185,&
 & 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185/)  
YDEAERC_TEGEN%RAERV(: ,1966,12) = (/&
 & 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170,&
 & 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170/)  
YDEAERC_TEGEN%RAERV(: ,1967, 1) = (/&
 & 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155,&
 & 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155/)  
YDEAERC_TEGEN%RAERV(: ,1967, 2) = (/&
 & 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140,&
 & 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140/)  
YDEAERC_TEGEN%RAERV(: ,1967, 3) = (/&
 & 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125,&
 & 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125/)  
YDEAERC_TEGEN%RAERV(: ,1967, 4) = (/&
 & 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,&
 & 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110/)  
YDEAERC_TEGEN%RAERV(: ,1967, 5) = (/&
 & 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130,&
 & 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130/)  
YDEAERC_TEGEN%RAERV(: ,1967, 6) = (/&
 & 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150,&
 & 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150/)  
YDEAERC_TEGEN%RAERV(: ,1967, 7) = (/&
 & 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170,&
 & 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170/)  
YDEAERC_TEGEN%RAERV(: ,1967, 8) = (/&
 & 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190,&
 & 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190/)  
YDEAERC_TEGEN%RAERV(: ,1967, 9) = (/&
 & 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210,&
 & 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210/)  
YDEAERC_TEGEN%RAERV(: ,1967,10) = (/&
 & 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230,&
 & 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230/)  
YDEAERC_TEGEN%RAERV(: ,1967,11) = (/&
 & 207, 207, 207, 207, 207, 207, 207, 207, 207, 207, 207, 207,&
 & 207, 207, 207, 207, 207, 207, 207, 207, 207, 207, 207, 207/)  
YDEAERC_TEGEN%RAERV(: ,1967,12) = (/&
 & 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183,&
 & 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183/)  
YDEAERC_TEGEN%RAERV(: ,1968, 1) = (/&
 & 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160,&
 & 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160/)  
YDEAERC_TEGEN%RAERV(: ,1968, 2) = (/&
 & 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137,&
 & 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137/)  
YDEAERC_TEGEN%RAERV(: ,1968, 3) = (/&
 & 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113,&
 & 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113/)  
YDEAERC_TEGEN%RAERV(: ,1968, 4) = (/&
 & 90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,&
 & 90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90/)  
YDEAERC_TEGEN%RAERV(: ,1968, 5) = (/&
 & 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170,&
 & 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170/)  
YDEAERC_TEGEN%RAERV(: ,1968, 6) = (/&
 & 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250,&
 & 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250/)  
YDEAERC_TEGEN%RAERV(: ,1968, 7) = (/&
 & 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330,&
 & 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330, 330/)  
YDEAERC_TEGEN%RAERV(: ,1968, 8) = (/&
 & 410, 410, 410, 410, 410, 410, 410, 410, 410, 410, 410, 410,&
 & 410, 410, 410, 410, 410, 410, 410, 410, 410, 410, 410, 410/)  
YDEAERC_TEGEN%RAERV(: ,1968, 9) = (/&
 & 490, 490, 490, 490, 490, 490, 490, 490, 490, 490, 490, 490,&
 & 490, 490, 490, 490, 490, 490, 490, 490, 490, 490, 490, 490/)  
YDEAERC_TEGEN%RAERV(: ,1968,10) = (/&
 & 570, 570, 570, 570, 570, 570, 570, 570, 570, 570, 570, 570,&
 & 570, 570, 570, 570, 570, 570, 570, 570, 570, 570, 570, 570/)  
YDEAERC_TEGEN%RAERV(: ,1968,11) = (/&
 & 551, 551, 551, 551, 551, 551, 551, 551, 551, 551, 551, 551,&
 & 551, 551, 551, 551, 551, 551, 551, 551, 551, 551, 551, 551/)  
YDEAERC_TEGEN%RAERV(: ,1968,12) = (/&
 & 533, 533, 533, 533, 533, 533, 533, 533, 533, 533, 533, 533,&
 & 533, 533, 533, 533, 533, 533, 533, 533, 533, 533, 533, 533/)  
YDEAERC_TEGEN%RAERV(: ,1969, 1) = (/&
 & 514, 514, 514, 514, 514, 514, 514, 514, 514, 514, 514, 514,&
 & 514, 514, 514, 514, 514, 514, 514, 514, 514, 514, 514, 514/)  
YDEAERC_TEGEN%RAERV(: ,1969, 2) = (/&
 & 496, 496, 496, 496, 496, 496, 496, 496, 496, 496, 496, 496,&
 & 496, 496, 496, 496, 496, 496, 496, 496, 496, 496, 496, 496/)  
YDEAERC_TEGEN%RAERV(: ,1969, 3) = (/&
 & 477, 477, 477, 477, 477, 477, 477, 477, 477, 477, 477, 477,&
 & 477, 477, 477, 477, 477, 477, 477, 477, 477, 477, 477, 477/)  
YDEAERC_TEGEN%RAERV(: ,1969, 4) = (/&
 & 459, 459, 459, 459, 459, 459, 459, 459, 459, 459, 459, 459,&
 & 459, 459, 459, 459, 459, 459, 459, 459, 459, 459, 459, 459/)  
YDEAERC_TEGEN%RAERV(: ,1969, 5) = (/&
 & 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440,&
 & 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440, 440/)  
YDEAERC_TEGEN%RAERV(: ,1969, 6) = (/&
 & 421, 421, 421, 421, 421, 421, 421, 421, 421, 421, 421, 421,&
 & 421, 421, 421, 421, 421, 421, 421, 421, 421, 421, 421, 421/)  
YDEAERC_TEGEN%RAERV(: ,1969, 7) = (/&
 & 403, 403, 403, 403, 403, 403, 403, 403, 403, 403, 403, 403,&
 & 403, 403, 403, 403, 403, 403, 403, 403, 403, 403, 403, 403/)  
YDEAERC_TEGEN%RAERV(: ,1969, 8) = (/&
 & 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384,&
 & 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384, 384/)  
YDEAERC_TEGEN%RAERV(: ,1969, 9) = (/&
 & 366, 366, 366, 366, 366, 366, 366, 366, 366, 366, 366, 366,&
 & 366, 366, 366, 366, 366, 366, 366, 366, 366, 366, 366, 366/)  
YDEAERC_TEGEN%RAERV(: ,1969,10) = (/&
 & 347, 347, 347, 347, 347, 347, 347, 347, 347, 347, 347, 347,&
 & 347, 347, 347, 347, 347, 347, 347, 347, 347, 347, 347, 347/)  
YDEAERC_TEGEN%RAERV(: ,1969,11) = (/&
 & 329, 329, 329, 329, 329, 329, 329, 329, 329, 329, 329, 329,&
 & 329, 329, 329, 329, 329, 329, 329, 329, 329, 329, 329, 329/)  
YDEAERC_TEGEN%RAERV(: ,1969,12) = (/&
 & 310, 310, 310, 310, 310, 310, 310, 310, 310, 310, 310, 310,&
 & 310, 310, 310, 310, 310, 310, 310, 310, 310, 310, 310, 310/)  
YDEAERC_TEGEN%RAERV(: ,1970, 1) = (/&
 & 291, 291, 291, 291, 291, 291, 291, 291, 291, 291, 291, 291,&
 & 291, 291, 291, 291, 291, 291, 291, 291, 291, 291, 291, 291/)  
YDEAERC_TEGEN%RAERV(: ,1970, 2) = (/&
 & 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273,&
 & 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273/)  
YDEAERC_TEGEN%RAERV(: ,1970, 3) = (/&
 & 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254,&
 & 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254/)  
YDEAERC_TEGEN%RAERV(: ,1970, 4) = (/&
 & 236, 236, 236, 236, 236, 236, 236, 236, 236, 236, 236, 236,&
 & 236, 236, 236, 236, 236, 236, 236, 236, 236, 236, 236, 236/)  
YDEAERC_TEGEN%RAERV(: ,1970, 5) = (/&
 & 217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 217,&
 & 217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 217, 217/)  
YDEAERC_TEGEN%RAERV(: ,1970, 6) = (/&
 & 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199,&
 & 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199/)  
YDEAERC_TEGEN%RAERV(: ,1970, 7) = (/&
 & 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180,&
 & 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180/)  
YDEAERC_TEGEN%RAERV(: ,1970, 8) = (/&
 & 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161,&
 & 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161/)  
YDEAERC_TEGEN%RAERV(: ,1970, 9) = (/&
 & 143, 143, 143, 143, 143, 143, 143, 143, 143, 143, 143, 143,&
 & 143, 143, 143, 143, 143, 143, 143, 143, 143, 143, 143, 143/)  
YDEAERC_TEGEN%RAERV(: ,1970,10) = (/&
 & 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124,&
 & 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124/)  
YDEAERC_TEGEN%RAERV(: ,1970,11) = (/&
 & 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106,&
 & 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106/)  
YDEAERC_TEGEN%RAERV(: ,1970,12) = (/&
 & 87,  87,  87,  87,  87,  87,  87,  87,  87,  87,  87,  87,&
 & 87,  87,  87,  87,  87,  87,  87,  87,  87,  87,  87,  87/)  
YDEAERC_TEGEN%RAERV(: ,1971, 1) = (/&
 & 69,  69,  69,  69,  69,  69,  69,  69,  69,  69,  69,  69,&
 & 69,  69,  69,  69,  69,  69,  69,  69,  69,  69,  69,  69/)  
YDEAERC_TEGEN%RAERV(: ,1971, 2) = (/&
 & 50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,&
 & 50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50/)  
YDEAERC_TEGEN%RAERV(: ,1971, 3) = (/&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42/)  
YDEAERC_TEGEN%RAERV(: ,1971, 4) = (/&
 & 33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,&
 & 33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33/)  
YDEAERC_TEGEN%RAERV(: ,1971, 5) = (/&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25/)  
YDEAERC_TEGEN%RAERV(: ,1971, 6) = (/&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17/)  
YDEAERC_TEGEN%RAERV(: ,1971, 7) = (/&
 & 8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,&
 & 8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8/)  
YDEAERC_TEGEN%RAERV(: ,1971, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,1971, 9) = (/&
 & 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,&
 & 2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2,   2/)  
YDEAERC_TEGEN%RAERV(: ,1971,10) = (/&
 & 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,&
 & 4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4/)  
YDEAERC_TEGEN%RAERV(: ,1971,11) = (/&
 & 6,   6,   6,   6,   6,   6,   6,   6,   6,   6,   6,   6,&
 & 6,   6,   6,   6,   6,   6,   6,   6,   6,   6,   6,   6/)  
YDEAERC_TEGEN%RAERV(: ,1971,12) = (/&
 & 8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,&
 & 8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8/)  
YDEAERC_TEGEN%RAERV(: ,1972, 1) = (/&
 & 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,&
 & 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10/)  
YDEAERC_TEGEN%RAERV(: ,1972, 2) = (/&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14/)  
YDEAERC_TEGEN%RAERV(: ,1972, 3) = (/&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18/)  
YDEAERC_TEGEN%RAERV(: ,1972, 4) = (/&
 & 21,  21,  21,  21,  21,  21,  21,  21,  21,  21,  21,  21,&
 & 21,  21,  21,  21,  21,  21,  21,  21,  21,  21,  21,  21/)  
YDEAERC_TEGEN%RAERV(: ,1972, 5) = (/&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25/)  
YDEAERC_TEGEN%RAERV(: ,1972, 6) = (/&
 & 29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29,&
 & 29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29/)  

IF (LHOOK) CALL DR_HOOK('YOEAERC:PART2',1,ZHOOK_HANDLE)
END SUBROUTINE PART2


SUBROUTINE PART3

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOEAERC:PART3',0,ZHOOK_HANDLE)

YDEAERC_TEGEN%RAERV(: ,1972, 7) = (/&
 & 33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,&
 & 33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33/)  
YDEAERC_TEGEN%RAERV(: ,1972, 8) = (/&
 & 37,  37,  37,  37,  37,  37,  37,  37,  37,  37,  37,  37,&
 & 37,  37,  37,  37,  37,  37,  37,  37,  37,  37,  37,  37/)  
YDEAERC_TEGEN%RAERV(: ,1972, 9) = (/&
 & 41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41,&
 & 41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41/)  
YDEAERC_TEGEN%RAERV(: ,1972,10) = (/&
 & 44,  44,  44,  44,  44,  44,  44,  44,  44,  44,  44,  44,&
 & 44,  44,  44,  44,  44,  44,  44,  44,  44,  44,  44,  44/)  
YDEAERC_TEGEN%RAERV(: ,1972,11) = (/&
 & 48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48,&
 & 48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48/)  
YDEAERC_TEGEN%RAERV(: ,1972,12) = (/&
 & 52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52,&
 & 52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52/)  
YDEAERC_TEGEN%RAERV(: ,1973, 1) = (/&
 & 56,  56,  56,  56,  56,  56,  56,  56,  56,  56,  56,  56,&
 & 56,  56,  56,  56,  56,  56,  56,  56,  56,  56,  56,  56/)  
YDEAERC_TEGEN%RAERV(: ,1973, 2) = (/&
 & 60,  60,  60,  60,  60,  60,  60,  60,  60,  60,  60,  60,&
 & 60,  60,  60,  60,  60,  60,  60,  60,  60,  60,  60,  60/)  
YDEAERC_TEGEN%RAERV(: ,1973, 3) = (/&
 & 64,  64,  64,  64,  64,  64,  64,  64,  64,  64,  64,  64,&
 & 64,  64,  64,  64,  64,  64,  64,  64,  64,  64,  64,  64/)  
YDEAERC_TEGEN%RAERV(: ,1973, 4) = (/&
 & 67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67,&
 & 67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67,  67/)  
YDEAERC_TEGEN%RAERV(: ,1973, 5) = (/&
 & 71,  71,  71,  71,  71,  71,  71,  71,  71,  71,  71,  71,&
 & 71,  71,  71,  71,  71,  71,  71,  71,  71,  71,  71,  71/)  
YDEAERC_TEGEN%RAERV(: ,1973, 6) = (/&
 & 75,  75,  75,  75,  75,  75,  75,  75,  75,  75,  75,  75,&
 & 75,  75,  75,  75,  75,  75,  75,  75,  75,  75,  75,  75/)  
YDEAERC_TEGEN%RAERV(: ,1973, 7) = (/&
 & 79,  79,  79,  79,  79,  79,  79,  79,  79,  79,  79,  79,&
 & 79,  79,  79,  79,  79,  79,  79,  79,  79,  79,  79,  79/)  
YDEAERC_TEGEN%RAERV(: ,1973, 8) = (/&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83/)  
YDEAERC_TEGEN%RAERV(: ,1973, 9) = (/&
 & 86,  86,  86,  86,  86,  86,  86,  86,  86,  86,  86,  86,&
 & 86,  86,  86,  86,  86,  86,  86,  86,  86,  86,  86,  86/)  
YDEAERC_TEGEN%RAERV(: ,1973,10) = (/&
 & 90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,&
 & 90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90,  90/)  
YDEAERC_TEGEN%RAERV(: ,1973,11) = (/&
 & 94,  94,  94,  94,  94,  94,  94,  94,  94,  94,  94,  94,&
 & 94,  94,  94,  94,  94,  94,  94,  94,  94,  94,  94,  94/)  
YDEAERC_TEGEN%RAERV(: ,1973,12) = (/&
 & 98,  98,  98,  98,  98,  98,  98,  98,  98,  98,  98,  98,&
 & 98,  98,  98,  98,  98,  98,  98,  98,  98,  98,  98,  98/)  
YDEAERC_TEGEN%RAERV(: ,1974, 1) = (/&
 & 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102,&
 & 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102/)  
YDEAERC_TEGEN%RAERV(: ,1974, 2) = (/&
 & 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106,&
 & 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106/)  
YDEAERC_TEGEN%RAERV(: ,1974, 3) = (/&
 & 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109,&
 & 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109/)  
YDEAERC_TEGEN%RAERV(: ,1974, 4) = (/&
 & 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113,&
 & 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113/)  
YDEAERC_TEGEN%RAERV(: ,1974, 5) = (/&
 & 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117,&
 & 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117/)  
YDEAERC_TEGEN%RAERV(: ,1974, 6) = (/&
 & 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121,&
 & 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121/)  
YDEAERC_TEGEN%RAERV(: ,1974, 7) = (/&
 & 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125,&
 & 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125/)  
YDEAERC_TEGEN%RAERV(: ,1974, 8) = (/&
 & 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129,&
 & 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129/)  
YDEAERC_TEGEN%RAERV(: ,1974, 9) = (/&
 & 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132,&
 & 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132/)  
YDEAERC_TEGEN%RAERV(: ,1974,10) = (/&
 & 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136,&
 & 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136/)  
YDEAERC_TEGEN%RAERV(: ,1974,11) = (/&
 & 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140,&
 & 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140/)  
YDEAERC_TEGEN%RAERV(: ,1974,12) = (/&
 & 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183,&
 & 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183/)  
YDEAERC_TEGEN%RAERV(: ,1975, 1) = (/&
 & 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227,&
 & 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227/)  
YDEAERC_TEGEN%RAERV(: ,1975, 2) = (/&
 & 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270,&
 & 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270/)  
YDEAERC_TEGEN%RAERV(: ,1975, 3) = (/&
 & 313, 313, 313, 313, 313, 313, 313, 313, 313, 313, 313, 313,&
 & 313, 313, 313, 313, 313, 313, 313, 313, 313, 313, 313, 313/)  
YDEAERC_TEGEN%RAERV(: ,1975, 4) = (/&
 & 357, 357, 357, 357, 357, 357, 357, 357, 357, 357, 357, 357,&
 & 357, 357, 357, 357, 357, 357, 357, 357, 357, 357, 357, 357/)  
YDEAERC_TEGEN%RAERV(: ,1975, 5) = (/&
 & 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400,&
 & 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400/)  
YDEAERC_TEGEN%RAERV(: ,1975, 6) = (/&
 & 337, 337, 337, 337, 337, 337, 337, 337, 337, 337, 337, 337,&
 & 337, 337, 337, 337, 337, 337, 337, 337, 337, 337, 337, 337/)  
YDEAERC_TEGEN%RAERV(: ,1975, 7) = (/&
 & 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273,&
 & 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273, 273/)  
YDEAERC_TEGEN%RAERV(: ,1975, 8) = (/&
 & 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210,&
 & 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210, 210/)  
YDEAERC_TEGEN%RAERV(: ,1975, 9) = (/&
 & 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147,&
 & 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147/)  
YDEAERC_TEGEN%RAERV(: ,1975,10) = (/&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83/)  
YDEAERC_TEGEN%RAERV(: ,1975,11) = (/&
 & 20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,&
 & 20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20/)  
YDEAERC_TEGEN%RAERV(: ,1975,12) = (/&
 & 20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,&
 & 20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20,  20/)  
YDEAERC_TEGEN%RAERV(: ,1976, 1) = (/&
 & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,&
 & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19/)  
YDEAERC_TEGEN%RAERV(: ,1976, 2) = (/&
 & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,&
 & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19/)  
YDEAERC_TEGEN%RAERV(: ,1976, 3) = (/&
 & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,&
 & 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19,  19/)  
YDEAERC_TEGEN%RAERV(: ,1976, 4) = (/&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18/)  
YDEAERC_TEGEN%RAERV(: ,1976, 5) = (/&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18/)  
YDEAERC_TEGEN%RAERV(: ,1976, 6) = (/&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,&
 & 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18,  18/)  
YDEAERC_TEGEN%RAERV(: ,1976, 7) = (/&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17/)  
YDEAERC_TEGEN%RAERV(: ,1976, 8) = (/&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,&
 & 17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17,  17/)  
YDEAERC_TEGEN%RAERV(: ,1976, 9) = (/&
 & 16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,&
 & 16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16/)  
YDEAERC_TEGEN%RAERV(: ,1976,10) = (/&
 & 16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,&
 & 16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16/)  
YDEAERC_TEGEN%RAERV(: ,1976,11) = (/&
 & 16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,&
 & 16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16,  16/)  
YDEAERC_TEGEN%RAERV(: ,1976,12) = (/&
 & 15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,&
 & 15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15/)  
YDEAERC_TEGEN%RAERV(: ,1977, 1) = (/&
 & 15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,&
 & 15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15/)  
YDEAERC_TEGEN%RAERV(: ,1977, 2) = (/&
 & 15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,&
 & 15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15/)  
YDEAERC_TEGEN%RAERV(: ,1977, 3) = (/&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14/)  
YDEAERC_TEGEN%RAERV(: ,1977, 4) = (/&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14/)  
YDEAERC_TEGEN%RAERV(: ,1977, 5) = (/&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,&
 & 14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14/)  
YDEAERC_TEGEN%RAERV(: ,1977, 6) = (/&
 & 13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,&
 & 13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13/)  
YDEAERC_TEGEN%RAERV(: ,1977, 7) = (/&
 & 13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,&
 & 13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13/)  
YDEAERC_TEGEN%RAERV(: ,1977, 8) = (/&
 & 13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,&
 & 13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13,  13/)  
YDEAERC_TEGEN%RAERV(: ,1977, 9) = (/&
 & 12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,&
 & 12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12/)  
YDEAERC_TEGEN%RAERV(: ,1977,10) = (/&
 & 12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,&
 & 12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12,  12/)  
YDEAERC_TEGEN%RAERV(: ,1977,11) = (/&
 & 11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,&
 & 11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11/)  
YDEAERC_TEGEN%RAERV(: ,1977,12) = (/&
 & 11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,&
 & 11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11/)  
YDEAERC_TEGEN%RAERV(: ,1978, 1) = (/&
 & 11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,&
 & 11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11,  11/)  
YDEAERC_TEGEN%RAERV(: ,1978, 2) = (/&
 & 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,&
 & 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10/)  
YDEAERC_TEGEN%RAERV(: ,1978, 3) = (/&
 & 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,&
 & 10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10/)  
YDEAERC_TEGEN%RAERV(: ,1978, 4) = (/&
 & 28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,&
 & 28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28/)  
YDEAERC_TEGEN%RAERV(: ,1978, 5) = (/&
 & 47,  47,  47,  47,  47,  47,  47,  47,  47,  47,  47,  47,&
 & 47,  47,  47,  47,  47,  47,  47,  47,  47,  47,  47,  47/)  
YDEAERC_TEGEN%RAERV(: ,1978, 6) = (/&
 & 65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65,&
 & 65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65/)  
YDEAERC_TEGEN%RAERV(: ,1978, 7) = (/&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,&
 & 83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83,  83/)  
YDEAERC_TEGEN%RAERV(: ,1978, 8) = (/&
 & 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102,&
 & 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102/)  
YDEAERC_TEGEN%RAERV(: ,1978, 9) = (/&
 & 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120,&
 & 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120/)  
YDEAERC_TEGEN%RAERV(: ,1978,10) = (/&
 & 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122,&
 & 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122/)  
YDEAERC_TEGEN%RAERV(: ,1978,11) = (/&
 & 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123,&
 & 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123/)  
YDEAERC_TEGEN%RAERV(: ,1978,12) = (/&
 & 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125,&
 & 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125/)  
YDEAERC_TEGEN%RAERV(: ,1979, 1) = (/&
 & 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,&
 & 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127/)  
YDEAERC_TEGEN%RAERV(: ,1979, 2) = (/&
 & 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,&
 & 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128/)  
YDEAERC_TEGEN%RAERV(: ,1979, 3) = (/&
 & 27,  27,  27,  27,  27,  27,  27,  27,  27,  27,  27,  27,&
 & 22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22/)  
YDEAERC_TEGEN%RAERV(: ,1979, 4) = (/&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,&
 & 23,  23,  23,  23,  23,  23,  23,  23,  23,  23,  23,  23/)  
YDEAERC_TEGEN%RAERV(: ,1979, 5) = (/&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,&
 & 24,  24,  24,  24,  24,  24,  24,  24,  24,  24,  24,  24/)  
YDEAERC_TEGEN%RAERV(: ,1979, 6) = (/&
 & 23,  23,  23,  23,  23,  23,  23,  23,  23,  23,  23,  23,&
 & 24,  24,  24,  24,  24,  24,  24,  24,  24,  24,  24,  24/)  
YDEAERC_TEGEN%RAERV(: ,1979, 7) = (/&
 & 22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,&
 & 26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26/)  
YDEAERC_TEGEN%RAERV(: ,1979, 8) = (/&
 & 22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,&
 & 27,  27,  27,  27,  27,  27,  27,  27,  27,  27,  27,  27/)  
YDEAERC_TEGEN%RAERV(: ,1979, 9) = (/&
 & 23,  23,  23,  23,  23,  23,  23,  23,  23,  23,  23,  23,&
 & 29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29/)  
YDEAERC_TEGEN%RAERV(: ,1979,10) = (/&
 & 24,  24,  24,  24,  24,  24,  24,  24,  24,  24,  24,  24,&
 & 26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26/)  
YDEAERC_TEGEN%RAERV(: ,1979,11) = (/&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25/)  
YDEAERC_TEGEN%RAERV(: ,1979,12) = (/&
 & 30,  30,  30,  30,  30,  30,  30,  30,  30,  30,  30,  30,&
 & 25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25/)  
YDEAERC_TEGEN%RAERV(: ,1980, 1) = (/&
 & 33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,  33,&
 & 27,  27,  27,  27,  27,  27,  27,  27,  27,  27,  27,  27/)  
YDEAERC_TEGEN%RAERV(: ,1980, 2) = (/&
 & 22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,&
 & 29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29/)  
YDEAERC_TEGEN%RAERV(: ,1980, 3) = (/&
 & 34,  34,  34,  34,  34,  34,  34,  34,  34,  34,  34,  34,&
 & 29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29/)  
YDEAERC_TEGEN%RAERV(: ,1980, 4) = (/&
 & 31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31,&
 & 29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29,  29/)  
YDEAERC_TEGEN%RAERV(: ,1980, 5) = (/&
 & 30,  30,  30,  30,  30,  30,  30,  30,  30,  30,  30,  30,&
 & 30,  30,  30,  30,  30,  30,  30,  30,  30,  30,  30,  30/)  
YDEAERC_TEGEN%RAERV(: ,1980, 6) = (/&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,&
 & 31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31/)  
YDEAERC_TEGEN%RAERV(: ,1980, 7) = (/&
 & 55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,&
 & 31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31/)  

IF (LHOOK) CALL DR_HOOK('YOEAERC:PART3',1,ZHOOK_HANDLE)
END SUBROUTINE PART3


SUBROUTINE PART4

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOEAERC:PART4',0,ZHOOK_HANDLE)

YDEAERC_TEGEN%RAERV(: ,1980, 8) = (/&
 & 65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65,&
 & 31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31,  31/)  
YDEAERC_TEGEN%RAERV(: ,1980, 9) = (/&
 & 68,  68,  68,  68,  68,  68,  68,  68,  68,  68,  68,  68,&
 & 32,  32,  32,  32,  32,  32,  32,  32,  32,  32,  32,  32/)  
YDEAERC_TEGEN%RAERV(: ,1980,10) = (/&
 & 68,  68,  68,  68,  68,  68,  68,  68,  68,  68,  68,  68,&
 & 34,  34,  34,  34,  34,  34,  34,  34,  34,  34,  34,  34/)  
YDEAERC_TEGEN%RAERV(: ,1980,11) = (/&
 & 65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65,  65,&
 & 41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41/)  
YDEAERC_TEGEN%RAERV(: ,1980,12) = (/&
 & 61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,&
 & 41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41,  41/)  
YDEAERC_TEGEN%RAERV(: ,1981, 1) = (/&
 & 59,  59,  59,  59,  59,  59,  59,  59,  59,  59,  59,  59,&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42/)  
YDEAERC_TEGEN%RAERV(: ,1981, 2) = (/&
 & 55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42/)  
YDEAERC_TEGEN%RAERV(: ,1981, 3) = (/&
 & 49,  49,  49,  49,  49,  49,  49,  49,  49,  49,  49,  49,&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42/)  
YDEAERC_TEGEN%RAERV(: ,1981, 4) = (/&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42/)  
YDEAERC_TEGEN%RAERV(: ,1981, 5) = (/&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,&
 & 42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42/)  
YDEAERC_TEGEN%RAERV(: ,1981, 6) = (/&
 & 52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52,&
 & 44,  44,  44,  44,  44,  44,  44,  44,  44,  44,  44,  44/)  
YDEAERC_TEGEN%RAERV(: ,1981, 7) = (/&
 & 61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,&
 & 48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48/)  
YDEAERC_TEGEN%RAERV(: ,1981, 8) = (/&
 & 57,  57,  57,  57,  57,  57,  57,  57,  57,  57,  57,  57,&
 & 51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51/)  
YDEAERC_TEGEN%RAERV(: ,1981, 9) = (/&
 & 53,  53,  53,  53,  53,  53,  53,  53,  53,  53,  53,  53,&
 & 48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48,  48/)  
YDEAERC_TEGEN%RAERV(: ,1981,10) = (/&
 & 54,  54,  54,  54,  54,  54,  54,  54,  54,  54,  54,  54,&
 & 46,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46/)  
YDEAERC_TEGEN%RAERV(: ,1981,11) = (/&
 & 55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,&
 & 46,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46/)  
YDEAERC_TEGEN%RAERV(: ,1981,12) = (/&
 & 55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,&
 & 46,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46/)  
YDEAERC_TEGEN%RAERV(: ,1982, 1) = (/&
 & 61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,&
 & 51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51/)  
YDEAERC_TEGEN%RAERV(: ,1982, 2) = (/&
 & 61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,&
 & 51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51/)  
YDEAERC_TEGEN%RAERV(: ,1982, 3) = (/&
 & 61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,  61,&
 & 51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51,  51/)  
YDEAERC_TEGEN%RAERV(: ,1982, 4) = (/&
 & 101, 101, 101, 101, 101,  23,  20,  13, 123,2823, 124,   8,&
 & 6,   6,   4,   3,   2,   2,   1,  66,  66,  66,  66,  66/)  
YDEAERC_TEGEN%RAERV(: ,1982, 5) = (/&
 & 185, 185, 185, 185, 185, 224, 194, 146,1158,3001,1271, 102,&
 & 77,  58,  44,  35,  29,  23,  17,  55,  55,  55,  55,  55/)  
YDEAERC_TEGEN%RAERV(: ,1982, 6) = (/&
 & 253, 253, 253, 253, 253, 274, 238, 279,1266,2642,1938, 443,&
 & 267, 163,  99,  78,  73,  68,  62,  66,  66,  66,  66,  66/)  
YDEAERC_TEGEN%RAERV(: ,1982, 7) = (/&
 & 284, 284, 284, 284, 284, 299, 263, 437,1184,2470,2386,1285,&
 & 664, 354, 187, 144, 151, 160, 163,  96,  96,  96,  96,  96/)  
YDEAERC_TEGEN%RAERV(: ,1982, 8) = (/&
 & 347, 347, 347, 347, 347, 377, 352, 504,1242,2180,2270,1575,&
 & 829, 468, 260, 197, 209, 227, 221, 143, 143, 143, 143, 143/)  
YDEAERC_TEGEN%RAERV(: ,1982, 9) = (/&
 & 387, 387, 387, 387, 387, 487, 485, 557,1319,1916,2074,1719,&
 & 938, 578, 345, 256, 273, 301, 275,  96,  96,  96,  96,  96/)  
YDEAERC_TEGEN%RAERV(: ,1982,10) = (/&
 & 476, 476, 476, 476, 476, 630, 670, 618,1400,1681,1895,1873,&
 & 1060, 715, 459, 337, 359, 402, 342,  46,  46,  46,  46,  46/)  
YDEAERC_TEGEN%RAERV(: ,1982,11) = (/&
 & 689, 689, 689, 689, 689, 882, 900, 773,1218,1404,1738,1834,&
 & 1188, 826, 531, 343, 402, 473, 430, 195, 195, 195, 195, 195/)  
YDEAERC_TEGEN%RAERV(: ,1982,12) = (/&
 & 1085,1085,1085,1085,1085,1240,1126,1016, 890,1134,1575,1601,&
 & 1255, 865, 532, 287, 397, 497, 526, 461, 461, 461, 461, 461/)  
YDEAERC_TEGEN%RAERV(: ,1983, 1) = (/&
 & 1214,1214,1214,1214,1214,1300,1143,1018, 803,1015,1364,1344,&
 & 1053, 755, 494, 310, 449, 519, 545, 549, 549, 549, 549, 549/)  
YDEAERC_TEGEN%RAERV(: ,1983, 2) = (/&
 & 1377,1377,1377,1377,1377,1315,1132, 976, 747, 923,1176,1126,&
 & 853, 647, 452, 352, 514, 534, 545, 575, 575, 575, 575, 575/)  
YDEAERC_TEGEN%RAERV(: ,1983, 3) = (/&
 & 2013,2013,2013,2013,2013,1341,1128, 932, 688, 832,1010, 941,&
 & 693, 558, 413, 400, 589, 545, 540, 487, 487, 487, 487, 487/)  
YDEAERC_TEGEN%RAERV(: ,1983, 4) = (/&
 & 1681,1681,1681,1681,1681,1365,1124, 892, 636, 749, 864, 780,&
 & 558, 476, 372, 458, 675, 553, 536, 410, 410, 410, 410, 410/)  
YDEAERC_TEGEN%RAERV(: ,1983, 5) = (/&
 & 1460,1460,1460,1460,1460,1359,1093, 846, 590, 678, 746, 656,&
 & 457, 409, 337, 502, 740, 552, 526, 331, 331, 331, 331, 331/)  
YDEAERC_TEGEN%RAERV(: ,1983, 6) = (/&
 & 1572,1572,1572,1572,1572,1200, 943, 748, 542, 620, 674, 598,&
 & 419, 367, 296, 447, 644, 491, 485, 279, 279, 279, 279, 279/)  
YDEAERC_TEGEN%RAERV(: ,1983, 7) = (/&
 & 1294,1294,1294,1294,1294,1035, 793, 652, 498, 568, 615, 554,&
 & 389, 333, 261, 381, 536, 428, 440, 290, 290, 290, 290, 290/)  
YDEAERC_TEGEN%RAERV(: ,1983, 8) = (/&
 & 1038,1038,1038,1038,1038, 890, 663, 569, 459, 518, 560, 513,&
 & 365, 301, 229, 327, 444, 372, 402, 329, 329, 329, 329, 329/)  
YDEAERC_TEGEN%RAERV(: ,1983, 9) = (/&
 & 1012,1012,1012,1012,1012, 768, 558, 496, 420, 473, 510, 475,&
 & 341, 273, 200, 279, 371, 323, 366, 188, 188, 188, 188, 188/)  
YDEAERC_TEGEN%RAERV(: ,1983,10) = (/&
 & 818, 818, 818, 818, 818, 662, 468, 433, 385, 432, 465, 440,&
 & 319, 246, 176, 238, 308, 282, 333, 122, 122, 122, 122, 122/)  
YDEAERC_TEGEN%RAERV(: ,1983,11) = (/&
 & 667, 667, 667, 667, 667, 571, 393, 378, 354, 396, 424, 407,&
 & 296, 223, 154, 204, 256, 245, 303, 186, 186, 186, 186, 186/)  
YDEAERC_TEGEN%RAERV(: ,1983,12) = (/&
 & 567, 567, 567, 567, 567, 492, 330, 330, 326, 362, 386, 377,&
 & 277, 204, 135, 174, 213, 213, 275, 294, 294, 294, 294, 294/)  
YDEAERC_TEGEN%RAERV(: ,1984, 1) = (/&
 & 540, 540, 540, 540, 540, 424, 276, 287, 297, 329, 351, 349,&
 & 259, 183, 119, 147, 177, 185, 251, 310, 310, 310, 310, 310/)  
YDEAERC_TEGEN%RAERV(: ,1984, 2) = (/&
 & 503, 503, 503, 503, 503, 375, 242, 256, 268, 298, 318, 318,&
 & 237, 165, 106, 131, 155, 165, 227, 293, 293, 293, 293, 293/)  
YDEAERC_TEGEN%RAERV(: ,1984, 3) = (/&
 & 385, 385, 385, 385, 385, 338, 218, 230, 243, 268, 286, 286,&
 & 212, 147,  96, 118, 140, 149, 205, 263, 263, 263, 263, 263/)  
YDEAERC_TEGEN%RAERV(: ,1984, 4) = (/&
 & 433, 433, 433, 433, 433, 303, 197, 207, 217, 242, 257, 257,&
 & 191, 133,  86, 107, 125, 133, 184, 246, 246, 246, 246, 246/)  
YDEAERC_TEGEN%RAERV(: ,1984, 5) = (/&
 & 388, 388, 388, 388, 388, 273, 177, 186, 196, 217, 231, 231,&
 & 173, 120,  77,  95, 112, 119, 165, 217, 217, 217, 217, 217/)  
YDEAERC_TEGEN%RAERV(: ,1984, 6) = (/&
 & 321, 321, 321, 321, 321, 244, 158, 167, 176, 196, 208, 208,&
 & 154, 108,  69,  85, 101, 107, 147, 249, 249, 249, 249, 249/)  
YDEAERC_TEGEN%RAERV(: ,1984, 7) = (/&
 & 284, 284, 284, 284, 284, 220, 142, 150, 157, 175, 187, 187,&
 & 140,  98,  62,  76,  90,  96, 133, 276, 276, 276, 276, 276/)  
YDEAERC_TEGEN%RAERV(: ,1984, 8) = (/&
 & 274, 274, 274, 274, 274, 197, 127, 134, 142, 157, 167, 167,&
 & 124,  86,  56,  67,  81,  86, 119, 227, 227, 227, 227, 227/)  
YDEAERC_TEGEN%RAERV(: ,1984, 9) = (/&
 & 264, 264, 264, 264, 264, 177, 114, 121, 127, 141, 151, 151,&
 & 112,  78,  51,  61,  73,  77, 107, 114, 114, 114, 114, 114/)  
YDEAERC_TEGEN%RAERV(: ,1984,10) = (/&
 & 256, 256, 256, 256, 256, 158, 102, 109, 114, 127, 135, 135,&
 & 101,  70,  44,  55,  66,  69,  95,  84,  84,  84,  84,  84/)  
YDEAERC_TEGEN%RAERV(: ,1984,11) = (/&
 & 249, 249, 249, 249, 249, 142,  92,  98, 102, 114, 122, 122,&
 & 90,  63,  40,  50,  59,  62,  86, 139, 139, 139, 139, 139/)  
YDEAERC_TEGEN%RAERV(: ,1984,12) = (/&
 & 241, 241, 241, 241, 241, 128,  83,  88,  92, 102, 110, 110,&
 & 81,  56,  35,  44,  53,  55,  77, 199, 199, 199, 199, 199/)  
YDEAERC_TEGEN%RAERV(: ,1985, 1) = (/&
 & 229, 229, 229, 229, 229, 256, 251, 186, 134, 120, 133, 134,&
 & 140, 145, 128, 104, 112, 145, 182, 213, 213, 213, 213, 213/)  
YDEAERC_TEGEN%RAERV(: ,1985, 2) = (/&
 & 237, 237, 237, 237, 237, 268, 262, 178, 127, 115, 119, 117,&
 & 125, 127, 102, 107, 111, 133, 212, 213, 211, 213, 215, 215/)  
YDEAERC_TEGEN%RAERV(: ,1985, 3) = (/&
 & 222, 222, 222, 222, 223, 220, 202, 122,  88,  90, 103, 115,&
 & 116, 116, 112, 109, 110, 121, 157, 208, 205, 202, 202, 202/)  
YDEAERC_TEGEN%RAERV(: ,1985, 4) = (/&
 & 234, 234, 234, 240, 238, 191, 155, 115,  93,  92,  97, 117,&
 & 123, 115, 105, 107, 111, 134, 165, 192, 192, 192, 192, 192/)  
YDEAERC_TEGEN%RAERV(: ,1985, 5) = (/&
 & 210, 210, 210, 210, 210, 155, 138, 108,  88,  84,  98, 109,&
 & 112, 105,  96, 100, 115, 149, 161, 176, 176, 176, 176, 176/)  
YDEAERC_TEGEN%RAERV(: ,1985, 6) = (/&
 & 172, 172, 172, 172, 173, 151, 118,  94,  83,  86, 103, 112,&
 & 116, 105,  92,  95, 121, 171, 178, 182, 182, 182, 182, 182/)  
YDEAERC_TEGEN%RAERV(: ,1985, 7) = (/&
 & 150, 150, 150, 150, 150, 122,  98,  79,  79,  88, 109, 115,&
 & 122,  98,  79,  86, 120, 161, 176, 205, 205, 205, 205, 205/)  
YDEAERC_TEGEN%RAERV(: ,1985, 8) = (/&
 & 152, 152, 152, 152, 152, 119,  97,  82,  78,  86, 103, 113,&
 & 118,  93,  73,  78, 107, 144, 165, 191, 191, 191, 191, 191/)  
YDEAERC_TEGEN%RAERV(: ,1985, 9) = (/&
 & 155, 155, 155, 155, 156, 116,  96,  86,  77,  84,  97, 110,&
 & 115,  88,  67,  69,  93, 127, 154,  85,  89,  89,  89,  89/)  
YDEAERC_TEGEN%RAERV(: ,1985,10) = (/&
 & 157, 157, 157, 157, 157, 125,  97,  77,  76,  84,  93, 107,&
 & 108,  86,  73,  80, 105, 126, 130,  75,  76,  73,  73,  73/)  
YDEAERC_TEGEN%RAERV(: ,1985,11) = (/&
 & 146, 146, 146, 146, 146, 133, 111,  83,  77,  78,  94, 105,&
 & 104,  81,  70,  63,  76, 107, 132,  87,  85,  85,  85,  85/)  
YDEAERC_TEGEN%RAERV(: ,1985,12) = (/&
 & 150, 150, 150, 150, 150, 151, 134, 102,  77,  79, 139, 254,&
 & 188,  98,  71,  62,  71,  93, 113, 116, 116, 116, 116, 116/)  
YDEAERC_TEGEN%RAERV(: ,1986, 1) = (/&
 & 178, 178, 178, 178, 178, 168, 157, 121,  76,  80, 185, 404,&
 & 271, 116,  71,  60,  66,  79,  94, 120, 120, 120, 120, 120/)  
YDEAERC_TEGEN%RAERV(: ,1986, 2) = (/&
 & 224, 224, 224, 224, 224, 184, 182, 138, 104, 111, 239, 389,&
 & 284, 115,  32,  64,  67,  83,  98, 114, 114, 115, 118, 118/)  
YDEAERC_TEGEN%RAERV(: ,1986, 3) = (/&
 & 152, 152, 152, 152, 150, 160, 140,  86,  81,  94, 197, 364,&
 & 290, 140,  81,  67,  67,  86,  97, 121, 118, 117, 117, 117/)  
YDEAERC_TEGEN%RAERV(: ,1986, 4) = (/&
 & 173, 173, 173, 181, 179, 149, 140, 110,  94, 120, 178, 354,&
 & 273, 138,  81,  68,  69,  89, 107, 116, 116, 116, 116, 116/)  
YDEAERC_TEGEN%RAERV(: ,1986, 5) = (/&
 & 172, 172, 172, 172, 176, 130, 113,  94,  91, 110, 186, 301,&
 & 271, 136,  79,  67,  76,  96, 117, 116, 116, 116, 116, 116/)  
YDEAERC_TEGEN%RAERV(: ,1986, 6) = (/&
 & 149, 149, 149, 149, 150, 130,  99,  79,  76,  96, 163, 251,&
 & 240, 125,  70,  71,  92, 119, 131, 126, 126, 126, 126, 126/)  
YDEAERC_TEGEN%RAERV(: ,1986, 7) = (/&
 & 135, 135, 135, 135, 135, 116,  85,  63,  62,  83, 140, 202,&
 & 208, 138,  87,  78, 107, 122, 134, 201, 201, 201, 201, 201/)  
YDEAERC_TEGEN%RAERV(: ,1986, 8) = (/&
 & 141, 141, 141, 140, 139, 101,  64,  71,  78,  95, 117, 188,&
 & 194, 145, 110,  86, 122, 144, 150, 303, 303, 303, 303, 303/)  
YDEAERC_TEGEN%RAERV(: ,1986, 9) = (/&
 & 149, 149, 149, 151, 151, 112,  94,  80,  93, 109, 144, 181,&
 & 178, 155,  79,  70, 104, 129, 140, 139, 144, 144, 144, 144/)  
YDEAERC_TEGEN%RAERV(: ,1986,10) = (/&
 & 148, 148, 148, 148, 148, 118, 100,  82,  87, 114, 169, 207,&
 & 186, 153, 105,  86, 100, 136, 138,  59,  63,  62,  62,  62/)  
YDEAERC_TEGEN%RAERV(: ,1986,11) = (/&
 & 144, 144, 144, 144, 144, 141, 118,  97,  84, 113, 159, 184,&
 & 181, 151,  95,  73,  85,  97, 111,  92,  94,  94,  94,  94/)  
YDEAERC_TEGEN%RAERV(: ,1986,12) = (/&
 & 161, 161, 161, 161, 161, 152, 136, 107,  87, 106, 150, 162,&
 & 160, 135,  85,  51,  60,  86, 105, 119, 120, 120, 120, 120/)  
YDEAERC_TEGEN%RAERV(: ,1987, 1) = (/&
 & 247, 247, 247, 247, 247, 159, 105, 113,  91, 112, 136, 141,&
 & 138, 119,  74,  28,  35,  76,  93, 118, 118, 118, 118, 118/)  
YDEAERC_TEGEN%RAERV(: ,1987, 2) = (/&
 & 187, 187, 187, 187, 187, 157, 155, 119, 104, 104, 119, 121,&
 & 116,  82,  76,  49,  51,  77, 101, 112, 110, 113, 112, 112/)  
YDEAERC_TEGEN%RAERV(: ,1987, 3) = (/&
 & 176, 176, 176, 176, 176, 162, 145,  98,  83,  87, 110, 113,&
 & 108, 100,  77,  69,  68,  78,  82, 128, 128, 128, 128, 128/)  
YDEAERC_TEGEN%RAERV(: ,1987, 4) = (/&
 & 164, 164, 164, 167, 170, 135, 124, 104,  86,  89,  96, 105,&
 & 105,  97,  76,  70,  69,  83,  93, 115, 115, 115, 115, 115/)  
YDEAERC_TEGEN%RAERV(: ,1987, 5) = (/&
 & 158, 158, 158, 158, 152, 148, 118, 100,  92,  90, 102,  95,&
 & 95,  91,  80,  65,  73,  82,  93, 105, 105, 105, 105, 105/)  
YDEAERC_TEGEN%RAERV(: ,1987, 6) = (/&
 & 151, 151, 151, 151, 153, 119, 104,  90,  85,  86,  95,  91,&
 & 89,  87,  91,  83,  95, 115, 115, 142, 142, 142, 142, 142/)  
YDEAERC_TEGEN%RAERV(: ,1987, 7) = (/&
 & 124, 124, 124, 124, 126, 108,  89,  80,  79,  82,  90,  86,&
 & 83,  87,  83,  88, 102, 114, 118, 230, 230, 230, 230, 230/)  
YDEAERC_TEGEN%RAERV(: ,1987, 8) = (/&
 & 132, 132, 129, 129, 130, 102,  91,  81,  80,  83,  80,  85,&
 & 87,  84,  85,  86, 109, 125, 119, 230, 230, 230, 230, 230/)  
YDEAERC_TEGEN%RAERV(: ,1987, 9) = (/&
 & 149, 149, 149, 149, 149, 110,  93,  82,  82,  83,  82,  88,&
 & 85,  79,  71,  80, 100, 121, 114, 180, 179, 179, 179, 179/)  
YDEAERC_TEGEN%RAERV(: ,1987,10) = (/&
 & 140, 140, 140, 140, 140, 119,  99,  81,  82,  76,  77,  86,&
 & 86,  73,  68,  72, 102, 122, 120,  86,  91,  82,  82,  82/)  
YDEAERC_TEGEN%RAERV(: ,1987,11) = (/&
 & 144, 144, 144, 144, 144, 125, 111,  83,  70,  67,  74,  83,&
 & 83,  68,  66,  74,  89, 104, 109,  70,  73,  73,  73,  73/)  
YDEAERC_TEGEN%RAERV(: ,1987,12) = (/&
 & 139, 139, 139, 139, 139, 125, 115,  82,  63,  59,  68,  82,&
 & 79,  72,  63,  65,  76, 102,  99, 110, 114, 114, 114, 114/)  
YDEAERC_TEGEN%RAERV(: ,1988, 1) = (/&
 & 139, 139, 139, 139, 139, 136, 113,  80,  54,  57,  67,  80,&
 & 75,  75,  61,  56,  63,  74,  87, 114, 115, 115, 115, 115/)  
YDEAERC_TEGEN%RAERV(: ,1988, 2) = (/&
 & 177, 177, 177, 177, 177, 135, 123,  88,  52,  47,  66,  76,&
 & 73,  64,  63,  60,  64,  72,  81, 113, 111, 112, 112, 112/)  
YDEAERC_TEGEN%RAERV(: ,1988, 3) = (/&
 & 127, 127, 127, 127, 126, 112,  94,  66,  45,  41,  57,  74,&
 & 73,  68,  64,  63,  66,  71,  75, 113, 113, 113, 113, 113/)  
YDEAERC_TEGEN%RAERV(: ,1988, 4) = (/&
 & 122, 122, 122, 124, 129,  91,  82,  61,  42,  44,  58,  66,&
 & 68,  54,  60,  58,  62,  73,  88, 108, 108, 108, 108, 108/)  
YDEAERC_TEGEN%RAERV(: ,1988, 5) = (/&
 & 103, 103, 103, 103, 102,  84,  74,  54,  41,  44,  55,  66,&
 & 67,  63,  57,  56,  67,  82,  94,  97,  97,  97,  97,  97/)  
YDEAERC_TEGEN%RAERV(: ,1988, 6) = (/&
 & 81,  81,  81,  81,  83,  70,  49,  46,  40,  44,  56,  65,&
 & 67,  70,  61,  67,  79,  96, 101, 110, 110, 110, 110, 110/)  
YDEAERC_TEGEN%RAERV(: ,1988, 7) = (/&
 & 68,  68,  68,  68,  69,  54,  46,  39,  40,  43,  57,  64,&
 & 67,  65,  58,  58,  79,  98, 104, 243, 243, 243, 243, 243/)  

IF (LHOOK) CALL DR_HOOK('YOEAERC:PART4',1,ZHOOK_HANDLE)
END SUBROUTINE PART4


SUBROUTINE PART5

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOEAERC:PART5',0,ZHOOK_HANDLE)

YDEAERC_TEGEN%RAERV(: ,1988, 8) = (/&
 & 77,  77,  74,  72,  73,  62,  50,  42,  41,  45,  53,  65,&
 & 68,  66,  58,  65,  90, 118, 113, 378, 378, 378, 378, 378/)  
YDEAERC_TEGEN%RAERV(: ,1988, 9) = (/&
 & 94,  94,  94,  94,  96,  65,  54,  45,  43,  47,  57,  66,&
 & 70,  61,  51,  56,  79,  94, 104,  95,  95,  95,  95,  95/)  
YDEAERC_TEGEN%RAERV(: ,1988,10) = (/&
 & 90,  90,  90,  90,  90,  67,  59,  45,  41,  44,  53,  60,&
 & 63,  62,  56,  57,  73, 100, 105,  67,  63,  62,  62,  62/)  
YDEAERC_TEGEN%RAERV(: ,1988,11) = (/&
 & 88,  88,  88,  88,  88,  76,  67,  50,  43,  45,  53,  59,&
 & 63,  63,  56,  54,  65,  80,  91, 100, 100, 100, 100, 100/)  
YDEAERC_TEGEN%RAERV(: ,1988,12) = (/&
 & 84,  84,  84,  84,  84,  82,  79,  57,  46,  49,  50,  60,&
 & 60,  65,  54,  51,  59,  70,  86, 114, 112, 112, 112, 112/)  
YDEAERC_TEGEN%RAERV(: ,1989, 1) = (/&
 & 160, 160, 160, 160, 160,  88,  77,  57,  41,  44,  49,  61,&
 & 58,  67,  51,  49,  53,  63,  76, 100, 100, 100, 100, 100/)  
YDEAERC_TEGEN%RAERV(: ,1989, 2) = (/&
 & 83,  83,  83,  83,  83,  91,  83,  61,  43,  42,  49,  53,&
 & 54,  46,  51,  51,  54,  63,  76,  98,  96,  96,  94,  94/)  
YDEAERC_TEGEN%RAERV(: ,1989, 3) = (/&
 & 74,  74,  74,  73,  76,  73,  63,  53,  37,  40,  47,  55,&
 & 60,  57,  52,  53,  55,  64,  75, 100,  98,  95,  95,  95/)  
YDEAERC_TEGEN%RAERV(: ,1989, 4) = (/&
 & 79,  79,  79,  81,  82,  72,  67,  49,  44,  42,  46,  51,&
 & 54,  51,  48,  47,  49,  61,  71,  93,  93,  93,  93,  93/)  
YDEAERC_TEGEN%RAERV(: ,1989, 5) = (/&
 & 67,  67,  67,  67,  66,  61,  57,  42,  38,  37,  44,  49,&
 & 2,  12,   1,  22,  49,  68,  79, 102, 102, 102, 102, 102/)  
YDEAERC_TEGEN%RAERV(: ,1989, 6) = (/&
 & 80,  80,  80,  80,  81,  56,  48,  38,  36,  37,  43,  46,&
 & 24,  55,  51,  49,  63,  76,  79,  97,  97,  97,  97,  97/)  
YDEAERC_TEGEN%RAERV(: ,1989, 7) = (/&
 & 80,  80,  80,  80,  78,  43,  40,  34,  34,  37,  43,  44,&
 & 46,  50,  49,  52,  65,  80,  86, 203, 203, 203, 203, 203/)  
YDEAERC_TEGEN%RAERV(: ,1989, 8) = (/&
 & 69,  69,  67,  65,  65,  42,  42,  36,  38,  42,  37,  40,&
 & 47,  50,  52,  57,  73,  89,  90, 369, 369, 369, 369, 369/)  
YDEAERC_TEGEN%RAERV(: ,1989, 9) = (/&
 & 71,  71,  71,  71,  71,  40,  45,  39,  42,  46,  45,  43,&
 & 42,  46,  46,  42,  53,  79,  79, 161, 162, 162, 162, 162/)  
YDEAERC_TEGEN%RAERV(: ,1989,10) = (/&
 & 71,  71,  71,  71,  71,  55,  48,  40,  42,  43,  42,  41,&
 & 41,  46,  48,  50,  65,  76,  81,  63,  60,  58,  58,  58/)  
YDEAERC_TEGEN%RAERV(: ,1989,11) = (/&
 & 73,  73,  73,  73,  73,  64,  57,  47,  42,  42,  43,  45,&
 & 45,  47,  47,  48,  58,  63,  75,  63,  66,  66,  66,  66/)  
YDEAERC_TEGEN%RAERV(: ,1989,12) = (/&
 & 74,  74,  74,  74,  74,  73,  71,  47,  42,  36,  36,  40,&
 & 40,  44,  44,  43,  48,  56,  68,  77,  77,  77,  77,  77/)  
YDEAERC_TEGEN%RAERV(: ,1990, 1) = (/&
 & 157, 157, 157, 157, 157,  81,  72,  57,  38,  35,  35,  34,&
 & 36,  41,  41,  39,  40,  45,  59,  74,  77,  77,  77,  77/)  
YDEAERC_TEGEN%RAERV(: ,1990, 2) = (/&
 & 313, 313, 313, 313, 313,  76,  69,  56,  41,  33,  31,  33,&
 & 49,  28,  58,  44,  42,  49,  59,  75,  77,  77,  77,  77/)  
YDEAERC_TEGEN%RAERV(: ,1990, 3) = (/&
 & 77,  77,  77,  67,  68,  68,  59,  48,  37,  33,  41,  58,&
 & 70, 107,  74,  48,  45,  52,  59,  84,  82,  82,  82,  82/)  
YDEAERC_TEGEN%RAERV(: ,1990, 4) = (/&
 & 77,  77,  77,  75,  75,  71,  64,  50,  39,   7,   1,   1,&
 & 1,  61,  78,  56,  25,  28,  49,  73,  73,  73,  73,  73/)  
YDEAERC_TEGEN%RAERV(: ,1990, 5) = (/&
 & 74,  74,  74,  74,  72,  55,  50,  39,  36,  37,  39,  39,&
 & 1,  13,  22,  28,  32,  75,  84,  65,  65,  65,  65,  65/)  
YDEAERC_TEGEN%RAERV(: ,1990, 6) = (/&
 & 74,  74,  74,  74,  74,  51,  40,  35,  34,  37,  41,  46,&
 & 32,  72,  61,  56,  76,  91,  94,  93,  93,  93,  93,  93/)  
YDEAERC_TEGEN%RAERV(: ,1990, 7) = (/&
 & 68,  68,  68,  68,  65,  62,  36,  30,  33,  38,  43,  53,&
 & 64,  64,  55,  62,  82,  95, 105, 209, 209, 209, 209, 209/)  
YDEAERC_TEGEN%RAERV(: ,1990, 8) = (/&
 & 63,  63,  63,  62,  61,  60,  40,  34,  35,  40,  47,  56,&
 & 69,  66,  60,  62,  84, 101, 107, 203, 203, 203, 203, 203/)  
YDEAERC_TEGEN%RAERV(: ,1990, 9) = (/&
 & 64,  64,  64,  64,  66,  58,  44,  37,  38,  42,  51,  57,&
 & 66,  58,  47,  51,  64,  87,  92, 162, 171, 171, 171, 171/)  
YDEAERC_TEGEN%RAERV(: ,1990,10) = (/&
 & 73,  73,  73,  73,  73,  57,  47,  40,  38,  40,  47,  57,&
 & 63,  57,  49,  55,  68,  86,  94,  52,  54,  50,  50,  50/)  
YDEAERC_TEGEN%RAERV(: ,1990,11) = (/&
 & 73,  73,  73,  73,  73,  60,  55,  42,  39,  43,  48,  59,&
 & 65,  58,  51,  54,  65,  78,  87,  61,  66,  66,  66,  66/)  
YDEAERC_TEGEN%RAERV(: ,1990,12) = (/&
 & 79,  79,  79,  79,  79,  67,  62,  48,  35,  37,  39,  59,&
 & 65,  58,  51,  48,  51,  66,  81,  81,  82,  82,  82,  82/)  
YDEAERC_TEGEN%RAERV(: ,1991, 1) = (/&
 & 102, 103, 110, 115,  87,  66,  64,  52,  38,  36,  43,  51,&
 & 53,  50,  44,  47,  49,  56,  71,  98,  91,  83,  89,  89/)  
YDEAERC_TEGEN%RAERV(: ,1991, 2) = (/&
 & 85,  86,  83,  75,  69,  68,  65,  54,  39,  36,  43,  49,&
 & 48,  45,  44,  44,  49,  56,  67,  91, 101,  90,  85,  84/)  
YDEAERC_TEGEN%RAERV(: ,1991, 3) = (/&
 & 73,  73,  72,  68,  65,  59,  51,  42,  34,  34,  41,  49,&
 & 51,  49,  47,  48,  51,  58,  68,  83, 111, 109,  92,  90/)  
YDEAERC_TEGEN%RAERV(: ,1991, 4) = (/&
 & 77,  76,  71,  70,  66,  59,  52,  43,  36,  34,  38,  43,&
 & 45,  44,  43,  43,  47,  55,  65,  78, 105, 105,  90,  88/)  
YDEAERC_TEGEN%RAERV(: ,1991, 5) = (/&
 & 118, 114,  63,  59,  57,  52,  45,  38,  33,  33,  40,  47,&
 & 47,  47,  44,  48,  56,  64,  68,  87, 100, 101,  86,  84/)  
YDEAERC_TEGEN%RAERV(: ,1991, 6) = (/&
 & 168, 160,  66,  57,  55,  52,  51,  45,  59, 164, 378, 720,&
 & 558, 175,  60,  51,  63,  70,  74,  98,  99, 110,  97,  95/)  
YDEAERC_TEGEN%RAERV(: ,1991, 7) = (/&
 & 228, 224, 161, 138, 126, 118, 103,  88, 213, 498, 673,1202,&
 & 1218, 780, 146,  78,  93,  92,  89, 134, 151, 171, 111, 105/)  
YDEAERC_TEGEN%RAERV(: ,1991, 8) = (/&
 & 342, 338, 274, 235, 225, 216, 193, 200, 480, 956,1083,1729,&
 & 1938,1826, 721, 332, 329, 272, 191, 152, 168, 183, 121, 115/)  
YDEAERC_TEGEN%RAERV(: ,1991, 9) = (/&
 & 366, 366, 368, 330, 336, 325, 312, 398, 883,1489,1634,2008,&
 & 2263,2023, 990, 543, 571, 565, 450, 217, 186, 256, 157, 149/)  
YDEAERC_TEGEN%RAERV(: ,1991,10) = (/&
 & 547, 559, 748, 531, 515, 529, 552, 613,1068,1773,1947,2163,&
 & 2292,2090,1168, 780, 900,1019, 871, 609, 547, 442, 206, 189/)  
YDEAERC_TEGEN%RAERV(: ,1991,11) = (/&
 & 715, 732,1010, 693, 685, 768, 804, 777,1125,1716,1959,2124,&
 & 2301,2558,1457, 826,1020,1242,1268,1240,1057, 944, 446, 406/)  
YDEAERC_TEGEN%RAERV(: ,1991,12) = (/&
 & 833, 854,1201,1085,1239,1635,1578,1356,1226,1606,1612,1748,&
 & 1977,2073,1161, 720, 785, 924,1173,1413,1220,1210,1236,1236/)  
YDEAERC_TEGEN%RAERV(: ,1992, 1) = (/&
 & 924, 949,1358,1392,1564,1731,1761,1665,1330,1539,1837,1873,&
 & 2006,2067,1269, 708, 717, 854,1158,1423,1462,1600,1668,1668/)  
YDEAERC_TEGEN%RAERV(: ,1992, 2) = (/&
 & 1003,1031,1454,1573,1677,1801,1837,1788,1444,1520,1846,1880,&
 & 1863,1911,1213, 740, 729, 871,1162,1418,1553,1559,1492,1489/)  
YDEAERC_TEGEN%RAERV(: ,1992, 3) = (/&
 & 1062,1089,1449,1524,1571,1674,1693,1581,1370,1485,1705,1760,&
 & 1765,1699,1224, 826, 758, 886,1134,1350,1602,1635,1472,1459/)  
YDEAERC_TEGEN%RAERV(: ,1992, 4) = (/&
 & 1120,1145,1468,1496,1461,1480,1465,1327,1285,1427,1556,1637,&
 & 1681,1672,1310, 911, 885,1025,1242,1361,1598,1710,1442,1420/)  
YDEAERC_TEGEN%RAERV(: ,1992, 5) = (/&
 & 1179,1193,1392,1422,1428,1404,1365,1298,1256,1435,1560,1624,&
 & 1613,1566,1214,1035,1117,1238,1294,1369,1472,1527,1378,1366/)  
YDEAERC_TEGEN%RAERV(: ,1992, 6) = (/&
 & 1233,1241,1334,1334,1304,1237,1110, 971, 995,1160,1294,1394,&
 & 1414,1383,1149,1081,1295,1464,1436,1395,1479,1454,1323,1316/)  
YDEAERC_TEGEN%RAERV(: ,1992, 7) = (/&
 & 1269,1270,1292,1299,1180,1117, 957, 836, 942,1175,1369,1415,&
 & 1306,1257,1099,1067,1317,1532,1533,1366,1427,1347,1272,1270/)  
YDEAERC_TEGEN%RAERV(: ,1992, 8) = (/&
 & 1272,1271,1251,1261,1203,1069, 956, 828, 912,1059,1230,1096,&
 & 1031,1098,1098,1078,1368,1532,1418,1263,1255,1177,1212,1216/)  
YDEAERC_TEGEN%RAERV(: ,1992, 9) = (/&
 & 1254,1253,1233,1231,1167,1041, 913, 779, 864,1066,1179, 969,&
 & 888, 980, 965, 933,1185,1338,1283,1175,1037,1023,1156,1166/)  
YDEAERC_TEGEN%RAERV(: ,1992,10) = (/&
 & 1193,1193,1183,1161,1147,1085, 942, 796, 829, 994,1012, 889,&
 & 824, 866, 832, 871,1070,1208,1228,1113,1018,1012,1111,1119/)  
YDEAERC_TEGEN%RAERV(: ,1992,11) = (/&
 & 1117,1117,1115,1120,1132,1100, 979, 820, 721, 783, 884, 862,&
 & 785, 821, 797, 803, 935,1063,1138,1136,1039, 989,1060,1066/)  
YDEAERC_TEGEN%RAERV(: ,1992,12) = (/&
 & 1044,1044,1050,1056,1143,1081, 984, 833, 663, 620, 667, 660,&
 & 647, 663, 580, 577, 686, 865,1027,1091,1000, 995,1015,1017/)  
YDEAERC_TEGEN%RAERV(: ,1993, 1) = (/&
 & 981, 982, 995, 996,1165,1123,1003, 823, 601, 520, 552, 610,&
 & 628, 599, 535, 533, 577, 751, 915,1026,1041,1013, 995, 994/)  
YDEAERC_TEGEN%RAERV(: ,1993, 2) = (/&
 & 904, 906, 926, 926, 985,1002, 927, 690, 482, 430, 496, 572,&
 & 563, 521, 462, 478, 519, 638, 789, 907,1006,1054,1042,1041/)  
YDEAERC_TEGEN%RAERV(: ,1993, 3) = (/&
 & 827, 831, 883, 887, 902, 886, 798, 610, 427, 388, 448, 518,&
 & 521, 485, 437, 443, 483, 607, 758, 853, 913,1021,1054,1052/)  
YDEAERC_TEGEN%RAERV(: ,1993, 4) = (/&
 & 749, 754, 817, 825, 833, 811, 725, 530, 383, 347, 407, 474,&
 & 465, 431, 410, 412, 466, 589, 721, 808, 839, 934,1006,1005/)  
YDEAERC_TEGEN%RAERV(: ,1993, 5) = (/&
 & 664, 664, 672, 704, 717, 702, 603, 436, 335, 323, 393, 457,&
 & 447, 424, 398, 418, 527, 649, 713, 744, 811, 881, 945, 945/)  
YDEAERC_TEGEN%RAERV(: ,1993, 6) = (/&
 & 579, 584, 641, 656, 644, 587, 478, 337, 285, 300, 348, 388,&
 & 402, 409, 356, 389, 527, 674, 728, 752, 780, 826, 883, 884/)  
YDEAERC_TEGEN%RAERV(: ,1993, 7) = (/&
 & 508, 509, 521, 512, 506, 470, 351, 275, 271, 290, 332, 369,&
 & 385, 372, 338, 375, 523, 704, 753, 757, 747, 771, 825, 826/)  
YDEAERC_TEGEN%RAERV(: ,1993, 8) = (/&
 & 505, 505, 506, 497, 466, 416, 331, 270, 262, 296, 357, 361,&
 & 378, 380, 322, 333, 463, 654, 720, 719, 696, 711, 756, 758/)  
YDEAERC_TEGEN%RAERV(: ,1993, 9) = (/&
 & 509, 509, 506, 459, 423, 372, 305, 258, 260, 299, 345, 357,&
 & 348, 323, 296, 319, 421, 554, 638, 660, 658, 667, 679, 679/)  
YDEAERC_TEGEN%RAERV(: ,1993,10) = (/&
 & 488, 489, 488, 432, 407, 373, 315, 254, 244, 263, 303, 343,&
 & 344, 322, 292, 310, 409, 521, 599, 614, 617, 624, 595, 593/)  
YDEAERC_TEGEN%RAERV(: ,1993,11) = (/&
 & 467, 467, 470, 422, 402, 381, 325, 260, 230, 237, 258, 301,&
 & 342, 319, 284, 295, 373, 450, 537, 538, 532, 535, 539, 539/)  
YDEAERC_TEGEN%RAERV(: ,1993,12) = (/&
 & 446, 446, 452, 405, 395, 382, 342, 255, 212, 215, 239, 265,&
 & 273, 261, 236, 232, 263, 349, 455, 503, 537, 537, 512, 510/)  
YDEAERC_TEGEN%RAERV(: ,1994, 1) = (/&
 & 422, 422, 431, 384, 384, 384, 341, 254, 193, 190, 216, 242,&
 & 260, 250, 219, 206, 223, 287, 376, 435, 450, 461, 454, 453/)  
YDEAERC_TEGEN%RAERV(: ,1994, 2) = (/&
 & 392, 392, 405, 375, 369, 368, 328, 241, 167, 172, 208, 220,&
 & 248, 234, 203, 188, 204, 252, 325, 390, 422, 454, 470, 471/)  
YDEAERC_TEGEN%RAERV(: ,1994, 3) = (/&
 & 360, 361, 379, 382, 363, 323, 267, 202, 153, 151, 175, 204,&
 & 215, 206, 182, 166, 180, 226, 286, 343, 381, 449, 475, 476/)  
YDEAERC_TEGEN%RAERV(: ,1994, 4) = (/&
 & 330, 330, 341, 350, 332, 303, 258, 195, 149, 145, 160, 178,&
 & 180, 175, 157, 150, 173, 219, 267, 307, 338, 419, 447, 448/)  
YDEAERC_TEGEN%RAERV(: ,1994, 5) = (/&
 & 323, 322, 316, 299, 264, 244, 211, 167, 135, 138, 155, 168,&
 & 170, 174, 158, 156, 195, 249, 270, 283, 323, 386, 414, 414/)  
YDEAERC_TEGEN%RAERV(: ,1994, 6) = (/&
 & 270, 269, 257, 253, 227, 195, 155, 128, 122, 137, 153, 145,&
 & 145, 159, 144, 142, 186, 253, 280, 281, 299, 354, 383, 383/)  
YDEAERC_TEGEN%RAERV(: ,1994, 7) = (/&
 & 205, 205, 208, 206, 196, 168, 126, 109, 116, 137, 159, 150,&
 & 146, 171, 148, 150, 196, 265, 291, 280, 275, 319, 349, 349/)  
YDEAERC_TEGEN%RAERV(: ,1994, 8) = (/&
 & 211, 211, 205, 195, 184, 162, 125, 110, 119, 141, 184, 154,&
 & 142, 151, 142, 136, 179, 240, 267, 266, 254, 282, 313, 313/)  
YDEAERC_TEGEN%RAERV(: ,1994, 9) = (/&
 & 218, 218, 209, 186, 175, 149, 124, 107, 113, 135, 159, 137,&
 & 133, 131, 125, 131, 174, 223, 250, 249, 239, 256, 273, 273/)  
YDEAERC_TEGEN%RAERV(: ,1994,10) = (/&
 & 250, 250, 240, 182, 173, 159, 129, 111, 109, 125, 154, 154,&
 & 154, 146, 125, 134, 185, 239, 256, 246, 238, 238, 231, 231/)  
YDEAERC_TEGEN%RAERV(: ,1994,11) = (/&
 & 294, 293, 285, 198, 188, 189, 160, 122, 106, 110, 123, 134,&
 & 149, 140, 118, 118, 152, 198, 232, 237, 238, 228, 219, 219/)  
YDEAERC_TEGEN%RAERV(: ,1994,12) = (/&
 & 367, 367, 363, 309, 261, 201, 177, 126,  95,  90, 100, 124,&
 & 148, 132, 114, 107, 127, 164, 211, 233, 241, 246, 267, 268/)  
YDEAERC_TEGEN%RAERV(: ,1995, 1) = (/&
 & 338, 338, 334, 284, 240, 185, 162, 116,  87,  83,  92, 114,&
 & 136, 121, 105,  99, 117, 151, 194, 214, 222, 226, 245, 247/)  
YDEAERC_TEGEN%RAERV(: ,1995, 2) = (/&
 & 311, 311, 308, 262, 221, 170, 149, 107,  80,  76,  85, 105,&
 & 125, 111,  96,  91, 108, 139, 179, 197, 204, 208, 226, 227/)  
YDEAERC_TEGEN%RAERV(: ,1995, 3) = (/&
 & 286, 286, 283, 241, 204, 156, 138,  98,  74,  70,  78,  96,&
 & 115, 102,  89,  83,  99, 128, 165, 181, 188, 192, 208, 209/)  
YDEAERC_TEGEN%RAERV(: ,1995, 4) = (/&
 & 263, 263, 260, 222, 187, 144, 127,  90,  68,  64,  72,  89,&
 & 106,  94,  81,  77,  91, 118, 151, 167, 173, 176, 191, 192/)  
YDEAERC_TEGEN%RAERV(: ,1995, 5) = (/&
 & 242, 242, 240, 204, 172, 132, 116,  83,  63,  59,  66,  82,&
 & 97,  87,  75,  71,  84, 108, 139, 153, 159, 162, 176, 177/)  
YDEAERC_TEGEN%RAERV(: ,1995, 6) = (/&
 & 223, 223, 220, 187, 158, 122, 107,  76,  58,  54,  61,  75,&
 & 90,  80,  69,  65,  77, 100, 128, 141, 146, 149, 162, 163/)  
YDEAERC_TEGEN%RAERV(: ,1995, 7) = (/&
 & 205, 205, 203, 173, 146, 112,  99,  70,  53,  50,  56,  69,&
 & 82,  73,  63,  60,  71,  92, 118, 130, 135, 137, 149, 150/)  
YDEAERC_TEGEN%RAERV(: ,1995, 8) = (/&
 & 189, 189, 187, 159, 134, 103,  91,  65,  49,  46,  51,  64,&
 & 76,  68,  58,  55,  65,  84, 108, 120, 124, 126, 137, 138/)  
YDEAERC_TEGEN%RAERV(: ,1995, 9) = (/&
 & 174, 174, 172, 146, 123,  95,  83,  59,  45,  42,  47,  58,&
 & 70,  62,  54,  51,  60,  78, 100, 110, 114, 116, 126, 127/)  
YDEAERC_TEGEN%RAERV(: ,1995,10) = (/&
 & 160, 160, 158, 134, 114,  87,  77,  55,  41,  39,  43,  54,&
 & 64,  57,  49,  47,  55,  71,  92, 101, 105, 107, 116, 117/)  
YDEAERC_TEGEN%RAERV(: ,1995,11) = (/&
 & 147, 147, 145, 124, 104,  80,  71,  50,  38,  36,  40,  50,&
 & 59,  53,  45,  43,  51,  66,  84,  93,  96,  98, 107, 107/)  
YDEAERC_TEGEN%RAERV(: ,1995,12) = (/&
 & 135, 135, 134, 114,  96,  74,  65,  46,  35,  33,  37,  46,&
 & 54,  48,  42,  39,  47,  61,  78,  86,  89,  91,  98,  99/)  
YDEAERC_TEGEN%RAERV(: ,1996, 1) = (/&
 & 124, 124, 123, 105,  88,  68,  60,  43,  32,  30,  34,  42,&
 & 50,  45,  38,  36,  43,  56,  72,  79,  82,  83,  90,  91/)  
YDEAERC_TEGEN%RAERV(: ,1996, 2) = (/&
 & 114, 114, 113,  96,  81,  62,  55,  39,  30,  28,  31,  39,&
 & 46,  41,  35,  33,  40,  51,  66,  72,  75,  77,  83,  83/)  
YDEAERC_TEGEN%RAERV(: ,1996, 3) = (/&
 & 105, 105, 104,  89,  75,  57,  51,  36,  27,  26,  29,  35,&
 & 42,  38,  33,  31,  36,  47,  61,  67,  69,  70,  76,  77/)  
YDEAERC_TEGEN%RAERV(: ,1996, 4) = (/&
 & 97,  97,  96,  81,  69,  53,  47,  33,  25,  24,  26,  33,&
 & 39,  35,  30,  28,  33,  43,  56,  61,  64,  65,  70,  71/)  
YDEAERC_TEGEN%RAERV(: ,1996, 5) = (/&
 & 89,  89,  88,  75,  63,  49,  43,  31,  23,  22,  24,  30,&
 & 36,  32,  28,  26,  31,  40,  51,  56,  58,  60,  65,  65/)  
YDEAERC_TEGEN%RAERV(: ,1996, 6) = (/&
 & 82,  82,  81,  69,  58,  45,  39,  28,  21,  20,  22,  28,&
 & 33,  29,  25,  24,  28,  37,  47,  52,  54,  55,  59,  60/)  
YDEAERC_TEGEN%RAERV(: ,1996, 7) = (/&
 & 75,  75,  75,  63,  54,  41,  36,  26,  19,  18,  21,  25,&
 & 30,  27,  23,  22,  26,  34,  43,  48,  50,  51,  55,  55/)  
YDEAERC_TEGEN%RAERV(: ,1996, 8) = (/&
 & 69,  69,  69,  58,  49,  38,  33,  24,  18,  17,  19,  23,&
 & 28,  25,  21,  20,  24,  31,  40,  44,  46,  46,  50,  51/)  

IF (LHOOK) CALL DR_HOOK('YOEAERC:PART5',1,ZHOOK_HANDLE)
END SUBROUTINE PART5


SUBROUTINE PART6

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOEAERC:PART6',0,ZHOOK_HANDLE)

YDEAERC_TEGEN%RAERV(: ,1996, 9) = (/&
 & 64,  64,  63,  54,  45,  35,  31,  22,  17,  16,  17,  22,&
 & 26,  23,  20,  19,  22,  29,  37,  40,  42,  43,  46,  47/)  
YDEAERC_TEGEN%RAERV(: ,1996,10) = (/&
 & 59,  59,  58,  49,  42,  32,  28,  20,  15,  14,  16,  20,&
 & 24,  21,  18,  17,  20,  26,  34,  37,  39,  39,  43,  43/)  
YDEAERC_TEGEN%RAERV(: ,1996,11) = (/&
 & 54,  54,  53,  45,  38,  30,  26,  19,  14,  13,  15,  18,&
 & 22,  19,  17,  16,  19,  24,  31,  34,  35,  36,  39,  39/)  
YDEAERC_TEGEN%RAERV(: ,1996,12) = (/&
 & 50,  50,  49,  42,  35,  27,  24,  17,  13,  12,  14,  17,&
 & 20,  18,  15,  14,  17,  22,  29,  32,  33,  33,  36,  36/)  
YDEAERC_TEGEN%RAERV(: ,1997, 1) = (/&
 & 46,  46,  45,  38,  33,  25,  22,  16,  12,  11,  12,  15,&
 & 18,  16,  14,  13,  16,  20,  26,  29,  30,  31,  33,  33/)  
YDEAERC_TEGEN%RAERV(: ,1997, 2) = (/&
 & 42,  42,  42,  35,  30,  23,  20,  14,  11,  10,  11,  14,&
 & 17,  15,  13,  12,  15,  19,  24,  27,  28,  28,  31,  31/)  
YDEAERC_TEGEN%RAERV(: ,1997, 3) = (/&
 & 39,  39,  38,  33,  28,  21,  19,  13,  10,   9,  11,  13,&
 & 16,  14,  12,  11,  13,  17,  22,  25,  25,  26,  28,  28/)  
YDEAERC_TEGEN%RAERV(: ,1997, 4) = (/&
 & 36,  36,  35,  30,  25,  19,  17,  12,   9,   9,  10,  12,&
 & 14,  13,  11,  10,  12,  16,  20,  23,  23,  24,  26,  26/)  
YDEAERC_TEGEN%RAERV(: ,1997, 5) = (/&
 & 33,  33,  32,  28,  23,  18,  16,  11,   8,   8,   9,  11,&
 & 13,  12,  10,  10,  11,  15,  19,  21,  22,  22,  24,  24/)  
YDEAERC_TEGEN%RAERV(: ,1997, 6) = (/&
 & 30,  30,  30,  25,  21,  16,  14,  10,   8,   7,   8,  10,&
 & 12,  11,   9,   9,  10,  13,  17,  19,  20,  20,  22,  22/)  
YDEAERC_TEGEN%RAERV(: ,1997, 7) = (/&
 & 28,  28,  27,  23,  20,  15,  13,  10,   7,   7,   8,   9,&
 & 11,  10,   9,   8,  10,  12,  16,  18,  18,  19,  20,  20/)  
YDEAERC_TEGEN%RAERV(: ,1997, 8) = (/&
 & 26,  26,  25,  21,  18,  14,  12,   9,   7,   6,   7,   9,&
 & 10,   9,   8,   7,   9,  11,  15,  16,  17,  17,  19,  19/)  
YDEAERC_TEGEN%RAERV(: ,1997, 9) = (/&
 & 23,  23,  23,  20,  17,  13,  11,   8,   6,   6,   6,   8,&
 & 9,   8,   7,   7,   8,  11,  14,  15,  15,  16,  17,  17/)  
YDEAERC_TEGEN%RAERV(: ,1997,10) = (/&
 & 22,  22,  21,  18,  15,  12,  10,   7,   6,   5,   6,   7,&
 & 9,   8,   7,   6,   7,  10,  12,  14,  14,  14,  16,  16/)  
YDEAERC_TEGEN%RAERV(: ,1997,11) = (/&
 & 20,  20,  20,  17,  14,  11,  10,   7,   5,   5,   5,   7,&
 & 8,   7,   6,   6,   7,   9,  11,  13,  13,  13,  14,  15/)  
YDEAERC_TEGEN%RAERV(: ,1997,12) = (/&
 & 18,  18,  18,  15,  13,  10,   9,   6,   5,   4,   5,   6,&
 & 7,   7,   6,   5,   6,   8,  11,  12,  12,  12,  13,  13/)  
YDEAERC_TEGEN%RAERV(: ,1998, 1) = (/&
 & 17,  17,  17,  14,  12,   9,   8,   6,   4,   4,   5,   6,&
 & 7,   6,   5,   5,   6,   8,  10,  11,  11,  11,  12,  12/)  
YDEAERC_TEGEN%RAERV(: ,1998, 2) = (/&
 & 15,  15,  15,  13,  11,   8,   7,   5,   4,   4,   4,   5,&
 & 6,   6,   5,   5,   5,   7,   9,  10,  10,  10,  11,  11/)  
YDEAERC_TEGEN%RAERV(: ,1998, 3) = (/&
 & 14,  14,  14,  12,  10,   8,   7,   5,   4,   3,   4,   5,&
 & 6,   5,   4,   4,   5,   6,   8,   9,   9,  10,  10,  10/)  
YDEAERC_TEGEN%RAERV(: ,1998, 4) = (/&
 & 13,  13,  13,  11,   9,   7,   6,   4,   3,   3,   4,   4,&
 & 5,   5,   4,   4,   5,   6,   8,   8,   9,   9,  10,  10/)  
YDEAERC_TEGEN%RAERV(: ,1998, 5) = (/&
 & 12,  12,  12,  10,   9,   7,   6,   4,   3,   3,   3,   4,&
 & 5,   4,   4,   4,   4,   5,   7,   8,   8,   8,   9,   9/)  
YDEAERC_TEGEN%RAERV(: ,1998, 6) = (/&
 & 11,  11,  11,   9,   8,   6,   5,   4,   3,   3,   3,   4,&
 & 4,   4,   3,   3,   4,   5,   6,   7,   7,   7,   8,   8/)  
YDEAERC_TEGEN%RAERV(: ,1998, 7) = (/&
 & 10,  10,  10,   9,   7,   6,   5,   3,   3,   2,   3,   3,&
 & 4,   4,   3,   3,   4,   5,   6,   6,   7,   7,   7,   7/)  
YDEAERC_TEGEN%RAERV(: ,1998, 8) = (/&
 & 9,   9,   9,   8,   7,   5,   5,   3,   2,   2,   3,   3,&
 & 4,   3,   3,   3,   3,   4,   5,   6,   6,   6,   7,   7/)  
YDEAERC_TEGEN%RAERV(: ,1998, 9) = (/&
 & 9,   9,   9,   7,   6,   5,   4,   3,   2,   2,   2,   3,&
 & 3,   3,   3,   3,   3,   4,   5,   5,   6,   6,   6,   6/)  
YDEAERC_TEGEN%RAERV(: ,1998,10) = (/&
 & 8,   8,   8,   7,   6,   4,   4,   3,   2,   2,   2,   3,&
 & 3,   3,   2,   2,   3,   4,   5,   5,   5,   5,   6,   6/)  
YDEAERC_TEGEN%RAERV(: ,1998,11) = (/&
 & 7,   7,   7,   6,   5,   4,   4,   3,   2,   2,   2,   2,&
 & 3,   3,   2,   2,   3,   3,   4,   5,   5,   5,   5,   5/)  
YDEAERC_TEGEN%RAERV(: ,1998,12) = (/&
 & 7,   7,   7,   6,   5,   4,   3,   2,   2,   2,   2,   2,&
 & 3,   2,   2,   2,   2,   3,   4,   4,   4,   5,   5,   5/)  
YDEAERC_TEGEN%RAERV(: ,1999, 1) = (/&
 & 6,   6,   6,   5,   4,   3,   3,   2,   2,   2,   2,   2,&
 & 2,   2,   2,   2,   2,   3,   4,   4,   4,   4,   4,   5/)  
YDEAERC_TEGEN%RAERV(: ,1999, 2) = (/&
 & 6,   6,   6,   5,   4,   3,   3,   2,   1,   1,   2,   2,&
 & 2,   2,   2,   2,   2,   3,   3,   4,   4,   4,   4,   4/)  
YDEAERC_TEGEN%RAERV(: ,1999, 3) = (/&
 & 5,   5,   5,   4,   4,   3,   3,   2,   1,   1,   1,   2,&
 & 2,   2,   2,   2,   2,   2,   3,   3,   3,   4,   4,   4/)  
YDEAERC_TEGEN%RAERV(: ,1999, 4) = (/&
 & 5,   5,   5,   4,   3,   3,   2,   2,   1,   1,   1,   2,&
 & 2,   2,   1,   1,   2,   2,   3,   3,   3,   3,   3,   4/)  
YDEAERC_TEGEN%RAERV(: ,1999, 5) = (/&
 & 4,   4,   4,   4,   3,   2,   2,   2,   1,   1,   1,   1,&
 & 2,   2,   1,   1,   2,   2,   3,   3,   3,   3,   3,   3/)  
YDEAERC_TEGEN%RAERV(: ,1999, 6) = (/&
 & 4,   4,   4,   3,   3,   2,   2,   1,   1,   1,   1,   1,&
 & 2,   1,   1,   1,   1,   2,   2,   3,   3,   3,   3,   3/)  
YDEAERC_TEGEN%RAERV(: ,1999, 7) = (/&
 & 4,   4,   4,   3,   3,   2,   2,   1,   1,   1,   1,   1,&
 & 2,   1,   1,   1,   1,   2,   2,   2,   2,   3,   3,   3/)  
YDEAERC_TEGEN%RAERV(: ,1999, 8) = (/&
 & 3,   3,   3,   3,   2,   2,   2,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3/)  
YDEAERC_TEGEN%RAERV(: ,1999, 9) = (/&
 & 3,   3,   3,   3,   2,   2,   2,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   2/)  
YDEAERC_TEGEN%RAERV(: ,1999,10) = (/&
 & 3,   3,   3,   2,   2,   2,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   2/)  
YDEAERC_TEGEN%RAERV(: ,1999,11) = (/&
 & 3,   3,   3,   2,   2,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2,   2/)  
YDEAERC_TEGEN%RAERV(: ,1999,12) = (/&
 & 2,   2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   2,   2,   2,   2,   2/)  
YDEAERC_TEGEN%RAERV(: ,2000, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2000, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2000,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2001, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2001,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2002, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2002, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2002,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2003, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2003, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2003,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2004, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2004, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

IF (LHOOK) CALL DR_HOOK('YOEAERC:PART6',1,ZHOOK_HANDLE)
END SUBROUTINE PART6


SUBROUTINE PART7

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOEAERC:PART7',0,ZHOOK_HANDLE)

YDEAERC_TEGEN%RAERV(: ,2004, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2004,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2005, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2005, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2005,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2006, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2006, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2006,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2007, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2007, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2007,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2008, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2008, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2008,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2009, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2009, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2009,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

YDEAERC_TEGEN%RAERV(: ,2010, 1) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   2/)  
YDEAERC_TEGEN%RAERV(: ,2010, 2) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010, 3) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010, 4) = (/&
 & 2,   2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010, 5) = (/&
 & 2,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010, 6) = (/&
 & 2,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010, 7) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010, 8) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010, 9) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010,10) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010,11) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  
YDEAERC_TEGEN%RAERV(: ,2010,12) = (/&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,&
 & 1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/)  

IF (LHOOK) CALL DR_HOOK('YOEAERC:PART7',1,ZHOOK_HANDLE)
END SUBROUTINE PART7

END SUBROUTINE SUECAEC


SUBROUTINE SU_AER_CLIMATOLOGY(SELF_MACC)

! PURPOSE
! -------
! Load vertically integrated mass of aerosol species in the form of
! monthly averages with 3x3 degree lon-lat resolution, from CAMS
! free-running control run (experiment gbst), bias corrected with
! total aerosol optical depth CAMS-Interim reanalysis 2003-2011
! (experiment eac4).
!
! INTERFACE
! ---------
! SU_AER_CLIMATOLOGY is called from TAERC_MACC%SETUP()
!
!
! AUTHORS
! -------
! R. Hogan and A. Bozzo, ECMWF 2016-09-28
!
! MODIFICATIONS
! -------------
!   2019-02-02  R. Hogan  Open/close file once, rather than for each variable
!
!-----------------------------------------------------------------------
USE PARKIND1  , ONLY : JPIM, JPRB
USE EASY_NETCDF_READ_MPI,ONLY : NETCDF_FILE
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULERR

IMPLICIT NONE

! ARGUMENTS
! ------------

TYPE(TEAERC_MACC), INTENT(INOUT) :: SELF_MACC



! Expected dimension lengths for month, latitude and longitude
INTEGER(KIND=JPIM), PARAMETER :: IMONTH = 12
INTEGER(KIND=JPIM), PARAMETER :: ILAT   = 61
INTEGER(KIND=JPIM), PARAMETER :: ILON   = 120

! Full and base names of the aerosol climatology file
CHARACTER(len=:),ALLOCATABLE :: CL_FILE_NAME, CL_AER_CLIM_FILE
! "DATA" directory
CHARACTER(LEN=512)  :: CLDIRECTORY

! The NetCDF file containing the input climatology
TYPE(NETCDF_FILE)  :: FILE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_AER_CLIMATOLOGY',0,ZHOOK_HANDLE)

  CL_AER_CLIM_FILE = "aerosol_cams_climatology_43R3a.nc"
  CALL GET_ENVIRONMENT_VARIABLE("DATA",CLDIRECTORY)
  IF(CLDIRECTORY /= " ") THEN
    CL_FILE_NAME = TRIM(CLDIRECTORY) // "/ifsdata/" // CL_AER_CLIM_FILE
  ELSE
    CL_FILE_NAME = CL_AER_CLIM_FILE
  ENDIF

  CALL FILE%OPEN(TRIM(CL_FILE_NAME), IVERBOSE=4)

  CALL FILE%GET('Black_Carbon_hydrophilic', SELF_MACC%RMACBC1)

  ! Check dimensions
  IF (    SIZE(SELF_MACC%RMACBC1,1) /= ILON .OR. &
       &  SIZE(SELF_MACC%RMACBC1,2) /= ILAT .OR. &
       &  SIZE(SELF_MACC%RMACBC1,3) /= IMONTH) THEN
    WRITE(NULERR,*) 'Aerosols in', CL_FILE_NAME, 'must be dimensioned (', &
         &  ILON, ",", ILAT, ",", IMONTH, ")"
    CALL ABOR1("Error reading aerosol climatology NetCDF file")
  ENDIF

  CALL FILE%GET('Black_Carbon_hydrophobic', SELF_MACC%RMACBC2)
  CALL FILE%GET('Organic_Matter_hydrophilic', SELF_MACC%RMACOR1)
  CALL FILE%GET('Organic_Matter_hydrophobic', SELF_MACC%RMACOR2)
  CALL FILE%GET('Mineral_Dust_bin1', SELF_MACC%RMACSD1)
  CALL FILE%GET('Mineral_Dust_bin2', SELF_MACC%RMACSD2)
  CALL FILE%GET('Mineral_Dust_bin3', SELF_MACC%RMACSD3)
  CALL FILE%GET('Sea_Salt_bin1', SELF_MACC%RMACSS1)
  CALL FILE%GET('Sea_Salt_bin2', SELF_MACC%RMACSS2)
  CALL FILE%GET('Sea_Salt_bin3', SELF_MACC%RMACSS3)
  CALL FILE%GET('Sulfates', SELF_MACC%RMACSU1)

  CALL FILE%CLOSE()

  ! Allocate the interpolation arrays
  IF (ANY( (/ ALLOCATED(SELF_MACC%RMAEBC1), ALLOCATED(SELF_MACC%RMAEBC2), ALLOCATED(SELF_MACC%RMAEOR1), ALLOCATED(SELF_MACC%RMAEOR2), &
 &            ALLOCATED(SELF_MACC%RMAESD1), ALLOCATED(SELF_MACC%RMAESD2), ALLOCATED(SELF_MACC%RMAESD3), &
 &            ALLOCATED(SELF_MACC%RMAESS1), ALLOCATED(SELF_MACC%RMAESS2), ALLOCATED(SELF_MACC%RMAESS3), ALLOCATED(SELF_MACC%RMAESU1)/) )) THEN
    WRITE(NULERR,*) "Aerosol climatology interpolation arrays already allocated, are we expecting this? "
  ELSE
    ALLOCATE(SELF_MACC%RMAEBC1(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAEBC2(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAEOR1(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAEOR2(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAESD1(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAESD2(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAESD3(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAESS1(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAESS2(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAESS3(ILON,ILAT))
    ALLOCATE(SELF_MACC%RMAESU1(ILON,ILAT))
  ENDIF


IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_AER_CLIMATOLOGY',1,ZHOOK_HANDLE)

END SUBROUTINE SU_AER_CLIMATOLOGY

SUBROUTINE SU_AER_CLIMATOLOGY3D(SELF_MACC)

! PURPOSE
! -------
! Load 3D climatologies of mass of aerosol species in the form of
! monthly averages with 3x3 degree lon-lat resolution, 60 levels from CAMS
! free-running control run (experiment gbst), bias corrected with
! total aerosol optical depth CAMS-Interim reanalysis 2003-2011
! (experiment eac4).
!
! INTERFACE
! ---------
! SU_AER_CLIMATOLOGY3D is called from TAERC_MACC%SETUP()
!
!
! AUTHORS
! -------
! A. Bozzo, ECMWF 2018-01-22
!
! MODIFICATIONS
! -------------
!   2019-02-02  R. Hogan   Open/close file once, rather than for each variable
!   2021-03-15  O. Marsden Allow climatology arrays to be read in simultaneously by multiple MPI ranks. 
!
!-----------------------------------------------------------------------
USE PARKIND1            , ONLY : JPIM, JPRB
USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN              , ONLY : NULERR, NULOUT
USE EASY_NETCDF_READ_MPI, ONLY : NETCDF_FILE
USE YOMMP0              , ONLY : MYPROC, NPROC
USE MPL_MODULE          , ONLY : MPL_WAIT, JP_BLOCKING_STANDARD, JP_BLOCKING_BUFFERED, &
 &                               JP_NON_BLOCKING_STANDARD, JP_NON_BLOCKING_BUFFERED
USE OMP_LIB

IMPLICIT NONE

! ARGUMENTS
! ------------

TYPE(TEAERC_MACC), INTENT(INOUT) :: SELF_MACC

! LOCAL VARIABLES
! ---------------

! Expected dimension lengths for month, latitude and longitude
INTEGER(KIND=JPIM), PARAMETER :: IMONTH = 12
INTEGER(KIND=JPIM), PARAMETER :: ILAT   = 61
INTEGER(KIND=JPIM), PARAMETER :: ILON   = 120
INTEGER(KIND=JPIM), PARAMETER :: ILEV   = 60

! Full and base names of the aerosol climatology file
CHARACTER(LEN=512) :: CL_FILE_NAME, CL_AER_CLIM_FILE
! "DATA" directory
CHARACTER(LEN=512)  :: CLDIRECTORY

! The NetCDF file containing the input climatology
TYPE(NETCDF_FILE)  :: FILE
LOGICAL :: LCHECK_NC_ARRAYS = .FALSE.

! NetCDF read timings
REAL(KIND=JPRB) :: T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15
LOGICAL :: LPRINT_TIMINGS = .FALSE.

!! set LACTIVE to true if multiple tasks should cooperate on reading in aerosol climatology arrays
LOGICAL :: LACTIVE = .TRUE.
!!INTEGER(KIND=JPIM)   :: COMM_TYPE = JP_BLOCKING_STANDARD 
!!INTEGER(KIND=JPIM)   :: COMM_TYPE = JP_NON_BLOCKING_BUFFERED 
INTEGER(KIND=JPIM)   :: COMM_TYPE = JP_NON_BLOCKING_STANDARD 
INTEGER(KIND=JPIM)   :: II,JJ,KK
!! list of MPI tasks that should cooperate on input if LACTIVE=.TRUE.
INTEGER(KIND=JPIM)   :: IACTIVE_PROCS(13), IPROC_SKIP
! For non-blocking MPI mode, the list of comm requests generated
INTEGER(KIND=JPIM)   :: IREQ_LIST(13)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_AER_CLIMATOLOGY3D',0,ZHOOK_HANDLE)

!! check whether the main storage arrays have already been initialised
IF (SMART_MACC3D_POINTER_COPIES > 0) THEN
  SMART_MACC3D_POINTER_COPIES = SMART_MACC3D_POINTER_COPIES + 1 !! this is used to be able to deallocate cleanly only if number reaches 0
ELSE

  CL_AER_CLIM_FILE = "aerosol_cams_climatology_43R3_3D.nc"
  CALL GET_ENVIRONMENT_VARIABLE("DATA",CLDIRECTORY)
  IF(CLDIRECTORY /= " ") THEN
    CL_FILE_NAME = TRIM(CLDIRECTORY) // "/ifsdata/" // CL_AER_CLIM_FILE
  ELSE
    CL_FILE_NAME = CL_AER_CLIM_FILE
  ENDIF
  
  IF (.NOT. LACTIVE) then 
    CALL FILE%OPEN(TRIM(CL_FILE_NAME), IVERBOSE=4)
  ELSE
    CALL FILE%OPEN_ACTIVE(TRIM(CL_FILE_NAME), IVERBOSE=4)
    IPROC_SKIP = NPROC / 13
    IACTIVE_PROCS = [ (1+II*IPROC_SKIP,II=0,12) ]
  ENDIF
  
  
  IF(LPRINT_TIMINGS) T1 = OMP_GET_WTIME()
  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Black_Carbon_hydrophilic', RMAC3D_BC1)
  ELSE
    CALL FILE%GET_ACTIVE('Black_Carbon_hydrophilic', RMAC3D_BC1, IACTIVE_RANK=IACTIVE_PROCS(1), IREQUEST=IREQ_LIST(1), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_BC1, 'BC1')
  ! Check dimensions
  IF (    SIZE(RMAC3D_BC1,1) /= ILON .OR. &
       &  SIZE(RMAC3D_BC1,2) /= ILAT .OR. &
       &  SIZE(RMAC3D_BC1,3) /= ILEV .OR. &
       &  SIZE(RMAC3D_BC1,4) /= IMONTH) THEN
    WRITE(NULERR,*) 'Aerosols in', CL_FILE_NAME, 'must be dimensioned (', &
         &  ILON, ",", ILAT, ",", ILEV, ",", IMONTH, ")"
    CALL ABOR1("Error reading 3D aerosol climatology NetCDF file")
  ENDIF
  IF(LPRINT_TIMINGS) T2 = OMP_GET_WTIME()
    
  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Black_Carbon_hydrophobic', RMAC3D_BC2)
  ELSE
    CALL FILE%GET_ACTIVE('Black_Carbon_hydrophobic', RMAC3D_BC2, IACTIVE_RANK=IACTIVE_PROCS(2), IREQUEST=IREQ_LIST(2), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_BC2 , 'BC2')
  IF(LPRINT_TIMINGS) T3 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Organic_Matter_hydrophilic', RMAC3D_OR1)
  ELSE
    CALL FILE%GET_ACTIVE('Organic_Matter_hydrophilic', RMAC3D_OR1, IACTIVE_RANK=IACTIVE_PROCS(3), IREQUEST=IREQ_LIST(3), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_OR1 , 'OR1')
  IF(LPRINT_TIMINGS) T4 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Organic_Matter_hydrophobic', RMAC3D_OR2)
  ELSE
    CALL FILE%GET_ACTIVE('Organic_Matter_hydrophobic', RMAC3D_OR2, IACTIVE_RANK=IACTIVE_PROCS(4), IREQUEST=IREQ_LIST(4), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_OR2 , 'OR2')
  IF(LPRINT_TIMINGS) T5 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Mineral_Dust_bin1', RMAC3D_SD1)
  ELSE
    CALL FILE%GET_ACTIVE('Mineral_Dust_bin1', RMAC3D_SD1, IACTIVE_RANK=IACTIVE_PROCS(5), IREQUEST=IREQ_LIST(5), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_SD1 , 'SD1')
  IF(LPRINT_TIMINGS) T6 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Mineral_Dust_bin2', RMAC3D_SD2)
  ELSE
    CALL FILE%GET_ACTIVE('Mineral_Dust_bin2', RMAC3D_SD2, IACTIVE_RANK=IACTIVE_PROCS(6), IREQUEST=IREQ_LIST(6), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_SD2 , 'SD2')
  IF(LPRINT_TIMINGS) T7 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Mineral_Dust_bin3', RMAC3D_SD3)
  ELSE
    CALL FILE%GET_ACTIVE('Mineral_Dust_bin3', RMAC3D_SD3, IACTIVE_RANK=IACTIVE_PROCS(7), IREQUEST=IREQ_LIST(7), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_SD3, 'SD3')
  IF(LPRINT_TIMINGS) T8 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Sea_Salt_bin1', RMAC3D_SS1)
  ELSE
    CALL FILE%GET_ACTIVE('Sea_Salt_bin1', RMAC3D_SS1, IACTIVE_RANK=IACTIVE_PROCS(8), IREQUEST=IREQ_LIST(8), IMP_TYPE=COMM_TYPE)
  endif
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_SS1 , 'SS1')
  IF(LPRINT_TIMINGS) T9 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Sea_Salt_bin2', RMAC3D_SS2)
  ELSE
    CALL FILE%GET_ACTIVE('Sea_Salt_bin2', RMAC3D_SS2, IACTIVE_RANK=IACTIVE_PROCS(9), IREQUEST=IREQ_LIST(9), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_SS2 , 'SS2')
  IF(LPRINT_TIMINGS) T10 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Sea_Salt_bin3', RMAC3D_SS3)
  ELSE
    CALL FILE%GET_ACTIVE('Sea_Salt_bin3', RMAC3D_SS3, IACTIVE_RANK=IACTIVE_PROCS(10), IREQUEST=IREQ_LIST(10), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_SS3,  'SS3')
  IF(LPRINT_TIMINGS) T11 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('Sulfates', RMAC3D_SU1)
  ELSE
    CALL FILE%GET_ACTIVE('Sulfates', RMAC3D_SU1, IACTIVE_RANK=IACTIVE_PROCS(11), IREQUEST=IREQ_LIST(11), IMP_TYPE=COMM_TYPE)
  ENDIF
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_SU1, 'SU1')
  IF(LPRINT_TIMINGS) T12 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('half_level_pressure', RMAC3D_PRS)
  ELSE
    CALL FILE%GET_ACTIVE('half_level_pressure', RMAC3D_PRS, IACTIVE_RANK=IACTIVE_PROCS(12), IREQUEST=IREQ_LIST(12), IMP_TYPE=COMM_TYPE)
  endif
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_PRS, 'PRS')
  IF(LPRINT_TIMINGS) T13 = OMP_GET_WTIME()

  IF (.NOT. LACTIVE) THEN 
  CALL FILE%GET('half_level_delta_pressure', RMAC3D_DPRS)
  ELSE
    CALL FILE%GET_ACTIVE('half_level_delta_pressure', RMAC3D_DPRS, IACTIVE_RANK=IACTIVE_PROCS(13), IREQUEST=IREQ_LIST(13), IMP_TYPE=COMM_TYPE)
  endif
  IF (LCHECK_NC_ARRAYS) CALL CHECK_ARRAY(RMAC3D_DPRS , 'DPRS')
  IF(LPRINT_TIMINGS) T14 = OMP_GET_WTIME()

  IF (LACTIVE .AND. (COMM_TYPE .EQ. JP_NON_BLOCKING_STANDARD .OR. COMM_TYPE .EQ. JP_NON_BLOCKING_BUFFERED)) THEN
    CALL MPL_WAIT( IREQ_LIST(1:13) )
  END IF
  IF(LPRINT_TIMINGS) T15 = OMP_GET_WTIME()
  IF (MYPROC == 1 .AND. LPRINT_TIMINGS) WRITE(NULOUT,*) 'NetCDF read times in SU_AER_CLIMATOLOGY3D: ', & 
 &        T2-T1, T3-T2, T4-T3, T5-T4, T6-T5, T7-T6, T8-T7, T9-T8, T10-T9, T11-T10, T12-T11, T13-T12, T14-T13, T15-T14, ' total time : ',T15-T1
  
  CALL FILE%CLOSE()
  
  !! this is the first time the storage arrays are being set up, so set to 1
  SMART_MACC3D_POINTER_COPIES = 1 !! this is used to be able to deallocate cleanly only if number reaches 0
  
ENDIF

! Allocate the interpolation arrays
IF (ANY( (/ ALLOCATED(SELF_MACC%RMAEBC13D), ALLOCATED(SELF_MACC%RMAEBC23D),  &
 &          ALLOCATED(SELF_MACC%RMAEOR13D), ALLOCATED(SELF_MACC%RMAEOR23D),  &
 &          ALLOCATED(SELF_MACC%RMAESD13D), ALLOCATED(SELF_MACC%RMAESD23D),  ALLOCATED(SELF_MACC%RMAESD33D), &
 &          ALLOCATED(SELF_MACC%RMAESS13D), ALLOCATED(SELF_MACC%RMAESS23D),  ALLOCATED(SELF_MACC%RMAESS33D), &
 &          ALLOCATED(SELF_MACC%RMAESU13D), ALLOCATED(SELF_MACC%REF_AER_PRS),ALLOCATED(SELF_MACC%REF_AER_DPRS)  /) ) ) THEN
  WRITE(NULERR,*) "Aerosol 3D climatology interpolation arrays already allocated, are we expecting this? "
ELSE
  ALLOCATE(SELF_MACC%RMAEBC13D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAEBC23D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAEOR13D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAEOR23D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAESD13D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAESD23D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAESD33D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAESS13D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAESS23D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAESS33D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%RMAESU13D(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%REF_AER_PRS(ILON,ILAT,ILEV))
  ALLOCATE(SELF_MACC%REF_AER_DPRS(ILON,ILAT,ILEV))

  IF (LACTIVE) THEN
    !!first touch
    !$OMP PARALLEL DO PRIVATE(II,JJ,KK) COLLAPSE(2) SCHEDULE(STATIC) 
    DO KK=1,ILEV
      DO JJ=1,ILAT
      !$OMP SIMD
        DO II=1,ILON
          SELF_MACC%RMAEBC13D(II,JJ,KK) = 0.0_JPRB        
          SELF_MACC%RMAEBC23D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAEOR13D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAEOR23D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAESD13D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAESD23D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAESD33D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAESS13D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAESS23D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAESS33D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%RMAESU13D(II,JJ,KK) = 0.0_JPRB
          SELF_MACC%REF_AER_PRS(II,JJ,KK)=0.0_JPRB
          SELF_MACC%REF_AER_DPRS(II,JJ,KK)=0.0_JPRB 
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  ENDIF

ENDIF

SELF_MACC%RMACBC13D => RMAC3D_BC1(:,:,:,:)
SELF_MACC%RMACBC23D => RMAC3D_BC2(:,:,:,:)  
SELF_MACC%RMACOR13D => RMAC3D_OR1(:,:,:,:)
SELF_MACC%RMACOR23D => RMAC3D_OR2(:,:,:,:)
SELF_MACC%RMACSD13D => RMAC3D_SD1(:,:,:,:)
SELF_MACC%RMACSD23D => RMAC3D_SD2(:,:,:,:)
SELF_MACC%RMACSD33D => RMAC3D_SD3(:,:,:,:)
SELF_MACC%RMACSS13D => RMAC3D_SS1(:,:,:,:)
SELF_MACC%RMACSS23D => RMAC3D_SS2(:,:,:,:)
SELF_MACC%RMACSS33D => RMAC3D_SS3(:,:,:,:)
SELF_MACC%RMACSU13D => RMAC3D_SU1(:,:,:,:)
SELF_MACC%REF_MON_AER_PRS  => RMAC3D_PRS(:,:,:,:)
SELF_MACC%REF_MON_AER_DPRS => RMAC3D_DPRS(:,:,:,:)


IF (LHOOK) CALL DR_HOOK('YOEAERC:SU_AER_CLIMATOLOGY3D',1,ZHOOK_HANDLE)

END SUBROUTINE SU_AER_CLIMATOLOGY3D

SUBROUTINE CHECK_ARRAY( ARR, NAME_STR )
USE YOMLUN, ONLY : NULOUT
USE, INTRINSIC :: IEEE_ARITHMETIC
  REAL(KIND=JPRB), INTENT(IN) :: ARR(:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN) :: NAME_STR

  REAL(KIND=JPRB) :: BIG_VAL = HUGE(1.0_JPRB), SMALL_VAL = TINY(1.0_JPRB)
  LOGICAL :: LPROBLEM

  LPROBLEM = .FALSE.
  
  IF ( ANY( .NOT. IEEE_IS_FINITE(ARR ) ) ) THEN
    PRINT *, 'ARRAY RMAC '//NAME_STR, ' HAS NON-FINITE ELEMENTS'
    LPROBLEM = .TRUE.
  ENDIF 
  IF ( ANY( .NOT. IEEE_IS_NORMAL(ARR) ) ) THEN
    PRINT *, 'ARRAY RMAC '//NAME_STR, ' HAS NON-NORMAL ELEMENTS'
    LPROBLEM = .TRUE.
  ENDIF 
  IF ( ANY( (ABS(ARR) .GE. BIG_VAL) ) ) THEN
    PRINT *, 'ARRAY RMAC '//NAME_STR, ' HAS OVERFLOWING ELEMENTS'
    LPROBLEM = .TRUE.
  ENDIF 
  IF ( ANY( (ABS(ARR) .LE. SMALL_VAL) .AND. ABS(ARR) .GT. 0.0_JPRB ) ) THEN
    PRINT *, 'ARRAY RMAC '//NAME_STR, ' HAS UNDERFLOWING ELEMENTS'
    LPROBLEM = .TRUE.
  ENDIF 

  IF (LPROBLEM) THEN
    WRITE(NULOUT,*) 'ARRAY RMAC '//NAME_STR,' HAS A NUMERICAL INITIALISATION ISSUE'
  ELSE
    WRITE(NULOUT,*) 'ARRAY RMAC '//NAME_STR,' IS CLEANLY INITIALISED'
  ENDIF

END SUBROUTINE CHECK_ARRAY

END MODULE YOEAERC

