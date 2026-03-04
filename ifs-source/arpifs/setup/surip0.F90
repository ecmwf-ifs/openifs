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

SUBROUTINE SURIP0(KULOUT, KDAT, KSSS)

!**** *SURIP0* - Routine to initialize the module YOMRIP0

!     Purpose.
!     --------
!        Initialize and print the module YOMRIP0

!**   Interface.
!     ----------
!        *CALL* *SURIP0(...)

!        Explicit arguments :
!        --------------------
!        KULOUT  - Logical unit of the output
!        KDAT    - Date in the form AAAAMMDD
!        KSSS    - Number of seconds in the day

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified : 2011-Mar A.Alias LASTRF added to prevent drift in insolation
!      Y. Bouteloup (Feb 2011) : Add RCODECF  ,RSIDECF  ,RCOVSRF  ,RSIVSRF
!      K. Yessad (July 2014) : Move model-independent variables from YOMRIP/SURIP.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMARG   , ONLY : NSUPERSEDE, NUDATE, NUSSSS
USE YOMLUN   , ONLY : NULNAM
USE YOMRIP0  , ONLY : NINDAT, NSSSSS, RTIMST, LASTRF
USE YOMCST   , ONLY : RDAY, REA, REPSM
USE YOMDYNCORE,ONLY : LAPE

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)           :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KDAT 
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KSSS 

INTEGER(KIND=JPIM) :: IA, ID, IM, IA2

REAL(KIND=JPRD) :: ZJU
REAL(KIND=JPRB) :: ZDE, ZET, ZRS, ZRSREL, ZTETA, ZTI, ZTI2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "posnam.intfb.h"

#include "namrip0.nam.h"

#include "fctast.func.h"
#include "fcttim.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURIP0',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

NINDAT=19880826
NSSSSS=43200

!     ------------------------------------------------------------------

IF (PRESENT(KDAT) .AND. PRESENT(KSSS)) THEN
  
  NINDAT=KDAT
  NSSSSS=KSSS

ELSE

!*       2.    READ NAMELIST.
!              --------------


  IF (NSUPERSEDE==0) THEN
    
    LASTRF=.FALSE.

    !*     2.1   READ NAMELIST.

    CALL POSNAM(NULNAM,'NAMRIP0')
    READ(NULNAM,NAMRIP0)

  ELSE

    !*     2.2   OVERWRITE NAMELIST WITH VALUES COMPUTED IN SUARG

    NINDAT=NUDATE
    NSSSSS=NUSSSS

  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       3     CHECKINGS AND ADDITIONAL SETUP
!              ------------------------------

ID=NDD(NINDAT)
IM=NMM(NINDAT)
IA=NCCAA(NINDAT)
ZJU=RJUDAT(IA,IM,ID)
ZTI=RTIME(IA,IM,ID,NSSSSS,RDAY)
RTIMST=ZTI
IF(LASTRF) THEN
  IA2 = MOD(IA,4)+2000
  ZTI2 = RTIME(IA2,IM,ID,NSSSSS,RDAY)
  ZTETA = RTETA(ZTI2)
ELSE
  ZTETA=RTETA(ZTI)
ENDIF
IF( LAPE ) THEN
  ZRS=RRSAQUA(ZTETA)
  ZDE=RDSAQUA(ZTETA)
  ZET=RETAQUA(ZTETA)
ELSE
  ZRS=RRS(ZTETA)
  ZDE=RDS(ZTETA)
  ZET=RET(ZTETA)
ENDIF
ZRSREL=ZRS/REA

!      ----------------------------------------------------------------

!*       4.    PRINT NAMELIST VARIABLES.
!              -------------------------

WRITE(KULOUT,'('' Printings for module YOMRIP0 '')')
WRITE(KULOUT,'('' NINDAT='',I10,'' NSSSSS='',I6)') NINDAT,NSSSSS
WRITE(KULOUT,'('' Correction of solar drift = '',L7)') LASTRF
WRITE(KULOUT,'('' The initial date of the run is : '',I4,1X,I2,1X,I2)') IA,IM,ID
WRITE(KULOUT,'('' The Julian date is : '',F11.2)') ZJU
WRITE(KULOUT,'('' Time of the model  : '',F15.2,'' s'')') ZTI
WRITE(KULOUT,'('' Distance Earth-Sun : '',E13.7,'' m'')') ZRS
WRITE(KULOUT,'('' Relative Dist. E-S : '',E13.7,'' m'')') ZRSREL
WRITE(KULOUT,'('' Declination        : '',F12.5)') ZDE
WRITE(KULOUT,'('' Eq. of time        : '',F12.5,'' s'')') ZET

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURIP0',1,ZHOOK_HANDLE)
END SUBROUTINE SURIP0
