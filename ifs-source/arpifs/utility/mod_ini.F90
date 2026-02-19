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

SUBROUTINE MOD_INI(KUDATE,KUSSSS,KSTEP,KUCONF,KUSTOP,&
 & KSTEPINI,PTSTEP,KULOUT)  

!**** *MOD_INI*

!     Purpose.
!     --------
!           Initialize date of model
!**   Interface.
!     ----------

!        Explicit arguments :
!        --------------------
!          KUDATE: date in model           
!          KUSSSS: time in model           
!          KUCONF: configuration            
!          KUSTOP: number of time step for the model
!          KSTEPINI: initial step in hours
!          PTSTEP: time step lenght in model 
!          KULOUT: output unit number                

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by SUARG.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Ryad El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 96-01-04
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Fisher : 2004-08-18 Increment date correctly when KSTEPINI>=48
!        M.Hamrud      01-May-2006 Generalized IO scheme
!        Carla Cardinali    13 January 2009 KUDATE not update for 801
!        M.Janiskova   11-Apr-2016 Update for TL/AD diag.run as 4D-var
!        B.Ingleby     14-Jan-2019 Remove obsolete section 3 (fsobs)
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LECMWF
USE YOETLDIAG, ONLY : LTL4DREP
!USE YOMCST   , ONLY : RDAY

IMPLICIT NONE

!     EXTERNAL INTEGER FUNCTIONS
INTEGER(KIND=JPIM),INTENT(INOUT)   :: KUDATE 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: KUSSSS 
INTEGER(KIND=JPIM),INTENT(IN)      :: KSTEP
INTEGER(KIND=JPIM),INTENT(IN)      :: KUCONF 
INTEGER(KIND=JPIM),INTENT(IN)      :: KUSTOP 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: KSTEPINI 
REAL(KIND=JPRB)   ,INTENT(IN)      :: PTSTEP 
INTEGER(KIND=JPIM),INTENT(IN)      :: KULOUT 
INTEGER(KIND=JPIM), EXTERNAL :: ICD2YMD
INTEGER(KIND=JPIM), EXTERNAL :: IYMD2CD

INTEGER(KIND=JPIM) :: ILMO(12)

INTEGER(KIND=JPIM) :: ICD, ID1, IDD, IHOUR, IHOUR1, IM1, IMM, INSSS,&
 & IRET, ISECRUN, ISHOUR, ISSS, &
 & IYMD, IYYYY, IYYYY1, INC

REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttim.func.h"

!#include "abor1.intfb.h"
#include "updcal.intfb.h"

!     ------------------------------------------------------------------

!*        1. Preparations

IF (LHOOK) CALL DR_HOOK('MOD_INI',0,ZHOOK_HANDLE)

ZEPS=1.E-10_JPRB
ISHOUR=3600

!*        2. Update the current date and time with forecast time

IF(KUCONF == 302 .OR. KUCONF == 2.OR. (KUCONF/100) == 1 &
 & .OR. ((KUCONF == 401 .OR. KUCONF == 501) .AND. LTL4DREP)) THEN
  IYYYY=NCCAA(KUDATE)
  IMM=NMM(KUDATE)
  IDD=NDD(KUDATE)
  IHOUR=NCTH(KUSSSS)
  IHOUR1=MOD(IHOUR+NCTH(KSTEP),24)
  KSTEPINI=NCTH(KSTEP)
  WRITE(KULOUT,FMT='('' GRIB IY IM ID IH:'',I6,3I3,&
   & '' NEW HOUR:'',I3)')&
   & IYYYY,IMM,IDD,IHOUR,IHOUR1  
  WRITE(KULOUT,FMT='('' KSTEPINI :'',I6)') KSTEPINI
  IF(IHOUR+ NCTH(KSTEP)>= 24) THEN
    INC=(IHOUR+NCTH(KSTEP))/24
    CALL UPDCAL(IDD,IMM,IYYYY,INC,ID1,IM1,IYYYY1,ILMO,KULOUT)
  ELSE
    IYYYY1=IYYYY
    IM1=IMM
    ID1=IDD
  ENDIF
  WRITE(KULOUT,FMT='('' NEW DATE:'',I6,2I3)') IYYYY1,IM1,ID1
  KUDATE=IYYYY1*10**4+IM1*10**2+ID1
  KUSSSS=ISHOUR*IHOUR1
ENDIF

IF (LHOOK) CALL DR_HOOK('MOD_INI',1,ZHOOK_HANDLE)
END SUBROUTINE MOD_INI
