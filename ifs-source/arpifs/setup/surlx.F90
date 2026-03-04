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

SUBROUTINE SURLX(YDDIM,YDDIMV,YDRIP,KULOUT)

!**** *SURLX* * - ROUTINE TO INITIALIZE COEFS FOR RELAXATION

!     PURPOSE.
!     --------
!        SET DEFAULT VALUES, THEN READS NAMELIST NAMRLX

!**   INTERFACE.
!     ----------
!        *CALL* *SURLX(KULOUT)*

!     EXPLICIT ARGUMENTS :  KULOUT
!     --------------------

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMRLX

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        Thomas Jung
!        ORIGINAL : 06-07-28 

!     MODIFICATIONS.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      J Flemming (Jan 2014): Relaxation for ozone
!      MO Koehler (May 2025): File path variables for relaxation files
!-----------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMRIP   , ONLY : TRIP
USE YOMRLX   , ONLY : LRLXG   ,NRLXAN   ,&
 & NFRLXG     ,NFRLXU    ,NPRLX3D   ,NPRLX2D   ,&
 & LRLXVO     ,LRLXDI    ,LRLXTE    ,LRLXQ     ,LRLXLP   ,&
 & LRLXQI     ,LRLXQL    ,LRLXQC    ,LRLXO3    ,&
 & XRLXVO     ,XRLXDI    ,XRLXTE    ,XRLXQ     ,XRLXLP   ,XRLXO3  ,&
 & ALATRLX1   ,ALATRLX2  ,ALONRLX1  ,ALONRLX2  ,&
 & NRLXLMIN   ,NRLXLMAX  ,NRLXLMINU ,NRLXLMAXU ,AXRLX ,AYRLX ,AZRLX, NRLXSMAX, &
 & CRLXPATHGG ,CRLXPATHSH
USE YOMSRLX   ,ONLY : TRLXVO   ,TRLXDI   ,TRLXTE,TRLXQ   ,TRLXLP, TRLXO3
USE YOMCST   , ONLY : RPI

!      ----------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM),         INTENT(IN) :: YDDIM
TYPE(TDIMV),        INTENT(IN) :: YDDIMV
TYPE(TRIP),         INTENT(INOUT):: YDRIP
INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT 

!      ----------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "namrlx.nam.h"

!      ----------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURLX',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NSMAX=>YDDIM%NSMAX, &
 & NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP)
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

!        1.1 Set implicit default values

LRLXG= .FALSE.

LRLXVO=.FALSE.
LRLXDI=.FALSE.
LRLXTE=.FALSE.
LRLXQ=.FALSE.
LRLXQI=.FALSE.
LRLXQL=.FALSE.
LRLXQC=.FALSE.
LRLXLP=.FALSE.
LRLXO3=.FALSE.

! default values for relaxation file paths
CRLXPATHGG='./'
CRLXPATHSH='./'

NFRLXG=2
NFRLXU=6

XRLXVO =0._JPRB
XRLXDI =0._JPRB
XRLXTE =0._JPRB
XRLXQ =0._JPRB
XRLXO3 =0._JPRB
XRLXLP =0._JPRB

IF(ALLOCATED(TRLXVO)) THEN
  TRLXVO =0._JPRB
ENDIF
IF(ALLOCATED(TRLXDI)) THEN
  TRLXDI =0._JPRB
ENDIF
IF(ALLOCATED(TRLXTE)) THEN
  TRLXTE =0._JPRB
ENDIF
IF(ALLOCATED(TRLXQ)) THEN
  TRLXQ =0._JPRB
ENDIF
IF(ALLOCATED(TRLXLP)) THEN
  TRLXLP =0._JPRB
ENDIF
IF(ALLOCATED(TRLXO3)) THEN
  TRLXO3 =0._JPRB
ENDIF

NPRLX3D = 0

! default values in degrees
ALATRLX1=90._JPRB
ALATRLX2=-90._JPRB
ALONRLX1=0._JPRB
ALONRLX2=360._JPRB

AXRLX=-0.5_JPRB
AYRLX=-0.5_JPRB
AZRLX=1.0_JPRB

NRLXLMIN=1
NRLXLMAX=NFLEVG
! default for velocity is same levels as for temperature
NRLXLMINU=-1
NRLXLMAXU=-1

NRLXSMAX=NSMAX

!      ----------------------------------------------------------------

!*       2.    MODIFIES DEFAULT VALUES.
!              ------------------------

CALL POSNAM(NULNAM,'NAMRLX')
READ       (NULNAM, NAMRLX)

LRLXG=(LRLXG.AND.(LRLXVO.OR.LRLXDI.OR.LRLXTE&
 & .OR.LRLXLP ))   

IF(NRLXLMINU == -1) NRLXLMINU=NRLXLMIN
IF(NRLXLMAXU == -1) NRLXLMAXU=NRLXLMAX

! count the number of 3D and 2D fields, respectively, being relaxed
IF(LRLXVO) NPRLX3D = NPRLX3D + 1
IF(LRLXDI) NPRLX3D = NPRLX3D + 1
IF(LRLXTE) NPRLX3D = NPRLX3D + 1
IF(LRLXQ)  NPRLX3D = NPRLX3D + 1
IF(LRLXQI)  NPRLX3D = NPRLX3D + 1
IF(LRLXQL)  NPRLX3D = NPRLX3D + 1
IF(LRLXQC)  NPRLX3D = NPRLX3D + 1
IF(LRLXO3)  NPRLX3D = NPRLX3D + 1
IF(LRLXLP) NPRLX2D = NPRLX2D + 1

! total number of events for which reference fields are available
WRITE(KULOUT,*) NSTOP,TSTEP,NFRLXU
NRLXAN=(NSTOP*TSTEP)/(NFRLXU*3600.0_JPRB)
NRLXAN=NRLXAN+1

! convert to radians
ALATRLX1=ALATRLX1*RPI/180._JPRB
ALATRLX2=ALATRLX2*RPI/180._JPRB
ALONRLX1=ALONRLX1*RPI/180._JPRB
ALONRLX2=ALONRLX2*RPI/180._JPRB

! convert relaxation time-scales from hours to time-steps
XRLXVO=XRLXVO*TSTEP/3600._JPRB
XRLXDI=XRLXDI*TSTEP/3600._JPRB
XRLXTE=XRLXTE*TSTEP/3600._JPRB
XRLXQ=XRLXQ*TSTEP/3600._JPRB
XRLXO3=XRLXO3*TSTEP/3600._JPRB
XRLXLP=XRLXLP*TSTEP/3600._JPRB

! set surface pressure relaxation to zero if necessary
IF (NRLXLMAX /= NFLEVG) THEN
  XRLXLP=0._JPRB
  WRITE(UNIT=KULOUT,FMT='('' WARNING: XRLXLP SET TO ZERO '')')
ENDIF

!      -----------------------------------------------------------

!*       3.    PRINTS FINAL VALUES.
!              --------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMRLX '')')
WRITE(UNIT=KULOUT,FMT='(&
 & '' LRLXG  ='',L2,9X,/,&
 & '' XRLXVO ='',E11.4,'' XRLXDI ='',E11.4,'' XRLXTE ='',E11.4,&
 & '' XRLXQ ='',E11.4,/,'' XRLXLP ='',E11.4,'' XRLXO3 ='',E11.4,/,&
 & '' LRLXVO ='',L2,9X,'' LRLXDI ='',L2,9X,'' LRLXTE ='',L2,9X,&
 & '' LRLXQ ='',L2,9X,'' LRLXLP ='',L2,&
 & '' LRLXQI ='',L2,9X,'' LRLXQL ='',L2,9X,'' LRLXQC ='',L2,9X&
 & ,'' LRLXO3 ='',L2,9X )')&
 & LRLXG,&
 & XRLXVO, XRLXDI, XRLXTE,XRLXQ, XRLXLP, XRLXO3,&
 & LRLXVO, LRLXDI, LRLXTE,LRLXQ, LRLXLP,&
 & LRLXQI,LRLXQL,LRLXQC,LRLXO3
WRITE(UNIT=KULOUT,FMT='(&
 & '' CRLXPATHGG = '',A,1X )')&
 & TRIM(CRLXPATHGG)
WRITE(UNIT=KULOUT,FMT='(&
 & '' CRLXPATHSH = '',A,1X )')&
 & TRIM(CRLXPATHSH)
WRITE(UNIT=KULOUT,FMT='('' NFRLXU ='',I3,1X,'' NRLXAN ='',I3,1X)')&
 & NFRLXU, NRLXAN
WRITE(UNIT=KULOUT,FMT='('' NPRLX3D ='',I3,1X,'' NPRLX2D ='',I3,1X)')&
 & NPRLX3D, NPRLX2D
WRITE(KULOUT,*)'ALATRLX1 (deg): ',ALATRLX1*180._JPRB/RPI
WRITE(KULOUT,*)'ALATRLX2 (deg): ',ALATRLX2*180._JPRB/RPI
WRITE(KULOUT,*)'ALONRLX1 (deg): ',ALONRLX1*180._JPRB/RPI
WRITE(KULOUT,*)'ALONRLX2 (deg): ',ALONRLX2*180._JPRB/RPI

WRITE(KULOUT,*)'AXRLX: ',AXRLX
WRITE(KULOUT,*)'AYRLX: ',AYRLX
WRITE(KULOUT,*)'AZRLX: ',AZRLX

WRITE(KULOUT,*)'NRLXLMIN: ',NRLXLMIN
WRITE(KULOUT,*)'NRLXLMAX: ',NRLXLMAX
WRITE(KULOUT,*)'NRLXLMINU: ',NRLXLMINU
WRITE(KULOUT,*)'NRLXLMAXU: ',NRLXLMAXU

WRITE(KULOUT,*)'NRLXSMAX: ',NRLXSMAX

!      -----------------------------------------------------------

!*       4.    VARIOUS CHECKS.
!              ------------------

IF (LRLXG.AND.(NFRLXG /= 2)) THEN
  WRITE(KULOUT,FMT='('' WHEN RELAXATION, NFRLXG MUST BE =2 '')')
  WRITE(KULOUT,FMT='('' LINEAR INTERPOLATION ONLY '')')
  CALL ABOR1('SURLX: ABOR1 CALLED')
ENDIF

!IF(LRLXG.AND.(NFRLXG == 0)) THEN
!  WRITE(KULOUT,FMT='('' WHEN RELAXATION, NFRLXG MUST BE >0 '')')
!  CALL ABOR1('SURLX: ABOR1 CALLED')
!ENDIF

!      ----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURLX',1,ZHOOK_HANDLE)
END SUBROUTINE SURLX
