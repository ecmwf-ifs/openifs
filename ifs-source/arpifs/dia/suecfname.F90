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

SUBROUTINE SUECFNAME(LDINC,CDLEVTYPE,KSTEP,PTSTEP,KSTOP,CDFNGG,CDFNSH)

!**** *SUECFNAME* - setup name of output files

!     Purpose.
!     --------
!     CREATE FILENAMES BEFORE WRITING

!**   Interface.
!     ----------
!        *CALL* *SUECFNAME*(...)

!        Explicit arguments :
!        ------------------
!            CDLEVTYPE - type of level for fields

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
!        Mats Hamrud  *ECMWF*
!       Original : 01-12-14 Adapted from WR****

!     Modifications.
!     --------------
!      R. El Khatib 25-Apr-2017 from wroutgpgb
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMCT0   , ONLY : NCONF, LECFPOS, CNMEXP
USE YOMVAR   , ONLY : MCGLVEC, MBGVEC, LTRREF, LTWANA, LTWGRA, LTWINC, LTWBGV, LTWCGL  
USE YOMOPH0  , ONLY : CFNAN, CFNBGV, CFNCGL, CFNGR, CFNINSH, CFNRF, CFPATH, CFNGG, CFNINGG, CFNUA
USE ALGORITHM_STATE_MOD  , ONLY : GET_NSIM4D

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL,            INTENT(IN)  :: LDINC
CHARACTER(LEN=1),   INTENT(IN)  :: CDLEVTYPE 
INTEGER(KIND=JPIM), INTENT(IN)  :: KSTEP
REAL(KIND=JPRB),    INTENT(IN)  :: PTSTEP
INTEGER(KIND=JPIM), INTENT(IN)  :: KSTOP
CHARACTER(LEN=*),   INTENT(OUT) :: CDFNGG
CHARACTER(LEN=*),   INTENT(OUT) :: CDFNSH

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IINC, IMTS
CHARACTER(LEN=5) :: CLR1, CLR4, CLR5
CHARACTER(LEN=4) :: CLR2, CLR1INC
CHARACTER(LEN=LEN(CDFNGG)) :: CLFNGG, CLFNUA, CLFNOC, CLFBGV
CHARACTER(LEN=LEN(CDFNSH)) :: CLFNRF, CLFANA, CLFINC, CLFGRA, CLFCGL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "suofname.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUECFNAME',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------



CLFNGG = ' '
CLR1 = CFNGG(1:5)
CLR4 = CFNUA(1:5)
IF (LTWINC) CLR1INC = CFNINGG(1:4)
CLR2 = CNMEXP
CLR5 = 'ICMOC'   !KPP
IF (LDINC) THEN
  IINC = NINT(REAL(KSTEP,JPRB)*PTSTEP/3600._JPRB)
ELSE
  IINC = KSTEP
ENDIF

IF (KSTEP < 0) THEN
  IF (LECFPOS) THEN
    WRITE(CLFNGG,'(A5,A4,I5.4)') CLR1,CLR2,KSTEP
    WRITE(CLFNUA,'(A5,A4,I5.4)') CLR4,CLR2,KSTEP
  ELSE
    WRITE(CLFNGG,'(A5,A4,I5.4)') CLR1,CLR2,KSTEP
  ENDIF
  WRITE(CLFNOC,'(A5,A4,I5.4)') CLR5,CLR2,KSTEP      !KPP  
ELSE
  IF(LTWINC) THEN
    IMTS = NINT(REAL(KSTEP,JPRB)*PTSTEP/60._JPRB)
    IINC = (IMTS/60)*100+MOD(IMTS,60)
    WRITE(CLFNGG,'(A4,A4,I3.3,''+'',I6.6)') CLR1INC,CLR2,GET_NSIM4D(),IINC
  ELSE
    WRITE(CLFNGG,'(A5,A4,''+'',I6.6)') CLR1,CLR2,IINC
  ENDIF
  WRITE(CLFNUA,'(A5,A4,''+'',I6.6)') CLR4,CLR2,IINC
  WRITE(CLFNOC,'(A5,A4,''+'',I6.6)') CLR5,CLR2,IINC      !KPP  
ENDIF

!IF( LECFPOS .AND. NFPOS == 2 )  THEN
! ecmwf uses nfpos=0 (=>lecfpos=F) or nfpos=2 only
IF( LECFPOS)  THEN
  IF(CDLEVTYPE == 's') THEN
    CDFNGG = CLFNGG
  ELSE
    CDFNGG = CLFNUA
  ENDIF
ELSE
  CDFNGG = CLFNGG
ENDIF

IF(CDLEVTYPE == 'o') THEN    !KPP
  CDFNGG = CLFNOC            !KPP
ENDIF                      !KPP

IF(LTWBGV) THEN
  CLFBGV=CFNBGV
  WRITE(CLFBGV(7:),'(I3.3)') MBGVEC
  WRITE(CLFBGV(1:2),'(A2)') 'UA'
  CDFNGG=CLFBGV
ENDIF




IF (LTRREF) THEN
  
  CLR1 = CFNRF(1:5)
  CLR2 = CFNRF(6:9)
  IF (KSTEP < 0) THEN
    WRITE(CLFNRF,'(A5,A4,I7.6)') CLR1,CLR2,KSTEP
  ELSE
    IF (LDINC) THEN
      IINC = NINT(REAL(KSTEP,JPRB)*PTSTEP/3600._JPRB)
    ELSE
      IINC = KSTEP
    ENDIF
    WRITE(CLFNRF,'(A5,A4,''+'',I6.6)') CLR1,CLR2,IINC
  ENDIF
  CDFNSH=CLFNRF

ELSEIF(LTWANA) THEN
    
  CLFANA = CFNAN
  IF (KSTEP < 0) THEN
    CLFANA(10:10) = '-'
    IINC = ABS(KSTEP)
  ELSEIF (NCONF/100 == 1) THEN
    IMTS = NINT(REAL(KSTEP,JPRB)*PTSTEP/60._JPRB)
    IINC = (IMTS/60)*100+MOD(IMTS,60)
  ELSE
    IINC = KSTEP
  ENDIF
  WRITE(CLFANA(7:9),'(I3.3)') GET_NSIM4D()
  WRITE(CLFANA(11:),'(I6.6)') IINC
  CDFNSH=CLFANA
    
ELSEIF(LTWINC) THEN
  
  CLFINC = CFNINSH
  ! Output files labelled with time step in configuration 501 (TL evolution)
  IF (NCONF == 501) THEN
    IINC = KSTEP
  ELSE
    IMTS = NINT(REAL(KSTEP,JPRB)*PTSTEP/60._JPRB)
    IINC = (IMTS/60)*100+MOD(IMTS,60)
  ENDIF
  WRITE(CLFINC(9:11),'(I3.3)') GET_NSIM4D()
  IF (IINC >= 0) THEN
    WRITE(CLFINC(13:),'(I6.6)') IINC
  ELSE
    WRITE(CLFINC(12:12),'(A1)') '-'
    WRITE(CLFINC(13:),'(I6.6)') -IINC
  ENDIF
  CDFNSH=CLFINC
  
ELSEIF(LTWGRA) THEN
    
  CLFGRA = CFNGR
  IF (KSTEP < 0) THEN
    CLFGRA(10:10) = '-'
    IINC = ABS(KSTEP)
  ELSEIF (NCONF/100 == 1) THEN
    IMTS = NINT(REAL(KSTEP,JPRB)*PTSTEP/60._JPRB)
    IINC = (IMTS/60)*100+MOD(IMTS,60)
  ELSEIF (IABS(NCONF)/100 == 8) THEN
    CLFGRA(10:10) = '-'
    IINC = KSTOP-KSTEP
  ELSE
    IINC = KSTEP
  ENDIF
  WRITE(CLFGRA(7:9),'(I3.3)') GET_NSIM4D()
  WRITE(CLFGRA(11:),'(I6.6)') IINC
  CDFNSH=CLFGRA
  
ELSEIF(LTWBGV) THEN
  
  CLFBGV = CFNBGV
  WRITE(CLFBGV(7:),'(I3.3)') MBGVEC
  CDFNSH=CLFBGV
    
ELSEIF(LTWCGL) THEN
  
  CLFCGL = CFNCGL
  WRITE(CLFCGL(7:),'(I3.3)') MCGLVEC
  CDFNSH=CLFCGL
  
ELSE

  IF (LDINC) THEN
    CALL SUOFNAME(KSTEP=KSTEP,PTSTEP=PTSTEP,KDIGITS=6,CDPATH=CFPATH,CDLABEL='ICMSH',CDOFNAME=CDFNSH)
  ELSE
    CALL SUOFNAME(KSTEP=KSTEP,KDIGITS=6,CDPATH=CFPATH,CDLABEL='ICMSH',CDOFNAME=CDFNSH)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUECFNAME',1,ZHOOK_HANDLE)
END SUBROUTINE SUECFNAME
