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

SUBROUTINE SUMCUF(YDMCUF,YDDIM,YDRIP)

! Purpose:
! --------  
!   *SUMCUF*  - Setup of the parameters for the monitoring of 
!               update frequency of the coupling files for ALADIN

! Interface:
! ----------
!   *CALL* *SUMCUF*

! Externals:
! ----------
!   None.

! Method:
! -------
!   See documentation.

! Reference:
! ----------
!   Mon. Wea. Rev., 132, 2130-2141

! Author:
! -------
!   31-Aug-2004 Piet Termonia                    *RMI Belgium*

! Modicifations:
! --------------
!  Rachida El Ouaraini & Ryad El Khatib: 05-07-27 Add YDMCUF%LREACUF to read CUF fields. 
!  Ryad El Khatib 27-Feb-2008 Bugfix
!  O.Spaniel 27-Sep-2010 Bugfix - size for the declaration of arrays 
!  K. Yessad (Aug 2013): merge PARMCUF and YOMMCUF into YOMMCUF.
!  K. Yessad (July 2014): Move some variables.
! End Modifications
!------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMRIP   , ONLY : TRIP
USE YOMMCUF   ,ONLY : TMCUF  ,JPMFNR   ,JPMFOR  
USE YOMLUN    ,ONLY : NULOUT   ,NULNAM
USE YOMCST    ,ONLY : RPI

!------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TRIP)  ,INTENT(IN):: YDRIP
TYPE(TMCUF),TARGET , INTENT(INOUT) :: YDMCUF

INTEGER(KIND=JPIM)              :: JOR      ,JNR      ,JA
INTEGER(KIND=JPIM)              :: IN       ,ICI      ,ICII
INTEGER(KIND=JPIM), ALLOCATABLE :: ICFC(:)  ,IBI(:)
REAL(KIND=JPHOOK)                 :: ZHOOK_HANDLE
REAL(KIND=JPRB)                 :: ZAMAX    ,ZF       ,ZOMC     ,ZARG     ,&
                                &  ZD       ,ZTR      ,ZTI      ,ZTTR     ,&
                                &  ZQ
REAL(KIND=JPRB),    ALLOCATABLE :: ZPOLER(:),ZPOLEI(:),ZPOLEN(:),ZFAC(:)  ,&
                                &  ZBR(:)   ,ZBI(:)
LOGICAL                         :: LLCONT   ,LLJ

!-----------------------------------------------------------------
LOGICAL, POINTER             :: LMCUF, LREACUF                             
REAL(KIND=JPRB), POINTER     :: RMCUFI(:)           
REAL(KIND=JPRB)              :: RMCUFA   ,RMCUFB  
INTEGER(KIND=JPIM), POINTER  :: NCUFOR   ,NCUFNR
!------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "nammcuf.nam.h"

!------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUMCUF',0,ZHOOK_HANDLE)
ASSOCIATE(NSPEC2=>YDDIM%NSPEC2, &
 & TSTEP=>YDRIP%TSTEP)
!------------------------------------------------------------------

! 1. Defaults and namelist:
! -------------------------

LMCUF   => YDMCUF%LMCUF
LREACUF => YDMCUF%LREACUF
NCUFNR  => YDMCUF%NCUFNR
NCUFOR  => YDMCUF%NCUFOR
RMCUFI  => YDMCUF%RMCUFI

!YDMCUF%LMCUF   =.FALSE.
!YDMCUF%LREACUF =.FALSE.
!YDMCUF%NCUFNR  =1
!YDMCUF%NCUFOR  =2
!YDMCUF%RMCUFI  =10800._JPRB


LMCUF   =.FALSE.
LREACUF =.FALSE.
NCUFNR  =1
NCUFOR  =2
RMCUFI  =10800._JPRB



ZAMAX=99999._JPRB
DO JNR=1,JPMFNR
  YDMCUF%RMCUFA(0,JNR)=2._JPRB*ZAMAX
  YDMCUF%RMCUFA(1:JPMFOR,JNR)=0._JPRB
  YDMCUF%RMCUFB(1:JPMFOR,JNR)=0._JPRB
ENDDO
CALL POSNAM(NULNAM,'NAMMCUF')
READ(NULNAM,NAMMCUF)

IF (.NOT.YDMCUF%LMCUF) THEN
! Reset number of filters to zero -REK-
  YDMCUF%NCUFNR=0
ENDIF

WRITE(UNIT=NULOUT,FMT=*) 'YDMCUF%NCUFNR=',YDMCUF%NCUFNR,'YDMCUF%NCUFOR=',YDMCUF%NCUFOR,'YDMCUF%LMCUF=', &
&                          YDMCUF%LMCUF,'YDMCUF%LREACUF=',YDMCUF%LREACUF
WRITE(UNIT=NULOUT,FMT=*) 'YDMCUF%LMCUF =',YDMCUF%LMCUF,'YDMCUF%LREACUF=',YDMCUF%LREACUF
WRITE(UNIT=NULOUT,FMT=*) 'YDMCUF%RMCUFI=',YDMCUF%RMCUFI

IF (YDMCUF%NCUFOR>JPMFOR) CALL ABOR1('SUMCUF: maximum order of filter exceeded!')
IF (YDMCUF%NCUFNR>JPMFNR) CALL ABOR1('SUMCUF: maximum nummer of filters exceeded!')
IF (YDMCUF%LMCUF .AND. YDMCUF%NCUFNR <= 0) CALL ABOR1('SUMCUF:YDMCUF%LMCUF=T BUT NO FILTERS !')

! 2. Compute coefficients:
! ------------------------

ALLOCATE(ZPOLER(1:2*YDMCUF%NCUFOR),ZPOLEI(1:2*YDMCUF%NCUFOR),ZPOLEN(1:2*YDMCUF%NCUFOR))
ALLOCATE(ZFAC(0:YDMCUF%NCUFOR),ICFC(1:YDMCUF%NCUFOR),IBI(0:YDMCUF%NCUFOR))
ALLOCATE(ZBR(1:YDMCUF%NCUFOR),ZBI(1:YDMCUF%NCUFOR))
ZFAC(0)=1.0_JPRB
DO JOR=1,YDMCUF%NCUFOR
  ZFAC(JOR)=REAL(JOR,JPRB)*ZFAC(JOR-1)
ENDDO
DO JNR=1,YDMCUF%NCUFNR
  IF (YDMCUF%RMCUFA(0,JNR) > ZAMAX) THEN
    WRITE(NULOUT,*) 'FILTER COEFFICIENTS COMPUTED FOR FILTER',JNR
    ZOMC=0.9_JPRB*RPI/YDMCUF%RMCUFI(JNR)
    IN=1
    ZF=TAN(ZOMC*TSTEP*0.5_JPRB)
    DO JA=-YDMCUF%NCUFOR+1,YDMCUF%NCUFOR
      ZARG=RPI*REAL(2*JA-1,JPRB)/REAL(2*YDMCUF%NCUFOR,JPRB)
      ZD=1._JPRB+ZF*ZF-2._JPRB*ZF*SIN(ZARG)
      ZPOLER(IN)=1._JPRB-ZF*ZF
      ZPOLER(IN)=ZPOLER(IN)/ZD
      ZPOLEI(IN)=2._JPRB*ZF*COS(ZARG)
      ZPOLEI(IN)=ZPOLEI(IN)/ZD
      ZPOLEN(IN)=ZPOLER(IN)**2+ZPOLEI(IN)**2 
      IF (ZPOLEN(IN) > 1._JPRB) IN=IN+1
    ENDDO
    DO JOR=1,YDMCUF%NCUFOR
      ZBR(JOR)=0.0_JPRB
      ZBI(JOR)=0.0_JPRB
      ICFC(JOR)=0
      IBI(0)=YDMCUF%NCUFOR
      DO ICI=1,JOR
        IBI(ICI)=ICI
      ENDDO
      LLCONT=.TRUE.
      DO WHILE(LLCONT)
        ZTR=1.0_JPRB
        ZTI=0.0_JPRB
        DO ICI=1,JOR
          ZTTR=-ZTR*ZPOLER(IBI(ICI))-ZTI*ZPOLEI(IBI(ICI))
          ZTTR=ZTTR/ZPOLEN(IBI(ICI))
          ZTI=-ZTI*ZPOLER(IBI(ICI))+ZTR*ZPOLEI(IBI(ICI))
          ZTI=ZTI/ZPOLEN(IBI(ICI))
          ZTR=ZTTR
        ENDDO
        ZBR(JOR)=ZBR(JOR)+ZTR
        ZBI(JOR)=ZBI(JOR)+ZTI
        ICFC(JOR)=ICFC(JOR)+1
        LLJ=.FALSE.
        DO ICI=1,JOR
          IF (IBI(ICI-1) < IBI(ICI)-1) THEN
            IBI(ICI-1)=IBI(ICI-1)+1
            DO ICII=1,ICI-2
              IBI(ICII)=ICII
            ENDDO
            LLJ=.TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT.LLJ) THEN
          IF (IBI(JOR) < YDMCUF%NCUFOR) THEN
            IBI(JOR)=IBI(JOR)+1
            DO ICI=1,JOR-1
              IBI(ICI)=ICI
            ENDDO
          ELSE
            LLCONT=.FALSE.
          ENDIF
       ENDIF
      ENDDO
    ENDDO
    DO JOR=1,YDMCUF%NCUFOR
      YDMCUF%RMCUFB(JOR,JNR)=ZBR(JOR)
    ENDDO
    ZQ=1._JPRB
    DO JOR=1,YDMCUF%NCUFOR
      ZQ=ZQ+(-1._JPRB)**JOR*ZBR(JOR)
    ENDDO
    ZQ=-ZQ/REAL(2**YDMCUF%NCUFOR,JPRB)
    DO JOR=0,YDMCUF%NCUFOR
      YDMCUF%RMCUFA(JOR,JNR)=ZQ*(-1._JPRB)**JOR*ZFAC(YDMCUF%NCUFOR)
      YDMCUF%RMCUFA(JOR,JNR)=YDMCUF%RMCUFA(JOR,JNR)/ZFAC(JOR)
      YDMCUF%RMCUFA(JOR,JNR)=YDMCUF%RMCUFA(JOR,JNR)/ZFAC(YDMCUF%NCUFOR-JOR)
    ENDDO
  ENDIF
ENDDO
DEALLOCATE(ZPOLER,ZPOLEI,ZPOLEN)
DEALLOCATE(ZFAC,ICFC,IBI)
DEALLOCATE(ZBR,ZBI)

! 3. Initialisation:
! ------------------

IF (YDMCUF%LMCUF) THEN
  DO JNR=1,YDMCUF%NCUFNR
    WRITE(UNIT=NULOUT,FMT=*) 'YDMCUF%RMCUFA(',JNR,')=',YDMCUF%RMCUFA(0:YDMCUF%NCUFOR,JNR)
    WRITE(UNIT=NULOUT,FMT=*) 'YDMCUF%RMCUFB(',JNR,')=',YDMCUF%RMCUFB(1:YDMCUF%NCUFOR,JNR)
  ENDDO
  ALLOCATE(YDMCUF%RMCUFSP(1:NSPEC2,0:YDMCUF%NCUFOR))
  ALLOCATE(YDMCUF%RMCUFFP(1:NSPEC2,0:YDMCUF%NCUFOR,1:YDMCUF%NCUFNR))
  YDMCUF%RMCUFSP(:,:)=0._JPRB
  YDMCUF%RMCUFFP(:,:,:)=0._JPRB
ENDIF

!------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUMCUF',1,ZHOOK_HANDLE)
END SUBROUTINE SUMCUF
