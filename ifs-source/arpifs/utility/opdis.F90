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

SUBROUTINE OPDIS(CDCONF,CDCALLER,PCP,PVCP,PRT,PACP,PART,KSTEP,LDNHDYN)

!***  *OPDIS* - SUBROUTINE TO GENERATE AN OPERATOR DISPLAY

!     PURPOSE.
!     -------

!     *OPDIS* PREPARES A TEXT ARRAY CONTAINING RELEVANT INFORMATION
!     ABOUT THE CURRENT STATE OF MODEL EXECUTION.
!     THE TEXT IS WRITTEN TO A FILE SO THAT THE BUFFER IS IMMEDIATELY
!     FLUSHED.

!     INTERFACE
!     ---------

!     *CALL OPDIS*

!     CHARACTER*9 CDCONF
!   1 : configuration of IOPACK     IO handling
!   2 : configuration of LTINV      inverse Legendre transform
!   3 : configuration of FTINV      inverse Fourier transform (under SCAN2M*)
!   4 : configuration of CPG        grid point computations   (under SCAN2M*)
!   5 : configuration of POS        post processing           (under SCAN2M*)
!   6 : configuration of OBS        comparison to observations(under SCAN2M*)
!   7 : configuration of FTDIR      direct Fourier transform  (under SCAN2M*)
!   8 : configuration of LTDIR      direct Legendre transform
!   9 : configuration of SPC,SPHDG  spectral space computations

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        D.Dent    ECMWF

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!        G.Mozdzynski  09-Dec-2004 Remove ifs.disp
!      F. Vana  05-Mar-2015  Support for single precision
!      R. El khatib 22-Feb-2016 argument KSTEP instead of module variable NSTEP
!      J. Vivoda July 2018 add NH variables

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPIB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : CNDISPP
USE YOMSPNRM , ONLY : AVNRMDIV, AVNRMPD, AVNRMVD
USE YOMLUN   , ONLY : NULSTAT
USE YOMMP0   , ONLY : MYPROC

IMPLICIT NONE

CHARACTER(LEN=9)  ,INTENT(IN)    :: CDCONF
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDCALLER
REAL(KIND=JPRD)   ,INTENT(IN)    :: PCP
REAL(KIND=JPRD)   ,INTENT(IN)    :: PVCP
REAL(KIND=JPRD)   ,INTENT(IN)    :: PRT
REAL(KIND=JPRD)   ,INTENT(IN)    :: PACP
REAL(KIND=JPRD)   ,INTENT(IN)    :: PART
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
LOGICAL           ,INTENT(IN)    :: LDNHDYN
!     CDCONF: current configuration flags
!     CDCALLER: name of calling routine
!     PCP:   cp time of last step
!     PVCP:  vector cp time of last step
!     PRT:   real time of last step
!     PACP:  accumulated cp time
!     PART:  accumulated real time
!     KSTEP : model time step

CHARACTER(LEN=18) :: CLCALLER

CHARACTER (LEN = 200) ::  CLMESS
CHARACTER (LEN = 10) ::   CLCP,CLTIME
INTEGER(KIND=JPIM) :: IVALUES(8)

CHARACTER (LEN = 10) ::  CLDATEOD,CLTIMEOD,CLZONEOD
INTEGER(KIND=JPIM),SAVE :: IFIRST=0

INTEGER(KIND=JPIM) :: IB, IMCP, IMIT, ISCP, ISIT, JSTEP, ILEN

INTEGER(KIND=JPIB), EXTERNAL :: GETMAXHWM, GETMAXSTK, GETVMPEAK
INTEGER(KIND=JPIB) :: IGETHWM, IGETSTK, IGETVMP
CHARACTER(LEN=2) :: CLHWM, CLSTK, CLVMP

INTEGER(KIND=JPIB) :: IENERGY, IPOWER

REAL(KIND=JPRD) :: ZCP, ZIT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "ec_pmon.intfb.h"

!*         1.     CHECK SWITCHES AND BLANK DISPLAY ARRAY
!                 ----- -------- --- ----- ------- -----

IF (LHOOK) CALL DR_HOOK('OPDIS',0,ZHOOK_HANDLE)

IF (MYPROC == 1) THEN

  CALL GSTATS(1930,0)

  IF(IFIRST == 0) THEN
    IFIRST=1
    IB=INDEX(CNDISPP,' ') -1
    IF(IB >= 1) THEN
      OPEN(NULSTAT,FILE=CNDISPP(1:IB)//'ifs.stat',ACTION='WRITE',STATUS='REPLACE')
    ELSE
      OPEN(NULSTAT,FILE='ifs.stat',ACTION='WRITE',STATUS='REPLACE')
    ENDIF
    JSTEP=0
  ENDIF

  IF(IFIRST /= 0) THEN
    JSTEP=KSTEP
  ELSE
    JSTEP=0
  ENDIF

!     ---------------------------------------------------------------------

!*         0.5     FIND STEP TYPE

  CLCALLER = ' '
  CLCALLER = ADJUSTL(TRIM(CDCALLER))
  CLMESS=' '
  CALL DATE_AND_TIME(CLDATEOD,CLTIMEOD,CLZONEOD,IVALUES)
  WRITE(CLMESS(2:9),'(A,'':'',A,'':'',A)')&
   & CLTIMEOD(1:2),CLTIMEOD(3:4),CLTIMEOD(5:6)
  CLMESS(11:19)=CDCONF

  WRITE(CLMESS(21:),'(A,I6)') CLCALLER,JSTEP
  ILEN = LEN(TRIM(CLMESS))

!*         WRITE TIMING STATISTICS

  ZIT=PART
  IMIT=INT(ZIT/60._JPRD)
  ISIT=INT(ZIT-60._JPRD*IMIT)
  WRITE(CLTIME,'(I7,":",I2.2)') IMIT,ISIT
  ZCP=PACP
  IMCP=INT(ZCP/60._JPRD)
  ISCP=INT(ZCP-60._JPRD*IMCP)
  WRITE(CLCP  ,'(I7,":",I2.2)') IMCP,ISCP
!-- Approximate max heap usage on *this* MPI task
  IGETHWM = GETMAXHWM()
  IF (IGETHWM < 10_JPIB * 1048576_JPIB) THEN
     CLHWM = "KB" ! Display in KBytes
     IGETHWM = IGETHWM / 1024_JPIB
  ELSE
     CLHWM = "MB" ! Display in MBytes
     IGETHWM = IGETHWM / 1048576_JPIB
  ENDIF
!-- Approximate max stack size on *this* MPI task's master thread
  IGETSTK = GETMAXSTK()
  IF (IGETSTK < 10_JPIB * 1048576_JPIB) THEN
     CLSTK = "KB" ! Display in KBytes
     IGETSTK = IGETSTK / 1024_JPIB
  ELSE
     CLSTK = "MB" ! Display in MBytes
     IGETSTK = IGETSTK / 1048576_JPIB
  ENDIF
!-- Virtual memory peak on *this* MPI-task 
  IGETVMP = GETVMPEAK()
  IF (IGETVMP < 10_JPIB * 1073741824_JPIB) THEN
     CLVMP = "MB" ! Display in MBytes
     IGETVMP = IGETVMP / 1048576_JPIB
  ELSE
     CLVMP = "GB" ! Display in GBytes
     IGETVMP = IGETVMP / 1073741824_JPIB
  ENDIF
  CALL EC_PMON(IENERGY,IPOWER)
  IENERGY = IENERGY/1000_JPIB ! In kJoules (values are per node -- not per task [Cray])
  IF (LDNHDYN) THEN
    WRITE(CLMESS(ILEN+1:),207) PVCP,PCP,PRT,CLCP,CLTIME,AVNRMDIV,AVNRMVD,AVNRMPD
  ELSE
    WRITE(CLMESS(ILEN+1:),208) PVCP,PCP,PRT,CLCP,CLTIME,AVNRMDIV
  ENDIF
  ILEN = LEN(TRIM(CLMESS))
  WRITE(CLMESS(ILEN+2:),209) IGETHWM,CLHWM,&
      & IGETSTK,CLSTK,&
      & IGETVMP,CLVMP,&
      & IENERGY,IPOWER
207 FORMAT(3F9.3,2A10,3E21.14,1X)
208 FORMAT(3F9.3,2A10,E21.14,43X)
209 FORMAT(I7,A2,1X,I7,A2,1X,I7,A2,1X,I6,"kJ",1X,I4,"W")
  WRITE(NULSTAT,'(A)') TRIM(CLMESS)
  CALL FLUSH(NULSTAT)

  CALL GSTATS(1930,1)
ENDIF

IF (LHOOK) CALL DR_HOOK('OPDIS',1,ZHOOK_HANDLE)
END SUBROUTINE OPDIS
