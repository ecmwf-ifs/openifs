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

SUBROUTINE INICOU(YDML_AOC,YDEPHY,KASTP,KEXCH,KSTEP)

!**** *INICOU*  - Initialize coupled mode communication

!     Purpose.
!     --------
!     Exchange process id and timestep information
!     between AGCM, OGCM and COUPLER.

!**   Interface.
!     ----------
!       *CALL*  *INICOU(...)*

!     Input:
!     -----
!       KASTP  : total number of timesteps in atmospheric model
!       KEXCH  : frequency of exchange for the fluxes (in time steps)
!       KSTEP  : time step

!     Output:
!     ------

!     Method:
!     ------
!     Depends on value of NOACOMM

!     Externals:
!     ---------
!     GETPID, PBOPEN, PBREAD, PBWRITE,
!     SVIPC_OPEN, SVIPC_READ, SVIPC_WRITE

!     Reference:
!     ---------
!     See Epicoa 0803 (1992)

!     Author:
!     -------
!     Laurent Terray  92-09-01

!     Modifications:
!     --------------
!     T. Stockdale  07-02-05  Add initial sea-ice, change synchronization
!     -----------------------------------------------------------

USE MODEL_ATMOS_OCEAN_COUPLING_MOD , ONLY : MODEL_ATMOS_OCEAN_COUPLING_TYPE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMSC   , ONLY : NINTLEN
USE YOEPHY   , ONLY : TEPHY
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : MYPROC

IMPLICIT NONE

TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_ATMOS_OCEAN_COUPLING_TYPE),INTENT(INOUT):: YDML_AOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KASTP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEXCH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
EXTERNAL GETPID
LOGICAL :: LLFILE, LLDONE

INTEGER(KIND=JPIM) :: IMESS(4)
!        SORRY DOCTOR, THE NEXT LINE IS A CRAY/FUJITSU FUNCTION
#ifdef WITH_NEMO
INTEGER(KIND=JPIM) :: GETPID
#endif

INTEGER(KIND=JPIM) :: IA, ICOUNT, IRET, IRET1, IRET2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INICOU',0,ZHOOK_HANDLE)
ASSOCIATE(NCULMR=>YDML_AOC%YRCOU%NCULMR, NCULMW=>YDML_AOC%YRCOU%NCULMW, &
 & NPIAT=>YDML_AOC%YRCOU%NPIAT, &
 & LECURR=>YDEPHY%LECURR, &
 & LMCCICEIC=>YDML_AOC%YRMCC%LMCCICEIC, NOACOMM=>YDML_AOC%YRMCC%NOACOMM, &
 & YDCOU=>YDML_AOC%YRCOU)

IF(MYPROC == 1) THEN

!*    1. Open statements with correct read/write attributes
!        --------------------------------------------------

!  With either pipe method, the pipes must be created by the script 
!  before the coupled model is started, to guarantee their existence.

  IF(NOACOMM == 1) THEN

!  Use unit numbers reserved for atmosphere DA, for full Fortran 90 
!  compatibility. Fortran pipes often not properly supported.
      WRITE(NULOUT,*) 'FORTRAN pipes no longer supported by IFS'
      WRITE(NULOUT,*) 'Please use C pipes (NOACOMM=2) or another method'
      CALL ABOR1('ABORT in INICOU')

  ELSEIF((NOACOMM == 2).OR.(NOACOMM == 4)) THEN

    CALL PBOPEN(NCULMR,'Preadm01','w',IRET1)
    CALL PBOPEN(NCULMW,'Pwritm01','r',IRET2)
    IF((IRET1 < 0).OR.(IRET2 < 0)) THEN
      WRITE(NULOUT,*) 'Problems opening C pipes'
      CALL ABOR1('ABORT in INICOU')
    ENDIF

  ELSEIF(NOACOMM == 3) THEN

!  Wait for OASIS to signal creation of shared memory pool before we 
!  attempt to attach to it.
    LLFILE = .FALSE.
    ICOUNT = 0
    DO WHILE (.NOT.LLFILE )
      INQUIRE(FILE='DUMMY_SIPC',EXIST=LLFILE)
      CALL SLEEP(2)
      ICOUNT=ICOUNT+1
      IF(ICOUNT  ==  120) THEN
        WRITE(NULOUT,*) 'DUMMY_SIPC file not found'
        WRITE(NULOUT,*) 'Time out after waiting 240 seconds'
        CALL ABOR1('ABORT in INICOU')
      ENDIF
    ENDDO
    CALL SVIPC_OPEN(NCULMR,'ifs   R',0,IRET1)
    CALL SVIPC_OPEN(NCULMW,'ifs   W',0,IRET2)
    IF((IRET1 < 0).OR.(IRET2 < 0)) THEN
      WRITE(NULOUT,*) 'Problems attaching to memory pool'
      CALL ABOR1('ABORT in INICOU')
    ENDIF

  ELSE
    WRITE(NULOUT,*) 'Ocean-atmosphere communication method not known'
    WRITE(NULOUT,*) 'NOACOMM = ',NOACOMM
    CALL ABOR1('ABORT in INICOU')
  ENDIF

!*    2. Send timestep and pid info to coupler
!        -------------------------------------

#ifdef WITH_NEMO
  NPIAT = GETPID()
#endif
  IMESS(1) = KASTP
  IMESS(2) = KEXCH
  IMESS(3) = KSTEP
  IMESS(4) = NPIAT

  IF((NOACOMM == 2).OR.(NOACOMM == 4)) THEN
    CALL PBWRITE(NCULMR,IMESS,4*NINTLEN,IRET)
    CALL PBFLUSH(NCULMR)
  ELSEIF(NOACOMM == 3) THEN
    CALL SVIPC_WRITE(NCULMR,IMESS,4*NINTLEN,IRET)
  ENDIF
  IF(IRET < 0) THEN
    WRITE(NULOUT,*) 'Problems sending timestep to coupler'
    CALL ABOR1('ABORT in INICOU')
  ENDIF

!*    3. Read timestep and pid info from coupler
!        ---------------------------------------

  IF(NOACOMM == 2) THEN
    CALL PBREAD(NCULMW,IMESS,4*NINTLEN,IRET)
  ELSEIF(NOACOMM == 3) THEN
    CALL SVIPC_READ(NCULMW,IMESS,4*NINTLEN,IRET)
  ELSEIF(NOACOMM == 4) THEN
    LLDONE=.FALSE.
    IA=0
    DO WHILE (.NOT.LLDONE)
      CALL PBREAD(NCULMW,IMESS,4*NINTLEN,IRET)
      IF(IRET <= 0) THEN
        IA=IA+1
        WRITE(NULOUT,*) 'IRET is ',IRET
        WRITE(NULOUT,*) 'Sleep ',IA,' seconds'
        CALL FLUSH(NULOUT)
        CALL SLEEP(IA)
        IF(IA > 40) THEN
          WRITE(NULOUT,*) 'Timed out waiting for coupler'
          CALL ABOR1('ABORT in INICOU')
        ENDIF
      ELSE
        LLDONE=.TRUE.
      ENDIF
    ENDDO
  ENDIF
  IF(IRET < 0) THEN
    WRITE(NULOUT,*) 'Problems receiving info from coupler'
    CALL ABOR1('ABORT in INICOU')
  ENDIF

  WRITE(NULOUT,FMT='('' INICOU printout: '')')
  WRITE(NULOUT,FMT='('' Communication established with OASIS '')')
  WRITE(NULOUT,FMT='('' OASIS pid is : '',I6)')  IMESS(4)
  CALL FLUSH(NULOUT)

!*    4. OPEN PIPES FOR SIGNALLING FIELD TRANSFERS
!     --------------------------------------------

!  In the case NOACOMM=2, the order of the open statements is
!  important: it must match the order of the fields given in
!  the OASIS namelist input.

  IF((NOACOMM == 2).OR.(NOACOMM == 4)) THEN
    IF(LMCCICEIC) THEN
      CALL PBOPEN(YDCOU%NCULF(0),'IICEATMO','w',IRET)
    ENDIF
    CALL PBOPEN(YDCOU%NCULF(1),'SSTATMOS','r',IRET)
    CALL PBOPEN(YDCOU%NCULF(2),'ICEATMOS','r',IRET)
    CALL PBOPEN(YDCOU%NCULF(3),'QNET_ATM','w',IRET)
    CALL PBOPEN(YDCOU%NCULF(4),'QSOL_ATM','w',IRET)
    CALL PBOPEN(YDCOU%NCULF(5),'TAUX_ATM','w',IRET)
    CALL PBOPEN(YDCOU%NCULF(6),'TAUY_ATM','w',IRET)
    CALL PBOPEN(YDCOU%NCULF(7),'PMEV_ATM','w',IRET)
    IF(LECURR) THEN
      CALL PBOPEN(YDCOU%NCULF(8),'UV_ATMOS','r',IRET)
      CALL PBOPEN(YDCOU%NCULF(9),'VV_ATMOS','r',IRET)
    ENDIF
  ELSEIF(NOACOMM == 3) THEN
    IF(LMCCICEIC) THEN
      CALL SVIPC_OPEN(YDCOU%NCULF(0),'PIICEATMO',0,IRET)
    ENDIF
    CALL SVIPC_OPEN(YDCOU%NCULF(1),'PSSTATMOS',0,IRET)
    CALL SVIPC_OPEN(YDCOU%NCULF(2),'PICEATMOS',0,IRET)
    CALL SVIPC_OPEN(YDCOU%NCULF(3),'PQNET_ATM',0,IRET)
    CALL SVIPC_OPEN(YDCOU%NCULF(4),'PQSOL_ATM',0,IRET)
    CALL SVIPC_OPEN(YDCOU%NCULF(5),'PTAUX_ATM',0,IRET)
    CALL SVIPC_OPEN(YDCOU%NCULF(6),'PTAUY_ATM',0,IRET)
    CALL SVIPC_OPEN(YDCOU%NCULF(7),'PPMEV_ATM',0,IRET)
    IF(LECURR) THEN
      CALL SVIPC_OPEN(YDCOU%NCULF(8),'PUV_ATMOS',0,IRET)
      CALL SVIPC_OPEN(YDCOU%NCULF(9),'PVV_ATMOS',0,IRET)
    ENDIF
  ENDIF
  IF(IRET < 0) THEN
    WRITE(NULOUT,*) 'Problems opening field communication'
    CALL ABOR1('ABORT in INICOU')
  ENDIF

!*    5. Wait for OASIS2 signal that SST is ready
!     -------------------------------------------

! moved to routine UPDCLIE

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INICOU',1,ZHOOK_HANDLE)
END SUBROUTINE INICOU
