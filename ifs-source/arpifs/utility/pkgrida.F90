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

SUBROUTINE PKGRIDA(YDGEOMETRY,KFIELDS,KBIT,PGP)

!**** *PKGRIDA*  - PACK GRIDPOINT DATA THROUGH ARPEGE FILE

!     PURPOSE.
!     --------
!        To pack model gridpoint fields by writing to, then reading from a
!         file ARPEGE
!        In DM version : Fields are distributed between the processors ; 
!        Each processor performs its own loop on fields to write & read back
!        on a local file.

!**   INTERFACE.
!     ----------
!       *CALL* *PKGRIDA*

!        EXPLICIT ARGUMENTS
!        --------------------
!        PGPBUF    - memory buffer (used in case of no real IO)

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 95-07-18

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-03-16 No I/O in packing
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      Apr 2008  K. Yessad: use DISGRID instead of DISGRID_C + cleanings
!      R. El Khatib : 01-Apr-2010 Overhead reduction
!      P.Marguinaud : 28-05-2010 Change SUMPIOH interface
!      P.Marguinaud : 11-09-2012 Initialize INBARI
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El khatib 16-May-2014 Optimization of in-line/off-line post-processing reproducibility
!      P.Marguinaud 04-Oct-2016 Port to single precision
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NSCRTCH
USE YOMOPH0      , ONLY : CNMCA
USE YOMMP0       , ONLY : NPROC    ,MYPROC
USE DIWRGRID_MOD , ONLY : DIWRGRID_SEND, DIWRGRID_RECV
USE DISGRID_MOD  , ONLY : DISGRID_SEND, DISGRID_RECV
USE FA_MOD   , ONLY : JPPRCM

IMPLICIT NONE

TYPE(GEOMETRY) , INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN) :: KBIT(KFIELDS)
REAL(KIND=JPRB), INTENT(INOUT) :: PGP(YDGEOMETRY%YRGEM%NGPTOT,KFIELDS)

!     ------------------------------------------------------------------

REAL(KIND=JPRB), ALLOCATABLE :: ZREALG (:,:)

REAL(KIND=JPRB), ALLOCATABLE :: ZVALCO (:)

INTEGER(KIND=JPIM) :: INFD(NPROC), IFLDOFF(NPROC)
INTEGER(KIND=JPIM) :: IEND, IFIELDS, IST, JROC, ICH
INTEGER(KIND=JPIM) :: JF, INBARI, INOMA, ILONGA, INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL,IBITGP
INTEGER(KIND=JPIM) :: IREP, IPFAOVSZ

CHARACTER(LEN=16) :: CLSCRTCH, CLNOMA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     IFIELDS  : Number of fields packed by MYPROC
!     ZREALG   : global array of fields

!     ------------------------------------------------------------------

#include "sumpioh.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PKGRIDA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG)
!     ------------------------------------------------------------------

!*       1.1 DISTRIBUTE FIELDS
!            -----------------

CALL SUMPIOH(NPROC,NPROC,KFIELDS,INFD,IFLDOFF)
ICH=1

!*       1.3 PACK DATA THROUGH WRITE AND READ
!            --------------------------------

DO JROC=1,NPROC
  IF (INFD(JROC) > 0 .AND. JROC /= MYPROC) THEN
    IST=IFLDOFF(JROC)+1
    IEND=IFLDOFF(JROC)+INFD(JROC)
    CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,JROC,INFD(JROC),PGP(:,IST:IEND),ICH)
  ENDIF
ENDDO

!*       1.3.2 Receive, pack and send back the fields I pack

IFIELDS=INFD(MYPROC)
IF (IFIELDS > 0) THEN
 
  CALL FASGRA (IREP, CNMCA, IPFAOVSZ)

  ALLOCATE (ZVALCO (YDGEM%NGPTOTG+JPPRCM*IPFAOVSZ))
  ALLOCATE (ZREALG (YDGEM%NGPTOTG,IFIELDS))

  IST=IFLDOFF(MYPROC)+1
  IEND=IFLDOFF(MYPROC)+IFIELDS
  CALL DIWRGRID_RECV(YDGEOMETRY,IFIELDS,PGP(:,IST:IEND),ICH,ZREALG)
  INBARI=0
  CALL FANOUV(IREP,NSCRTCH,.FALSE.,CLSCRTCH,'NEW',.TRUE.,.TRUE.,1,1,INBARI,CNMCA)  
  ! Get default value
  CALL FAVEUR(IREP,NSCRTCH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  DO JF=1,IFIELDS
    ! set number of bits for packing
    IF (KBIT(IFLDOFF(MYPROC)+JF) > 0) THEN
      IBITGP=KBIT(IFLDOFF(MYPROC)+JF)
    ELSE
      IBITGP=INBPDG
    ENDIF
    CALL FAGOTE(IREP,NSCRTCH,INGRIB,IBITGP,INBCSP,ISTRON,IPUILA,IDMOPL)
    ! pack/unpack
    CALL FACOND(IREP,NSCRTCH,'S',0,'FIELD',ZREALG(1,JF),.FALSE.,CLNOMA,INOMA,&
     & ZVALCO,ILONGA)  
    CALL FADECO(IREP,NSCRTCH,'S',0,'FIELD',.FALSE.,CLNOMA,INOMA,ZVALCO,&
     & ILONGA,ZREALG(1,JF))  
  ENDDO
  CALL FAIRNO(IREP,NSCRTCH,'DELETE')
  CALL DISGRID_SEND(YDGEOMETRY,IFIELDS,ZREALG,ICH,PGP(:,IST:IEND))

  DEALLOCATE (ZREALG)
  DEALLOCATE (ZVALCO)

ENDIF


!*       1.3.3 Receive fields I have not packed and include back all fields

DO JROC=1,NPROC
  IF (INFD(JROC) > 0 .AND. JROC /= MYPROC) THEN
    IST=IFLDOFF(JROC)+1
    IEND=IFLDOFF(JROC)+INFD(JROC)
    CALL DISGRID_RECV(YDGEOMETRY,JROC,INFD(JROC),PGP(:,IST:IEND),ICH)
  ENDIF
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PKGRIDA',1,ZHOOK_HANDLE)
END SUBROUTINE PKGRIDA

