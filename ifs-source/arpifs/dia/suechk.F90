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

SUBROUTINE SUECHK(YDGEOMETRY,YDDIMF)

!     PURPOSE
!     -------
!       INITIALIZE CONSTANTS FOR EVOLUTION DIAGNOSTICS

!     INTERFACE
!     ---------
!       *CALL* *SUECHK*

!     EXPLICIT ARGUMENTS : None
!     ------------------

!     IMPLICIT ARGUMENTS : COMMON YEMCHK
!     ------------------

!     REFERENCE
!     ---------
!        ARPEGE/ALADIN documentation

!     AUTHOR
!     ------
!      Gabor Radnoti GMAP/MICECO
!      Original: 1992.12.24

!     MODIFICATIONS
!     -------------
!      C.Fisher, L. Gaytandjieva : 01-04-02 provisional initialization of nchktend
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT, NULNAM
USE YOMCT0       , ONLY : LECMWF
USE YOMDIMF      , ONLY : TDIMF
USE YOMCHK       , ONLY : JPGPCHK, JPFLDCHK, LECHKEVO, LECHKTND, LECHKPS,&
 &                        NXCHK, NYCHK, NNFCHK, NFRQCHK, NFLDCHK, NGPCHK, NLENCHK, NCHKTEND  

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDIMF)    ,INTENT(IN) :: YDDIMF
INTEGER(KIND=JPIM) :: J, JF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "posnam.intfb.h"

#include "namchk.nam.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUECHK',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDLON=>YDDIM%NDLON, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

!       1. SET UP DEFAULT VALUES FOR NAMCHK
!          --------------------------------

!        1.1 Set implicit default values

LECHKEVO= .FALSE.
LECHKTND= .FALSE.
LECHKPS = .FALSE.
NFRQCHK= 1
NFLDCHK= 1
NGPCHK= 0

! DEFAULT-FIELD: SURFACE PRESSURE
! !!! Caution: GFL and extra-GFL variables are currently not taken in
!     account in NNFCHK(1).
NNFCHK(1)= NFLEVG*(NFTHER+4)+1

DO JF=2,JPFLDCHK
  NNFCHK(JF)= -999
ENDDO
DO J=1,JPGPCHK
  NXCHK(J)= -999
  NYCHK(J)= -999
ENDDO

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
ENDIF

!     ------------------------------------------------------------------

!       2. READ AND CHECK NAMCHK
!          ---------------------

!        2.1 Read the namelist

CALL POSNAM(NULNAM,'NAMCHK')
READ(NULNAM,NAMCHK)

!        2.2 Check values

IF (LECHKEVO) THEN
  NFLDCHK= MIN(JPFLDCHK,MAX(1,NFLDCHK))
  NFRQCHK= MAX(1,ABS(NFRQCHK))

  ! !!! Caution: GFL and extra-GFL variables are currently not taken in
  !     account in array NNFCHK.
  DO JF=1,NFLDCHK
    NNFCHK(JF)= MIN(NFLEVG*(NFTHER+4)+1,MAX(1,NNFCHK(JF)))
  ENDDO

  NGPCHK= MIN(JPGPCHK,MAX(0,NGPCHK))
  IF (NGPCHK > 0) THEN
    DO J=1,NGPCHK
      NXCHK(J)= MIN(NDLON,MAX(1,NXCHK(J)))
      NYCHK(J)= MIN(NDGLG,MAX(1,NYCHK(J)))
    ENDDO
  ELSEIF (.NOT.LECHKTND) THEN
    WRITE(NULOUT,*) ' ! NO EVOLUTION DIAGNOSTICS REQUIRED !!'
    LECHKEVO=.FALSE.
  ENDIF
ENDIF

!        2.3 Compute additional values

! Provisional initialization of NCHKTEND
! NCHKTEND will be reset in chkevo

IF (LECHKTND) THEN
  NCHKTEND= NFLDCHK*NDLON*NDGLG
ELSE
  NCHKTEND= 1
ENDIF
NLENCHK = NFLDCHK*(NGPCHK+2)

!        2.4 Print

IF (LECHKEVO) THEN
  WRITE(NULOUT,*) '  '
  WRITE(NULOUT,*) ' **** EVOLUTION DIAGNOSTICS REQUIRED ****'
  WRITE(NULOUT,*) ' **** COMMON /YEMCHK/ ****'
  WRITE(NULOUT,*) ' NUMBER OF FIELDS IN DIAGNOSTICS : ',NFLDCHK
  WRITE(NULOUT,*) '  REQUIRED FIELDS : ',(NNFCHK(JF),JF=1,NFLDCHK)
  IF (LECHKPS)WRITE(NULOUT,*) '  CONVERSION ln(Ps) -> Ps WHENEVER REQUIRED'
  IF (LECHKTND)WRITE(NULOUT,*) ' GLOBAL DIAGNOSTICS FOR TENDENCIES REQUIRED'
  WRITE(NULOUT,*) ' FREQUENCY OF DIAGNOSTICS : ',NFRQCHK,' s'
  IF (NGPCHK > 0) THEN
    WRITE(NULOUT,*) ' LOCAL DIAGNOSTICS REQUIRED'
    WRITE(NULOUT,*) ' NUMBER OF GRIDPOINTS : ',NGPCHK
    WRITE(NULOUT,*) '   LATITUDE  INDEX : ',(NYCHK(J),J=1,NGPCHK)
    WRITE(NULOUT,*) '   LONGITUDE INDEX : ',(NXCHK(J),J=1,NGPCHK)
  ENDIF
  WRITE(NULOUT,*) '  '
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUECHK',1,ZHOOK_HANDLE)
END SUBROUTINE SUECHK
