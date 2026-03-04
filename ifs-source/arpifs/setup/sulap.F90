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

SUBROUTINE SULAP(YDLAP,YDDIM)

!**** *SULAP * - Routine to initialize the Laplace space working common

!     Purpose.
!     --------
!           Initialize COMMON YOMLAP , initialize pointer arrays.

!**   Interface.
!     ----------
!        *CALL* *SULAP *

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        COMMON YOMLAP, corresponding arrays of the area 1 of the memory manager

!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
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
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      R. El Khatib : 04-08-05 use YEMDIM instead of calling ellips
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      0. Spaniel 15-Nov-2007: cleaning
!      R. El Khatib 27-Sep-2013 Cleaning for LAM
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMLAP   , ONLY : TLAP
USE YOMDIM   , ONLY : TDIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LALLOPR, LELAM
USE YOMCST   , ONLY : RA
USE YOMMP0   , ONLY : NPRINTLEV, LOUTPUT, NPRTRW, MYSETW
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TLAP) , INTENT(INOUT) :: YDLAP
TYPE(TDIM) , INTENT(IN)    :: YDDIM
INTEGER(KIND=JPIM), PARAMETER :: JPKD=JPRD

INTEGER(KIND=JPIM) :: IC,  IM, INSM, IU,  JM, JMLOC, JN, JJ  

LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "trans_inq.h"
#include "ini_spec_dist.h"

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SULAP',0,ZHOOK_HANDLE)
ASSOCIATE(NCMAX=>YDDIM%NCMAX, NRESOL=>YDDIM%NRESOL, NSMAX=>YDDIM%NSMAX, &
 & NSPEC=>YDDIM%NSPEC, NSPEC2=>YDDIM%NSPEC2, NUMCP=>YDDIM%NUMCP, &
 & NUMP=>YDDIM%NUMP)
!     ------------------------------------------------------------------

!*       1.    Initialize Laplace space positioning.
!              -------------------------------------

IU = NULOUT
LLP = NPRINTLEV >= 1.OR. LALLOPR

ALLOCATE(YDLAP%MYMS(NUMCP))
IF(LLP)WRITE(IU,9) 'YDLAP%MYMS     ',SIZE(YDLAP%MYMS     ),SHAPE(YDLAP%MYMS     )
ALLOCATE(YDLAP%RLAPDI(-1:NSMAX+2))
IF(LLP)WRITE(IU,9) 'YDLAP%RLAPDI   ',SIZE(YDLAP%RLAPDI   ),SHAPE(YDLAP%RLAPDI   )
ALLOCATE(YDLAP%RLAPIN(-1:NSMAX+2))
IF(LLP)WRITE(IU,9) 'YDLAP%RLAPIN   ',SIZE(YDLAP%RLAPIN   ),SHAPE(YDLAP%RLAPIN   )
ALLOCATE(YDLAP%NASN0 (0:NSMAX))
IF(LLP)WRITE(IU,9) 'YDLAP%NASN0    ',SIZE(YDLAP%NASN0    ),SHAPE(YDLAP%NASN0    )
ALLOCATE(YDLAP%NASM0 (0:NSMAX))
IF(LLP)WRITE(IU,9) 'YDLAP%NASM0    ',SIZE(YDLAP%NASM0    ),SHAPE(YDLAP%NASM0    )
ALLOCATE(YDLAP%NASM0G(0:NSMAX))
IF(LLP)WRITE(IU,9) 'YDLAP%NASM0G   ',SIZE(YDLAP%NASM0G   ),SHAPE(YDLAP%NASM0G   )
ALLOCATE(YDLAP%NSPZERO(0:NSMAX))
IF(LLP)WRITE(IU,9) 'YDLAP%NSPZERO  ',SIZE(YDLAP%NSPZERO  ),SHAPE(YDLAP%NSPZERO  )
ALLOCATE(YDLAP%NVALUE (NSPEC2))
IF(LLP)WRITE(IU,9) 'YDLAP%NVALUE   ',SIZE(YDLAP%NVALUE   ),SHAPE(YDLAP%NVALUE   )

IF (LELAM) THEN
  ! nasm0 is not used in LAM models
ELSE
  CALL TRANS_INQ(KRESOL=NRESOL,KASM0=YDLAP%NASM0)
ENDIF
IF(LELAM .OR. NCMAX==NSMAX) THEN
  CALL TRANS_INQ(KRESOL=NRESOL,KMYMS=YDLAP%MYMS)
ELSE
  CALL INI_SPEC_DIST(NCMAX,NCMAX,NPRTRW,MYSETW,KMYMS=YDLAP%MYMS)
ENDIF


DO JN=0,NSMAX
  YDLAP%NASN0(JN)=JN*JN+JN+1
ENDDO

IF (.NOT.LELAM) THEN
  ALLOCATE(YDLAP%NSE0L(NUMP))
  IF(LLP)WRITE(IU,9) 'YDLAP%NSE0L    ',SIZE(YDLAP%NSE0L    ),SHAPE(YDLAP%NSE0L    )
  IC=0
  DO JMLOC=1,NUMP
    IM=YDLAP%MYMS(JMLOC)
    YDLAP%NSE0L(JMLOC)=IC
    IC=IC+(NSMAX+1-IM)
  ENDDO
  IF (IC /= NSPEC) CALL ABOR1('SULAP: IC /= NSPEC')
ENDIF


IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  WRITE(NULOUT,'('' YDLAP%MYMS '')')
  WRITE(NULOUT,'(30(1X,I3))')(YDLAP%MYMS(JJ),JJ=1,NUMP)
  IF (.NOT.LELAM) THEN
    WRITE(NULOUT,'('' NASM0 '')')
    WRITE(NULOUT,'(16(1X,I6))')(YDLAP%NASM0(JJ),JJ=0,NSMAX)
  ENDIF
  WRITE(UNIT=NULOUT,FMT='('' YDLAP%NASN0'')')
  WRITE(UNIT=NULOUT,FMT='(16(1X,I6))')YDLAP%NASN0
  IF (.NOT.LELAM) THEN
    WRITE(UNIT=NULOUT,FMT='('' YDLAP%NSE0L '')')
    WRITE(UNIT=NULOUT,FMT='(16(1X,I6))')YDLAP%NSE0L
  ENDIF
ENDIF

!      -----------------------------------------------------------------

!*       2.    Initialize Laplacian operator and its inverse.
!              ----------------------------------------------

DO JN=1,NSMAX+2
  YDLAP%RLAPDI(JN)=REAL(-REAL(JN*(JN+1),JPKD)/(REAL(RA,JPKD)*REAL(RA,JPKD)),JPRB)
  YDLAP%RLAPIN(JN)=REAL(-(REAL(RA,JPKD)*REAL(RA,JPKD))/REAL(JN*(JN+1),JPKD),JPRB)
ENDDO
YDLAP%RLAPDI(0)=0._JPRB
YDLAP%RLAPIN(0)=0._JPRB
YDLAP%RLAPDI(-1)=0.0_JPRB
YDLAP%RLAPIN(-1)=0.0_JPRB

IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  WRITE(UNIT=NULOUT,FMT='('' EIGEN-VALUES OF THE LAPLACIAN'')')
  WRITE(UNIT=NULOUT,FMT='('' TRANSFORMED SPHERE'')')
  WRITE(UNIT=NULOUT,FMT='(14(1X,E8.2))')YDLAP%RLAPDI
  WRITE(UNIT=NULOUT,FMT='('' EIGEN-VALUES OF ITS INVERSE'')')
  WRITE(UNIT=NULOUT,FMT='('' TRANSFORMED SPHERE'')')
  WRITE(UNIT=NULOUT,FMT='(14(1X,E8.2))')YDLAP%RLAPIN
ENDIF

!      -----------------------------------------------------------------

!*       3.    Initialize arrays used in the Legendre transforms.
!              --------------------------------------------------

IF (.NOT.LELAM) THEN
  IC=1
  DO JMLOC=1,NUMP
    IM=YDLAP%MYMS(JMLOC)
    DO JN=IM,NSMAX
      YDLAP%NVALUE(IC  )=JN
      YDLAP%NVALUE(IC+1)=JN
      IC=IC+2
    ENDDO
  ENDDO

  !      -----------------------------------------------------------------

!*       4.    Spectral control variable and state variable relationship
!              ---------------------------------------------------------

! Global version is needed for control space calculations

  INSM=1
  DO JM=0,NSMAX
    YDLAP%NASM0G(JM)=INSM
    INSM=INSM+(NSMAX+1-JM)*2
  ENDDO

  IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
    WRITE(UNIT=NULOUT,FMT='('' YDLAP%NASM0G '')')
    WRITE(UNIT=NULOUT,FMT='(14(1X,I6))')YDLAP%NASM0G
  ENDIF


!  Define imaginary m=0 positions.

  DO JN=0,NSMAX
    YDLAP%NSPZERO(JN)=1+JN*2+1
  ENDDO

ENDIF

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)


!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SULAP',1,ZHOOK_HANDLE)
END SUBROUTINE SULAP
