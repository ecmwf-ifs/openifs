! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUMPFPOS_DEP(YDFPGEO_DEP)

!**** *SUMPFPOS_DEP*  - FULL-POS HORIZONTAL INTERPOLATIONS

!     PURPOSE.
!     --------
!       Statistics on the lengths of non-empty output rows

!**   INTERFACE.
!     ----------
!       *CALL* *SUMPFPOS_DEP*

!        EXPLICIT ARGUMENTS
!        ------------------
!        NONE.

!        IMPLICIT ARGUMENTS
!        --------------------
!        See modules below.

!     METHOD.
!     -------
!        SEE DOCUMENTATION ABOUT FULL-POS

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE* (SUFPSC2)

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 03-01-29 Remove SM aspects/cleanings
!      R. El Khatib : 03-04-17 Fullpos improvments
!      M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 23-may-2005 NFPBLOFF
!      M. Jidane    : 19-04-2006  Correction of a bug in allocation of NFPBLOFF
!      K. Yessad    : 27-Feb-2007 Rename into SUMPFPOS_DEP, adapt, clean.
!      R. El Khatib : 17-Sep-2007 Proper default value for NFPROMA_DEP if the default is negative
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2 (differentiation
!      of interpolation grid and output grid)
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 08-Aug-2017 recycled
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMFPGEO, ONLY : TFPGEO
USE YOMMP0   , ONLY : MYPROC   ,NPROC
USE YOMTAG   , ONLY : MTAGPART
USE MPL_MODULE, ONLY : MPL_BROADCAST

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE (TFPGEO),  INTENT(IN) :: YDFPGEO_DEP

INTEGER(KIND=JPIM) :: ILENROWS(NPROC)
! Number of active/idle proc. for pp. work
INTEGER(KIND=JPIM) :: IACTIVG(1), IIDLEG(1)
! Proc. of the smallest/biggest pp. row
INTEGER(KIND=JPIM) :: IPROCN(1), IPROCX(1)

! Size of the smallest/biggest pp. row 
INTEGER(KIND=JPIM) :: ISMALLG, IBIGG
! Mean size of active pp. rows
INTEGER(KIND=JPIM) :: IMEAN

INTEGER(KIND=JPIM) :: ITAG, JROC, IPROC

REAL(KIND=JPRB) :: ZFILL(NPROC) ! Buffer fill-up ratio
REAL(KIND=JPRB) :: ZWORK(NPROC) ! Working ratio
REAL(KIND=JPRB) :: ZNUME(NPROC) ! For more statistics : numerator
REAL(KIND=JPRB) :: ZDENO(NPROC) ! For more statistics : denominator

REAL(KIND=JPRB) :: ZRFPRGPG     ! 1/NFPRGPG_DEP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUMPFPOS_DEP',0,ZHOOK_HANDLE)
ASSOCIATE(NFPRGPL_DEP=>YDFPGEO_DEP%NFPRGPL, NFPRGPG_DEP=>YDFPGEO_DEP%NFPRGPG, &
 & NFPROMA_DEP=>YDFPGEO_DEP%NFPROMA)

!     ------------------------------------------------------------------


  ! Preinitialisation :
  ILENROWS(:)=-999

  ! Fill my own row length :
  ILENROWS(MYPROC)=NFPRGPL_DEP

  ! Broadcast my row length / gather the other rows lengths :
  SELECT CASE (NPROC)
  CASE (2:)
    DO JROC=1,NPROC
      ITAG=MTAGPART+NPROC+JROC
      IPROC=JROC
      CALL MPL_BROADCAST(ILENROWS(IPROC:IPROC),KROOT=IPROC,KTAG=ITAG, &
       & CDSTRING='SUMPFPOS_DEP')
    ENDDO
  CASE(1)
    CONTINUE
  CASE DEFAULT
    CALL ABOR1('SUMPFPOS_DEP:ILLEGAL VALUE NPROC')
  END SELECT   

  ! Number of active/idle processors :
  IIDLEG(1)=COUNT(ILENROWS==0)
  IACTIVG(1)=COUNT(ILENROWS>0)
  IF (IIDLEG(1)+IACTIVG(1) /= NPROC) THEN
    CALL ABOR1('SUMPFPOS_DEP:ILLEGAL VALUE ILENROWS')
  ENDIF

  ! Smallest and biggest rows :
  IPROCN=MINLOC(ILENROWS)
  IPROCX=MAXLOC(ILENROWS)
  ISMALLG=ILENROWS(IPROCN(1))
  IBIGG=ILENROWS(IPROCX(1))

  ! Mean value of active row lengths :
  IMEAN=NINT(REAL(NFPRGPG_DEP,JPRB)/REAL(IACTIVG(1),JPRB))

    WRITE(UNIT=NULOUT,FMT='(A)') &
     & ' SUMPFPOS_DEP : GLOBAL STATISTICS FOR ALL POST-PROCESSED POINTS :'  
    WRITE(UNIT=NULOUT,FMT='('' MINIMUM SIZE OF ROWS = '',I6, &
     & '' PROC # '',I4)') ISMALLG,IPROCN  
    WRITE(UNIT=NULOUT,FMT='('' MAXIMUM SIZE OF ROWS = '',I6, &
     & '' PROC # '',I4)') IBIGG,IPROCX  
    WRITE(UNIT=NULOUT,FMT='('' MEAN SIZE OF ROWS '',I6)') IMEAN
    WRITE(UNIT=NULOUT,FMT='('' NUMBER OF IDLE PROCESSORS = '',I3)') IIDLEG(1)

    ! More statistics : buffers fill-up ratio :
    ZRFPRGPG=1.0_JPRB/REAL(NFPRGPG_DEP,JPRB)
    DO JROC=1,NPROC
      ZNUME(JROC)=100._JPRB*REAL(ILENROWS(JROC),JPRB)
      ZWORK(JROC)=ZNUME(JROC)*ZRFPRGPG
      ZDENO(JROC)=REAL(((ILENROWS(JROC)-1)/NFPROMA_DEP +1)*NFPROMA_DEP,JPRB)
    ENDDO
    DO JROC=1,NPROC
      SELECT CASE (ILENROWS(JROC))
      CASE (1:)
        ZFILL(JROC)=ZNUME(JROC)/ZDENO(JROC)
      CASE DEFAULT
        ZFILL(JROC)=-999._JPRB
      END SELECT
    ENDDO
    WRITE(NULOUT,*) ' (PROCESSOR  LENGHT  WORKING RATIO(%)  FILL-UP RATIO(%) )'
    WRITE(NULOUT,'('' # '',4I8)') &
     & (JROC,ILENROWS(JROC),NINT(ZWORK(JROC)),NINT(ZFILL(JROC)),JROC=1,NPROC)  

! -------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUMPFPOS_DEP',1,ZHOOK_HANDLE)
END SUBROUTINE SUMPFPOS_DEP
