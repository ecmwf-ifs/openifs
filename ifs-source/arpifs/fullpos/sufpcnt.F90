! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPCNT(KFRFPOS,KFPOSTS,KFPOSTSMIN,CDIRLST,CDPATH,CDMONIPATH_IN,CDMONIPATH_OUT,CDFPNCF,KFPCONF,YDFPCNT,CDNAM)

!**** *SUFPCNT*  - INITIALIZE FULL POST PROCESSING (control options)

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPCNT

!        EXPLICIT ARGUMENTS     
!        --------------------
!           KFRFPOS,KFPOSTS,KFPOSTSMIN : control of compute frequency
!           CDIRLST  : basename of the file containing the list of pp request namelists
!           CDPATH  : object-relative dirname of the file containing the list of pp request namelists
!           CDMONIPATH_IN : dirname of the control files in input (allows to support multiple objects)
!           CDMONIPATH_OUT : dirname of the control files in output (allows to support multiple objects)
!           CDFPNCF : basename of the control file used to monitor output files achievement     
!           KFPCONF : configuration of the post-processing :
!                     0 : vertical interpolation only (<CFPFMT='MODEL'>)
!                     1 : gridpoint post-processing, possibly with spectral filters (<NFPOS=1>)
!                     2 : gridpoint/spectral post-processing (spectral outputs possible) (<NFPOS=2>)
!           CDNAM    : namelist file attached to this object

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
!      ORIGINAL : 06-Jun-2017 from sufpc

!     MODIFICATIONS.
!     --------------

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT  ,NPPPSH  ,NPDIRL
USE YOMCT0   , ONLY : LECMWF
USE YOMMP0   , ONLY : NPROC, MYPROC
USE YOMFPCNT , ONLY : TFPCNT
USE MPL_MODULE, ONLY : MPL_ALLGATHERV, MPL_BARRIER



!     ------------------------------------------------------------------

IMPLICIT NONE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), INTENT(IN) :: KFRFPOS
INTEGER(KIND=JPIM), INTENT(IN) :: KFPOSTS(0:)
INTEGER(KIND=JPIM), INTENT(IN) :: KFPOSTSMIN(0:)
CHARACTER(LEN=*),   INTENT(IN) :: CDIRLST
CHARACTER(LEN=*),   INTENT(IN) :: CDPATH
CHARACTER(LEN=*),   INTENT(IN) :: CDMONIPATH_IN
CHARACTER(LEN=*),   INTENT(IN) :: CDMONIPATH_OUT
CHARACTER(LEN=*),   INTENT(IN) :: CDFPNCF
INTEGER(KIND=JPIM), INTENT(IN) :: KFPCONF
TYPE(TFPCNT),       INTENT(OUT) :: YDFPCNT
CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: CDNAM

! CLIST : the actual name of the file containing the list of pp request namelists (it depends of the object)
! CLDIR : the actual name of the directory from which CLIST can be made
! CLCAT : the actual name of the file containing the concatenated lists of pp request namelists (it depends of the object)
CHARACTER(LEN=180) :: CLIST, CLDIR, CLCAT
INTEGER(KIND=JPIM) :: IEXIST, IIEXIST(NPROC), IERR, J
LOGICAL :: LLEXIST
CHARACTER(LEN=4)   :: CIERR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPCNT',0,ZHOOK_HANDLE)
ASSOCIATE(LFPCNT=>YDFPCNT%LFPCNT, LFPNAMELIST=>YDFPCNT%LFPNAMELIST, NFPCONF=>YDFPCNT%NFPCONF, CFPNCF=>YDFPCNT%CFPNCF)

! Make CLIST (or at least try to ...)
IF (CDPATH == ' ' .AND. CDIRLST == ' ') THEN
  CLIST=' '
ELSEIF(CDPATH == ' ' .AND. CDIRLST /= ' ') THEN
  IF (CDMONIPATH_IN /= ' ' .AND. CDMONIPATH_IN /= '.') THEN
    CLIST=TRIM(CDMONIPATH_IN)//'/'//CDIRLST
  ELSE
    CLIST=CDIRLST
  ENDIF
ELSEIF(CDPATH /= ' ' .AND. CDIRLST == ' ') THEN
  IF (CDMONIPATH_IN /= ' ' .AND. CDMONIPATH_IN /= '.') THEN
    CLIST=TRIM(CDMONIPATH_IN)//'/'//TRIM(CDPATH)//'dirlst'
  ELSE
    CLIST=TRIM(CDPATH)//'dirlst'
  ENDIF
ELSE
  IF (CDMONIPATH_IN /= ' ' .AND. CDMONIPATH_IN /= '.') THEN
    CLIST=TRIM(CDMONIPATH_IN)//'/'//TRIM(CDPATH)//TRIM(CDIRLST)
  ELSE
    CLIST=TRIM(CDPATH)//'/'//TRIM(CDIRLST)
  ENDIF
ENDIF
! Wait until all task gives the same answer for the existence of the file CLIST :
DO
  INQUIRE(FILE=CLIST,EXIST=LLEXIST)
  IF (LLEXIST) THEN
    IEXIST=1
  ELSE
    IEXIST=0
  ENDIF
  IF (NPROC > 1) THEN
    CALL MPL_ALLGATHERV(IEXIST,IIEXIST,CDSTRING='SUFPCNT:')
  ELSE
    IIEXIST(:)=IEXIST
  ENDIF
  IF ( ALL(IIEXIST==1) .OR. ALL(IIEXIST==0) ) THEN
    EXIT
  ENDIF
ENDDO
IF(.NOT.LLEXIST .AND. CDPATH /= ' ') THEN
  ! The file CLIST does not exist, make it from what is in the directory (objdir)/CDPATH 
  ! Keep CDMONIPATH_IN readonly. We may write on CDMONIPATH_OUT
  IF (CDMONIPATH_IN /= ' ' .AND. CDMONIPATH_IN /= '.') THEN
    CLDIR=TRIM(CDMONIPATH_IN)//'/'//TRIM(CDPATH)
  ELSE
    CLDIR=TRIM(CDPATH)
  ENDIF
  IF (CDIRLST == ' ') THEN
    IF (CDMONIPATH_OUT /= ' ' .AND. CDMONIPATH_OUT /= '.') THEN
      CLIST=TRIM(CDMONIPATH_OUT)//'/'//'dirlst'
    ELSE
      CLIST='dirlst'
    ENDIF
  ELSE
    IF (CDMONIPATH_OUT /= ' ' .AND. CDMONIPATH_OUT /= '.') THEN
      CLIST=TRIM(CDMONIPATH_OUT)//'/'//TRIM(CDIRLST)
    ELSE
      CLIST=TRIM(CDIRLST)
    ENDIF
  ENDIF
  IF (MYPROC==1) THEN
    CALL EXECUTE_COMMAND_LINE('/bin/ls '//TRIM(CLDIR)//'/* > '//TRIM(CLIST),EXITSTAT=IERR)
    IF (IERR /= 0) THEN
      WRITE(CIERR,'(A4)') IERR
      CALL ABOR1('SUFPCNT : FAILED IN EXECUTE_COMMAND_LINE /bin/ls EXITSTAT='//CIERR)
    ENDIF
  ENDIF
  CALL MPL_BARRIER(CDSTRING='SUFPCNT:')
! Wait until all task gives the same answer for the existence of the file CLIST :
  DO
    INQUIRE(FILE=CLIST,EXIST=LLEXIST)
    IF (LLEXIST) THEN
      IEXIST=1
    ELSE
      IEXIST=0
    ENDIF
    IF (NPROC > 1) THEN
      CALL MPL_ALLGATHERV(IEXIST,IIEXIST,CDSTRING='SUFPCNT:')
    ELSE
      IIEXIST(:)=IEXIST
    ENDIF
    IF ( ALL(IIEXIST==1) .OR. ALL(IIEXIST==0) ) THEN
      EXIT
    ENDIF
  ENDDO
  IF (.NOT.LLEXIST) THEN
    CALL ABOR1(' ERROR : SUFPCNT CANNOT MAKE CLIST')
  ENDIF
ENDIF


! NOTICE : the content of cdirlst could be stored in objects in memory

LFPCNT=LLEXIST
IF (LFPCNT) THEN
  WRITE(UNIT=NULOUT,FMT='('' LFPCNT = '',L2,'' CLIST = '',A)') LFPCNT, TRIM(CLIST)
  OPEN(NPDIRL,FILE=CLIST,ACTION='READ')
  IF (LECMWF) THEN
    IF (CDMONIPATH_IN /= ' ' .AND. CDMONIPATH_IN /= '.') THEN
      CLCAT=TRIM(CDMONIPATH_IN)//'/ppnamelist'
    ELSE
      CLCAT='ppnamelist'
    ENDIF
    INQUIRE(FILE=TRIM(CLCAT),EXIST=LFPNAMELIST)
    IF (LFPNAMELIST) THEN
      WRITE(UNIT=NULOUT,FMT='('' LFPNAMELIST = '',L2, '' CLCAT = '',A)') LFPNAMELIST, TRIM(CLCAT)
      OPEN(UNIT=NPPPSH,FILE=TRIM(CLCAT),ACTION='READ')
    ELSE
      WRITE(UNIT=NULOUT,FMT='('' LFPNAMELIST = '',L2)') LFPNAMELIST
    ENDIF
  ENDIF
ELSE
  WRITE(UNIT=NULOUT,FMT='('' LFPCNT = '',L2)') LFPCNT
  LFPNAMELIST=.FALSE.
ENDIF

NFPCONF=KFPCONF
IF (CDFPNCF /= ' ') THEN
  CFPNCF=TRIM(CDMONIPATH_OUT)//'/'//TRIM(CDFPNCF)
ELSE
  CFPNCF=' '
ENDIF
IF (PRESENT(CDNAM)) THEN
  YDFPCNT%CNAM=CDNAM
ENDIF
YDFPCNT%NFRFPOS=KFRFPOS
ALLOCATE(YDFPCNT%NFPOSTS(0:ABS(KFPOSTS(0))))
YDFPCNT%NFPOSTS(0:)=KFPOSTS(0:ABS(KFPOSTS(0)))
ALLOCATE(YDFPCNT%NFPOSTSMIN(0:ABS(KFPOSTSMIN(0))))
YDFPCNT%NFPOSTSMIN(0:)=KFPOSTSMIN(0:ABS(KFPOSTSMIN(0)))

WRITE(UNIT=NULOUT,FMT='('' NFPCONF = '', I2, '' CFPNCF = '',A, '' CNAM = '',A)') NFPCONF, TRIM(CFPNCF), TRIM(CDNAM)
WRITE(UNIT=NULOUT,FMT='('' NFPOSTS(0) = '',I5,'' NFRFPOS = '',I3)') YDFPCNT%NFPOSTS(0), YDFPCNT%NFRFPOS
IF (YDFPCNT%NFPOSTS(0) /= 0) THEN
  WRITE(UNIT=NULOUT,FMT='(6('' NFPOSTS('',I3,'') = '',I4))') (J,YDFPCNT%NFPOSTS(J),J=1,ABS(YDFPCNT%NFPOSTS(0)))
ENDIF
IF (YDFPCNT%NFPOSTSMIN(0) /= 0) THEN
  WRITE(UNIT=NULOUT,FMT='(5('' NFPOSTSMIN('',I3,'') = '',I4))') (J,YDFPCNT%NFPOSTSMIN(J),J=1,ABS(YDFPCNT%NFPOSTSMIN(0)))
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPCNT',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPCNT
