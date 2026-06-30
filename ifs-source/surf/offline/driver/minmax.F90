SUBROUTINE MINMAX ( CNAME, PA, KLON, KLAT, LMASK, KULOUT )

#ifdef DOC
! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *MINMAX * - SPECIFYING STATISTICS OF ARRAY

!     Purpose.
!     --------
!            DISPLAY OF ARRAY STATISTICS

!**   Interface.
!     ----------
!        *CALL* *MINMAX*

!     Explicit arguments :
!     --------------------
!     CNAME	CHARACTER*8	FIELD NAME
!     PA	REAL()		ARRAY
!     KLON	INTEGER		X-DIMENSION
!     KLAT	INTEGER		Y-DIMENSION
!     LMASK	LOGICAL		MASK ARRAY
!     KULOUT	INTEGER		LOGICAL UNIT OUTPUT

!     Implicit arguments :
!     --------------------


!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd Hurk, KNMI

!     Modifications.
!     --------------
!        Original : 2000-07-07

!     ------------------------------------------------------------------

#endif
USE PARKIND1  ,ONLY : JPIM     ,JPRB      ,JPRD
IMPLICIT NONE

!*    DECLARATION OF GLOBAL PARAMETERS                                  
!     --------------------------------                                  

CHARACTER(LEN=*),INTENT(IN) :: CNAME
REAL(KIND=JPRB),INTENT(IN) :: PA(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KLON,KLAT,KULOUT
LOGICAL,INTENT(IN) :: LMASK(:)

!*      local variables
!       ---------------
INTEGER(KIND=JPIM) :: ILNLT,NMSK,IMINL,IMAXL,JL,IYMIN,IYMAX,IXMIN,IXMAX
REAL(KIND=JPRB) :: ZSUM,ZMIN,ZMAX

!---------------------------------------------------------------------- 


ILNLT = KLON*KLAT
NMSK=COUNT(LMASK)
IF(NMSK == 0)RETURN

ZSUM=0.
ZMIN=1.e30
ZMAX=-1.E30
IMINL=0
IMAXL=0
DO JL=1,ILNLT
  IF(LMASK(JL))THEN
    ZSUM = ZSUM+PA(JL)
    IF(PA(JL) < ZMIN)THEN
      ZMIN=PA(JL)
      IMINL=JL
    ENDIF
    IF(PA(JL) > ZMAX)THEN
      ZMAX=PA(JL)
      IMAXL=JL
    ENDIF
  ENDIF
ENDDO
ZSUM = ZSUM/REAL(NMSK,KIND=JPRB)

IYMIN = (IMINL-1)/KLON + 1
IXMIN =  IMINL-(IYMIN-1)*KLON
IYMAX = (IMAXL-1)/KLON + 1
IXMAX =  IMAXL-(IYMAX-1)*KLON

IF ( CNAME == 'LATITUDE' .OR. CNAME == 'LONGITUD' ) THEN
!to convert the lat/lon written value to deg
WRITE(KULOUT,'(1X,/,1X,A8)') CNAME
WRITE(KULOUT,'(1X,''MEAN-VALUE: '',G20.8,'' NUMBER: '',I10)')ZSUM*57.3248_JPRB,NMSK
WRITE(KULOUT,'(1X,''MIN -VALUE: '',G20.8,'' AT (X,Y) = '',2I5)')&
      &ZMIN*57.3248_JPRB,IXMIN,IYMIN                                            
WRITE(KULOUT,'(1X,''MAX -VALUE: '',G20.8,'' AT (X,Y) = '',2I5)')&
      &ZMAX*57.3248_JPRB,IXMAX,IYMAX                                            
ELSE

WRITE(KULOUT,'(1X,/,1X,A8)') CNAME
WRITE(KULOUT,'(1X,''MEAN-VALUE: '',G20.8,'' NUMBER: '',I10)')ZSUM,NMSK
WRITE(KULOUT,'(1X,''MIN -VALUE: '',G20.8,'' AT (X,Y) = '',2I5)')&
      &ZMIN,IXMIN,IYMIN                                            
WRITE(KULOUT,'(1X,''MAX -VALUE: '',G20.8,'' AT (X,Y) = '',2I5)')&
      &ZMAX,IXMAX,IYMAX                                            
ENDIF

RETURN
END SUBROUTINE
