! (C) Copyright 1989- Meteo-France.

SUBROUTINE RDECCLIMO(YDAFN,CDFILE,KMONTH,YDFPUSERGEO,PFIELD,KGP,KFIELDS,KCOD)

!**** *RDECCLIMO*  - READ OUTPUT CLIMATOLOGY

!     PURPOSE.
!     --------
!        Open output climatology files and read the needed fields

!**   INTERFACE.
!     ----------
!       *CALL* *RDECCLIMO*

!        EXPLICIT ARGUMENTS
!        --------------------
!        CDFILE : filename
!        KMONTH : month
!        PFIELD : OUTPUT FIELDS ARRAY
!        KGP    : NUMBER OF OUTPUT POINTS
!        KFIELDS: NUMBER OF OUTPUT FIELDS
!        KCOD   : FIELDS CODES

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
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      Nils Wedi
!      ORIGINAL : 95-12-15

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 05-03-15 Cleanings
!      T. Wilhelmsson 11-01-24 Switch to GRIB_API
!      R. El Khatib : 01-Aug-2013 Remove LDMASK
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 22-mar-2016 flexible climatology filename with SUFPCLIFNAME
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMAFN, ONLY : TAFN
USE GRIB_API_INTERFACE, ONLY : IGRIB_OPEN_FILE, JPGRIB_SUCCESS, JPGRIB_END_OF_FILE, &
 & IGRIB_NEW_FROM_FILE, IGRIB_GET_VALUE, IGRIB_RELEASE, IGRIB_CLOSE_FILE, &
 & IGRIB_ERROR_MESSAGE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TAFN),        INTENT(IN) :: YDAFN
CHARACTER(LEN=*),  INTENT(IN) :: CDFILE
INTEGER(KIND=JPIM),INTENT(IN) :: KMONTH
TYPE (TFPUSERGEO) ,INTENT(IN) :: YDFPUSERGEO
INTEGER(KIND=JPIM),INTENT(IN) :: KGP 
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PFIELD(KGP,KFIELDS) 
INTEGER(KIND=JPIM),INTENT(IN) :: KCOD(KFIELDS) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZREAL(KGP)
INTEGER(KIND=JPIM) :: ILONS(YDFPUSERGEO%NLAT)
INTEGER(KIND=JPIM) :: IFLD, IGRIB
INTEGER(KIND=JPIM) :: IPARAM, IRET, IUNITGG, J, JFLD, JGL, JJ  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RDECCLIMO',0,ZHOOK_HANDLE)
ASSOCIATE(GFP=>YDAFN%GFP, GFP_PHYDS=>YDAFN%GFP_PHYDS, NLAT=>YDFPUSERGEO%NLAT)
!     ------------------------------------------------------------------

!*       2. LOOP ON SUBDOMAINS
!           ------------------


!*       2.1 OPEN FILE

  CALL IGRIB_OPEN_FILE(IUNITGG,CDFILE,'r')
  WRITE(NULOUT,*) 'READ FROM CLIMATE FILE ',CDFILE

!*       2.2 READ FIELDS

  FIELD_LOOP: DO JFLD=1,KFIELDS

    CALL IGRIB_NEW_FROM_FILE(IUNITGG,IGRIB,IRET)
    IF(IRET /= JPGRIB_SUCCESS) THEN
      CALL IGRIB_ERROR_MESSAGE(IRET)
      CALL ABOR1('RDECCLIMO: PROBLEM IN IGRIB_NEW_FROM_FILE')
    ELSEIF(IRET == JPGRIB_END_OF_FILE) THEN
      EXIT FIELD_LOOP
    ENDIF
    CALL GSTATS(1703,0)
    CALL IGRIB_GET_VALUE(IGRIB,'values',ZREAL)
    CALL GSTATS(1703,1)

    CALL IGRIB_GET_VALUE(IGRIB,'paramId',IPARAM)
    DO JJ=1,KFIELDS
      IF( GFP_PHYDS(KCOD(JJ))%IGRIB == IPARAM) IFLD=JJ
    ENDDO

    IF (IPARAM /= GFP%LSM%IGRIB .AND. IPARAM /= GFP%GFIS%IGRIB) THEN
      WRITE(NULOUT,*) 'CLIMATE FIELD NOT USED ... ',IPARAM
      CYCLE
    ELSE
      WRITE(NULOUT,*) 'PARAMETER ',IPARAM,' READ FROM CLIMATE ','FILE'
    ENDIF

!   *     2.3    CHECK PARAMETER.

    IF(YDFPUSERGEO%NFPHTYP /= 0) THEN
      CALL IGRIB_GET_VALUE(IGRIB,'pl',ILONS)
      DO JGL=1,NLAT
        IF(YDFPUSERGEO%NFPRGRI(JGL) /= ILONS(JGL)) THEN
          WRITE(NULOUT,*) ' INCONSISTENT REDUCED GRID'
          WRITE(NULOUT,*) ' IN MODEL ',(YDFPUSERGEO%NFPRGRI(J),J=1,NLAT)
          WRITE(NULOUT,*) ' IN FILE  ',(ILONS(J),J=1,NLAT)
          CALL ABOR1('RDECCLIMO')
        ENDIF
      ENDDO
    ENDIF

    IF (IPARAM == GFP%LSM%IGRIB) THEN
      DO J=1,KGP
        ZREAL(J)=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZREAL(J)-0.5_JPRB))
      ENDDO
    ENDIF

    PFIELD(:,IFLD)=ZREAL(:)

  ENDDO FIELD_LOOP

!*       2.4 CLOSE FILE

  CALL IGRIB_RELEASE(IGRIB)
  CALL IGRIB_CLOSE_FILE(IUNITGG)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RDECCLIMO',1,ZHOOK_HANDLE)
END SUBROUTINE RDECCLIMO
