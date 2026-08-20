! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WGRIBENCODE_MODEL (IU06, ITEST, I1, I2, FIELD,   &
 &                            ITABLE, IPARAM, KLEV, IK, IM, &
 &                            CDATE, IFCST, MARSTYPE,       &
 &                            IGRIB_HANDLE)

! ----------------------------------------------------------------------

!****  *WGRIBENCODE_MODEL*  ENCODES WAM MODEL FIELD INTO GRIB CODE AND OUTPUT

!       J. BIDLOT    ECMWF JULY 2009: USE GRIB API 

!       PURPOSE.
!       --------
!         SUBROUTINE PACKS WAVE FIELDS INTO THE GRIB CODE
!         !!!! PRESET_WGRIB_TEMPLATE NEEDS to BE CALLED BEFORE TO
!         !!!! INITIALISE THE GRIB TEMPLATES.

!**    INTERFACE.
!      ----------
!        *CALL* *WGRIBENCODE_MODEL (IU06, ITEST, I1, I2, FIELD,
!                             ITABLE, IPARAM, KLEV, IK, IM,
!                             CDATE, IFCST, MARSTYPE,
!                             IGRIB_HANDLE)
!          *IU06*    LOGFILE OUTPUT UNIT.
!          *ITEST*   TEST OUTPUT GIVEN IF ITEST GT 2.
!          *I1*      FIRST DIMENSION OF FIELD.
!          *I2*      SECOND DIMENSION OF FIELD.
!          *FIELD*   FIELD TO BE PACKED.
!          *ITABLE*  GRIB TABLE NUMBER.
!          *IPARAM*  PARAMETER IDENTIFIER.
!          *KLEV*    REFERENCE LEVEL IN full METER
!                    (SHOULD BE 0 EXCEPT FOR 233 AND 245)
!          *IK*      DIRECTION INDEX,
!                    ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *IM*      FREQUENCY INDEX,
!                    ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *CDATE*   DATE AND TIME.
!          *IFCST*   FORECAST STEP IN HOURS.
!          *MARSTYPE* TYPE OF CURRENT FIELD
!          *IGRIB_HANDLE  GRIB HANDLE CONTAINING THE DATA.

!      METHOD.
!      -------

!      EXTERNALS.
!      ----------

!      REFERENCES.
!      -----------

! ----------------------------------------------------------------------

 USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

 USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_I2, NGRIB_HANDLE_WAM_S         ! To clone
 USE YOWGRIBHD, ONLY : PPMISS   , PPEPS    ,PPREC    ,NTENCODE ,  &
            NGRBRESS ,PPRESOL  ,LGRHDIFS ,LNEWLVTP  ,PPMIN_RESET, &
            LPADPOLES, NDATE_TIME_WINDOW_END,NWINOFF
 USE YOWMAP   , ONLY : IRGG, AMONOP, AMOSOP, XDELLA, NLONRGG, CLDOMAIN
 USE YOWCOUP  , ONLY : KCOUSTEP
 USE YOWCOUT  , ONLY : LRSTST0
 USE YOWPCONS , ONLY : ZMISS
 USE YOWGRIB  , ONLY : IGRIB_CLONE
 USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
 IMPLICIT NONE
#include "preset_wgrib_template.intfb.h"
#include "wgrib_edition.intfb.h"
#include "wgribencode.intfb.h"

 INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ITEST, I1, I2
 INTEGER(KIND=JWIM), INTENT(IN) :: ITABLE, IPARAM, KLEV, IK, IM, IFCST
 INTEGER(KIND=JWIM), INTENT(OUT) :: IGRIB_HANDLE
 CHARACTER(LEN=2), INTENT(IN) :: MARSTYPE
 CHARACTER(LEN=14), INTENT(IN) :: CDATE
 REAL(KIND=JWRB), INTENT(INOUT) :: FIELD(I1,I2)


 INTEGER(KIND=JWIM) :: IPARAMID
 REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_MODEL',0,ZHOOK_HANDLE)

CALL GSTATS(1709,0)
IF (ITEST > 0) THEN
   WRITE(IU06,*) '   SUB. WGRIBENCODE_MODEL CALLED FOR ',IPARAM
   CALL FLUSH(IU06)
ENDIF

!     CLONE GRIB TEMPLATE:

IPARAMID=1000*ITABLE+IPARAM

IF (IPARAM == 251) THEN
   IF (LGRHDIFS) THEN
     CALL PRESET_WGRIB_TEMPLATE("S",IGRIB_HANDLE)
   ELSE
     IGRIB_HANDLE=-99
     CALL IGRIB_CLONE(NGRIB_HANDLE_WAM_S,IGRIB_HANDLE)
   ENDIF
ELSE
   IF (LGRHDIFS) THEN
     SELECT CASE(WGRIB_EDITION(IPARAMID))
     CASE (2) 
       CALL PRESET_WGRIB_TEMPLATE("I",IGRIB_HANDLE, NGRIBV=2)
     CASE DEFAULT
       CALL PRESET_WGRIB_TEMPLATE("I",IGRIB_HANDLE)
     END SELECT
   ELSE
     IGRIB_HANDLE=-99
     SELECT CASE(WGRIB_EDITION(IPARAMID))
     CASE (2) 
       CALL IGRIB_CLONE(NGRIB_HANDLE_WAM_I2,IGRIB_HANDLE)
     CASE DEFAULT
       CALL IGRIB_CLONE(NGRIB_HANDLE_WAM_I,IGRIB_HANDLE)
     END SELECT
   ENDIF
ENDIF

CALL WGRIBENCODE( IU06, ITEST, &
 &                I1, I2, &
 &                FIELD, &
 &                ITABLE, IPARAM, &
 &                KLEV, &
 &                IK, IM, &
 &                CDATE, IFCST, MARSTYPE, &
 &                PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, &
 &                LGRHDIFS, &
 &                NDATE_TIME_WINDOW_END, &
 &                NGRBRESS, LNEWLVTP, LPADPOLES, &
 &                SIZE(NLONRGG), NLONRGG, IRGG, &
 &                AMONOP, AMOSOP, XDELLA, CLDOMAIN, &
 &                KCOUSTEP, LRSTST0, &
 &                ZMISS, &
 &                IGRIB_HANDLE)

CALL GSTATS(1709,1)
IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE WGRIBENCODE_MODEL

