! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTSPEC (FL1, FF_NOW)

!----------------------------------------------------------------------

!**** *OUTSPEC*  ENCODES SPECTRA AS PARAMETER 251 USING GRIB API
!                AND WRITES TO FILE OR TO FDB.

!     J. BIDLOT   ECMWF  APRIL 2010 

!*    PURPOSE.
!     --------

!       ENCODES SPECTRA AND WRITES TO FDB OR TO A FILE.

!**   INTERFACE.
!     ----------

!     SUBROUTINE OUTSPEC (FL1, FF_NOW)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *FL1*       REAL      LOCAL SPECTRA OF CURRENT PE.
!      *FF_NOW*    REAL      FOR SEA ICE COVER.

!     METHOD.
!     -------

!      ENCODE SPECTRA PER FREQUENCY AND DIRECTION
!      INTO GRIB AND WRITE TO FDB OR TO A SINGLE FILE.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCOUP  , ONLY : LWCOU, LIFS_IO_SERV_ENABLED,                &
     &                      OUTWSPEC_IO_SERV
      USE YOWCOUT  , ONLY : LWAM_USE_IO_SERV
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWICE   , ONLY : LICERUN  ,CITHRSH  ,FLMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDATEA   ,    &
     &                      MARSTYPE ,LLSOURCE

      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
      USE YOWABORT, ONLY : WAM_ABORT

!-----------------------------------------------------------------------

      IMPLICIT NONE

#include "outwspec.intfb.h"
#include "difdate.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: IJSG, IJLG, IJSB, IJLB, KIJS, KIJL, ICHNK
      INTEGER(KIND=JWIM) :: IFCST

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM) :: ZMASK
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:, :, :) :: SPEC

      CHARACTER(LEN=14) :: CDATE, CDATED 

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('OUTSPEC',0,ZHOOK_HANDLE)

      CALL GSTATS(2080,0)

      IJSG = IJFROMCHNK(1,1)
      IJLG = IJSG + SUM(KIJL4CHNK) - 1
      ALLOCATE(SPEC(IJSG:IJLG, NANG, NFRE))

!*    APPLY SEA ICE MASK TO THE OUTPUT SPECTRA (IF NEEDED)
      IF (LICERUN .AND. LLSOURCE) THEN
!$OMP    PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, KIJS, KIJL, IJSB, IJLB, ZMASK, K, M)
         DO ICHNK = 1, NCHNK

           DO IJ = 1, KIJL4CHNK(ICHNK)
             IF (FF_NOW%CICOVER(IJ,ICHNK) > CITHRSH) THEN
               ZMASK(IJ) = 1.0_JWRB
             ELSE
               ZMASK(IJ) = 0.0_JWRB
             ENDIF
           ENDDO

           KIJS = 1
           IJSB = IJFROMCHNK(KIJS, ICHNK)
           KIJL = KIJL4CHNK(ICHNK)
           IJLB = IJFROMCHNK(KIJL, ICHNK)

           DO M = 1, NFRE
             DO K = 1, NANG
               SPEC(IJSB:IJLB,K,M) = (1.0_JWRB - ZMASK(KIJS:KIJL)) * FL1(KIJS:KIJL,K,M,ICHNK) + ZMASK(KIJS:KIJL)*FLMIN
             ENDDO
           ENDDO
         ENDDO
!$OMP    END PARALLEL DO

      ELSE

!$OMP    PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IJSB, IJLB)
         DO ICHNK = 1, NCHNK
           KIJS = 1
           IJSB = IJFROMCHNK(KIJS, ICHNK)
           KIJL = KIJL4CHNK(ICHNK)
           IJLB = IJFROMCHNK(KIJL, ICHNK)

           SPEC(IJSB:IJLB, :, :) = FL1(KIJS:KIJL, :, :, ICHNK)
         ENDDO
!$OMP    END PARALLEL DO
      ENDIF

      IF (CDTPRO <= CDATEF) THEN
!*    0.1.  THIS IS AN ANALYSIS DATE.
        IF (LWCOU .AND. MARSTYPE == 'fg') THEN
          CDATE=CDATEA
          CDATED=CDATEA
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSEIF (LWCOU .AND. MARSTYPE == '4v') THEN
          CDATE=CDATEA
          CDATED=CDATEA
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSE
          CDATE=CDTPRO
          CDATED=CDTPRO
          IFCST = 0
        ENDIF
      ELSE
!*    0.2.  THIS IS A  FORECAST DATE.
        CDATE=CDATEF
        CDATED=CDTPRO
        CALL DIFDATE (CDATEF, CDTPRO, IFCST)
        IFCST = IFCST/3600
      ENDIF
!-----------------------------------------------------------------------

!*    3. OUTPUT GRIB DATA
!     -------------------
      
      ! Use IFS IO server?
      IF (LWCOU .AND. LIFS_IO_SERV_ENABLED .AND. LWAM_USE_IO_SERV) THEN
        CALL OUTWSPEC_IO_SERV(IJSG, IJLG, SPEC, MARSTYPE, CDATE, IFCST)
      ELSE
        CALL OUTWSPEC(IJSG, IJLG, SPEC, MARSTYPE, CDATE, CDATED, IFCST)
      ENDIF

      CALL GSTATS(2080,1)

      DEALLOCATE(SPEC)

IF (LHOOK) CALL DR_HOOK('OUTWSPEC',1,ZHOOK_HANDLE)

END SUBROUTINE OUTSPEC
