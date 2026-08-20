! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTWPSP (FL1, FF_NOW)
! ----------------------------------------------------------------------

!**** *OUTWPSP* - MODEL OUTPUT OF SPECTRA AT GIVEN LOCATIONS 

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTWPSP (FL1, FF_NOW) 
!      *FL1*    - INPUT SPECTRUM.
!      *FF_NOW* - FORCING FIELDS

!     EXTERNALS.
!     ----------

!       *OUT_ONEGRDPT_SP* - ONE GRID POINT OUTPUT
!       *OUTERS*    - OUTPUT OF SATELLITE COLOCATION SPECTRA.
!   
!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK
      USE YOWMAP   , ONLY : CLDOMAIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDTPRO

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "out_onegrdpt_sp.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      TYPE(FORCING_FIELDS), INTENT(IN) :: FF_NOW


      INTEGER(KIND=JWIM) :: ICHNK, IJ
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTWPSP',0,ZHOOK_HANDLE)


!*    OUTPUT OF SPECTRA FOR ONE GRID POINT SIMULATION
!     -----------------------------------------------

      IF (CLDOMAIN == 's') THEN
        DO ICHNK = 1, NCHNK
          DO IJ = 1, KIJL4CHNK(ICHNK)
            CALL OUT_ONEGRDPT_SP(FL1(IJ:IJ,:,:,ICHNK), FF_NOW%UFRIC(IJ,ICHNK), CDTPRO)
          ENDDO
        ENDDO
      ENDIF

IF (LHOOK) CALL DR_HOOK('OUTWPSP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTWPSP
