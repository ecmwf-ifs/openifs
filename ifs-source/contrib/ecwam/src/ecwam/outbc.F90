! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTBC (FL1, BLK2GLO, IU19)

! ----------------------------------------------------------------------

!**** *OUTBC* - OUTPUT OF THE COARSE GRID BOUNDARY VALUES.

!     R. PORTZ     MPI          JANUARY 1991
!     J. BIDLOT    ECMWF        FEBRARY 1996  MESSAGE PASSING

!*    PURPOSE.
!     --------

!        WRITE THE BOUNDARY VALUE OUTPUT FILE(s).

!**   INTERFACE.
!     ----------

!    *CALL* *OUTBC (FL1, BLK2GLO, IU19)*
!      *FL1*     - BLOCK OF SPECTRA.
!      *BLK2GLO* - BLOCK TO GRID TRANSFORMATION
!      *IU19*    - OUTPUT UNIT OF BOUNDARY VALUES.


!     METHOD.
!     -------

!       SEQUENCIAL UNFORMATED WRITE TO UNIT(s).

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LL1D
      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC    ,GBOUNC   ,IPOGBO
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : NPROC  , IRANK
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSPEC  , ONLY : IJ2NEWIJ

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "femean.intfb.h"
#include "mpgatherbc.intfb.h"
#include "sthq.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      TYPE(WVGRIDGLO), INTENT(IN)                         :: BLK2GLO
      INTEGER(KIND=JWIM),DIMENSION(GBOUNC), INTENT(IN) :: IU19


      INTEGER(KIND=JWIM), PARAMETER :: NSCFLD=3
      INTEGER(KIND=JWIM) :: IJ, II, NGOU, IX, KX, K ,M 
      INTEGER(KIND=JWIM) :: ICHNK 
      INTEGER(KIND=JWIM) :: IRECV

      REAL(KIND=JWRB) :: XLON, XLAT
      REAL(KIND=JWRB),DIMENSION(NPROMA_WAM, NCHNK) :: EM, FM, TQ
      REAL(KIND=JWRB),DIMENSION(NBOUNC) :: EMPTS, FMPTS, TQPTS
      REAL(KIND=JWRB),DIMENSION(NBOUNC, NANG, NFRE) :: FLPTS


! ----------------------------------------------------------------------



!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
      DO ICHNK = 1, NCHNK
        CALL FEMEAN (1, NPROMA_WAM, FL1(:,:,:, ICHNK), EM(:, ICHNK), FM(:, ICHNK) )
        CALL STHQ (1, NPROMA_WAM, FL1(:,:,:, ICHNK), TQ(:, ICHNK))
      ENDDO
!$OMP END PARALLEL DO

      IRECV=1
      CALL MPGATHERBC(IRECV, NSCFLD,                   &
     &                FL1, EM, TQ, FM,                 &
     &                FLPTS, EMPTS, TQPTS, FMPTS)


!*    2. WRITE INFORMATION TO FILE(s) IU19.
!        ----------------------------------

      IF (IRANK == IRECV) THEN
        DO II=1,GBOUNC
          DO NGOU=IPOGBO(II-1)+1,IPOGBO(II)
            IF (LL1D .OR. NPROC == 1) THEN
              IJ = IJARC(NGOU)
            ELSE
              IJ = IJ2NEWIJ(IJARC(NGOU))
            ENDIF
            IX = BLK2GLO%IXLG(IJ)
            KX = BLK2GLO%KXLT(IJ)
            XLON = AMOWEP + REAL(IX-1)*ZDELLO(KX)
            XLAT = AMOSOP + REAL(KX-1)*XDELLA

            WRITE(IU19(II)) XLON, XLAT, CDTPRO,                         &
     &          EMPTS(NGOU), TQPTS(NGOU), FMPTS(NGOU) 
            WRITE(IU19(II)) ((FLPTS(NGOU,K,M),K=1,NANG),M=1,NFRE)
          ENDDO
        ENDDO
      ENDIF

END SUBROUTINE OUTBC
