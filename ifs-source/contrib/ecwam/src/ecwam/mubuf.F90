! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MUBUF (IU01, IU08, NPROPAGS)
! ----------------------------------------------------------------------

!**** *MUBUF* - ROUTINE TO ARRANGE YOWMON UBUF FOR ONE BLOCK.

!     H.GUNTHER            ECMWF       04/04/1990
!     J. BIDLOT            ECMWF       APRIL 2000: add second closest
!                                                  grid points.
!     J. BIDLOT            ECMWF       CLOSEST AND SECOND CLOSEST GRID
!                                      POINT FOR THE ROTATED CELL.
!                               IN ORDER TO SAVE MEMORY THE OBSTRUCTION
!                               COEFFICIENTS ARE READ IN AND PROCESSED
!                               SEQUENTIALLY.
!                          ECMWF       MODIFIED TO COMPUTE AND WRITE ALL
!                                      ARRAYS SEQUENTIALLY.

!*    PURPOSE.
!     -------

!       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A BLOCK

!**   INTERFACE.
!     ----------

!       *CALL* *MUBUF (IU08)*
!          *IU01*  -  LOGICAL INPUT UNIT OF TOPOGRAPHIC DATA.
!          *IU08*    - LOGICAL UNITS FOR OUTPUT OF GRID BLOCKING
!                      COMMON UBUF (UNFORMATED)

!     METHOD.
!     -------

!       THE INDICES OF THE NEXT POINTS ON LAT. AND LONG. ARE
!       COMPUTED. ZERO INDICATES A LAND POINT IS NEIGHBOUR.
!       THE FINAL COMMON UBUF IS WRITTEN OUT.

!     EXTERNALS.
!     ----------

!       *OUTUBUF*   - WRITE OUT COMMON UBUF.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NFRE_RED
      USE YOWMAP   , ONLY : BLK2GLO  ,NGX      ,NGY      ,NIBLO     ,  &
     &                      IPER     ,IRGG     ,NLONRGG  ,LLOBSTRCT
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU01
      INTEGER(KIND=JWIM), INTENT(IN) :: IU08(0:NPROPAGS)
      INTEGER(KIND=JWIM), INTENT(IN) :: NPROPAGS


      INTEGER(KIND=JWIM) :: NFREMAX, IX, IXLP
      INTEGER(KIND=JWIM) :: IJ, IJP, I, K, IP, IH, IS, M
      INTEGER(KIND=JWIM) :: IC, ICP, ICL, ICR
      INTEGER(KIND=JWIM), DIMENSION(NGX, NGY) :: IDUM
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: KDUM

      CHARACTER(LEN=5) :: CX
      CHARACTER(LEN=11) :: FORMAT

! ----------------------------------------------------------------------


!     DETERMINE THE REAL REDUCTION FACTORS TO BE USED IN THE PROPAGATION

!     CHECK INOUT IS CONSISTENT WITH CURRENT SETUP
      IF (LLOBSTRCT) THEN
        READ (IU01,'(I4)') NFREMAX
        IF (NFREMAX /= NFRE_RED ) THEN
          WRITE (IU06,*) ' *******************************************'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. MUBUF         *'
          WRITE (IU06,*) ' * NFREMAX MUST BE = NFRE_RED              *'
          WRITE (IU06,*) ' * NFREMAX = ',NFREMAX
          WRITE (IU06,*) ' * NFRE_RED = ',NFRE_RED 
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                 *'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF
      ENDIF

      ALLOCATE(KDUM(NIBLO))

      WRITE(CX,'(I5.5)') NLONRGG(1)
      FORMAT='('//CX//'I4)'

      DO M=1,NFRE_RED  ! loop over frequencies used for the propagation

!     KOBSLAT
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                READ(IU01,FORMAT) (IDUM(IX,K),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
            DO IP = 1,NIBLO
              I = BLK2GLO%IXLG(IP)
              K = BLK2GLO%KXLT(IP)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          DO IJP=0,NPROPAGS
            WRITE(IU08(IJP)) KDUM 
          ENDDO
        ENDDO

!     KOBSLON
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                READ(IU01,FORMAT) (IDUM(IX,K),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
            DO IP = 1,NIBLO
              I = BLK2GLO%IXLG(IP)
              K = BLK2GLO%KXLT(IP)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          DO IJP=0,NPROPAGS
            WRITE(IU08(IJP)) KDUM 
          ENDDO
        ENDDO

!     KOBSRLAT
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                READ(IU01,FORMAT)(IDUM(IX,K),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
            DO IP = 1,NIBLO
              I = BLK2GLO%IXLG(IP)
              K = BLK2GLO%KXLT(IP)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          WRITE(IU08(1)) KDUM 
        ENDDO

!     KOBSRLON
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                READ(IU01,FORMAT)(IDUM(IX,K),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
            DO IP = 1,NIBLO
              I = BLK2GLO%IXLG(IP)
              K = BLK2GLO%KXLT(IP)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          WRITE(IU08(1)) KDUM 
        ENDDO

!     KOBSCOR
        DO IS=1,4
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                READ(IU01,FORMAT) (IDUM(IX,K),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
            DO IP = 1,NIBLO
              I = BLK2GLO%IXLG(IP)
              K = BLK2GLO%KXLT(IP)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          WRITE(IU08(2)) KDUM 
        ENDDO

      ENDDO  ! end loop over frequencies

      DEALLOCATE(KDUM)

      END SUBROUTINE MUBUF
