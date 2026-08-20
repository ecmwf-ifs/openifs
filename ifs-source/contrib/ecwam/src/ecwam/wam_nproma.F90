! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WAM_NPROMA(IJS, IJL, MTHREADS, NPROMA)
! ----------------------------------------------------------------------

!**** *WAM_NPROMA* - ADJUSTS DOWN NPROMA BASED ON IJS, IJL AND MTHREADS
!                    TO INSURE THAT NTOT IS DISTRIBUTED AS EQUALLY
!                    AS POSSIBLE ON MTHREADS THREADS

!**   INTERFACE.
!     ----------


!       *IJS*     - INDEX OF FIRST GRIDPOINT
!       *IJL*     - INDEX OF LAST GRIDPOINT
!       *MTHREADS*- NUMBER OF OPENMP THREADS
!       *NPROMA*  - MAXIMUM NUMBER OF POINTS PER THREAD PER CALL

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, MTHREADS
      INTEGER(KIND=JWIM), INTENT(INOUT) :: NPROMA
      INTEGER(KIND=JWIM) :: NLOOP, NREMAIN, NPROMA_MAX


      INTEGER(KIND=JWIM) :: NTOT
! ----------------------------------------------------------------------

      NTOT=IJL-IJS+1

      NPROMA=MIN(NPROMA,NTOT)
      NPROMA_MAX=NPROMA
      NLOOP=NTOT/(MTHREADS*NPROMA)
      NREMAIN=NTOT-NLOOP*MTHREADS*NPROMA
      IF (NREMAIN > 0 .AND. MTHREADS > 1 .AND. NTOT > 0) THEN
        DO WHILE (NREMAIN < (NPROMA*MTHREADS) .AND. NPROMA > 1)
          NPROMA=NPROMA-1
          NREMAIN=NTOT-NLOOP*MTHREADS*NPROMA
        ENDDO
        NPROMA=MIN(NPROMA+1,NPROMA_MAX)
      ENDIF

      END SUBROUTINE WAM_NPROMA
