! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE  BASCOE_READ_COLDAT( CD_FILENM, KN_OUT, PZ, PVAL, CD_EXT_VALS )
!----------------------------------------------------------------------
!     General function to READ 2-cols (ZF, ZVALF) ASCII input file
!     and interpolate it  to (PZ, PVAL) of dimension KN_OUT using
!     subroutine BASCOE_INTERP8.
!     CD_FILENM is the input file name
!     Expected format for input file *.dat :
!         - (comment lines, starting with '!', to skip)
!         - an integer = number of data points -> nf   (free format)
!         - (comment lines, starting with '!', to skip)
!         - nf lines of 2 reals: x (or z), f(x) (or f(z)) (free format)
!     CD_EXT_VALS: see BASCOE_INTERP8
!
!----------------------------------------------------------------------
USE YOMLUN             , ONLY : NULOUT
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

!----------------------------------------------------------------------
!       ... Dummy args
!----------------------------------------------------------------------
  CHARACTER(LEN=*), intent(in) :: CD_FILENM
  INTEGER(KIND=JPIM), intent(in) :: KN_OUT
  REAL(KIND=JPRB), dimension(KN_OUT), intent(in) :: PZ
  REAL(KIND=JPRB), dimension(KN_OUT), intent(inout) :: PVAL
  CHARACTER(LEN=*), intent(in) :: CD_EXT_VALS

!-----------------------------------------------------------------------
!       ... parameters
!-----------------------------------------------------------------------
  CHARACTER(LEN=*), PARAMETER    :: CL_MY_NAME   = 'BASCOE_READ_COLDAT'

  INTEGER(KIND=JPIM), parameter :: INMAX = 10000
  INTEGER(KIND=JPIM), PARAMETER :: IUTMP=77        ! default LUN for external files ?

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM) :: ios, iz, I_NF
  REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE
  REAL(KIND=JPRB), dimension(INMAX) :: ZF, ZVALF
  LOGICAL :: LL_OK

#include "abor1.intfb.h"
#include "bascoe_interp8.intfb.h"
#include "bascoe_cskip.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_READ_COLDAT',0,ZHOOK_HANDLE )

!-----------------------------------------------------------------------
!       ... READ file
!-----------------------------------------------------------------------
      WRITE(NULOUT,*) CL_MY_NAME//': reading >'// trim(CD_FILENM) // '<'
      OPEN( unit = IUTMP, &
     &      file = trim(CD_FILENM), &
     &      form = 'FORMATTED', &
     &      status = 'OLD', &
     &      iostat = ios )
      IF( ios /= 0 ) then
         WRITE(NULOUT,*) CL_MY_NAME//': Error opening >'//TRIM(CD_FILENM)//'< , iostat=',ios
         CALL ABOR1(CL_MY_NAME//': Error opening >'//TRIM(CD_FILENM))
      ENDIF
      CALL BASCOE_CSKIP( '!',IUTMP )
      READ( IUTMP, * ) I_NF
      IF( i_NF > INMAX ) then
         WRITE(NULOUT,*) CL_MY_NAME//': Error number of lines:',i_NF,' > INMAX:',INMAX
         CALL ABOR1(CL_MY_NAME//': Error too many lines')
      ENDIF
      CALL BASCOE_CSKIP( '!',IUTMP )
      DO iz = 1, i_NF
            READ(IUTMP,*) ZF(iz), ZVALF(iz)
      ENDDO
      CLOSE( IUTMP )
!-----------------------------------------------------------------------
!       ... Put data from file on output grid "PZ"
!-----------------------------------------------------------------------
       CALL BASCOE_INTERP8( KN_OUT, PZ, PVAL, I_NF, ZF(1:I_NF), ZVALF(1:I_NF), &
     &               LL_OK,  CD_EXT_VALS)
       IF( .NOT. LL_OK ) then
         WRITE(NULOUT,*) CL_MY_NAME//': Error interpolation to target grid failed'
         CALL ABOR1(CL_MY_NAME//': Error interpolation to target grid failed')
       ENDIF

IF (LHOOK) CALL DR_HOOK('BASCOE_READ_COLDAT',1,ZHOOK_HANDLE )

END SUBROUTINE BASCOE_READ_COLDAT
