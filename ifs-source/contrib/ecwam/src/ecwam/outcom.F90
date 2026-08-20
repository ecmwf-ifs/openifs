! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUTCOM (IU07, BATHY)

! ----------------------------------------------------------------------

!**** *OUTCOM* - ROUTINE TO WRITE MODULES TO DISK

!     H.GUNTHER            ECMWF       04/04/1990

!     J. BIDLOT            ECMWF       10/1998.
!!!!!!                     COMMON BLOCKS HAVE BEEN CONVERTED TO MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN, 
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE READPRE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     -------

!       TO WRITE OUT THE COMPUTED MODULES
!       (MODULE UBUF IS WRITTEN IN MUBUF)

!**   INTERFACE.
!     ----------

!       *CALL* *OUTCOM (IU07, BATHY)*
!          *IU07*   - LOGICAL UNIT FOR  UNFORMATED WRITE.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIBHD, ONLY : IMDLGRBID_G
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWMAP   , ONLY : NGX      ,NGY      , IPER     ,IRGG    ,    &
     &                      AMOWEP   ,AMOSOP   ,AMOEAP    ,AMONOP  ,    &
     &                      XDELLA   ,XDELLO   ,NLONRGG
      USE YOWABORT, ONLY : WAM_ABORT

#ifdef WAM_HAVE_UNWAM
      USE UNWAM , ONLY : UNWAM_OUT
      USE YOWUNPOOL ,ONLY : LPREPROC
#endif

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "outnam.intfb.h"
 
      INTEGER(KIND=JWIM), INTENT(IN) :: IU07
      REAL(KIND=JWRB), INTENT(IN) :: BATHY(NGX, NGY)

      INTEGER(KIND=JWIM) :: IDUM, K, M, L
      INTEGER(KIND=JWIM) :: NKIND !Precision used when writing

! ----------------------------------------------------------------------

!     WRITE IDENTIFIERS (MAKE SURE TO UPDATE THE VALUES IF YOU CHANGE
!     ANYTHING TO THE MODEL). 
      
      NKIND = KIND(AMOSOP)

      WRITE(IU07) NKIND, IMDLGRBID_G 


! ----------------------------------------------------------------------

!*    GRID INFORMATION
!     ----------------

      WRITE(IU07) NGX, NGY

      WRITE(IU07) NLONRGG

      WRITE(IU07) IPER, IRGG, AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO

      WRITE(IU07) BATHY

! ----------------------------------------------------------------------

!*    WRITE NAMELIST PARWAM.
!     ----------------------
      
      CALL OUTNAM  (NGX, NGY)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        IF (LPREPROC) THEN
          CALL UNWAM_OUT(IU07)
        ENDIF
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      END IF

      END SUBROUTINE OUTCOM
