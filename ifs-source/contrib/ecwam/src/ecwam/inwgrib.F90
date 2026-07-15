! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INWGRIB (FILNM, IREAD, CDATE, IPARAM, KZLEV,  & 
     &                    NXS, NXE, NYS, NYE, FIELDG, FIELD)

! -----------------------------------------------------------------     

!***  *INWGRIB* - READS IN AND UNPACKS WAVE MODEL GRIB FIELDS
!                 IT WILL ALSO INTERPOLATE THE DATA TO THE MODEL
!                 GRID (IF NECESSARY).
!                 (Except for wave spectra. See *GETSPEC*)
!                 MESSAGE PASSING WILL BE USE TO DISTRIBUTED THE
!                 CODED DATA TO ALL PE'S BEFORE DECODING.

!      J. BIDLOT    ECMWF    APRIL 2010. 

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL INWGRIB* (FILNM, IREAD, CDATE, IPARAM, KZLEV,
!    &                 NXS, NXE, NYS, NYE, FIELDG, FIELD)

!        *FILNM*  - DATA INPUT FILENAME. 
!        *IREAD*  - ACCESS TO FILE ONLY FOR PE=IREAD.
!        *CDATE*  - DATE/TIME OF THE DATA READ.         
!        *IPARAM* - PARAMETER NUMBER. 
!        *KZLEV*  - REFERENCE LEVEL IN FULL METER
!                   SHOULD BE 0 EXCEPT FOR 233, 245 AND 249 WHERE IT
!                   MIGHT BE DIFFERENT THAN ZERO. 
!        *NXS:NXE*  FIRST DIMENSION OF FIELDG
!        *NYS:NYE*  SECOND DIMENSION OF FIELDG
!        *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!        *FIELD*  - UNPACKED DATA.

!     EXTERNALS.                                                        
!     ----------                                                        

!     GRIBAPI 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWMAP   , ONLY : NGY      ,NIBLO    ,IRGG     ,NLONRGG
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI 
      USE YOWPCONS , ONLY : ZMISS
      USE YOWTEST  , ONLY : IU06
#ifdef WAM_HAVE_UNWAM
      USE YOWPD, ONLY : MNP => npa
#endif

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BARRIER, MPL_BROADCAST
      USE YOWGRIB, ONLY : IGRIB_OPEN_FILE, IGRIB_CLOSE_FILE, IGRIB_RELEASE, &
                        & IGRIB_NEW_FROM_MESSAGE, IGRIB_READ_FROM_FILE, &
                        & JPKSIZE_T, JPGRIB_BUFFER_TOO_SMALL, &
                        & JPGRIB_END_OF_FILE, JPGRIB_SUCCESS
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "grib2wgrid.intfb.h"
#include "kgribsize.intfb.h"

      CHARACTER(LEN=24), INTENT(IN) :: FILNM
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      CHARACTER(LEN=14), INTENT(OUT) :: CDATE
      INTEGER(KIND=JWIM), INTENT(OUT) :: IPARAM, KZLEV
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), INTENT(OUT) :: FIELD(NXS:NXE, NYS:NYE)


      INTEGER(KIND=JWIM) :: NBIT

      INTEGER(KIND=JWIM) :: IFORP
      INTEGER(KIND=JWIM) :: LFILE, KFILE_HANDLE, KGRIB_HANDLE
      INTEGER(KIND=JWIM) :: IRET, ISIZE
      INTEGER(KIND=JWIM) :: KK, MM
      INTEGER(KIND=JWIM) :: IBUF(2) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: INGRIB(:)
      INTEGER(KIND=JWIM) :: NLONRGG_LOC(NGY)

      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LLEXIST

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INWGRIB',0,ZHOOK_HANDLE)

      NBIT=NIBLO

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        NLONRGG_LOC(:)=MNP
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
      ENDIF

!     READ DATA ON PE IREAD
      IF (IRANK == IREAD) THEN
        LLEXIST=.FALSE.
        LFILE = LEN_TRIM(FILNM)
        INQUIRE(FILE=FILNM(1:LFILE),EXIST=LLEXIST)
        IF (.NOT. LLEXIST) THEN
          WRITE (IU06,*) '*************************************'
          WRITE (IU06,*) '*                                   *'
          WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
          WRITE (IU06,*) '*  IN INWGRIB:                      *'
          WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILNM
          WRITE (IU06,*) '*                                   *'
          WRITE (IU06,*) '*************************************'
          WRITE (*,*) '*************************************'
          WRITE (*,*) '*                                   *'
          WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
          WRITE (*,*) '*  IN INWGRIB:                      *'
          WRITE (*,*) '*  COULD NOT FIND FILE ',FILNM
          WRITE (*,*) '*                                   *'
          WRITE (*,*) '*************************************'
          CALL ABORT1
        ENDIF

        CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILNM(1:LFILE),'r')

1021    ISIZE=NBIT
        KBYTES=ISIZE*NPRECI
        IF (.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
          CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,KBYTES,IRET)
        IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
!!!       *IGRIB_READ_FROM_FILE* does not read through the file if
!!!       the size is too small, so figure out the size and read again.
          CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'INWGRIB')
          DEALLOCATE(INGRIB)
          GOTO 1021
        ELSEIF (IRET == JPGRIB_END_OF_FILE) THEN
          WRITE(IU06,*) '**********************************'
          WRITE(IU06,*) '* INWGRIB: END OF FILE ENCOUNTED'
          WRITE(IU06,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(IU06,*) '**********************************'
          CALL ABORT1
          ELSEIF (IRET /= JPGRIB_SUCCESS) THEN
          WRITE(IU06,*) '**********************************'
          WRITE(IU06,*) '* INWGRIB: FILE HANDLING ERROR'
          WRITE(IU06,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(IU06,*) '**********************************'
          CALL ABORT1
        ENDIF
      ENDIF

      WRITE(IU06,*) ' SUB. INWGRIB - READ FROM ',FILNM

      CALL MPL_BARRIER(CDSTRING='INWGRIB: DATA READ IN')

!     SEND GRIB DATA TO THE OTHER PE'S
      IF (NPROC > 1) THEN
        CALL GSTATS(619,0)
        IF (IRANK == IREAD) THEN
          IBUF(1)=ISIZE
          IBUF(2)=KBYTES
        ENDIF
        CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=1,                &
     &                     CDSTRING='INWGRIB IBUF:')
        IF (IRANK /= IREAD) THEN
          ISIZE=IBUF(1)
          KBYTES=IBUF(2)
          ALLOCATE(INGRIB(ISIZE))
        ENDIF

        CALL MPL_BROADCAST(INGRIB(1:ISIZE),KROOT=IREAD,KTAG=2,          &
     &                     CDSTRING='INWGRIB: INGRIB')
        CALL GSTATS(619,1)
      ENDIF

!     DECODE THE GRIB DATA
!     (and interpolate to model grid if necessary)

      KGRIB_HANDLE=-99
      CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

      CALL GRIB2WGRID (IU06, NPROMA_WAM,                                &
     &                 KGRIB_HANDLE, INGRIB, ISIZE,                     &
     &                 LLUNSTR,                                         &
     &                 NGY, IRGG, NLONRGG_LOC,                          &
     &                 NXS, NXE, NYS, NYE,                              &
     &                 FIELDG%XLON, FIELDG%YLAT,                        &
     &                 ZMISS, PPREC, PPEPS,                             &
     &                 CDATE, IFORP, IPARAM, KZLEV, KK, MM, FIELD)

      CALL IGRIB_RELEASE(KGRIB_HANDLE)

      IF (ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

      IF (LHOOK) CALL DR_HOOK('INWGRIB',1,ZHOOK_HANDLE)

      END SUBROUTINE INWGRIB
