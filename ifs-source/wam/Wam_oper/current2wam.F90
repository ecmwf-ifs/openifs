! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE CURRENT2WAM (FILNM, IREAD, CDATEIN,        &
     &                        IFROMIJ, JFROMIJ,             &
     &                        NXS, NXE, NYS, NYE, FIELDG,   &
     &                        UCUR, VCUR)

!--------------------------------------------------------------------

!**** *CURRENT2WAM* 

!     J. BIDLOT    ECMWF  NOVEMBER 2002.
!     J. BIDLOT    ECMWF  APRIL 2010: USE GRIBAPI 

!     PURPOSE.
!     --------

!     TO READ THE OCEAN CURRENT IN AND TO TRANSFER THEM ONTO THE
!     WAVE MODEL GRID.

!**   INTERFACE
!     ---------
!     *CALL* *CURRENT2WAM(FILNM, IREAD, CDATEIN,
!                         IFROMIJ, JFROMIJ,
!                         NXS, NXE, NYS, NYE, FIELDG,
!                         UCUR, VCUR)

!     *FILNM*     DATA INPUT FILENAME.
!     *IREAD*     RANK OF THE PROCESS WHICH INPUTS THE DATA. 
!     *CDATEIN*   DATE OF THE DECODED DATA. 
!     *IFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *JFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *NXS:NXE*  FIRST DIMENSION OF FIELDG
!     *NYS:NYE*  SECOND DIMENSION OF FIELDG
!     *FIELDG*   INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!     *UCUR*      U-COMPONENT OF SURFACE CURRENT
!     *VCUR*      V-COMPONENT OF SURFACE CURRENT


!     METHOD.
!     -------


!    EXTERNAL.
!    ---------

!---------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCURR  , ONLY : CURRENT_MAX
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWMAP   , ONLY : IRGG     ,NLONRGG
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI
      USE YOWPARAM , ONLY : NGY
      USE YOWPCONS , ONLY : ZMISS    ,EPSMIN
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06
      USE YOWPD, ONLY : MNP => npa
      USE YOWUNPOOL ,ONLY : LLUNSTR

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BARRIER, MPL_BROADCAST
      USE YOWGRIB , ONLY : IGRIB_NEW_FROM_MESSAGE, IGRIB_OPEN_FILE, &
                         & IGRIB_READ_FROM_FILE, IGRIB_RELEASE, &
                         & JPGRIB_BUFFER_TOO_SMALL, JPGRIB_END_OF_FILE, &
                         & JPGRIB_SUCCESS, JPKSIZE_T

!-----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "grib2wgrid.intfb.h"
#include "kgribsize.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      CHARACTER(LEN=24), INTENT(IN) :: FILNM
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEIN
      INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), DIMENSION(NXS:NXE, NYS:NYE), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(OUT) :: UCUR, VCUR


      INTEGER(KIND=JWIM) :: NBIT = 1000000

      INTEGER(KIND=JWIM) :: KFILE_HANDLE1
      INTEGER(KIND=JWIM) :: LFILE, KGRIB_HANDLE
      INTEGER(KIND=JWIM) :: KRET, IVAR, KPLENG, KLEN, KLENG, IDM, KDM, MDM
      INTEGER(KIND=JWIM) :: IRET, ISIZE
      INTEGER(KIND=JWIM) :: KK, MM
      INTEGER(KIND=JWIM) :: IFORP, IPARAM, KZLEV, ICHNK, IJ, IC, IX, JY
      INTEGER(KIND=JWIM) :: IBUF(2)
      INTEGER(KIND=JWIM), ALLOCATABLE :: INGRIB(:)
      INTEGER(KIND=JWIM) :: NLONRGG_LOC(NGY)

      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JWRB), PARAMETER :: WLOWEST=0.0001_JWRB
      REAL(KIND=JWRB) :: ZDUM1, ZDUM2
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: FIELD(NXS:NXE, NYS:NYE)

      CHARACTER(LEN=14) :: CDATEIN_OLD 

      LOGICAL :: FRSTIME

      DATA FRSTIME / .TRUE. /
      SAVE KFILE_HANDLE1
      SAVE FRSTIME

! --------------------------------------------------------------------  

      IF (LHOOK) CALL DR_HOOK('CURRENT2WAM',0,ZHOOK_HANDLE)

      IF (LLUNSTR) THEN
        NLONRGG_LOC(:)=MNP
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
      ENDIF

      LFILE = LEN_TRIM(FILNM)
      ZDUM1 = 0._JWRB
      ZDUM2 = 0._JWRB

      IF (IRANK == IREAD) THEN

        IF (FRSTIME) THEN
          KFILE_HANDLE1=-99
          CALL IGRIB_OPEN_FILE(KFILE_HANDLE1,FILNM(1:LFILE),'r')
          FRSTIME = .FALSE.
        ENDIF  

        ISIZE=NBIT
      ENDIF


      READCURRENT: DO IVAR=1,2

        IF (IRANK == IREAD) THEN
1021        ISIZE=NBIT
            KBYTES=ISIZE*NPRECI
            IF (.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
            CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE1,INGRIB,KBYTES,IRET)
            IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
!!!           *IGRIB_READ_FROM_FILE* does not read through the file if
!!!            the size is too small, so figure out the size and read again.
              CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'CURRENT2WAM')
              DEALLOCATE(INGRIB)
              GOTO 1021
            ELSEIF (IRET == JPGRIB_END_OF_FILE) THEN
              WRITE(IU06,*) '**************************************'
              WRITE(IU06,*) '* CURRENT2WAM: END OF FILE ENCOUNTED *'
              WRITE(IU06,*) '**************************************'
              CALL ABORT1
            ELSEIF (IRET /= JPGRIB_SUCCESS) THEN
               WRITE(IU06,*) '************************************'
               WRITE(IU06,*) '* CURRENT2WAM: FILE HANDLING ERROR *'
               WRITE(IU06,*) '************************************'
              CALL ABORT1
            ENDIF

          WRITE(IU06,*) ' SUB. CURRENT2WAM - READ FROM ',FILNM
          WRITE(IU06,*) ' IVAR = ',IVAR 

        ENDIF

!       BROADCAST GRIB DATA TO OTHER PE'S
        CALL GSTATS(622,0)
        IF (NPROC > 1) THEN

          CALL MPL_BARRIER(CDSTRING='CURRENT2WAM: INGRIB ')

          IF (IRANK == IREAD) THEN
            IBUF(1)=ISIZE
            IBUF(2)=KLEN
          ENDIF
          CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=IVAR,           &
     &                       CDSTRING='CURRENT2WAM IBUF:')
          IF (IRANK /= IREAD) THEN
            ISIZE=IBUF(1)
            KLEN=IBUF(2)
            ALLOCATE(INGRIB(ISIZE))
          ENDIF

          CALL MPL_BROADCAST(INGRIB(1:ISIZE),KROOT=IREAD,KTAG=IVAR,     &
     &                       CDSTRING='CURRENT2WAM INGRIB:')

        ENDIF
        CALL GSTATS(622,1)
! ----------------------------------------------------------------------

!*      UNPACK GRIB FIELDS.
!       -------------------
        CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

        KK=0
        MM=0
        CALL GRIB2WGRID (IU06, NPROMA_WAM,                              &
     &                   KGRIB_HANDLE, INGRIB, ISIZE,                   &
     &                   LLUNSTR,                                       &
     &                   NGY, IRGG, NLONRGG_LOC,                        &
     &                   NXS, NXE, NYS, NYE,                            &
     &                   FIELDG%XLON, FIELDG%YLAT,                      &
     &                   ZMISS, ZDUM1, ZDUM2,                           &
     &                   CDATEIN, IFORP, IPARAM, KZLEV, KK, MM, FIELD)


        CALL IGRIB_RELEASE(KGRIB_HANDLE)

        IF (IVAR == 2 .AND. CDATEIN /= CDATEIN_OLD) THEN
          WRITE(IU06,*) ' *****************************************'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. CURRENT2WAM    *'
          WRITE(IU06,*) ' *      ============================     *'
          WRITE(IU06,*) ' * DATE FOR PARAMETER ',IPARAM 
          WRITE(IU06,*) ' * IS ',CDATEIN
          WRITE(IU06,*) ' * IT SHOULD BE ',CDATEIN_OLD
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *****************************************'
          CALL FLUSH(IU06)
          CALL ABORT1
       ENDIF

       CDATEIN_OLD=CDATEIN


        IF (IPARAM == 131) THEN
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
          DO ICHNK = 1, NCHNK
            DO IJ = 1, NPROMA_WAM 
              IX = IFROMIJ(IJ,ICHNK)
              JY = JFROMIJ(IJ,ICHNK)
              UCUR(IJ,ICHNK) = FIELD(IX,JY)
!             SOME WAM MODEL GRID POINTS MAY HAVE A MISSING DATA FROM
!             OCEAN MODEL. THEY ARE SET TO 0.
!             0. WILL BE USED TO DETECT THE INABILITY TO COMPUTE THE GRADIANT
              IF (ABS(UCUR(IJ,ICHNK)) <= WLOWEST) UCUR(IJ,ICHNK)=0.0_JWRB
              IF (UCUR(IJ,ICHNK) <= ZMISS) UCUR(IJ,ICHNK)=0.0_JWRB
              UCUR(IJ,ICHNK)=SIGN(MIN(ABS(UCUR(IJ,ICHNK)),CURRENT_MAX),UCUR(IJ,ICHNK))
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO

        ELSEIF (IPARAM == 132) THEN
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
          DO ICHNK = 1, NCHNK
            DO IJ = 1, NPROMA_WAM 
              IX = IFROMIJ(IJ,ICHNK)
              JY = JFROMIJ(IJ,ICHNK)
              VCUR(IJ,ICHNK) = FIELD(IX,JY)
!             SOME WAM MODEL GRID POINTS MAY HAVE A MISSING DATA FROM
!             OCEAN MODEL. THEY ARE SET TO 0.
!             0. WILL BE USED TO DETECT THE INABILITY TO COMPUTE THE GRADIANT
              IF (ABS(VCUR(IJ,ICHNK)) <= WLOWEST) VCUR(IJ,ICHNK)=0.0_JWRB
              IF (VCUR(IJ,ICHNK) <= ZMISS) VCUR(IJ,ICHNK)=0.0_JWRB
              VCUR(IJ,ICHNK)=SIGN(MIN(ABS(VCUR(IJ,ICHNK)),CURRENT_MAX),VCUR(IJ,ICHNK))
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO

        ELSE
          WRITE(IU06,*) ' *****************************************'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. CURRENT2WAM    *'
          WRITE(IU06,*) ' *      ============================     *'
          WRITE(IU06,*) ' * THE INPUT PARAMETER ',IPARAM 
          WRITE(IU06,*) ' * IS NOT RECOGNISED AS BEING A CURRENT !*' 
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *****************************************'
          CALL ABORT1
        ENDIF

        DEALLOCATE(INGRIB)

      ENDDO READCURRENT 

      IF (LHOOK) CALL DR_HOOK('CURRENT2WAM',1,ZHOOK_HANDLE)

      END SUBROUTINE CURRENT2WAM
