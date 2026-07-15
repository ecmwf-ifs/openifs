! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUTINT(CDATE, CDATED, IFCST, BOUT)

! ----------------------------------------------------------------------

!**** *OUTINT* - OUTPUT OF INTEGRATED FIELDS WITHOUT THE I/O SERVER

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *OUTINT(CDATE, CDATED, IFCST, BOUT)
!          *CDATE*   - ACTUAL DATE AND TIME. IF NOT A FORECAST (IFCST<=0), OTHERWISE
!                      DATE AND TIME OF THE START OF THE FOERCAST.
!          *CDATED*  - DATE USED FOR NAMING THE OUTPUT FILE (if data written to files)
!          *IFCST*   - FORECAST STEP IN HOURS.
!          *BOUT*    - OUTPUT PARAMETERS BUFFER


!     EXTERNALS.
!     ----------
!       *OUTGRID*
!       *WGRIBENOUT*

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : FFLAG    ,FFLAG20  ,GFLAG    ,GFLAG20  ,    &
     &                      JPPFLAG  ,LFDB     ,IPFGTBL  ,ITOBOUT  ,    &
     &                      INFOBOUT , LOUTINT  ,NIPRMOUT
      USE YOWGRID  , ONLY : DELPHI   ,NPROMA_WAM, NCHNK
      USE YOWINTP  , ONLY : GOUT
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &                      NGX      ,NGY      ,CLDOMAIN,               &
     &                      AMONOP   ,ZDELLO   ,NLONRGG
      USE YOWMPP   , ONLY : IRANK
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,CDTINTT  ,MARSTYPE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH
      USE YOWUNIT  , ONLY : IU20

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR
      USE YOWGRIB  , ONLY : IGRIB_OPEN_FILE, IGRIB_CLOSE_FILE

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "grstname.intfb.h"
#include "outgrid.intfb.h"
#include "out_onegrdpt.intfb.h"
#include "wgribenout.intfb.h"

      CHARACTER(LEN=14), INTENT(IN) :: CDATE, CDATED
      INTEGER(KIND=JWIM), INTENT(IN) :: IFCST
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT


      INTEGER(KIND=JWIM) :: I, J, ITG, IFLAG, IT
      INTEGER(KIND=JWIM) :: ICOUNT      ! Field counter
      INTEGER(KIND=JWIM) :: IPARAM, ITABLE, IZLEV
      INTEGER(KIND=JWIM) :: IERR
      INTEGER(KIND=JWIM) :: LFILE, IUOUT 
  
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=3) :: FILEID
      CHARACTER(LEN=296) :: OUTFILEN

      LOGICAL :: LLIUOUTOPEN

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTINT',0,ZHOOK_HANDLE)

      LOUTINT=.TRUE.

!     1. COLLECT INTEGRATED PARAMETERS FOR OUTPUT ON SELECTED PE's 
!        i.e BOUT => GOUT
!        ---------------------------------------------------------
      CALL OUTGRID(BOUT)

      CALL GSTATS(1704,0)

!     2. ONE GRID POINT OUTPUT
!        ---------------------
      IF (CLDOMAIN == 's') CALL OUT_ONEGRDPT(IU06)


!*    3. OUTPUT IN PURE BINARY FORM (obsolete)
!        --------------------------
      IF (FFLAG20 .AND. CDTINTT == CDTPRO ) THEN
        DO IFLAG = 1, JPPFLAG
          IF (FFLAG(IFLAG) .AND. IRANK == IPFGTBL(IFLAG)) THEN
            WRITE (IU20) CDTPRO, NGX, NGY, IRGG
            WRITE (IU20) AMOWEP,AMOSOP,AMOEAP,AMONOP
            WRITE (IU20) ZDELLO, NLONRGG, DELPHI
            WRITE (IU20) GOUT(IFLAG,:,:)
          ENDIF
        ENDDO
      ENDIF

!*    4.  PACK INTEGRATED PARAMETERS INTO GRIB AND OUTPUT
!         -----------------------------------------------
      IF (GFLAG20 .AND. CDTINTT == CDTPRO ) THEN

!       DEFINE OUTPUT FILE FOR GRIB DATA (if not written to FDB)
        IF (.NOT.LFDB .AND. IRANK == 1 ) THEN
       !  output to file should only take place on one PE
         IF (MARSTYPE == 'fg') THEN
           ! fg only exist in uncoupled anlysis experiments
           FILEID = 'MFG'
         ELSE
           FILEID = 'MPP'
         ENDIF
         CALL GRSTNAME(CDATED,CDATEF,IFCST,FILEID,ICPLEN,CPATH,OUTFILEN)
         LFILE=LEN_TRIM(OUTFILEN)
         CALL IGRIB_OPEN_FILE(IUOUT,OUTFILEN(1:LFILE),'w')
         WRITE(IU06,*) '  INTEGRATED PARAMETERS WRITTEN TO FILE ',OUTFILEN(1:LFILE)
         LLIUOUTOPEN=.TRUE.
        ELSE
          LLIUOUTOPEN=.FALSE.
          IUOUT=0
        ENDIF

!       OUTPUT:
        ICOUNT=0
        DO IFLAG = 1, JPPFLAG
          IF (GFLAG(IFLAG)) THEN
            IF (IRANK == IPFGTBL(IFLAG)) THEN
              ICOUNT=ICOUNT+1
              IF (ICOUNT > SIZE(GOUT,1)) THEN
                WRITE(NULERR,*) ' -------------------------------------'
                WRITE(NULERR,*) ' ERROR in OUTINT '
                WRITE(NULERR,*) ' ACCESSING MORE FIELDS THAN AVAILABLE'
                WRITE(NULERR,*) ' SIZE(GOUT,1) = ',SIZE(GOUT,1) 
                WRITE(NULERR,*) ' -------------------------------------'
                CALL ABORT1
              ENDIF

              IT=ITOBOUT(IFLAG)
              ITABLE=INFOBOUT(IT,1)
              IPARAM=INFOBOUT(IT,2)
              IZLEV=INFOBOUT(IT,3)

              CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, GOUT(ICOUNT,:,:),  &
     &                        ITABLE, IPARAM, IZLEV, 0 , 0,             &
     &                        CDATE, IFCST, MARSTYPE, LFDB, IUOUT)
            ENDIF
          ENDIF
        ENDDO

        IF (LLIUOUTOPEN) THEN
          CALL IGRIB_CLOSE_FILE(IUOUT)
        ENDIF

      ENDIF ! end grib output

!     CLEAN-UP
      IF (ALLOCATED(GOUT)) DEALLOCATE(GOUT)

      CALL GSTATS(1704,1)

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTINT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTINT
