! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE OUTINT_IO_SERV (NIPRMOUT, BOUT, INFOBOUT, MARSTYPE, CDATE, IFCST)

!----------------------------------------------------------------------

!**** *OUTINT*  WRITES INTEGRATED PARAMETERS USING IO SERVER

!     J. HAWKES   ECMWF  APRIL 2018

!     See also: wv_io_serv_rec
!     See also: outwspec_io_serv

!-------------------------------------------------------------------

! Wave Model
USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, NTOTIJ, KIJL4CHNK
USE YOWPARAM , ONLY : NIBLO
USE YOWCOUP  , ONLY : KCOUSTEP
USE YOWGRIBHD, ONLY : NDATE_TIME_WINDOW_END

! IFS/Coupling
USE MPL_MODULE
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOWCOUT  , ONLY : LRSTST0
USE YOWCOUP  , ONLY : IFSNSTEP, IFSTSTEP

! IO Server
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE YOMIO_SERV, ONLY : IO_SERV_C001
USE YOMIO_SERV_MAP_PLAN, ONLY : IO_SERV_SEND_PLAN
USE YOMIO_SERV_HDR, ONLY : NIO_SERV_HDR_IDOM_WV
USE IOMULTIBUF_MOD, ONLY : IOMULTIBUF_INCR_IDX, IOMULTIBUF_INIT_IDX

!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN)  :: NIPRMOUT                              ! Number of integrated params
REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT             ! Data
INTEGER(KIND=JWIM), DIMENSION(NIPRMOUT, 3), INTENT(IN) :: INFOBOUT      ! grib table num, grib paramid, grib ref level
CHARACTER(LEN=2), INTENT(IN)    :: MARSTYPE 
CHARACTER(LEN=14), INTENT(IN)   :: CDATE
INTEGER (KIND=JWIM), INTENT(IN) :: IFCST

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

TYPE (IO_SERV_SEND_PLAN)      :: YLIOSMPP          ! Send map plan
TYPE (IOFLDDESC), ALLOCATABLE :: YLFLDSC (:)       ! Field descriptor array
INTEGER(KIND=JWIM)            :: IFNUM, JFNUM      ! Number of fields
INTEGER(KIND=JWIM)            :: IGPTOTG, IGPTOT   ! Number of points (global & local)
INTEGER (KIND=JWIM), PARAMETER :: WAM_INTEGRATED_IO_TAG = 334
INTEGER (KIND=JWIM)           :: I, IOSERVNUM, JSPBUFY
INTEGER (KIND=JWIM)           :: IPRM, ICHNK

#include "io_serv_map_send_part1.intfb.h"
#include "io_serv_map_send_part2.intfb.h"
#include "io_serv_log.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('OUTINT_IO_SERV',0,ZHOOK_HANDLE)

      ! -- Create field descriptors and populate with meta-data about the fields
      JFNUM=NIPRMOUT
      ALLOCATE (YLFLDSC (JFNUM))

      IGPTOTG = NIBLO 
      IGPTOT  = NTOTIJ 
      
      DO IFNUM = 1, JFNUM
          
        ! Generic field header data
        YLFLDSC(IFNUM)%CPREF  = 'wv_int'
        YLFLDSC(IFNUM)%NFLSZ  = IGPTOTG
        
        ! Wave model field header data
        ASSOCIATE( WAMFLD => YLFLDSC(IFNUM)%YWAM )
        WAMFLD%IANGLE = -1
        WAMFLD%IFREQ  = -1
        WAMFLD%NDATE_TIME_WINDOW_END = NDATE_TIME_WINDOW_END
        WAMFLD%KCOUSTEP = KCOUSTEP
        WAMFLD%LRSTST0 = LRSTST0
        WAMFLD%ITABLE = INFOBOUT(IFNUM,1)
        WAMFLD%IPARAM = INFOBOUT(IFNUM,2)
        WAMFLD%KLEV =   INFOBOUT(IFNUM,3)
        WAMFLD%CDATE = CDATE
        WAMFLD%IFCST = IFCST
        WAMFLD%MARSTYPE = MARSTYPE
        WAMFLD%NSTEP  = IFSNSTEP
        END ASSOCIATE

      ENDDO

      ! -- Part 1, distributes fields between IO procs and allocates send-buffers in send-plan (YLIOSMPP)
      CALL IO_SERV_MAP_SEND_PART1(  IO_SERV_C001, &
                                 &  IGPTOT, &
                                 &  IGPTOTG, &
                                 &  YLFLDSC, &
                                 &  WAM_INTEGRATED_IO_TAG, &
                                 &  YLIOSMPP, &
                                 &  NIO_SERV_HDR_IDOM_WV, &
                                 &  PTSTEP=IFSTSTEP)

      ! -- Populate buffers
      CALL IOMULTIBUF_INIT_IDX (YLIOSMPP%YLBUFD, IOSERVNUM, JSPBUFY)

      DO IFNUM = 1, JFNUM
        
        CALL IOMULTIBUF_INCR_IDX (YLIOSMPP%YLBUFD, IOSERVNUM, JSPBUFY)
       
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IPRM, I)
        DO ICHNK = 1, NCHNK
          DO IPRM = 1, KIJL4CHNK(ICHNK)
            I = IPRM + NPROMA_WAM*(ICHNK-1)
            YLIOSMPP%YLBUFD( IOSERVNUM)%P(I, JSPBUFY) = BOUT( IPRM, IFNUM, ICHNK )
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

      ENDDO

      IF (IO_SERV_C001%LDBUG) CALL IO_SERV_LOG (IO_SERV_C001, 'outint_io_serv', 3_JWIM)

      ! -- Part 2, finalizes send to IO server
      CALL IO_SERV_MAP_SEND_PART2 (IO_SERV_C001, YLIOSMPP)

      DEALLOCATE (YLFLDSC)

      IF (LHOOK) CALL DR_HOOK('OUTINT_IO_SERV',1,ZHOOK_HANDLE)

END SUBROUTINE OUTINT_IO_SERV
