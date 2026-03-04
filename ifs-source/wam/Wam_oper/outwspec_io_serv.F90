! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE OUTWSPEC_IO_SERV (IJS, IJL, SPEC, MARSTYPE, CDATE, IFCST)

!----------------------------------------------------------------------

!**** *OUTWSPEC*  WRITES SPECTRA AS PARAMETER 251 USING IO SERVER

!     J. HAWKES   ECMWF  OCTOBER 2017 

!     See also: wv_io_serv_rec
!     See also: outint_io_serv

!-------------------------------------------------------------------

! Wave Model
USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWMPP   , ONLY : IRANK    ,NPROC
USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NGX      ,NGY
USE YOWSPEC  , ONLY : NSTART   ,NEND
USE YOWTEST  , ONLY : IU06     ,ITEST
USE YOWGRIBHD, ONLY : NDATE_TIME_WINDOW_END

! IFS/Coupling
USE MPL_MODULE
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOWCOUT  , ONLY : LRSTST0
USE YOWCOUP  , ONLY : KCOUSTEP, IFSNSTEP, IFSTSTEP

! IO Server
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE YOMIO_SERV, ONLY : IO_SERV_C001
USE YOMIO_SERV_MAP_PLAN, ONLY : IO_SERV_SEND_PLAN
USE YOMIO_SERV_HDR, ONLY : NIO_SERV_HDR_IDOM_WV
USE IOMULTIBUF_MOD, ONLY : IOMULTIBUF_INCR_IDX, IOMULTIBUF_INIT_IDX

!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: SPEC
CHARACTER(LEN=2), INTENT(IN)    :: MARSTYPE
CHARACTER(LEN=14), INTENT(IN)   :: CDATE
INTEGER, INTENT(IN)             :: IFCST

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

TYPE (IO_SERV_SEND_PLAN)      :: YLIOSMPP          ! Send map plan
TYPE (IOFLDDESC), ALLOCATABLE :: YLFLDSC (:)       ! Field descriptor array
INTEGER (KIND=JWIM)           :: IFNUM, JFNUM      ! Number of fields
INTEGER (KIND=JWIM)           :: IGPTOTG, IGPTOT   ! Number of points (global & local)
INTEGER (KIND=JWIM), PARAMETER :: WAM_SPECTRAL_IO_TAG = 333
INTEGER (KIND=JWIM), PARAMETER :: DUMMY = 0
INTEGER (KIND=JWIM)           :: I, IANGLE, IFREQ, IOSERVNUM, JSPBUFY
INTEGER (KIND=JWIM)           :: IPARAM, ITABLE, KLEV

#include "io_serv_map_send_part1.intfb.h"
#include "io_serv_map_send_part2.intfb.h"
#include "io_serv_log.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV',0,ZHOOK_HANDLE)

      IPARAM = 251
      ITABLE = 140
      KLEV = 0
      
      ! -- Create field descriptors and populate with meta-data about the fields
      JFNUM=NFRE_RED*NANG
      ALLOCATE (YLFLDSC (JFNUM))

      IGPTOTG = NEND(NPROC) - NSTART(1) + 1
      IGPTOT  = NEND(IRANK) - NSTART(IRANK) + 1

      IFNUM = 1
      DO IFREQ = 1, NFRE_RED
        DO IANGLE = 1, NANG
          
          ! Generic field header data
          YLFLDSC(IFNUM)%CPREF  = 'wv_spec'
          YLFLDSC(IFNUM)%NFLSZ  = IGPTOTG
          
          ! Wave model field header data
          ASSOCIATE( WAMFLD => YLFLDSC(IFNUM)%YWAM )
          WAMFLD%IANGLE = IANGLE
          WAMFLD%IFREQ  = IFREQ
          WAMFLD%NDATE_TIME_WINDOW_END = NDATE_TIME_WINDOW_END
          WAMFLD%KCOUSTEP = KCOUSTEP
          WAMFLD%LRSTST0 = LRSTST0
          WAMFLD%ITABLE = ITABLE
          WAMFLD%IPARAM = IPARAM
          WAMFLD%KLEV = KLEV
          WAMFLD%CDATE = CDATE
          WAMFLD%IFCST = IFCST
          WAMFLD%MARSTYPE = MARSTYPE
          WAMFLD%NSTEP  = IFSNSTEP
          END ASSOCIATE
          IFNUM = IFNUM + 1
        ENDDO
      ENDDO

      ! -- Part 1, distributes fields between IO procs and allocates send-buffers in send-plan (YLIOSMPP)
      CALL IO_SERV_MAP_SEND_PART1(  IO_SERV_C001, &
                                 &  IGPTOT, &
                                 &  IGPTOTG, &
                                 &  YLFLDSC, &
                                 &  WAM_SPECTRAL_IO_TAG, &
                                 &  YLIOSMPP, &
                                 &  NIO_SERV_HDR_IDOM_WV, &
                                 &  PTSTEP=IFSTSTEP)

      ! -- Populate buffers
      CALL IOMULTIBUF_INIT_IDX (YLIOSMPP%YLBUFD, IOSERVNUM, JSPBUFY)

      DO IFNUM = 1, JFNUM
        
        CALL IOMULTIBUF_INCR_IDX (YLIOSMPP%YLBUFD, IOSERVNUM, JSPBUFY)
        
        IFREQ  = YLFLDSC (IFNUM)%YWAM%IFREQ
        IANGLE = YLFLDSC (IFNUM)%YWAM%IANGLE

        DO I = 1, IGPTOT
            YLIOSMPP%YLBUFD(IOSERVNUM)%P(I,JSPBUFY) = SPEC( NSTART(IRANK) + I - 1, IANGLE, IFREQ)
        ENDDO

      ENDDO

      ! -- Part 2, finalizes send to IO server
      CALL IO_SERV_MAP_SEND_PART2 (IO_SERV_C001, YLIOSMPP)

      IF (IO_SERV_C001%LDBUG) CALL IO_SERV_LOG (IO_SERV_C001, 'outwspec_io_serv', 3_JWIM)

      DEALLOCATE (YLFLDSC)

      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV',1,ZHOOK_HANDLE)

END SUBROUTINE OUTWSPEC_IO_SERV
