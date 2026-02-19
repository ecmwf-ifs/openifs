! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE DIWRGRID_MOD

!     Purpose.
!     --------
!       To gather a set of global gridpoint fields from processors.

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*
!        Original : 31-Mar-2010 from DIWGRID

!     Modifications.
!     --------------
!        P. Marguinaud : 01-Jan-2011; add DIWRGRID_RECVX to handle fields
!                        with a header 
!        0. Vignes : 08-Jun-2010 Added GATHERV option (from S. Saarinen)
!        P. Marguinaud : 05-Jul-2011 Cleaning
!        P. Marguinaud : 11-Sep-2012 Use IOMULTIBUF_MOD to avoid unnecessary 
!                                    copies of arrays
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!      R. El Khatib 16-May-2019 optimize memory access in NGPSET2PE
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : NPRINTLEV, MYPROC, NPRCIDS, NGPSET2PE, N_REGIONS, N_REGIONS_NS
USE YOMTAG   , ONLY : MTAGDISTGP
USE YOMGEM   , ONLY : TGEM
USE MPL_MODULE, ONLY: MPL_SEND, MPL_RECV, MPL_GATHERV
USE IOMULTIBUF_MOD, ONLY : IOMULTIBUF,&
                         & IOMULTIBUF_SIZE_IDX,&
                         & IOMULTIBUF_COMP_IDX

PRIVATE
PUBLIC DIWRGRID_SEND, DIWRGRID_RECV, DIWRGRID_RECVX

INTERFACE DIWRGRID_RECV
  MODULE PROCEDURE DIWRGRID_RECV1, DIWRGRID_RECV2, DIWRGRID_RECVX, DIWRGRID_RECVY
END INTERFACE 

CONTAINS
!=============================================================
SUBROUTINE DIWRGRID_SEND(YDGEM,KIOPROC,KNUM,PREAL,KCH,LD_GATHERV)

!     Purpose.
!     --------
!       To distribute a set of global gridpoint fields among processors.

!        Explicit arguments :
!        --------------------
!           - KIOPROC : I/O processor involved
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - KCH     : chunk number
!           - LD_GATHERV : whether to use MPL_GATHERV (optional)
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGEM)        ,INTENT(IN) :: YDGEM
INTEGER(KIND=JPIM),INTENT(IN) :: KIOPROC 
INTEGER(KIND=JPIM),INTENT(IN) :: KNUM 
INTEGER(KIND=JPIM),INTENT(IN) :: KCH
REAL(KIND=JPRB),INTENT(IN)    :: PREAL(YDGEM%NGPTOT,KNUM)
LOGICAL,OPTIONAL,INTENT(IN)   :: LD_GATHERV

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILEN, ITAG
LOGICAL :: LL_GATHERV
REAL(KIND=JPRB), ALLOCATABLE :: ZBUF(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_SEND',0,ZHOOK_HANDLE)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT)
!     ------------------------------------------------------------------

LL_GATHERV = .FALSE.
IF (PRESENT(LD_GATHERV)) LL_GATHERV = LD_GATHERV

!     send fields to output node  :  KIOPROC

ILEN=NGPTOT*KNUM
IF (ILEN > 0 .AND. KIOPROC /= MYPROC) THEN
  IF (LL_GATHERV) THEN
    ! Should really make a version of MPL_GATHERV accepting 2D fields
    ALLOCATE(ZBUF(ILEN))
    ZBUF(:) = RESHAPE(PREAL,SHAPE(ZBUF))
    CALL MPL_GATHERV(ZBUF,NPRCIDS(KIOPROC),CDSTRING='DIWRGRID_MOD:DIWRGRID_SEND')
    IF (NPRINTLEV == 2) THEN
      WRITE(NULOUT,*) 'DIWRGRID_SEND : gatherv g.p. field to ',NPRCIDS(KIOPROC)
    ENDIF
    DEALLOCATE(ZBUF)
  ELSE
    ITAG=MTAGDISTGP+KCH
    CALL MPL_SEND(PREAL,KDEST=NPRCIDS(KIOPROC),KTAG=ITAG,CDSTRING='DIWRGRID_MOD:DIWRGRID_SEND')  
    IF (NPRINTLEV == 2) THEN
      WRITE(NULOUT,*) 'DIWRGRID_SEND :  send gridpoint field to ',NPRCIDS(KIOPROC),&
       & ' with tag ',ITAG, ' KNUM = ', KNUM, ' NGPTOT = ', NGPTOT
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_SEND',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DIWRGRID_SEND

!=============================================================
SUBROUTINE DIWRGRID_RECV1(YDGEOMETRY,KNUM,PREAL,KCH,PREALG,LD_GATHERV)

!     Purpose.
!     --------
!       To distribute a set of global gridpoint fields among processors.

!        Explicit arguments :
!        --------------------
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - KCH     : chunk number
!           - PREALG  : global gridpoint fields
!           - LD_GATHERV : whether to use MPL_GATHERV (optional)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUM 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREAL(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCH
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREALG(YDGEOMETRY%YRGEM%NGPTOTG)
LOGICAL,OPTIONAL,INTENT(IN)   :: LD_GATHERV

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECV1',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

CALL DIWRGRID_RECVX(YDGEOMETRY,KNUM,PREAL,KCH,PREALG,0,LD_GATHERV)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECV1',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DIWRGRID_RECV1

!=============================================================
SUBROUTINE DIWRGRID_RECV2(YDGEOMETRY,KNUM,PREAL,KCH,PREALG,LD_GATHERV)

!     Purpose.
!     --------
!       To distribute a set of global gridpoint fields among processors.

!        Explicit arguments :
!        --------------------
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - KCH     : chunk number
!           - PREALG  : global gridpoint fields
!           - LD_GATHERV : whether to use MPL_GATHERV (optional)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUM 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREAL(YDGEOMETRY%YRGEM%NGPTOT,KNUM)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCH
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREALG(YDGEOMETRY%YRGEM%NGPTOTG,KNUM)
LOGICAL,OPTIONAL,INTENT(IN)   :: LD_GATHERV

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECV2',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

CALL DIWRGRID_RECVX(YDGEOMETRY,KNUM,PREAL,KCH,PREALG,0,LD_GATHERV)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECV2',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DIWRGRID_RECV2

!=============================================================
SUBROUTINE DIWRGRID_RECVX(YDGEOMETRY,KNUM,PREAL,KCH,PREALG,KHDRS,LD_GATHERV)

!     Purpose.
!     --------
!       To distribute a set of global gridpoint fields among processors.

!        Explicit arguments :
!        --------------------
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - KCH     : chunk number
!           - PREALG  : global gridpoint fields
!           - KHDRS   : size of header if any (optional)
!           - LD_GATHERV : whether to use MPL_GATHERV (optional)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)          :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)          :: KNUM 
REAL(KIND=JPRB)   ,INTENT(IN)          :: PREAL(YDGEOMETRY%YRGEM%NGPTOT,KNUM)
INTEGER(KIND=JPIM),INTENT(IN)          :: KCH
INTEGER(KIND=JPIM),INTENT(IN)          :: KHDRS
REAL(KIND=JPRB)   ,INTENT(OUT), TARGET :: PREALG(KHDRS+YDGEOMETRY%YRGEM%NGPTOTG,KNUM)
LOGICAL,OPTIONAL  ,INTENT(IN)          :: LD_GATHERV

TYPE (IOMULTIBUF) :: YLGPL (1)

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECVX',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

YLGPL(1)%P => PREALG (KHDRS+1:KHDRS+YDGEOMETRY%YRGEM%NGPTOTG,:)

CALL DIWRGRID_RECVY(YDGEOMETRY,KNUM, PREAL, KCH, YLGPL, LD_GATHERV)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECVX',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DIWRGRID_RECVX

!=============================================================
SUBROUTINE DIWRGRID_RECVY(YDGEOMETRY,KNUM,PREAL,KCH,YDGPL,LD_GATHERV)

!     Purpose.
!     --------
!       To distribute a set of global gridpoint fields among processors.

!        Explicit arguments :
!        --------------------
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - KCH     : chunk number
!           - PREALG  : global gridpoint fields
!           - LD_GATHERV : whether to use MPL_GATHERV (optional)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUM 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREAL(YDGEOMETRY%YRGEM%NGPTOT,KNUM)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCH
TYPE (IOMULTIBUF) ,INTENT(INOUT) :: YDGPL (:)
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LD_GATHERV

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) ::  IGL1, IGL2, IGLOFF, ILAST, ILEN, ILENB, IOFF, IRCV, IFNUM
INTEGER(KIND=JPIM) ::  IBUFLEN, ISENDER, ITAG, ITAGR, IOFFZ, JA, JB, JGL, J
INTEGER(KIND=JPIM), ALLOCATABLE :: IRECVCOUNTS(:)
LOGICAL :: LL_GATHERV
REAL(KIND=JPRB), ALLOCATABLE :: ZFLD(:) ! communication buffer
REAL(KIND=JPRB), ALLOCATABLE :: ZBUF(:)
INTEGER (KIND=JPIM) :: IY (KNUM), JY (KNUM)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECVY',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPTOTG=>YDGEM%NGPTOTG, NLOENG=>YDGEM%NLOENG, NGPTOTL=>YDGEM%NGPTOTL, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & NFRSTLAT=>YDMP%NFRSTLAT, NONL=>YDMP%NONL, &
 & NPTRFRSTLAT=>YDMP%NPTRFRSTLAT, NLSTLAT=>YDMP%NLSTLAT)
!     ------------------------------------------------------------------

CALL IOMULTIBUF_SIZE_IDX (YDGPL, IFNUM)

IF (IFNUM /= KNUM)&
  & CALL ABOR1 ('DIWRGRID_RECVY: FIELD NUMBER MISMATCH')

CALL IOMULTIBUF_COMP_IDX (YDGPL, IY, JY)

LL_GATHERV = .FALSE.
IF (PRESENT(LD_GATHERV)) LL_GATHERV = LD_GATHERV

IF (LL_GATHERV) THEN
  ALLOCATE(IRECVCOUNTS(1:SUM(N_REGIONS(1:N_REGIONS_NS))))
  ALLOCATE(ZFLD(NGPTOTG*KNUM))
  ALLOCATE(ZBUF(NGPTOT*KNUM))
  ZBUF(:) = RESHAPE(PREAL,SHAPE(ZBUF))
  DO JA=1,N_REGIONS_NS
    DO JB=1,N_REGIONS(JA)
      IRCV=NGPSET2PE(JB,JA)
      IRECVCOUNTS(IRCV) = NGPTOTL(JA,JB)*KNUM
    ENDDO
  ENDDO
  CALL MPL_GATHERV(ZBUF,NPRCIDS(MYPROC),PRECVBUF=ZFLD,&
   &               KRECVCOUNTS=IRECVCOUNTS,CDSTRING='DIWRGRID_MOD:DIWRGRID_RECVY')
  IF (NPRINTLEV == 2) THEN
    WRITE(NULOUT,*) 'DIWRGRID_RECV : gatherv recv g.p. fields'
  ENDIF
ENDIF

DO JA = 1, N_REGIONS_NS

  IGLOFF = NPTRFRSTLAT(JA)
  IGL1   = NFRSTLAT(JA)
  IGL2   = NLSTLAT(JA)
  IOFF   = 0

  IF (JA > 1) THEN
    IF (NLSTLAT(JA-1) == NFRSTLAT(JA)) THEN
      ILAST=NLSTLAT(JA-1)-1
    ELSE
      ILAST=NLSTLAT(JA-1)
    ENDIF
    DO J=NFRSTLAT(1),ILAST
      IOFF=IOFF+NLOENG(J)
    ENDDO
  ENDIF

  DO JB=1,N_REGIONS(JA)

    ILEN=0
    DO JGL=IGL1,IGL2
      ILEN=ILEN+NONL(IGLOFF+JGL-IGL1,JB)
    ENDDO

    IRCV=NGPSET2PE(JB,JA)
    ILENB=ILEN*KNUM

    IF (LL_GATHERV) THEN
      IOFFZ = SUM(IRECVCOUNTS(1:IRCV-1)) + 1
      IF (ILENB /= IRECVCOUNTS(IRCV)) THEN
        WRITE(NULOUT,*) 'IRCV = ',IRCV, ' ILENB = ',ILENB,&
         &              ' IBUFLEN = ',IRECVCOUNTS(IRCV)
        CALL ABOR1('DIWRGRID_RECV: INVALID GATHERV RECEIVE MESSAGE LENGTH')
      ENDIF
      CALL DIWRGRID_SETTLE(YDGEOMETRY,IGLOFF, IGL1, IGL2, JB, IOFF, KNUM, ILEN,&
       &                   ZFLD(IOFFZ), YDGPL, IY, JY)
    ELSE
      ITAG=MTAGDISTGP+KCH
      ALLOCATE(ZFLD(ILENB))
      IF (ILENB > 0 .AND. IRCV /= MYPROC) THEN
        CALL MPL_RECV(ZFLD,KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
         & KOUNT=IBUFLEN,KFROM=ISENDER,KRECVTAG=ITAGR,CDSTRING='DIWRGRID_MOD:DIWRGRID_RECVY')  
        IF (ILENB /= IBUFLEN) THEN
          WRITE(NULOUT,*) 'KSOURCE = ', NPRCIDS(IRCV), 'KFROM = ', ISENDER
          WRITE(NULOUT,*) 'ILENB = ',ILENB,' IBUFLEN = ',IBUFLEN
          CALL ABOR1('DIWRGRID_RECV: INVALID RECEIVE MESSAGE LENGTH')
        ENDIF 
        IF (NPRINTLEV == 2) THEN
          WRITE(NULOUT,*) 'DIWRGRID_RECV : recv gridpoint field from ',ISENDER,&
         & ' with tag ',ITAGR, ' KNUM = ', KNUM
        ENDIF
        CALL DIWRGRID_SETTLE(YDGEOMETRY,IGLOFF, IGL1, IGL2, JB, IOFF, KNUM, ILEN, ZFLD, YDGPL, IY, JY)
      ELSE
        IF (ILEN /= NGPTOT) THEN
          WRITE(NULOUT,*) 'ILEN = ',ILEN,' NGPTOT = ',NGPTOT
          CALL ABOR1('DIWRGRID_RECV: INVALID COMPUTED LOCAL LENGTH')
        ENDIF
        CALL DIWRGRID_SETTLE(YDGEOMETRY,IGLOFF, IGL1, IGL2, JB, IOFF, KNUM, NGPTOT, PREAL, YDGPL, IY, JY)
      ENDIF
      DEALLOCATE(ZFLD)
    ENDIF

  ENDDO
ENDDO

IF (LL_GATHERV) THEN
  DEALLOCATE(IRECVCOUNTS)
  DEALLOCATE(ZFLD)
  DEALLOCATE(ZBUF)
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_RECVY',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DIWRGRID_RECVY

!=============================================================
SUBROUTINE DIWRGRID_SETTLE(YDGEOMETRY,KGLOFF, KGL1, KGL2, KB, KOFF, KNUM,&
                          & KGPTOT, PREAL, YDGPL, KIY, KJY)

!     Purpose.
!     --------
!       To settle into order the gridpoint contribution from each proc.

!        Explicit arguments :
!        --------------------
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - PREALG  : global gridpoint fields
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)     , INTENT(IN)    :: YDGEOMETRY
INTEGER (KIND=JPIM), INTENT(IN)    :: KGLOFF,KGL1,KGL2,KB,KOFF,KNUM,KGPTOT
REAL (KIND=JPRB)   , INTENT(IN)    :: PREAL(KGPTOT,KNUM)
TYPE (IOMULTIBUF)  , INTENT(INOUT) :: YDGPL (:)
INTEGER (KIND=JPIM), INTENT (IN)   :: KIY (:), KJY (:)

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JFLD, JGL, JLON, II, ILOFF
INTEGER(KIND=JPIM) :: IY, JY

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_SETTLE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NLOENG=>YDGEM%NLOENG, NONL=>YDMP%NONL, NSTA=>YDMP%NSTA)
!     ------------------------------------------------------------------

DO JFLD =1, KNUM
 
  IY = KIY (JFLD)
  JY = KJY (JFLD)

  II=0
  ILOFF=0

  DO JGL=KGL1,KGL2
    DO JLON=1,NONL(KGLOFF+JGL-KGL1,KB)
      YDGPL(IY)%P(KOFF+ILOFF+NSTA(KGLOFF+JGL-KGL1,KB)+JLON-1,JY) =&
       & PREAL(II+JLON,JFLD)  
    ENDDO
    II=II+NONL(KGLOFF+JGL-KGL1,KB)
    ILOFF=ILOFF + NLOENG(JGL)
  ENDDO

ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DIWRGRID_MOD:DIWRGRID_SETTLE',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DIWRGRID_SETTLE

!=============================================================
END MODULE DIWRGRID_MOD
