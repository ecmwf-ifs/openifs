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

MODULE DISGRID_MOD

!     Purpose.
!     --------
!       To scatter a set of global gridpoint fields among processors.

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*
!        Original : 31-Mar-2010 from DISGRID

!     Modifications.
!     --------------
!        17/10/2011 : P.Marguinaud : Fix DISGRID_RECV intent
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!      R. El Khatib 16-May-2019 optimize memory access in NGPSET2PE
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULERR
USE YOMMP0   , ONLY : NPRINTLEV, MYPROC, NPRCIDS, NGPSET2PE, N_REGIONS, N_REGIONS_NS
USE YOMTAG   , ONLY : MTAGDISTGP
USE MPL_MODULE, ONLY: MPL_SEND, MPL_RECV

PRIVATE
PUBLIC DISGRID_SEND, DISGRID_RECV

CONTAINS
!=============================================================
SUBROUTINE DISGRID_SEND(YDGEOMETRY,KNUM,PREALG,KCH,PREAL)

!     Purpose.
!     --------
!       To scatter a set of global gridpoint fields among processors.

!        Explicit arguments :
!        --------------------
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - KCH     : chunk number
!           - PREALG  : global gridpoint fields
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUM 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREALG(YDGEOMETRY%YRGEM%NGPTOTG,KNUM)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCH
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREAL(YDGEOMETRY%YRGEM%NGPTOT,KNUM)

!     ------------------------------------------------------------------


INTEGER(KIND=JPIM) ::  IGL1, IGL2, IGLOFF, ILAST, ILEN, ILENB, IOFF, ISND
INTEGER(KIND=JPIM) ::  ITAG, JA, JB, JGL, J
REAL(KIND=JPRB), ALLOCATABLE :: ZFLD(:,:) ! communication buffer
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DISGRID_MOD:DISGRID_SEND',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPTOTG=>YDGEM%NGPTOTG, NLOENG=>YDGEM%NLOENG, NGPTOT=>YDGEM%NGPTOT, &
 & NFRSTLAT=>YDMP%NFRSTLAT, NONL=>YDMP%NONL, &
 & NPTRFRSTLAT=>YDMP%NPTRFRSTLAT, NLSTLAT=>YDMP%NLSTLAT)
!     ------------------------------------------------------------------

DO JA=1,N_REGIONS_NS
  IGLOFF=NPTRFRSTLAT(JA)
  IGL1  = NFRSTLAT(JA)
  IGL2  = NLSTLAT(JA)
  IOFF=0
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
    ISND=NGPSET2PE(JB,JA)
    ITAG=MTAGDISTGP+KCH
    ILENB=ILEN*KNUM
    ALLOCATE(ZFLD(ILEN,KNUM))
    IF (ILENB > 0 .AND. ISND /= MYPROC) THEN
      CALL DISGRID_SETTLE(YDGEOMETRY,IGLOFF,IGL1,IGL2,JB,IOFF,KNUM,ILEN,PREALG,ZFLD)
      CALL MPL_SEND(ZFLD,KDEST=NPRCIDS(ISND),KTAG=ITAG,CDSTRING='DISGRID_MOD:DISGRID_SEND')  
      IF (NPRINTLEV == 2) THEN
        WRITE(NULOUT,*) 'DISGRID_SEND : ', MYPROC,  ' send gridpoint field to ',NPRCIDS(ISND),&
       & ' with tag ',ITAG
      ENDIF
    ELSE
      IF (ILEN /= NGPTOT) THEN
        WRITE(NULOUT,*) 'ILEN = ',ILEN,' NGPTOT = ',NGPTOT
        CALL ABOR1('DISGRID_SEND: INVALID COMPUTED LOCAL LENGHT')
      ENDIF
      CALL DISGRID_SETTLE(YDGEOMETRY,IGLOFF,IGL1,IGL2,JB,IOFF,KNUM,NGPTOT,PREALG,PREAL)
    ENDIF
    DEALLOCATE(ZFLD)
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DISGRID_MOD:DISGRID_SEND',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DISGRID_SEND

!=============================================================
SUBROUTINE DISGRID_RECV(YDGEOMETRY,KIOPROC,KNUM,PREAL,KCH)

!     Purpose.
!     --------
!       To scatter a set of global gridpoint fields among processors.

!        Explicit arguments :
!        --------------------
!           - KIOPROC : I/O processor involved
!           - KNUM    : number of fields to send/receive
!           - PREAL   : local gridpoint fields
!           - KCH     : chunk number
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KIOPROC 
INTEGER(KIND=JPIM),INTENT(IN) :: KNUM 
INTEGER(KIND=JPIM),INTENT(IN) :: KCH
REAL(KIND=JPRB),INTENT(INOUT) :: PREAL(YDGEOMETRY%YRGEM%NGPTOT,KNUM)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILEN, ITAG, ISENDER, ITAGR, ILENB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DISGRID_MOD:DISGRID_RECV',0,ZHOOK_HANDLE)
ASSOCIATE(NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT)
!     ------------------------------------------------------------------

!     send fields to output node  :  KIOPROC

ILEN=NGPTOT*KNUM
IF (ILEN > 0 .AND. KIOPROC /= MYPROC) THEN
  ITAG=MTAGDISTGP+KCH
  CALL MPL_RECV(PREAL,KSOURCE=NPRCIDS(KIOPROC),KTAG=ITAG,KOUNT=ILENB,&
   & KFROM=ISENDER,KRECVTAG=ITAGR,CDSTRING='DISGRID_MOD:DISGRID_RECV')  
  IF (NPRINTLEV == 2) THEN
    WRITE(NULOUT,*) 'DISGRID_RECV : ', MYPROC, ' recv gridpoint field from ',ISENDER,&
     & ' with tag ',ITAGR
  ENDIF
  IF (ILENB /= ILEN) THEN
    WRITE (NULERR, *) " ILENB = ", ILENB, " ILEN = ", ILEN,&
                     &" MYPROC = ", MYPROC, " KSOURCE = ", NPRCIDS(KIOPROC),&
                     &" KFROM = ", ISENDER, " KRECVTAG = ", ITAGR, " KTAG = ", ITAG
    CALL FLUSH (NULERR)
    CALL ABOR1('DISGRID: INVALID RECEIVE MESSAGE LENGTH')
  ENDIF
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DISGRID_MOD:DISGRID_RECV',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DISGRID_RECV

!=============================================================
SUBROUTINE DISGRID_SETTLE(YDGEOMETRY,KGLOFF,KGL1,KGL2,KB,KOFF,KNUM,KGPTOT,PREALG,PREAL)

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

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KGLOFF,KGL1,KGL2,KB,KOFF,KNUM,KGPTOT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREALG(YDGEOMETRY%YRGEM%NGPTOTG,KNUM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREAL(KGPTOT,KNUM)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) ::  JFLD, JGL, JLON, II, ILOFF
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DISGRID_MOD:DISGRID_SETTLE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NLOENG=>YDGEM%NLOENG, NONL=>YDMP%NONL, NSTA=>YDMP%NSTA)
!     ------------------------------------------------------------------

DO JFLD=1,KNUM
  II=0
  ILOFF=0
  DO JGL=KGL1,KGL2
    DO JLON=1,NONL(KGLOFF+JGL-KGL1,KB)
      PREAL(II+JLON,JFLD)=PREALG(KOFF+ILOFF+NSTA(KGLOFF+JGL-KGL1,KB)+JLON-1,JFLD)
    ENDDO
    II=II+NONL(KGLOFF+JGL-KGL1,KB)
    ILOFF=ILOFF + NLOENG(JGL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DISGRID_MOD:DISGRID_SETTLE',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE DISGRID_SETTLE

!=============================================================
END MODULE DISGRID_MOD
