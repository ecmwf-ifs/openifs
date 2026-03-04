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

SUBROUTINE DISSPEC0(YDGEOMETRY,KIOPROC,KSETV,KCH,KKIND,KSPEC2,KFIELDS,KPOSSP,KPRTRW, &
 & PSPEC,LDC,KNUM)

!**** *DISSPEC0*  - DISTRIBUTION OF A SET OF GLOBAL SPECTRAL FIELDS          

!     Purpose.
!     --------
!       To distribute a set of global spectral fields among processors. 
!       N.B. : all fields in the set are expected to belong to the same V-set   

!**   Interface.
!     ----------
!        *CALL* *DISSPEC0(.....)*

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           - KIOPROC : I/O processor involved
!           - KSETV   : Set-V involved
!           - KCH     : chunk number
!           - KKIND    : 0 -> SEND ; 1 -> RECV
!           - KSPEC2  : spectral dimension (send=>global ; recv=>local)
!           - KFIELDS : number of fields to send/receive
!           - KPOSSP  : location of first spectral coefficient for each set-W
!                       (+ 1 extra value locating the "last + 1" set-W)
!           - KPRTRW  : number of processors for set-W
!         * IN-OUT:
!           - PSPEC   : spectral fields to send/receive
!         * OPTIONAL INPUT:
!           - LDC     : T/F: use this routine for climate/general purpose.
!                       (T: cf. former DISSPEC)
!           - KNUM    : number of the field (required if LDC=T).

!        Implicit arguments :
!        --------------------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      MPP Group *METEO-FRANCE*
!      95-xx-xx : Original

!     Modifications.
!     --------------
!      Modified : 01-07-06 by R. El Khatib : Fix on MPE to MPL conversion
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Ouaraini & R. El Khatib: 03-Oct-2006: partionning of global
!                 spectral data among PEs 
!      Modified : 08-Nov-07 by K. Yessad: merged version of DISSPEC+DISSPEC0
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOMMP0       , ONLY : NPRINTLEV, NPRCIDS, MYPROC, MYSETV, NPROC
USE YOMCT0       , ONLY : LELAM
USE YOMTAG       , ONLY : MTAGDISTSP
USE MPL_MODULE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRTRW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIOPROC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSETV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KKIND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPOSSP(KPRTRW+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPEC(KSPEC2,KFIELDS) 
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LDC
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KNUM

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM)               :: IDIM0GG1(0:YDGEOMETRY%YRDIM%NMSMAX)
INTEGER(KIND=JPIM)               :: IDIM0GG2(0:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB)                  :: ZSPEC2(KSPEC2) 

LOGICAL :: LLBSET,LLREORD,LLC

! ISND    : processor receiving data
! ZBUF    : communication buffer
! IBUFLEN : lenght ofcommunication buffer
INTEGER(KIND=JPIM) :: ILEN, ISENDER, ITAG, ITAGR, JSETW, IBUFLEN, JFLD, IP, ISND
INTEGER(KIND=JPIM) :: JM, JN, IPOS, ISP, ISPG, ISEND

REAL(KIND=JPRB), ALLOCATABLE :: ZBUF(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "set2pe.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DISSPEC0',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDEDIM=>YDGEOMETRY%YREDIM, YDELAP=>YDGEOMETRY%YRELAP)
ASSOCIATE(NMSMAX=>YDDIM%NMSMAX, NSMAX=>YDDIM%NSMAX, &
 & NISNAX=>YDEDIM%NISNAX, &
 & NDIM0G=>YDMP%NDIM0G)

!     ------------------------------------------------------------------

LLC=.FALSE.
IF (PRESENT(LDC)) LLC=LDC

IF (LLC) THEN
  LLBSET=(KSETV==MYSETV)
  LLREORD=.FALSE.
ELSE
  LLBSET=.TRUE.
  LLREORD=.TRUE.
ENDIF

IF (KKIND == 0) THEN

  ! * Reordering (like in SPREORD) when required:

  IF (LLREORD) THEN
    IF (LELAM)  THEN
      IPOS=1
      DO JM=0,NMSMAX
        IDIM0GG1(JM)=IPOS
        IPOS=IPOS+YDELAP%NCPL4M(JM)
      ENDDO
    ELSE
      IPOS=1
      DO JM=0,NSMAX
        IDIM0GG2(JM)=IPOS
        IPOS=IPOS+(NSMAX+1-JM)*2
      ENDDO
    ENDIF   
    IF (LELAM)  THEN
      DO JFLD=1,KFIELDS
        ZSPEC2=0.0_JPRB
        DO JM=0,NMSMAX
          DO JN=0,NISNAX(JM)
            ISP=NDIM0G(JM)+4*JN
            ISPG=IDIM0GG1(JM)+4*JN
            ZSPEC2(ISP:ISP+3)=PSPEC(ISPG:ISPG+3,JFLD)
          ENDDO
        ENDDO
        PSPEC(:,JFLD)=ZSPEC2(:)
      ENDDO
    ELSE
      DO JFLD=1,KFIELDS
        ZSPEC2=0.0_JPRB
        DO JM=0,NSMAX
          DO JN=JM,NSMAX
            ISP=NDIM0G(JM)+2*(JN-JM)
            ISPG=IDIM0GG2(JM)+2*(JN-JM)
            ZSPEC2(ISP)=PSPEC(ISPG,JFLD)
            ZSPEC2(ISP+1)=PSPEC(ISPG+1,JFLD)
          ENDDO
        ENDDO
        PSPEC(:,JFLD)=ZSPEC2(:)
      ENDDO
    ENDIF
  ENDIF

  ! * Distribute fields to nodes:

  DO JSETW=1,KPRTRW
    CALL SET2PE(ISND,0,0,JSETW,KSETV)
    ILEN = KPOSSP(JSETW+1)-KPOSSP(JSETW)
    IF ( ILEN > 0 .AND. NPRCIDS(ISND) /= MYPROC) THEN
      ! There is something to send to somebody
      IBUFLEN=ILEN*KFIELDS
      ALLOCATE(ZBUF(IBUFLEN))
      IP=1
      DO JFLD=1,KFIELDS
        ZBUF(IP:IP+ILEN-1)=PSPEC(KPOSSP(JSETW):KPOSSP(JSETW)+ILEN-1,JFLD)
        IP=IP+ILEN
      ENDDO
      IF (LLC) THEN
        ITAG=MTAGDISTSP+(KCH-1)+KNUM*NPROC+ISND
        ISEND=NPRCIDS(ISND)
      ELSE 
        ITAG=MTAGDISTSP+KCH
        ISEND=ISND
      ENDIF
      CALL MPL_SEND(ZBUF(1:IBUFLEN),KDEST=NPRCIDS(ISEND),KTAG=ITAG, &
       & CDSTRING='DISSPEC0:')  
      IF (NPRINTLEV == 2) THEN
        WRITE(NULOUT,*) 'DISSPEC0 : send spectral fields to ',&
         & NPRCIDS(ISND),' with tag ',ITAG  
      ENDIF
      DEALLOCATE(ZBUF)
    ENDIF
  ENDDO

ELSE

  ! * Receive fields:

  IF (LLBSET .AND. KSPEC2 > 0 .AND. NPRCIDS(KIOPROC) /= MYPROC) THEN
    ! There is something to receive from somebody
    IBUFLEN=KSPEC2*KFIELDS
    ALLOCATE(ZBUF(IBUFLEN))
    IF (LLC) THEN
      ITAG=MTAGDISTSP+(KCH-1)+KNUM*NPROC+MYPROC
    ELSE
      ITAG=MTAGDISTSP+KCH
    ENDIF
    CALL MPL_RECV(ZBUF(1:IBUFLEN),KSOURCE=NPRCIDS(KIOPROC),KTAG=ITAG, &
     & KOUNT=ILEN,KFROM=ISENDER,KRECVTAG=ITAGR,CDSTRING='DISSPEC0:')  
    IF (NPRINTLEV == 2) THEN
      WRITE(NULOUT,*) 'DISSPEC0 : recv spectral field from ',ISENDER,&
       & ' with tag ',ITAGR  
    ENDIF
    IF (ILEN /= IBUFLEN) THEN
      WRITE(NULOUT,*) 'ILEN = ',ILEN,' IBUFLEN = ',IBUFLEN
      CALL ABOR1('DISSPEC0: INVALID RECEIVE MESSAGE LENGHT')
    ELSE
      IP=1
      DO JFLD=1,KFIELDS
        PSPEC(1:KSPEC2,JFLD)=ZBUF(IP:IP+KSPEC2-1)
        IP=IP+KSPEC2
      ENDDO
    ENDIF
    DEALLOCATE(ZBUF)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DISSPEC0',1,ZHOOK_HANDLE)
END SUBROUTINE DISSPEC0
