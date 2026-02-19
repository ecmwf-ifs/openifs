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

SUBROUTINE GATHERGPF(YDGEOMETRY,PLOCGP,PGLOBGP,KFLDS,KROOT,KSCALE)

!**** *GATHERGPF* - Gather global gridpoint fields

!     Purpose.
!     --------
!     Gather global gridpoint fields

!**   Interface.
!     ----------
!        *CALL* *GATHERGPF*(...)

!        Explicit arguments : PLOCGP - local gridpoint fields
!        -------------------- PGLOBGP global gridpoint fields (on KROOT)
!                             KFLDS - number of fields
!                             KROOT - master process

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    - MPL_GATHERV
!                     MPL_ALLGATHERV 

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!       Original : 10-10-02
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!        M.Fisher      17-Mar-2004 Wavelet Jb
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : MYPROC, NPROC
USE MPL_MODULE   , ONLY : MPL_ALLGATHERV, MPL_GATHERV
USE YOMJG        , ONLY : JB_STRUCT
USE YOMCOSJO     , ONLY : OBSCOR_GRID_DEF
IMPLICIT NONE

TYPE(GEOMETRY)    ,TARGET,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDS 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCGP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGLOBGP(:,:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KROOT 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KSCALE 
REAL(KIND=JPRB), ALLOCATABLE :: ZRECVBUF(:)
 
INTEGER(KIND=JPIM) :: IRECVCOUNTS(NPROC),IOFFPROC(NPROC+1)
INTEGER(KIND=JPIM) :: JROC,JFLD,JROF,JGL,JLON
INTEGER(KIND=JPIM) :: IA,IB,IDUM1,IDUM2,IGLOFF,IGL1,IGL2,IOFF,ILAST
INTEGER(KIND=JPIM) :: ILEN,ILOFF,IOFFR,IGPTOTG,IGPTOT

INTEGER(KIND=JPIM), POINTER :: IPONL(:,:), IPSTA(:,:), IPTRFRSTLAT(:), IPFRSTLAT(:), IPLSTLAT(:),&
 & IPLOENG(:), IPGPTOTL(:,:)  
INTEGER(KIND=JPIM) :: II,JJ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "pe2set.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GATHERGPF',0,ZHOOK_HANDLE)


IF (PRESENT(KSCALE)) THEN
  IF (KSCALE /= 9999) THEN
    IPONL       => JB_STRUCT%GRID_DEFINITION(KSCALE)%NONL
    IPSTA       => JB_STRUCT%GRID_DEFINITION(KSCALE)%NSTA
    IPTRFRSTLAT => JB_STRUCT%GRID_DEFINITION(KSCALE)%NPTRFRSTLAT
    IPFRSTLAT   => JB_STRUCT%GRID_DEFINITION(KSCALE)%NFRSTLAT
    IPLSTLAT    => JB_STRUCT%GRID_DEFINITION(KSCALE)%NLSTLAT
    IPLOENG     => JB_STRUCT%GRID_DEFINITION(KSCALE)%NLOENG
    IGPTOT      =  JB_STRUCT%GRID_DEFINITION(KSCALE)%NGPTOT
    IGPTOTG     =  JB_STRUCT%GRID_DEFINITION(KSCALE)%NGPTOTG
    IPGPTOTL    => JB_STRUCT%GRID_DEFINITION(KSCALE)%NGPTOTL
  ELSE
    IPONL       => OBSCOR_GRID_DEF%NONL
    IPSTA       => OBSCOR_GRID_DEF%NSTA
    IPTRFRSTLAT => OBSCOR_GRID_DEF%NPTRFRSTLAT
    IPFRSTLAT   => OBSCOR_GRID_DEF%NFRSTLAT
    IPLSTLAT    => OBSCOR_GRID_DEF%NLSTLAT
    IPLOENG     => OBSCOR_GRID_DEF%NLOENG
    IGPTOT      =  OBSCOR_GRID_DEF%NGPTOT
    IGPTOTG     =  OBSCOR_GRID_DEF%NGPTOTG
    IPGPTOTL    => OBSCOR_GRID_DEF%NGPTOTL
  ENDIF 
ELSE
  IPONL       => YDGEOMETRY%YRMP%NONL
  IPSTA       => YDGEOMETRY%YRMP%NSTA
  IPTRFRSTLAT => YDGEOMETRY%YRMP%NPTRFRSTLAT
  IPFRSTLAT   => YDGEOMETRY%YRMP%NFRSTLAT
  IPLSTLAT    => YDGEOMETRY%YRMP%NLSTLAT
  IPLOENG     => YDGEOMETRY%YRGEM%NLOENG
  IGPTOT      =  YDGEOMETRY%YRGEM%NGPTOT
  IGPTOTG     =  YDGEOMETRY%YRGEM%NGPTOTG
  IPGPTOTL    => YDGEOMETRY%YRGEM%NGPTOTL
ENDIF

ALLOCATE (ZRECVBUF(IGPTOTG*KFLDS))

IF(NPROC > 1) THEN

! Prepare for message passing

  IOFFPROC(1) = 0
  DO JROC=1,NPROC
    CALL PE2SET(JROC,IA,IB,IDUM1,IDUM2)
    IRECVCOUNTS(JROC)= IPGPTOTL(IA,IB)*KFLDS
    IOFFPROC(JROC+1) = IOFFPROC(JROC)+IRECVCOUNTS(JROC)
  ENDDO
  
! Gather field
  CALL GSTATS_BARRIER(769)
  CALL GSTATS(604,0)
  IF(KROOT == -1) THEN
    CALL MPL_ALLGATHERV(PLOCGP,&
     & PRECVBUF=ZRECVBUF,KRECVCOUNTS=IRECVCOUNTS, &
     & CDSTRING='GATHERGPF:')  
  ELSE
    CALL MPL_GATHERV(PLOCGP,KROOT=KROOT,&
     & PRECVBUF=ZRECVBUF,KRECVCOUNTS=IRECVCOUNTS, &
     & CDSTRING='GATHERGPF:')  
  ENDIF
  CALL GSTATS(604,1)
  CALL GSTATS_BARRIER2(769)

! Unpack recieve buffer

  IF(KROOT == -1 .OR. MYPROC == KROOT) THEN
    CALL GSTATS(1428,0)
!$OMP PARALLEL DO SCHEDULE(STATIC)&
!$OMP&PRIVATE(JROC,IA,IB,IDUM1,IDUM2,IGLOFF,IGL1,IGL2,IOFF,ILAST,JJ,&
!$OMP&ILEN,ILOFF,JFLD,IOFFR,JGL,JLON)
    DO JROC=1,NPROC
      CALL PE2SET(JROC,IA,IB,IDUM1,IDUM2)
      IGLOFF = IPTRFRSTLAT(IA)
      IGL1 = IPFRSTLAT(IA)
      IGL2 = IPLSTLAT(IA)
      IOFF = 0
      IF(IA > 1) THEN
        IF( IPLSTLAT(IA-1) == IPFRSTLAT(IA) )THEN
          ILAST = IPLSTLAT(IA-1)-1
        ELSE
          ILAST = IPLSTLAT(IA-1)
        ENDIF
        DO JJ=IPFRSTLAT(1),ILAST
          IOFF = IOFF+IPLOENG(JJ)
        ENDDO
      ENDIF
      DO JFLD=1,KFLDS
        ILEN = 0
        ILOFF = 0
        IOFFR = IOFFPROC(JROC)+(JFLD-1)*IPGPTOTL(IA,IB)
        DO JGL=IGL1,IGL2
          DO JLON=1,IPONL(IGLOFF+JGL-IGL1,IB)
            PGLOBGP(IOFF+ILOFF+IPSTA(IGLOFF+JGL-IGL1,IB)+JLON-1,JFLD) &
             & = ZRECVBUF(IOFFR+ILEN+JLON)  
          ENDDO
          ILEN = ILEN + IPONL(IGLOFF+JGL-IGL1,IB)
          ILOFF = ILOFF + IPLOENG(JGL)
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    CALL GSTATS(1428,1)
  ENDIF

ELSE ! If only one PE copy

  CALL GSTATS(1428,0)
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JFLD,JROF,II)
  DO JFLD=1,KFLDS
    DO JROF=1,IGPTOT
      II=IGPTOT*(JFLD-1)+JROF
      PGLOBGP(JROF,JFLD) = PLOCGP(II)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1428,1)
ENDIF

DEALLOCATE (ZRECVBUF)
IF (LHOOK) CALL DR_HOOK('GATHERGPF',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE GATHERGPF
