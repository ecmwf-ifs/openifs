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

SUBROUTINE SUWAM(YDGEOMETRY,YDEWCOU)
#ifdef WITH_WAVE

!     PURPOSE.
!     --------

!           THIS ROUTINE CALLS THE WAVE MODEL TO OBTAIN THE CONTRIBUTION

!     REFERENCES.
!     -----------

!           NONE.

!  ADD EXTRA FIELDS CODE TO COMMON BLOCK, CHANGE CALLPAR AND THIS 
!  ROUTINE

!  3 FLAVOURS FOR ZBETA, LABELLED ZBETAXX, WHERE XX STANDS FOR

!          WN  - NON-EXTENDED WAVE MODEL GRID
!          WX  - POLE-EXTENDED WAVE MODEL GRID
!          RG  - ATM REDUCED GAUSSIAN GRID

!     AUTHOR.
!     -------
!      P. VITERBO               E.C.M.W.F.      06/09/88

!     MODIFICATIONS.
!     --------------
!      O. la (ECMWF) 27/11/01: Generalized ATM-WAVE interface
!     ------------------------------------------------------------------


USE GEOMETRY_MOD, ONLY : GEOMETRY
USE PARKIND1    , ONLY : JPIM, JPRB
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0      , ONLY : LFDBOP
USE YOMLUN      , ONLY : NULOUT
USE YOEWCOU     , ONLY : TEWCOU
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TEWCOU)      ,INTENT(INOUT) :: YDEWCOU
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IGPTOTG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL :: LLRNL

!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUWAM',0,ZHOOK_HANDLE)
ASSOCIATE(NDGLG=>YDGEOMETRY%YRDIM%NDGLG, &
 & LWCOU=>YDEWCOU%LWCOU, LWCOU2W=>YDEWCOU%LWCOU2W, &
 & LWCOURNW=>YDEWCOU%LWCOURNW, LWCOUHMF=>YDEWCOU%LWCOUHMF, &
 & LWCOUNORMS=>YDEWCOU%LWCOUNORMS, LWFLUX=>YDEWCOU%LWFLUX, &
 & LWFRSTIME=>YDEWCOU%LWFRSTIME, &
 & LWVIN_MASK_NOT_SET=>YDEWCOU%LWVIN_MASK_NOT_SET, &
 & LWVIN_UNINITIALISED=>YDEWCOU%LWVIN_UNINITIALISED, &
 & LWVOUT_MASK_NOT_SET=>YDEWCOU%LWVOUT_MASK_NOT_SET, &
 & NLAT1W=>YDEWCOU%NLAT1W, &
 & NLATW=>YDEWCOU%NLATW, NLON1W=>YDEWCOU%NLON1W, NLONW=>YDEWCOU%NLONW, &
 & NNORXW=>YDEWCOU%NNORXW, NWV_W2IWGHT=>YDEWCOU%NWV_W2IWGHT, &
 & RDEGREW=>YDEWCOU%RDEGREW, RMISSW=>YDEWCOU%RMISSW, RNORTW=>YDEWCOU%RNORTW, &
 & RSOUTW=>YDEWCOU%RSOUTW, &
 & NGPTOTG=>YDGEOMETRY%YRGEM%NGPTOTG)
!     ------------------------------------------------------------------

IGPTOTG=NGPTOTG


WRITE(NULOUT,*) 'SETTING UP THE WAVE MODEL FOR TWO-WAY INTERACTION BETWEEN WIND AND WAVES '

!*    ------------------------------------------------------------------
!*       0.     WAVE MODEL SETUP 
!               ----------------

IF (.NOT.ALLOCATED(YDEWCOU%MASK_WAVE_OUT) ) THEN

  LLRNL=.TRUE.
  CALL WVWAMINIT(LWCOU,NULOUT,LLRNL,NLONW,NLATW,RSOUTW,RNORTW)

  RDEGREW=(RNORTW-RSOUTW)/(NLATW-1)

! Allocate and initialise masks for optimised coupling communications
  ALLOCATE(YDEWCOU%MASK_WAVE_IN(IGPTOTG))
  ALLOCATE(YDEWCOU%MASK_WAVE_OUT(NLONW,NLATW))
  YDEWCOU%MASK_WAVE_IN(:)=0
  YDEWCOU%MASK_WAVE_OUT(:,:)=0
! Set flags to indicate mask and comms data structures are not initialised
  LWVIN_MASK_NOT_SET=.TRUE.
  LWVOUT_MASK_NOT_SET=.TRUE.
  LWVIN_UNINITIALISED=.TRUE.
! Inactive until WV_W2IWGHT is defined.
  NWV_W2IWGHT=0
! Set default for producing global norms
  LWCOUNORMS=.FALSE.

ENDIF


!IF(LWFRSTIME) THEN
   CALL WVWAMINIT1(LWCOU, LWCOU2W, LWCOURNW, LWCOUHMF, LWFLUX, LFDBOP)
   CALL WVWAMDECOMP
!ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUWAM',1,ZHOOK_HANDLE)
#endif
END SUBROUTINE SUWAM
