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

SUBROUTINE SPAREORD(YDDIM,KFLEV,PSPFILE,PSPBUF,LD_FILE_TO_MODEL)

USE YOMDIM   , ONLY : TDIM
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK



!**** *SPAREORD*  - Reorder spectral data from old to new (grib-like) ordering
!                  or the reverse.

!     Purpose.
!     --------
!           Reorder spectral data structure from file ordering to model ordering 
!           or vice-versa.
!           SM : model ordering is the global spectrum ordering used in the model.

!**   Interface.
!     ----------
!        *CALL* *SPAREORD

!        Explicit arguments :
!        --------------------
!        PSPFILE          : Spectral array ready to be used for file
!        PSPBUF           : Spectral (buffer) array ready to be used in model
!        KFLEV            : Number of fields
!        LD_FILE_TO_MODEL : .TRUE. to order from file to model data structure

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!     None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Lars Isaksen *ECMWF*
!      Original : 95-06-10

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 04-04-14 : optimization
!      R. El Ouaraini : 03-Oct-2006 : use IDIM0GG instead of NDIM0G !
!      R. El Khatib : 30-Mar-2010 : Optimisation (reverse arrays ordering)
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM)        ,INTENT(IN)    :: YDDIM
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPFILE(YDDIM%NSEFRE,KFLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPBUF(YDDIM%NSPEC2G,KFLEV) 
LOGICAL           ,INTENT(IN)    :: LD_FILE_TO_MODEL 
INTEGER(KIND=JPIM)               :: IDIM0GG(0:YDDIM%NSMAX)
INTEGER(KIND=JPIM) :: II, IM, ISP, JLEV, JM, JN, IPOS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPAREORD',0,ZHOOK_HANDLE)
ASSOCIATE(NSEFRE=>YDDIM%NSEFRE, NSMAX=>YDDIM%NSMAX, NSPEC2G=>YDDIM%NSPEC2G)
!     ------------------------------------------------------------------

!*       1.   Rearrange spectral data structure 
!             ---------------------------------

IPOS=1
DO JN=0,NSMAX
   IDIM0GG(JN)=IPOS
   IPOS=IPOS+(NSMAX+1-JN)*2
ENDDO

IF (LD_FILE_TO_MODEL) THEN
  DO JLEV=1,KFLEV
    II=0
    DO JN=0,NSMAX
      DO JM=-JN,-1
        ISP=IDIM0GG(-JM)+(JN+JM)*2 +1
        II = II + 1
        PSPBUF(ISP,JLEV)=PSPFILE(II,JLEV)
      ENDDO
      ISP=IDIM0GG(0)+JN*2
      II = II + 1
      PSPBUF(ISP,JLEV)=PSPFILE(II,JLEV)
      PSPBUF(ISP+1,JLEV)=0.0_JPRB
      DO JM=1,JN
        ISP=IDIM0GG(JM)+(JN-JM)*2
        II = II + 1
        PSPBUF(ISP,JLEV)=PSPFILE(II,JLEV)
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    II=0
    DO JN=0,NSMAX
      DO JM=-JN,JN
        IM=ABS(JM)
        IF (JM < 0) THEN
          ISP=IDIM0GG(IM)+(JN-IM)*2 +1
        ELSE
          ISP=IDIM0GG(IM)+(JN-IM)*2
        ENDIF
        II = II + 1
        PSPFILE(II,JLEV)=PSPBUF(ISP,JLEV)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPAREORD',1,ZHOOK_HANDLE)
END SUBROUTINE SPAREORD
