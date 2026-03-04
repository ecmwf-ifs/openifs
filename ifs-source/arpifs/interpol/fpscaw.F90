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

!ocl list_copy(32,YDFPSTRUCT%NSLOFF)
!ocl list_copy(32,YDFPSTRUCT%NSLEXT)
SUBROUTINE FPSCAW(YDFPSTRUCT,KPROMA,KSTART,KEND,KFPROW,KLAT,KLONN,KLON,KLOS,KLOSS,KL0)

!**** *FPSCAW  -  FULL-POS:
!                 Numbering of the points (I is the interpolation point):
!                     13       5       6      14

!                      7       1       2       8
!                                 (I)
!                      9       3       4      10

!                     15      11      12      16

!                 Once known the coordinates 
!                 (DM-local ilon, DM-local ilat) of the points numbered 
!                 13, 7, 9, 15 (respectively (KLONN,KLAT-1), (KLON,KLAT),
!                 (KLOS,KLAT+1), (KLOSS,KLAT+2)), one computes the position
!                 irof=KL0(inumlat,ilen) in the interpolation buffer
!                 dimensioned with NPROMA. KL0(inumlat =1,2,3,4,.) correspond 
!                 respectively with (KLONN,KLAT-1), (KLON,KLAT),
!                 (KLOS,KLAT+1), (KLOSS,KLAT+2), 
!                 i.e. resp. points numbered 13, 7, 9
!                 and 15 (for consistency with routine LASCAW).

!                 This routine has been written after LASCAW.
!                 SUHOW1+SUHOW2+FPSCAW for FULL-POS plays the same part as
!                 LASCAW for semi-Lagrangian or observation interpolations.

!                 Computations are DM-local if distributed memory.

!**   Interface.
!     ----------
!        *CALL* *FPSCAW(...)

!        Explicit arguments :
!        --------------------
!        All the other following quantities involved in DM computations 
!         are DM-local.

!        INPUT:
!          KPROMA       - horizontal dimension
!          KSTART     - first element of arrays where computations are 
!                       performed.
!          KEND       - depth of work.
!          KOFF      - pp packet offset

!        OUTPUT:
!          KL0        - index of the four western points
!                       of the 16 points interpolation grid.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        No external.

!     Reference.
!     ----------
!        Documentation about FULL-POS.

!     Author.
!     -------
!      K. YESSAD (MARCH 1997).

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 28-Jul-2006 Porting to NEC
!      K. Yessad: 28-02-2007 DM-features optimisations in FULL-POS
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      P.Marguinaud (Oct 2014): Allow a wider halo
!      R. El Khatib 27-Jul-2016 interpolations over C+I+E
!     ----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE EINT_MOD, ONLY : SL_STRUCT

!     ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),   INTENT(IN)  :: YDFPSTRUCT
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)  :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPROW
INTEGER(KIND=JPIM),INTENT(IN)  :: KLAT(KPROMA)
INTEGER(KIND=JPIM),INTENT(IN)  :: KLONN(KPROMA)
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON(KPROMA)
INTEGER(KIND=JPIM),INTENT(IN)  :: KLOS(KPROMA)
INTEGER(KIND=JPIM),INTENT(IN)  :: KLOSS(KPROMA)
INTEGER(KIND=JPIM),INTENT(OUT) :: KL0(KPROMA,KFPROW)

!     ----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILA0, ILA1, ILA1G, ILA2, ILA2G, ILA3, ILA3G, ILAG, ILO
INTEGER(KIND=JPIM) :: ILO1, ILO2, ILO3, JLEN

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPSCAW',0,ZHOOK_HANDLE)

!     ----------------------------------------------------------------------

!*       2.    COMPUTE ADRESSES
!              ----------------

IF (YDFPSTRUCT%NSLWIDE  == 0) THEN

  !      2.1  Pointer for  buffer transfer

  IF (KFPROW < 1) CALL ABOR1('FPSCAW : KFPROW TOO SMALL')

  DO JLEN=KSTART,KEND

    KL0(JLEN,1)=JLEN

  ENDDO

ELSE

  !      2.2  Pointer for 12 points interpolations and bilinear interpolation.

  IF (KFPROW < 4) CALL ABOR1('FPSCAW : KFPROW TOO SMALL')

  DO JLEN=KSTART,KEND

    ! * Latitude of point 1:
    ILA1G=KLAT(JLEN)+YDFPSTRUCT%NFRSTLOFF
    ILA1=ILA1G-YDFPSTRUCT%NFRSTLOFF
    ! * Latitude of point 3:
    ILA2G =ILA1G+1
    ILA2=ILA2G-YDFPSTRUCT%NFRSTLOFF
    ! * Latitude of point 5:
    ILAG  =ILA1G-1
    ILA0=ILAG-YDFPSTRUCT%NFRSTLOFF
    ! * Latitude of point 11:
    ILA3G =ILA2G+1
    ILA3=ILA3G-YDFPSTRUCT%NFRSTLOFF

    ILO =KLONN(JLEN)
    ILO1=KLON (JLEN)
    ILO2=KLOS (JLEN)
    ILO3=KLOSS(JLEN)

    KL0(JLEN,1)=YDFPSTRUCT%NSLOFF(ILA0)+YDFPSTRUCT%NSLEXT(ILO ,ILA0)
    KL0(JLEN,2)=YDFPSTRUCT%NSLOFF(ILA1)+YDFPSTRUCT%NSLEXT(ILO1,ILA1)
    KL0(JLEN,3)=YDFPSTRUCT%NSLOFF(ILA2)+YDFPSTRUCT%NSLEXT(ILO2,ILA2)
    KL0(JLEN,4)=YDFPSTRUCT%NSLOFF(ILA3)+YDFPSTRUCT%NSLEXT(ILO3,ILA3)

  ENDDO

ENDIF

!     ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPSCAW',1,ZHOOK_HANDLE)

END SUBROUTINE FPSCAW
