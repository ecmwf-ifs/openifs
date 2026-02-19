! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GEMS_TEND(YGFL,KIDIA,KFDIA,KLEV,KLON,PTENGFL,PTENC)

!**   INTERFACE.
!     ----------
!          *GEMS_TEND* IS CALLED FROM *CALLPAR*.

! INPUTS:
!  -------

! INPUTS/OUTPUTS:
!---------------

!-----------------------------------------------------------------------

!     Externals.
!     ---------


!     Author
!    --------
!         2008-10-15, R. Engelen 

!     Modifications :
!    ----------------
!        A. Inness     28-Mar-2012  Changes for CHEM fields
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD

IMPLICIT NONE

TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENGFL(KLON,KLEV,YGFL%NDIM1)

! In the declaration below the INTENT attribute has been removed to comply 
! strict f95 standards. The attribute INTENT(INOUT) should be put back after all 
! compiler (especially NEC) support this extension. - R. El Khatib 04-Jun-2009
REAL(KIND=JPRB), POINTER :: PTENC(:,:,:)

INTEGER(KIND=JPIM) :: ITRC
INTEGER(KIND=JPIM) :: JK, JL, JEXT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GEMS_TEND',0,ZHOOK_HANDLE)
ASSOCIATE(NAERO=>YGFL%NAERO, NCHEM=>YGFL%NCHEM, NDIM1=>YGFL%NDIM1, &
 & NGHG=>YGFL%NGHG, YAERO=>YGFL%YAERO, YCHEM=>YGFL%YCHEM, YGHG=>YGFL%YGHG)

ITRC=0
DO JEXT=1,NGHG
  ITRC=ITRC+1 
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PTENGFL(JL,JK,YGHG(JEXT)%MP1)=PTENC(JL,JK,ITRC)
    ENDDO
  ENDDO
ENDDO
DO JEXT=1,NAERO
  ITRC=ITRC+1 
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PTENGFL(JL,JK,YAERO(JEXT)%MP1)=PTENC(JL,JK,ITRC)
    ENDDO
  ENDDO
ENDDO
DO JEXT=1,NCHEM
  ITRC=ITRC+1 
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PTENGFL(JL,JK,YCHEM(JEXT)%MP1)=PTENC(JL,JK,ITRC)
    ENDDO
  ENDDO
ENDDO



!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GEMS_TEND',1,ZHOOK_HANDLE)
END SUBROUTINE GEMS_TEND
