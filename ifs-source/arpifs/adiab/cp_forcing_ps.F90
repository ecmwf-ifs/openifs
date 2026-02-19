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

SUBROUTINE CP_FORCING_PS( &
 ! - INPUT --------------------------------------------------------------
 & YDGEOMETRY,KST,KEND,&
 ! - OUTPUT -------------------------------------------------------------
 & PAPRS)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMCT3   , ONLY : NSTEP

!     ------------------------------------------------------------------
!**** *CP_FORCING_PS* - 

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *CP_FORCING_PS(...)*

!        Explicit arguments :
!        --------------------
!        * INPUT:
!        KST                : first element of work
!        KEND               : last element of work

!        * OUTPUT:
!        PAPRS              : half-level pressure.

!        Implicit arguments : None.
!        --------------------

!     Method.
!     -------
!

!     Externals.
!     ----------

!     Reference.
!     ----------
       
!     Author.
!     -------
!     Romain ROEHRIG

! Modifications
! -------------
! End Modifications
!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAPRS(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF, JC
REAL(KIND=JPRB) :: ZPS
CHARACTER :: CLSTEP*20
CHARACTER :: CLFILE*200

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CP_FORCING_PS',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------

WRITE(CLSTEP,FMT='(I5)') NSTEP
DO JC=1,5
  IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
ENDDO
WRITE(CLFILE,FMT='(3A)') 'files/Ps_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
read(82,*) ZPS
CLOSE(82)

DO JROF=KST,KEND
  PAPRS(JROF,NFLEVG)=ZPS
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CP_FORCING_PS',1,ZHOOK_HANDLE)
END SUBROUTINE CP_FORCING_PS
