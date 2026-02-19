! (C) Copyright 1989- Meteo-France.

SUBROUTINE HPOS_DYN(YDQTYPE,YDNAMFPSCI,YDAFN,LDFPOSHOR,YDGEOMETRY,KEND,YDGSGEOM,KGPP,KAUX,PGPP,PAUX,KFLDIN,PARFP1,PFPBUF1)

!**** *HPOS_DYN*  - HORIZONTAL POST-PROCESSING

!     PURPOSE.
!     --------
!        FILL THE "SEMI-LAGRANGIAN" BUFFERS WITH THE FIELDS TO TRANSFER
!        OR INTERPOLATE

!**   INTERFACE.
!     ----------
!       *CALL* *HPOS_DYN*

!        EXPLICIT ARGUMENTS
!        --------------------
!         INPUT:
!          LDFPOSHOR  : .true. if actual horizontal interpolations
!          KEND       : last adress in in horizontal dimension
!          YDGSGEOM   : Grid point geometry

!         OUTPUT:
!          PFPBUF1    : interpolation buffer containing the fields to interpolate

!         OUTPUT:
!          PARFP1     : parity for interpolation buffer PFPBUF1 1=scalar; -1=vector.

!        IMPLICIT ARGUMENTS
!        --------------------
!        se #include below.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*
!        ORIGINAL : 16-Sep-2016 from HPOS

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMAFN, ONLY : TAFN
USE YOMFPC, ONLY : TNAMFPSCI
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN
USE YOMGSGEOM, ONLY : TGSGEOM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE (TNAMFPSCI),  INTENT(IN)    :: YDNAMFPSCI
TYPE(TAFN)        ,INTENT(IN)    :: YDAFN
LOGICAL           ,INTENT(IN)    :: LDFPOSHOR
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
TYPE(TGSGEOM)     ,INTENT(IN)    :: YDGSGEOM
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPP
INTEGER(KIND=JPIM),INTENT(IN)    :: KAUX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGPP(YDGEOMETRY%YRDIM%NPROMA,KGPP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAUX(YDGEOMETRY%YRDIM%NPROMA,KAUX)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDIN
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PARFP1(KFLDIN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPBUF1(YDGEOMETRY%YRDIM%NPROMA,KFLDIN) 

REAL(KIND=JPRB) :: ZGM(YDGEOMETRY%YRDIM%NPROMA)

INTEGER(KIND=JPIM) :: IORDER(KFLDIN)
INTEGER(KIND=JPIM) :: INC, IFLD, IST, IJ, IOFF, II, J, JFLD, JI

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fptensor.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HPOS_DYN',0,ZHOOK_HANDLE)
ASSOCIATE(LFPLOSP=>YDNAMFPSCI%LFPLOSP, &
 & TFP_DYNDS=>YDAFN%TFP_DYNDS, TFP=>YDAFN%TFP, NPROMA=>YDGEOMETRY%YRDIM%NPROMA)
!     ------------------------------------------------------------------

!*       1. PRELIMINARY INITIALISATIONS.
!           ----------------------------

!     IST     : first adress in post-processing buffers to read

IST=1

ZGM(IST:KEND)=YDGSGEOM%GM(IST:KEND)

!     ------------------------------------------------------------------

!*       3. DYNAMICS
!           --------

IOFF=0
DO J=1, YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  INC=SUM(YDQTYPE%ISPD(1:YDQTYPE%ILEV(IJ),IJ))
  DO IFLD=1,INC
    JFLD=IFLD+IOFF
    IF (YDQTYPE%ISF(IJ) >= 1) THEN
    ! Spectrally fitted data
      II=YDQTYPE%IGT1(IJ,0)+IFLD-1
      DO JI=IST,KEND
        PFPBUF1(JI,JFLD)=PGPP(JI,II)
      ENDDO
    ELSE
    ! Non-spectrally fitted data
      II=YDQTYPE%IGPX(IJ)+IFLD-1
      DO JI=IST,KEND
        PFPBUF1(JI,JFLD)=PAUX(JI,II)
      ENDDO
    ENDIF
  ENDDO
  IOFF=IOFF+INC
ENDDO

IF (LDFPOSHOR) THEN
  ! Tensors aspects :
  IOFF=0
  DO J=1, YDQTYPE%NFPOSDYN
    IJ=YDQTYPE%NFPTRDYN(J)
    INC=SUM(YDQTYPE%ISPD(1:YDQTYPE%ILEV(IJ),IJ))
    DO IFLD=1,INC
      JFLD=IFLD+IOFF
      IORDER(JFLD)=ABS(TFP_DYNDS(IJ)%IORDR)
    ENDDO
    IOFF=IOFF+INC
  ENDDO
  CALL FPTENSOR(IST,KEND,NPROMA,KFLDIN,IORDER,ZGM,PFPBUF1,PARFP1)
ELSE
  ! For security
  PARFP1(:)=HUGE(1._JPRB)
ENDIF

! Particular cases :
IOFF=0
DO J=1, YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  INC=SUM(YDQTYPE%ISPD(1:YDQTYPE%ILEV(IJ),IJ))
  DO IFLD=1,INC
    JFLD=IFLD+IOFF
    IF (IJ==TFP%ABS%ICOD .AND. LDFPOSHOR) THEN
      ! Absolute vorticity => add Coriolis parameter
      DO JI=IST,KEND
        PFPBUF1(JI,JFLD)=PFPBUF1(JI,JFLD)+YDGSGEOM%RCORI(JI)
      ENDDO
    ELSEIF (IJ==TFP%ABS%ICOD .AND. .NOT.LDFPOSHOR) THEN
      ! (Reduced) Absolute vorticity => add (Reduced) Coriolis parameter
      DO JI=IST,KEND
        PFPBUF1(JI,JFLD)=PFPBUF1(JI,JFLD)+YDGSGEOM%RCORI(JI)/(ZGM(JI)**2)
      ENDDO
    ELSEIF (IJ==TFP%SP%ICOD .AND. LDFPOSHOR .AND.(.NOT.LFPLOSP) ) THEN
      ! Ps => Ln(Ps)
      DO JI=IST,KEND
        PFPBUF1(JI,JFLD)=LOG(PFPBUF1(JI,JFLD))
      ENDDO
    ENDIF

  ENDDO
  IOFF=IOFF+INC
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('HPOS_DYN',1,ZHOOK_HANDLE)
END SUBROUTINE HPOS_DYN
