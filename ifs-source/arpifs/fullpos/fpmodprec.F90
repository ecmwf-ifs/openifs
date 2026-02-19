! (C) Copyright 1989- Meteo-France.

!OPTION! -pvctl noouterunroll -pvctl noouterstrip
SUBROUTINE FPMODPREC(YDFPSTRUCT,YDGEOMETRY,KGP,PGPBUF,KFLDPTP)

!**** *FPMODPREC*  - MODIFY PRECIPITATION FIELDS FOR POST-PROCESSING

!     PURPOSE.
!     --------
!        To modify precipitation fields before horizontal post-processing.

!**   INTERFACE.
!     ----------
!       *CALL* *FPMODPREC*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!        KGP       - number of fields to modify
!        KFLDPTP   - array of field pointers for extraction

!        INPUT/OUTPUT:
!        PGPBUF - memory buffer

!        IMPLICIT ARGUMENTS
!        --------------------
!        See above 'use yom...'.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!      * Is called by FPMODCFU and FPMODXFU.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 96-08-06

!     MODIFICATIONS.
!     --------------
!      G.Mozdzynski : 02-10-01 support for radiation on-demand comms
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      R. El Khatib  20-May-2005 NFPWIDE moved to YOMWFPDS
!      R. El Khatib  17-Oct-2008 bf for ALADIN fullpos with equal distribution
!                             of C+I+E (RDISTR_E=1.)
!      K. Yessad Dec 2008: merge SLCOMM+SLCOMM1 -> SLCOMM.
!      K. Yessad Dec 2008: merge the different (E)SLEXTPOL.. -> (E)SLEXTPOL.
!      E. Sevault 18-May-2009 : SX9 compiler workaround : optimisations
!            "outerunroll" and "outerstrip" manually switched off to prevent
!            unexpected numerical variations. Note that this workaround may
!            be irrelevent when the source code of this subroutine is deeply 
!            re-written.
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (May 2012): further cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 27-Jul-2016 interface to eslextpol
!      R. El Khatib 11-Dec-2018 Bugfix : latitude rows can extend over LAM E-zone
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0       , ONLY : LELAM
USE YOMMP0       , ONLY : MY_REGION_EW, NPROC
USE EINT_MOD, ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),   INTENT(IN) :: YDFPSTRUCT
TYPE (GEOMETRY)    , INTENT (IN)    :: YDGEOMETRY
REAL (KIND=JPRB)   , INTENT (INOUT) :: PGPBUF (:,:,:)
INTEGER (KIND=JPIM), INTENT (IN)    :: KGP 
INTEGER (KIND=JPIM), INTENT (IN)    :: KFLDPTP(SIZE (PGPBUF, 2))

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INBTOPG(YDGEOMETRY%YRDIM%NDGLG), INBBOTG(YDGEOMETRY%YRDIM%NDGLG)
INTEGER(KIND=JPIM) :: ILAT(YDGEOMETRY%YRGEM%NGPTOT),ILON(YDGEOMETRY%YRGEM%NGPTOT),ILOTOP(YDGEOMETRY%YRGEM%NGPTOT),&
 & ILOBOT(YDGEOMETRY%YRGEM%NGPTOT)
INTEGER(KIND=JPIM) :: IL0(YDGEOMETRY%YRDIM%NPROMA,0:2)
REAL(KIND=JPRB) :: ZLON(YDGEOMETRY%YRGEM%NGPTOT), ZWTOP(YDGEOMETRY%YRGEM%NGPTOT), ZWBOT(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZPARFP1(KGP)
REAL(KIND=JPRB) :: ZDATAI(YDFPSTRUCT%NASLB1,KGP)

INTEGER(KIND=JPIM) :: IB, IEND, IFLDPTP, IGLG, IGLO, ILEN, IROF, IST, IOFF, IFLDCORE
INTEGER(KIND=JPIM) :: JF, JGL, ILA2G, ILO, IBLK
INTEGER(KIND=JPIM) :: JKGLO, JLON, JROF, ILA0, ILA0G, ILA1, ILA1G, ILA2, ILAMAX
INTEGER(KIND=JPIM) :: ILO1, ILO2, IROF0, IDUMARR(2)

REAL(KIND=JPRB) :: ZRATIO
REAL(KIND=JPRB) :: Z0, ZVAL0, ZVALE, ZVALN, ZVALNE, ZVALNW, ZVALS, ZVALSE
REAL(KIND=JPRB) :: ZVALW, ZVALX, ZVALSW
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     INBTOPG: number of neighbour points on top for each latitude
!     INBBOTG: number of neighbour points on bottom for each latitude
!     ILOTOP : local longitude index of nearest top left point
!     ILOBOT : local longitude index of nearest bottom left point
!     IWTOP  : weight of nearest top left point
!     IWBOT  : weight of nearest bottom left point
!     ZLON   : (work array) longitude abcissa (real) 
!               of the current point on top (bottom) row 
!     ILAT   : current local latitude index
!     ILON   : current local longitude index

!     ------------------------------------------------------------------

#include "fphalo.intfb.h"
#include "slcomm.intfb.h"
#include "slextpol.intfb.h"
#include "eslextpol.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPMODPREC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGENL=>YDDIM%NDGENL, NDGLG=>YDDIM%NDGLG, NDGUXG=>YDDIM%NDGUXG, &
 & NDLON=>YDDIM%NDLON, NPROMA=>YDDIM%NPROMA, &
 & NGPTOT=>YDGEM%NGPTOT, NLOENG=>YDGEM%NLOENG, &
 & NFRSTLOFF=>YDMP%NFRSTLOFF, NONL=>YDMP%NONL, NPTRFLOFF=>YDMP%NPTRFLOFF, &
 & NSTA=>YDMP%NSTA)
!     ------------------------------------------------------------------

!*       1. FIND NEIGHBOUR POINTS
!           ---------------------

INBTOPG(1)=1
DO JGL=2,NDGLG
  IF (NLOENG(JGL-1) == NLOENG(JGL)) THEN
    INBTOPG(JGL)=1
  ELSE
    INBTOPG(JGL)=2
  ENDIF
ENDDO
DO JGL=1,NDGLG-1
  IF (NLOENG(JGL+1) == NLOENG(JGL)) THEN
    INBBOTG(JGL)=1
  ELSE
    INBBOTG(JGL)=2
  ENDIF
ENDDO
INBBOTG(NDGLG)=1

IB=MY_REGION_EW
IROF0=0
! Loop on all the local latitudes + the possible split ones
DO JGL=1,NDGENL
  IGLG=JGL+NFRSTLOFF ! global adressing for latitudes
  IGLO=JGL+NPTRFLOFF ! global adressing for NSTA, NONL including lat. splitting
  ILEN=NONL(IGLO,IB)
  IROF=IROF0
  DO JLON=NSTA(IGLO,IB),NSTA(IGLO,IB)+ILEN-1
    IROF=IROF+1
    ILAT(IROF)=JGL
    ILON(IROF)=JLON
  ENDDO
  IROF=IROF0
  IF (INBTOPG(IGLG) == 1) THEN
    DO JLON=NSTA(IGLO,IB),NSTA(IGLO,IB)+ILEN-1
      IROF=IROF+1
      ILOTOP(IROF)=JLON
      ZWTOP (IROF)=1.0_JPRB
    ENDDO
  ELSE
    ZRATIO=REAL(NLOENG(MAX(IGLG-1,1)),JPRB)/REAL(NLOENG(IGLG),JPRB)
    DO JLON=NSTA(IGLO,IB),NSTA(IGLO,IB)+ILEN-1
      IROF=IROF+1
      ZLON  (IROF)=REAL(JLON-1,JPRB)*ZRATIO + 1.0_JPRB
      ILOTOP(IROF)=INT(ZLON(IROF))
      ZWTOP (IROF)=1.0_JPRB-(ZLON(IROF)-REAL(ILOTOP(IROF),JPRB))
    ENDDO
  ENDIF
  IROF=IROF0
  IF (INBBOTG(IGLG) == 1) THEN
    DO JLON=NSTA(IGLO,IB),NSTA(IGLO,IB)+ILEN-1
      IROF=IROF+1
      ILOBOT(IROF)=JLON
      ZWBOT (IROF)=1.0_JPRB
    ENDDO
  ELSE
    ZRATIO=REAL(NLOENG(MIN(IGLG+1,NDGLG)),JPRB)/REAL(NLOENG(IGLG),JPRB)
    DO JLON=NSTA(IGLO,IB),NSTA(IGLO,IB)+ILEN-1
      IROF=IROF+1
      ZLON  (IROF)=REAL(JLON-1,JPRB)*ZRATIO + 1.0_JPRB
      ILOBOT(IROF)=INT(ZLON(IROF))
      ZWBOT (IROF)=1.0_JPRB-(ZLON(IROF)-REAL(ILOBOT(IROF),JPRB))
    ENDDO
  ENDIF
  IROF0=IROF0+ILEN
ENDDO

!     ------------------------------------------------------------------

!*       2. FIELD MODIFICATION.
!           -------------------

ZPARFP1(:)=1.0_JPRB
IST=1
IFLDCORE=SIZE(PGPBUF,DIM=2)

!          *** transfer data in a halo ZDATAI.

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,IEND,IBLK,JROF,JF,IFLDPTP,IOFF)
DO JKGLO=1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBLK=1+(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA
  CALL FPHALO(YDFPSTRUCT,NPROMA,KGP,JKGLO,IEND,IFLDCORE,PGPBUF(:,:,IBLK),ZDATAI,KFLDPTP)
ENDDO
!$OMP END PARALLEL DO

!          *** Complete buffer for extra-longitudes and
!              compute halo (call to SLCOMM)
!              and extra-polar latitudes (call to SLEXTPOL).

IF (YDFPSTRUCT%NSLWIDE > 0) THEN
  IF (NPROC > 1) THEN
    CALL SLCOMM(YDFPSTRUCT,IDUMARR,KGP,.FALSE.,0,ZDATAI)
  ENDIF
  IF (LELAM) THEN
    CALL ESLEXTPOL(YDGEOMETRY,YDFPSTRUCT,KGP,IDUMARR,1,ZDATAI) 
  ELSE
    CALL SLEXTPOL(YDGEOMETRY%YRDIM,YDFPSTRUCT,KGP,IDUMARR,1,ZPARFP1,ZDATAI)  
  ENDIF
ENDIF

!          *** Store modified field

DO JKGLO=1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBLK=1+(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA

  DO JROF=IST,IEND
!         * Latitude of current row:
    ILA1=ILAT(JROF+JKGLO-1)
    ILA1G=ILAT(JROF+JKGLO-1)+NFRSTLOFF
!         * Latitude of southern row:
    ILA2G =ILA1G+1
    IF (LELAM .AND. ILA2G > NDGUXG) THEN
      ILAMAX=NDGLG ! task fully in E-zone
    ELSE
      ILAMAX=NDGUXG
    ENDIF
    ILA2=MIN(MAX(ILA2G,1),ILAMAX)-NFRSTLOFF
!         * Latitude of northern row:
    ILA0G =ILA1G-1
    IF (LELAM .AND. ILA0G > NDGUXG) THEN
      ILAMAX=NDGLG ! task fully in E-zone
    ELSE
      ILAMAX=NDGUXG
    ENDIF
    ILA0=MIN(MAX(ILA0G,1),ILAMAX)-NFRSTLOFF
!         * longitudes :
    ILO =MIN(MAX(ILOTOP(JROF+JKGLO-1),1),NDLON)-1
    ILO1=MIN(MAX(ILON(JROF+JKGLO-1),1),NDLON)-1
    ILO2=MIN(MAX(ILOBOT(JROF+JKGLO-1),1),NDLON)-1
    IL0(JROF,0)=YDFPSTRUCT%NSLOFF(ILA0)+YDFPSTRUCT%NSLEXT(ILO ,ILA0)
    IL0(JROF,1)=YDFPSTRUCT%NSLOFF(ILA1)+YDFPSTRUCT%NSLEXT(ILO1,ILA1)
    IL0(JROF,2)=YDFPSTRUCT%NSLOFF(ILA2)+YDFPSTRUCT%NSLEXT(ILO2,ILA2)
  ENDDO

  DO JF=1,KGP
    IFLDPTP=KFLDPTP(JF)
    DO JROF=IST,IEND

      ZVAL0 =ZDATAI(IL0(JROF,1)+1,JF)
      ZVALE =ZDATAI(IL0(JROF,1)+2,JF)
      ZVALSW=ZDATAI(IL0(JROF,2)+1,JF)
      ZVALSE=ZDATAI(IL0(JROF,2)+2,JF)
      ZVALNW=ZDATAI(IL0(JROF,0)+1,JF)
      ZVALNE=ZDATAI(IL0(JROF,0)+2,JF)
      ZVALW =ZDATAI(IL0(JROF,1)  ,JF)

      ZVALS=(1.0_JPRB-ZWBOT(JROF+JKGLO-1))*ZVALSE + ZWBOT(JROF+JKGLO-1)*ZVALSW
      ZVALN=(1.0_JPRB-ZWTOP(JROF+JKGLO-1))*ZVALNE + ZWTOP(JROF+JKGLO-1)*ZVALNW

      ZVALX=MAX(ZVALN,ZVALS,ZVALW,ZVALE)
      Z0=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-ZVAL0))
      PGPBUF(JROF,IFLDPTP,IBLK) = Z0*ZVAL0 -(1.0_JPRB-Z0)*ZVALX

    ENDDO
  ENDDO

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPMODPREC',1,ZHOOK_HANDLE)
END SUBROUTINE FPMODPREC

