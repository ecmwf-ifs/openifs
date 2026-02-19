! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUPROCFP_DEP(YDGEOMETRY,KGP,PLA,PLO,KPROCFP)

!**** *SUPROCFP_DEP*  - Setup processor number - Fullpos

!     PURPOSE.
!     --------
!        To initialize the array telling for each output point which 
!        processor treats this point, in DM-distribution linked to
!        the departure geometry.

!**   INTERFACE.
!     ----------
!       *CALL* *SUPROCFP_DEP*

!        EXPLICIT ARGUMENTS
!        --------------------
!         INPUT : 
!          KGP      : number of points on the interpolation grid
!          PLA, PLO : global latitudes, longitudes of all output points
!         OUTPUT : 
!          KPROCFP  : processor number each point belongs

!        IMPLICIT ARGUMENTS
!        --------------------
!        See module below

!     METHOD.
!     -------
!         ** Associate to interpolation point a model grid point, the rule is
!            the following:
!            - if the Gaussian/Lobatto latitude of the north-western model
!              grid-point is between 1 and YRDIM%NDGLG, THEN take it.
!            - if the Gaussian/Lobatto latitude of the north-western model
!              grid-point is zero (Northern pole or high resolution pole),
!              add 1 to latitude (SW point instead of NW point).
!            Remarks:
!            - LELAM=.T. coded according to similar relations used in
!              the routines SUEHOW1 and ELASCAW (compute IFPLANW).
!              the differences to Arpege come from the simplified geometry.
!            - a more complicated rule has to be studied later to minimize the
!              inter-processor variance of number of interpolation points
!              => the trick would be to define a larger halo so each processor
!                 could treat the same number of points (in that case, several
!                 points would be interpolated without the core!
!                 To be tried ... but the distribution of the output point will
!                 be difficult and the communication may increase ?? (Ryad. El Khatib)
!                 Then define the processor treating the interpolation point 
!                 (KPROCFP) as the processor treating the associated model 
!                 grid point (NGLOBALPROC).
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!       None

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE* (SUPROCFP)
!        K. Yessad (renaming into SUPROCFP_DEP and adapting/cleaning)

!     MODIFICATIONS.
!     --------------
!        Original : 98-10-05 from sufpg
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K.Yessad : 27-Feb-2007: rename into SUPROCFP_DEP, adapt, clean.
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!        R. El Khatib 27-Sep-2013 Differentiation between interpolation grid and
!        output grid
!        G.Mozdzynski (Nov 2013): delete call to suprocgp, as we have already computed
!                                 what it is doing (and far more efficiently)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCT0   , ONLY : LELAM
USE YOMCST   , ONLY : RPI

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KGP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLA(KGP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLO(KGP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPROCFP(KGP)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IASGPLAG, IASGPLOG, IAXNW, IFPLANW, ILAG,&
 & IROF, IZLATG, JGLGLO, JJPG  
INTEGER(KIND=JPIM) :: ISTAGPG(YDGEOMETRY%YRDIM%NDGLG)

REAL(KIND=JPRB) :: ZFRL, ZLOO, ZLSDEPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPROCFP_DEP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDEGEO=>YDGEOMETRY%YREGEO)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDGUNG=>YDDIM%NDGUNG, NDLUNG=>YDDIM%NDLUNG, &
 & NDLUXG=>YDDIM%NDLUXG, &
 & EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, &
 & NGPTOTG=>YDGEM%NGPTOTG, NLOENG=>YDGEM%NLOENG, R4JP=>YDGEM%R4JP, &
 & NGLOBALPROC=>YDMP%NGLOBALPROC)
!-----------------------------------------------------------------------

!*     2. COMPUTE KPROCFP
!         ---------------

!     1.1 Starting adress of each global latitude

IROF=1
DO JGLGLO=1,NDGLG
  ISTAGPG(JGLGLO)=IROF
  IROF=IROF+NLOENG(JGLGLO)
ENDDO

IF (LELAM) THEN
  DO JJPG=1,KGP

    ZFRL =PLA(JJPG)/EDELY
    ILAG =INT(ZFRL)+NDGUNG
    ! * IFPLANW: DM-global regular y-coordinate number situated 
    !   immediately south of the interpolation point.
    IFPLANW=ILAG
    ! * IASGPLAG: DM-global regular y-coordinate number of the
    !   model grid-point which is associated to the interpolation point.
    IASGPLAG=MAX(1,MIN(NDGLG,IFPLANW))
    ! * IASGPLOG: x-coordinate number of the point of y-coord. number 
    !   IASGPLAG, immediately west to the interpolation point.
    ZFRL  =PLO(JJPG)/EDELX
    IAXNW =INT(ZFRL)+NDLUNG
    IASGPLOG=MAX(NDLUNG,MIN(IAXNW,NDLUXG))

    KPROCFP(JJPG)=NGLOBALPROC(ISTAGPG(IASGPLAG)+IASGPLOG-1)

  ENDDO
ELSE
  DO JJPG=1,KGP

    IZLATG=INT(R4JP*(0.5_JPRB*RPI-PLA(JJPG))+0.75_JPRB)
    ILAG  =INT(R4JP*(0.5_JPRB*RPI-PLA(JJPG))+0.75_JPRB)&
     & +NINT(SIGN(0.5_JPRB,YDCSGLEG%RLATIG(IZLATG)-PLA(JJPG))-1.5_JPRB)  
    ! * IFPLANW: DM-global Gaussian/Lobatto latitude number situated 
    !   immediately north of the interpolation point.
    IFPLANW=ILAG+1
    ! * IASGPLAG: DM-global Gaussian/Lobatto latitude number of the
    !   model grid-point which is associated to the interpolation point.
    IASGPLAG=MAX(1,MIN(NDGLG,IFPLANW))
    ! * IASGPLOG: Longitude number of the point of latitude number 
    !   IASGPLAG, immediately west to the interpolation point.
    ZLOO = MOD(PLO(JJPG),2.0_JPRB*RPI)
    ZLSDEPI=REAL(NLOENG(IASGPLAG),JPRB)/(2.0_JPRB*RPI)
    IASGPLOG=INT(ZLOO*ZLSDEPI)
    IASGPLOG=IASGPLOG+NLOENG(IASGPLAG)*MAX(0,1-ABS(IASGPLOG))

    KPROCFP(JJPG)=NGLOBALPROC(ISTAGPG(IASGPLAG)+IASGPLOG-1)

  ENDDO
ENDIF

!-----------------------------------------------------------------------


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPROCFP_DEP',1,ZHOOK_HANDLE)
END SUBROUTINE SUPROCFP_DEP
