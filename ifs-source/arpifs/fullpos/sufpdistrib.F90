! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPDISTRIB(YDNAMFPSCI,YDFPGEO_DEP,YDFPGIND,YDFPGEO,YDFPUSERGEO,KNUMPROCFP_DEP,KNUMPROCFP)

!**** *SUFPDISTRIB*  - SET UP MESSAGE PASSINF FOR FULL-POS

!     PURPOSE.
!     --------
!        To initialize control array for the distribution and 
!        to distribute some global fields (DM distribution linked
!        to arrival geometry).

!        Computes the following variables:
!        NFPRGPL,NFPRGPLX
!        NFPRGPNUM,NFPRGPIND
!        NUMPROCFP

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPDISTRIB*

!        EXPLICIT ARGUMENTS
!        ------------------
!         INPUT : 
!         KNUMPROCFP_DEP : MPI task number for each point on the interpolation grid
!         Output :
!         KNUMPROCFP     : MPI task number for each point on the output grid

!        IMPLICIT ARGUMENTS
!        ------------------
!        See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!       SUPROCFP

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      Original : 98-10-05 from sufpg

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad: 27-Feb-2007 Move old code in SUFPDISTRIB_DEP, adapt to DM-arrival geometry.
!      R. El Khatib  24-Jul-2012 LFPDISTRIB replaced by IDISTRIB
!      R. El Khatib & Tayfun Dalkilic 13-Sep-2012 IDISTRIB=2
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2 (differentiation
!      between interpolation grid and output grid)
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 27-Jul-2016 bugfix for the case NFPBW* > 0 + recode LWIDER_DOM + optimize for Boyd
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMFPGEO , ONLY : TFPGEO
USE YOMFPGIND, ONLY : TFPGIND
USE YOMMP0   , ONLY : MYPROC, NPROC
USE YOMFPC   , ONLY : TNAMFPSCI, LALLOFP

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
TYPE (TFPGEO),  INTENT(IN) :: YDFPGEO_DEP
TYPE (TFPGIND),  INTENT(INOUT) :: YDFPGIND
TYPE (TFPGEO),  INTENT(IN) :: YDFPGEO
TYPE (TFPUSERGEO) ,INTENT(IN)    :: YDFPUSERGEO(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUMPROCFP_DEP(YDFPGEO_DEP%NFPRGPG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUMPROCFP(YDFPGEO%NFPRGPG)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JI, IGP, IBWXO, IBWYO, IPROC_DEP, IPROC
INTEGER(KIND=JPIM) :: JGL, JLON, ILONF, ILATF, IND, IND1
!     Arrays for data transposition :
INTEGER(KIND=JPIM) :: IDEP(YDFPGEO%NFPRGPG), IARR(YDFPGEO_DEP%NFPRGPG)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPDISTRIB',0,ZHOOK_HANDLE)
ASSOCIATE(NFPBOYD=>YDNAMFPSCI%NFPBOYD, &
 & NFPRGPL_DEP=>YDFPGEO_DEP%NFPRGPL, NFPRGPG_DEP=>YDFPGEO_DEP%NFPRGPG, &
 & NFPRGPIND_DEP=>YDFPGEO_DEP%NFPRGPIND, &
 & NFPRGPL=>YDFPGEO%NFPRGPL, NFPRGPG=>YDFPGEO%NFPRGPG, NFPRGPIND=>YDFPGEO%NFPRGPIND, &
 & NLON=>YDFPUSERGEO%NLON, NLAT=>YDFPUSERGEO%NLAT, NFPLUX=>YDFPUSERGEO%NFPLUX, &
 & NFPRESOL=>YDFPUSERGEO%NFPRESOL, NFPSIZEG=>YDFPUSERGEO%NFPSIZEG, &
 & NFPGUX=>YDFPUSERGEO%NFPGUX, NFPBWX=>YDFPUSERGEO%NFPBWX, NFPBWY=>YDFPUSERGEO%NFPBWY, &
 & CFPGRID=>YDFPUSERGEO%CFPGRID)

!-----------------------------------------------------------------------


!*     2. TRANSPOSITION
!         -------------

  ALLOCATE(YDFPGIND%NFP2SEND(NPROC))
  IF (LALLOFP) THEN
    WRITE(NULOUT,9990) 'NFP2SEND',SIZE(YDFPGIND%NFP2SEND),SHAPE(YDFPGIND%NFP2SEND)
  ENDIF
  ALLOCATE(YDFPGIND%NFP2RECV(NPROC))
  IF (LALLOFP) THEN
    WRITE(NULOUT,9990) 'NFP2RECV',SIZE(YDFPGIND%NFP2RECV),SHAPE(YDFPGIND%NFP2RECV)
  ENDIF

  IF (NFPBOYD /= 0) THEN
!   Global addresses of the interpolation gridpoints which are not in the output grid are set to zero :
    IARR(:)=0
    IDEP(:)=0
    IBWXO=NFPBWX(1)+(NLON(1)-NFPLUX(1))
    IBWYO=NFPBWY(1)+(NLAT(1)-NFPGUX(1))
    DO JLON=1+IBWXO,NLON(1)+IBWXO
      ILONF=MODULO(JLON-1-NFPBWX(1)-(NLON(1)-NFPLUX(1)),NLON(1))+1
      DO JGL=1+IBWYO,NLAT(1)+IBWYO
        ILATF=MODULO(JGL-1-NFPBWY(1)-(NLAT(1)-NFPGUX(1)),NLAT(1))+1
        IND1=(NFPLUX(1)+2*NFPBWX(1)+2*(NLON(1)-NFPLUX(1)))*(JGL-1)+JLON
        IND=NLON(1)*(ILATF-1)+ILONF
!       Global addresses of the output gridpoints in the interpolation grid
        IDEP(IND)=IND1
!       Global addresses of the interpolation gridpoints in the output grid 
        IARR(IND1)=IND
      ENDDO
    ENDDO
    ALLOCATE(YDFPGIND%NFP2SENDG(NFPRGPG_DEP))
    IF (LALLOFP) THEN
      WRITE(NULOUT,9990) 'NFP2SENDG',SIZE(YDFPGIND%NFP2SENDG),SHAPE(YDFPGIND%NFP2SENDG)
    ENDIF
  ELSE
    DO JI=1,NFPRGPG
      IDEP(JI)=JI
      IARR(JI)=JI
    ENDDO
  ENDIF

  YDFPGIND%NFP2SEND(:)=0
  DO JI=1,NFPRGPL_DEP
    IGP=IARR(NFPRGPIND_DEP(JI,MYPROC))
!   Skip gridpoints of the "Boyd window"
    IF (IGP > 0) THEN
      IPROC=KNUMPROCFP(IGP)
      YDFPGIND%NFP2SEND(IPROC)=YDFPGIND%NFP2SEND(IPROC)+1
    ENDIF
  ENDDO
! local indexes of the gridpoints in the source task to be distributed to the target distribution
  ALLOCATE(YDFPGIND%NFPSOURCE(MAXVAL(YDFPGIND%NFP2SEND(:)),NPROC))
  IF (LALLOFP) THEN
    WRITE(NULOUT,9990) 'NFPSOURCE',SIZE(YDFPGIND%NFPSOURCE),SHAPE(YDFPGIND%NFPSOURCE)
  ENDIF
  YDFPGIND%NFP2SEND(:)=0
  DO JI=1,NFPRGPL_DEP
    IGP=IARR(NFPRGPIND_DEP(JI,MYPROC))
!   Skip gridpoints of the "Boyd window"
    IF (IGP > 0) THEN
      IPROC=KNUMPROCFP(IGP)
      YDFPGIND%NFP2SEND(IPROC)=YDFPGIND%NFP2SEND(IPROC)+1
      YDFPGIND%NFPSOURCE(YDFPGIND%NFP2SEND(IPROC),IPROC)=JI
    ENDIF
  ENDDO

  YDFPGIND%NFP2RECV(:)=0
  DO JI=1,NFPRGPL
    IGP=IDEP(NFPRGPIND(JI,MYPROC))
    IF (IGP > 0) THEN
      IPROC_DEP=KNUMPROCFP_DEP(IGP)
      YDFPGIND%NFP2RECV(IPROC_DEP)=YDFPGIND%NFP2RECV(IPROC_DEP)+1
    ENDIF
  ENDDO
! Local indexes of the gridpoints in the target task to be received from the source distribution
  ALLOCATE(YDFPGIND%NFPTARGET(MAXVAL(YDFPGIND%NFP2RECV(:)),NPROC))
  IF (LALLOFP) THEN
    WRITE(NULOUT,9990) 'NFPTARGET',SIZE(YDFPGIND%NFPTARGET),SHAPE(YDFPGIND%NFPTARGET)
  ENDIF
  YDFPGIND%NFP2RECV(:)=0
  DO JI=1,NFPRGPL
    IGP=IDEP(NFPRGPIND(JI,MYPROC))
    IF (IGP > 0) THEN
      IPROC_DEP=KNUMPROCFP_DEP(IGP)
      YDFPGIND%NFP2RECV(IPROC_DEP)=YDFPGIND%NFP2RECV(IPROC_DEP)+1
      YDFPGIND%NFPTARGET(YDFPGIND%NFP2RECV(IPROC_DEP),IPROC_DEP)=JI
    ENDIF
  ENDDO

  IF (NFPBOYD /= 0) THEN

    YDFPGIND%NFP2SENDG(:)=0
    DO JI=1,NFPRGPG_DEP
      IGP=IARR(JI)
!     Skip gridpoints of the "Boyd window"
      IF (IGP > 0) THEN
        IPROC=KNUMPROCFP(IGP)
        YDFPGIND%NFP2SENDG(IPROC)=YDFPGIND%NFP2SENDG(IPROC)+1
      ENDIF
    ENDDO
!   local indexes of the gridpoints in the source task to be distributed to the target distribution
    ALLOCATE(YDFPGIND%NFPSOURCEG(MAXVAL(YDFPGIND%NFP2SENDG(:)),NPROC))
    IF (LALLOFP) THEN
      WRITE(NULOUT,9990) 'NFPSOURCEG',SIZE(YDFPGIND%NFPSOURCEG),SHAPE(YDFPGIND%NFPSOURCEG)
    ENDIF
    YDFPGIND%NFP2SENDG(:)=0
    DO JI=1,NFPRGPG_DEP
      IGP=IARR(JI)
!     Skip gridpoints of the "Boyd window"
      IF (IGP > 0) THEN
        IPROC=KNUMPROCFP(IGP)
        YDFPGIND%NFP2SENDG(IPROC)=YDFPGIND%NFP2SENDG(IPROC)+1
        YDFPGIND%NFPSOURCEG(YDFPGIND%NFP2SENDG(IPROC),IPROC)=JI
      ENDIF
    ENDDO
  ENDIF

! ---------------------------------------------------------------------
9990 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPDISTRIB',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPDISTRIB
