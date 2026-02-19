! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUBFPOS(YDFPOS,YDHOSTGEOMETRY,YDMODEL,KFPOS,YDPOSGEOMETRY,YDXPOSGEOMETRY,YDPOSVAB,CDNAM)

!**** *SUBFPOS*  - CONSTRUCTOR OF FULLPOS

!     PURPOSE.
!     --------
!        To interface and call the constructor of 'Fullpos' object (geometries, filters, interpolators ...
!        ... but not the fields list request which may vary from one invokation of the post-processor to another)
!        There are 3 possible methods of construction for the geometries 
!        1/ with a set of namelists to define the geometries
!        2/ with a model geometry object to define the geometry
!        3/ with a set of model geometry objects to define the geometries.
!         That last construction method is here more because it is possible than because it is usefull ;-)

!**   INTERFACE.
!     ----------
!        1/ CALL SUBFPOS(YDFPOS,YDHOSTGEOMETRY,YDMODEL,KFPOS,YDPOSVAB,CDNAM)
!        2/ CALL SUBFPOS(YDFPOS,YDHOSTGEOMETRY,YDMODEL,KFPOS,YDPOSGEOMETRY,YDPOSVAB,CDNAM)
!        3/ CALL SUBFPOS(YDFPOS,YDHOSTGEOMETRY,YDMODEL,KFPOS,YDXPOSGEOMETRY,YDPOSVAB,CDNAM)

!        EXPLICIT ARGUMENTS
!        --------------------
!           YDHOSTGEOMETRY : host model geometry
!           KFPOS : configuration of the post-processing :
!                   1 : gridpoint post-processing, possibly with spectral filters
!                   2 : gridpoint/spectral post-processing (spectral outputs possible)
!           YDPOSGEOMETRY : requested post-processing geometry
!           YDXPOSGEOMETRY : requested post-processing multiple geometries
!           NB : YDPOSGEOMETRY and YDXPOSGEOMETRY are optional and exclusive
!           YDPOSVAB : vertical eta coordinates object for the post-processing
!           CDNAM : alternative namelist filename

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib     20-May-2005 NFPWIDE moved to YOMFPDS
!      K. Yessad: 27-02-2007 Optimisation of distributed memory in FULL-POS
!      K. Yessad Aug-2007: allow to relax the C1 constraint (NDLNPR=0) in NH
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      R. El Khatib : 01-Mar-2012 LFPOS => NFPOS
!      R. El Khatib  24-Jul-2012 move out geometry setup + modularization of distribution
!      R. El Khatib  31-Jul-2012 Setup spectral transforms
!      R. El Khatib 27-Jul-2016 dynamic computation of halo width for interpolations over C+I+E
!     ------------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM    ,JPRB
USE YOMHOOK      , ONLY : LHOOK   ,DR_HOOK, JPHOOK

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMLUN       , ONLY : NULNAM, NULOUT
USE YOMCT0       , ONLY : CFPNCF, CFDIRLST, CNPPATH
USE YOMFPD       , ONLY : TNAMFPD
USE YOMFPG       , ONLY : TNAMFPG, TNAMFPV
USE YOMFPF       , ONLY : TNAMFPF
USE YOMFPC       , ONLY : TNAMFPOBJ, LTRACEFP
USE YOMFPIOS     , ONLY : TNAMFPIOS
USE FULLPOS      , ONLY : TFPOS
USE YOMVERT      , ONLY : TVAB
USE YOMFPGEOMETRY, ONLY : LFPOSBUF
USE YOMMP0       , ONLY : MYPROC
                        
IMPLICIT NONE           
                        
TYPE(TFPOS)       ,INTENT(OUT) :: YDFPOS
TYPE(GEOMETRY)    ,INTENT(IN) :: YDHOSTGEOMETRY
TYPE(MODEL)       ,INTENT(IN) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN) :: KFPOS
TYPE(GEOMETRY)    ,INTENT(IN), OPTIONAL :: YDPOSGEOMETRY
TYPE(GEOMETRY)    ,INTENT(IN), OPTIONAL :: YDXPOSGEOMETRY(:)
TYPE(TVAB)        ,INTENT(IN), OPTIONAL :: YDPOSVAB
CHARACTER(LEN=*)  ,INTENT(IN), OPTIONAL :: CDNAM

!IFPCONF : configuration of the post-processing :
!  0 : vertical interpolation only (<CFPFMT='MODEL'>)
!  1 : gridpoint post-processing, possibly with spectral filters (<NFPOS=1>)
!  2 : gridpoint/spectral post-processing (spectral outputs possible) (<NFPOS=2>)
INTEGER(KIND=JPIM) :: IFPCONF, J, IFPDOM, IFPROMA, IFPROMA_DEP, IFPWIDE_NAM, IOS, IERR, ICMDSTAT
CHARACTER(LEN=180), ALLOCATABLE :: CLFPFN(:), CLFPCLIFNAME(:), CLFPSFXFNAME(:)
REAL(KIND=JPRB) :: ZCO_EZO
LOGICAL :: LLNEWCL(2), LLELAM, LLOPENED
TYPE(TNAMFPD) :: YLNAMFPD
TYPE(TNAMFPG) :: YLNAMFPG
TYPE(TNAMFPV) :: YLNAMFPV
TYPE(TNAMFPF) :: YLNAMFPF
TYPE(TNAMFPIOS) :: YLNAMFPIOS
TYPE(TNAMFPOBJ) :: YLNAMFPOBJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM), EXTERNAL :: I_SYSTEM
!     ------------------------------------------------------------------

#include "sufpc.intfb.h"
#include "updtrans.intfb.h"
#include "sufpmodelgeo.intfb.h"
#include "sufpusergeo.intfb.h"
#include "sumpfpos.intfb.h"
#include "sufpvert.intfb.h"
#include "suafn.intfb.h"
#include "sufpcnt.intfb.h"
#include "sufpd.intfb.h"
#include "sufpg.intfb.h"
#include "sufpv.intfb.h"
#include "sufpf.intfb.h"
#include "sufpfilters.intfb.h"
#include "sufpios.intfb.h"
#include "sufpioh.intfb.h"
#include "sufpgeometry.intfb.h"
#include "sufpsc2_dep.intfb.h"
#include "sufpsc2.intfb.h"
#include "sufpezo.intfb.h"
#include "sufpwfpbuf.intfb.h"
#include "sufpwide.intfb.h"

#include "abor1.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUBFPOS',0,ZHOOK_HANDLE)
ASSOCIATE(CFPFMT=>YLNAMFPOBJ%CFPFMT, CFPDOM=>YLNAMFPOBJ%CFPDOM, CFPDIR=>YLNAMFPOBJ%CFPDIR, NFPGRIB=>YLNAMFPOBJ%NFPGRIB, &
 & NFRFPOS=>YLNAMFPOBJ%NFRFPOS, NFPOSTS=>YLNAMFPOBJ%NFPOSTS, NFPOSTSMIN=>YLNAMFPOBJ%NFPOSTSMIN, &
 & CFPMONIPATH_IN=>YLNAMFPOBJ%CFPMONIPATH_IN, CFPMONIPATH_OUT=>YLNAMFPOBJ%CFPMONIPATH_OUT)

CALL GSTATS(2001,0)

!     ------------------------------------------------------------------

IF (KFPOS == 0) THEN
  CALL ABOR1('SUBFPOS : YOU ARE NOT SUPPOSED TO CONSTRUCT A POST-PROCESSOR WITHOUT POST-PROCESSING INTENTION !!')
ENDIF

!            0. CHANGE TO ALTERNATIVE NAMELIST FILE
!               ===================================

IF (PRESENT(CDNAM)) THEN
  INQUIRE(NULNAM,OPENED=LLOPENED)
  IF (LLOPENED) CLOSE(NULNAM)
  OPEN(NULNAM,FILE=CDNAM,ACTION='READ',IOSTAT=IOS)
  IF (IOS /= 0) THEN
    CALL ABOR1('SUBFPOS FAILED TO OPEN NAMELIST FILE '//TRIM(CDNAM))
  ENDIF
ENDIF

!            1. GENERAL CONFIGURATION OF POST-PROCESSOR AND ASSOCIATED DOMAINS LABELS
!               =====================================================================

! Partial reading of namfpc :
CALL SUFPC(LDPRINT=.TRUE.,YDNAMFPOBJ=YLNAMFPOBJ,YDNAMFPSCI=YDFPOS%YNAMFPSCI,YDNAMFPINT=YDFPOS%YNAMFPINT,KFPCONF=KFPOS)

IF (PRESENT(YDPOSGEOMETRY).OR.PRESENT(YDXPOSGEOMETRY)) THEN
  IFPCONF=MAX(1,MIN(2,KFPOS))
ELSEIF (CFPFMT == 'MODEL') THEN
  IFPCONF=0
ELSE
  IFPCONF=MAX(1,MIN(2,KFPOS))
ENDIF

DO J=1,SIZE(CFPDOM)
  IF (CFPDOM(J) == ' ') EXIT
ENDDO
IF (J == 1 .AND. IFPCONF /= 0) THEN
  CALL ABOR1('SUBFPOS : CFPDOM CANNOT BE AN ARRAY OF BLANK STRINGS !')
ELSE
  IFPDOM=MAX(1,J-1)
ENDIF

WRITE(NULOUT,'('' == Full-Pos constructor : setup control variables == '')')
CALL SUFPCNT(NFRFPOS,NFPOSTS,NFPOSTSMIN,CFDIRLST,CNPPATH,CFPMONIPATH_IN,CFPMONIPATH_OUT,CFPNCF,IFPCONF,YDFPOS%YFPCNT,CDNAM)


!            2. HORIZONTAL GEOMETRIES ASPECTS
!               =============================

WRITE(NULOUT,'('' == Full-Pos constructor : setup geometries == '')')
IF (IFPCONF == 0) THEN
  IF (IFPDOM > 1) THEN
    CALL ABOR1('SUBFPOS : ONLY 1 DOMAIN POSSIBLE WITH CFPFMT=MODEL !')
  ELSE
    CALL UPDTRANS(YDHOSTGEOMETRY%YRDIM%NRESOL,LLELAM)
    WRITE(UNIT=NULOUT,FMT='('' CFPFMT = '', A5,'' NFPDOM = '',I2)') CFPFMT,IFPDOM
    ALLOCATE(YDFPOS%YFPGEOMETRY%YFPUSERGEO(IFPDOM))
    WRITE(UNIT=NULOUT,FMT='('' DOMAIN STRUCTURE : '')')
    CALL SUFPMODELGEO(LLELAM,YDHOSTGEOMETRY,CFPDOM(IFPDOM),YDFPOS%YFPGEOMETRY%YFPUSERGEO(IFPDOM),YDHOSTGEOMETRY,IFPCONF, &
     & YDFPOS%YNAMFPSCI%NFPBOYD)
    ! use the negative value to be sure the potentially optimized value will be used
    IFPROMA=-ABS(YDHOSTGEOMETRY%YRDIM%NPROMA)
    WRITE(UNIT=NULOUT,FMT='('' '')')
    ! Initialize Fullpos horizontal output geometry
    YDFPOS%YFPGEOMETRY%LFPOSHOR=.FALSE.
    YDFPOS%YFPGEOMETRY%YFPGEO_DEP%NFPRGPG=YDHOSTGEOMETRY%YRGEM%NGPTOTG
    YDFPOS%YFPGEOMETRY%YFPGEO%NFPRGPG=YDHOSTGEOMETRY%YRGEM%NGPTOTG
    CALL SUMPFPOS(' ',YDHOSTGEOMETRY%YRGEM%NGPTOTG,YDHOSTGEOMETRY%YRMP%NGLOBALPROC,IFPROMA,YDFPOS%YFPGEOMETRY%YFPGEO%NFPRGPNUM, &
     & YDFPOS%YFPGEOMETRY%YFPGEO%NFPRGPL,YDFPOS%YFPGEOMETRY%YFPGEO%NFPRGPLX,YDFPOS%YFPGEOMETRY%YFPGEO%NFPRGPIND, &
     & YDFPOS%YFPGEOMETRY%YFPGEO%NFPROMA,YDFPOS%YFPGEOMETRY%YFPGEO%NFPBLOCS,YDFPOS%YFPGEOMETRY%YFPGEO%NFPEND)
    ! Map factor can be of some help ...
    ALLOCATE(YDFPOS%YFPGEOMETRY%YFPGEO%RFPGM(YDFPOS%YFPGEOMETRY%YFPGEO%NFPROMA,YDFPOS%YFPGEOMETRY%YFPGEO%NFPBLOCS))
    DO J=1,YDFPOS%YFPGEOMETRY%YFPGEO%NFPBLOCS
      YDFPOS%YFPGEOMETRY%YFPGEO%RFPGM(1:YDFPOS%YFPGEOMETRY%YFPGEO%NFPEND(J),J)= &
       & YDHOSTGEOMETRY%YRGSGEOM(J)%GM(1:YDFPOS%YFPGEOMETRY%YFPGEO%NFPEND(J)) ! array sized ngptot, not nproma*ngpblks :-(
    ENDDO
  ENDIF
ELSE
  IF (PRESENT(YDPOSGEOMETRY).AND..NOT.PRESENT(YDXPOSGEOMETRY)) THEN
    IF (IFPDOM > 1) THEN
      CALL ABOR1('SUBFPOS : ONLY 1 DOMAIN POSSIBLE WITH YDPOSGEOMETRY !')
    ELSE
      CALL UPDTRANS(YDPOSGEOMETRY%YRDIM%NRESOL,LLELAM)
      IF (LLELAM) THEN
        CFPFMT='LELAM'
      ELSE
        CFPFMT='GAUSS'
      ENDIF
      WRITE(UNIT=NULOUT,FMT='('' CFPFMT = '', A5,'' NFPDOM = '',I2)') CFPFMT,IFPDOM
      ALLOCATE(YDFPOS%YFPGEOMETRY%YFPUSERGEO(IFPDOM))
      WRITE(UNIT=NULOUT,FMT='('' DOMAIN STRUCTURE : '')')
      CALL SUFPMODELGEO(LLELAM,YDPOSGEOMETRY,CFPDOM(IFPDOM),YDFPOS%YFPGEOMETRY%YFPUSERGEO(IFPDOM),YDHOSTGEOMETRY,IFPCONF, &
       & YDFPOS%YNAMFPSCI%NFPBOYD)
      ! must keep to the provided geometries (though not always needed for intermediate grid - IFPROMA_DEP)
      ! use the negative value to be sure the potentially optimized value will be used
      IFPROMA_DEP=-ABS(YDHOSTGEOMETRY%YRDIM%NPROMA)
      IFPROMA=-ABS(YDPOSGEOMETRY%YRDIM%NPROMA)
      IFPWIDE_NAM=0
    ENDIF
  ELSEIF (PRESENT(YDXPOSGEOMETRY).AND..NOT.PRESENT(YDPOSGEOMETRY)) THEN
    IF (IFPDOM  > SIZE(YDXPOSGEOMETRY)) THEN
      WRITE(UNIT=NULOUT,FMT='('' SIZE(YDXPOSGEOMETRY) = '',I2)') SIZE(YDXPOSGEOMETRY)
      CALL ABOR1('SUBFPOS : SIZE OF YDXPOSGEOMETRY LESS THAN EXPECTED NUMBER OF POST-PROCESSING GEOMETRIES !')
    ELSE
      CALL UPDTRANS(YDXPOSGEOMETRY(1)%YRDIM%NRESOL,LLELAM)
      IF (LLELAM) THEN
        CFPFMT='LELAM'
      ELSE
        CFPFMT='GAUSS'
      ENDIF
      WRITE(UNIT=NULOUT,FMT='('' CFPFMT = '', A5,'' NFPDOM = '',I2)') CFPFMT,IFPDOM
      DO J=2,SIZE(YDXPOSGEOMETRY)
        CALL UPDTRANS(YDXPOSGEOMETRY(J)%YRDIM%NRESOL,LLELAM)
        IF (LLELAM .AND. CFPFMT /= 'LELAM') THEN
          CALL ABOR1('SUBFPOS : MIXING GLOBAL AND LAM GEOMETRIES IS NOT YET POSSIBLE !')
        ENDIF
      ENDDO
      ALLOCATE(YDFPOS%YFPGEOMETRY%YFPUSERGEO(IFPDOM))
      DO J=1,IFPDOM
        WRITE(UNIT=NULOUT,FMT='('' DOMAIN STRUCTURE '',I2,'' : '')') J
        CALL SUFPMODELGEO(LLELAM,YDXPOSGEOMETRY(J),CFPDOM(J),YDFPOS%YFPGEOMETRY%YFPUSERGEO(J),YDHOSTGEOMETRY,IFPCONF, &
         & YDFPOS%YNAMFPSCI%NFPBOYD)
        WRITE(UNIT=NULOUT,FMT='('' '')')
      ENDDO
      ! Obviously fullpos buffer shape
      IFPROMA_DEP=-ABS(YDHOSTGEOMETRY%YRDIM%NPROMA)
      IFPROMA=-ABS(YDHOSTGEOMETRY%YRDIM%NPROMA)
      CALL SUFPSC2_DEP(IFPROMA_DEP)
      CALL SUFPSC2(IFPROMA,IFPWIDE_NAM)
    ENDIF
  ELSEIF (.NOT.(PRESENT(YDPOSGEOMETRY).OR.PRESENT(YDXPOSGEOMETRY))) THEN
    WRITE(UNIT=NULOUT,FMT='('' CFPFMT = '', A5,'' NFPDOM = '',I2)') CFPFMT,IFPDOM
    ! horizontal geometries-related namelists
    WRITE(NULOUT,'(''- Set up F-post processing, Horizontal geometries -'')')
    CALL SUFPD(KFPOS,YDHOSTGEOMETRY,CFPFMT,CFPDOM(1:IFPDOM),YLNAMFPD)
    CALL SUFPG(YDHOSTGEOMETRY,CFPFMT,CFPDOM(1:IFPDOM),YLNAMFPD,YLNAMFPG)
    ! user geoemtries
    ALLOCATE(YDFPOS%YFPGEOMETRY%YFPUSERGEO(IFPDOM))
    DO J=1,IFPDOM
      WRITE(UNIT=NULOUT,FMT='('' DOMAIN STRUCTURE '',I2,'' : '')') J
      CALL SUFPUSERGEO(YDFPOS%YNAMFPSCI,CFPFMT,CFPDOM(J),YLNAMFPD,YLNAMFPG,KFPOS,YDHOSTGEOMETRY,J, &
       & YDFPOS%YFPGEOMETRY%YFPUSERGEO(J))
      WRITE(UNIT=NULOUT,FMT='('' '')')
    ENDDO
    ! use the negative value to be sure the potentially optimized value will be used
    IFPROMA_DEP=-ABS(YDHOSTGEOMETRY%YRDIM%NPROMA)
    IFPROMA=-ABS(YDHOSTGEOMETRY%YRDIM%NPROMA)
    IF (LFPOSBUF(YDFPOS%YFPGEOMETRY)) THEN
      CALL SUFPSC2_DEP(IFPROMA_DEP)
      CALL SUFPSC2(IFPROMA,IFPWIDE_NAM)
    ELSE
      IFPWIDE_NAM=0
    ENDIF
  ELSE
    CALL ABOR1('SUBFPOS : YDPOSGEOMETRY AND YDXPOSGEOMETRY ARE EXCLUSIVE !')
  ENDIF
  IF (ANY(YDFPOS%YFPGEOMETRY%YFPUSERGEO(:)%LFPBIPER)) THEN
    CALL SUFPEZO(ZCO_EZO)
  ENDIF
  ! Initialize Fullpos horizontal output geometries
  CALL SUFPGEOMETRY(YDFPOS%YNAMFPINT,YDFPOS%YNAMFPSCI,YDFPOS%YFPGEOMETRY,YDHOSTGEOMETRY,IFPROMA_DEP,IFPROMA,ZCO_EZO)
  ! Initialize Fullpos horizontal interpolator
  IF (LFPOSBUF(YDFPOS%YFPGEOMETRY)) THEN
    WRITE(NULOUT,'('' == Full-Pos constructor : setup halo for horizontal interpolations == '')')
    CALL SUFPWIDE(YDFPOS%YFPGEOMETRY,YDFPOS%YFPSTRUCT,YDHOSTGEOMETRY,IFPWIDE_NAM,YDFPOS%YFPGEOMETRY%YFPGEO_DEP)
    WRITE(NULOUT,'('' == Full-Pos constructor : setup standard interpolation weights == '')')
    CALL SUFPWFPBUF(YDFPOS%YNAMFPINT,YDFPOS%YFPGEOMETRY%YFPGEO_DEP,YDFPOS%YFPWSTD,YDFPOS%YFPSTRUCT,YDHOSTGEOMETRY)
  ENDIF
ENDIF
! Tag of the input horizontal resolution :
YDFPOS%YFPGEOMETRY%NMDLRESOL=YDHOSTGEOMETRY%YRDIM%NRESOL

IF(LTRACEFP) THEN
  WRITE(NULOUT,'(A)') ' '
  WRITE(UNIT=NULOUT,FMT='('' LFPOSBUF = '',L2,'' LFPOSHOR = '',L2,'' NFPDISTRIB = '',I2,'' NFPRGPG_DEP = '',I6, &
   & '' NFPRGPG = '',I6)') LFPOSBUF(YDFPOS%YFPGEOMETRY), YDFPOS%YFPGEOMETRY%LFPOSHOR, &
   & YDFPOS%YFPGEOMETRY%YFPUSERGEO(1)%NFPDIST, YDFPOS%YFPGEOMETRY%YFPGEO_DEP%NFPRGPG, YDFPOS%YFPGEOMETRY%YFPGEO%NFPRGPG 
ENDIF

!            3. ASSOCIATED VERTICAL ETA COORDINATE
!               ==================================

WRITE(NULOUT,'('' == Full-Pos constructor : setup output vertical eta coordinates == '')')
! either a specific coordinate is specified in argument, or a namelist is read with the model coordinate as default
IF (PRESENT(YDPOSVAB)) THEN
  CALL SUFPVERT(YDFPOS%YFPVAB,YDVAB=YDPOSVAB)
ELSE
  CALL SUFPV(YDHOSTGEOMETRY,YLNAMFPV)
  CALL SUFPVERT(YDFPOS%YFPVAB,YDNAMFPV=YLNAMFPV)
ENDIF

!            4. I/O HANDLINGS
!               =============

ALLOCATE(CLFPFN(IFPDOM))
ALLOCATE(CLFPCLIFNAME(IFPDOM))
ALLOCATE(CLFPSFXFNAME(IFPDOM))
WRITE(UNIT=NULOUT,FMT='('' CFPDIR = '',A)') TRIM(CFPDIR)
CALL SUFPIOS(NFPGRIB,YDFPOS%YNAMFPSCI%NFPSURFEX,CFPDIR,YDFPOS%YFPGEOMETRY%YFPUSERGEO(:)%CFPDOM, &
 & CLFPFN,CLFPCLIFNAME,CLFPSFXFNAME,YLNAMFPIOS)
WRITE(NULOUT,'('' == Full-Pos constructor : setup I/O handlings == '')')
CALL SUFPIOH(YDFPOS%YNAMFPSCI%NFPSURFEX,YLNAMFPIOS,CLFPFN,CLFPCLIFNAME,CLFPSFXFNAME,YDFPOS%YFPGEOMETRY,YDFPOS%YFPVAB,YDFPOS%YFPIOH)


!            5. FIELDS DESCRIPTORS
!               ==================

WRITE(NULOUT,'(''--- Set up Fullpos fields descriptors '')')
!*    Initialize Fullpos field names and descriptors
LLNEWCL(1)= .NOT.YDMODEL%YRML_GCONF%YGFL%YL%LGP &
  & .AND..NOT.YDMODEL%YRML_GCONF%YGFL%YI%LGP &
  & .AND..NOT.YDMODEL%YRML_GCONF%YGFL%YA%LGP &
  & .AND..NOT.YDMODEL%YRML_PHY_SLIN%YRPHNC%LEPCLD2 &
  & .AND.YDMODEL%YRML_PHY_SLIN%YRPHNC%LENCLD2
LLNEWCL(2)=.NOT.YDMODEL%YRML_GCONF%YGFL%YA%LGP &
  & .AND.YDMODEL%YRML_PHY_SLIN%YRPHNC%LEPCLD2
CALL SUAFN(YDFPOS%YNAMFPSCI,YDFPOS%YAFN,LLNEWCL,YDMODEL%YRML_GCONF%YGFL,YDMODEL%YRML_DYN%YRDYNA%LNHDYN, &
  &        IFPCONF,YDMODEL%YRML_CHEM%YRCOMPO,YDMODEL%YRML_PHY_RAD%YREAERATM,CLFPSFXFNAME)
! -------- this could be changed later (moving it down to cprep3/fullpos_drv) -------------


!            6. SPECTRAL FILTERS
!               ================

CALL SUFPF(YLNAMFPF,IFPCONF,YDFPOS%YFPGEOMETRY%YFPUSERGEO)
WRITE(NULOUT,'('' == Full-Pos constructor : setup spectral filters == '')')
CALL SUFPFILTERS(YDFPOS%YNAMFPSCI,YLNAMFPF,YDHOSTGEOMETRY,YDMODEL%YRML_PHY_G%YRDPHY%NCSNEC,IFPDOM,IFPCONF,YDFPOS%YAFN%TFP_DYNDS(:)%ISF, &
 & YDFPOS%YAFN%TFP_DYNDS(:)%CLNAME,YDFPOS%YFPFILTERS)
! Remark : later, the filtering matrixes for Arpege output lat-lon domains could be replaced
! by pointers, allocated as target in a third-party object, and computed in an associated external method
! on the basis of a multi-objects post-processor.
! That would allow to make one post-processor object per post-processing time step, where the list of output subdomain
! could change : all that to spare cpu time.


!            7. CHANGE BACK TO DEFAULT NAMELIST FILE
!               ====================================

IF (PRESENT(CDNAM)) THEN
  CLOSE(NULNAM)
ENDIF


!            8. PREPARE MONITORING
!               ==================

IF (YDFPOS%YFPIOH%YNAMFPIOS%NFPWRITE==1 .AND. YDFPOS%YFPCNT%CFPNCF /= ' ' .AND. &
  & TRIM(CFPMONIPATH_OUT) /= '.' .AND. MYPROC == 1) THEN
  CALL EXECUTE_COMMAND_LINE('/bin/mkdir -p '//TRIM(CFPMONIPATH_OUT),EXITSTAT=IERR,CMDSTAT=ICMDSTAT)
  IF (ICMDSTAT /= 0) THEN
    WRITE(NULOUT,'(''SUBFPOS:SYSTEM CALL FOR mkdir RETURNED STATUS '',I2)') IERR
    CALL ABOR1('SUBFPOS : FAILED IN EXECUTE_COMMAND_LINE /bin/mkdir')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

CALL GSTATS(2001,1)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUBFPOS',1,ZHOOK_HANDLE)

END SUBROUTINE SUBFPOS
