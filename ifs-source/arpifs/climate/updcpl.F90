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

SUBROUTINE UPDCPL(YDGEOMETRY,YDSURF,YDML_AOC,YDRIP,YDPHY1,KGP)

!**** *UPDCPL*

!     PURPOSE.
!     --------

!     Updates the climatological data on coupled mode

!**   INTERFACE.
!     ----------

!     CALL UPDCPL(KGP)

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     UPDCAL
!     CHIEN
!     FAITOU
!     FADIES
!     FACILE
!     FAIRME

!     AUTHORS.
!     --------

!        CNRM/GMGEC/EAC/Jean-Philippe Piedelievre 09-1999

!     MODIFICATIONS.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!        original 99-01-28 from updcli
!     30-10-2002 by P. Marquet : transcoded to *.f90
!        M.Hamrud      01-Jul-2006 Revised surface fields
!     30-08-2006 by A. Alias   : RZHMER replaced by RZHZ0M*RZ0MER
!     05-09-2006 by A. Alias   : RZHGLA replaced by RZHZ0G*RZ0GLA
!     Apr 2008  K. Yessad: use DISGRID instead of DISGRID_C + cleanings
!     Aug 2009 by A.Alias      : ocean currents added (J.F. Gueremy)
!                                No Correction of SST by the orography
!                                check pole in case of LELAM (S. Calmanti)
!                                set in agreement with external updcli (M.Deque)
!        R. El Khatib 03-Jun-2010 Remove dead call to clim_import
!     K. Yessad (Jan 2010): remove useless variables.
!     K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!     R. El Khatib : 23-Apr-2010 use disgrid_mod instead of disgrid
!     G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     K. Yessad (July 2014): Move some variables.
!     E.Dutra/G.Arduini (Jan 2018): Snow multi-layer, changes in SP_SG
!     ------------------------------------------------------------------

USE MODEL_ATMOS_OCEAN_COUPLING_MOD , ONLY : MODEL_ATMOS_OCEAN_COUPLING_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT   
USE YOMVERT  , ONLY : VP00
USE YOMCST   , ONLY : RG
USE YOMPHY1  , ONLY : TPHY1
USE YOMCT0   , ONLY : NQUAD, LELAM
USE YOMCT3   , ONLY : NSTEP
USE YOMRIP0  , ONLY : NINDAT
USE YOMRIP   , ONLY : TRIP
USE YOMMP0   , ONLY : MYPROC, NPROC
USE DISGRID_MOD, ONLY : DISGRID_SEND, DISGRID_RECV
USE YOMCPL   , ONLY : FCSTSU, FCSTSV, FCCHAS, FCRSOS, FCCHLL, FCCHLN, FCCHSS
USE YOM_OAS  , ONLY : CPHAN
USE PAR_COU  , ONLY : JPFLDO2A

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)     ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)        ,INTENT(INOUT) :: YDSURF
TYPE(MODEL_ATMOS_OCEAN_COUPLING_TYPE),INTENT(INOUT):: YDML_AOC
TYPE(TPHY1)        ,INTENT(INOUT) :: YDPHY1
TYPE(TRIP)         ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KGP
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JPCCL
PARAMETER(JPCCL=20)

REAL(KIND=JPRB) ::      ZREALG(YDGEOMETRY%YRDIM%NDLON*(YDGEOMETRY%YRDIM%NDGLG+2))
REAL(KIND=JPRB) ::      ZFIELDBUF(YDGEOMETRY%YRGEM%NGPTOTG),ZLSMG(YDGEOMETRY%YRGEM%NGPTOTG),ZLSM(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ::      ZFIELD(YDGEOMETRY%YRGEM%NGPTOT,KGP)
REAL(KIND=JPRB) ::      ZGEOP(YDGEOMETRY%YRGEM%NGPTOT)

INTEGER(KIND=JPIM) ::   ILMOIS(12),IDATEF(11)
INTEGER(KIND=JPIM) ::   INFOS(3)

CHARACTER (LEN = 3)  :: CLJOBNAM_R
CHARACTER (LEN = 16) :: CLNOMC
CHARACTER (LEN = 20) :: CLNOMF,CLOST(JPCCL)

!CLIM
 
INTEGER(KIND=JPIM) ::      IA0
INTEGER(KIND=JPIM) ::      IAN
INTEGER(KIND=JPIM) ::      IEND
INTEGER(KIND=JPIM) ::      IINF
INTEGER(KIND=JPIM) ::      IJ0
INTEGER(KIND=JPIM) ::      IJOUR
INTEGER(KIND=JPIM) ::      IJT1
INTEGER(KIND=JPIM) ::      IJT2
INTEGER(KIND=JPIM) ::      IM0
INTEGER(KIND=JPIM) ::      IMOIS
INTEGER(KIND=JPIM) ::      IMT1
INTEGER(KIND=JPIM) ::      IMT2
INTEGER(KIND=JPIM) ::      INBARI
INTEGER(KIND=JPIM) ::      INBARP
INTEGER(KIND=JPIM) ::      INDEX
INTEGER(KIND=JPIM) ::      INIMES
INTEGER(KIND=JPIM) ::      INULCL
INTEGER(KIND=JPIM) ::      INULCL1
INTEGER(KIND=JPIM) ::      INULCL2
INTEGER(KIND=JPIM) ::      INUM
INTEGER(KIND=JPIM) ::      IBL
INTEGER(KIND=JPIM) ::      IREG
INTEGER(KIND=JPIM) ::      IREP
INTEGER(KIND=JPIM) ::      IROF
INTEGER(KIND=JPIM) ::      IST
INTEGER(KIND=JPIM) ::      ISTEPCPL
INTEGER(KIND=JPIM) ::      ITIM
INTEGER(KIND=JPIM) ::      IZOA
INTEGER(KIND=JPIM) ::      IZOI
INTEGER(KIND=JPIM) ::      IZOT
INTEGER(KIND=JPIM) ::      IZOU
INTEGER(KIND=JPIM) ::      IZOV
INTEGER(KIND=JPIM) ::      IZST

INTEGER(KIND=JPIM) ::      J
INTEGER(KIND=JPIM) ::      JCSS
INTEGER(KIND=JPIM) ::      JKGLO
INTEGER(KIND=JPIM) ::      JM
INTEGER(KIND=JPIM) ::      JNULCL
INTEGER(KIND=JPIM) ::      JROF
INTEGER(KIND=JPIM) ::      JV

REAL(KIND=JPRB) ::      ZEPS
REAL(KIND=JPRB) ::      Z_CLIM_OK

REAL(KIND=JPRB) ::      ZPOID1
REAL(KIND=JPRB) ::      ZPOID2

REAL(KIND=JPRB) ::      ZWXMER
REAL(KIND=JPRB) ::      ZWXSUR
REAL(KIND=JPRB) ::      ZWXGLA
REAL(KIND=JPRB) :: ZRANO , ZRDNO
 
LOGICAL ::      LLERFA
LOGICAL ::      LLEXIST
LOGICAL ::      LLFICP
LOGICAL ::      LLIMST
LOGICAL ::      LLNOMM
LOGICAL ::      LLOPENED
LOGICAL ::      LLMEGL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "fcttim.func.h"

!-----------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "sipc_read_model.intfb.h"
#include "slab.intfb.h"
#include "updcal.intfb.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDCPL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  YDCSGLEG=>YDGEOMETRY%YRCSGLEG, YDVAB=>YDGEOMETRY%YRVAB, &
 & YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,  YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
  & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDMCC=>YDML_AOC%YRMCC,  &
  & YDCOM=>YDML_AOC%YRCOM)

ASSOCIATE(NTVGLA=>YDPHY1%NTVGLA, NTVMER=>YDPHY1%NTVMER, RZ0GLA=>YDPHY1%RZ0GLA, &
 & RZHZ0M=>YDPHY1%RZHZ0M, EMMMER=>YDPHY1%EMMMER, ALBMER=>YDPHY1%ALBMER, &
 & RZHZ0G=>YDPHY1%RZHZ0G, RD2GLA=>YDPHY1%RD2GLA, RD1=>YDPHY1%RD1, &
 & RD2MER=>YDPHY1%RD2MER, TMERGL=>YDPHY1%TMERGL, ALBGLA=>YDPHY1%ALBGLA, &
 & RZ0MER=>YDPHY1%RZ0MER, EMMGLA=>YDPHY1%EMMGLA, &
 & GTTLIN=>YDCOM%GTTLIN, LOMLDTH=>YDCOM%LOMLDTH, NVCOM=>YDCOM%NVCOM, &
 & OMLDTH=>YDCOM%OMLDTH, SSTMSK=>YDCOM%SSTMSK, SSTPRE=>YDCOM%SSTPRE, &
 & TRAFLX=>YDCOM%TRAFLX, &
 & NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG, &
 & NDLON=>YDDIM%NDLON, NPROMA=>YDDIM%NPROMA, NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, NHTYP=>YDGEM%NHTYP, &
 & NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG, NSTTYP=>YDGEM%NSTTYP, &
 & RLOCEN=>YDGEM%RLOCEN, RMUCEN=>YDGEM%RMUCEN, RSTRET=>YDGEM%RSTRET, &
 & LCURR=>YDMCC%LCURR, LMCC02=>YDMCC%LMCC02, LMCC03=>YDMCC%LMCC03, &
 & LMCC05=>YDMCC%LMCC05, NOACOMM=>YDMCC%NOACOMM, &
 & NSTADD=>YDRIP%NSTADD, TSTEP=>YDRIP%TSTEP, &
 & SD_VF=>YDSURF%SD_VF, SD_VP=>YDSURF%SD_VP, SD_VV=>YDSURF%SD_VV, &
 & SP_RR=>YDSURF%SP_RR, SP_SB=>YDSURF%SP_SB, SP_SG=>YDSURF%SP_SG, &
 & YSD_VF=>YDSURF%YSD_VF, YSD_VP=>YDSURF%YSD_VP, YSD_VV=>YDSURF%YSD_VV, &
 & YSD_VVD=>YDSURF%YSD_VVD, YSP_RR=>YDSURF%YSP_RR, YSP_SB=>YDSURF%YSP_SB, &
 & YSP_SBD=>YDSURF%YSP_SBD, YSP_SG=>YDSURF%YSP_SG,YSP_SGD=>YDSURF%YSP_SGD)

!-----------------------------------------------------------------------

WRITE(NULOUT,*)  ' ---------------------- '
WRITE(NULOUT,*)  ' UPDCPL ON PROC ',MYPROC
WRITE(NULOUT,*)  ' VERSION CORRIGEE GELATO XCLIM'
WRITE(NULOUT,*)  ' ---------------------- '

!     1. SETTING CONSTANT VALUES.
!     ---------------------------

!*    1.1 Arpege files

!SLAB MODEL (LMCC05=.TRUE.)
CLOST(1)   =  'PROFMIXLAY.DEPTH' ! epaisseur de la couche de melange oceanique.
CLOST(2)   =  'PROFGRT.THERMOCL' ! gradient thermique de la thermocline.

LLNOMM=.TRUE.
LLERFA=.TRUE.
LLIMST=.FALSE.
INIMES=1
INBARP=14
CLNOMC='Cadre.Clim'
LLFICP=.FALSE.
INULCL1=0
INULCL2=0
LLOPENED=.FALSE.
LLEXIST =.FALSE.

!* Calendar

IJ0=NDD(NINDAT)
IM0=NMM(NINDAT)
IA0=NCCAA(NINDAT)
IJOUR=1
IMOIS=1
IAN  =1
CALL UPDCAL (IJ0,IM0,IA0,NSTADD,IJOUR,IMOIS,IAN,ILMOIS,NULOUT)
IF (IJOUR > 15) THEN
  IMT1=IMOIS
  IMT2=1+MOD(IMOIS,12)
  IJT1=15
  IJT2=15+ILMOIS(IMT1)
ELSE
  IMT1=1+MOD(IMOIS+10,12)
  IMT2=IMOIS
  IJT1=15-ILMOIS(IMT1)
  IJT2=15
ENDIF
ZPOID1=REAL(IJT2-IJOUR,JPRB)/REAL(IJT2-IJT1,JPRB)
ZPOID2=1.0_JPRB-ZPOID1

!*    1.4 Physics


INUM=0
ZFIELD=0.0_JPRB

IF (LMCC03) THEN

! SST
  IZOT=1
! SEA ICE COVER
  IZOI=2
! ALBEDO
  IZOA=3
! U-CURRENT
IZOU=4
! V-CURRENT
IZOV=5
!    Define the number of field exchanged from ocean to atmosphere
IF (LCURR) THEN
JPFLDO2A=5
ELSE
JPFLDO2A=3
ENDIF

!*        Couplage OASIS

!*    2.8 Reading surface temperature and sea-ice mask from coupler

!CLIM com (LMCC03)
  Z_CLIM_OK=0
  ISTEPCPL=NSTEP+1
  IF(NSTEP == 0) ISTEPCPL=0
!ENDCLIM
!MPI1
ITIM=(NSTEP+1)*TSTEP
IF(NSTEP == 0) ITIM=0
IF (NOACOMM == 5) THEN
  WRITE(NULOUT,'(A,I8)')'Getting variables from OASIS3 for time : ',ITIM
ENDIF
!ENDMPI1

!* 1) Sea-ice surface temperature
!* 2) sea-ice index or ice fraction
!* 3) sea-ice surface albedo

! COUPLAGE A OASIS

! RECEPTION DES FLUX VENANT D'OASIS

IF (NOACOMM == 5) THEN

   CALL ABOR1('UPDCPL: ABOR1 CALLED OASIS3 NOT IMPLEMENTED')

ELSE

  DO JV=1,JPFLDO2A
    INUM=INUM+1
    IF (MYPROC == 1) THEN

!* Reading of input field sea-surface-temperature SISUTESU
!* or sea-ice cover SIICECOV in shared-memory segment
!* jv is the index of the field in total number of fields jpfldo2a

      IF    ( CPHAN == 'SIPC') THEN

        CALL SIPC_READ_MODEL(JV,NGPTOTG,CLJOBNAM_R,INFOS,ZFIELDBUF)

      ENDIF

      CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,ZFIELD(:,JV))
    ELSE
      CALL DISGRID_RECV(YDGEOMETRY,1,1,ZFIELD(1,JV),INUM)
    ENDIF
  ENDDO

ENDIF

!  CORRECTION OF SURFACE TEMPERATURE BY GIBBS OROGRAPHY

  INUM=INUM+1
  IF (MYPROC == 1) THEN

!*  Opening the local climatic boundary conditions file
!*  used by the external module updcli at the beginning
!*  of the run.

    DO JNULCL=1,99
      INQUIRE(UNIT=JNULCL,OPENED=LLOPENED)
      IF (.NOT.LLOPENED) EXIT
    ENDDO
    INULCL=JNULCL
    CALL FAITOU (IREP,INULCL,LLNOMM,'ICMSH0123INIT','OLD',&
     & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  
    ZEPS=1.E-10_JPRB
    IINF=0
    IF (.NOT.LELAM) THEN
      CALL CHIEN (CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,&
       & NSMAX,NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,&
       & VP00,YDVAB%VALH,YDVAB%VBH,NQUAD,IINF,NDGSAG,NDGENG,&
       & ZEPS,LLFICP,NULOUT)
    ELSE
      LLFICP=.FALSE. 
    ENDIF 

!*    Indices for horizontal loops
    IF (LLFICP) THEN
      INDEX=NLOENG(1)+1
    ELSE
      INDEX=1
    ENDIF
    CALL FACILE(IREP,INULCL,'SURF',1,'GEOPOTENTIEL',ZREALG,.FALSE.)
    CALL FAIRME (IREP,INULCL,'KEEP')
    ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
    CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,ZGEOP)
  ELSE
    CALL DISGRID_RECV(YDGEOMETRY,1,1,ZGEOP,INUM)
  ENDIF

ELSEIF (LMCC05) THEN

!  SST COM
  IZOT=1
!  SST STATISTIQUE LUE DANS LE FICHIER BCOND
  IZST=2

!* Couplage COM

!* Reads the buffers

! On  recupere TSTAT 

  DO JKGLO=1,NGPTOT,NPROMA
    IST=1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    DO JROF=IST,IEND
      IROF=JROF+JKGLO-1
      ZFIELD(IROF,IZST)=SD_VP(JROF,YSD_VP%YTPC%MP,IBL)
      ZLSM(IROF)=SD_VF(JROF,YSD_VF%YLSM%MP,IBL)
    ENDDO
  ENDDO

  IF (NSTEP == 0) THEN

! Initialisation

! au premier appel, la profondeur de melange et le gradT thermocline
! ------------------------------------------------------------------
    INQUIRE(FILE='ICMCM0123INIT',EXIST=LLEXIST)

    IF (LLEXIST) THEN
! sont lus dans le fichier historique s'il existe.
! (copie de sauvegarde d'un passge anterieur)
! Attention : le fichier sera detruit apres lecture.
! (STATUS='DELETE' dans FAIRME)
! ------------------------------------------------
      IF (MYPROC == 1) THEN
!*    Opening Arpege files : month 1 for the stat-dyn mask
        WRITE(UNIT=CLNOMF,FMT='(''Const.Clim.'',I2.2)') IMT1
        DO JNULCL=1,99
          INQUIRE(UNIT=JNULCL,OPENED=LLOPENED)
          IF (.NOT.LLOPENED) EXIT
        ENDDO
        INULCL=JNULCL
        CALL FAITOU (IREP,INULCL,LLNOMM,CLNOMF,'OLD',&
         & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  

        WRITE (NULOUT,*) ' FICHIER ',CLNOMF,' OUVERT'

        ZEPS=1.E-10_JPRB
        IINF=0
        CALL CHIEN (CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,&
         & NSMAX,NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,&
         & VP00,YDVAB%VALH,YDVAB%VBH,NQUAD,IINF,NDGSAG,NDGENG,&
         & ZEPS,LLFICP,NULOUT)  
        IF (.NOT.LLFICP) THEN
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED: LLFICP=.F.')
        ENDIF
        CALL FADIES (IREP,INULCL,IDATEF)
        IF (IDATEF(2) /= IMT1) THEN
          WRITE(NULOUT,'('' ERROR IDATEF '',11I6/'' IMT1 '',I4)')IDATEF,IMT1
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED')
        ENDIF
!*    Indices for horizontal loops
        IF (LLFICP) THEN
          INDEX=NLOENG(1)+1
        ELSE
          INDEX=1
        ENDIF

! masque dynamique-statistique (dyn=1./stat=0./+ transitions)
        CALL FACILE (IREP,INULCL,'SURF',1,'COMMASK',ZREALG,.FALSE.)

        WRITE (NULOUT,*) ' LECTURE DE COMMASK'

        ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
        WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
         & '' FIELD:'')') 'SURFCOMMASK'  
        INUM=INUM+1
        CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,SSTMSK)

        WRITE (NULOUT,*) ' ENVOIE DE COMMASK DEPUIS PROC 1'

! correction de flux du au transport
        CALL FACILE (IREP,INULCL,'SURF',1,'TRANSFLUX',ZREALG,.FALSE.)

        WRITE (NULOUT,*) ' LECTURE DE TRANSFLUX'

        ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
        WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
         & '' FIELD:'')') 'SURFTRANSFLUX'  
        INUM=INUM+1
        CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,TRAFLX)

        WRITE (NULOUT,*) ' ENVOIE DE TRANSFLUX DEPUIS PROC 1'

        CALL FAIRME (IREG,INULCL,'KEEP')

        DO JNULCL=1,99
          INQUIRE(UNIT=JNULCL,OPENED=LLOPENED)
          IF (.NOT.LLOPENED) EXIT
        ENDDO
        INULCL=JNULCL
        CALL FAITOU (IREP,INULCL,LLNOMM,'ICMCM0123INIT','OLD',&
         & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  

        WRITE (NULOUT,*) ' OUVERTURE DU FICHIER ICMCM0123INIT'

! valeurs des parametres du modele archivees a la fin du passage precedent
        DO JV=1,NVCOM
          CALL FACILE (IREP,INULCL,CLOST(JV)(1:4),1,&
           & CLOST(JV)(5:),ZREALG,.FALSE.)  

          WRITE (NULOUT,*) ' LECTURE DE ',CLOST(JV)(5:)

          ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
          WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
           & '' FIELD:'')') CLOST(JV)  
          INUM=INUM+1
          IF (JV == 1) THEN
            CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,OMLDTH)
          ELSEIF (JV == 2) THEN
            CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,GTTLIN)
          ENDIF

          WRITE (NULOUT,*) ' ENVOIE DE ',CLOST(JV)(5:),' DEPUIS PROC 1'

        ENDDO
! temperature de surface archivee a la fin du passage precedent
        CALL FACILE (IREP,INULCL,'SURF',1,'TEMPERATURE',ZREALG,.FALSE.)

        WRITE (NULOUT,*) ' LECTURE DE LA TEMPERATURE'

        ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
        WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
         & '' FIELD:'')') 'SURFTEMPERATURE'  
        INUM=INUM+1
        CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,SSTPRE)

        WRITE (NULOUT,*) ' ENVOIE DE LA TEMPERATURE DEPUIS PROC 1'

! Initialisation de la temperature
        CALL FAIRME (IREG,INULCL,'DELETE')
      ELSE
        INUM=INUM+1
        CALL DISGRID_RECV(YDGEOMETRY,1,1,SSTMSK,INUM)

        WRITE (NULOUT,*) ' RECEPTION DE SSTMSK SUR PROC ',NPROC

        DO JV=1,NVCOM
          INUM=INUM+1
          IF (JV == 1) THEN
            CALL DISGRID_RECV(YDGEOMETRY,1,1,OMLDTH,INUM)

            WRITE (NULOUT,*) ' RECEPTION DE OMLDTH SUR PROC ',NPROC

          ELSEIF (JV == 2) THEN
            CALL DISGRID_RECV(YDGEOMETRY,1,1,GTTLIN,INUM)

            WRITE (NULOUT,*) ' RECEPTION DE GTTLIN SUR PROC ',NPROC

          ENDIF
        ENDDO
        INUM=INUM+1
        CALL DISGRID_RECV(YDGEOMETRY,1,1,SSTPRE,INUM)

        WRITE (NULOUT,*) ' RECEPTION DE SSTPRE SUR PROC ',NPROC

      ENDIF

    ELSE

! sinon on applique la methode updcli (interpolation lineaire entre deux mois)
! pour les parametres du modele et on recupere la temperature calculee par updcli
! ---------------------------------------------------------------------------
! On  recupere TSURF dans GPPBUF ecrite par updcli

      DO JKGLO=1,NGPTOT,NPROMA
        IST=1
        IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        DO JROF=IST,IEND
          IROF=JROF+JKGLO-1
          SSTPRE(IROF)=SP_RR(JROF,YSP_RR%YT%MP,IBL)
        ENDDO
      ENDDO

!*    Opening Arpege files : month 1

      IF (MYPROC == 1) THEN
        WRITE(UNIT=CLNOMF,FMT='(''Const.Clim.'',I2.2)') IMT1
        DO JNULCL=1,99
          INQUIRE(UNIT=JNULCL,OPENED=LLOPENED)
          IF (.NOT.LLOPENED) EXIT
        ENDDO
        INULCL1=JNULCL
        CALL FAITOU (IREP,INULCL1,LLNOMM,CLNOMF,'OLD',&
         & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  
        ZEPS=1.E-10_JPRB
        IINF=0
        CALL CHIEN (CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,&
         & NSMAX,NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,&
         & VP00,YDVAB%VALH,YDVAB%VBH,NQUAD,IINF,NDGSAG,NDGENG,&
         & ZEPS,LLFICP,NULOUT)  
        IF (.NOT.LLFICP) THEN
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED: LLFICP=.F.')
        ENDIF
        CALL FADIES (IREP,INULCL1,IDATEF)
        IF (IDATEF(2) /= IMT1) THEN
          WRITE(NULOUT,'('' ERROR IDATEF '',11I6/'' IMT1 '',I4)')IDATEF,IMT1
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED')
        ENDIF

!*    Opening Arpege files : month 2

        WRITE(UNIT=CLNOMF,FMT='(''Const.Clim.'',I2.2)') IMT2
        DO JNULCL=1,99
          INQUIRE(UNIT=JNULCL,OPENED=LLOPENED)
          IF (.NOT.LLOPENED) EXIT
        ENDDO
        INULCL2=JNULCL
        CALL FAITOU (IREP,INULCL2,LLNOMM,CLNOMF,'OLD',&
         & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  
        CALL FADIES (IREP,INULCL2,IDATEF)
        IF (IDATEF(2) /= IMT2) THEN
          WRITE(NULOUT,'('' ERROR IDATEF '',11I6/'' IMT2 '',I4)')IDATEF,IMT2
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED')
        ENDIF

!*    Indices for horizontal loops

        IF (LLFICP) THEN
          INDEX=NLOENG(1)+1
        ELSE
          INDEX=1
        ENDIF

! indice terre mer pour limiter les modifications aux surfaces maritimes
        CALL FACILE (IREP,INULCL1,'SURF',1,'IND.TERREMER',ZREALG,.FALSE.)
        ZLSMG=ZREALG(INDEX:NGPTOTG+INDEX-1)

        DO JV=1,NVCOM
          DO JM=1,2
            IF (JM == 1) THEN
              INULCL=INULCL1
            ELSE
              INULCL=INULCL2
            ENDIF
            CALL FACILE (IREP,INULCL,CLOST(JV)(1:4),1,&
             & CLOST(JV)(5:),ZREALG,.FALSE.)  
            WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
             & '' FIELD:'')') CLOST(JV)  
            IF (JM == 1) THEN
              ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
            ELSE
              DO JROF=1,NGPTOTG
                IF (ZLSMG(JROF) < 0.5_JPRB) THEN
                  ZFIELDBUF(JROF)=ZFIELDBUF(JROF)*ZPOID1+&
                   & ZREALG(INDEX+JROF-1)*ZPOID2  
                ENDIF
              ENDDO
            ENDIF
          ENDDO ! JM=1,2
          INUM=INUM+1
          IF (JV == 1) THEN
            CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,OMLDTH)
          ELSEIF (JV == 2) THEN
            CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,GTTLIN)
          ENDIF
        ENDDO

        CALL FACILE(IREP,INULCL,'SURF',1,'COMMASK',ZREALG,.FALSE.)
        WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
         & '' FIELD:'')') 'COMMASK'  
        ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
        INUM=INUM+1
        CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,SSTMSK)
        WRITE (NULOUT,*) ' LECTURE INITIALE SSTMSK'

        CALL FACILE (IREP,INULCL,'SURF',1,'TRANSFLUX',ZREALG,.FALSE.)
        ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
        WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
         & '' FIELD:'')') 'TRANSFLUX'  
        INUM=INUM+1
        CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,TRAFLX)

        CALL FAIRME (IREG,INULCL1,'KEEP')
        CALL FAIRME (IREG,INULCL2,'KEEP')

      ELSE
! Get  data
        DO JV=1,NVCOM
          INUM=INUM+1
          IF (JV == 1) THEN
            CALL DISGRID_RECV(YDGEOMETRY,1,1,OMLDTH,INUM)
          ELSEIF (JV == 2) THEN
            CALL DISGRID_RECV(YDGEOMETRY,1,1,GTTLIN,INUM)
          ENDIF
        ENDDO
        INUM=INUM+1
        CALL DISGRID_RECV(YDGEOMETRY,1,1,SSTMSK,INUM)
      ENDIF  !MYPROC
    ENDIF  !LEXIST

  ELSEIF (NSTEP > 0) THEN

! INPUTS du modele.

! Flux cumules sur NFRCPL*TSTEP secondes

! FCSTSU, FCSTSV : Somme (stress U et V)*TSTEP
! FCCHAS         : Somme (RTHER + CHALEAU + CHALNEI + CHASEN)*TSTEP
! FCRSOS         : Somme (RSOL)*TSTEP
! FCCHLL         : Somme (CHALEAU)*TSTEP
! FCCHLN         : Somme (CHALNEI)*TSTEP
! FCCHSS         : Somme (CHASEN)*TSTEP
! OMLDTH         : Profondeur de la couche de melange oceanique
! GTTLIN         : Gradient thermique de la thermocline
! SSTPRE         : Temperature de surface

! OUTPUTS.

! IF ( LOMLDTH ) THEN
! OMLDTH         : Profondeur de la couche de melange oceanique
! ENDIF
! GTTLIN         : Gradient thermique de la thermocline
! SSTPRE         : Temperature de surface prevue par le modele dynamique 1D

    IF ( .NOT. LOMLDTH ) THEN

! L'evolution de la profondeur de melange n'est pas calculee par le modele
! de couche de melange

! On applique l'interpolation lineaire entre deux mois
! pour la profondeur de melange
! -----------------------------

!*    Opening Arpege files : month 1

      IF (MYPROC == 1) THEN
        WRITE(UNIT=CLNOMF,FMT='(''Const.Clim.'',I2.2)') IMT1
        DO JNULCL=1,99
          INQUIRE(UNIT=JNULCL,OPENED=LLOPENED)
          IF (.NOT.LLOPENED) EXIT
        ENDDO
        INULCL1=JNULCL
        CALL FAITOU (IREP,INULCL1,LLNOMM,CLNOMF,'OLD',&
         & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  
        ZEPS=1.E-10_JPRB
        IINF=0
        CALL CHIEN (CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,&
         & NSMAX,NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,&
         & VP00,YDVAB%VALH,YDVAB%VBH,NQUAD,IINF,NDGSAG,NDGENG,&
         & ZEPS,LLFICP,NULOUT)  
        IF (.NOT.LLFICP) THEN
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED: LLFICP=.F.')
        ENDIF
        CALL FADIES (IREP,INULCL1,IDATEF)
        IF (IDATEF(2) /= IMT1) THEN
          WRITE(NULOUT,'('' ERROR IDATEF '',11I6/'' IMT1 '',I4)')IDATEF,IMT1
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED')
        ENDIF

!*    Opening Arpege files : month 2

        WRITE(UNIT=CLNOMF,FMT='(''Const.Clim.'',I2.2)') IMT2
        DO JNULCL=1,99
          INQUIRE(UNIT=JNULCL,OPENED=LLOPENED)
          IF (.NOT.LLOPENED) EXIT
        ENDDO
        INULCL2=JNULCL
        CALL FAITOU (IREP,INULCL2,LLNOMM,CLNOMF,'OLD',&
         & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  
        CALL FADIES (IREP,INULCL2,IDATEF)
        IF (IDATEF(2) /= IMT2) THEN
          WRITE(NULOUT,'('' ERROR IDATEF '',11I6/'' IMT2 '',I4)')IDATEF,IMT2
          CALL ABOR1 ('UPDCPL: ABOR1 CALLED')
        ENDIF

!*    Indices for horizontal loops

        IF (LLFICP) THEN
          INDEX=NLOENG(1)+1
        ELSE
          INDEX=1
        ENDIF

        DO JM=1,2
          IF (JM == 1) THEN
            INULCL=INULCL1
          ELSE
            INULCL=INULCL2
          ENDIF
          CALL FACILE (IREP,INULCL,CLOST(1)(1:4),1,&
           & CLOST(1)(5:),ZREALG,.FALSE.)  
          WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
           & '' FIELD:'')') CLOST(JV)  
          IF (JM == 1) THEN
            ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
          ELSE
            DO JROF=1,NGPTOTG
              IF (ZLSMG(JROF) < 0.5_JPRB) THEN
                ZFIELDBUF(JROF)=ZFIELDBUF(JROF)*ZPOID1+&
                 & ZREALG(INDEX+JROF-1)*ZPOID2  
              ENDIF
            ENDDO
          ENDIF
        ENDDO ! JM=1,2
        INUM=INUM+1
        CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,INUM,OMLDTH)

        CALL FAIRME (IREG,INULCL1,'KEEP')
        CALL FAIRME (IREG,INULCL2,'KEEP')

      ELSE
! Get  data
        INUM=INUM+1
        CALL DISGRID_RECV(YDGEOMETRY,1,1,OMLDTH,INUM)
      ENDIF  !MYPROC

    ENDIF

    WRITE (NULOUT,*) ' NSTEP=',NSTEP,' DANS UPDCPL.F AVANT L''APPEL AU SLAB'

    CALL SLAB(YDGEOMETRY%YRGEM,YDCOM,ZLSM)

    WRITE (NULOUT,*) ' DANS UPDCPL.F APRET L''APPEL AU SLAB'

! Remise a zero des cumuls de flux

    FCSTSU=0.0_JPRB
    FCSTSV=0.0_JPRB
    FCCHAS=0.0_JPRB
    FCRSOS=0.0_JPRB
    FCCHLL=0.0_JPRB
    FCCHLN=0.0_JPRB
    FCCHSS=0.0_JPRB

  ENDIF !NSTEP == 0

  DO J=1,NGPTOT

! MIXAGE STATISTIQUE-DYNAMIQUE

! WARNING !!!  LE MODELE DYNAMIQUE NE GERE PAS LA BANQUISE.

    IF ( ZFIELD(J,IZST) < TMERGL ) THEN
      SSTPRE(J)=ZFIELD(J,IZST)
      ZFIELD(J,IZOT)=SSTPRE(J)
    ELSE
      SSTPRE(J)=MAX(TMERGL,SSTPRE(J))
      ZFIELD(J,IZOT)=SSTMSK(J)*SSTPRE(J)+(1.0_JPRB-SSTMSK(J))*ZFIELD(J,IZST)
      SSTPRE(J)=ZFIELD(J,IZOT)
    ENDIF

  ENDDO

ENDIF !LMCC03-LMCC05

!*    3.1   Reads the buffers

! This section is now in agreement with updsst.F90 and external updcli
ZWXMER=1000._JPRB*RD2MER
ZWXGLA=1000._JPRB*RD2GLA
ZWXSUR=1000._JPRB*RD1
ZRANO=0.70_JPRB
ZRDNO=0.30_JPRB
DO JKGLO=1,NGPTOT,NPROMA

  IST=1
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1

  DO JROF=IST,IEND
    IROF=JROF+JKGLO-1

    SD_VF(JROF,YSD_VF%YUCUR%MP,IBL)=0._JPRB
    SD_VF(JROF,YSD_VF%YVCUR%MP,IBL)=0._JPRB

    LLMEGL=(NINT(SD_VV(JROF,YSD_VV%YIVEG%MP,IBL)) == NTVMER.OR.&
       & (NINT(10*SD_VV(JROF,YSD_VV%YIVEG%MP,IBL))) == (10*NTVGLA+1))  
    IF(LLMEGL)THEN

      IF (LMCC03) THEN

        IF (ZFIELD(IROF,IZOI) >= 0.5_JPRB) THEN

!     - - - - -
!     Sea Ice 
!     - - - - -
          IF (NINT(SD_VV(JROF,YSD_VV%YIVEG%MP,IBL)) == NTVMER) THEN
!       New formed Sea Ice
            SP_RR(JROF,YSP_RR%YT%MP,IBL)=ZFIELD(IROF,IZOT)
            DO JCSS=1,YSP_SBD%NLEVS
              SP_SB(JROF,JCSS,YSP_SB%YT%MP,IBL)=TMERGL
            ENDDO
            SP_SG(JROF,1,YSP_SG%YA%MP,IBL)=ZRANO
            SP_SG(JROF,1:YSP_SGD%NLEVS,YSP_SG%YR%MP,IBL)=ZRDNO
            IF (LMCC02) THEN
              SD_VF(JROF,YSD_VF%YLSM%MP,IBL)=1.
            ELSE
              SD_VF(JROF,YSD_VF%YLSM%MP,IBL)=0.
            ENDIF
            SD_VF(JROF,YSD_VF%YZ0F%MP,IBL)=RZ0GLA*RG
            IF(JPFLDO2A >= 3) THEN
              SD_VF(JROF,YSD_VF%YALBF%MP,IBL)=ZFIELD(IROF,IZOA)
            ELSE
              SD_VF(JROF,YSD_VF%YALBF%MP,IBL)=ALBGLA
            ENDIF
            SD_VF(JROF,YSD_VF%YEMISF%MP,IBL)=EMMGLA
            SP_RR(JROF,YSP_RR%YW%MP,IBL)=0.
            SP_RR(JROF,YSP_RR%YIC%MP,IBL)=ZWXSUR
            SP_SB(JROF,1,YSP_SB%YQ%MP,IBL)=0.
            SP_SB(JROF,1,YSP_SB%YTL%MP,IBL)=ZWXGLA
            SP_SG(JROF,1:YSP_SGD%NLEVS,YSP_SG%YF%MP,IBL)=0.
            SD_VV(JROF,YSD_VV%YIVEG%MP,IBL)=REAL(NTVGLA)+.1
            SD_VV(JROF,YSD_VV%YD2%MP,IBL)=RD2GLA
            IF (YSD_VVD%NUMFLDS >= 8) THEN
              SD_VV(JROF,YSD_VV%YZ0H%MP,IBL)=RZHZ0G*RZ0GLA*RG
            ENDIF
          ELSE
!       Old Sea Ice
            IF (LMCC02) THEN
              SD_VF(JROF,YSD_VF%YLSM%MP,IBL)=1.
            ELSE
              SP_RR(JROF,YSP_RR%YT%MP,IBL)=ZFIELD(IROF,IZOT)
              DO JCSS=1,YSP_SBD%NLEVS
                SP_SB(JROF,JCSS,YSP_SB%YT%MP,IBL)=SP_RR(JROF,YSP_RR%YT%MP,IBL)
              ENDDO
              SD_VF(JROF,YSD_VF%YLSM%MP,IBL)=0.
            ENDIF
            SD_VF(JROF,YSD_VF%YZ0F%MP,IBL)=RZ0GLA*RG
            IF(JPFLDO2A >= 3) THEN
              SD_VF(JROF,YSD_VF%YALBF%MP,IBL)=ZFIELD(IROF,IZOA)
            ELSE
              SD_VF(JROF,YSD_VF%YALBF%MP,IBL)=ALBGLA
            ENDIF
            SD_VF(JROF,YSD_VF%YEMISF%MP,IBL)=EMMGLA
            SP_RR(JROF,YSP_RR%YW%MP,IBL)=0.
            SP_RR(JROF,YSP_RR%YIC%MP,IBL)=ZWXSUR
            SP_SB(JROF,1,YSP_SB%YQ%MP,IBL)=0.
            SP_SB(JROF,1,YSP_SB%YTL%MP,IBL)=ZWXGLA
            SD_VV(JROF,YSD_VV%YIVEG%MP,IBL)=REAL(NTVGLA)+.1
            SD_VV(JROF,YSD_VV%YD2%MP,IBL)=RD2GLA
            IF (YSD_VVD%NUMFLDS >= 8) THEN
              SD_VV(JROF,YSD_VV%YZ0H%MP,IBL)=RZHZ0G*RZ0GLA*RG
            ENDIF
          ENDIF
        ELSE
!     - - - - -
!     Open Sea
!     - - - - -

          SD_VF(JROF,YSD_VF%YLSM%MP,IBL)=0.
  ! should be updated by Charnock formula
          SD_VF(JROF,YSD_VF%YZ0F%MP,IBL)=RZ0MER*RG
          IF(JPFLDO2A >= 3) THEN
            SD_VF(JROF,YSD_VF%YALBF%MP,IBL)=ZFIELD(IROF,IZOA)
          ELSE
            SD_VF(JROF,YSD_VF%YALBF%MP,IBL)=ALBMER
          ENDIF
          SD_VF(JROF,YSD_VF%YEMISF%MP,IBL)=EMMMER
          SD_VF(JROF,YSD_VF%YUCUR%MP,IBL)=ZFIELD(IROF,IZOU)
          SD_VF(JROF,YSD_VF%YVCUR%MP,IBL)=ZFIELD(IROF,IZOV)
          SP_RR(JROF,YSP_RR%YT%MP,IBL)=ZFIELD(IROF,IZOT)
          DO JCSS=1,YSP_SBD%NLEVS
            SP_SB(JROF,JCSS,YSP_SB%YT%MP,IBL)=SP_RR(JROF,YSP_RR%YT%MP,IBL)
          ENDDO
          SP_SG(JROF,1,YSP_SG%YA%MP,IBL)=ZRANO
          SP_SG(JROF,1:YSP_SGD%NLEVS,YSP_SG%YR%MP,IBL)=ZRDNO
          SP_RR(JROF,YSP_RR%YW%MP,IBL)=ZWXSUR
          SP_RR(JROF,YSP_RR%YIC%MP,IBL)=0.
          SP_SB(JROF,1,YSP_SB%YQ%MP,IBL)=ZWXMER
          SP_SB(JROF,1,YSP_SB%YTL%MP,IBL)=0.
          SP_SG(JROF,1:YSP_SGD%NLEVS,YSP_SG%YF%MP,IBL)=0.
          SD_VV(JROF,YSD_VV%YIVEG%MP,IBL)=REAL(NTVMER)
          SD_VV(JROF,YSD_VV%YD2%MP,IBL)=RD2MER
          IF (YSD_VVD%NUMFLDS >= 8) THEN
            SD_VV(JROF,YSD_VV%YZ0H%MP,IBL)=RZHZ0M*RZ0MER*RG
          ENDIF
        ENDIF

      ELSEIF (LMCC05) THEN

        SP_RR(JROF,YSP_RR%YT%MP,IBL)=ZFIELD(IROF,IZOT)
        DO JCSS=1,YSP_SBD%NLEVS
          SP_SB(JROF,JCSS,YSP_SB%YT%MP,IBL)=SP_RR(JROF,YSP_RR%YT%MP,IBL)
        ENDDO

      ENDIF  !LMCC03/LMCC05

    ENDIF ! no land

  ENDDO !END OF JROF=IST,IEND

ENDDO  !END OF JKGLO=1,NGPTOT,NPROMA

!-----------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDCPL',1,ZHOOK_HANDLE)
END SUBROUTINE UPDCPL

