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

SUBROUTINE CUCONVCA(YDGEOMETRY,YDGMV,YDMODEL,KSTEP)

!**** *CUCONVCA*  - Control Cellular Automaton for convective forcing

!     Purpose.
!     --------
!           Cellular automaton for convection control routine

!**   Interface.
!     ----------
!        *CALL* *CUCONVCA*

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        Module YOE_CUCONVCA

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------

!     Author.
!     -------
!        P. Bechtold Feb. 2009 copied on J Berner previous routines surand* 

!     Modifications.
!     --------------
!        M. Steinheimer Mar-Aug 2009   changed to probabilistic CA,
!           added optional wind dependency, bilinear interpolation, ... 
!        L. Bengtsson-Sedlar Mar-Dec 2010 Added deterministic rules according to GOL
!                                         Adapted to LAM version...
!        F. Vana  Jan-2011 : optimization - reordered arrays, exchanged loops,
!                            NPROMA slicing, more OpenMP... 
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type SL_STRUCT
!        L. Bengtsson (Aug 2014): Spin-up correction if LCA_GLOBAL, correction if GOL=TRUE
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE YOE_CUCONVCA       , ONLY : NIJH, INITIALIZE_CELLS, WRITE_FIELD, UPDCELAUT_RGG
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE RANDOM_NUMBERS_MIX , ONLY : UNIFORM_DISTRIBUTION
USE YOMMP0             , ONLY : MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP ! model time step
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM)              :: IT, ISPINUP
INTEGER(KIND=JPIM)              :: JROF, JK

INTEGER(KIND=JPIM), ALLOCATABLE :: ILIVES(:),IINI(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IDX(:), IDY(:)

REAL(KIND=JPRB), ALLOCATABLE    :: ZU(:,:), ZV(:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZELAMG(:,:), ZELATG(:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZELAM(:,:), ZELAT(:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZRAND1DLS(:),ZRAND1DSS(:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZRAND1DLSG(:),ZRAND1DSSG(:)

REAL(KIND=JPHOOK)                 :: ZHOOK_HANDLE

INTEGER(KIND=JPIM)              :: IMAINPROC
CHARACTER(LEN=25)               :: CLFNAME, CLFNAME3, CLFNAME4

INTEGER(KIND=JPIM)              :: JSS, JLS, JLSG
INTEGER(KIND=JPIM)              :: JI, ICEND, IBL, IOFF, IGP
INTEGER(KIND=JPIM)              :: INIJH2

!stuff for SL-halo
LOGICAL                         :: LLINC
INTEGER(KIND=JPIM)              :: IDUMARR(2)
INTEGER(KIND=JPIM)              :: IFIXSFLD(2)
INTEGER(KIND=JPIM)              :: IFLDSLB1
REAL(KIND=JPRB), ALLOCATABLE    :: ZPB1(:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZFERT(:,:),ZCELL(:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZCUCONVCAG(:,:),ZCUCONVCA(:,:)

!     ------------------------------------------------------------------

#include "gathergpf.intfb.h"
#include "slcomm.intfb.h"

!     ------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('CUCONVCA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDECUMF=>YDMODEL%YRML_PHY_EC%YRECUMF, &
 & YDECUCONVCA=>YDMODEL%YRML_PHY_EC%YRECUCONVCA, &
 & YDML_PHY_STOCH=>YDMODEL%YRML_PHY_STOCH, YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & CA_PROB=>YDECUCONVCA%CA_PROB, CA_WIND=>YDECUCONVCA%CA_WIND, &
 & LCA_ADVECT=>YDECUCONVCA%LCA_ADVECT, LCA_EXTRACT=>YDECUCONVCA%LCA_EXTRACT, &
 & LCA_GLOBAL=>YDECUCONVCA%LCA_GLOBAL, LCA_RANTROP=>YDECUCONVCA%LCA_RANTROP, &
 & LCA_TEST=>YDECUCONVCA%LCA_TEST, NCELLCU=>YDECUCONVCA%NCELLCU, &
 & NDXUNREAL=>YDECUCONVCA%NDXUNREAL, NDYUNREAL=>YDECUCONVCA%NDYUNREAL, &
 & NFERTCU=>YDECUCONVCA%NFERTCU, NFRCASEED=>YDECUCONVCA%NFRCASEED, &
 & NLIVES=>YDECUCONVCA%NLIVES, NSPINUP=>YDECUCONVCA%NSPINUP, &
 & NTESTGP=>YDECUCONVCA%NTESTGP, NTESTPROC=>YDECUCONVCA%NTESTPROC, &
 & RCA_SEEDPROB=>YDECUCONVCA%RCA_SEEDPROB, RCUCONVCA=>YDECUCONVCA%RCUCONVCA, &
 & RNLCONVCA=>YDECUCONVCA%RNLCONVCA, RWASALIVE=>YDECUCONVCA%RWASALIVE, &
 & RWGHTCU=>YDECUCONVCA%RWGHTCU, &
 & YD_RANDOM_STREAM_CA=>YDECUCONVCA%YD_RANDOM_STREAM_CA, &
 & NJKT4=>YDECUMF%NJKT4, NJKT5=>YDECUMF%NJKT5, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, &
 & GMV=>YDGMV%GMV, YPH9=>YDGMV%YPH9, &
 & NGLOBALINDEX=>YDMP%NGLOBALINDEX, &
 & LSTOPH_CASBS=>YDSTOPH%LSTOPH_CASBS, RSTOPHCA=>YDSTOPH%RSTOPHCA, YDSL=>YDMODEL%YRML_DYN%YRSL)
!     ------------------------------------------------------------------

!  WRITE(NULOUT,'("xxxx in CUCONVCA at step ",I4)') kstep
!  CALL FLUSH(NULOUT)
!     ------------------------------------------------------
!*       1.0  Setup constants and allocate arrays
!     ------------------------------------------------------

INIJH2=NIJH*NIJH

ALLOCATE(ILIVES(NGPTOT),IINI(NGPTOT))

ALLOCATE(IDX(NGPTOT),IDY(NGPTOT))

ALLOCATE(ZRAND1DSS(INIJH2*NGPTOT),ZRAND1DLS(NGPTOT))
ALLOCATE(ZRAND1DSSG(INIJH2*NGPTOTG),ZRAND1DLSG(NGPTOTG))

ILIVES=0
IINI=0

!stuff for SL halo
LLINC=.FALSE.
IFIXSFLD(:)=0
IFLDSLB1=2*INIJH2
ALLOCATE(ZPB1(YDSL%NASLB1,IFLDSLB1))
ALLOCATE(ZFERT(YDSL%NASLB1,INIJH2))
ALLOCATE(ZCELL(YDSL%NASLB1,INIJH2))

!     ------------------------------------------------------
!*       2.0  prepare data
!     ------------------------------------------------------

!*       2.1  get U and V from GMV array to local array
!     ------------------------------------------------------

ALLOCATE(ZU(NGPTOT,1), ZV(NGPTOT,1))

IF( CA_WIND == "REAL" ) THEN
  ! calculate mean wind of 850hPa (model level NJKT4) and 500hPa (ml=NJKT5)
  ! and reshape wind arrays from GMV(NPROMA,NFLEVG,NDIMGMV,NGPBLKS) to
  ! ZU(NGPTOT,1)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,ICEND,IBL,IOFF,JI,IGP)
  DO JK=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JK+1)
    IBL=(JK-1)/NPROMA+1
    IOFF=JK
    DO JI=1,ICEND
      IGP=IOFF+JI-1
      ZU(IGP,1)=  0.5*(GMV(JI,NJKT4,YPH9%MU,IBL) + GMV(JI,NJKT5,YPH9%MU,IBL) )
      ZV(IGP,1)=  0.5*(GMV(JI,NJKT4,YPH9%MV,IBL) + GMV(JI,NJKT5,YPH9%MV,IBL) )

      IDX(IGP)=NINT(SIGN(1._JPRB,ZU(IGP,1)))
      IDY(IGP)=NINT(SIGN(1._JPRB,ZV(IGP,1)))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

ELSE

  !option of idealized wind for testing
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JROF)
  DO JK=1,NGPTOT,NPROMA
    DO JROF=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
      IDX(JROF)= NDXUNREAL
      IDY(JROF)= NDYUNREAL
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!     ------------------------------------------------------
!*       3.0  Evolve the CA
!     ------------------------------------------------------

IF (LCA_GLOBAL) THEN
  IF (KSTEP==0) THEN
    ISPINUP=NSPINUP
  ELSE
    ISPINUP=1
  ENDIF
ELSE
  ISPINUP=1           !spin-up currently not used for convective CA
  IF (LCA_RANTROP .AND. KSTEP==1) ISPINUP=NSPINUP
ENDIF


  IF (LCA_GLOBAL) THEN
    IF((MOD(KSTEP,NFRCASEED) == 0 .AND. KSTEP > 0) .OR. KSTEP == 0) THEN
      CALL UNIFORM_DISTRIBUTION (ZRAND1DLSG,YD_RANDOM_STREAM_CA)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JROF)
      DO JK=1,NGPTOT,NPROMA
        DO JROF=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
          IF (ZRAND1DLSG(NGLOBALINDEX(JROF)) >= RCA_SEEDPROB ) THEN
            IINI(JROF)= 1
          ELSE
            IINI(JROF)=0
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JROF)
    DO JK=1,NGPTOT,NPROMA
      DO JROF=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
        ILIVES(JROF)=NLIVES
      ENDDO
    ENDDO
!$OMP END PARALLEL DO  
  ELSE
    IF (LCA_TEST) THEN
      !these 3 lines can be used to initialize only a single
      !gaussian gridpoint for testing
      RCUCONVCA(:)=0
      IF (MYPROC == NTESTPROC ) RCUCONVCA(NTESTGP)=1
      RNLCONVCA(:)=NLIVES
    ENDIF
    IF (LCA_RANTROP) THEN
      CALL UNIFORM_DISTRIBUTION (ZRAND1DLSG,YD_RANDOM_STREAM_CA)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JROF)
      DO JK=1,NGPTOT,NPROMA
        DO JROF=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
          ILIVES(JROF)= INT(RNLCONVCA(JROF))
          IF (ZRAND1DLSG(NGLOBALINDEX(JROF)) >= RCA_SEEDPROB .AND.&
           &  RCUCONVCA(JROF)==1._JPRB) THEN
            IINI(JROF)=1
          ELSE
            IINI(JROF)=0
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ELSE
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JROF)
      DO JK=1,NGPTOT,NPROMA
        DO JROF=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
          ILIVES(JROF)= INT(RNLCONVCA(JROF))
          IINI(JROF)= INT(RCUCONVCA(JROF))
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF
  ENDIF

!     ------------------------------------------------------
!*       3.1  Initialise the cellular automaton from convection
!     ------------------------------------------------------

  IF (LCA_GLOBAL) THEN
    IF(KSTEP==0) CALL INITIALIZE_CELLS(YDMODEL%YRML_PHY_EC%YRECUCONVCA,YDGEOMETRY,IINI,ILIVES,NCELLCU,NFERTCU)
  ELSE
    !Initialize CA at first step (at step 0 physics field not yet filled)
    IF(KSTEP==1) THEN
      CALL INITIALIZE_CELLS(YDMODEL%YRML_PHY_EC%YRECUCONVCA,YDGEOMETRY,IINI,ILIVES,NCELLCU,NFERTCU)
    ENDIF
  ENDIF
  !Initialize CA at each NFRCASEED'th step
  IF(MOD(KSTEP,NFRCASEED) == 0 .AND. KSTEP > 0) THEN
    CALL INITIALIZE_CELLS(YDMODEL%YRML_PHY_EC%YRECUCONVCA,YDGEOMETRY,IINI,ILIVES,NCELLCU,NFERTCU)
  ENDIF

!  WRITE(NULOUT,'("xxxx CA initialized")')
!  CALL FLUSH(NULOUT)


!*       3.2  get halo for CA
!     ------------------------------------------------------

DO IT=1,ISPINUP

  ZPB1=0
  ZFERT=0
  ZCELL=0

!fill buffer
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JLS,JSS)
  DO JK=1,NGPTOT,NPROMA
    DO JSS=1,INIJH2
      DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
        ZPB1(YDSL%NSLCORE(JLS),JSS)=NCELLCU(JLS,JSS)
        ZPB1(YDSL%NSLCORE(JLS),JSS+INIJH2)=NFERTCU(JLS,JSS)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

!  WRITE(NULOUT,'("xxxx SL buffer set up")')
!  CALL FLUSH(NULOUT)

!get halo
  CALL SLCOMM(YDSL,IDUMARR,IFLDSLB1,LLINC,0,ZPB1)

!  WRITE(NULOUT,'("xxxx SL-Halo retrieved")')
!  CALL FLUSH(NULOUT)

!extract from buffer
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JLS,JSS)
  DO JK=1,YDSL%NASLB1,NPROMA
    DO JSS=1,INIJH2
      DO JLS=JK,JK+MIN(NPROMA,YDSL%NASLB1-JK+1)-1
        ZCELL(JLS,JSS)=ZPB1(JLS,JSS)
        ZFERT(JLS,JSS)=ZPB1(JLS,JSS+INIJH2)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

!  WRITE(NULOUT,'("xxxx buffer extracted")')
!  CALL FLUSH(NULOUT)

!     ------------------------------------------------------
!*       3.3  Update cellular automaton 
!             and transform fields back to reduced Gaussian grid
!     ------------------------------------------------------

  ! create random numbers
  IF (CA_PROB=='GOL')THEN
    ZRAND1DSSG(:)=1._JPRB
  ELSE
    CALL UNIFORM_DISTRIBUTION (ZRAND1DSSG,YD_RANDOM_STREAM_CA)
  ENDIF

  !extract local random numbers from global field
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JLS,JSS,JLSG)
  DO JK=1,NGPTOT,NPROMA
    DO JSS=1,INIJH2
      DO JLS=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
        JLSG=NGLOBALINDEX(JLS)
        ZRAND1DSS((JLS-1)*INIJH2+JSS)=ZRAND1DSSG((JLSG-1)*INIJH2+JSS)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  CALL UPDCELAUT_RGG(YDGEOMETRY,YDECUCONVCA,YDECUMF,YDSL,ILIVES,IDX,IDY,NFERTCU,NCELLCU,ZFERT,ZCELL,&
   & RWGHTCU,ZRAND1DSS,LCA_ADVECT)

ENDDO  !end spin-up loop

!  WRITE(NULOUT,'("xxxx CA updated")')
!  CALL FLUSH(NULOUT)

!set WASALIVE field to one at gridpoints where CA was active
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JK,JROF)
DO JK=1,NGPTOT,NPROMA
  DO JROF=JK,JK+MIN(NPROMA,NGPTOT-JK+1)-1
    RCUCONVCA(JROF)=RWGHTCU(JROF)
    IF (LSTOPH_CASBS) RSTOPHCA(JROF)=RWGHTCU(JROF)
    IF (RCUCONVCA(JROF) > 0) RWASALIVE(JROF)=1
    RNLCONVCA(JROF)=RWASALIVE(JROF)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

IF (LCA_EXTRACT) THEN
  ALLOCATE(ZCUCONVCA(NGPTOT,1),ZCUCONVCAG(NGPTOTG,1))
  ALLOCATE(ZELAMG(NGPTOTG,1),ZELATG(NGPTOTG,1))
  ALLOCATE(ZELAM(NGPTOT,1),ZELAT(NGPTOT,1))

  IMAINPROC=1

  ZCUCONVCA(:,1)=RCUCONVCA

!get whole fields on 1 processor

  CALL GATHERGPF(YDGEOMETRY,ZCUCONVCA(:,1),ZCUCONVCAG,1,IMAINPROC)
  CALL GATHERGPF(YDGEOMETRY,ZELAM(:,1),ZELAMG,1,IMAINPROC)
  CALL GATHERGPF(YDGEOMETRY,ZELAT(:,1),ZELATG,1,IMAINPROC)

!binary output of wind for debugging
  IF (IMAINPROC == MYPROC) THEN
    WRITE(CLFNAME,'(A,I4.4,A)') 'CA',KSTEP,'.dat'
    CALL WRITE_FIELD(CLFNAME, REAL(ZCUCONVCAG(:,1),JPRB), NGPTOTG, 1)
    IF (KSTEP == 0) THEN
      WRITE(CLFNAME3,'(A,I4.4,A)') 'LO',KSTEP,'.dat'
      WRITE(CLFNAME4,'(A,I4.4,A)') 'LA',KSTEP,'.dat'
      CALL WRITE_FIELD(CLFNAME3, REAL(ZELAMG(:,1),JPRB), NGPTOTG, 1)
      CALL WRITE_FIELD(CLFNAME4, REAL(ZELATG(:,1),JPRB), NGPTOTG, 1)
    ENDIF
  ENDIF

  DEALLOCATE(ZELAMG,ZELATG,ZELAM,ZELAT)
  DEALLOCATE(ZCUCONVCA,ZCUCONVCAG)

ENDIF

!     ------------------------------------------------------
!*       4.0 Cleanup
!     ------------------------------------------------------

DEALLOCATE(ILIVES)
DEALLOCATE(IINI)

DEALLOCATE(ZU, ZV)
DEALLOCATE(IDX, IDY)
DEALLOCATE(ZRAND1DSS,ZRAND1DLS)
DEALLOCATE(ZRAND1DSSG,ZRAND1DLSG)

DEALLOCATE(ZPB1)
DEALLOCATE(ZFERT)
DEALLOCATE(ZCELL)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUCONVCA',1,ZHOOK_HANDLE)
END SUBROUTINE CUCONVCA
