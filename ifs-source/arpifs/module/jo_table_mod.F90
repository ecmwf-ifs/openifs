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

MODULE JO_TABLE_MOD

USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
USE PARDIMO            , ONLY : JPNOTP
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMANCS            , ONLY : RMDI  
USE VARNO_MODULE       , ONLY : VARNO
USE NVAR_CLASS         , ONLY : CLASS_NVAR
USE YOMCOCTP           , ONLY : NSATEM, NSATOB, NLIMB, NALLSKY, NSCATT, MSQBYCTP
USE YOMCOSJO           , ONLY : JPST
USE OBSOP_SETS         , ONLY : TYPE_SET_INFO
USE YOMCT0             , ONLY : LSCREEN, NCONF, LGUESS, L_OOPS
USE YOMMP0             , ONLY : MYPROC, NPROC
USE YOMCMA             , ONLY : NMXOTP
USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA,GET_NSIM4D,L_OBS_IN_FC
USE YOMNMEV            , ONLY : NODEV1
USE YOMLUN             , ONLY : NULERR, NULOUT, RESERVE_LUN, FREE_LUN
USE YOMDIAGVAR         , ONLY : DIAG_4DVAR
USE YOMCHEV            , ONLY : CH1EVENT
USE MPL_MODULE


IMPLICIT NONE

!       JOT JO table, derived type
TYPE JOTSUB_T
  CHARACTER(LEN=32)               :: SNAME         ! name of sub-table
  REAL(KIND=JPRB)   , ALLOCATABLE :: COST(:,:)     ! cost function value
  INTEGER(KIND=JPIM), ALLOCATABLE :: COUNT(:,:)    ! data count (where fg_depar /= RMDI)
  INTEGER(KIND=JPIM), ALLOCATABLE :: COUNTALL(:,:) ! data count (ALL data)
  REAL(KIND=JPRB)   , ALLOCATABLE :: OBIMPACT(:,:) ! obs impact
  REAL(KIND=JPRB)   , ALLOCATABLE :: OBSERR(:)     ! Sum of (obs error)**2
  REAL(KIND=JPRB)   , ALLOCATABLE :: BGERR(:)      ! Sum of (bg error)**2
  INTEGER(KIND=JPIM)              :: CODETYPE      ! obs codetype from ODB
  INTEGER(KIND=JPIM), ALLOCATABLE :: STATUS(:,:)   ! datum status summary 
  INTEGER(KIND=JPIM), ALLOCATABLE :: EVENT1(:,:)   ! datum event1 summary 

CONTAINS
  PROCEDURE :: CREATE
END TYPE JOTSUB_T

TYPE JOT_T
  CHARACTER(LEN=32)           :: NAME              ! name of table entry
  INTEGER(KIND=JPIM)          :: SIZE              ! size of sub-table
  TYPE(JOTSUB_T), ALLOCATABLE :: JOST(:)           ! Jo sub-table
END TYPE JOT_T

TYPE JO_TABLE
  TYPE(JOT_T)                 :: JOT(JPNOTP)
  INTEGER(KIND=JPIM)          :: NPRTBINS          ! Number of bins for Jo prints
  INTEGER(KIND=JPIM)          :: NJOTTSZ           ! Total size of cost and count
  INTEGER(KIND=JPIM)          :: NJOTISZ           ! Size of table per item
  TYPE(CLASS_NVAR)            :: NVAR              ! Mapping of VARNOs into the Jo table

CONTAINS
  PROCEDURE, PUBLIC       :: SETUP
  PROCEDURE, PUBLIC       :: RESET
  PROCEDURE, PUBLIC       :: UPDATE
  PROCEDURE, PUBLIC       :: GATHER
  PROCEDURE, PUBLIC       :: PRINT_JO
  PROCEDURE, PUBLIC       :: PRINT_SCREEN_STATS
END TYPE JO_TABLE


!     ------------------------------------------------------------------
CONTAINS

SUBROUTINE SETUP(YDJOT)
   CLASS(JO_TABLE), INTENT(INOUT)   :: YDJOT
   WRITE(0,*) "JO_TABLE % SETUP : NOT IMPLEMENTED YET"
   ! A place holder for allocations and category description setups done under suobs in suamv, suscat, surad etc
   ! And the NVARS are set up in suobs.F90 as well.

END SUBROUTINE SETUP

SUBROUTINE RESET(YDJOT)

   CLASS(JO_TABLE), INTENT(INOUT)   :: YDJOT
   INTEGER(KIND=JPIM)               :: I
   INTEGER(KIND=JPIM)               :: J
   REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

   !     ------------------------------------------------------------------
   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:RESET',0,ZHOOK_HANDLE)
   !     ------------------------------------------------------------------

   DO I = 1,SIZE(YDJOT%JOT,DIM=1)
      DO J = 1,YDJOT%JOT(I)%SIZE
         YDJOT%JOT(I)%JOST(J)%COST(:,:)     = 0.0_JPRB
         YDJOT%JOT(I)%JOST(J)%OBSERR(:)     = 0.0_JPRB
         YDJOT%JOT(I)%JOST(J)%BGERR(:)      = 0.0_JPRB
         YDJOT%JOT(I)%JOST(J)%OBIMPACT(:,:) = 0.0_JPRB
         YDJOT%JOT(I)%JOST(J)%COUNT(:,:)    = 0
         YDJOT%JOT(I)%JOST(J)%COUNTALL(:,:) = 0
         YDJOT%JOT(I)%JOST(J)%STATUS(:,:)   = 0
         YDJOT%JOT(I)%JOST(J)%EVENT1(:,:)   = 0
      ENDDO
   ENDDO

   !     ------------------------------------------------------------------
   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:RESET',1,ZHOOK_HANDLE)
   !     ------------------------------------------------------------------

END SUBROUTINE RESET

SUBROUTINE PRINT_JO(YDJOT,KNSMAX)
   
   !**** *PRINT_JO*- PRINT DIFFERENT TERMS OF THE COST FUNCTION FROM THE
   !****            DIAGNOSTIC JO-TABLE
   
   !     Purpose.
   !     -------
   !           PRINT DIFFERENT TERMS OF THE COST FUNCTION IN ORDER TO
   !     FIND OUT THE RESPECTIVE CONTRIBUTIONS OF EACH OBS TYPE, EACH
   
   !     SUB OBS TYPE, EACH OBSERVED VARIABLE....
   !**   Interface.
   !**   ---------
   !**   *CALL* * PRTJO*
   !**         The common YOMCOSJO is an implicit INPUT of SUDIMO  AS IT
   !**         CONTAINS THE COST FUNCTION TABLE JOT.
   
   !     AUTHOR.
   !     -------
   !      Author: Jean Pailleux
   !      Date  : 89-12-29
   
   !     Modifications
   !     -------------
   !      M.Hamrud      01-Oct-2003 CY28 Cleaning
   !      Y.Tremolet    29-May-2008 Print Jo table by bins
   !      A. Geer          22-Nov-2010 Event summary table 
   !      N. Bormann     1-Aug-2016 CVarBC
   !      P. Lean       24-Mar-2107 Moved inside JO_TABLE object for OOPS
   !      S. Massart    19-Feb-2019 Parameter optimisation
   !      B. Ingleby    24-Apr-2020 Tidy call to map_varno_to_nvar
   !     ------------------------------------------------------------------
   
   IMPLICIT NONE
   
   CLASS(JO_TABLE), INTENT(IN)           :: YDJOT
   INTEGER(KIND=JPIM), INTENT(IN)        :: KNSMAX
   
   REAL(KIND=JPRB)    :: ZJOSUM,ZCOST,ZRMS_OBSERR,ZRMS_BGERR,ZJOTOT
   REAL(KIND=JPRB)    :: ZSTATUS(JPST), ZEVENT1(NODEV1)
   INTEGER(KIND=JPIM) :: JOTP, JCTP, JVAR, IJOSUM,ICOUNT,IJOTOT,JBIN,ISTAT
   INTEGER(KIND=JPIM) :: IULTMP
   CHARACTER(LEN=16)  :: CLJOB
   CHARACTER(LEN=10)  :: CLV
   REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
   
   !     ------------------------------------------------------------------
   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:PRINT_JO',0,ZHOOK_HANDLE)
   !     ------------------------------------------------------------------
   
   ! Print the JO table
   IF(LSCREEN) THEN
     CLJOB='SCREENING JOB   '
   ELSEIF(L_OBS_IN_FC()) THEN
     CLJOB='TRAJECTORY JOB  '
   ELSE
     CLJOB='MINIMISATION JOB'
   ENDIF
   
   WRITE(NULOUT,'(/,A,A,A,I4.4,3(A,I5))') &
    & ' Diagnostic JO-table (JOT) ',CLJOB,' T',KNSMAX,&
    & ' NCONF= ',NCONF,' NSIM4D= ',GET_NSIM4D(),' NUPTRA= ',GET_NUPTRA()
   WRITE(NULOUT,*) &
    & '====================================================================='//&
    & '===================='
   
   ZRMS_BGERR=0.0_JPRB
   ZJOTOT=0.0_JPRB
   IJOTOT=0
   DO JOTP=1,NMXOTP
     ZJOSUM=0.0_JPRB
     IJOSUM=0
     WRITE(NULOUT,*)
     WRITE(NULOUT,'(6X,A,I5,A,A32)') &
      & 'Obstype ',JOTP,' === ',YDJOT%JOT(JOTP)%NAME  
     WRITE(NULOUT,'(6X,A)') &
      & '--------------------------------------------------'  
   
     DO JCTP=1,YDJOT%JOT(JOTP)%SIZE
       IF(ANY(YDJOT%JOT(JOTP)%JOST(JCTP)%COUNT(:,:) /= 0)) THEN
         WRITE(NULOUT,'(8X,A,I5,A,A32)') &
          & 'Codetype ',YDJOT%JOT(JOTP)%JOST(JCTP)%CODETYPE,' === ',&
          & YDJOT%JOT(JOTP)%JOST(JCTP)%SNAME  
         WRITE(NULOUT,'(10X,2A)') &
          & 'Variable          DataCount          Jo_Costfunction'&
          & ,'         JO/n       ObsErr      BgErr'  
       ENDIF
       DO JVAR=1,YDJOT%NVAR%JPXVAR
         ICOUNT  = SUM(YDJOT%JOT(JOTP)%JOST(JCTP)%COUNT(JVAR,:))
         ZCOST   = SUM(YDJOT%JOT(JOTP)%JOST(JCTP)%COST(JVAR,:))
         ZRMS_OBSERR = YDJOT%JOT(JOTP)%JOST(JCTP)%OBSERR(JVAR)
         CALL YDJOT%NVAR%GET_NAME(JVAR,CLV)
         IF(LGUESS) ZRMS_BGERR  = YDJOT%JOT(JOTP)%JOST(JCTP)%BGERR(JVAR)
         IF(ICOUNT /= 0) THEN
           IF(ZRMS_OBSERR /= 0) &
            & ZRMS_OBSERR = SQRT(ZRMS_OBSERR/ICOUNT)
           IF(LGUESS .AND. ZRMS_BGERR > 0.0_JPRB .AND. ZRMS_BGERR /= RMDI) &
            & ZRMS_BGERR  = SQRT(ZRMS_BGERR /ICOUNT)
           IF(ZRMS_BGERR /= RMDI) THEN
             WRITE(NULOUT,'(12X,A10,3X,I10,3X,G27.13,3X,F8.2,1X,2(2x,E10.3))') &
            & CLV,ICOUNT,ZCOST,ZCOST/ICOUNT,ZRMS_OBSERR,ZRMS_BGERR
           ELSE
             WRITE(NULOUT,'(12X,A10,3X,I10,3X,G27.13,3X,F8.2,1X,3X,E10.3,4X,A4)') &
            & CLV,ICOUNT,ZCOST,ZCOST/ICOUNT,ZRMS_OBSERR,'RMDI'
           ENDIF
           IJOSUM=IJOSUM + ICOUNT
           ZJOSUM=ZJOSUM + ZCOST
         ENDIF
       ENDDO
     ENDDO
     ZJOTOT=ZJOTOT+ZJOSUM
     IJOTOT=IJOTOT+IJOSUM
     IF(IJOSUM /= 0) THEN
       WRITE(NULOUT,'(25X,A)') &
        & '----------   ---------------------------   --------'  
       WRITE(NULOUT,'(5X,A,I2,A,I10,3X,G27.13,3X,F8.2)') &
        & 'ObsType ',JOTP,' Total: ',IJOSUM,ZJOSUM,ZJOSUM/IJOSUM  
     ENDIF
   ENDDO
   WRITE(NULOUT,*)
   WRITE(NULOUT,*) &
    & ' ---------------------------------------------------------------------------'  
   IF (IJOTOT > 0) THEN
     WRITE(NULOUT,'(8X,A,2x,I10,3X,G27.13,3X,F8.2)') &
      & 'Jo Global : ',IJOTOT,ZJOTOT,ZJOTOT/IJOTOT  
   ELSE
     WRITE(NULOUT,'(8X,A,2x,I10,3X,G27.13,3X,A8)') &
      & 'Jo Global : ',IJOTOT,ZJOTOT,'********'  
   ENDIF
   DIAG_4DVAR%NOBS=IJOTOT
   IF (MYPROC == 1 .AND. .NOT. L_OOPS) THEN
     IULTMP=RESERVE_LUN()
     OPEN(UNIT=IULTMP,FILE='diagvar.txt')
     WRITE(IULTMP,'(F20.5,I10,I2,7F20.8)') DIAG_4DVAR%CONDNUM,DIAG_4DVAR%NOBS,DIAG_4DVAR%NCOST, &
    &          DIAG_4DVAR%JO,DIAG_4DVAR%JB,DIAG_4DVAR%JC,DIAG_4DVAR%JQ,DIAG_4DVAR%JP,DIAG_4DVAR%JH,DIAG_4DVAR%JCVARBC 
     CLOSE(IULTMP)
     CALL FREE_LUN(IULTMP)
   ENDIF
   WRITE(NULOUT,*) &
    & ' ==========================================================================='  
   WRITE(NULOUT,*)
   WRITE(NULOUT,*) ' End of JO-table (JOT)'
   
   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:PRINT_JO',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_JO

SUBROUTINE PRINT_SCREEN_STATS(YDJOT)

   !**** Print screening statistics 
   CLASS(JO_TABLE), INTENT(IN)           :: YDJOT

   REAL(KIND=JPRB)    :: ZJOSUM,ZCOST,ZRMS_OBSERR,ZRMS_BGERR,ZJOTOT
   REAL(KIND=JPRB)    :: ZSTATUS(JPST), ZEVENT1(NODEV1)
   INTEGER(KIND=JPIM) :: JOTP, JCTP, JVAR, ICOUNT,ISTAT
   INTEGER(KIND=JPIM) :: IEV(NODEV1)
   CHARACTER(LEN=10)  :: CLV
   LOGICAL            :: LLEMPTY
   REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:PRINT_SCREEN_STATS',0,ZHOOK_HANDLE)

   ! Prepare event number list for printing later
   DO ISTAT=1,NODEV1
     IEV(ISTAT)=ISTAT
   ENDDO

   WRITE(NULOUT,*)
   WRITE(NULOUT,*) &
    & ' Observation usage summary '
   WRITE(NULOUT,*) &
    & '====================================================================='//&
    & '===================='

   DO JOTP=1,NMXOTP

     LLEMPTY = .TRUE.
     DO JCTP=1,YDJOT%JOT(JOTP)%SIZE
       IF(ANY(YDJOT%JOT(JOTP)%JOST(JCTP)%COUNT(:,:) /= 0)) THEN
         LLEMPTY = .FALSE.
       ENDIF
     ENDDO
     IF (LLEMPTY) CYCLE

     WRITE(NULOUT,*)
     WRITE(NULOUT,'(8X,A,I5,A,A32)') &
      & 'Obstype ',JOTP,' === ',YDJOT%JOT(JOTP)%NAME
     WRITE(NULOUT,'(8X,A)') &
      & '--------------------------------------------------'
     WRITE(NULOUT,'(2A)') &
      & '                                                         datum_status %       ',&
      & '       datum_event1 %    (rounded UP to nearest integer, see key below)'
     WRITE(NULOUT,'(A,32I4)') &
      & 'Codetype                         Variable      DataCount    Activ  Pass   Rej Black   ',&
      & IEV

     DO JCTP=1,YDJOT%JOT(JOTP)%SIZE
       IF(ANY(YDJOT%JOT(JOTP)%JOST(JCTP)%COUNT(:,:) /= 0)) THEN
         WRITE(NULOUT,*)
       ENDIF
       DO JVAR=1,YDJOT%NVAR%JPXVAR
         CALL YDJOT%NVAR%GET_NAME(JVAR,CLV)
         ICOUNT  = SUM(YDJOT%JOT(JOTP)%JOST(JCTP)%COUNTALL(JVAR,:))
         IF(ICOUNT /= 0) THEN
           DO ISTAT = 1, JPST
             ZSTATUS(ISTAT) = 100.0_JPRB * REAL(YDJOT%JOT(JOTP)%JOST(JCTP)%STATUS(JVAR,ISTAT),JPRB) &
                          & / REAL(ICOUNT,JPRB)
           ENDDO
           DO ISTAT = 1, NODEV1
             ZEVENT1(ISTAT) = 100.0_JPRB * REAL(YDJOT%JOT(JOTP)%JOST(JCTP)%EVENT1(JVAR,ISTAT),JPRB) &
                          & / REAL(ICOUNT,JPRB)
           ENDDO
           WRITE(NULOUT,'(I4,X,A30,X,A10,I10,3X,4F6.1,3X,32I4)') &
            & YDJOT%JOT(JOTP)%JOST(JCTP)%CODETYPE, YDJOT%JOT(JOTP)%JOST(JCTP)%SNAME, &
            & CLV,ICOUNT,ZSTATUS,CEILING(ZEVENT1)
         ENDIF
       ENDDO
     ENDDO
   ENDDO
   WRITE(NULOUT,*)

   WRITE(NULOUT,*) &
    & ' ==========================================================================='
   WRITE(NULOUT,*) 'datum_event1 key'
   WRITE(NULOUT,*)
   DO ISTAT = 1,NODEV1
     WRITE(NULOUT,*) CH1EVENT(ISTAT)
   ENDDO
   WRITE(NULOUT,*) &
    & ' ==========================================================================='
   WRITE(NULOUT,*)
   WRITE(NULOUT,*) ' End of observation usage summary'
   WRITE(NULOUT,*)

   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:PRINT_SCREEN_STATS',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_SCREEN_STATS

SUBROUTINE UPDATE(THIS,ZJO_PER_DATUM,YDSET,ZVARNOS,ZFINAL_OBS_ERRORS,ZFG_ERRORS,ZCODETYPES,LLCOMBINED_UV)
   ! Updates Jo table entries given Jo values for all observations in a given set
   ! (code extracted from ifs/module/hjo.F90)
   CLASS(JO_TABLE), INTENT(INOUT)   :: THIS
   REAL(KIND=JPRB), INTENT(IN)      :: ZJO_PER_DATUM(:)
   TYPE(TYPE_SET_INFO), INTENT(IN)  :: YDSET
   REAL(KIND=JPRD), INTENT(IN)      :: ZVARNOS(:)
   REAL(KIND=JPRD), INTENT(IN)      :: ZFINAL_OBS_ERRORS(:)
   REAL(KIND=JPRD), INTENT(IN)      :: ZFG_ERRORS(:)
   REAL(KIND=JPRB), INTENT(IN)      :: ZCODETYPES(:)
   LOGICAL, INTENT(IN)              :: LLCOMBINED_UV
   INTEGER(KIND=JPIM)               :: IBODY
   INTEGER(KIND=JPIM)               :: IIVAR
   INTEGER(KIND=JPIM)               :: IIVNM
   INTEGER(KIND=JPIM)               :: ISTATUS
   INTEGER(KIND=JPIM)               :: IGRP
   INTEGER(KIND=JPIM)               :: ICTPSQ
   INTEGER(KIND=JPIM)               :: ITIME
   INTEGER(KIND=JPIM)               :: IADD
   LOGICAL                          :: LLU
   REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:UPDATE',0,ZHOOK_HANDLE)

   DO IBODY=1,YDSET%SUMBDY
        IF(ZJO_PER_DATUM(IBODY) == RMDI) CYCLE

        ! Find nvar corresponding to this varno
        IIVNM = NINT(ZVARNOS(IBODY),JPIM)
        CALL THIS%NVAR%MAP_VARNO(IIVNM,IIVAR,ISTATUS)
        IF (ISTATUS/=0) THEN
          ! This is not an error - some varnos do not map to nvar
          !WRITE(NULERR,*)'JO_TABLE % UPDATE: ibody,iivar,NVAR%JPXVAR,zvarno,obstype=',IBODY,IIVAR,NVAR%JPXVAR,ZVARNOS(IBODY),&
          !               &  YDSET%OBSTYPE, MAXVAL(ZCODETYPES),MINVAL(ZCODETYPES)
          CYCLE
        ENDIF

        ! Support combined u,v Jo values
        IADD = 1
        IF(LLCOMBINED_UV)THEN
          LLU=IIVNM == VARNO%U .OR. IIVNM == VARNO%U10M .OR. IIVNM == VARNO%SCATU
          IF(LLU) IADD = 2   ! If combined u,v and this is a wind observation
        ENDIF

        ! Pick up obs codes associated with this set
        IF(YDSET%OBSTYPE == NSATEM .OR. YDSET%OBSTYPE == NSATOB .OR. &
           & YDSET%OBSTYPE == NLIMB .OR. YDSET%OBSTYPE == NALLSKY .OR. YDSET%OBSTYPE == NSCATT) THEN
          IGRP = YDSET%GROUP
          ICTPSQ  = IGRP
        ELSE ! Conventional observations can have more than one codetype per set
          ICTPSQ=MSQBYCTP(NINT(ZCODETYPES(IBODY),JPIM))
        ENDIF

        ! Assign time bin
        ITIME = 1
!        ZTIME=PTIME(JOBS) * REAL(THIS%NPRTBINS,JPRB)
!        ITIME=INT(ZTIME)+1
!        IF (ITIME<1) ITIME=1
!        IF (ITIME>THIS%NPRTBINS) ITIME=THIS%NPRTBINS

        ! Check JO_TABLE structure is initialised as required
        IF (ALLOCATED(THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COST)) THEN
          IF (SIZE(THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COST)==0) THEN
            WRITE(NULERR,*)'JO_TABLE % UPDATE: YDSET%OBSTYPE,ICTPSQ,size(cost),NVAR%JPXVAR,NPRTBINS=', &
             & YDSET%OBSTYPE,ICTPSQ,SHAPE(THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COST),THIS%NVAR%JPXVAR,THIS%NPRTBINS
            CALL ABOR1('JO_TABLE % UDPATE: cost size 0')
          ENDIF
        ELSE
          WRITE(NULERR,*)'JO_TABLE % UPDATE: YDSET%OBSTYPE,ICTPSQ,NVAR%JPXVAR,NPRTBINS=', &
                             & YDSET%OBSTYPE,ICTPSQ,THIS%NVAR%JPXVAR,THIS%NPRTBINS
          CALL ABOR1('JO_TABLE % UPDATE: cost not allocated')
        ENDIF
        IF (ALLOCATED(THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COUNT)) THEN
          IF (SIZE(THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COUNT)==0) THEN
            WRITE(NULERR,*)'JO_TABLE % UDPATE: YDSET%OBSTYPE,ICTPSQ,size(count),NVAR%JPXVAR,NPRTBINS=', &
             & YDSET%OBSTYPE,ICTPSQ,SHAPE(THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COUNT),THIS%NVAR%JPXVAR,THIS%NPRTBINS
            CALL ABOR1('JO_TABLE % UPDATE: count size 0')
          ENDIF
        ELSE
          WRITE(NULERR,*)'JO_TABLE % UPDATE: YDSET%OBSTYPE,ICTPSQ,NVAR%JPXVAR,NPRTBINS=', &
                             & YDSET%OBSTYPE,ICTPSQ,THIS%NVAR%JPXVAR,THIS%NPRTBINS
          CALL ABOR1('JO_TABLE % UDPATE: count not allocated')
        ENDIF

        ! Update JO-table using data from this set
        THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COST(IIVAR,ITIME)= &
         & THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COST(IIVAR,ITIME)+ZJO_PER_DATUM(IBODY)

        THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COUNT(IIVAR,ITIME)= &
         & THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%COUNT(IIVAR,ITIME) + IADD

        THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%OBSERR(IIVAR)= &
         & THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%OBSERR(IIVAR) + SUM(ZFINAL_OBS_ERRORS(IBODY:IBODY+IADD-1)**2)

        IF(THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%BGERR(IIVAR) == RMDI) THEN
          ! DO Nothing
        ELSEIF(ZFG_ERRORS(IBODY) == RMDI) THEN
          THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%BGERR(IIVAR) = RMDI
        ELSE
          THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%BGERR(IIVAR)= &
           & THIS%JOT(YDSET%OBSTYPE)%JOST(ICTPSQ)%BGERR(IIVAR) + SUM(ZFG_ERRORS(IBODY:IBODY+IADD-1)**2)
        ENDIF

    ENDDO

    IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:UPDATE',1,ZHOOK_HANDLE)

END SUBROUTINE UPDATE

SUBROUTINE GATHER(YDJOT)

   !**** Gather JO_TABLE contributions
   
   !     Purpose.
   !     --------
   !           Build up global diagnostic Jo cost function arrays
   
   !**   Interface.
   !     ----------
   !        *CALL* *YDJOT%GATHER(..)
   
   !        Explicit arguments :
   !        --------------------
   
   
   !        Implicit arguments :
   !        --------------------
   !        None
   
   !     Method.
   !     -------
   !        See documentation
   
   !     Externals.
   !     ----------
   
   !     Reference.
   !     ----------
   !        ECMWF Research Department documentation of the IFS
   
   !     Author.
   !     -------
   !      Mats Hamrud ECMWF
   !      Original  : 96-04-04
   
   !     Modifications.
   !     --------------
   !      M.Hamrud      01-Oct-2003 CY28 Cleaning
   !      Y.Tremolet    29-May-2008 Print Jo table by bins
   !      A. Geer          22-Nov-2010 Event summary table 
   !      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
   !      P. Lean       29-Mar-2017  Moved into jo_table_mod and renamed from gathercosto to gather
   !     ------------------------------------------------------------------
   
   IMPLICIT NONE
   
   CLASS(JO_TABLE), INTENT(INOUT)  :: YDJOT
   
   REAL(KIND=JPRB) :: ZCOMSEND(YDJOT%NJOTTSZ)
   REAL(KIND=JPRB),ALLOCATABLE :: ZCOMRECV(:)
   
   INTEGER(KIND=JPIM) :: ISZSEND,  IOFF
   INTEGER(KIND=JPIM) :: JROC, JOTP, JCTP, JBIN
   INTEGER(KIND=JPIM) :: IRECVCOUNTS(NPROC)
   REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
   
#include "abor1.intfb.h"
   
   !     ------------------------------------------------------------------
   
   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:GATHER',0,ZHOOK_HANDLE)
   ISZSEND = SIZE(ZCOMSEND)
   
   !     1.    Store this PE's cost function contributions in global array
   !           -----------------------------------------------------------
   
   IF(NPROC == 1) THEN
   
   !     2.    Make global view
   !           ----------------
   
   ELSE
     IF(MYPROC == 1) THEN
       ALLOCATE(ZCOMRECV(NPROC*YDJOT%NJOTTSZ))
     ELSE
       ALLOCATE(ZCOMRECV(1))
     ENDIF
       
     IOFF=0
     DO JOTP=1,NMXOTP
       DO JCTP=1,YDJOT%JOT(JOTP)%SIZE
         DO JBIN=1,YDJOT%NPRTBINS
           ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%COST(:,JBIN)
           IOFF=IOFF+YDJOT%NVAR%JPXVAR
           ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%COUNT(:,JBIN)
           IOFF=IOFF+YDJOT%NVAR%JPXVAR
           ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%COUNTALL(:,JBIN)
           IOFF=IOFF+YDJOT%NVAR%JPXVAR
           ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%OBIMPACT(:,JBIN)
           IOFF=IOFF+YDJOT%NVAR%JPXVAR
         ENDDO
         IF (LSCREEN) THEN
           DO JBIN=1,JPST
             ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%STATUS(:,JBIN)
             IOFF=IOFF+YDJOT%NVAR%JPXVAR
           ENDDO
           DO JBIN=1,NODEV1
             ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%EVENT1(:,JBIN)
             IOFF=IOFF+YDJOT%NVAR%JPXVAR
           ENDDO
         ENDIF
         ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%OBSERR(:)
         IOFF=IOFF+YDJOT%NVAR%JPXVAR
         ZCOMSEND(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)=YDJOT%JOT(JOTP)%JOST(JCTP)%BGERR(:)
         IOFF=IOFF+YDJOT%NVAR%JPXVAR
       ENDDO
     ENDDO
     IF(IOFF /= ISZSEND) CALL ABOR1('GATHERCOSTO: ISZSEND /= IOFF')
   
     DO JROC=1,NPROC
       IRECVCOUNTS(JROC)=YDJOT%NJOTTSZ
     ENDDO
     
     CALL GSTATS_BARRIER(712)
     CALL GSTATS(610,0)
     CALL MPL_GATHERV(PSENDBUF=ZCOMSEND,PRECVBUF=ZCOMRECV,&
      & KRECVCOUNTS=IRECVCOUNTS,KROOT=1,CDSTRING='JO_TABLE_MOD:GATHER')  
     CALL GSTATS(610,1)
     CALL GSTATS_BARRIER2(712)
     
     IF(MYPROC == 1) THEN
       IOFF=0
       DO JROC=1,NPROC
         IF(JROC /= MYPROC ) THEN
           CALL GSTATS(1826,0)
           DO JOTP=1,NMXOTP
             DO JCTP=1,YDJOT%JOT(JOTP)%SIZE
               DO JBIN=1,YDJOT%NPRTBINS
                 YDJOT%JOT(JOTP)%JOST(JCTP)%COST(:,JBIN)  = &
                  & YDJOT%JOT(JOTP)%JOST(JCTP)%COST(:,JBIN) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)
                 IOFF=IOFF+YDJOT%NVAR%JPXVAR
                 YDJOT%JOT(JOTP)%JOST(JCTP)%COUNT(:,JBIN) = &
                  & YDJOT%JOT(JOTP)%JOST(JCTP)%COUNT(:,JBIN) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)
                 IOFF=IOFF+YDJOT%NVAR%JPXVAR
                 YDJOT%JOT(JOTP)%JOST(JCTP)%COUNTALL(:,JBIN) = &
                  & YDJOT%JOT(JOTP)%JOST(JCTP)%COUNTALL(:,JBIN) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)
                 IOFF=IOFF+YDJOT%NVAR%JPXVAR
                 YDJOT%JOT(JOTP)%JOST(JCTP)%OBIMPACT(:,JBIN) = &
                  & YDJOT%JOT(JOTP)%JOST(JCTP)%OBIMPACT(:,JBIN) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)  
                 IOFF=IOFF+YDJOT%NVAR%JPXVAR
               ENDDO
               IF (LSCREEN) THEN
                 DO JBIN=1,JPST
                   YDJOT%JOT(JOTP)%JOST(JCTP)%STATUS(:,JBIN) = &
                    & YDJOT%JOT(JOTP)%JOST(JCTP)%STATUS(:,JBIN) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)
                   IOFF=IOFF+YDJOT%NVAR%JPXVAR
                 ENDDO
                 DO JBIN=1,NODEV1
                   YDJOT%JOT(JOTP)%JOST(JCTP)%EVENT1(:,JBIN) = &
                    & YDJOT%JOT(JOTP)%JOST(JCTP)%EVENT1(:,JBIN) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)
                   IOFF=IOFF+YDJOT%NVAR%JPXVAR
                 ENDDO
               ENDIF
               YDJOT%JOT(JOTP)%JOST(JCTP)%OBSERR(:) = &
                & YDJOT%JOT(JOTP)%JOST(JCTP)%OBSERR(:) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)
               IOFF=IOFF+YDJOT%NVAR%JPXVAR
               WHERE(ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR) == RMDI)
                  YDJOT%JOT(JOTP)%JOST(JCTP)%BGERR(:) = RMDI
               ELSEWHERE
                  YDJOT%JOT(JOTP)%JOST(JCTP)%BGERR(:) = &
                &   YDJOT%JOT(JOTP)%JOST(JCTP)%BGERR(:) + ZCOMRECV(IOFF+1:IOFF+YDJOT%NVAR%JPXVAR)
               ENDWHERE
               IOFF=IOFF+YDJOT%NVAR%JPXVAR
             ENDDO
           ENDDO
           CALL GSTATS(1826,1)
         ELSE
           IOFF=IOFF+YDJOT%NJOTTSZ
         ENDIF

       ENDDO
     ENDIF
     DEALLOCATE(ZCOMRECV)
   ENDIF

   IF (LHOOK) CALL DR_HOOK('JO_TABLE_MOD:GATHER',1,ZHOOK_HANDLE)
END SUBROUTINE GATHER

SUBROUTINE CREATE(YDJOTSUB,KPXVAR,KPRTBINS)
   CLASS(JOTSUB_T),    INTENT(INOUT) :: YDJOTSUB
   INTEGER(KIND=JPIM), INTENT(IN)    :: KPXVAR,KPRTBINS

   REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

   IF (LHOOK) CALL DR_HOOK('JO_TABLE:JOTSUB_T:CREATE',0,ZHOOK_HANDLE)

   ALLOCATE(YDJOTSUB%COST    (KPXVAR,KPRTBINS))
   ALLOCATE(YDJOTSUB%COUNT   (KPXVAR,KPRTBINS))
   ALLOCATE(YDJOTSUB%COUNTALL(KPXVAR,KPRTBINS))
   ALLOCATE(YDJOTSUB%EVENT1  (KPXVAR,NODEV1))
   ALLOCATE(YDJOTSUB%OBIMPACT(KPXVAR,KPRTBINS))
   ALLOCATE(YDJOTSUB%OBSERR  (KPXVAR))
   ALLOCATE(YDJOTSUB%BGERR   (KPXVAR))
   ALLOCATE(YDJOTSUB%STATUS  (KPXVAR,JPST))

   YDJOTSUB%COST(:,:)= 0.0_JPRB
   YDJOTSUB%OBSERR(:)= 0.0_JPRB
   YDJOTSUB%BGERR(:) = 0.0_JPRB
   YDJOTSUB%OBIMPACT(:,:) = 0.0_JPRB
   YDJOTSUB%COUNT(:,:) = 0
   YDJOTSUB%COUNTALL(:,:) = 0
   YDJOTSUB%STATUS(:,:) = 0
   YDJOTSUB%EVENT1(:,:) = 0
     
   IF (LHOOK) CALL DR_HOOK('JO_TABLE:JOTSUB_T:CREATE',1,ZHOOK_HANDLE)
END SUBROUTINE CREATE

END MODULE JO_TABLE_MOD
