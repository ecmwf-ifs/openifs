! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE IFSTOWAM (BLK2LOC,                      &
 &                   NFIELDS, NGPTOTG, NCA, NRA,   &
 &                   FIELDS, LWCUR, MASK_IN,       &
 &                   NXS, NXE, NYS, NYE, FIELDG)

!***  *IFSTOWAM* - INTERPOLATES FORCING FIELDS FROM IFS TO WAM GRID. 


!     PURPOSE                                                       
!     -------                                                      

!     *IFSTOWAM* - INTERPOLATES FORCING FIELDS FROM IFS TO WAM GRID. 
!                  ONLY FOR COUPLED RUNS.

!     INTERFACE                                                         
!     ---------                                                         

!     *CALL* *IFSTOWAM (BLK2LOC,
!    &                  NFIELDS, NGPTOTG, NCA, NRA,
!    &                  FIELDS, LWCUR, MASK_IN,
!    &                  NXS, NXE, NYS, NYE, FIELDG)
!
!        *BLK2LOC*- POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!        *NFIELDS*- NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!        *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *FIELDS* - ATMOSPHERIC FIELDS AS FOLLOWS:
!                   FIELDS(:,1) = U COMPONENT OF WIND SPEED
!                   FIELDS(:,2) = V COMPONENT OF WIND SPEED
!                   FIELDS(:,3) = AIR DENSITY
!                   FIELDS(:,4) = w* USED FOR GUSTINESS
!                   FIELDS(:,5) = SEA ICE FRACTION 
!                   FIELDS(:,6) = LAKE FRACTION 
!                   FIELDS(:,7) = U COMPONENT OF SURFACE CURRENT 
!                   FIELDS(:,8) = V COMPONENT OF SURFACE CURRENT 
!        *LWCUR* - LOGICAL CONTROLLING THE PRESENCE OF MEANINGFUL
!                   SURFACE CURRENTS (I.E. NOT ALL ZEROS).
!        *MASK_IN*- MASK FOR FIELDS TO ONLY POINTS TO THE PART
!                   OF FIELDS THAT ARE NEEDED ON EACH PE.
!        *NXS:NXE*  FIRST DIMENSION OF FIELDG
!        *NYS:NYE*  SECOND DIMENSION OF FIELDG
!        *FIELDG*   INPUT FORCING FIELDS ON THE WAVE MODEL GRID

!     EXTERNALS                                                         
!     ---------                                                         

!     *ABORT1*         TERMINATES PROCESSING
!     *INCDATE*        INCREMENTS DATE 
!     *INITIALINT*     COMPUTES INTERPOLATION COEFFICIENTS (FOR FLDINTER)
!     *FLDINTER*       INTERPOLATES ATMOSPHERIC FIELDS TO WAM GRID

! --------------------------------------------------------------------- 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, FORCING_FIELDS

      USE YOWCOUP  , ONLY : LWCOU    ,LWCOUNORMS, I_MASK_IN, N_MASK_IN, &
     &            LLNORMIFS2WAM      ,LWCOUSAMEGRID
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK
      USE YOWICE   , ONLY : IPARAMCI
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,ZDELLO   ,NLONRGG
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOMMP0   , ONLY : MYSETW   ,NPRCIDS  ,NPRTRW   ,NPRTRV
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO    ,LL1D
      USE YOWPCONS , ONLY : ZMISS    ,ROAIR    ,WSTAR0
      USE YOWSTAT  , ONLY : LADEN    ,LGUST
      USE YOWTEST  , ONLY : IU06
      USE YOWWNDG  , ONLY : ICODE_CPL,IWPER    ,ICOORD 
      USE YOWWIND  , ONLY : NC       ,NR       , CWDFILE  ,LLNEWCURR

      USE YOMTAG   , ONLY : MTAGWAMNORM
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, MPL_GATHERV, &
                           & JP_NON_BLOCKING_STANDARD, &
                           & MPL_ALL_LEVS_COMM

! --------------------------------------------------------------------- 

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "fldinter.intfb.h"
#include "initialint.intfb.h"

      TYPE(WVGRIDLOC), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: BLK2LOC
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS, NGPTOTG, NCA, NRA
      REAL(KIND=JWRB), INTENT(IN) :: FIELDS(NGPTOTG,NFIELDS)
      LOGICAL, INTENT(IN) :: LWCUR
      INTEGER(KIND=JWIM), INTENT(INOUT) :: MASK_IN(NGPTOTG)
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), DIMENSION(NXS:NXE, NYS:NYE), INTENT(INOUT) :: FIELDG


      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL
      INTEGER(KIND=JWIM) :: IJ, I, J, JSN, JF, L, IC, IFLD, ICOUNT, IST, IP
      INTEGER(KIND=JWIM) :: NCOMLOC, NCOMBUF, NMASK, NTOT
      INTEGER(KIND=JWIM) :: IMASTER, IR, ITOT, IINC, IRECV, JSETW, ITAG
      INTEGER(KIND=JWIM) :: IBEG, IEND, IBEGOFF, IENDOFF, IA, IB, ISETW, ISETV 
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS
      INTEGER(KIND=JWIM), SAVE :: IPERIODIC
      INTEGER(KIND=JWIM) :: IREQ(NPROC+1)
      INTEGER(KIND=JWIM) :: IRECVCOUNTS(NPROC)
      INTEGER(KIND=JWIM), ALLOCATABLE, SAVE :: JJ(:), II(:,:), IIP(:,:)
      INTEGER(KIND=JWIM), ALLOCATABLE, SAVE :: ISTART(:)
      INTEGER(KIND=JWIM), ALLOCATABLE, SAVE :: ILONRGG(:), IJBLOCK(:,:)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZDUM(2)
      REAL(KIND=JWRB) :: VAL 
      REAL(KIND=JWRB), DIMENSION(NFIELDS,4) :: NORMS ! 1:average, 2:min, 3:max, 4:variance 
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUF(:),ZCOMBUF1(:)
      REAL(KIND=JWRB), ALLOCATABLE :: ZBUFW(:)
      REAL(KIND=JWRB), ALLOCATABLE, SAVE :: DJ1(:), DII1(:,:), DIIP1(:,:)

      CHARACTER(LEN=16) :: AFLABEL(NFIELDS)

      LOGICAL, SAVE :: LLINTERPOL
      LOGICAL, SAVE :: FRSTMATM
      LOGICAL, SAVE :: LFRSTIME
      LOGICAL, SAVE :: LFRST

      LOGICAL :: LLKC

      DATA FRSTMATM / .TRUE. /
      DATA LFRSTIME / .TRUE. /
      DATA LFRST /.TRUE./

! --------------------------------------------------------------------  

IF (LHOOK) CALL DR_HOOK('IFSTOWAM',0,ZHOOK_HANDLE)

ASSOCIATE(IFROMIJ => BLK2LOC%IFROMIJ, &
 &        JFROMIJ => BLK2LOC%JFROMIJ)

      CALL GSTATS_BARRIER(796)

!     ONLY ACTIVE IF COUPLED TO IFS.
!     -----------------------------
      IF (LWCOU) THEN
        IF (FRSTMATM) THEN
!        
!         1.1 FIRST TIME ONLY, COMPUTE INTERPOLATION COEFFICIENTS:
!             ----------------------------------------------------
!
          IF (.NOT.ALLOCATED(ILONRGG)) ALLOCATE (ILONRGG(NRA))
          IF (.NOT.ALLOCATED(JJ)) ALLOCATE (JJ(NGY))
          IF (.NOT.ALLOCATED(II)) ALLOCATE (II(NGX,NGY))
          IF (.NOT.ALLOCATED(IIP)) ALLOCATE (IIP(NGX,NGY))
          IF (.NOT.ALLOCATED(DJ1)) ALLOCATE (DJ1(NGY))
          IF (.NOT.ALLOCATED(DII1)) ALLOCATE (DII1(NGX,NGY))
          IF (.NOT.ALLOCATED(DIIP1)) ALLOCATE (DIIP1(NGX,NGY))

          
          CALL INITIALINT(IU06, NCA, NRA,                               &
     &                    NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,      &
     &                    AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,    &
     &                    ILONRGG,                                      &
     &                    LLINTERPOL,DJ1, DII1, DIIP1, JJ, II, IIP)

          LWCOUSAMEGRID = .NOT. LLINTERPOL
          WRITE(IU06,*) ' SUB. IFSTOWAM - LWCOUSAMEGRID = ',LWCOUSAMEGRID

          IF (.NOT.ALLOCATED(IJBLOCK)) ALLOCATE(IJBLOCK(0:NCA+1,NRA))

          L = 0
          DO J=1,NRA
            JSN=NRA-J+1
            DO I=1,ILONRGG(JSN)
              L = L+1
              IJBLOCK(I,J) = L
            ENDDO
          ENDDO
          IF (NGPTOTG /= L)THEN
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ####  WARNING IN  IFSTOWAM WARNING  ###'
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ##  NGPTOTG  IS  NOT  EQUAL TO  L  ###'
            WRITE(IU06,*)' ##  NGPTOTG = ', NGPTOTG
            WRITE(IU06,*)' ##        L = ', L
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ##  NCA = ', NCA, '   NRA = ', NRA
            WRITE(IU06,*)' ##  SUM ILONRGG = ',SUM(ILONRGG(1:NRA))
            WRITE(IU06,*)' ##  ILONRGG = ', ILONRGG(1:NRA)
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' #######################################'
            CALL FLUSH(IU06)
          ENDIF

          IF (IPERIODIC == 1) THEN
            DO J=1,NRA
            JSN=NRA-J+1
            IJBLOCK(0,J)= IJBLOCK(ILONRGG(JSN),J)
            IJBLOCK(ILONRGG(JSN)+1,J)= IJBLOCK(1,J)
            ENDDO
          ENDIF
     
          FRSTMATM = .FALSE.
          
          IF (LADEN .AND. NFIELDS < 3) THEN 
            LADEN = .FALSE.
            WRITE(IU06,*) ' WAM_IFSTOWAM - WARNING: '
            WRITE(IU06,*) '      AIR DENSITY RUN WAS REQUESTED BUT '    &
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 3'
            WRITE(IU06,*) '      RESET  LADEN  TO  .FALSE.'
          ENDIF
          
          IF (LADEN) THEN 
            WRITE(IU06,*) ' WAM_IFSTOWAM - VARIABLE AIR DENSITY RUN'
          ELSE
            WRITE(IU06,*) ' WAM_IFSTOWAM - AIR DENSITY IS '             &
     &                     //' CONSTANT = ', ROAIR
          ENDIF
          
          IF (LGUST .AND. NFIELDS < 4) THEN 
            LGUST = .FALSE.
            WRITE(IU06,*) ' SUB. IFSTOWAM - WARNING: '
            WRITE(IU06,*) '      GUSTINESS RUN WAS REQUESTED BUT '      &
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 4'
            WRITE(IU06,*) '      RESET  LGUST  TO  .FALSE.'
          ENDIF
          
          IF (LGUST) THEN 
            WRITE(IU06,*) ' WAM_IFSTOWAM - GUSTINESS RUN'
          ELSE
            WRITE(IU06,*) ' WAM_IFSTOWAM - NO GUSTINESS EFFECT.'
          ENDIF

        ENDIF  ! FRSTMATM

!        
!       COMPUTATION OF THE NORMS OF INPUT FIELDS
!
        IF ( LFRSTIME .OR. LWCOUNORMS .OR. LLNORMIFS2WAM )THEN
          SELECT CASE (ICODE_CPL)
          CASE(0)
            AFLABEL(1)='PROBLEM !!!!!!! '
            AFLABEL(2)='PROBLEM !!!!!!! '
          CASE(1)
            AFLABEL(1)='u* x-component  '
            AFLABEL(2)='u* y-component  '
          CASE(2)
            AFLABEL(1)='tau x-component '
            AFLABEL(2)='tau y-component '
          CASE(3)
            AFLABEL(1)='Neutral 10m U   '
            AFLABEL(2)='Neutral 10m V   '
          END SELECT 

          IF (NFIELDS >= 3) AFLABEL(3)='Air density     '
          IF (NFIELDS >= 4) AFLABEL(4)='w*              '
          IF (NFIELDS >= 5) AFLABEL(5)='Sea ice fraction'
          IF (NFIELDS >= 6) AFLABEL(6)='Lake fraction'
          IF (NFIELDS >= 7) AFLABEL(7)='U-ocean-current '
          IF (NFIELDS >= 8) AFLABEL(8)='V-ocean-current '
        ENDIF

!       COMPUTATION OF THE NORMS OF INPUT FIELDS
!       ----------------------------------------

!       *GLOBAL* ONLY AVAILABLE THE FIRST TIME OR
!       IF REQUESTED BY SETTING LWCOUNORMS=TRUE IN INPUT NAMELIST.
        IF ( LFRSTIME .OR. LWCOUNORMS )THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' *GLOBAL* NORMS OF INPUT ATMOSPHERIC FIELDS :'
          CALL GSTATS(1436,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(IFLD,JF,IC,I)
          DO IFLD=1,NFIELDS
            JF=IFLD
            DO IC=1,3
              NORMS(JF,IC)=FIELDS(1,JF)
            ENDDO
            NORMS(JF,4)=FIELDS(1,JF)**2
            DO I=2,NGPTOTG 
              NORMS(JF,1)=NORMS(JF,1)+FIELDS(I,JF)
              NORMS(JF,2)=MIN(NORMS(JF,2),FIELDS(I,JF))
              NORMS(JF,3)=MAX(NORMS(JF,3),FIELDS(I,JF))
              NORMS(JF,4)=NORMS(JF,4)+FIELDS(I,JF)**2
            ENDDO
            NORMS(JF,1)=NORMS(JF,1)/NGPTOTG
            NORMS(JF,4)=NORMS(JF,4)/NGPTOTG
            NORMS(JF,4)=NORMS(JF,4)-NORMS(JF,1)**2
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1436,1)
          LFRSTIME=.FALSE.
          IF (.NOT.ALLOCATED(ISTART)) ALLOCATE (ISTART(NFIELDS))
          ISTART(:)=0
          N_MASK_IN=NGPTOTG

          DO IFLD=1,NFIELDS
            WRITE(IU06,*) ' ',AFLABEL(IFLD),                            &
     &       NORMS(IFLD,1),NORMS(IFLD,2),NORMS(IFLD,3),N_MASK_IN
            WRITE(IU06,111)                                             &
     &       NORMS(IFLD,1),NORMS(IFLD,2),NORMS(IFLD,3)
          ENDDO

        ELSEIF (LLNORMIFS2WAM) THEN
!       *LOCAL* NORM

          IMASTER=1

          IF (IRANK == IMASTER) THEN
            WRITE(IU06,*) ' NORMS OF GLOBAL INPUT ATMOSPHERIC FIELDS :'
          ELSE
            WRITE(IU06,*) ' *Local* NORMS OF INPUT ATMOSPHERIC FIELDS : '
          ENDIF
            WRITE(IU06,*) ' *Can only be compared if the log files'
            WRITE(IU06,*) '  are for processor: ',IRANK 
            WRITE(IU06,*) '  have the same number of processors: ',NPROC 
          IF (LL1D) THEN
            WRITE(IU06,*) '  and have the same model decomposition: ONE D'
          ELSE
            WRITE(IU06,*) '  and have the same model decomposition: TWO D'
          ENDIF

!         INITIALISE I POINTER TO MASK_IN
          IF (LFRST) THEN
            LFRST=.FALSE.
            JF=1
            N_MASK_IN=0
            DO I=1,NGPTOTG 
              IF (MASK_IN(I) == 1) THEN
                N_MASK_IN=N_MASK_IN+1
              ENDIF
            ENDDO
            ALLOCATE(I_MASK_IN(MAX(1,N_MASK_IN)))
            I_MASK_IN(1)=1

            IC=0
            DO I=1,NGPTOTG 
              IF (MASK_IN(I) == 1) THEN
                IC=IC+1
                I_MASK_IN(IC)=I
              ENDIF
            ENDDO
          ENDIF

!         COMPUTE NORMS
          CALL GSTATS(1436,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(IFLD,JF,VAL,IC)
          DO IFLD=1,NFIELDS
            JF=IFLD
            VAL=FIELDS(I_MASK_IN(1),JF)
            NORMS(JF,1)=VAL
            NORMS(JF,2)=VAL
            NORMS(JF,3)=VAL
            NORMS(JF,4)=VAL**2
            DO IC=2,N_MASK_IN
                VAL=FIELDS(I_MASK_IN(IC),JF)
                NORMS(JF,1)=NORMS(JF,1)+VAL
                NORMS(JF,2)=MIN(NORMS(JF,2),VAL)
                NORMS(JF,3)=MAX(NORMS(JF,3),VAL)
                NORMS(JF,4)=NORMS(JF,4)+VAL**2
            ENDDO
            NORMS(JF,1)=NORMS(JF,1)/MAX(1,N_MASK_IN)
            NORMS(JF,4)=NORMS(JF,4)/MAX(1,N_MASK_IN)
            NORMS(JF,4)=NORMS(JF,4)-NORMS(JF,1)**2
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1436,1)

!         FOR PRIMARY PE (IMASTER), COLLECT ALL THE NORMS TO PRODUCE A
!         PSEUDO GLOBAL NORM.
 
          NCOMLOC=1+4*NFIELDS
          NCOMBUF=NCOMLOC*NPROC
          ALLOCATE(ZCOMBUF(NCOMBUF))
          ALLOCATE(ZCOMBUF1(NCOMLOC))

          ICOUNT=1
          ZCOMBUF1(ICOUNT)=N_MASK_IN
          DO IFLD=1,NFIELDS
            DO IC=1,4
              ICOUNT=ICOUNT+1
              ZCOMBUF1(ICOUNT)=NORMS(IFLD,IC)
            ENDDO
          ENDDO

          CALL GSTATS(647,0)
          IRECVCOUNTS(:)=NCOMLOC

!         CALL MPL_GATHERV(ZCOMBUF1(:),KROOT=IMASTER,
!    &       PRECVBUF=ZCOMBUF,KRECVCOUNTS=IRECVCOUNTS,
!    &       CDSTRING='IFSTOWAM:')

!         Optimise comms to use a 2D gather scheme

          IMASTER=1
          IBEG=(MYSETW-1)*NPRTRV+1
          IEND=IBEG+NPRTRV-1
          ALLOCATE(ZBUFW(SUM(IRECVCOUNTS(IBEG:IEND))))
          CALL MPL_GATHERV(ZCOMBUF1(1:NCOMLOC),KROOT=IMASTER,           &
     &      PRECVBUF=ZBUFW,KRECVCOUNTS=IRECVCOUNTS(IBEG:IEND),          &
     &      KCOMM=MPL_ALL_LEVS_COMM,CDSTRING='IFSTOWAM:')
          CALL PE2SET(IRANK,IA,IB,ISETW,ISETV)
          ITAG=MTAGWAMNORM
          IR=0
          IF ( IRANK == IMASTER )THEN
            ITOT=0
            DO JSETW=1,NPRTRW
              CALL SET2PE(IRECV,0,0,JSETW,1)
              IR=IR+1
              IBEG=(JSETW-1)*NPRTRV+1
              IEND=IBEG+NPRTRV-1
              IINC=SUM(IRECVCOUNTS(IBEG:IEND))
              IBEGOFF=1+ITOT
              IENDOFF=IBEGOFF+IINC-1
              CALL MPL_RECV(ZCOMBUF(IBEGOFF:IENDOFF),                   &
     &          KSOURCE=NPRCIDS(IRECV),                                 &
     &          KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),    &
     &          KTAG=ITAG,CDSTRING='IFSTOWAM:' )
              ITOT=ITOT+IINC
            ENDDO
          ENDIF
          IF ( ISETV == 1 )THEN
            IR=IR+1
            CALL MPL_SEND(ZBUFW,KDEST=NPRCIDS(IMASTER),                 &
     &        KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),      &
     &        KTAG=ITAG,CDSTRING='IFSTOWAM:' )
          ENDIF
          IF (IR > 0) THEN
            CALL MPL_WAIT(KREQUEST=IREQ(1:IR),                          &
     &        CDSTRING='IFSTOWAM: WAIT FOR SENDS AND RECEIVES')
          ENDIF
          DEALLOCATE(ZBUFW)
          CALL GSTATS(647,1)


!         COMPUTE PSEUDO GLOBAL NORM
          IF (IRANK == IMASTER) THEN
            IST=1+(IRANK-1)*NCOMLOC
            ICOUNT=IST
            NMASK=ZCOMBUF(ICOUNT)
            NTOT=NMASK
            DO JF=1,NFIELDS
              ICOUNT=ICOUNT+1
              NORMS(JF,1)=NMASK*ZCOMBUF(ICOUNT)
              ICOUNT=ICOUNT+3
              NORMS(JF,4)=NMASK*ZCOMBUF(ICOUNT)
            ENDDO
            DO IP=1,NPROC
              IF (IP /= IMASTER) THEN
                IST=1+(IP-1)*NCOMLOC
                ICOUNT=IST
                NMASK=ZCOMBUF(ICOUNT)
                NTOT=NTOT+NMASK
                DO JF=1,NFIELDS
                  ICOUNT=ICOUNT+1
                  NORMS(JF,1)=NORMS(JF,1)+NMASK*ZCOMBUF(ICOUNT)
                  ICOUNT=ICOUNT+1
                  NORMS(JF,2)=MIN(NORMS(JF,2),ZCOMBUF(ICOUNT))
                  ICOUNT=ICOUNT+1
                  NORMS(JF,3)=MAX(NORMS(JF,3),ZCOMBUF(ICOUNT))
                  ICOUNT=ICOUNT+1
                  NORMS(JF,4)=NORMS(JF,4)+NMASK*ZCOMBUF(ICOUNT)
                ENDDO
              ENDIF
            ENDDO
            DO JF=1,NFIELDS
              NORMS(JF,1)=NORMS(JF,1)/NTOT
              NORMS(JF,4)=NORMS(JF,4)/NTOT
            ENDDO
          ELSE
            NTOT=N_MASK_IN
          ENDIF

          DO IFLD=1,NFIELDS
            WRITE(IU06,*) ' ',AFLABEL(IFLD),                            &
     &         NORMS(IFLD,1),NORMS(IFLD,2),NORMS(IFLD,3),NTOT,          &
     &                    IRANK, NPROC, LL1D
            WRITE(IU06,111)                                             &
     &         NORMS(IFLD,1),NORMS(IFLD,2),NORMS(IFLD,3)
          ENDDO

          DEALLOCATE(ZCOMBUF)
          DEALLOCATE(ZCOMBUF1)

        ENDIF  ! end production of norms
111     FORMAT(19x,'HEX: ',3(Z16.16,2x))


        IWPER = 1
        ICOORD = 1

        LLNEWCURR=LWCUR

!!!!!!!!! for now make use of lake cover all the time
!!!debile
        LLKC=.true.

        CALL GSTATS(1437,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, KIJL)
        DO ICHNK = 1, NCHNK
          KIJS = 1
          KIJL = NPROMA_WAM 

          CALL FLDINTER (IU06, NGPTOTG, NCA, NRA, NFIELDS, FIELDS,                         &
     &                   NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,                          &
     &                   IFROMIJ(:,ICHNK), JFROMIJ(:,ICHNK), KIJS, KIJL, KIJL4CHNK(ICHNK), &
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,                        &
     &                   ILONRGG, IJBLOCK, ZMISS,                                          &
     &                   LADEN, ROAIR, LGUST, WSTAR0, LLKC, LWCUR,                         &
     &                   LLINTERPOL,                                                       &
     &                   DJ1, DII1, DIIP1, JJ, II, IIP, MASK_IN,                           &
     &                   NXS, NXE, NYS, NYE, FIELDG)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1437,1)

        IPARAMCI=31

      ELSE
        LWCOUSAMEGRID = .FALSE. 
      ENDIF

IF (LHOOK) CALL DR_HOOK('IFSTOWAM',1,ZHOOK_HANDLE)

END ASSOCIATE

END SUBROUTINE IFSTOWAM
