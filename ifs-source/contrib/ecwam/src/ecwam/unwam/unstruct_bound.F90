! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UNSTRUCT_BOUND
  !
  ! Rule of conversion of angle in WAM.
  ! theta_{WAM} = 90 - theta_{trigonometric}
  !
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWSTAT,  ONLY : IREFRA
      USE yowpd, only: MNE=>ne, INE, MNP=>npa, NP_RES => np
      USE yowpd, only: XP=>x, YP=>y, DEP=>z
      USE yowpd, only: exchange
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      IMPLICIT NONE

      REAL(KIND=JWRU), ALLOCATABLE :: SPSIG(:)
      REAL(KIND=JWRU), ALLOCATABLE :: WBAC (:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: WBAC1(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: WBAC2(:,:,:)
      INTEGER(KIND=JWIM) :: nbDirWWM, nbFreqWWM
      REAL(KIND=JWRU), ALLOCATABLE :: SPDIR_WWM(:)
      REAL(KIND=JWRU), ALLOCATABLE :: SPSIG_WWM(:)
      REAL(KIND=JWRU), ALLOCATABLE :: SPDIR_WAM(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: WWM_ID1(:), WWM_ID2(:)
      REAL(KIND=JWRU), ALLOCATABLE :: WWM_WD1(:), WWM_WD2(:)
      INTEGER(KIND=JWIM) :: recTime1, recTime2
      INTEGER(KIND=JWIM), ALLOCATABLE :: Indexes_boundary(:)
      REAL(KIND=JWRU) :: WAV_BoucTime = 0._JWRU
      CHARACTER(LEN=*), PARAMETER :: eFileBnd = 'wwm_bouc_format.nc'
      INTEGER(KIND=JWIM)                :: IWBMNP ! number of wave boundary points
      INTEGER(KIND=JWIM)                :: IWBMNPGL
      INTEGER(KIND=JWIM), ALLOCATABLE   :: IWBNDLC(:) ! local wave boundary index
      INTEGER(KIND=JWIM), ALLOCATABLE   :: IWBNDLC_REV(:) ! local wave boundary index
      INTEGER(KIND=JWIM), ALLOCATABLE   :: IOBPD(:,:) ! boundary direction pointer
      INTEGER(KIND=JWIM), ALLOCATABLE   :: IOBWB(:)   ! gl. wave boundary index ... will vanish in the decomp.
      INTEGER(KIND=JWIM), ALLOCATABLE   :: IOBP(:)    ! boundary points index
      LOGICAL                :: LBCWA = .FALSE.
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BOUNDARY
      IMPLICIT NONE
      ALLOCATE(IOBP(MNP)  ); IOBP = 0
      ALLOCATE(IOBPD(NANG,MNP)); IOBPD = 0
      ALLOCATE(IOBWB(MNP) ); IOBWB = 1 ! for boundary nodes we set this to 0 since in this case we omit advection at these nodes
      CALL SET_IOBP          ! boundary point marker
      CALL SET_IOBPD         ! boundary directional marker
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_NODE(IFILE,NKR,NKG)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     XFN header nodes ...

!     Externals. 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IFILE
      INTEGER(KIND=JWIM), INTENT(OUT):: NKR, NKG

      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*) NKR
      READ(IFILE,*)
      READ(IFILE,*) NKG
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)

      END SUBROUTINE RHEADER_NODE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_ELEMENT(IFILE,NELEM)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     XFN header elements

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IFILE
      INTEGER(KIND=JWIM), INTENT(OUT):: NELEM
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*) NELEM
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      END SUBROUTINE RHEADER_ELEMENT
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_STATUS(STATUS)
      IMPLICIT NONE
      INTEGER(KIND=JWIM), intent(out) :: STATUS(MNP)
      INTEGER(KIND=JWIM) :: COLLECTED(MNP)
      INTEGER(KIND=JWIM) :: PREVVERT(MNP)
      INTEGER(KIND=JWIM) :: NEXTVERT(MNP)
      INTEGER(KIND=JWIM) :: IE, I, IPREV, INEXT, IP
      INTEGER(KIND=JWIM) :: IPPREV, IPNEXT
      INTEGER(KIND=JWIM) :: ZNEXT
      LOGICAL IsFinished
      STATUS(:) = 0
      DO IE=1,MNE
        DO I=1,3
          IF (I.EQ.1) THEN
            IPREV=3
          ELSE
            IPREV=I-1
          END IF
          IF (I.EQ.3) THEN
            INEXT=1
          ELSE
            INEXT=I+1
          END IF
          IP=INE(I,IE)
          IPNEXT=INE(INEXT,IE)
          IPPREV=INE(IPREV,IE)
          IF (STATUS(IP).EQ.0) THEN
            STATUS(IP)=1
            PREVVERT(IP)=IPPREV
            NEXTVERT(IP)=IPNEXT
          END IF
        END DO
      END DO
      STATUS(:)=0
      DO
        COLLECTED(:)=0
        DO IE=1,MNE
          DO I=1,3
            IF (I.EQ.1) THEN
              IPREV=3
            ELSE
              IPREV=I-1
            END IF
            IF (I.EQ.3) THEN
              INEXT=1
            ELSE
              INEXT=I+1
            END IF
            IP=INE(I,IE)
            IPNEXT=INE(INEXT,IE)
            IPPREV=INE(IPREV,IE)
            IF (STATUS(IP).eq.0) THEN
              ZNEXT=NEXTVERT(IP)
              IF (ZNEXT.eq.IPPREV) THEN
                COLLECTED(IP)=1
                NEXTVERT(IP)=IPNEXT
                IF (NEXTVERT(IP).eq.PREVVERT(IP)) THEN
                  STATUS(IP)=1
                END IF
              END IF
            END IF
          END DO
        END DO
        ISFINISHED=.TRUE.
        DO IP=1,MNP
          IF ((COLLECTED(IP).eq.0).and.(STATUS(IP).eq.0)) THEN
            STATUS(IP)=-1
          END IF
          IF (STATUS(IP).eq.0) THEN
            ISFINISHED=.FALSE.
          END IF
        END DO
        IF (ISFINISHED) THEN
          EXIT
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBP
      USE YOWUNPOOL
      USE yowpd, only: ipgl, iplg, np_global
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: I, IWILD(MNP)
      INTEGER(KIND=JWIM) :: I1, I2, I3, IE, IP, ID, IFSTAT, IPglob
      INTEGER(KIND=JWIM) :: ZNEXT, ITMP
      INTEGER(KIND=JWIM) :: IOBPcopy(MNP), eDiff, nbDiff
      INTEGER(KIND=JWIM) :: STATUS(MNP)
      INTEGER(KIND=JWIM) :: istat
      INTEGER(KIND=JWIM), ALLOCATABLE :: Indexes(:)
      REAL(KIND=JWRB) :: BNDTMP
      REAL(KIND=JWRU) :: x1, y1, x2, y2
      REAL(KIND=JWRU) :: EVX, EVY
      REAL(KIND=JWRU) :: eDet1, eDet2
      REAL(KIND=JWRU) :: ATMP, BTMP
      REAL(KIND=JWRU) :: rtemp(MNP)
      LOGICAL :: LFLIVE
      INTEGER(KIND=JWIM) :: IOBPglobal(np_global)
      INTEGER(KIND=JWIM) :: idx
      INTEGER(KIND=JWIM) :: eIOBP
!
! open and read boundary nodes file ...
!
      IOBP    = 0
      IOBPD   = 0


      OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
      CALL RHEADER_NODE(BND%FHNDL,ITMP,ITMP)
      DO IPglob = 1, np_global
        READ(BND%FHNDL,*) ITMP, BNDTMP, BNDTMP, BNDTMP
        eIOBP=INT(BNDTMP)
        IF (.NOT.  LBCWA) THEN
          IF ((eIOBP .EQ. 2).or.(eIOBP .EQ. 3)) eIOBP = 1
        END IF
        IOBPglobal(IPglob) = eIOBP
        IP = ipgl(IPglob)
        IF(IP /= 0) THEN
          IOBP(IP) = eIOBP
        END IF
      END DO
      CLOSE(BND%FHNDL)
      IOBPcopy=IOBP
!
! find islands and domain boundary ....
!
      CALL COMPUTE_STATUS(STATUS)
      DO IP=1,MNP
        IF (STATUS(IP).eq.-1 .AND. IOBP(IP) .EQ. 0) THEN
          IOBP(IP)=1
        END IF
      END DO
!
! reporting differences found
!
! Determining number of boundary nodes
!
      IWBMNP = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2) IWBMNP = IWBMNP + 1 ! Local number of boundary nodes ...
      END DO
!
! map boundary nodes ... needed later for the decomposition ...
!
      ALLOCATE(IWBNDLC(IWBMNP), IWBNDLC_REV(MNP), stat=istat)
      IWBNDLC_REV = 0
      idx = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2) THEN
          idx = idx + 1
          IWBNDLC(idx)    = IP
          IWBNDLC_REV(IP) = idx
        END IF
      END DO
!
! allocate wave boundary arrays ... 
!
      ALLOCATE(Indexes(np_global), stat=istat)
      Indexes = 0
      idx = 0
      DO IPglob = 1, np_global
        IF (IOBPglobal(IPglob) == 2) THEN
          idx = idx + 1
          Indexes(IPglob)=idx
        END IF
      END DO
      IWBMNPGL = idx
      allocate(Indexes_boundary(IWBMNP), stat=istat)
      DO idx=1,IWBMNP
        IP=IWBNDLC(idx)
        IPglob=iplg(IP)
        Indexes_boundary(idx) = Indexes(IPglob)
      END DO
      IF (LBCWA) THEN
        CALL INIT_FILE_BOUNDARY
      END IF      
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Estimate spectral direction that are pointing into the domain from the boundary ...
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL
      USE YOWFRED, ONLY : DFIM, COSTH, SINTH

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: I1, I2, I3, IE, IP, ID
      INTEGER(KIND=JWIM) :: I, IWILD(MNP)
      REAL(KIND=JWRU) :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
      REAL(KIND=JWRU) :: x1, y1, x2, y2
      REAL(KIND=JWRU) :: EVX, EVY
      REAL(KIND=JWRU) :: eDet1, eDet2
      REAL(KIND=JWRU) :: rtemp(MNP)
!
! SET IOBPD ...
!
        DO IE=1,MNE
          I1   =   INE(1,IE)
          I2   =   INE(2,IE)
          I3   =   INE(3,IE)
          DXP1 =   IEN(6,IE)
          DYP1 = - IEN(5,IE)
          DXP2 =   IEN(2,IE)
          DYP2 = - IEN(1,IE)
          DXP3 =   IEN(4,IE)
          DYP3 = - IEN(3,IE)
!2do ... modifly wave direction by currents ...
          DO ID=1,MDC
            EVX=SINTH(ID)
            EVY=COSTH(ID)
            DO I=1,3
              IF (I.eq.1) THEN
                x1=   DXP1
                y1=   DYP1
                x2= - DXP3
                y2= - DYP3
                IP=   I1
              END IF
              IF (I.eq.2) THEN
                x1 =   DXP2
                y1 =   DYP2
                x2 = - DXP1
                y2 = - DYP1
                IP =   I2
              END IF
              IF (I.eq.3) THEN
                x1 =   DXP3
                y1 =   DYP3
                x2 = - DXP2
                y2 = - DYP2
                IP =   I3
              END IF
              eDet1 = SMALL-x1*EVY+y1*EVX
              eDet2 = SMALL+x2*EVY-y2*EVX
              IF ((eDet1.gt.0.0_JWRU).and.(eDet2.gt.0.0_JWRU)) THEN
                IOBPD(ID,IP)=1
              END IF
            END DO
          END DO
        END DO

        DO IP = 1, MNP
          IF ( LBCWA ) THEN
            IF ( IOBP(IP) == 2 ) THEN
              IOBWB(IP) = 0
              IOBPD(:,IP) = 1
            ENDIF
          END IF
          IF ( IOBP(IP) == 3 ) THEN ! If Neumann boundary condition is given set IOBP to 4
            IOBPD(:,IP) = 1 ! Update Neumann nodes ...
          END IF
        END DO

        DO ID=1, MDC
          rtemp = IOBPD(ID,:)
          call exchange(rtemp)
          IOBPD(ID,:) = INT(rtemp)
        END DO
      END SUBROUTINE SET_IOBPD
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_UP_WBAC
      USE YOWSTAT  , ONLY : IDELT
      IMPLICIT NONE
      REAL(KIND=JWRU) :: eTimeDay
      REAL(KIND=JWRU) w1, w2
      INTEGER(KIND=JWIM) :: iTime1, iTime2
      IF (iTime1 .ne. recTime1) THEN
        recTime1 = iTime1
      END IF
      IF (iTime2 .ne. recTime2) THEN
        recTime2 = iTime2
      END IF
      WBAC = w1*WBAC1 + w2*WBAC2
      WAV_BoucTime = WAV_BoucTime + DBLE(IDELT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE APPLY_BOUNDARY_CONDITION(FL)
      USE YOWMPP   , ONLY : NINF, NSUP
      IMPLICIT NONE
      REAL(KIND=JWRB), INTENT(INOUT) :: FL(NINF-1:NSUP,NANG,NFRE)
      INTEGER(KIND=JWIM) :: ID, IS, IP, idx
      DO IS=1,NFRE
        DO ID=1,NANG
          DO idx=1,IWBMNP
            IP=IWBNDLC(idx)
            FL(IP,ID,IS) = REAL(WBAC(ID,IS,idx),JWRB)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FILE_BOUNDARY
      USE YOWPCONS, ONLY : PI, ZPI
      USE YOWFRED, ONLY : FR, TH
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE
