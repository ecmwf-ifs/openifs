! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MAKEGRID (BLOCK, GRID, PMISS)

! ----------------------------------------------------------------------

!**** *MAKEGRID* - MAKE GRIDDED WAM MODEL FIELDS FROM BLOCKED FIELDS.

!     H. GUNTHER       ECMWF    NOVEMBER 1989

!*    PURPOSE.
!     --------

!       GRIDDED WAVE FIELDS ARE CREATED FROM BLOCKED WAVE FIELDS.

!**   INTERFACE.
!     ----------

!        *CALL MAKEGRID*
!           *BLOCK*   REAL   DATA IN BLOCKED FORMAT
!           *GRID*    REAL   DATA IN GRID FORMAT
!           *PMISS*   REAL  VALUE GIVEN FOR ALL NON SEA POINTS. 

!     METHOD.
!     -------

!       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
!       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
!       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
!       INITIALISED WITH PMISS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : BLK2GLO  ,AMOWEP   ,AMONOP   ,    &
     &            XDELLA   ,ZDELLO   ,IPER     ,              &
     &            NGX      ,NGY      ,NIBLO
      USE YOWPARAM , ONLY : LLUNSTR
#ifdef WAM_HAVE_UNWAM
      USE YOWUNPOOL, ONLY : OUT_METHOD, IE_OUTPTS
      USE OUTPUT_STRUCT, ONLY : IXarr, IYarr
      USE YOWPD,     ONLY : NODES=>nodes_global,INE_GLOBAL
      USE OUTPUT_STRUCT, ONLY : INTELEMENT_IPOL
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif
      USE YOWABORT  ,ONLY : WAM_ABORT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(INOUT) :: BLOCK(:)
      REAL(KIND=JWRB), INTENT(IN) :: PMISS
      REAL(KIND=JWRB), INTENT(OUT) :: GRID(NGX,NGY)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRU) :: PMS8
      INTEGER(KIND=JWIM) :: I, J

      INTEGER(KIND=JWIM) :: IY, IX, IJ
#ifdef WAM_HAVE_UNWAM
      INTEGER(KIND=JWIM) :: IE, KI, IR, IP
      INTEGER(KIND=JWIM) :: NI(3)

      REAL(KIND=JWRB), ALLOCATABLE :: BLOCK_G(:)

      REAL(KIND=JWRU) :: XLO, XLA
      REAL(KIND=JWRU) :: GRID8
      REAL(KIND=JWRU) :: XYELE(2,3), SKALAR(3)

      LOGICAL :: LLCLST

      INTEGER(KIND=JWIM) :: iNode
      INTEGER(KIND=JWIM) :: NIBLO_OUT
#endif

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MAKEGRID',0,ZHOOK_HANDLE)

      CALL GSTATS(1996,0)

!*    1. INITIALIZATION
!     -----------------

      PMS8=PMISS

      DO J = 1,NGY
        DO I = 1,NGX
          GRID(I,J) = PMISS
        ENDDO
      ENDDO

!*    2. MAKE GRIDDED FIELD
!     ---------------------

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        NIBLO_OUT=SIZE(BLOCK)
        IF (OUT_METHOD == 1) THEN
          ALLOCATE(BLOCK_G(NIBLO_OUT))
          CALL UNBLKRORD(1, 1, NIBLO_OUT, 1, 1, 1, 1, BLOCK, BLOCK_G)
!!! I believe one could have an openmp loop here !!!!
          DO IY=1,NGY
            XLA = REAL(AMONOP-REAL(IY-1,JWRB)*XDELLA,JWRU)
            DO IX=1,NGX
              XLO = REAL(AMOWEP+REAL(IX-1,JWRB)*ZDELLO(IY),JWRU)
!!debile
!!! need to find out unstructured grid left longitude (I assume -180 for now)
              IF (XLO >= 180.0_JWRU) XLO=MAX(XLO-IPER*360.0_JWRU,-180.0_JWRU)

              IE=IE_OUTPTS(IX,IY)
              IF (IE /= -1) THEN
!!! this is not very efficient
                LLCLST=.FALSE.
                NI = INE_GLOBAL(:,IE)
                XYELE(1,:)=NODES(NI)%X
                XYELE(2,:)=NODES(NI)%Y
                SKALAR(:)=BLOCK_G(NI)
                IF (ANY(BLOCK_G(NI) == PMISS)) LLCLST=.TRUE.
                CALL INTELEMENT_IPOL(LLCLST, XYELE, SKALAR, XLO, XLA, IE, PMS8, GRID8)
                GRID(IX,IY)=GRID8
              ENDIF
            ENDDO
          ENDDO
          DEALLOCATE(BLOCK_G)
        ENDIF
        IF (OUT_METHOD == 2) THEN
          DO iNode=1,NIBLO_OUT
            IX=IXarr(iNode)
            IY=IYarr(iNode)
            GRID(IX,IY) = BLOCK(iNode)
          ENDDO
        ENDIF
#else
      CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSE
        DO IJ = 1, NIBLO 
          IX = BLK2GLO%IXLG(IJ)
          IY = NGY-BLK2GLO%KXLT(IJ)+1
          GRID(IX,IY) = BLOCK(IJ)
        ENDDO
      ENDIF ! LLUNSTR

      CALL GSTATS(1996,1)

      IF (LHOOK) CALL DR_HOOK('MAKEGRID',1,ZHOOK_HANDLE)

      END SUBROUTINE MAKEGRID
