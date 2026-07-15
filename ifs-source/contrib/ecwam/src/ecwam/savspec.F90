! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SAVSPEC(FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 : use PBIO SOFTWARE
!     B. HANSEN    ECMWF      FEBRUARY  1998 : use FDB SOFTWARE

!*    PURPOSE.
!     --------
!     WRITE BINARY SPECTRA TO DISK.

!**   INTERFACE.
!     ----------
!     *CALL* *SAVSPEC(FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)
!     *FL1*       ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *CDTPRO*    CHAR*14   END DATE OF PROPAGATION.
!     *CDATEF*    CHAR*14   END DATE OF ANALYSIS RUN (YYMMDDHHMM).
!     *CDATER*    CHAR*14   DATE FOR OUTPUT OF BOTH RESTART FILES

!     METHOD.
!     -------
!     IN CASE OF MESSAGE PASSING, THE OUTPUT IS DIRECTED STRAIGHT TO ITS
!     PLACE ON DISK, THE CORRESPONDING NAME AND DIRECTORY IS DETERMINED
!     BY GRSTNAME OR BY THE FDB SOFTWARE IF FDB IS USED.
!     THE CONTRIBUTION FROM ALL PE'S NEED TO BE GATHERED ON THE OUTPUT
!     PE BY CALLING MPGATHERFL.
!     IN ORDER TO USE MPL_GATHERV IN MPGATHERFL THE SPECTRA WILL BE
!     GATHERED MDEL FREQUENCIES AT A TIME
!     AND KDEL DIRECTIONS AT A TIME.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,IPFGTBL  ,KDEL    ,MDEL      ,    &
     &                      LRSTPARALW
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET,        &
     &                      NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWMAP   , ONLY : NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "expand_string.intfb.h"
#include "grstname.intfb.h"
#include "mpgatherfl.intfb.h"
#include "writefl.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE
      CHARACTER(LEN=14), INTENT(IN) :: CDTPRO, CDATEF, CDATER


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: IJSG, IJLG, IJSB, IJLB, KIJS, KIJL, IPRM, ICHNK
      INTEGER(KIND=JWIM) :: IUNIT
      INTEGER(KIND=JWIM) :: IFCST
      INTEGER(KIND=JWIM) :: IRECV, MLOOP, KLOOP, MINF, MSUP, KINF, KSUP 
      INTEGER(KIND=JWIM) :: LNAME

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:) :: RFL

      CHARACTER(LEN=296) :: FILENAME

      LOGICAL :: LOUNIT

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SAVSPEC',0,ZHOOK_HANDLE)

      LOUNIT = .TRUE.

      IFCST = 0
      CALL GRSTNAME(CDTPRO,CDATEF,IFCST,'BLS',ICPLEN,CPATH,FILENAME)

      IF (LRSTPARALW) THEN
!        RESTART FILES FROM ALL PS's
         LNAME = LEN_TRIM(FILENAME)
         FILENAME=FILENAME(1:LNAME)//'.%p_%n'
         CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAME,1)

         IJSG = IJFROMCHNK(1,1)
         IJLG = IJSG + SUM(KIJL4CHNK) - 1
         ALLOCATE(RFL(IJSG:IJLG, NANG, NFRE))
!$OMP    PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IJSB, IJLB)
         DO ICHNK = 1, NCHNK
           KIJS = 1
           IJSB = IJFROMCHNK(KIJS, ICHNK)
           KIJL = KIJL4CHNK(ICHNK)
           IJLB = IJFROMCHNK(KIJL, ICHNK)

           RFL(IJSB:IJLB, :, :) = FL1(KIJS:KIJL, :, :, ICHNK)
         ENDDO
!$OMP    END PARALLEL DO

         CALL WRITEFL(RFL, IJSG, IJLG, 1, NANG, 1, NFRE,               &
     &                FILENAME, IUNIT, LOUNIT, LRSTPARALW)

         DEALLOCATE(RFL)
      ELSE

!       COLLECT THE DIFFERENT CONTRIBUTIONS AND PROCEED TO THE OUTPUT OF
!       SPECTRUM FROM PE IPFGTBL(JPPFLAG+1)
        IRECV=IPFGTBL(JPPFLAG+1)

        DO MLOOP= 1, NFRE, MDEL
          MINF=MLOOP
          MSUP=MIN(MLOOP+MDEL-1,NFRE)
          DO KLOOP = 1, NANG, KDEL
            KINF=KLOOP
            KSUP=MIN(KLOOP+KDEL-1,NANG)

            LOUNIT = .FALSE.
            IF (MINF == 1 .AND. KINF == 1) LOUNIT = .TRUE.

            ALLOCATE(RFL(NIBLO, KINF:KSUP, MINF:MSUP))

!$OMP       PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, M, K, IPRM, IJ)
            DO ICHNK = 1, NCHNK
              DO M = MINF, MSUP
                DO K = KINF, KSUP 
                  DO IPRM = 1, KIJL4CHNK(ICHNK)
                    IJ = IJFROMCHNK(IPRM,ICHNK)
                    IF (IJ >= IJSLOC .AND. IJ <= IJLLOC) THEN
                      IJ = IJ + IJGLOBAL_OFFSET 
                      RFL(IJ, K, M) = FL1(IPRM, K, M, ICHNK)
                    ENDIF
                  ENDDO
                 ENDDO
               ENDDO
             ENDDO
!$OMP        END PARALLEL DO

            CALL MPGATHERFL(IRECV, NBLKS, NBLKE, KINF, KSUP, MINF, MSUP, RFL)


            IF (IRANK == IPFGTBL(JPPFLAG+1) .OR. NPROC == 1) THEN
              CALL WRITEFL(RFL, 1, NIBLO, KINF, KSUP, MINF, MSUP,       &
     &                     FILENAME, IUNIT, LOUNIT, LRSTPARALW)
            ENDIF
            DEALLOCATE(RFL)
          ENDDO
        ENDDO
      ENDIF

      WRITE(IU06,*) ' SPECTRUM FILE DISPOSED AT........ CDTPRO  = ', CDTPRO

      IF (LHOOK) CALL DR_HOOK('SAVSPEC',1,ZHOOK_HANDLE)

      END SUBROUTINE SAVSPEC
