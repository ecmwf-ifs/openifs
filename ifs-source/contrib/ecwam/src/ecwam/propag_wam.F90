! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PROPAG_WAM (BLK2GLO, WVENVI, WVPRPT, FL1)

! ----------------------------------------------------------------------

!**** *PROPAG_WAM* - WAVE PROPGATION

!*    PURPOSE.
!     --------

!     PROPAGATION

!**   INTERFACE.
!     ----------

!     *CALL* *PROPAG_WAM (BLK2GLO, WVENVI, WVPRPT, FL1)
!          *BLK2GLO*   - BLOCK TO GRID TRANSFORMATION
!          *WVENVI*    - WAVE ENVIRONMENT FIELDS
!          *WVPRPT*    - WAVE PROPERTIES FIELDS
!          *FL1*       - SPECTRUM


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY

      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWMAP   , ONLY : NIBLO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,LLUNSTR
      USE YOWREFD  , ONLY : LLUPDTTD ,THDD     ,THDC     ,SDOT
      USE YOWSTAT  , ONLY : IPROPAGS ,IFRELFMAX, DELPRO_LF, IDELPRO
      USE YOWUBUF  , ONLY : LUPDTWGHT
#ifdef WAM_HAVE_UNWAM
      USE UNWAM    , ONLY : PROPAG_UNWAM
#endif

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "ctuwupdt.intfb.h"
#include "mpexchng.intfb.h"
#include "proenvhalo.intfb.h"
#include "propags.intfb.h"
#include "propags1.intfb.h"
#include "propags2.intfb.h"
#include "propdot.intfb.h"

      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(ENVIRONMENT), INTENT(IN) :: WVENVI
      TYPE(FREQUENCY), INTENT(IN) :: WVPRPT
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M, J
      INTEGER(KIND=JWIM) :: JKGLO, NPROMA, MTHREADS
      INTEGER(KIND=JWIM) :: NSTEP_LF, ISUBST
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS
      INTEGER(KIND=JWIM) :: IJSG, IJLG, ICHNK, KIJS, KIJL, IJSB, IJLB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     Spectra extended with the halo exchange for the propagation
!     But limited to NFRE_RED frequencies
!!! the advection schemes are still written in block structure
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NANG, NFRE_RED) :: FL1_EXT, FL3_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED) :: WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1) :: DELLAM1_EXT, COSPHM1_EXT 
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1) :: DEPTH_EXT, UCUR_EXT, VCUR_EXT

      LOGICAL :: L1STCALL

      DATA L1STCALL / .TRUE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',0,ZHOOK_HANDLE)


      IF (NIBLO > 1) THEN

        IJSG = IJFROMCHNK(1,1)
        IJLG = IJSG + SUM(KIJL4CHNK) - 1

        MTHREADS=1
!$      MTHREADS=OMP_GET_MAX_THREADS()
        NPROMA=(IJLG-IJSG+1)/MTHREADS + 1


!!! the advection schemes are still written in block structure
!!! mapping chuncks to block ONLY for actual grid points !!!!
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, M, K)
        DO ICHNK = 1, NCHNK
          KIJS = 1
          IJSB = IJFROMCHNK(KIJS, ICHNK)
          KIJL = KIJL4CHNK(ICHNK)
          IJLB = IJFROMCHNK(KIJL, ICHNK)
          DO M = 1, NFRE_RED
            DO K = 1, NANG
              FL1_EXT(IJSB:IJLB, K, M) = FL1(KIJS:KIJL, K, M, ICHNK)
            ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

!       SET THE DUMMY LAND POINT TO 0.
        FL1_EXT(NSUP+1,:,:) = 0.0_JWRB 


        IF (LLUNSTR) THEN 

        WRITE(NULERR,*) '!!! ********************************* !!'
        WRITE(NULERR,*) '!!! in PROPAG_WAM Not yet ready !!!' 
        WRITE(NULERR,*) '!!! PROPAG_UNWAM will need to be adapted to new FL1 structure !!!' 
        WRITE(NULERR,*) '!!! ********************************* !!'
        CALL ABORT1

!!!!!           CALL PROPAG_UNWAM(FL1_EXT, FL1)

        ELSE


!          OBTAIN INFORMATION AT NEIGHBORING GRID POINTS (HALO)
!          ----------------------------------------------------
           CALL MPEXCHNG(FL1_EXT, NANG, 1, NFRE_RED)


           CALL GSTATS(1430,0)

!          UPDATE DOT THETA TERM
!          ---------------------

           IF (LLUPDTTD) THEN
             IF (.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(IJSG:IJLG, NANG))
             IF (.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(IJSG:IJLG, NANG))
             IF (.NOT.ALLOCATED(SDOT)) ALLOCATE(SDOT(IJSG:IJLG, NANG, NFRE_RED))

!            NEED HALO VALUES
             CALL  PROENVHALO (NINF, NSUP,                            &
&                              WVPRPT,                                &
&                              WVENVI,                                &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DELLAM1_EXT, COSPHM1_EXT,              & 
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )


!            DOT THETA TERM:

!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPDOT(KIJS, KIJL, NINF, NSUP,                                     &
     &                      BLK2GLO,                                                    &
     &                      WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT,                      &
     &                      COSPHM1_EXT, DEPTH_EXT, UCUR_EXT, VCUR_EXT,                 &
     &                      THDC(KIJS:KIJL,:), THDD(KIJS:KIJL,:), SDOT(KIJS:KIJL,:,:))
             ENDDO
!$OMP        END PARALLEL DO

             LLUPDTTD = .FALSE.
           ENDIF


!          IPROPAGS = 2 is the default option (the other options are kept but usually not used)
!          -----------------------------------------------------------------------------------
           SELECT CASE (IPROPAGS)

           CASE(2)

             IF (LUPDTWGHT) THEN
!              NEED HALO VALUES
               CALL  PROENVHALO (NINF, NSUP,                            &
&                                WVPRPT,                                &
&                                WVENVI,                                &
&                                WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                                DELLAM1_EXT, COSPHM1_EXT,              & 
&                                DEPTH_EXT, UCUR_EXT, VCUR_EXT )

!              COMPUTES ADVECTION WEIGTHS AND CHECK CFL CRITERIA
               CALL CTUWUPDT(IJSG, IJLG, NINF, NSUP,                      &
&                            BLK2GLO,                                     &
&                            CGROUP_EXT, OMOSNH2KD_EXT,                   &
&                            COSPHM1_EXT, DEPTH_EXT, UCUR_EXT, VCUR_EXT )

               LUPDTWGHT=.FALSE.
             ENDIF

!$OMP        PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS2(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL, NANG, 1, NFRE_RED)
             ENDDO
!$OMP        END PARALLEL DO


!            SUB TIME STEPPING FOR FAST WAVES (only if IFRELFMAX > 0)
             IF (IFRELFMAX > 0 ) THEN
               NSTEP_LF = NINT(REAL(IDELPRO, JWRB)/DELPRO_LF)
               ISUBST = 2  ! The first step was done as part of the previous call to PROPAGS2

               DO WHILE (ISUBST <= NSTEP_LF)

!$OMP            PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL, M, K, IJ)
                 DO JKGLO = IJSG, IJLG, NPROMA
                   KIJS=JKGLO
                   KIJL=MIN(KIJS+NPROMA-1, IJLG)
                   DO M = 1, IFRELFMAX 
                     DO K = 1, NANG
                       DO IJ = KIJS, KIJL
                         FL1_EXT(IJ, K, M) = FL3_EXT(IJ, K, M)
                       ENDDO
                     ENDDO
                   ENDDO
                 ENDDO
!$OMP            END PARALLEL DO

                 CALL MPEXCHNG(FL1_EXT(:,:,1:IFRELFMAX), NANG, 1, IFRELFMAX)

!$OMP            PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL)
                 DO JKGLO = IJSG, IJLG, NPROMA
                   KIJS=JKGLO
                   KIJL=MIN(KIJS+NPROMA-1, IJLG)
                   CALL PROPAGS2(FL1_EXT(:,:,1:IFRELFMAX), FL3_EXT(:,:,1:IFRELFMAX), NINF, NSUP, KIJS, KIJL, NANG, 1, IFRELFMAX)
                 ENDDO
!$OMP            END PARALLEL DO

                 ISUBST = ISUBST + 1

               ENDDO
             ENDIF  ! end sub time steps (if needed)

           CASE(1)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO (NINF, NSUP,                            &
&                              WVPRPT,                                &
&                              WVENVI,                                &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DELLAM1_EXT, COSPHM1_EXT,              & 
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )

!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS1(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL,       &
&                            BLK2GLO,                                        &
&                            DEPTH_EXT,                                      &
&                            CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                            DELLAM1_EXT, COSPHM1_EXT,                       &
&                            UCUR_EXT, VCUR_EXT,                             &
&                            L1STCALL)
             ENDDO
!$OMP       END PARALLEL DO

           CASE(0)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO (NINF, NSUP,                            &
&                              WVPRPT,                                &
&                              WVENVI,                                &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DELLAM1_EXT, COSPHM1_EXT,              & 
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )

!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL,       &
&                           BLK2GLO,                                        &
&                           DEPTH_EXT,                                      &
&                           CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                           DELLAM1_EXT, COSPHM1_EXT,                       &
&                           UCUR_EXT, VCUR_EXT,                             &
&                           L1STCALL)
             ENDDO
!$OMP        END PARALLEL DO
           END SELECT 


!!! the advection schemes are still written in block structure
!!!  So need to convert back to the nproma_wam chuncks
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, M, K)
          DO ICHNK = 1, NCHNK
            KIJS = 1
            IJSB = IJFROMCHNK(KIJS, ICHNK)
            KIJL = KIJL4CHNK(ICHNK)
            IJLB = IJFROMCHNK(KIJL, ICHNK)
            DO M = 1, NFRE_RED
              DO K = 1, NANG
                FL1(KIJS:KIJL, K, M, ICHNK) = FL3_EXT(IJSB:IJLB, K, M)
              ENDDO
            ENDDO

            IF (KIJL < NPROMA_WAM) THEN
              !!! make sure fictious points keep values of the first point in the chunk
              DO M = 1, NFRE_RED
                DO K = 1, NANG
                  FL1(KIJL+1:NPROMA_WAM, K, M, ICHNK) = FL1(1, K, M, ICHNK)
                ENDDO
              ENDDO
            ENDIF

          ENDDO
!$OMP     END PARALLEL DO


           CALL GSTATS(1430,1)

        ENDIF  ! end propagation

      ENDIF ! more than one grid point

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.

IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAG_WAM
