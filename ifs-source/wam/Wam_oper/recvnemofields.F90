! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RECVNEMOFIELDS(BLK2LOC, WVENVI, NEMO2WAM,  &
 &                        NXS, NXE, NYS, NYE, FIELDG, &
 &                        FF_NOW, LREST, LINIT) 

!****  *RECVNEMOFIELDS* - UPDATE FIELDS WAVE FIELDS WITH NEMO INFORMATION

!      KRISTIAN MOGENSEN ECMWF    MARCH 2013

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES NEMO INFORMATION THROUGH TO
!          WAM VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!*    INTERFACE.                                                        
!     ----------                                                        


!     METHOD.                                                           
!     -------                                                           

!          PARALLEL INTERPOLATION BASED ON PREDETERMINED WEIGHTS

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_GET  -  UPDATE NEMO FIELDS IN WAM

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

! GRID POINTS CHUNKS
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, NTOTIJ, KIJL4CHNK
! MODULES NEEDED FOR LAKE MASK HANDLING
      USE YOWWIND  , ONLY : LLNEWCURR 
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_MODULE, ONLY : MPL_COMM
! COUPLING INFORMATION
      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOUCIC, LWNEMOCOUCIT, LWNEMOCOUCUR
! ICE AND CURRENT INFORMATION 
      USE YOWCURR  , ONLY : CURRENT_MAX
! OUTPUT FORTRAN UNIT
      USE YOWTEST  , ONLY : IU06
! NEMO FIELDS ON WAVE GRID
      USE YOWNEMOFLDS,ONLY: LNEMOCITHICK, LNEMOICEREST  
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! -------------------------------------------------------------------   

      IMPLICIT NONE

#include "abort1.intfb.h"

      TYPE(WVGRIDLOC), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: WVENVI
      TYPE(OCEAN2WAVE), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: NEMO2WAM
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), DIMENSION(NXS:NXE, NYS:NYE), INTENT(IN) :: FIELDG

      TYPE(FORCING_FIELDS), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: FF_NOW ! FORCING FIELDS
      LOGICAL, INTENT(IN) :: LREST ! RESTART SO UPDATE FROM RESTART VALUES
      LOGICAL, INTENT(IN) :: LINIT ! UPDATE CICOVER, CITHICK, UCUR, VCUR AT INITIAL TIME
                                   ! IF NEEDED.


      INTEGER(KIND=JWIM), PARAMETER :: NFIELD = 5
      INTEGER(KIND=JWIM) :: IX, JY, IJ
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IC, IFLD
      REAL(KIND=JWRO), DIMENSION(NTOTIJ, NFIELD) :: ZNEMOTOWAM
       REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LLFLDUPDT
! -------------------------------------------------------------------   

IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',0,ZHOOK_HANDLE)

ASSOCIATE(IFROMIJ => BLK2LOC%IFROMIJ, &
 &        JFROMIJ => BLK2LOC%JFROMIJ, &
 &        UCUR => WVENVI%UCUR, &
 &        VCUR => WVENVI%VCUR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        CITHICK => FF_NOW%CITHICK, &
 &        NEMOSST => NEMO2WAM%NEMOSST, &
 &        NEMOCICOVER => NEMO2WAM%NEMOCICOVER, &
 &        NEMOCITHICK => NEMO2WAM%NEMOCITHICK, &
 &        NEMOUCUR => NEMO2WAM%NEMOUCUR, &
 &        NEMOVCUR => NEMO2WAM%NEMOVCUR)


      ! IF WE ARE IN RESTART JUST COPY THE RESTART FILE INFO.
      IF (LREST) THEN
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)
        DO ICHNK = 1, NCHNK
          NEMOCICOVER(:, ICHNK)=CICOVER(:, ICHNK)
          NEMOCITHICK(:, ICHNK)=CITHICK(:, ICHNK)
          NEMOUCUR(:, ICHNK)=UCUR(:, ICHNK)
          NEMOVCUR(:, ICHNK)=VCUR(:, ICHNK)
        ENDDO
!$OMP   END PARALLEL DO

        LNEMOICEREST=.TRUE.
        LNEMOCITHICK=.TRUE.
      ELSE
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_GET( IRANK-1, NPROC, MPL_COMM,     &
     &                            NTOTIJ, NFIELD, ZNEMOTOWAM,   &
     &                            LNEMOCITHICK, LLFLDUPDT )

        LLNEWCURR=.TRUE.

        IF (LLFLDUPDT) THEN
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IC, IFLD)
          DO ICHNK = 1, NCHNK

            KIJS = 1
            DO IC = 1, ICHNK-1
               KIJS = KIJS + KIJL4CHNK(IC)
            ENDDO
            KIJL = KIJS + KIJL4CHNK(ICHNK) - 1

            IFLD = 1
            NEMOSST(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOCICOVER(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOCITHICK(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOUCUR(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOVCUR(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)

            IF ( KIJL4CHNK(ICHNK) < NPROMA_WAM ) THEN
!             values for fictious points
              NEMOSST(KIJL4CHNK(ICHNK)+1:NPROMA_WAM, ICHNK)     = NEMOSST(1, ICHNK)
              NEMOCICOVER(KIJL4CHNK(ICHNK)+1:NPROMA_WAM, ICHNK) = NEMOCICOVER(1, ICHNK)
              NEMOCITHICK(KIJL4CHNK(ICHNK)+1:NPROMA_WAM, ICHNK) = NEMOCITHICK(1, ICHNK)
              NEMOUCUR(KIJL4CHNK(ICHNK)+1:NPROMA_WAM, ICHNK)    = NEMOUCUR(1, ICHNK)
              NEMOVCUR(KIJL4CHNK(ICHNK)+1:NPROMA_WAM, ICHNK)    = NEMOVCUR(1, ICHNK)
            ENDIF
          ENDDO
!$OMP     END PARALLEL DO
        ENDIF

#endif
        LNEMOICEREST=.FALSE.

      ENDIF

!     UPDATE CICOVER, CITHICK UCUR AND VCUR AT INITIAL TIME ONLY !!!!
      IF (LINIT) THEN

        WRITE(IU06,*)' RECVNEMOFIELDS: INITIALISE OCEAN FIELDS'

        IF (LWCOU) THEN
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
          DO ICHNK = 1, NCHNK

           IF (LWNEMOCOUCIC) THEN
              DO IJ = 1, NPROMA_WAM
                IX = IFROMIJ(IJ,ICHNK)
                JY = JFROMIJ(IJ,ICHNK)
!              if lake cover = 0, we assume open ocean point, then get sea ice directly from NEMO
                IF (FIELDG(IX,JY)%LKFR <= 0.0_JWRB ) THEN
                  CICOVER(IJ,ICHNK)=NEMOCICOVER(IJ,ICHNK)
                ELSE
!              get ice information from atmopsheric model
                  CICOVER(IJ,ICHNK)=FIELDG(IX,JY)%CICOVER 
                ENDIF
              ENDDO
            ENDIF

            IF (LWNEMOCOUCIT) THEN
              DO IJ = 1, NPROMA_WAM
                IX = IFROMIJ(IJ, ICHNK)
                JY = JFROMIJ(IJ, ICHNK)
!              if lake cover = 0, we assume open ocean point, then get sea ice thickness directly from NEMO
                IF (FIELDG(IX,JY)%LKFR <= 0.0_JWRB ) THEN
                  CITHICK(IJ,ICHNK)=NEMOCICOVER(IJ,ICHNK)*NEMOCITHICK(IJ,ICHNK)
                ELSE
                  CICOVER(IJ,ICHNK)=0.5_JWRB*NEMOCICOVER(IJ,ICHNK)
                ENDIF
              ENDDO
            ENDIF

            IF (LWNEMOCOUCUR) THEN
              DO IJ = 1, NPROMA_WAM
                IX = IFROMIJ(IJ, ICHNK)
                JY = JFROMIJ(IJ, ICHNK)
!              if lake cover = 0, we assume open ocean point, then get currents directly from NEMO
                IF (FIELDG(IX,JY)%LKFR <= 0.0_JWRB ) THEN
                  UCUR(IJ,ICHNK) = SIGN(MIN(ABS(NEMOUCUR(IJ,ICHNK)),CURRENT_MAX),NEMOUCUR(IJ,ICHNK))
                  VCUR(IJ,ICHNK) = SIGN(MIN(ABS(NEMOVCUR(IJ,ICHNK)),CURRENT_MAX),NEMOVCUR(IJ,ICHNK))
                ELSE
                  UCUR(IJ,ICHNK)=0.0_JWRB
                  VCUR(IJ,ICHNK)=0.0_JWRB
                ENDIF
              ENDDO
            ENDIF

          ENDDO
!$OMP   END PARALLEL DO

        ELSE

!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)
          DO ICHNK = 1, NCHNK
            IF (LWNEMOCOUCIC) CICOVER(:,ICHNK)=NEMOCICOVER(:,ICHNK)
            IF (LWNEMOCOUCIT) CITHICK(:,ICHNK)=NEMOCICOVER(:,ICHNK)*NEMOCITHICK(:,ICHNK)
            IF (LWNEMOCOUCUR) THEN
             UCUR(:,ICHNK)=NEMOUCUR(:,ICHNK)
             VCUR(:,ICHNK)=NEMOVCUR(:,ICHNK)
            ENDIF
          ENDDO
!$OMP     END PARALLEL DO

        ENDIF

      ENDIF

      IF (LWNEMOCOUCIT.AND.(.NOT.LNEMOCITHICK)) THEN
        WRITE(IU06,*) ' --------------------------------'
        WRITE(IU06,*) ' LWNEMOCOUCIT ONLY MAKES SENSES  '
        WRITE(IU06,*) ' IF LIM IS ACTIVATED IN NEMO     '
        WRITE(IU06,*) ' --------------------------------'
        CALL ABORT1
      ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',1,ZHOOK_HANDLE)

END SUBROUTINE RECVNEMOFIELDS
