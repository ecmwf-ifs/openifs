! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PROPCONNECT(IJS, IJL, NEWIJ2IJ)
! ----------------------------------------------------------------------

!**** *PROPCONNECT* - ROUTINE TO DETERMINE THE NEIGHBOURING POINTS FOR THE STRUCTURED GRID ADVECTION

!*    PURPOSE.
!     -------

!     TO ARRANGE NEIGHBOUR GRID POINT INDICES

!**   INTERFACE.
!     ----------

!       *CALL* *PROPCONNECT*

!     METHOD.
!     -------

!       THE INDICES OF THE NEXT POINTS ON LAT. AND LONG. ARE
!       COMPUTED. ZERO INDICATES A LAND POINT IS NEIGHBOUR.
!       THE FINAL COMMON UBUF IS WRITTEN OUT.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : BLK2GLO  ,IPER     ,IRGG      ,   &
     &                      XDELLA   ,ZDELLO   ,NLONRGG   ,   &
     &                      NGX      ,NGY      ,NIBLO
      USE YOWSHAL  , ONLY : LLOCEANMASK 
      USE YOWSTAT  , ONLY : IPROPAGS 
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR     ,              &
     &                      KRLAT    ,KRLON    ,                        &
     &                      WLAT     ,WCOR     ,WRLAT    ,WRLON

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: NEWIJ2IJ

      INTEGER(KIND=JWIM) :: IJ, IJP, I, K, IP, IP1D, IH, ID1D, IS, M
      INTEGER(KIND=JWIM) :: IMIN, IPLUS, IMIN2, IPLUS2
      INTEGER(KIND=JWIM) :: IC, ICP, ICL, ICR

      REAL(KIND=JWRB) :: XMIN, XPLUS, D0, D1, D2, D3, D4, D5, D6
      REAL(KIND=JWRB) :: XLON, XL, XP, TWOXDELLA
      REAL(KIND=JWRB) :: XLL, XPL, XLR, XPR


! ----------------------------------------------------------------------

      TWOXDELLA=2.0_JWRB*XDELLA

!*    2. COMPUTE INDICES OF NEIGHBOUR SEA POINTS.
!        ----------------------------------------

!*    2.1 LATITUDE NEIGHBOURS (KLAT)
!         --------------------------

      DO ICL=1,2
        DO IJP=1,2
          DO IJ=IJS, IJL
             KLAT(IJ,IJP,ICL) = 0
          ENDDO
        ENDDO
      ENDDO

      DO IP = IJS, IJL
        I = BLK2GLO%IXLG(IP)
        K = BLK2GLO%KXLT(IP)
        IP1D = NEWIJ2IJ(IP)
        IF (K > 1) THEN
          XMIN = REAL(I-1,JWRB)*ZDELLO(K)/ZDELLO(K-1)
          IMIN = NINT(XMIN) + 1

!         CLOSEST GRID POINT
          IF (LLOCEANMASK(IMIN,K-1)) THEN
            DO ID1D = IP1D,1,-1
              IH = IJ2NEWIJ(ID1D)
              IF (BLK2GLO%IXLG(IH) == IMIN .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
            ENDDO
            KLAT(IP,1,1) = IJ2NEWIJ(ID1D) 
          ENDIF

!         SECOND CLOSEST GRID POINT
          IF (IRGG == 1) THEN
            IF (XMIN <= REAL(IMIN-1,JWRB)) THEN
              IF (IMIN <= 1) THEN
                IMIN2=1
              ELSE
                IMIN2=IMIN-1
              ENDIF
            ELSE
              IF (IMIN >= NLONRGG(K-1)) THEN
                IMIN2=NLONRGG(K-1)
              ELSE
                IMIN2=IMIN+1
              ENDIF
            ENDIF

            IF (LLOCEANMASK(IMIN2,K-1)) THEN
              DO ID1D = IP1D,1,-1
                IH = IJ2NEWIJ(ID1D)
                IF (BLK2GLO%IXLG(IH) == IMIN2 .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
              ENDDO
              KLAT(IP,1,2) = IJ2NEWIJ(ID1D)
            ENDIF
          ELSE
            KLAT(IP,1,2) = KLAT(IP,1,1)
          ENDIF

        ENDIF

        IF (K < NGY) THEN
          XPLUS = REAL(I-1,JWRB)*ZDELLO(K)/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (LLOCEANMASK(IPLUS,K+1)) THEN
            DO ID1D = IP1D,NIBLO
              IH = IJ2NEWIJ(ID1D)
              IF (BLK2GLO%IXLG(IH) == IPLUS .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
            ENDDO
            KLAT(IP,2,1) = IJ2NEWIJ(ID1D)
          ENDIF

          IF (IRGG == 1) THEN
            IF (XPLUS <= REAL(IPLUS-1,JWRB)) THEN
              IF (IPLUS <= 1) THEN
                IPLUS2=1
              ELSE
                IPLUS2=IPLUS-1
              ENDIF
            ELSE
              IF (IPLUS >= NLONRGG(K+1)) THEN
                IPLUS2=NLONRGG(K+1)
              ELSE
                IPLUS2=IPLUS+1
              ENDIF
            ENDIF

            IF (LLOCEANMASK(IPLUS2,K+1)) THEN
              DO ID1D = IP1D,NIBLO
                IH = IJ2NEWIJ(ID1D)
                IF (BLK2GLO%IXLG(IH) == IPLUS2 .AND.  BLK2GLO%KXLT(IH) == K+1) EXIT
              ENDDO
              KLAT(IP,2,2) = IJ2NEWIJ(MIN(ID1D,NIBLO))
            ENDIF
          ELSE
            KLAT(IP,2,2) = KLAT(IP,2,1)
          ENDIF
        ENDIF

      ENDDO


!*    2.2 LONGITUDE NEIGHBOURS (KLON)
!         ---------------------------

      DO IJP=1,2
        DO IJ=IJS, IJL
           KLON(IJ,IJP) = 0
        ENDDO
      ENDDO

      DO IP = IJS, IJL
        I = BLK2GLO%IXLG(IP)
        K = BLK2GLO%KXLT(IP)
        IP1D = NEWIJ2IJ(IP)
        IF (I > 1) THEN
          IF (LLOCEANMASK(I-1,K)) KLON(IP,1) = IJ2NEWIJ(IP1D-1)
        ELSE
          IF (IPER == 1 .AND. LLOCEANMASK(NLONRGG(K),K)) THEN
            KLON(IP,1) = IP1D
            DO IH=2,NLONRGG(K)
              IF (LLOCEANMASK(IH,K)) KLON(IP,1) = KLON(IP,1)+1
            ENDDO
            KLON(IP,1) = IJ2NEWIJ(KLON(IP,1)) 
          ENDIF
        ENDIF
        IF (I < NLONRGG(K)) THEN
          IF (LLOCEANMASK(I+1,K)) KLON(IP,2) = IJ2NEWIJ(IP1D+1)
        ELSE
          IF (IPER == 1 .AND. LLOCEANMASK(1,K)) THEN
            KLON(IP,2) = IP1D
            DO IH=NLONRGG(K)-1,1,-1
              IF (LLOCEANMASK(IH,K)) KLON(IP,2) = KLON(IP,2)-1
            ENDDO
            KLON(IP,2) = IJ2NEWIJ(KLON(IP,2)) 
          ENDIF
        ENDIF
      ENDDO


!     COMPUTE THE CORNER GRID POINT AND THE CLOSEST GRID POINT ON EITHER
!     SIDE (ON A GIVEN LATITUDE) (KCOR)
!     ------------------------------------------------------------------

      IF (IPROPAGS == 2) THEN

        DO IC=1,2
          DO ICR=1,4
            DO IJ=IJS, IJL
               KCOR(IJ,ICR,IC) = 0
            ENDDO
          ENDDO
        ENDDO

        DO IP = IJS, IJL
          I = BLK2GLO%IXLG(IP)
          K = BLK2GLO%KXLT(IP)
          IP1D = NEWIJ2IJ(IP)
          XLON = REAL(I-1)*ZDELLO(K)

          IF (K > 1) THEN
!           CLOSEST GRID POINT IN SW GRID CORNER POINT
            XL=XLON-ZDELLO(K)

            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1

            IF (IPER == 1 .AND. IMIN < 1)  THEN
              IMIN = IMIN + NLONRGG(K-1)
              XMIN = XMIN + REAL(NLONRGG(K-1),JWRB)
            ENDIF
            IF (IMIN >= 1)  THEN
              IF (LLOCEANMASK(IMIN,K-1)) THEN
                DO ID1D = IP1D,1,-1
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IMIN .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                ENDDO
               KCOR(IP,3,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN SW GRID CORNER POINT
            IF (IMIN >= 1)  THEN
              IF (XMIN <= REAL(IMIN-1,JWRB)) THEN
                IMIN2=IMIN-1
                IF (IMIN <= 1) THEN
!!test                  IMIN2=1
                  IMIN2=NLONRGG(K-1)
                ELSE
                  IMIN2=IMIN-1
                ENDIF
              ELSE
                 IMIN2=IMIN+1
                IF (IMIN >= NLONRGG(K-1)) THEN
!!test                  IMIN2=NLONRGG(K-1)
                  IMIN2=1
                ELSE
                  IMIN2=IMIN+1
                ENDIF
              ENDIF

              IF (LLOCEANMASK(IMIN2,K-1)) THEN
                DO ID1D = IP1D,1,-1
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IMIN2 .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                ENDDO
                KCOR(IP,3,2) = IJ2NEWIJ(ID1D)
              ENDIF

            ENDIF


!           CLOSEST GRID POINT IN SE GRID CORNER
            XL=XLON+ZDELLO(K)
            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            IF (IPER == 1 .AND. IMIN > NLONRGG(K-1))  THEN
              IMIN = IMIN - NLONRGG(K-1)
              XMIN = XMIN - REAL(NLONRGG(K-1),JWRB)
            ENDIF
            IF (IMIN <= NLONRGG(K-1))  THEN
              IF (LLOCEANMASK(IMIN,K-1)) THEN
                DO ID1D = IP1D,1,-1
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IMIN .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                ENDDO
               KCOR(IP,2,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN SE GRID CORNER
            IF (IMIN <= NLONRGG(K-1))  THEN
              IF (XMIN <= REAL(IMIN-1,JWRB)) THEN
                 IMIN2=IMIN-1
                IF (IMIN <= 1) THEN
!!test                  IMIN2=1
                IMIN2=NLONRGG(K-1)
                ELSE
                  IMIN2=IMIN-1
                ENDIF
              ELSE
                 IMIN2=IMIN+1
                IF (IMIN >= NLONRGG(K-1)) THEN
!!test                  IMIN2=NLONRGG(K-1)
                  IMIN2=1
                ELSE
                  IMIN2=IMIN+1
                ENDIF
              ENDIF

              IF (LLOCEANMASK(IMIN2,K-1)) THEN
                DO ID1D = IP1D,1,-1
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IMIN2 .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                ENDDO
                KCOR(IP,2,2) = IJ2NEWIJ(ID1D) 
              ENDIF
            ENDIF


          ENDIF ! END K > 1


          IF (K < NGY) THEN
!           CLOSEST GRID POINT IN NW GRID CORNER 
            XL=XLON-ZDELLO(K)
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            IF (IPER == 1 .AND. IPLUS < 1)  THEN
              IPLUS = NLONRGG(K+1) + IPLUS 
              XPLUS = XPLUS + REAL(NLONRGG(K+1),JWRB)
            ENDIF

            IF (IPLUS >= 1)  THEN
              IF (LLOCEANMASK(IPLUS,K+1)) THEN
                DO ID1D = IP1D,NIBLO
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IPLUS .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                ENDDO
                KCOR(IP,4,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN NW GRID CORNER 
            IF (IPLUS >= 1)  THEN
              IF (XPLUS <= REAL(IPLUS-1,JWRB)) THEN
                 IPLUS2=IPLUS-1
                IF (IPLUS <= 1) THEN
!!test                  IPLUS2=1
                  IPLUS2=NLONRGG(K+1)
                ELSE
                  IPLUS2=IPLUS-1
                ENDIF
              ELSE
               IPLUS2=IPLUS+1
                IF (IPLUS >= NLONRGG(K+1)) THEN
!!test                  IPLUS2=NLONRGG(K+1)
                  IPLUS2=1
                ELSE
                  IPLUS2=IPLUS+1
                ENDIF
              ENDIF
 
              IF (LLOCEANMASK(IPLUS2,K+1)) THEN
                DO ID1D = IP1D,NIBLO
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IPLUS2 .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                ENDDO
                KCOR(IP,4,2) = IJ2NEWIJ(MIN(ID1D,NIBLO))
              ENDIF
            ENDIF


!           CLOSEST GRID POINT IN NE GRID CORNER
            XL=XLON+ZDELLO(K)
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            IF (IPER == 1 .AND. IPLUS > NLONRGG(K+1))  THEN
              IPLUS = IPLUS - NLONRGG(K+1)
              XPLUS = XPLUS - REAL(NLONRGG(K+1),JWRB)
            ENDIF

            IF (IPLUS <= NLONRGG(K+1))  THEN
              IF (LLOCEANMASK(IPLUS,K+1)) THEN
                DO ID1D = IP1D,NIBLO
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IPLUS .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                ENDDO
                KCOR(IP,1,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN NE GRID CORNER 
            IF (IPLUS <= NLONRGG(K+1))  THEN
              IF (XPLUS <= REAL(IPLUS-1,JWRB)) THEN
                 IPLUS2=IPLUS-1
                IF (IPLUS <= 1) THEN
!!test                  IPLUS2=1
                  IPLUS2=NLONRGG(K+1)
                ELSE
                  IPLUS2=IPLUS-1
                ENDIF
              ELSE
                 IPLUS2=IPLUS+1
                IF (IPLUS >= NLONRGG(K+1)) THEN
!!test                  IPLUS2=NLONRGG(K+1)
                  IPLUS2=1
                ELSE
                  IPLUS2=IPLUS+1
                ENDIF
              ENDIF

              IF (LLOCEANMASK(IPLUS2,K+1)) THEN
                DO ID1D = IP1D,NIBLO
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IPLUS2 .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                ENDDO
                KCOR(IP,1,2) =  IJ2NEWIJ(MIN(ID1D,NIBLO))
              ENDIF
            ENDIF

          ENDIF  ! END K < NGY

        ENDDO
      ENDIF  ! IPROPAGS == 2



!     COMPUTE THE CLOSEST GRID POINTS IN THE SW-NE and SE-NW DIRECTIONS
!     (I.E. GOING AT 45 DEGREE) (KRLAT, KRLON)
!     -----------------------------------------------------------------
      IF (IPROPAGS == 1 ) THEN

        DO ICL=1,2
          DO IJP=1,2
            DO IJ=IJS, IJL
               KRLAT(IJ,IJP,ICL) = 0
               KRLON(IJ,IJP,ICL) = 0
            ENDDO
          ENDDO
        ENDDO

        DO IP = IJS, IJL
          I = BLK2GLO%IXLG(IP)
          K = BLK2GLO%KXLT(IP)
          IP1D = NEWIJ2IJ(IP)
          XLON = REAL(I-1)*ZDELLO(K)

          IF (K > 1) THEN
!           CLOSEST GRID POINT IN SW CORNER
            XL=XLON-XDELLA
            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            IF (IPER == 1 .AND. IMIN .LT. 1)  THEN
              IMIN = IMIN + NLONRGG(K-1)
              XMIN = XMIN + REAL(NLONRGG(K-1),JWRB)
            ENDIF
            IF (IMIN >= 1)  THEN
              IF (LLOCEANMASK(IMIN,K-1)) THEN
                DO ID1D = IP1D,1,-1
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IMIN .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                ENDDO
               KRLON(IP,1,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN SW CORNER
            IF (IRGG == 1) THEN
              IF (IMIN >= 1)  THEN
                IF (XMIN <= REAL(IMIN-1,JWRB)) THEN
                  IF (IMIN <= 1) THEN
                    IMIN2=1
                  ELSE
                    IMIN2=IMIN-1
                  ENDIF
                ELSE
                  IF (IMIN >= NLONRGG(K-1)) THEN
                    IMIN2=NLONRGG(K-1)
                  ELSE
                    IMIN2=IMIN+1
                  ENDIF
                ENDIF

                IF (LLOCEANMASK(IMIN2,K-1)) THEN
                  DO ID1D = IP1D,1,-1
                    IH = IJ2NEWIJ(ID1D)
                    IF (BLK2GLO%IXLG(IH) == IMIN2 .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                  ENDDO
                  KRLON(IP,1,2) = IJ2NEWIJ(ID1D)
                ENDIF
              ENDIF
            ELSE
              KRLON(IP,1,2) = KRLON(IP,1,1)
            ENDIF

!           CLOSEST GRID POINT IN SE CORNER
            XL=XLON+XDELLA
            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            IF (IPER == 1 .AND. IMIN .GT. NLONRGG(K-1))  THEN
              IMIN = IMIN - NLONRGG(K-1)
              XMIN = XMIN - REAL(NLONRGG(K-1),JWRB)
            ENDIF
            IF (IMIN <= NLONRGG(K-1))  THEN
              IF (LLOCEANMASK(IMIN,K-1)) THEN
                DO ID1D = IP1D,1,-1
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IMIN .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                ENDDO
               KRLAT(IP,1,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN SE CORNER
            IF (IRGG == 1) THEN
              IF (IMIN <= NLONRGG(K-1))  THEN
                IF (XMIN <= REAL(IMIN-1,JWRB)) THEN
                  IF (IMIN <= 1) THEN
                    IMIN2=1
                  ELSE
                    IMIN2=IMIN-1
                  ENDIF
                ELSE
                  IF (IMIN >= NLONRGG(K-1)) THEN
                    IMIN2=NLONRGG(K-1)
                  ELSE
                    IMIN2=IMIN+1
                  ENDIF
                ENDIF

                IF (LLOCEANMASK(IMIN2,K-1)) THEN
                  DO ID1D = IP1D,1,-1
                    IH = IJ2NEWIJ(ID1D)
                    IF (BLK2GLO%IXLG(IH) == IMIN2 .AND. BLK2GLO%KXLT(IH) == K-1) EXIT
                  ENDDO
                  KRLAT(IP,1,2) = IJ2NEWIJ(ID1D)
                ENDIF
              ENDIF
            ELSE
              KRLAT(IP,1,2) = KRLAT(IP,1,1)
            ENDIF
          ENDIF


          IF (K < NGY) THEN
!           CLOSEST GRID POINT IN NW CORNER 
            XL=XLON-XDELLA
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            IF (IPER == 1 .AND. IPLUS < 1)  THEN
              IPLUS = NLONRGG(K+1) + IPLUS 
              XPLUS = XPLUS + REAL(NLONRGG(K+1),JWRB)
            ENDIF

            IF (IPLUS >= 1)  THEN
              IF (LLOCEANMASK(IPLUS,K+1)) THEN
                DO ID1D = IP1D,NIBLO
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IPLUS .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                ENDDO
                KRLAT(IP,2,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN NW CORNER 
            IF (IRGG == 1) THEN
              IF (IPLUS >= 1)  THEN
                IF (XPLUS <= REAL(IPLUS-1,JWRB)) THEN
                  IF (IPLUS <= 1) THEN
                    IPLUS2=1
                  ELSE
                    IPLUS2=IPLUS-1
                  ENDIF
                ELSE
                  IF (IPLUS >= NLONRGG(K+1)) THEN
                    IPLUS2=NLONRGG(K+1)
                  ELSE
                    IPLUS2=IPLUS+1
                  ENDIF
                ENDIF
 
                IF (LLOCEANMASK(IPLUS2,K+1)) THEN
                  DO ID1D = IP1D,NIBLO
                    IH = IJ2NEWIJ(ID1D)
                    IF (BLK2GLO%IXLG(IH) == IPLUS2 .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                  ENDDO
                  KRLAT(IP,2,2) = IJ2NEWIJ(MIN(ID1D,NIBLO))
                ENDIF
              ENDIF
            ELSE
              KRLAT(IP,2,2) = KRLAT(IP,2,1)
            ENDIF

!           CLOSEST GRID POINT IN NE CORNER
            XL=XLON+XDELLA
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            IF (IPER == 1 .AND. IPLUS > NLONRGG(K+1))  THEN
              IPLUS = IPLUS - NLONRGG(K+1)
              XPLUS = XPLUS - REAL(NLONRGG(K+1),JWRB)
            ENDIF

            IF (IPLUS <= NLONRGG(K+1))  THEN
              IF (LLOCEANMASK(IPLUS,K+1)) THEN
                DO ID1D = IP1D,NIBLO
                  IH = IJ2NEWIJ(ID1D)
                  IF (BLK2GLO%IXLG(IH) == IPLUS .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                ENDDO
                KRLON(IP,2,1) = IJ2NEWIJ(ID1D)
              ENDIF
            ENDIF

!           SECOND CLOSEST GRID POINT IN NE CORNER 
            IF (IRGG == 1) THEN
              IF (IPLUS <= NLONRGG(K+1))  THEN
                IF (XPLUS <= REAL(IPLUS-1,JWRB)) THEN
                  IF (IPLUS <= 1) THEN
                  IPLUS2=1
                  ELSE
                    IPLUS2=IPLUS-1
                  ENDIF
                ELSE
                  IF (IPLUS >= NLONRGG(K+1)) THEN
                    IPLUS2=NLONRGG(K+1)
                  ELSE
                    IPLUS2=IPLUS+1
                  ENDIF
                ENDIF

                IF (LLOCEANMASK(IPLUS2,K+1)) THEN
                  DO ID1D = IP1D,NIBLO
                    IH = IJ2NEWIJ(ID1D)
                    IF (BLK2GLO%IXLG(IH) == IPLUS2 .AND. BLK2GLO%KXLT(IH) == K+1) EXIT
                  ENDDO
                  KRLON(IP,2,2) = IJ2NEWIJ(MIN(ID1D,NIBLO))
                ENDIF
              ENDIF
            ELSE
              KRLON(IP,2,2) = KRLON(IP,2,1)
            ENDIF

          ENDIF

        ENDDO
      ENDIF  ! IF IPROPAGS == 1


!     COMPUTE THE WEIGHT FOR THE ADVECTION SCHEME ON AN IRREGULAR GRID


      DO IJP=1,2
        DO IJ=IJS, IJL
           WLAT(IJ,IJP)=1.0_JWRB
        ENDDO
      ENDDO

      IF (IPROPAGS == 2) THEN

        DO ICR=1,4
          DO IJ=IJS, IJL
             WCOR(IJ,ICR)=1.0_JWRB
          ENDDO
        ENDDO

      ELSEIF (IPROPAGS == 1 ) THEN

        DO IJP=1,2
          DO IJ=IJS, IJL
             WRLAT(IJ,IJP)=1.0_JWRB
             WRLON(IJ,IJP)=1.0_JWRB
          ENDDO
        ENDDO

      ENDIF

      IF (IRGG == 1) THEN
        DO IP = IJS, IJL
          I = BLK2GLO%IXLG(IP)
          K = BLK2GLO%KXLT(IP)
          IP1D = NEWIJ2IJ(IP)
          D0 = FLOAT(I-1)*ZDELLO(K)
          D3=D0-0.5_JWRB*ZDELLO(K)
          D5=D0+0.5_JWRB*ZDELLO(K)

          IF (K > 1) THEN
!           SOUTHERN POINT 
            XMIN = D0/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=REAL(IMIN-1,JWRB)*ZDELLO(K-1)
            D4=XP-0.5_JWRB*ZDELLO(K-1)
            D6=XP+0.5_JWRB*ZDELLO(K-1)
            IF (D0 <= XP) THEN
              IF (D4 <= D3 .OR. D6 <= D5) THEN
                 WLAT(IP,1) = 1.0_JWRB
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ELSE
              IF (D4 >= D3 .OR. D6 >= D5) THEN
                 WLAT(IP,1) = 1.0_JWRB
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ENDIF

            IF (IPROPAGS == 2) THEN
!             FOR CORNER TRANSPORT SCHEME
!             SOUTH-WEST CORNER 
              XL=D0-ZDELLO(K)
              XLL=XL-0.5_JWRB*ZDELLO(K)
              XLR=XL+0.5_JWRB*ZDELLO(K)

              XMIN = XL/ZDELLO(K-1)
              IMIN = NINT(XMIN) + 1
              XP=REAL(IMIN-1,JWRB)*ZDELLO(K-1)
              XPL=XP-0.5_JWRB*ZDELLO(K-1)
              XPR=XP+0.5_JWRB*ZDELLO(K-1)

              IF (XPL > XLL .AND. XPR < XLR) THEN
                D1=ZDELLO(K)
              ELSE
                D1=MIN(XLR,XPR)-MAX(XLL,XPL)
              ENDIF
              WCOR(IP,3)=MIN(1.0_JWRB,D1/ZDELLO(K))

!             FOR CORNER TRANSPORT SCHEME
!             SOUTH-EAST CORNER 
              XL=D0+ZDELLO(K)
              XLL=XL-0.5_JWRB*ZDELLO(K)
              XLR=XL+0.5_JWRB*ZDELLO(K)

              XMIN = XL/ZDELLO(K-1)
              IMIN = NINT(XMIN) + 1
              XP=REAL(IMIN-1,JWRB)*ZDELLO(K-1)
              XPL=XP-0.5_JWRB*ZDELLO(K-1)
              XPR=XP+0.5_JWRB*ZDELLO(K-1)

              IF (XPL > XLL .AND. XPR < XLR) THEN
                D1=ZDELLO(K)
              ELSE
                D1=MIN(XLR,XPR)-MAX(XLL,XPL)
              ENDIF
              WCOR(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K))

            ELSEIF (IPROPAGS == 1 ) THEN
!             SOUTH-WEST CORNER 
!             FOR DUAL ROTATED SCHEME
              XL=D0-XDELLA
              XMIN = XL/ZDELLO(K-1)
              IMIN = NINT(XMIN) + 1
              XP=REAL(IMIN-1,JWRB)*ZDELLO(K-1)
              D4=XP-0.5_JWRB*ZDELLO(K-1)
              D6=XP+0.5_JWRB*ZDELLO(K-1)
              IF (XL <= XP) THEN
                IF (D4 <= D0-TWOXDELLA) THEN
                   WRLON(IP,1)=1.0_JWRB
                ELSEIF (D6 <= D0) THEN
                   D2=XP-XL
                   D1=ZDELLO(K-1)-D2
                   WRLON(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
                ELSE
                   D2=D4-(D0-TWOXDELLA)
                   D1=TWOXDELLA-D2
                   WRLON(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ELSE
                IF (D6 >= D0) THEN
                   WRLON(IP,1)=1.0_JWRB
                ELSE IF (D4 >= D0-TWOXDELLA) THEN
                   D2=XL-XP
                   D1=ZDELLO(K-1)-D2
                   WRLON(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
                ELSE
                   D2=D0-D6
                   D1=TWOXDELLA-D2
                   WRLON(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ENDIF

!             SOUTH-EAST CORNER 
!             FOR DUAL ROTATED SCHEME
              XL=D0+XDELLA
              XMIN = XL/ZDELLO(K-1)
              IMIN = NINT(XMIN) + 1
              XP=REAL(IMIN-1,JWRB)*ZDELLO(K-1)
              D4=XP-0.5_JWRB*ZDELLO(K-1)
              D6=XP+0.5_JWRB*ZDELLO(K-1)
              IF (XL <= XP) THEN
                IF (D4 <= D0) THEN
                   WRLAT(IP,1)=1.0_JWRB
                ELSE IF (D6 <= D0+TWOXDELLA) THEN
                   D2=XP-XL
                   D1=ZDELLO(K-1)-D2
                   WRLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
                ELSE
                   D2=D4-D0
                   D1=TWOXDELLA-D2
                   WRLAT(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ELSE
                IF (D6 >= D0+TWOXDELLA) THEN
                   WRLAT(IP,1)=1.0_JWRB
                ELSE IF (D4 >= D0) THEN
                   D2=XL-XP
                   D1=ZDELLO(K-1)-D2
                   WRLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
                ELSE
                   D2=D0+TWOXDELLA-D6
                   D1=TWOXDELLA-D2
                   WRLAT(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ENDIF

            ENDIF  ! IPROPAGS 1 OR 2

          ENDIF  ! K > 1


          IF (K < NGY) THEN
!           NORTHERN POINT 
            XPLUS = D0/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=REAL(IPLUS-1,JWRB)*ZDELLO(K+1)
            D4=XP-0.5_JWRB*ZDELLO(K+1)
            D6=XP+0.5_JWRB*ZDELLO(K+1)
            IF (D0 <= XP) THEN
              IF (D4 <= D3 .OR. D6 <= D5) THEN
                 WLAT(IP,2) = 1.0_JWRB
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ELSE
              IF (D4 >= D3 .OR. D6 >= D5) THEN
                 WLAT(IP,2) = 1.0_JWRB
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ENDIF

            IF (IPROPAGS == 2) THEN

!             FOR CORNER TRANSPORT SCHEME
!             NORTH-WEST CORNER 
              XL=D0-ZDELLO(K)
              XLL=XL-0.5_JWRB*ZDELLO(K)
              XLR=XL+0.5_JWRB*ZDELLO(K)

              XPLUS = XL/ZDELLO(K+1)
              IPLUS = NINT(XPLUS) + 1
              XP=REAL(IPLUS-1,JWRB)*ZDELLO(K+1)
              XPL=XP-0.5_JWRB*ZDELLO(K+1)
              XPR=XP+0.5_JWRB*ZDELLO(K+1)

              IF (XPL > XLL .AND. XPR < XLR) THEN
                D1=ZDELLO(K)
              ELSE
                D1=MIN(XLR,XPR)-MAX(XLL,XPL)
              ENDIF
              WCOR(IP,4)=MIN(1.0_JWRB,D1/ZDELLO(K))

!             FOR CORNER TRANSPORT SCHEME
!             NORTH-EAST CORNER 
              XL=D0+ZDELLO(K)
              XLL=XL-0.5_JWRB*ZDELLO(K)
              XLR=XL+0.5_JWRB*ZDELLO(K)

              XPLUS = XL/ZDELLO(K+1)
              IPLUS = NINT(XPLUS) + 1
              XP=REAL(IPLUS-1,JWRB)*ZDELLO(K+1)
              XPL=XP-0.5_JWRB*ZDELLO(K+1)
              XPR=XP+0.5_JWRB*ZDELLO(K+1)

              IF (XPL > XLL .AND. XPR < XLR) THEN
                D1=ZDELLO(K)
              ELSE
                D1=MIN(XLR,XPR)-MAX(XLL,XPL)
              ENDIF
              WCOR(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K))

            ELSEIF (IPROPAGS == 1 ) THEN

!             NORTH-WEST CORNER 
!             FOR DUAL ROTATED SCHEME
              XL=D0-XDELLA
              XPLUS = XL/ZDELLO(K+1)
              IPLUS = NINT(XPLUS) + 1
              XP=REAL(IPLUS-1,JWRB)*ZDELLO(K+1)
              D4=XP-0.5_JWRB*ZDELLO(K+1)
              D6=XP+0.5_JWRB*ZDELLO(K+1)
              IF (XL <= XP) THEN
                IF (D4 <= D0-TWOXDELLA) THEN
                   WRLAT(IP,2)=1.0_JWRB
                ELSE IF (D6 <= D0) THEN
                   D2=XP-XL
                   D1=ZDELLO(K+1)-D2
                   WRLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
                ELSE
                   D2=D4-(D0-TWOXDELLA)
                   D1=TWOXDELLA-D2
                   WRLAT(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ELSE
                IF (D6 >= D0) THEN
                   WRLAT(IP,2)=1.0_JWRB
                ELSE IF (D4 >= D0-TWOXDELLA) THEN
                   D2=XL-XP
                   D1=ZDELLO(K+1)-D2
                   WRLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
                ELSE
                   D2=D0-D6
                   D1=TWOXDELLA-D2
                   WRLAT(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ENDIF

!             NORTH-EAST CORNER 
!             FOR DUAL ROTATED SCHEME
              XL=D0+XDELLA
              XPLUS = XL/ZDELLO(K+1)
              IPLUS = NINT(XPLUS) + 1
              XP=REAL(IPLUS-1,JWRB)*ZDELLO(K+1)
              D4=XP-0.5_JWRB*ZDELLO(K+1)
              D6=XP+0.5_JWRB*ZDELLO(K+1)
              IF (XL <= XP) THEN
                IF (D4 <= D0) THEN
                   WRLON(IP,2)=1.0_JWRB
                ELSE IF (D6 <= D0+TWOXDELLA) THEN
                   D2=XP-XL
                   D1=ZDELLO(K+1)-D2
                   WRLON(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
                ELSE
                   D2=D4-D0
                   D1=TWOXDELLA-D2
                   WRLON(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ELSE
                IF (D6 >= D0+TWOXDELLA) THEN
                   WRLON(IP,2)=1.0_JWRB
                ELSE IF (D4 >= D0) THEN
                   D2=XL-XP
                   D1=ZDELLO(K+1)-D2
                   WRLON(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
                ELSE
                   D2=D0+TWOXDELLA-D6
                   D1=TWOXDELLA-D2
                   WRLON(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
                ENDIF
              ENDIF

            ENDIF  ! IPROPAGS 1 OR 2

          ENDIF  ! K < NGY

        ENDDO
      ENDIF


END SUBROUTINE PROPCONNECT
