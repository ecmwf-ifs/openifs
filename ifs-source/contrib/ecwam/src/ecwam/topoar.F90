! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE TOPOAR (IU01, BATHY)

! ----------------------------------------------------------------------

!**** *TOPOAR* - ARRANGE SUBGRID TOPOGRAPHY.

!     S. HASSELMANN     MPIFM           1/6/86.

!     MODIFIED BY       H. GUNTHER      1/4/90  -  REARANGEMENT OF CODE.
!                       J. BIDLOT    FEB 2003 - 
!                                    INTRODUCE OBSTRUCTION COEFFICIENTS
!                                    MOVED TO MUBUF NOW !
!                       J. BIDLOT    NOV 2005 - 
!                                    INTRODUCE NEW INPUT FORMAT FOR
!                                    BATHYMETRY (the old format can still
!                                    be used).
!*    PURPOSE.
!     --------

!       TO READ IN TOPOGRAPHY ON INPUT GRID AND CONVERT TO OUTPUT GRID.

!**   INTERFACE.
!     ----------

!       *CALL* *TOPOAR (IU01,BATHY)*
!          *IU01*  -  LOGICAL INPUT UNIT OF TOPOGRAPHIC DATA.
!          *BATHY*   -  BATHYMETRY DATA IN OUTPUT GRID.

!     METHOD.
!     -------

!       IF OLD TOPAGRAPHY INPUT:
!       THE TOPOGRAPHY MUST BE ON A REGULAR LATITUDE-LONGITUDE
!       GRID ARRANGED FROM SOUTH TO NORTH, AND FROM WEST TO EAST.
!       IT IS ASSUMED THAT NEGATIVE VALUES ARE SEA DEPTHS (WHICH
!       ARE CONVERTED TO POSITIVE) AND THAT POSITIVE VALUES ARE
!       LAND ELEVATIONS (WHICH ARE CONVERTED TO -999 IDENTIFIERS).

!       THE TOPOGRAPHIC DATA IS READ IN ON THE INPUT GRID AND IT
!       IS STORED ONLY FOR THOSE LATITUDES WITHIN THE REQUESTED
!       GRID. THEN THE TOPOGRAPHIC DATA IS FURTHER RESTRICTED
!       TO LIE WITHIN THE SUBGRID LONGITUDES. IT IS THEN PUT ON
!       THE REQUESTED SUBGRID LAT-LONG RESOLUTION, ALWAYS USING
!       THE NEAREST POINT.
!       FINALLY THE SUBGRID TOPOGRAPHY MAY BE MANUALLY ADJUSTED BY
!       MEANS OF THE CARD INPUT AND A PRINTER OUTPUT IS DONE.

!       FOR NEW TOPOGRAPHY INPUT:
!       THE INPUT IS ALREADY ON THE DESIRED GRID AND THE REDUCTION
!       FACTORS DUE TO OBSTRUCTIONS ARE GIVEN AS WELL
!       !!! THE REDUCTION FACTORS ARE READ AND PROCESSED IN MUBUF!

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *ADJUST*    - CORRECTS LONGITUDE INPUT.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCINP  , ONLY : NOUT     ,XOUTW    ,XOUTS    ,XOUTE    ,    &
     &            XOUTN    ,NOUTD
      USE YOWMAP   , ONLY : NGX      ,NGY      ,AMOWEP   ,AMOSOP   ,    &
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,    &
     &            NLONRGG  ,LLOBSTRCT
      USE YOWTEST  , ONLY : IU06, ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "adjust.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU01
      REAL(KIND=JWRB), DIMENSION(NGX,NGY), INTENT(OUT) :: BATHY


      INTEGER(KIND=JWIM) :: I, J, K, JH, L, IX, IAA, IS
      INTEGER(KIND=JWIM) :: KLONRGG  
      INTEGER(KIND=JWIM) :: MLON, NLATMAX, KMAX, NLAT, N1, N2, LAST
      INTEGER(KIND=JWIM) :: NMINADJT
      INTEGER(KIND=JWIM) :: ILW, NLON, IH, NJ, JJ, NL, JRGG, KXLO, KAMOEAP
      INTEGER(KIND=JWIM) :: ILEN, IPAGE, IA, IE
      INTEGER(KIND=JWIM), ALLOCATABLE :: IDUM(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IA2H(:), IA1(:,:)

      REAL(KIND=JWRB) :: XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      REAL(KIND=JWRB) :: XLAT, XLON, XLAG, XLW, XLA, XLO, XLOH, XLOG
      REAL(KIND=JWRB), ALLOCATABLE :: XA2H(:), XA1(:,:)

      CHARACTER(LEN=1), ALLOCATABLE :: AX(:), AXX(:)
      CHARACTER(LEN=5) :: CX
      CHARACTER(LEN=11) :: FORMT
      CHARACTER(LEN=14) :: CHEADER 

      LOGICAL :: LLREALIN

! ----------------------------------------------------------------------

!*    1. READING THE TOPOGRAPHY OF THE INPUT GRID AND STORING THOSE
!*       LATITUDES WITHIN THE OUTPUT SUBGRID AREA.
!        -----------------------------------------

!       MLON - NUMBER OF TOPOGRAPHIC POINTS PER GRID LATITUDE.
!       KMAX - NUMBER OF RECORDS OF INPUT GRID PER LATITUDE.
!       XLAT - LATITUDE OF CURRENT GRID DATA.
!       NLAT - NUMBER OF GRID LATITUDES STORED.
!       XLAG - LATITUDE OF FIRST GRID LATITUDE STORED.

      REWIND (UNIT=IU01)
!     DETERMINE WHICH TYPE OF FILE IS USED
      READ (IU01,'(a14)') CHEADER 
      IF (CHEADER == 'WAM BATHYMETRY') THEN
!       REAL BATHYMETRY INPUT
        LLREALIN=.TRUE.
        WRITE (IU06,*) ' BATHYMETRY FROM REAL INPUT DATA' 
      ELSE
!       INTEGER BATHYMETRY INPUT
        LLREALIN=.FALSE.
        WRITE (IU06,*) ' BATHYMETRY FROM INTEGER INPUT DATA' 
      ENDIF
      REWIND (UNIT=IU01)

      IF (LLREALIN) THEN
        READ (IU01,'(a14)') CHEADER 
        READ (IU01,'(6F13.8)') XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      ELSE
        READ (IU01,'(6F10.5)') XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      ENDIF
      CALL ADJUST (XLOW, XLOE)

      WRITE (IU06,'(1H1,'' INPUT GRID''/)')
      WRITE (IU06,'(3X,''RESOLUTION LAT-LON '',2F8.3)') XDELA, XDELO
      WRITE (IU06,'(3X,'' SOUTHERN LAT '','' NORTHERN LAT '',           &
     &                 '' WESTERN LONG '','' EASTERN LONG'',            &
     &                 /,2X,4F13.8)') XLAS, XLAN, XLOW, XLOE

      BATHY(:,:) = 999.0_JWRB

      IF (LLOBSTRCT) THEN
!     NEW TREATMENT OF BATHYMETRY INPUT
!     ---------------------------------
        WRITE (IU06,*) '  NEW TREATMENT OF BATHYMETRY INPUT '
        WRITE (IU06,*) ' '
!       TEST DOMAIN
        IF (NINT(10000*XDELA) /= NINT(10000*XDELLA) .OR.                 &
     &      NINT(10000*XDELO) /= NINT(10000*XDELLO) .OR.                 &
     &      NINT(10000*XLAS) /= NINT(10000*AMOSOP) .OR.                  &
     &      NINT(10000*XLAN) /= NINT(10000*AMONOP) .OR.                  &
     &      NINT(10000*XLOW) /= NINT(10000*AMOWEP) .OR.                  &
     &      NINT(1000*XLOE) /= NINT(1000*AMOEAP) ) THEN
          WRITE (IU06,*) ' *******************************************'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR        *'
          WRITE (IU06,*) ' *      ===========================        *'
          WRITE (IU06,*) ' * THE INPUT AND COMPUTED DOMAIN DEFINITION*'
          WRITE (IU06,*) ' * DO NOT AGREE :'
          WRITE (IU06,*) ' * XDELLA ',XDELA,XDELLA
          WRITE (IU06,*) ' * XDELLO ',XDELO,XDELLO
          WRITE (IU06,*) ' * AMOSOP ',XLAS,AMOSOP
          WRITE (IU06,*) ' * AMONOP ',XLAN,AMONOP
          WRITE (IU06,*) ' * AMOWEP ',XLOW,AMOWEP
          WRITE (IU06,*) ' * AMOEAP ',XLOE,AMOEAP
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                 *'
          WRITE (IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF

        DO K=1,NGY

!       READ THE NUMBER OF POINTS PER LATITUDE
          READ (IU01,*) KLONRGG  

          IF (KLONRGG /= NLONRGG(K)) THEN
            WRITE (IU06,*) ' ******************************************'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR       *'
            WRITE (IU06,*) ' *      ===========================       *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' * THE INPUT NUMBER OF POINTS PER LATITUDE '
            WRITE (IU06,*) ' * IS NOT THE SAME AS THE COMPUTED ONE !!  '
            WRITE (IU06,*) ' * KLONRGG = ',KLONRGG,K 
            WRITE (IU06,*) ' * NLONRGG = ',NLONRGG(K),K 
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                *'
            WRITE (IU06,*) ' ******************************************'
            CALL ABORT1
          ENDIF
        ENDDO

!       READ THE BATHYMETRY DATA

        ALLOCATE(IDUM(NGX))
        CX='     '
        FORMT='          ' 
        DO K=1,NGY
          IF (LLREALIN) THEN
            WRITE(CX,'(I5.5)') NLONRGG(1)
            FORMT='('//CX//'F9.2)'
            DO IS = 1,NLONRGG(K),NLONRGG(1)
              READ (IU01,FORMT) (BATHY(IX,K),IX=IS,MIN(IS+NLONRGG(1)-1,NLONRGG(K)))
            ENDDO
          ELSE
            WRITE(CX,'(I4.4)') NLONRGG(K)
            FORMT='('//CX//'I4)'
            READ (IU01,FORMT) (IDUM(IX),IX=1,NLONRGG(K))
            DO IX=1,NLONRGG(K)
              BATHY(IX,K)=REAL(IDUM(IX),JWRB)
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(IDUM)
  
!       THE REDUCTION FACTORS ARE NOW READ IN MUBUF (the rest of IU01)

      ELSE
!     OLD TREATMENT OF BATHYMETRY INPUT
!     ---------------------------------
      MLON = NINT((XLOE-XLOW)/XDELO+1.0_JWRB)

      NLATMAX=NINT((XLAN-XLAS)/XDELA+1.0_JWRB)

      ALLOCATE(AX(MLON))
      ALLOCATE(IA2H(MLON), IA1(MLON,NLATMAX))
      ALLOCATE(XA2H(MLON), XA1(MLON,NLATMAX))

      KMAX=(MLON+11)/12

      XLAT=XLAS-XDELA
      NLAT=0
1005  CONTINUE
      NLAT=NLAT+1
      IF (NLAT > NLATMAX) THEN
        WRITE (IU06,*) ' ******************************************'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR       *'
        WRITE (IU06,*) ' *      ===========================       *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' * NUMBER OF LATITUDES IN INPUT GRID      *'
        WRITE (IU06,*) ' * EXCEEDS DIMENSION                      *'
        WRITE (IU06,*) ' * DIMENSION IS      NLATMAX = ', NLATMAX
        WRITE (IU06,*) ' * LAST LATITUDE READ IS  XLAT  = ', XLAT
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                *'
        WRITE (IU06,*) ' ******************************************'
        CALL ABORT1
      ENDIF
1010  CONTINUE
      XLAT=XLAT+XDELA
      DO K=1,KMAX
        N1 = 12*(K-1)+1
        N2 = MIN(12*K,MLON)
        IF (LLREALIN) THEN
          READ (IU01, '(12(F9.2))', END=1017) (XA1(IAA,NLAT),IAA=N1,N2)
        ELSE
          READ (IU01, '(12(I5,A1))', END=1017) (IA1(IAA,NLAT),AX(IAA),IAA=N1,N2)
        ENDIF
      ENDDO

      IF (XLAT+0.5_JWRB*XDELA <= AMOSOP) GO TO 1010
      IF (NLAT == 1) XLAG=XLAT
      IF (.NOT.LLREALIN) THEN
      DO I=1,MLON
          IF (AX(I) == 'E' .AND. IA1(I,NLAT) < 0) IA1(I,NLAT)=-IA1(I,NLAT)
          IF (AX(I) == 'D' .AND. IA1(I,NLAT) > 0) IA1(I,NLAT)=-IA1(I,NLAT)
          IF (AX(I) == 'E' .AND. IA1(I,NLAT) == 0) IA1(I,NLAT)=1
        ENDDO
      ENDIF
      IF (XLAT+0.5_JWRB*XDELA > AMONOP) GO TO 1020
      GO TO 1005
1017  CONTINUE
      NLAT=NLAT-1
      WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                           +'
      WRITE (IU06,*) ' +     WARNING ERROR IN SUB. TOPOAR          +'
      WRITE (IU06,*) ' +     ============================          +'
      WRITE (IU06,*) ' +                                           +'
      WRITE (IU06,*) ' + END OF FILE ON INPUT GRID  UNIT :', IU01
      WRITE (IU06,*) ' + NORTH GRID BOUNDARY CHANGED TO XLAT = ',XLAT
      WRITE (IU06,*) ' +                                           +'
      WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
1020  CONTINUE

! ----------------------------------------------------------------------

!*    2. STORING TOPOGRAPHIC DATA AT LONGITUDES WITHIN SUBGRID AREA.
!        ----------------------------------------------------------

!       ILW      -  INDEX OF NEAREST GRID LONGITUDE EQUIVALENT TO
!                     WESTERN SUBGRID BOUNDARY.
!       NLON     -  NUMBER OF GRID LONGITUDES WITHIN SUBGRID AREA.
!       IH,IH1   -  GRID LONGITUDE NUMBER.


      XLW= MOD(AMOWEP-XLOW+720.0_JWRB,360.0_JWRB)
!!!!      ILW= NINT(XLW/XDELO-0.0001_JWRB)
      ILW= NINT(XLW/XDELO)

      NLON=NINT((AMOEAP-AMOWEP)/XDELO+1.0_JWRB)
      DO J=1,NLAT
        IH=ILW
        DO I=1,NLON
          IH=IH+1
          IF (IH <= 0) IH=IH+MLON
          IF (IH > MLON) IH=IH-MLON
          IF (LLREALIN) THEN
            XA2H(I) =XA1(IH,J)
          ELSE
            IA2H(I) =IA1(IH,J)
          ENDIF
        ENDDO
        DO I=1,NLON
          IF (LLREALIN) THEN
            XA1(I,J)=XA2H(I)
          ELSE
            XA1(I,J)=REAL(IA2H(I),JWRB)
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. PUT ON XDELLO BY XDELLA SUBGRID.
!        --------------------------------

!     XLA  - LATITUDE OF OUTPUT GRID.
!     XLAG - LATITUDE OF INPUT GRID.
!     XLO  - LONGITUDE OF OUTPUT GRID.
!     XLOG - LONGITUDE OF INPUT GRID.

      XLA=AMOSOP
      XLAG=XLAG-XDELA
      NJ=0
      JJ=0

!  LOOP THROUGH LATITUDES

      DO J=1,NLAT
        XLAG=XLAG+XDELA
 3010   CONTINUE
        IF (XLA < XLAG-0.5_JWRB*XDELA .OR. XLA >= XLAG+0.5_JWRB*XDELA) GO TO 3070
        NJ=NJ+1
        IF (NJ > NGY) THEN
          WRITE (IU06,*) ' ******************************************'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR       *'
          WRITE (IU06,*) ' *      ===========================       *'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' * NUMBER OF LATITUDES IN OUTPUT GRID     *'
          WRITE (IU06,*) ' * EXCEEDS DIMENSION.                     *'
          WRITE (IU06,*) ' * DIMENSION IS            NGY = ', NGY
          WRITE (IU06,*) ' *                         NJ  = ', NJ
          WRITE (IU06,*) ' * LAST LATITUDE USED IS   XLA = ', XLA
          WRITE (IU06,*) ' *                         XLAG= ', XLAG
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                *'
          WRITE (IU06,*) ' ******************************************'
          CALL ABORT1
        ENDIF

!  LOOP THROUGH LONGITUDES

        NL=0
        XLOH = XLOW + (REAL(ILW,JWRB)-1.5_JWRB)*XDELO+720.0_JWRB
        XLO=AMOWEP
        IF (XLO < 0.0_JWRB) XLO =XLO +360.0_JWRB
        DO I=1,NLON
          XLO = MOD(XLO+720.0_JWRB,360.0_JWRB)
          XLOG = MOD(XLOH + REAL(I,JWRB)*XDELO,360.0_JWRB)
 3030     CONTINUE
          IF (XLO < XLOG) XLO  = XLO  + 360.0_JWRB

          IF (XLO <= XLOG+XDELO) THEN
            NL=NL+1
            IF (NL > NGX) THEN
              WRITE (IU06,*) ' **************************************'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR   *'
              WRITE (IU06,*) ' *      ===========================   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * NUMBER OF LONG. IN OUTPUT GRID     *'
              WRITE (IU06,*) ' * EXCEEDS DIMENSION.                 *'
              WRITE (IU06,*) ' * DIMENSION IS            NGX = ', NGX
              WRITE (IU06,*) ' * LAST LONGITUDE READ IS  XLO = ', XLO
              WRITE (IU06,*) ' * LONGITUDE INCREMENT IS  ', ZDELLO(J)
              WRITE (IU06,*) ' * INDEX J  IS  ', J
              WRITE (IU06,*) ' * INDEX JRGG  IS  ', JRGG
              WRITE (IU06,*) ' * INDEX I  IS  ', I
              WRITE (IU06,*) ' * NLAT IS  ', NLAT
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED            *'
              WRITE (IU06,*) ' **************************************'
              CALL ABORT1
            ENDIF
            BATHY(NL,NJ)=XA1(I,J)
            JRGG = 1 +NINT((XLA-AMOSOP)/XDELLA)
            XLO = AMOWEP + REAL(NL,JWRB)*ZDELLO(JRGG)
            KXLO=NINT(XLO*100.0_JWRB)
            KAMOEAP=NINT(AMOEAP*100.0_JWRB)
            IF (KXLO <= KAMOEAP) THEN
              XLO = MOD(XLO+720.0_JWRB,360.0_JWRB)
              GOTO 3030
            ENDIF
          ENDIF
        ENDDO

        JJ=JJ+1
        XLA=AMOSOP+REAL(JJ,JWRB)*XDELLA

        IF ((XLA-0.5_JWRB*XDELLA) > AMONOP) THEN
          GOTO 3080
        ELSE
          GOTO 3010
        ENDIF
 3070   CONTINUE
      ENDDO
 3080 CONTINUE
      IF (NJ /= NGY .OR. NL > NGX) THEN
        WRITE (IU06,*) ' *****************************************'
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR      *'
        WRITE (IU06,*) ' *      ===========================      *'
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' * NUMBER OF LONGITUDES OR LATITUDES IN  *'
        WRITE (IU06,*) ' * IS NOT EQUAL TO EXPECTED NUMBER       *'
        WRITE (IU06,*) ' * LATITUDES  FOUND      NJ  = ', NJ
        WRITE (IU06,*) ' * LATITUDES  EXPECTED   NGY = ', NGY
        WRITE (IU06,*) ' * LONGITUDES FOUND      NL  = ', NL
        WRITE (IU06,*) ' * LONGITUDES EXPECTED   NGX = ', NGX
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED               *'
        WRITE (IU06,*) ' *****************************************'
        CALL ABORT1
      ENDIF

      ENDIF  ! end on new vs old type bathymetry


! ----------------------------------------------------------------------

!*    4. CONVERT INPUT DEPTH TO MODEL DEPTH
!*       POSITIVE SEA DEPTH IN METRES (-999  FOR LAND).
!        ----------------------------------------------

      DO J=1,NGY
        DO I=1,NGX
          IF (BATHY(I,J) < 0.0_JWRB) THEN
            BATHY(I,J) = -BATHY(I,J)
          ELSE
            BATHY(I,J) = -999.0_JWRB
          ENDIF
        ENDDO
      ENDDO


! ----------------------------------------------------------------------

!*    5. MANUAL ADJUSTMENT OF TOPOGRAPHY.
!        --------------------------------

      IF (NOUT /= 0) THEN
        XLAT=AMOSOP-XDELLA
        DO J=1,NGY
          XLAT=XLAT+XDELLA
          XLON=AMOWEP-ZDELLO(J)
          IF (XLON < 0.0_JWRB) XLON=360.0_JWRB+XLON
          DO I=1,NGX
            XLON=XLON+ZDELLO(J)
            IF (XLON >= 360.0_JWRB) XLON=XLON-360.0_JWRB
            DO JH = 1,NOUT
              IF (XLON < XOUTW(JH)) XLON=XLON+360.0_JWRB
              IF (XLON > XOUTE(JH)) XLON=XLON-360.0_JWRB
              IF (XLON >= XOUTW(JH) .AND. XLAT >= XOUTS(JH) .AND.       &
     &         XLON <= XOUTE(JH) .AND. XLAT <= XOUTN(JH))               &
     &         BATHY(I,J) = REAL(NOUTD(JH),JWRB)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    7. AID TO USERS - SIMPLE PLOT OF GRID.
!        ------------------------------------

      WRITE (IU06,'(''0NUMBER OF LATITUDES IS       NGY = '',I5)') NGY
      WRITE (IU06,'('' MOST SOUTHERN LATITUDE IS AMOSOP = '',F7.3)') AMOSOP
      WRITE (IU06,'('' MOST NORTHERN LATITUDE IS AMONOP = '',F7.3)') AMONOP 
      WRITE (IU06,'('' LATITUDE INCREMENT IS     XDELLA = '',F7.3)') XDELLA
      WRITE (IU06,'(''0MAX NUMBER OF LONGITUDES IS  NGX = '',I5)') NGX
      WRITE (IU06,'('' MOST WESTERN LONGITUDE IS AMOWEP = '',F7.3)') AMOWEP
      WRITE (IU06,'('' MOST EASTERN LONGITUDE IS AMOEAP = '',F7.3)') AMOEAP
      WRITE (IU06,                                                      &
     &     '('' LONGITUDE INCREMENT AS FUNCTION OF LATITUDE IS'')')
      WRITE (IU06,'(10F8.3)') ZDELLO

      IF (ALLOCATED(AX)) DEALLOCATE(AX)
      IF (ALLOCATED(AXX)) DEALLOCATE(AXX)
      IF (ALLOCATED(IA2H)) DEALLOCATE(IA2H)
      IF (ALLOCATED(XA2H)) DEALLOCATE(XA2H)
      IF (ALLOCATED(IA1)) DEALLOCATE(IA1)
      IF (ALLOCATED(XA1)) DEALLOCATE(XA1)

      END SUBROUTINE TOPOAR
