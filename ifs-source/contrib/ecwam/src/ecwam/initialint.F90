! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INITIALINT (IU06, NCA, NRA,                           &
     &                       NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,  &
     &                       AMOWEP, AMOSOP, AMOEAP, AMONOP,           &
     &                       NCAD, NRAD,                               &
     &                       RMONOP, RMOSOP, RMOWEP, RMOEAP,           &
     &                       ILONRGG, IPERIODIC,                       &
     &                       NWX, NWY,                                 &
     &                       DJ1, DII1, DIIP1, JJ, II, IIP)
! ----------------------------------------------------------------------    

!***  *INITIALINT* - INITIALIZES THE INTERPOLATION PARAMETERS TO BE USED
!                    IN INTERPOLATING FIELDS PASSED FROM THE ATMOSPHIRIC 
!                    SIDE TO THE WAVE SIDE OF THE MODEL.

!     PURPOSE.                                                          
!     --------                                                          

!     IT COMPUTES THE WEIGHT VALUES AND THE 4 GRID POINTS TO BE USED
!     FOR INTERPOLATING THE FIELDS PASSED FROM THE ATMOSPHIRIC SIDE TO  
!     THE WAVE SIDE OF THE MODEL. "INITIALINT" MUST BE CALLED BEFORE
!     ANY CALL TO "FLDINTER".

!     LIMITATIONS.
!     ------------
!     1. IT IS ASSUMED THAT THE FIELDS ARE NOT ANGULAR (DIRECTIONAL).
!     2. IT IS ASSUMED THAT THERE IS NO MISSING DATA POINTS.
!     3. IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!        (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!        *IU06*   - OUTPUT UNIT.
!       ATMOSPHERIC MODEL GRID:
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!       WAVE MODEL GRID SPECIFICATION:
!        *NGX  *  - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY  *  - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *KRGG*   - GRID DEFINITION PARAMETER (0=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELLO* - GRID POINT SPACING PER LATITUDES. 
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).
!        *NCAD*     - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR DECODED
!        *NRAD*     - NUMBER OF ATM. ROWS OF LATITUDES DECODED
!        *RMOWEP* - MOST WESTERN LONGITUDE IN ATM GRID
!        *RMOSOP* - MOST SOUTHERN LATITUDE IN ATM GRID
!        *RMOEAP* - MOST EASTERN LONGITUDE IN ATM GRID
!        *RMONOP* - MOST NORTHERN LATITUDE IN ATM GRID
!        *ILONRGG*- NUMBER OF GRID POINTS FOR EACH LATITUDE (ATM. GRID)
!        *IPERIODIC* - SPECIFIES IF ATM. GRID IS PERIODIC (=1) OR NOT (=0)
!        *NWX*    - FIRST DIMENSION FOR THE INTERPOLATION COEFFICIENTS 
!        *NWY*    - SECOND DIMENSION FOR THE INTERPOLATION COEFFICIENTS

!       OUTPUT: INTERPOLATION COEFFICIENTS & INDECES
!        *DJ1*    - COEFFICIENT
!        *DII1M*  - COEFFICIENT
!        *DIIP1M* - COEFFICIENT
!        *JJM*    - INDEX
!        *IIM*    - INDEX
!        *IIPM*   - INDEX

!     EXTERNALS.                                                        
!     ----------                                                        

!        none
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

      USE YOWGRIB   ,ONLY : IGRIB_GET_VALUE

! ----------------------------------------------------------------------

      IMPLICIT NONE 

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, NCA, NRA
      INTEGER(KIND=JWIM), INTENT(IN) :: NGX, NGY, KRGG
      INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(IN) :: KLONRGG
      REAL(KIND=JWRB), INTENT(IN) :: XDELLA
      REAL(KIND=JWRB), DIMENSION(NGY), INTENT(IN) :: ZDELLO
      REAL(KIND=JWRB), INTENT(IN) :: AMOWEP, AMOSOP, AMOEAP, AMONOP
      INTEGER(KIND=JWIM), INTENT(IN) :: NCAD, NRAD
      REAL(KIND=JWRB), INTENT(IN) :: RMONOP, RMOSOP, RMOWEP, RMOEAP
      INTEGER(KIND=JWIM), DIMENSION(NRA), INTENT(IN) :: ILONRGG
      INTEGER(KIND=JWIM), INTENT(IN) :: IPERIODIC
      INTEGER(KIND=JWIM), INTENT(IN) :: NWX, NWY

      INTEGER(KIND=JWIM), DIMENSION(NWY), INTENT(OUT) :: JJ
      INTEGER(KIND=JWIM), DIMENSION(NWX, NWY), INTENT(OUT) :: II, IIP
      REAL(KIND=JWRB), DIMENSION(NWY), INTENT(OUT) :: DJ1
      REAL(KIND=JWRB), DIMENSION(NWX, NWY), INTENT(OUT) :: DII1, DIIP1


      INTEGER(KIND=JWIM) :: I, J, K, JSN
      INTEGER(KIND=JWIM) :: KSN, JJ1, KSN1, III, IIIP

      REAL(KIND=JWRB) :: DELLA, DELLO
      REAL(KIND=JWRB) :: XK, XLAT, XI, XII, XIIP
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NRAD) :: RDELLO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INITIALINT',0,ZHOOK_HANDLE)

      DELLA=(RMONOP-RMOSOP)/MAX(1,NRAD-1)

      IF (IPERIODIC == 1) THEN
        DELLO=360.0_JWRB/MAX(1,NCAD)
        DO J=1,NRAD
          JSN=NRAD-J+1
          RDELLO(JSN) = 360.0_JWRB/MAX(1,ILONRGG(JSN)) 
        ENDDO
      ELSE
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NCAD-1)
        DO J=1,NRAD
          JSN=NRAD-J+1
          RDELLO(JSN) = (RMOEAP-RMOWEP)/MAX(1,ILONRGG(JSN)-1) 
        ENDDO
      ENDIF


!     INTERPOLATION WEIGHT TO WAVE MODEL GRID
!     ----------------------------------------

      CALL GSTATS(1436,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(K,JSN,XLAT,XK,KSN,JJ1,KSN1,I,XI,XII,III,XIIP,IIIP )
      DO K=1,NGY                                                

        JSN=NGY-K+1

        XLAT = AMOSOP + (JSN-1)*XDELLA
        XLAT=MIN(MAX(XLAT,RMOSOP),RMONOP)
        XK = RMONOP - XLAT 
        XK = XK/DELLA+1.0_JWRB
        JJ(K) = MAX(1,INT(XK))
        KSN = NRAD-JJ(K)+1
        JJ1 = MIN(JJ(K)+1,NRAD)
        KSN1 = NRAD-JJ1+1
        DJ1(K)=XK-REAL(JJ(K))

        DO I=1,KLONRGG(JSN)

          XI = AMOWEP+(I-1)*ZDELLO(JSN) - RMOWEP
          XI = MOD(XI+720.0_JWRB,360.0_JWRB)

          XII = XI/RDELLO(KSN)+1.0_JWRB
          III = MAX(1-IPERIODIC,INT(XII))
          II(I,K) = MIN(III,ILONRGG(KSN))
          DII1(I,K) = XII-REAL(II(I,K))

          XIIP = XI/RDELLO(KSN1)+1.0_JWRB
          IIIP = MAX(1-IPERIODIC,INT(XIIP))
          IIP(I,K) = MIN(IIIP,ILONRGG(KSN1)) 
          DIIP1(I,K)=XIIP-REAL(IIP(I,K))

        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1436,1)

      IF (LHOOK) CALL DR_HOOK('INITIALINT',1,ZHOOK_HANDLE)

      END SUBROUTINE INITIALINT
