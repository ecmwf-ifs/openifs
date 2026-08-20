! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE FLDINTER (IU06, NGPTOTG, NC, NR, NFIELDS, FIELDS,        &
     &                     NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,       &
     &                     IFROMIJ, JFROMIJ, KIJS, KIJL, KIJLMAX,         &
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,     &
     &                     ILONRGG, IJBLOCK, PMISS,                       &
     &                     LADEN, ROAIR, LGUST, WSTAR0, LLKC, LWCUR,      &
     &                     LLINTERPOL,                                    &
     &                     NWX, NWY, DJ1M, DII1M, DIIP1M, JJM, IIM, IIPM, &
     &                     MASK_IN,                                       &
     &                     NXS, NXE, NYS, NYE, FIELDG)
! ----------------------------------------------------------------------    

!***  *FLDINTER* - INTERPOLATION OF ATMOSPHERIC FIELDS OVER
!                  THE WAVE MODEL GRID USING BILINEAR
!                  INTERPOLATION.

!      S. ABDALLA   ECMWF   DECEMBER 2000. 
!      J. BIDLOT  ECMWF AUGUST 2006  ADDED SEA ICE FRACTION WHICH
!                                    CONTAINS MISSING DATA: NEAREST
!                                    SEA POINT IS TAKEN (IF AVAILABLE).

!     PURPOSE.                                                          
!     --------                                                          

!     IT INTERPOLATES ATMOSPHERIC DATA FIELDS TO PRODUCE A FIELD 
!     ON THE WAVE MODEL GRID USING BILINEAR INTERPOLATION.
!     IT MAKES USE OF THE INTERPOLATION COEFFICIENTS PRODUCED BY 
!     SUBROUTINE "INITIALINT" WHICH SHOULD BE CALLED BEFORE ANY
!     CALL TO "FLDINTER".
!     IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!     (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!        *IU06*   - OUTPUT UNIT.
!        ATMOSPHERIC MODEL GRID AND FIELD (INPUT):
!        *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *ILONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (INPUT FIELD)
!        *FIELDS* - ATMOSPHERIC FIELDS AS FOLLOWS:
!                   FIELDS(:,1) = U COMPONENT OF NEUTRAL WIND SPEED (U10)
!                   FIELDS(:,2) = V COMPONENT OF NEUTRAL WIND SPEED (V10)
!                   FIELDS(:,3) = AIR DENSITY
!                   FIELDS(:,4) = CONVECTIVE VELOCITY SCALE (w*)
!                   FIELDS(:,5) = SEA ICE FRACTION 
!                   FIELDS(:,6) = LAKE FRACTION 
!                   FIELDS(:,7) = U COMPONENT OF ATMOSPHERIC SURFACE STRESS
!                   FIELDS(:,8) = V COMPONENT OF ATMOSPHERIC SURFACE STRESS
!                   FIELDS(:,9) = U COMPONENT OF OCEAN SURFACE CURREMT 
!                   FIELDS(:,10)= V COMPONENT OF OCEAN SURFACE CURREMT 

!        WAVE MODEL GRID SPECIFICATION (INPUT):
!        *NGX*    - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY*    - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *KRGG*   - GRID DEFINITION PARAMETER (0=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (OUTPUT FIELD)
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELLO* - GRID POINT SPACING PER LATITUDES. 
!        *IFROMIJ*- INTEGER  !!! LOCAL !!! LONG. GRID INDEX.
!        *JFROMIJ*- INTEGER  !!! LOCAL !!! LAT. GRID INDEX (NORTH-SOUTH).
!        *KIJS:KIJL* DIMENSION OF IFROMIJ AND JFROMIJ
!        *KIJLMAX*- LARGEST WAM GRID POINT INDEX USED BY THREAD
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).           
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).           
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).           
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).           
!        *IPERIODIC* - SPECIFIES IF WAM GRID IS PERIODIC (=1) OR NOT (=0)
!        *PMISS*  - MISSING DATA VALUE
!        INTERPOLATION COEFFICIENTS & INDECES (INPUT):
!        *IJBLOCK* - BLOCK INDEX FOR THE ATMOSPHERIC FIELDS.
!        *LADEN*  - VARIABLE AIR DENSITY IS USED.
!        *ROAIR*  - DEFAULT VALUE FOR AIR DENSITY.
!        *LGUST*  - GUSTINESS IS USED.
!        *WSTAR0* - DEFAULT VALUE FOR w*. 
!        *LWCUR*  - SURFACE CURRENTS ARE USED.
!        *LLKC*   - GRAB LAKE COVER INFORMATION
!        *LLINTERPOL* - FLAG (TRUE=DO INTERPOLATION, FALSE=ASSIGN ONLY)
!        *NWX*    - FIRST DIMENSION FOR THE INTERPOLATION COEFFICIENTS 
!        *NWY*    - SECOND DIMENSION FOR THE INTERPOLATION COEFFICIENTS
!        *DJ1M*   - COEFFICIENT
!        *DII1M*  - COEFFICIENT
!        *DIIP1M* - COEFFICIENT
!        *JJM*    - INDEX
!        *IIM*    - INDEX
!        *IIPM*   - INDEX
!        *MASK_IN*  INTEGER  MASK TO INDICATE WHICH PART OF FIELDS IS RELEVANT.
!        *NXS:NXE*  FIRST DIMENSION OF FIELDG
!        *NYS:NYE*  SECOND DIMENSION OF FIELDG
!        *FIELDG*   INPUT FORCING FIELDS ON THE WAVE MODEL GRID



!     EXTERNALS.                                                        
!     ----------                                                        

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWWIND  , ONLY : LLNEWCURR

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, NGPTOTG, NC, NR, NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGX, NGY, KRGG, KIJS, KIJL, KIJLMAX
      INTEGER(KIND=JWIM), INTENT(IN) :: IPERIODIC
      INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(IN) :: KLONRGG, JJM
      INTEGER(KIND=JWIM), DIMENSION(NR), INTENT(IN) :: ILONRGG
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ, JFROMIJ
      INTEGER(KIND=JWIM), DIMENSION(0:NC+1,NR), INTENT(IN) :: IJBLOCK
      INTEGER(KIND=JWIM), INTENT(IN) :: NWX, NWY
      INTEGER(KIND=JWIM), DIMENSION(NWX,NWY), INTENT(IN) :: IIM, IIPM
      INTEGER(KIND=JWIM), DIMENSION(NGPTOTG), INTENT(INOUT) :: MASK_IN
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG


      REAL(KIND=JWRB), INTENT(IN) :: XDELLA, AMOWEP, AMOSOP, AMOEAP, AMONOP, PMISS
      REAL(KIND=JWRB), INTENT(IN) :: ROAIR, WSTAR0
      REAL(KIND=JWRB), DIMENSION(NGY), INTENT(IN) :: ZDELLO
      REAL(KIND=JWRB), DIMENSION(NWY), INTENT(IN) :: DJ1M
      REAL(KIND=JWRB), DIMENSION(NGPTOTG,NFIELDS), INTENT(IN) :: FIELDS
      REAL(KIND=JWRB), DIMENSION(NWX,NWY), INTENT(IN) :: DII1M, DIIP1M

      LOGICAL, INTENT(IN):: LADEN, LGUST, LWCUR, LLKC, LLINTERPOL


      INTEGER(KIND=JWIM) :: IJ, I, J
      INTEGER(KIND=JWIM) :: JJ, JSN, JJ1, JSN1, II, II1, IIP, IIP1, JCL, ICL
      INTEGER(KIND=JWIM) :: NCOUNT

      REAL(KIND=JWRB) :: DJ1, DJ2, DII1, DII2, DIIP1, DIIP2 
      REAL(KIND=JWRB) :: F00, F10, F01, F11, CI
      REAL(KIND=JWRB) :: ZLADEN, ZLGUST, ZLLKC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FLDINTER',0,ZHOOK_HANDLE)

      IF (LADEN) THEN
        ZLADEN = 1.0_JWRB
      ELSE
        ZLADEN = 0.0_JWRB
      ENDIF

      IF (LGUST) THEN
        ZLGUST=1.0_JWRB
      ELSE
        ZLGUST=0.0_JWRB
      ENDIF

      IF (LLKC) THEN
        ZLLKC=1.0_JWRB
      ELSE
        ZLLKC=0.0_JWRB
      ENDIF

      IF (.NOT.LLINTERPOL) THEN

!       REARRANGE DATA FIELD.
!       --------------------

          DO IJ = KIJS, KIJLMAX
            I = IFROMIJ(IJ)
            J = JFROMIJ(IJ)

            MASK_IN(IJBLOCK(I,J))=1

            FIELDG%UWND(I,J) = FIELDS(IJBLOCK(I,J),1)
            FIELDG%VWND(I,J) = FIELDS(IJBLOCK(I,J),2)
            FIELDG%AIRD(I,J) = ZLADEN*FIELDS(IJBLOCK(I,J),3) + (1.0_JWRB-ZLADEN)*ROAIR
            FIELDG%WSTAR(I,J) = ZLGUST*FIELDS(IJBLOCK(I,J),4) + (1.0_JWRB-ZLGUST)*WSTAR0
            FIELDG%CICOVER(I,J) = FIELDS(IJBLOCK(I,J),5)
            FIELDG%LKFR(I,J) = ZLLKC*FIELDS(IJBLOCK(I,J),6)
            FIELDG%USTRA(I,J) = FIELDS(IJBLOCK(I,J),7)
            FIELDG%VSTRA(I,J) = FIELDS(IJBLOCK(I,J),8)
!!!!!!!!!!! not yet in place to receive from IFS the sea ice thickness !!!!!!!!!!!
            FIELDG%CITHICK(I,J) = 0.0_JWRB
          ENDDO

          IF (LLNEWCURR) THEN
            IF (LWCUR) THEN
              DO IJ = KIJS, KIJLMAX
                I = IFROMIJ(IJ)
                J = JFROMIJ(IJ)
                FIELDG%UCUR(I,J) = FIELDS(IJBLOCK(I,J),9)
                FIELDG%VCUR(I,J) = FIELDS(IJBLOCK(I,J),10)
              ENDDO
            ELSE
              DO IJ = KIJS, KIJLMAX
                I = IFROMIJ(IJ)
                J = JFROMIJ(IJ)
                FIELDG%UCUR(I,J) = 0.0_JWRB
                FIELDG%VCUR(I,J) = 0.0_JWRB
              ENDDO
            ENDIF
          ENDIF

      ELSE

!       INTERPOLATE TO WAVE MODEL GRID
!       ------------------------------

          DO IJ = KIJS, KIJLMAX 
            I = IFROMIJ(IJ)
            J = JFROMIJ(IJ)


            JJ = JJM(J)
            JSN= NR-JJ+1
            JJ1= MIN(JJ+1,NR)
            JSN1=NR-JJ1+1
            DJ1= DJ1M(J) 
            DJ2=1.0_JWRB-DJ1

            II = IIM(I,J)
            II1 = MIN(II+1,ILONRGG(JSN)+IPERIODIC)
            DII1=DII1M(I,J) 
            DII2=1.0_JWRB-DII1

            IIP = IIPM(I,J) 
            IIP1 = MIN(IIP+1,ILONRGG(JSN1)+IPERIODIC)
            DIIP1= DIIP1M(I,J) 
            DIIP2=1.0_JWRB-DIIP1

            MASK_IN(IJBLOCK(II,JJ))=1
            MASK_IN(IJBLOCK(II1,JJ))=1
            MASK_IN(IJBLOCK(IIP,JJ1))=1
            MASK_IN(IJBLOCK(IIP1,JJ1))=1

            FIELDG%UWND(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),1) +      &
     &                      DII1*FIELDS(IJBLOCK(II1,JJ),1) ) +          &
     &                DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),1) +          &
     &                      DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),1) )

            FIELDG%VWND(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),2) +      &
     &                      DII1*FIELDS(IJBLOCK(II1,JJ),2) ) +          &
     &                DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),2) +          &
     &                      DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),2) )

            IF (LADEN) THEN
              FIELDG%AIRD(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),3) +    &
     &                        DII1*FIELDS(IJBLOCK(II1,JJ),3) ) +        &
     &                  DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),3) +        &
     &                        DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),3) )
            ELSE
              FIELDG%AIRD(I,J) = ROAIR
            ENDIF
            IF (LGUST) THEN
              FIELDG%WSTAR(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),4) +   &
     &                         DII1*FIELDS(IJBLOCK(II1,JJ),4) ) +       &
     &                   DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),4) +       &
     &                         DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),4) )
            ELSE
              FIELDG%WSTAR(I,J) = WSTAR0 
            ENDIF

!           FOR SEA ICE FRACTION
!           DETERMINE WHETHER ANY OF THE 4 CONERS HAS A MISSING DATA
!           IF SO USE THE CLOSEST GRID POINT VALUE IF IT IS NOT MISSING
            F00=FIELDS(IJBLOCK(II,JJ),5)
            F10=FIELDS(IJBLOCK(II1,JJ),5)
            F01=FIELDS(IJBLOCK(IIP,JJ1),5)
            F11=FIELDS(IJBLOCK(IIP1,JJ1),5)
            IF (F00 == PMISS .OR.                                       &
     &          F10 == PMISS .OR.                                       &
     &          F01 == PMISS .OR.                                       &
     &          F11 == PMISS     ) THEN
              IF (DJ1 <= 0.5_JWRB) THEN
                JCL=JJ
                ICL=II1
                IF (DII1 <= 0.5_JWRB) ICL=II
              ELSE
                JCL=JJ1
                ICL=IIP1
                IF (DIIP1 <= 0.5_JWRB) ICL=IIP
              ENDIF
              CI=FIELDS(IJBLOCK(ICL,JCL),5)

!             NON MISSING VALUE OVER SEA IS NEEDED
              IF (CI == PMISS) THEN
                NCOUNT=0
                CI=0.
                IF (F00 /= PMISS) THEN
                  CI=CI+F00
                  NCOUNT=NCOUNT+1
                ENDIF
                IF (F10 /= PMISS) THEN
                  CI=CI+F10
                  NCOUNT=NCOUNT+1
                ENDIF
                IF (F01 /= PMISS) THEN
                  CI=CI+F01
                  NCOUNT=NCOUNT+1
                ENDIF
                IF (F11 /= PMISS) THEN
                  CI=CI+F11
                  NCOUNT=NCOUNT+1
                ENDIF
                IF (NCOUNT > 0) THEN
                  CI=CI/NCOUNT
                ELSE
                  CI=PMISS
                ENDIF
              ENDIF
              FIELDG%CICOVER(I,J)=CI
            ELSE
              FIELDG%CICOVER(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),5) + &
     &                      DII1*FIELDS(IJBLOCK(II1,JJ),5) ) +          &
     &                DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),5) +          &
     &                      DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),5) )
            ENDIF

!!!!!!!!!!! not yet in place to receive from IFS the sea ice thickness !!!!!!!!!!!
            FIELDG%CITHICK(I,J) = 0.0_JWRB


            IF (LLKC) THEN
              FIELDG%LKFR(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),6) +    &
     &                         DII1*FIELDS(IJBLOCK(II1,JJ),6) ) +       &
     &                   DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),6) +       &
     &                         DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),6) )
            ELSE
              FIELDG%LKFR(I,J) = 0.0_JWRB
            ENDIF

            FIELDG%USTRA(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),7) +  &
     &                        DII1*FIELDS(IJBLOCK(II1,JJ),7) ) +     &
     &                  DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),7) +     &
     &                        DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),7) )

            FIELDG%VSTRA(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),8) +  &
     &                        DII1*FIELDS(IJBLOCK(II1,JJ),8) ) +     &
     &                  DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),8) +     &
     &                        DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),8) )

            IF (LLNEWCURR) THEN
              IF (LWCUR) THEN
                FIELDG%UCUR(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),9) +  &
     &                           DII1*FIELDS(IJBLOCK(II1,JJ),9) ) +     &
     &                     DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),9) +     &
     &                           DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),9) )

                FIELDG%VCUR(I,J)=DJ2*( DII2*FIELDS(IJBLOCK(II,JJ),10) +  &
     &                           DII1*FIELDS(IJBLOCK(II1,JJ),10) ) +     &
     &                     DJ1*( DIIP2*FIELDS(IJBLOCK(IIP,JJ1),10) +     &
     &                           DIIP1*FIELDS(IJBLOCK(IIP1,JJ1),10) )
              ELSE
                FIELDG%UCUR(I,J) = 0.0_JWRB
                FIELDG%VCUR(I,J) = 0.0_JWRB
              ENDIF
            ENDIF

          ENDDO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('FLDINTER',1,ZHOOK_HANDLE)

      END SUBROUTINE FLDINTER
