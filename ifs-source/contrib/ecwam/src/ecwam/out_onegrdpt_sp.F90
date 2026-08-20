! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUT_ONEGRDPT_SP(SPEC,USTAR,CDTPRO)
 
!--------------------------------------------------------------------
 
!*****OUT_ONEGRDPT_SP** OUTPUT OF A NUMBER OF INTEGRATED PARAMETERS AND
!                       OF THE ONE AND TWO-DIMENSIONAL SPECTRUM AT
!                       ONE GRID POINT.
!
!     P.JANSSEN JUNE 2005
 
!     PURPOSE
!     -------
!             TO OUTPUT INFORMATIONS AT A FIXED GRID POINT.
 
!     INTERFACE
!     ---------
!             *CALL* *OUT_ONEGRDPT_SP*
 
!           *SPEC*   -   TWO-DIMENSIONAL SPECTRA.
!           *USTAR*  -   FRICTION VELOCITY
!           *CDTPRO* -   DATE TIME GROUP
!--------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
 
      USE YOWFRED  , ONLY : FR       ,DELTH    ,TH       ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE     
      USE YOWPCONS , ONLY : G        ,ZPI      ,DEG      ,EPSMIN

! ----------------------------------------------------------------------
      IMPLICIT NONE
 
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN) :: SPEC
      REAL(KIND=JWRB), INTENT(IN) :: USTAR

      CHARACTER(LEN=12), INTENT(IN) :: CDTPRO

      INTEGER(KIND=JWIM) ::  K, M
      INTEGER(KIND=JWIM) :: IU_SPEC
      INTEGER(KIND=JWIM) :: MMAX
      REAL(KIND=JWRB) :: SUMMAX 
      REAL(KIND=JWRB) :: FPK, X1, X2, X3, X4, X5, X6, AF, ASPF, AF2, ASPF2
      REAL(KIND=JWRB), DIMENSION(NFRE) :: SPF

      CHARACTER(LEN=72) :: CPATH

      LOGICAL, SAVE :: FRSTIME

      DATA FRSTIME/.TRUE./   
 
! -------------------------------------------------------------------
 
!       THE SWAMP CASE HAS BEEN REDUCED TO A ONE GRID POINT MODEL
!       NECESSARY PARAMETERS ARE OUTPUT HERE

      IU_SPEC = 998

      IF (FRSTIME) THEN
        CPATH = 'out_spec'
        OPEN (IU_SPEC,FILE=CPATH,STATUS='UNKNOWN', FORM='FORMATTED')
        FRSTIME = .FALSE.
      ENDIF
         
      DO M=1,NFRE
        SPF(M) = 0._JWRB
        DO K=1,NANG 
          SPF(M) = SPF(M) + SPEC(K,M)*DELTH
        ENDDO
      ENDDO 

      SUMMAX = 0._JWRB
      DO M=1,NFRE
        IF(SPF(M).GT.SUMMAX)THEN
          SUMMAX = SPF(M)
          MMAX = M
        END IF
      ENDDO

      IF(SUMMAX.EQ.0._JWRB)THEN
         FPK=0._JWRB
      ELSEIF(MMAX.EQ.1) THEN
         FPK=FR(MMAX)
      ELSEIF(MMAX.EQ.NFRE) THEN
         FPK=FR(MMAX)
      ELSE
         X1 = FR(MMAX)**2-FR(MMAX+1)**2
         X2 = FR(MMAX)-FR(MMAX+1)
         X3 = SPF(MMAX)-SPF(MMAX+1)
         X4 = FR(MMAX-1)**2-FR(MMAX+1)**2
         X5 = FR(MMAX-1)-FR(MMAX+1)
         X6 = SPF(MMAX-1)-SPF(MMAX+1)
         FPK =0.5_JWRB*X1*(X3*X4-X1*X6)/(X1*X3*X5-X1*X2*X6)
      END IF 

      DO M=1,NFRE
         AF = FR(M)/FPK
         ASPF = (ZPI*FPK)**3/(G*USTAR)*SPF(M)
         AF2 = USTAR*FR(M)/G
         ASPF2 = ZPI**3*G**3/USTAR**5*SPF(M)
!!         WRITE(IU_SPEC,62) CDTPRO,FR(M),SPF(M),AF,ASPF,AF2,ASPF2
         WRITE(IU_SPEC,63) CDTPRO,FR(M),SPF(M)
      ENDDO

 62   FORMAT(A12,1x,2(1X,F8.4),1X,F10.6,E10.3)
 63   FORMAT(A12,1x,2(1x,F15.10))
 
! -------------------------------------------------------------------
 
      END SUBROUTINE OUT_ONEGRDPT_SP
