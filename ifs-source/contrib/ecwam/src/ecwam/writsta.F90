! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WRITSTA (IDIN, CBDT, CEDT, ANALPD, FOREPD, IS,         &
     &                    CRSDT, CLSDT, CABDT, CAEDT, IASS, NF, ISTAT,  &
     &                    CCURA, LRSTPARAL, NPROC_RST )

                                                                        
! --------------------------------------------------------------------- 

!***      *WRITSTA* - WRITES A NEW WAMINFO FILE                         

!          P. JANSSEN                                                   
!          B. HANSEN APRIL 97  : add CRSDT and CLSDT definition
!                                move COMMON CARD INTO INCLUDE FILE.

!     PURPOSE.                                                          
!     --------                                                          

!          WRITES A NEW WAMINFO FILE TO UNIT IDIN                       

!     INTERFACE.                                                        
!     ----------                                                        

!         *IDIN*     INTEGER     UNIT NUMBER.                           
!         *CBDT*     CHAR*14     BEGIN DATE TIME GROUP.                 
!         *CEDT*     CHAR*14     END DATE TIME GROUP.                   
!         *ANALPD*   INTEGER     ANALYSIS PERIOD IN SECONDS.
!         *FOREPD*   INTEGER     FORECAST PERIOD IN SECONDS.
!         *IS*       INTEGER     WIND INPUT TIME STEP IN SECONDS.
!         *CRSDT*    CHAR*14     DATE FOR RESTART FILES OUTPUT
!         *CLSDT*    CHAR*14     LAST DATE FOR SPECTRUM FILE OUTPUT
!         *CABDT*    CHAR*14     BEGIN ANALYSIS DATE TIME GROUP.            
!         *CAEDT*    CHAR*14     END ANALYSIS DATE TIME GROUP.              
!         *IASS*     INTEGER     ASSIMILATION CONDITION.                
!         *NF*       INTEGER     NEW FORECAST CONDITION.                
!         *ISTAT*    INTEGER     JOB STATUS INDICATORS.                 
!         *CCURA*    CHAR*14     BEGIN DATE FOR USE OF CURRENTS.                 
!         *LRSTPARAL LOGICAL     TRUE IF RESTART FILES WRITTEN IN PARALLEL
!         *NPROC_RST INTEGER     NUMBER OF PARALLEL MPI TASKS THAT HAVE BEEN USED


!     METHOD.                                                           
!     -------                                                           

!          NEW WAMINFO FILE IS WRITTEN TO UNIT IDNIN                    

!     EXTERNALS.                                                        
!     ----------                                                        

!          NONE                                                         

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

!---------------------------------------------------------------------- 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCARD  , ONLY : JPCL     ,CARD

! --------------------------------------------------------------------- 

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IDIN, IS, IASS, NF 
      INTEGER(KIND=JWIM), INTENT(IN) :: ANALPD, FOREPD
      INTEGER(KIND=JWIM), INTENT(IN) :: ISTAT(3)
      INTEGER(KIND=JWIM), INTENT(IN) :: NPROC_RST

      CHARACTER(LEN=14), INTENT(IN) :: CBDT,CEDT,CABDT,CAEDT,CRSDT,CLSDT
      CHARACTER(LEN=14), INTENT(IN) :: CCURA

      LOGICAL, INTENT(IN) :: LRSTPARAL

      INTEGER(KIND=JWIM) :: I
! --------------------------------------------------------------------- 

!*    1.0  FILL CARD ARRAY.                                             
!          ----------------                                             

      WRITE(CARD(1), 1) CBDT, CEDT
    1 FORMAT('RUN MODEL FROM ',A14,' TO ',A14)
                                                                        
      CARD(2)='PARAM=10U/10V,'

      CARD(3)='REPRES=GG,'

      WRITE(CARD(4), 2) ANALPD                                          
    2 FORMAT('ANALYSIS PERIOD = ',I7)        
                                                                        
      WRITE(CARD(5), 3) FOREPD                                          
    3 FORMAT('FORECAST PERIOD = ',I10)         
                                                                        
      WRITE(CARD(6), 4) IS                                              
    4 FORMAT('WIND TIME STEP IN SECONDS = ',I7)
                                                                        
      WRITE(CARD(7), 5) CABDT, CAEDT
    5 FORMAT('ANALYSIS FROM ',A14,' TO ',A14)
                                                                        
      IF (IASS == 0) THEN                                               
        CARD(8) = 'ASSIMILATION NO'        
      ELSE                                                              
        CARD(8) = 'ASSIMILATION YES'      
      ENDIF                                                            
                                                                        
      IF (NF == 0) THEN                                                 
        CARD(9) = 'NEW FORECAST NO'      
      ELSE                                                              
        CARD(9) = 'NEW FORECAST YES'    
      ENDIF                                                            
                                                                        
      DO I=1,3                                                     
        IF (ISTAT(I) == 0) THEN                                        
          WRITE(CARD(9+I), 6) I                                       
    6     FORMAT('STATUS STORM',I1,'= UNFINISHED')                    
        ELSE                                                           
          WRITE(CARD(9+I), 7) I                                       
    7     FORMAT('STATUS STORM',I1,'= FINISHED')                      
        ENDIF                                                         
      ENDDO
      WRITE(CARD(13),8) CRSDT
    8 FORMAT('DATE FOR OUTPUT OF BOTH RESTART FILES = ',A14)

      WRITE(CARD(14),9) CLSDT
    9 FORMAT('LAST DATE FOR SPECTRA FILE OUTPUT = ',A14)

      WRITE(CARD(15),10) CCURA
   10 FORMAT('BEGIN DATE FOR USING SURFACE CURRENT = ',A14)

      IF (LRSTPARAL) THEN 
        CARD(16) = 'PARALLEL RESTART = YES'        
      ELSE
        CARD(16) = 'PARALLEL RESTART = NO'        
      ENDIF

      WRITE(CARD(17),11) NPROC_RST
   11 FORMAT('NUMBER OF MPI TASKS USED = ',I10)

!*    2.0  WRITE CARD TO IDIN.                                          
!          -------------------                                          

      DO I=1,JPCL
        WRITE(IDIN,'(A72)') CARD(I)                                    
      ENDDO

      END SUBROUTINE WRITSTA
