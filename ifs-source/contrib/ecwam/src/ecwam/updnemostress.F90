! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE UPDNEMOSTRESS


!****  *UPDNEMOSTRESS* - UPDATE FIELDS PASSED TO NEMO WITH WAVE INFORMATION

!      KRISTIAN MOGENSEN ECMWF    JULY 2013                    

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES WAVE INFORMATION THROUGH TO
!          NEMO VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!          !!!!! APPLIES TO ACCUMULATED FIELDS !!!!!
!          !!!!! SEE UPDNEMOFILDS FOR NON-ACCUMULATED FIELDS !!!!!

!*    INTERFACE.                                                        
!     ----------                                                        


!          THE NEMO COUPLING IN nemo/coupled/src/nemointerface

!     METHOD.                                                           
!     -------                                                           

!          PARALLEL INTERPOLATION BASED ON PREDETERMINED WEIGHTS

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_TAU_UPDATE  -  UPDATE WAM STRESS IN NEMO

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, NTOTIJ, KIJL4CHNK
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! DATA FOR NEMO 
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOCOUDEBUG, NEMONTAU
      USE YOWNEMOFLDS, ONLY : WAM2NEMO
      USE YOWSTAT  , ONLY : CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: ICHNK, IC, KIJS, KIJL

      REAL(KIND=JWRB) :: ZNEMONTAUM1
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRO), DIMENSION(NTOTIJ) :: NEMOTAUX, NEMOTAUY, NEMOWSWAVE, NEMOPHIF

! -------------------------------------------------------------------   
IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',0,ZHOOK_HANDLE)


      IF (LWNEMOCOU .AND. ALLOCATED(WAM2NEMO%NEMOTAUX) ) THEN

        IF (NEMONTAU > 0) THEN
          ZNEMONTAUM1= 1.0_JWRB/NEMONTAU

!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IC, KIJL)
          DO ICHNK = 1, NCHNK

            KIJS = 1
            DO IC = 1, ICHNK-1
               KIJS = KIJS + KIJL4CHNK(IC)
            ENDDO
            KIJL = KIJS + KIJL4CHNK(ICHNK) - 1

            NEMOTAUX(KIJS:KIJL)   = WAM2NEMO%NEMOTAUX(1:KIJL4CHNK(ICHNK), ICHNK) * ZNEMONTAUM1
            NEMOTAUY(KIJS:KIJL)   = WAM2NEMO%NEMOTAUY(1:KIJL4CHNK(ICHNK), ICHNK) * ZNEMONTAUM1
            NEMOWSWAVE(KIJS:KIJL) = WAM2NEMO%NEMOWSWAVE(1:KIJL4CHNK(ICHNK), ICHNK) * ZNEMONTAUM1
            NEMOPHIF(KIJS:KIJL)   = WAM2NEMO%NEMOPHIF(1:KIJL4CHNK(ICHNK), ICHNK) * ZNEMONTAUM1

          ENDDO
!$OMP     END PARALLEL DO

        ELSE
          NEMOTAUX(1:NTOTIJ) = 0.0_JWRO
          NEMOTAUY(1:NTOTIJ) = 0.0_JWRO
          NEMOWSWAVE(1:NTOTIJ) = 0.0_JWRO
          NEMOPHIF(1:NTOTIJ) = 0.0_JWRO
        ENDIF

#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE_STRESS( IRANK-1, NPROC, MPL_COMM,   &
     &      NTOTIJ, NEMOTAUX, NEMOTAUY, NEMOWSWAVE, NEMOPHIF, &
     &      CDTPRO, LWNEMOCOUDEBUG )
#endif
        ! INITIALIZE STRESS FOR ACCUMULATION

        NEMONTAU = 0
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)
        DO ICHNK = 1, NCHNK
          WAM2NEMO%NEMOTAUX(1:NPROMA_WAM,ICHNK) = 0.0_JWRO
          WAM2NEMO%NEMOTAUY(1:NPROMA_WAM,ICHNK) = 0.0_JWRO
          WAM2NEMO%NEMOWSWAVE(1:NPROMA_WAM,ICHNK) = 0.0_JWRO
          WAM2NEMO%NEMOPHIF(1:NPROMA_WAM,ICHNK) = 0.0_JWRO
        ENDDO
!$OMP   END PARALLEL DO

      ENDIF

IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',1,ZHOOK_HANDLE)

END SUBROUTINE UPDNEMOSTRESS
