! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE UPDNEMOFIELDS


!****  *UPDNEMOFIELDS* - UPDATE FIELDS IN NEMO WITH WAVE INFORMATION

!      KRISTIAN MOGENSEN ECMWF    NOVEMBER 2012                    

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES WAVE INFORMATION THROUGH TO
!          NEMO VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!          !!!!! APPLIES TO NON-ACCUMULATED FIELDS !!!!!
!          !!!!! SEE UPDNEMOSTRESS FOR ACCUMULATED FIELDS !!!!!

!*    INTERFACE.                                                        
!     ----------                                                        

!          THE NEMO COUPLING IN nemo/coupled/src/nemointerface

!     METHOD.                                                           
!     -------                                                           

!          PARALLEL INTERPOLATION BASED ON PREDETERMINED WEIGHTS

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_UPDATE  -  UPDATE WAM FIELDS IN NEMO

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
! COUPLING INFORMATION
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOCOUDEBUG
      USE YOWNEMOFLDS, ONLY : WAM2NEMO 
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! INFORMATION FOR OPTIONAL DEBUGGING
      USE YOWSTAT  , ONLY : CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: ICHNK, IC, KIJS, KIJL
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRO), DIMENSION(NTOTIJ) :: NSWH, NMWP, NPHIEPS, NTAUOC
      REAL(KIND=JWRO), DIMENSION(NTOTIJ) :: NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES

! -------------------------------------------------------------------   

IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',0,ZHOOK_HANDLE)

      IF (LWNEMOCOU) THEN

!$OMP  PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IC, KIJL)
       DO ICHNK = 1, NCHNK

          KIJS = 1
          DO IC = 1, ICHNK-1
             KIJS = KIJS + KIJL4CHNK(IC)
          ENDDO
          KIJL = KIJS + KIJL4CHNK(ICHNK) - 1

          NSWH(KIJS:KIJL)        = WAM2NEMO%NSWH(1:KIJL4CHNK(ICHNK), ICHNK)
          NMWP(KIJS:KIJL)        = WAM2NEMO%NMWP(1:KIJL4CHNK(ICHNK), ICHNK)
          NPHIEPS(KIJS:KIJL)     = WAM2NEMO%NPHIEPS(1:KIJL4CHNK(ICHNK), ICHNK)
          NTAUOC(KIJS:KIJL)      = WAM2NEMO%NTAUOC(1:KIJL4CHNK(ICHNK), ICHNK)
          NEMOSTRN(KIJS:KIJL)    = WAM2NEMO%NEMOSTRN(1:KIJL4CHNK(ICHNK), ICHNK)
          NEMOUSTOKES(KIJS:KIJL) = WAM2NEMO%NEMOUSTOKES(1:KIJL4CHNK(ICHNK), ICHNK)
          NEMOVSTOKES(KIJS:KIJL) = WAM2NEMO%NEMOVSTOKES(1:KIJL4CHNK(ICHNK), ICHNK)

       ENDDO
!$OMP  END PARALLEL DO

#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE( IRANK-1, NPROC, MPL_COMM,             &
     &                               NTOTIJ,                               &
     &                               NSWH, NMWP, NPHIEPS, NTAUOC,          &
     &                               NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES,   &
     &                               CDTPRO, LWNEMOCOUDEBUG )
#endif

      ENDIF

IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',1,ZHOOK_HANDLE)

END SUBROUTINE UPDNEMOFIELDS
