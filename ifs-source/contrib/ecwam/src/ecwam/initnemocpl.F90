! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE INITNEMOCPL(BLK2LOC) 

!****  *INITNEMOCPL* - INITIALIZE THE COUPLING OF WAM TO NEMO

!      KRISTIAN MOGENSEN ECMWF    NOVEMBER 2012                    

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PREPARES FOR THE COUPLING OF WAM
!          TO NEMO BY DEFINING A GLOBAL GRID INDEX NEEDED TO
!          CALL NEMOGCMCOUP_WAM_COUPINIT         

!*    INTERFACE.                                                        
!     ----------                                                        

!          THE NEMO COUPLING IN nemo/coupled/src/nemointerface

!     METHOD.                                                           
!     -------                                                           

!          CREATE A GLOBAL INDEX AND NEMOGCMCOUP_WAM_COUPINIT

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_COUPINIT  -  INITIALIZES INTERNAL
!                      DATA STRUCTURES FOR THE COUPLING

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

! MODULES NEED FOR KIND DEFINITION
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC

! GRID POINTS CHUNKS
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, NTOTIJ, KIJL4CHNK
! MODULES NEED FOR GRID DEFINTION      
      USE YOWMAP   , ONLY : NGX, NGY, NLONRGG
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! WAVE FIELDS TO NEMO 
      USE YOWCOUP  , ONLY :  NEMONTAU
! FOR PRINTING
      USE YOWTEST  , ONLY : IU06
! DR. HOOK
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! NULERR
      USE EC_LUN   , ONLY : NULERR

! -------------------------------------------------------------------   

      IMPLICIT NONE
#include "abort1.intfb.h"

! LOCAL WAM GRID POINTS
! POINTERS FROM LOCAL GRID POINTS TO 2-D MAPS
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
! GLOBAL INDEX DEFINITON.      
      INTEGER(KIND=JWIM) :: IGLOBAL(NGX, NGY)
! INDICES
      INTEGER(KIND=JWIM) :: IJ, IX, JSN, I, J, K, IC, ICHNK
! TEMPORARY ARRAYS FOR LOCAL MASK AND GLOBAL INDICES
      INTEGER(KIND=JWIM), DIMENSION(NTOTIJ) :: NLOCMSK, NGLOIND
! MISC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------   

IF (LHOOK) CALL DR_HOOK('INITNEMOCPL',0,ZHOOK_HANDLE)

      WRITE(IU06,*)
      WRITE(IU06,*)' **************************************************'
      WRITE(IU06,*)' * INITIALIZING THE WAM <-> NEMO COUPLING         *'
      WRITE(IU06,*)' **************************************************'
      WRITE(IU06,*)
      CALL FLUSH(IU06)

! CONSTRUCT GLOBAL INDICES GRID.
! THE CONVENTION IS NORTH TO SOUTH AND WEST TO EAST 
! STARTING AT 0 DEG LONGITUDE

      IGLOBAL(:,:)=-1
      K=0
      DO J=1,NGY
         DO I=1,NLONRGG(NGY-J+1)
            K=K+1
            IGLOBAL(I,NGY-J+1)=K
         ENDDO
      ENDDO

! SETUP GLOBAL INDICES 1D ARRAY + LOCAL MASK
! SINCE WE ONLY HAVE OCEAN POINTS THE LOCAL MASK IS ONE FOR ALL POINTS

      IC = 0
      DO ICHNK = 1, NCHNK
         DO IJ = 1, KIJL4CHNK(ICHNK)
           IX   = BLK2LOC%IFROMIJ(IJ,ICHNK)
           JSN  = BLK2LOC%KFROMIJ(IJ,ICHNK)
           IC = IC + 1
           NLOCMSK(IC) = 1
           NGLOIND(IC) = IGLOBAL(IX,JSN)
        ENDDO
      ENDDO

! CHECK FOR GLOBAL INDICES LESS THAN ZERO
      IF (MINVAL(NGLOIND) < 1) THEN
         WRITE(NULERR,*)'INITNEMOCPL : GLOBAL INDEX LESS THAN 1 FOUND!!!'
         IC = 0
         DO ICHNK = 1, NCHNK
           DO IJ = 1, KIJL4CHNK(ICHNK)
             IX   = BLK2LOC%IFROMIJ(IJ,ICHNK)
             JSN  = BLK2LOC%KFROMIJ(IJ,ICHNK)
             IC = IC + 1
             WRITE(NULERR,*) ICHNK, IJ, IX, JSN, NLONRGG(JSN), NGLOIND(IC)
           ENDDO
         ENDDO
         CALL ABORT1
      ENDIF         

! FOR ACCUMULATED FIELDS (WAM to NEMO)
      NEMONTAU     = 0

! CALL COUPLING INTERFACE.
! THE TASKS RANK CONVENTION IS 0 TO NPROC-1 SO WE SUBTRACT ONE FOR IRANK
! SINCE WE ONLY HAVE LOCAL OCEAN POINTS WE SPECIFY GLOBAL NUMBER OF
! POINTS AS THE SUM OF NLOGRGG

#ifdef WITH_NEMO
      CALL NEMOGCMCOUP_WAM_COUPINIT( IRANK-1, NPROC, MPL_COMM,    &
     &                               NTOTIJ, SUM(NLONRGG(1:NGY)), &
     &                               NLOCMSK, NGLOIND, 0)
#endif


IF (LHOOK) CALL DR_HOOK('INITNEMOCPL',1,ZHOOK_HANDLE)

END SUBROUTINE INITNEMOCPL
