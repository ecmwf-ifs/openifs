! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE MPEXCHNGSARIN(NSPECPE, ISPECPE, ITAG) 

!****  *MPEXCHNGSARIN* - EXCHANGES MODEL DATA AT SAR DATA POINTS. 

!     J. BIDLOT    ECMWF   NOVEMBER 1999

!     PURPOSE.
!     --------

!     BROADCASTS MODEL DATA AT SAR DATA POINTS FROM ALL PE'S TO
!     ALL PE'S 

!*    INTERFACE.
!     ----------

!     CALL *MPEXCHNGSARIN*(NSPECPE,ISPECPE,ITAG)


!     *NSPECPE* - NUMBER OF SAR DATA POINT ON THAT PE.
!     *ISPECPE* - ARRAY INDEX OF EACH POINT BELONGING TO THAT PE.
!     *ITAG*    - TAG USED FOR MESSAGE PASSING.  
!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSARAS , ONLY : NSPEC    ,SPEC     ,IJSAR    ,              &
     &            LONG     ,LAT      ,U10      ,USSAR    ,THW
      USE MPL_MODULE, ONLY : MPL_BROADCAST

!----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: NSPECPE, KSPECPE, IP, ISPEC, IFRE, IANG
      INTEGER(KIND=JWIM) :: ITAG, MZC, KCOUNT
      INTEGER(KIND=JWIM) :: ISPECPE(NSPEC)

      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUF(:)

!----------------------------------------------------------------------

      IF (NPROC.EQ.1) THEN
        RETURN
      ELSE
        MZC=1+NSPEC*(8+NANG*NFRE)
        ALLOCATE(ZCOMBUF(MZC))

        DO IP=1,NPROC

          IF (IRANK.EQ.IP) THEN
            KCOUNT=1
            ZCOMBUF(KCOUNT)=FLOAT(NSPECPE)
            DO ISPEC = 1,NSPECPE
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=FLOAT(ISPECPE(ISPEC))
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=FLOAT(IJSAR(ISPECPE(ISPEC)))
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=LONG(ISPECPE(ISPEC),1) 
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=LAT(ISPECPE(ISPEC),1) 
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=U10(ISPECPE(ISPEC)) 
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=USSAR(ISPECPE(ISPEC)) 
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=THW(ISPECPE(ISPEC)) 
              DO IFRE = 1,NFRE
                DO IANG = 1,NANG
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=SPEC(ISPECPE(ISPEC),IANG,IFRE,1)
                ENDDO 
              ENDDO 
            ENDDO 
          ENDIF

          CALL MPL_BROADCAST(ZCOMBUF(1:MZC),KROOT=IP,KTAG=ITAG,         &
     &       CDSTRING='MPEXCHNGSARIN:')

          ITAG=ITAG+1

          IF (IRANK.NE.IP) THEN
            KCOUNT=1
            DO ISPEC = 1,NINT(ZCOMBUF(1))
              KCOUNT=KCOUNT+1
              KSPECPE=NINT(ZCOMBUF(KCOUNT))
              KCOUNT=KCOUNT+1
              IJSAR(KSPECPE)=NINT(ZCOMBUF(KCOUNT))
              KCOUNT=KCOUNT+1
              LONG(KSPECPE,1)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              LAT(KSPECPE,1)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              U10(KSPECPE)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              USSAR(KSPECPE)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              THW(KSPECPE)=ZCOMBUF(KCOUNT)
              DO IFRE = 1,NFRE
                DO IANG = 1,NANG
                  KCOUNT=KCOUNT+1
                  SPEC(KSPECPE,IANG,IFRE,1)=ZCOMBUF(KCOUNT)
                ENDDO 
              ENDDO 
            ENDDO 
          ENDIF
        
        ENDDO 

        DEALLOCATE(ZCOMBUF)

      ENDIF

      END SUBROUTINE MPEXCHNGSARIN
