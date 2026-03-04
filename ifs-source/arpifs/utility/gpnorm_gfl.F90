! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE GPNORM_GFL(YDGEOMETRY,YDGFL,LDPRINT_TL)

!**** *GPNORM_GFL*  - ROUTINE FOR CALCULATING GRID POINT NORM OF
!                     GFL VARIABLES

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        *CALL* *GPNORM_GFL*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!       G. Mozdzynski *ECMWF*  (based on gpnorm3 + gpnorm1)
!       Original : 02 Oct 2008

!     MODIFICATIONS.
!     --------------
!     A.Bogatchev: 12 Jun 2009 , call to egpnorm_trans
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     -----------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL       , ONLY : TGFL
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : NPROC, MYPROC
USE YOMCT0       , ONLY : NSPPR, LELAM
USE YOMLUN       , ONLY : NULOUT

IMPLICIT NONE

TYPE(GEOMETRY)  , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)      , INTENT(IN)    :: YDGFL
LOGICAL,OPTIONAL, INTENT(IN)    :: LDPRINT_TL
REAL(KIND=JPRB) , ALLOCATABLE   :: ZGP(:,:,:)
REAL(KIND=JPRB) , ALLOCATABLE   :: ZAVE(:), ZMAX(:), ZMIN(:)
REAL(KIND=JPRB) , ALLOCATABLE   :: ZMAXBLK(:,:), ZMINBLK(:,:)


INTEGER(KIND=JPIM) :: IEND, IST, JF, JL, JKGLO, IBL, JGFL
INTEGER(KIND=JPIM) :: IFLD,IFIELDS
LOGICAL            :: LLEVELS, LLAVE_ONLY, LLNOT_TL_AD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gpnorm_trans.h"
#include "egpnorm_trans.h"

#include "gpnorm1.intfb.h"


!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPNORM_GFL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YGFL=>YDGFL%YGFL)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GFL=>YDGFL%GFL)
LLNOT_TL_AD=.NOT.PRESENT(LDPRINT_TL)
LLEVELS = NSPPR>0

IF( LLEVELS )THEN

  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LGP) THEN
      ALLOCATE(ZGP(NPROMA,NFLEVG,NGPBLKS))

      CALL GSTATS(1223,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IST,IEND,IBL,JF,JL)
      DO JKGLO=1,NGPTOT,NPROMA
        IST=1
        IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        DO JF=1,NFLEVG
          DO JL=IST,IEND
            ZGP(JL,JF,IBL) = GFL(JL,JF,JGFL,IBL)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1223,1)

      CALL GPNORM1(YDGEOMETRY,ZGP,NFLEVG,LLEVELS,YCOMP(JGFL)%CNAME)

      DEALLOCATE(ZGP)
    ENDIF
  ENDDO

ELSE

  IFLD=0
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LGP)THEN
      IF(LLNOT_TL_AD&
                              & .OR.(YCOMP(JGFL)%CNAME==' LIQUID WATER'.OR.&
                              & YCOMP(JGFL)%CNAME==' ICE WATER'.OR.&
                              & YCOMP(JGFL)%CNAME==' CLOUD FRACTION'.OR.&
                              &      YCOMP(JGFL)%CNAME==' HUMIDITY')) THEN
        IFLD=IFLD+1
      ENDIF
    ENDIF
  ENDDO
  IFIELDS=IFLD

  ALLOCATE(ZGP(NPROMA,IFIELDS,NGPBLKS))

  ALLOCATE(ZAVE(IFIELDS))
  ALLOCATE(ZMIN(IFIELDS))
  ALLOCATE(ZMAX(IFIELDS))
  ALLOCATE(ZMINBLK(IFIELDS,NGPBLKS))
  ALLOCATE(ZMAXBLK(IFIELDS,NGPBLKS))
  
  IFLD=0
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LGP)THEN
      IF(LLNOT_TL_AD&
                              & .OR.(YCOMP(JGFL)%CNAME==' LIQUID WATER'.OR.&
                              & YCOMP(JGFL)%CNAME==' ICE WATER'.OR.&
                              & YCOMP(JGFL)%CNAME==' CLOUD FRACTION'.OR.&
                              &      YCOMP(JGFL)%CNAME==' HUMIDITY')) THEN
        IFLD=IFLD+1

        ZMINBLK(IFLD,:)=GFL(1,1,JGFL,1)
        ZMAXBLK(IFLD,:)=GFL(1,1,JGFL,1)

        CALL GSTATS(1427,0)
!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IST,IEND,IBL,JF,JL)
        DO JKGLO=1,NGPTOT,NPROMA
          IST=1
          IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
          IBL=(JKGLO-1)/NPROMA+1
          DO JL=IST,IEND
            ZGP(JL,IFLD,IBL) = 0.0_JPRB
          ENDDO
          DO JF=1,NFLEVG
            DO JL=IST,IEND
              ZGP(JL,IFLD,IBL) = ZGP(JL,IFLD,IBL)+GFL(JL,JF,JGFL,IBL)
              ZMINBLK(IFLD,IBL) = MIN(ZMINBLK(IFLD,IBL),GFL(JL,JF,JGFL,IBL))
              ZMAXBLK(IFLD,IBL) = MAX(ZMAXBLK(IFLD,IBL),GFL(JL,JF,JGFL,IBL))
            ENDDO
          ENDDO
        ENDDO
!!$OMP END PARALLEL DO
        CALL GSTATS(1427,1)
        ZMIN(IFLD)=MINVAL(ZMINBLK(IFLD,1:NGPBLKS))
        ZMAX(IFLD)=MAXVAL(ZMAXBLK(IFLD,1:NGPBLKS))
      ENDIF
    ENDIF
  ENDDO

  DEALLOCATE(ZMINBLK)
  DEALLOCATE(ZMAXBLK)

  LLAVE_ONLY=.TRUE. 
  IF ( .NOT. LELAM ) THEN
   CALL GPNORM_TRANS(ZGP,IFIELDS,NPROMA,ZAVE,ZMIN,ZMAX,LLAVE_ONLY,NRESOL)
  ELSE
   CALL EGPNORM_TRANS(ZGP,IFIELDS,NPROMA,ZAVE,ZMIN,ZMAX,LLAVE_ONLY,NRESOL)
  ENDIF
  IF (MYPROC == 1.OR.NPROC == 1) THEN
    IFLD=0
    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LGP)THEN
        IF(LLNOT_TL_AD&
                                & .OR.(YCOMP(JGFL)%CNAME==' LIQUID WATER'.OR.&
                                & YCOMP(JGFL)%CNAME==' ICE WATER'.OR.&
                                & YCOMP(JGFL)%CNAME==' CLOUD FRACTION'.OR.&
                                &      YCOMP(JGFL)%CNAME==' HUMIDITY')) THEN

          IFLD=IFLD+1
          WRITE(NULOUT,'(5(1X,A))') 'GPNORM',YCOMP(JGFL)%CNAME,'    AVERAGE         ',&
           & '     MINIMUM        ','      MAXIMUM        '
          WRITE(NULOUT,'(9X,A5,3(1X,E21.15))') 'AVE  ',ZAVE(IFLD)/NFLEVG,ZMIN(IFLD),ZMAX(IFLD)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  DEALLOCATE(ZGP)
  DEALLOCATE(ZAVE)
  DEALLOCATE(ZMIN)
  DEALLOCATE(ZMAX)

ENDIF
!     -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPNORM_GFL',1,ZHOOK_HANDLE)
END SUBROUTINE GPNORM_GFL
