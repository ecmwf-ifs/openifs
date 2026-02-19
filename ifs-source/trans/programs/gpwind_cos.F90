! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

PROGRAM GPWIND_COS
! Nils Wedi 2010, ECMWF, test transform of vor-div to u/v on lat-lon

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE GRIB_API
USE MPL_MODULE

IMPLICIT NONE

CHARACTER(LEN=13) :: OPTIONS
DATA OPTIONS/'g:G:l:t:sLOh;'/
CHARACTER(LEN=127) :: CLARG , OPTLET
CHARACTER(LEN=127) :: COUTGPF,  CINGPF, COUTSPF, CINSPEC,  CINTEMP, cltypeOfGrid
INTEGER(KIND=JPIM) ::  IFILE, IFILESP, IFILET, IGRIB_SP,IGRIB_H, IGRIB_OUT, IGRIB_T

INTEGER(KIND=JPIM) :: ISTACK
INTEGER(KIND=JPIM) :: getstackusage
INTEGER(KIND=JPIM) :: OPTVAL

INTEGER(KIND=JPIM) :: NERR,NSMAX,NDGLI, NDGLO, NSMAX_I, IDGL, NSIZE, ICOUNT, NDLON
INTEGER(KIND=JPIM) :: NOUT,NSPEC2,NGPTOT,NGPTOTG,IMAXFLD,IFLD, NSPEC2G, ISTRUNC
INTEGER(KIND=JPIM) :: ITAG,IRET, IOUT, ILAT, ILON
INTEGER(KIND=JPIM) :: IOUTGPF, I
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLOEN(:),ITO(:),NPRCIDS(:)
INTEGER(KIND=JPIM) :: GETOPT,JJ
REAL(KIND=JPRB),ALLOCATABLE :: ZSPEC(:,:),ZFPDAT(:),ZSPECG(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZG(:,:,:),&
 & ZGG(:,:),ZGGOUT(:)
LOGICAL :: LLSCAL, LLATLON_OUT,LLATLON_IN, LLSPECIN

INTEGER(KIND=JPIM) :: IPOI, JL, JGL, J, IR
REAL(KIND=JPRB)    :: ZDEG2RAD, ZLON, ZLONI
REAL(KIND=JPRB)    :: RPI, RA

INTEGER(KIND=JPIM) :: IPROC, NPROC, MYPROC, IMYSETV
INTEGER(KIND=JPIM) :: ILONS, IPLPRESENT, ITEST

!     ------------------------------------------------------------------

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "dist_grid.h"
#include "dist_spec.h"
#include "gath_grid.h"
#include "gath_spec.h"
#include "trans_inq.h"

! Initializations
NERR = 0
NOUT = 6
IMAXFLD = 2

! Set defaults for options

CINSPEC = 'insp.grib'
CINGPF  = 'in.grib'
CINTEMP = 'template.grib'
COUTGPF = 'out.grib'
COUTSPF = 'outspec.grib'
IPROC   = 2
LLSCAL   = .false.
NDGLO    = 0
NSMAX = 0
LLATLON_OUT=.false.
LLATLON_IN=.false.
LLSPECIN=.false.

RPI = 2.0_JPRB*ASIN(1.0_JPRB)
RA  = 6371229._JPRB

! Crack options

DO
  OPTVAL = GETOPT(OPTIONS,CLARG)
  IF(OPTVAL <= 0) EXIT
  OPTLET=CHAR(OPTVAL)
  IF(OPTVAL <= 0) EXIT
  IF(OPTLET /= 'g'.AND.&
   &OPTLET  /= 'l'.AND.&
   &OPTLET  /= 't'.AND.&
   &OPTLET  /= 's'.AND.&
   &OPTLET  /= 'h'.AND.&
   &OPTLET  /= 'L'.AND.&
   &OPTLET  /= 'O'.AND.&
   &OPTLET  /= 'G')  THEN
    CALL USE_GPT
    CALL ABOR1('GPSCAL_DER:ERROR IN ARGUMENT')
  ENDIF
  IF(OPTLET == 'h') THEN
    CALL USE_GPT
    STOP
  ELSEIF(OPTLET == 'l')  THEN
    READ(CLARG,*) NDGLO
  ELSEIF(OPTLET == 't')  THEN
    READ(CLARG,*) NSMAX
  ELSEIF(OPTLET == 'g')  THEN
    CINGPF = CLARG
  ELSEIF(OPTLET == 'G') THEN
    COUTGPF = CLARG
  ELSEIF(OPTLET == 's') THEN
    LLSPECIN=.true.
  ELSEIF(OPTLET == 'L') THEN 
    LLATLON_IN=.true.
  ELSEIF(OPTLET == 'O') THEN
    LLATLON_OUT=.true.
  ENDIF
ENDDO

IF( NDGLO == 0 ) THEN
  WRITE(NERR,*) 'NDGLO =',NDGLO
  CALL ABOR1('GPT: NUMBER OF GAUSSIAN LATITUDES NOT SET (use -l option)')
ENDIF
IF( NSMAX == 0 ) THEN
  WRITE(NERR,*) 'NSMAX =',NSMAX
  CALL ABOR1('GPT: TRUNCATION NOT SET (use -t option)')
ENDIF

! Message passing setup
! Participating processors limited by -P option
IF(IPROC > 1 ) THEN
  CALL MPL_INIT
  NPROC   = MPL_NPROC()
  MYPROC = MPL_MYRANK()
ELSE
  NPROC  = 1
  MYPROC = 1
  LLSCAL=.true.
ENDIF

ALLOCATE(NPRCIDS(NPROC))
DO JJ=1,NPROC
  NPRCIDS(JJ) = JJ
ENDDO

!-------------------------
ITAG = 191919

!------------------------------------------------------
! INPUT GRIDPOINT FIELD FOR LAT-LON OR GAUSSIAN INPUT
!------------------------------------------------------
! open file
CALL GRIB_OPEN_FILE(IFILE,CINGPF,'R')
! get handle
CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB_H, IRET)
! get key and check
CALL GRIB_GET(IGRIB_H,'typeOfGrid', cltypeOfGrid)
!!IF ( TRIM(cltypeOfGrid) /='reduced_gg' ) THEN
!  CALL ABOR1('GPT: Need gridpoint field in output resolution as input!')
!ENDIF
CALL GRIB_GET(IGRIB_H,'Nj',NDGLI)
CALL GRIB_GET_SIZE(IGRIB_H,'values',NSIZE)
!CALL GRIB_GET(KGRIB_HANDLE,'latitudeOfFirstGridPoint',INP)
!CALL GRIB_GET(KGRIB_HANDLE,'latitudeOfLastGridPoint',ISP)
!ZDLAT = ZDTOR*REAL(ISP-INP,JPRB)/REAL(ILATS-1,JPRB)
ALLOCATE(NLOEN(NDGLI))
CALL GRIB_GET(IGRIB_H,'PLPresent',IPLPresent)
IF (IPLPresent == 0) THEN
  CALL GRIB_GET(IGRIB_H,'numberOfPointsAlongAParallel',ILONS)
  NLOEN(:) = ILONS
ELSE
  CALL GRIB_GET(IGRIB_H,'pl',NLOEN)
ENDIF
NDLON=ILONS

!------------------------------------------------------
! TEMPLATE GRIDPOINT FIELD JUST TO SETUP OUTPUT (could link to in.grib)
!------------------------------------------------------
CALL GRIB_OPEN_FILE(IFILET,CINTEMP,'R')
CALL GRIB_NEW_FROM_FILE(IFILET,IGRIB_T, IRET)
CALL GRIB_CLOSE_FILE(IFILET)  

!------------------------------------------------------
! INPUT SPECTRAL FIELD IS OPTIONAL AND WILL THEN JUST PERFORM INV_TRANS
!------------------------------------------------------

IF( LLSPECIN ) THEN
  ! open file
  CALL GRIB_OPEN_FILE(IFILESP,CINSPEC,'R')
  ! get handle
  CALL GRIB_NEW_FROM_FILE(IFILESP,IGRIB_SP,IRET)
  ! get key and check
  CALL GRIB_GET(IGRIB_SP,'typeOfGrid', cltypeOfGrid)
  IF ( TRIM(cltypeOfGrid) /= 'sh' ) THEN
    CALL ABOR1('GPT: Need spectral field on input!')
  ENDIF
  CALL GRIB_GET(IGRIB_SP,'J',NSMAX_I)
  IF( NSMAX_I /= NSMAX ) THEN
    CALL ABOR1('GPT: NSMAX in file different from specified!')
  ENDIF
ENDIF

!------------------------------------------------------

! Prepare for transforms
CALL SETUP_TRANS0(KOUT=NOUT,KERR=0,KPRINTLEV=0,KMAX_RESOL=1,&
 &                KPRGPNS=NPROC,KPRGPEW=1,KPRTRW=NPROC,LDMPOFF=LLSCAL)

! setup regular Gaussian grid with lat/lon equivalent dual
IF( NDGLI > NDGLO+1 ) THEN
  CALL ABOR1('GPT: NDGLI > NDGLO+1 not possible !')
ENDIF

IF( NDGLI == NDGLO+1 ) THEN
  IDGL=NDGLO+2
  CALL SETUP_TRANS(KSMAX=NSMAX,KDGL=NDGLO,KDLON=NDLON,LDSPLIT=.FALSE.,&
   & LDLL=.TRUE.)
ELSE
  IDGL=NDGLO
  CALL SETUP_TRANS(KSMAX=NSMAX,KDGL=NDGLO,KDLON=NDLON,LDSPLIT=.FALSE.,&
   & LDLL=.TRUE.,LDSHIFTLL=.TRUE.)
ENDIF

CALL TRANS_INQ(KSPEC2=NSPEC2,KSPEC2G=NSPEC2G,KGPTOT=NGPTOT,KGPTOTG=NGPTOTG,KMYSETV=IMYSETV)

write(NERR,*) 'LATITUDES INPUT, LATITUDES OUPUT, RESOL', NDGLI, NDGLO, NSMAX
write(NERR,*) 'DIMS ', NSPEC2,NSPEC2G,NGPTOT,NGPTOTG

ALLOCATE(ZSPEC(IMAXFLD,NSPEC2))
ALLOCATE(ZG(NGPTOT,IMAXFLD,1))

IF( MYPROC == 1 ) THEN
  ALLOCATE(ZSPECG(2,NSPEC2G))
  ALLOCATE(ZGG(NGPTOTG,IMAXFLD))
ELSE
  ALLOCATE(ZSPECG(2,NSPEC2))
  ALLOCATE(ZGG(NGPTOT,IMAXFLD))
ENDIF

ALLOCATE(ITO(IMAXFLD))

! Gridpoint to spectral transform

ITO(:) = 1
IFLD = 0

ITEST = 0
IFLD=2

IF( .NOT.LLSPECIN ) THEN
  IF( MYPROC == 1 ) THEN
    IF( ITEST == 1 ) THEN
      ALLOCATE(ZFPDAT(NGPTOTG))
      ! test function COS(lambda)
      IPOI=0
      ZDEG2RAD=RPI/180._JPRB
      DO JGL=1,NDGLI
        ZLONI=360._JPRB/FLOAT(NLOEN(JGL))*ZDEG2RAD
        DO JL=1,NLOEN(JGL)
          IPOI=IPOI+1
          ZLON=FLOAT(JL-1)*ZLONI
          ZFPDAT(IPOI)=COS(ZLON)
        ENDDO
      ENDDO
      ZGG(:,1) = ZFPDAT(:)
      DEALLOCATE(ZFPDAT) 
    ELSE

      !needs mapping to field suitable for transforms
      ! U wind
      IF( NDGLI == NDGLO+1 ) THEN

        write(NERR,*) 'MAPPING INPUT FIELD TO EVEN BY DUPLICATING EQUATOR'

        ! odd number of latitudes, lat/lon field
        ALLOCATE(ZGGOUT(NSIZE))
        ! read input data
        CALL GRIB_GET(IGRIB_H,'values',ZGGOUT)
        ICOUNT=0
        IR=0
        DO ILAT=1,IDGL
          DO ILON=1,NDLON
            ! need equator 2x in input data, add one
            IF( ILAT==IDGL/2+1 .AND. IR==0 ) THEN
              ICOUNT=ICOUNT-NDLON
              IR=1
            ENDIF
            ICOUNT=ICOUNT+1
            ZGG(ILON+(ILAT-1)*NDLON,1) =  ZGGOUT(ICOUNT)
          ENDDO
        ENDDO
        CALL GRIB_RELEASE(IGRIB_H)
        CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB_H,IRET)
        ! V wind
        ! read input data
        CALL GRIB_GET(IGRIB_H,'values',ZGGOUT)
        ICOUNT=0
        IR=0
        DO ILAT=1,IDGL
          DO ILON=1,NDLON
            ! need equator 2x in input data, add one
            IF( ILAT==IDGL/2+1 .AND. IR==0 ) THEN
              ICOUNT=ICOUNT-NDLON
              IR=1
            ENDIF
            ICOUNT=ICOUNT+1
            ZGG(ILON+(ILAT-1)*NDLON,2) =  ZGGOUT(ICOUNT)
          ENDDO
        ENDDO
        DEALLOCATE(ZGGOUT)

      ELSE
        !  standard dimensions
        CALL GRIB_GET(IGRIB_H,'values',ZGG(1:NGPTOTG,1))
        CALL GRIB_RELEASE(IGRIB_H)
        CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB_H,IRET)
        CALL GRIB_GET(IGRIB_H,'values',ZGG(1:NGPTOTG,2))
      ENDIF
    ENDIF
  ENDIF

  ! Distribute gridpoint fields to processors
  CALL DIST_GRID(PGPG=ZGG,KFDISTG=IFLD,KFROM=ITO,PGP=ZG)
  
  ! Spectral transform
  CALL DIR_TRANS(PSPVOR=ZSPEC(1:1,:),PSPDIV=ZSPEC(2:2,:),PGP=ZG,LDLATLON=LLATLON_IN)
  
  ! Gather spectral fields to processor 1
  CALL GATH_SPEC(PSPECG=ZSPECG,KFGATHG=1,KTO=ITO,PSPEC=ZSPEC)

  IF(MYPROC == 1) THEN
    !  -----------------------------
    
    !*    WRITE THE SPECTRAL FIELD
    
    ALLOCATE(ZFPDAT(NSPEC2G))
    CALL GRIB_NEW_FROM_SAMPLES(IGRIB_OUT,'sh_sfc_grib2')
    CALL GRIB_OPEN_FILE(IOUT,COUTSPF,'w')
    CALL GRIB_SET(IGRIB_OUT,'gridType','sh')
    CALL GRIB_SET(IGRIB_OUT,'numberOfBitsContainingEachPackedValue',16)
    CALL GRIB_SET(IGRIB_OUT,'pentagonalResolutionParameterJ',NSMAX)
    CALL GRIB_SET(IGRIB_OUT,'pentagonalResolutionParameterK',NSMAX)
    CALL GRIB_SET(IGRIB_OUT,'pentagonalResolutionParameterM',NSMAX)
    CALL GRIB_SET(IGRIB_OUT,'laplacianOperator',0.5_JPRB)
    IF(NSMAX>= 213) THEN
      ISTRUNC = 20
    ELSE
      ISTRUNC = MIN(10,NSMAX)
    ENDIF
    CALL GRIB_SET(IGRIB_OUT,'subSetJ',ISTRUNC)
    CALL GRIB_SET(IGRIB_OUT,'subSetK',ISTRUNC)
    CALL GRIB_SET(IGRIB_OUT,'subSetM',ISTRUNC)
    CALL GRIB_SET(IGRIB_OUT,'paramId',80)
    CALL GRIB_SET(IGRIB_OUT,'level',0)
    CALL GRIB_SET(IGRIB_OUT,'values',ZSPECG(1,:))
    CALL GRIB_WRITE(IGRIB_OUT,IOUT)
    CALL GRIB_RELEASE(IGRIB_OUT)
    CALL GRIB_CLOSE_FILE(IOUT)
    DEALLOCATE(ZFPDAT)
    
  ENDIF
  
  ! Syncronize processors
  IF(NPROC > 1) THEN
    CALL MPL_BARRIER(CDSTRING='GPWIND_COS:')
  ENDIF

ELSE
  
  IF(MYPROC == 1) THEN
    ! read spectral input data vor/div
    ! VOR
    CALL GRIB_GET(IGRIB_SP,'values',ZSPECG(1,:))
    CALL GRIB_RELEASE(IGRIB_SP)
    ! DIV
    CALL GRIB_NEW_FROM_FILE(IFILESP,IGRIB_SP,IRET)
    CALL GRIB_GET(IGRIB_SP,'values',ZSPECG(2,:))
  ENDIF

  CALL DIST_SPEC(PSPECG=ZSPECG,KFDISTG=IFLD,KFROM=ITO,PSPEC=ZSPEC)

ENDIF

! inverse transform
CALL INV_TRANS(PSPVOR=ZSPEC(1:1,:),PSPDIV=ZSPEC(2:2,:),PGP=ZG,LDLATLON=LLATLON_OUT)

! Gather gridpoint fields to processor 1
CALL GATH_GRID(PGPG=ZGG,KFGATHG=IFLD,KTO=ITO,PGP=ZG)
IF(MYPROC == 1) THEN

  !*    WRITE THE GRIDPOINT FIELD

  ! open file
  CALL GRIB_OPEN_FILE(IOUTGPF,COUTGPF,'W')
  print *,'IFLD ',IFLD
  CALL GRIB_CLONE(IGRIB_T,IGRIB_OUT)

  ! U wind
  CALL GRIB_SET(IGRIB_OUT,'paramId',131)
  IF( NDGLI == NDGLO+1 ) THEN

    write(NERR,*) 'MAPPING OUTPUT FIELD BY REMOVING EQUATOR'

    ! odd number of latitudes
    ALLOCATE(ZGGOUT(NSIZE))
    ICOUNT=0
    DO ILAT=1,IDGL
      DO ILON=1,NDLON
        ! equator is 2x in output data, remove one or both
        IF( ILAT/=IDGL/2+1 .AND. .NOT.(ILAT==IDGL/2.AND..NOT.LLATLON_OUT)) THEN
          ICOUNT=ICOUNT+1
          ZGGOUT(ICOUNT) = ZGG(ILON+(ILAT-1)*NDLON,1)
        ENDIF
      ENDDO
    ENDDO
    write(NERR,*) 'OUTPUT DIMENSION: ', ICOUNT
    CALL GRIB_SET(IGRIB_OUT,'values',ZGGOUT(1: ICOUNT))
    DEALLOCATE(ZGGOUT) 
  ELSE
    !  standard dimensions
    CALL GRIB_SET(IGRIB_OUT,'values',ZGG(1:NGPTOTG,1))
  ENDIF
  CALL GRIB_WRITE(IGRIB_OUT,IOUTGPF)
  CALL GRIB_RELEASE(IGRIB_OUT)
  
  ! V wind
  CALL GRIB_CLONE(IGRIB_T,IGRIB_OUT)
  CALL GRIB_SET(IGRIB_OUT,'paramId',132)
  IF( NDGLI == NDGLO+1 ) THEN
    ! odd number of latitudes
    ALLOCATE(ZGGOUT(NSIZE))
    ICOUNT=0
    DO ILAT=1,IDGL
      DO ILON=1,NDLON
        ! equator is 2x in output data, remove one
        IF( ILAT/=IDGL/2+1 .AND. .NOT.(ILAT==IDGL/2.AND..NOT.LLATLON_OUT)) THEN
          ICOUNT=ICOUNT+1
          ZGGOUT(ICOUNT) = ZGG(ILON+(ILAT-1)*NDLON,2)
        ENDIF
      ENDDO
    ENDDO
    write(NERR,*) 'OUTPUT DIMENSION: ', ICOUNT
    CALL GRIB_SET(IGRIB_OUT,'values',ZGGOUT(1: ICOUNT))
    DEALLOCATE(ZGGOUT) 
  ELSE
    !  standard dimensions
    CALL GRIB_SET(IGRIB_OUT,'values',ZGG(1:NGPTOTG,2))
  ENDIF
  CALL GRIB_WRITE(IGRIB_OUT,IOUTGPF)
  CALL GRIB_RELEASE(IGRIB_OUT)

  ! close output file
  CALL GRIB_CLOSE_FILE(IOUTGPF)  
ENDIF

! Syncronize processors
IF(NPROC > 1) THEN
  CALL MPL_BARRIER(CDSTRING='GPWIND_COS:')
ENDIF
IFLD = 0
  
!           gather stack usage statistics
ISTACK = GETSTACKUSAGE()
IF(MYPROC == 1) THEN
  PRINT 9000, istack
9000 FORMAT("Stack Utilisation Information",/,&
      &"=============================",//,&
      &"Node           Size(Bytes)",/,&
      &"====           ===========",//,&
      &"   1",11x,I10)
  
  DO I=2,NPROC
    CALL MPL_RECV(ISTACK,KSOURCE=NPRCIDS(I),KTAG=I, &
     & CDSTRING='GPWIND_COS:')
    PRINT '(I4,11X,I10)', I,ISTACK
  ENDDO
ELSE
  CALL MPL_SEND(ISTACK,KDEST=NPRCIDS(1),KTAG=MYPROC, &
   &   CDSTRING='GPWIND_COS:')
ENDIF

!Close down message passing
!--------------------------
!write(0,*) "Calling barrier"
CALL MPL_BARRIER(CDSTRING='GPWIND_COS:')
CALL MPL_END()
!--------------------------
STOP

CONTAINS  

!     ------------------------------------------------------------------

SUBROUTINE USE_GPT

!   WRITE MESSAGE TO INFORM ABOUT CORRECT USE OF PROGRAM
WRITE(NERR,*) ' CORRECT OPTIONS AND ARGUMENTS FOR GPSCALAR_DER ARE :'
WRITE(NERR,*) '   -g    input filename (grid-point) (def: fort.11)'
WRITE(NERR,*) '   -G    output grid-point filename (def: fort.10)'
WRITE(NERR,*) '   -l    ndgl'
WRITE(NERR,*) '   -h    displays options and arguments'
END SUBROUTINE USE_GPT

!     ------------------------------------------------------------------

END PROGRAM GPWIND_COS

