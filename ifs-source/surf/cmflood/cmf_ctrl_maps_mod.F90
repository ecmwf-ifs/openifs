MODULE CMF_CTRL_MAPS_MOD
!==========================================================
!* PURPOSE: Control CaMa-Flood map/topography data
!
!* CONTAINS:
! -- CMF_MAPS_NMLIST   : configuration from namelist
! -- CMF_RIVMAP_INIT  : read & set river network map 
! -- CMF_TOPO_INIT    : read & set topography
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
! shared variables in module
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM
IMPLICIT NONE
SAVE
!*** NAMELIST/NMAP/ from inputnam
CHARACTER(LEN=256)              :: CNEXTXY         !! river network nextxy
CHARACTER(LEN=256)              :: CGRAREA         !! catchment area
CHARACTER(LEN=256)              :: CELEVTN         !! bank top elevation
CHARACTER(LEN=256)              :: CNXTDST         !! distance to next outlet
CHARACTER(LEN=256)              :: CRIVLEN         !! river channel length
CHARACTER(LEN=256)              :: CFLDHGT         !! floodplain elevation profile
!* river channel parameters
CHARACTER(LEN=256)              :: CRIVWTH         !! channel width
CHARACTER(LEN=256)              :: CRIVHGT         !! channel depth
CHARACTER(LEN=256)              :: CRIVMAN         !! river manning coefficient
!* optional maps
CHARACTER(LEN=256)              :: CPTHOUT         !! bifurcation channel table
CHARACTER(LEN=256)              :: CGDWDLY         !! Groundwater Delay Parameter
CHARACTER(LEN=256)              :: CMEANSL         !! mean sea level
!* netCDF map
LOGICAL                         :: LMAPCDF         !! true for netCDF map input
CHARACTER(LEN=256)              :: CRIVCLINC       !! river map netcdf
CHARACTER(LEN=256)              :: CRIVPARNC       !! river parameter netcdf (WIDTH,HEIGHT, Manning, ground wateer delay)
CHARACTER(LEN=256)              :: CMEANSLNC       !! mean sea level netCDF

NAMELIST/NMAP/     CNEXTXY,  CGRAREA,  CELEVTN,  CNXTDST, CRIVLEN, CFLDHGT, &
                   CRIVWTH,  CRIVHGT,  CRIVMAN,  CPTHOUT, CGDWDLY, CMEANSL,          &
                   LMAPCDF,  CRIVCLINC,CRIVPARNC,CMEANSLNC


CONTAINS
!####################################################################
! -- CMF_MAP_NMLIST   : configuration from namelist
! -- CMF_RIVMAP_INIT  : read & set river network map 
! -- CMF_TOPO_INIT    : read & set topography
!
!
!####################################################################
SUBROUTINE CMF_MAPS_NMLIST
! reed setting from namelist
! -- Called from CMF_DRV_NMLIST
USE YOS_CMF_INPUT,      ONLY: CSETFILE,NSETFILE,LMEANSL,LGDWDLY
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================
!*** 1. open namelist
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

NSETFILE=INQUIRE_FID()
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::MAP_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 2. default value
CNEXTXY="./nextxy.bin"
CGRAREA="./ctmare.bin"
CELEVTN="./elevtn.bin"
CNXTDST="./nxtdst.bin"
CRIVLEN="./rivlen.bin"
CFLDHGT="./fldhgt.bin"

CRIVWTH="./rivwth.bin"
CRIVHGT="./rivhgt.bin"
CRIVMAN="./rivman.bin"

CPTHOUT="./bifprm.txt"
CGDWDLY="NONE"
CMEANSL="NONE"

LMAPCDF=.FALSE.
CRIVCLINC="NONE"
CRIVPARNC="NONE"
CMEANSLNC="NONE"

!*** 3. read namelist
REWIND(NSETFILE)
READ(NSETFILE,NML=NMAP)

WRITE(LOGNAM,*)     "=== NAMELIST, NMAP ==="
WRITE(LOGNAM,*)     "LMAPCDF:   ", LMAPCDF
IF( LMAPCDF )THEN
  WRITE(LOGNAM,*)   "CRIVCLINC: ", TRIM(CRIVCLINC)
  WRITE(LOGNAM,*)   "CRIVPARNC: ", TRIM(CRIVPARNC)
  IF( LMEANSL ) THEN
    WRITE(LOGNAM,*) "CMEANSLNC: ", TRIM(CMEANSLNC)
  ENDIF
ELSE
  WRITE(LOGNAM,*)   "CNEXTXY:   ", TRIM(CNEXTXY)
  WRITE(LOGNAM,*)   "CGRAREA:   ", TRIM(CGRAREA)
  WRITE(LOGNAM,*)   "CELEVTN:   ", TRIM(CELEVTN)
  WRITE(LOGNAM,*)   "CNXTDST:   ", TRIM(CNXTDST)
  WRITE(LOGNAM,*)   "CRIVLEN:   ", TRIM(CRIVLEN)
  WRITE(LOGNAM,*)   "CFLDHGT:   ", TRIM(CFLDHGT)

  WRITE(LOGNAM,*)   "CRIVWTH:   ", TRIM(CRIVWTH)
  WRITE(LOGNAM,*)   "CRIVHGT:   ", TRIM(CRIVHGT)
  WRITE(LOGNAM,*)   "CRIVMAN:   ", TRIM(CRIVMAN)

  WRITE(LOGNAM,*)   "CPTHOUT:   ", TRIM(CPTHOUT)
  IF( LGDWDLY )THEN
    WRITE(LOGNAM,*) "CGDWDLY:    ",TRIM(CGDWDLY)
  ENDIF
  IF( LMEANSL )THEN
    WRITE(LOGNAM,*) "CMEANSL:   ", TRIM(CMEANSL)
  ENDIF
ENDIF

CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::MAP_NMLIST: end"

END SUBROUTINE CMF_MAPS_NMLIST
!####################################################################





!####################################################################
SUBROUTINE CMF_RIVMAP_INIT
! read & set river network map 
! -- call from CMF_DRV_INIT
USE YOS_CMF_INPUT,      ONLY: TMPNAM, NX,NY,NLFP, LPTHOUT
USE YOS_CMF_MAP,        ONLY: I2NEXTX,I2NEXTY, I2RIVSEQ,I2REGION, REGIONALL,REGIONTHIS, &
                            & I1SEQX, I1SEQY,  I1NEXT,  I2VECTOR, D1LON,    D1LAT,      &
                            & NSEQRIV,  NSEQALL,  NSEQMAX,  RIVSEQMAX
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: river network initialization'

! *** 1. ALLOCATE ARRAYS
ALLOCATE( I2NEXTX(NX,NY) )
ALLOCATE( I2NEXTY(NX,NY) )
ALLOCATE( I2RIVSEQ(NX,NY) )
ALLOCATE( I2REGION(NX,NY) )
ALLOCATE( D1LON(NX) )
ALLOCATE( D1LAT(NY) )

!============================
!*** 2a. read river network map
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: read nextXY & set lat lon'
IF( LMAPCDF )THEN
#ifdef UseCDF_CMF
  CALL READ_MAP_CDF
#endif
ELSE
  CALL READ_MAP_BIN
ENDIF

!*** 2b. calculate river sequence & regions
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: calc rivseq & region'
CALL CALC_RIVSEQ
CALL CALC_REGION

!============================
!*** 3. conversion 2D map -> 1D vector
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: calculate 1d river sequence'

CALL CALC_1D_SEQ                                  !! 2D map to 1D vector conversion. for faster calculation

WRITE(LOGNAM,*) '  NSEQRIV=',NSEQRIV
WRITE(LOGNAM,*) '  NSEQALL=',NSEQALL

! for MPI, nor used in v395
!IF( REGIONTHIS==1 )THEN
!  WRITE(LOGNAM,*) 'INIT_MAP: NSEQMAX=',NSEQMAX
!  WRITE(LOGNAM,*) 'INIT_MAP: RIVSEQMAX=',RIVSEQMAX
!ENDIF

!*** 3c. Write Map Data                                       !! used for combining mpi distributed output into one map
IF( REGIONTHIS==1 )THEN
  TMPNAM=INQUIRE_FID()
  OPEN(TMPNAM,FILE='./mapdata.txt',FORM='FORMATTED')
  WRITE(TMPNAM,*) 'NX',        NX
  WRITE(TMPNAM,*) 'NY',        NY
  WRITE(TMPNAM,*) 'NLFP',      NLFP
  WRITE(TMPNAM,*) 'REGIONALL', REGIONALL
  WRITE(TMPNAM,*) 'NSEQMAX',   NSEQMAX
  CLOSE(TMPNAM)
ENDIF

!============================
!*** 4.  bifurcation channel parameters
IF( LPTHOUT )THEN
  WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: read bifurcation channel setting'
  CALL READ_BIFPARAM
ENDIF

WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: end'

CONTAINS
!==========================================================
!+ READ_MAP_BIN
!+ READ_MAP_CDF
!+ CALC_RIVSEQ
!+ CALC_REGION
!+ CALC_1D_SEQ
!+ READ_BIFPRM
!==========================================================
SUBROUTINE READ_MAP_BIN
USE YOS_CMF_INPUT,      ONLY: TMPNAM, LMAPEND
USE YOS_CMF_INPUT,      ONLY: WEST,EAST,NORTH,SOUTH
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID, CONV_ENDI
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)         :: IX,IY
!==========================================================
!*** read river map
WRITE(LOGNAM,*)'RIVMAP_INIT: nextxy binary: ',TRIM(CNEXTXY)
TMPNAM=INQUIRE_FID()
OPEN(TMPNAM,FILE=CNEXTXY,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) I2NEXTX
READ(TMPNAM,REC=2) I2NEXTY
CLOSE(TMPNAM)

IF ( LMAPEND )THEN
  CALL CONV_ENDI(I2NEXTX,NX,NY)
  CALL CONV_ENDI(I2NEXTY,NX,NY)
ENDIF

!*** calculate lat, lon
IF( WEST>-180._JPRB .and. EAST<360._JPRB .and. SOUTH>=-180._JPRB .and. NORTH<=180._JPRB )THEN
  DO IX=1,NX
    D1LON(IX)=WEST +(DBLE(IX)-0.5D0)*(EAST-WEST)  /DBLE(NX)
  ENDDO
  DO IY=1,NY
    D1LAT(IY)=NORTH-(DBLE(IY)-0.5D0)*(NORTH-SOUTH)/DBLE(NY)
  ENDDO
ENDIF

END SUBROUTINE READ_MAP_BIN
!==========================================================
!+
!+
!+
!==========================================================
#ifdef UseCDF_CMF
SUBROUTINE READ_MAP_CDF
USE CMF_UTILS_MOD  ,ONLY: NCERROR
USE NETCDF
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)              :: NCID,VARID
!================================================
WRITE(LOGNAM,*)'RIVMAP_INIT: nextxy netCDF: ', TRIM(CRIVCLINC)

CALL NCERROR (NF90_OPEN(CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVCLINC) )

!*** next xy
CALL NCERROR ( NF90_INQ_VARID(NCID,'nextx',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTX),'reading data' ) 

CALL NCERROR ( NF90_INQ_VARID(NCID,'nexty',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTY),'reading data' )

!*** lat, lon
CALL NCERROR ( NF90_INQ_VARID(NCID,'lat',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LAT),'reading data' )

CALL NCERROR ( NF90_INQ_VARID(NCID,'lon',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LON),'reading data' )

CALL NCERROR( NF90_CLOSE(NCID))

END SUBROUTINE READ_MAP_CDF
#endif
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CALC_RIVSEQ
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM) :: IX, IY, JX, JY
INTEGER(KIND=JPIM) :: SEQNOW, NEXT
!================================================
WRITE(LOGNAM,*) 'RIVMAP_INIT: calculate river sequence'
I2RIVSEQ(:,:)=1
DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY)>0 )THEN
      JX=I2NEXTX(IX,IY)
      JY=I2NEXTY(IX,IY)
      I2RIVSEQ(JX,JY)=0
    ELSEIF( I2NEXTX(IX,IY)==-9999 )THEN
      I2RIVSEQ(IX,IY)=-9999
    ENDIF
  END DO
END DO

SEQNOW=0
NEXT=1
DO WHILE(NEXT>0 )
  SEQNOW=SEQNOW+1
  NEXT=0
  DO IY=1, NY
    DO IX=1, NX
      IF( I2RIVSEQ(IX,IY)==SEQNOW )THEN
        IF( I2NEXTX(IX,IY)>0 )THEN
          JX=I2NEXTX(IX,IY)
          JY=I2NEXTY(IX,IY)
          IF( I2RIVSEQ(JX,JY)<=SEQNOW )THEN
            NEXT=NEXT+1
            I2RIVSEQ(JX,JY)=SEQNOW+1
          ENDIF
        ENDIF
      ENDIF
    END DO
  ENDDO
END DO
RIVSEQMAX=SEQNOW
WRITE(LOGNAM,*) '  RIVSEQMAX = ', RIVSEQMAX

END SUBROUTINE CALC_RIVSEQ
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CALC_REGION    !! evenly allocate pixels to mpi nodes (not used in v395)
USE YOS_CMF_INPUT,           ONLY: IMIS
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM),ALLOCATABLE  :: I2BASIN(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE  :: I2UPGRID(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE  :: BASINGRID(:)
INTEGER(KIND=JPIM),ALLOCATABLE  :: REGIONGRID(:)
INTEGER(KIND=JPIM),ALLOCATABLE  :: BASINREGION(:)
!
INTEGER(KIND=JPIM)              :: IX,IY,   JX,JY,   KX,KY, ISEQ, UPG
INTEGER(KIND=JPIM)              :: IBASIN,  BASINMAX
INTEGER(KIND=JPIM)              :: IREGION, JREGION, GRIDMIN
!================================================
WRITE(LOGNAM,*) 'RIVMAP_INIT: region code'

ALLOCATE(I2BASIN(NX,NY))
ALLOCATE(I2UPGRID(NX,NY))

WRITE(LOGNAM,*)'  calculate upstream grid number'
I2UPGRID(:,:)=0
DO IY=1, NY
  DO IX=1, NX
    IF( I2RIVSEQ(IX,IY)==1 ) THEN
      JX=IX
      JY=IY
      I2UPGRID(JX,JY)=1
      UPG=I2UPGRID(JX,JY)
      DO WHILE( I2NEXTX(JX,JY)>0 )  !! IF RIVER REACHES MOUTH, END LOOP
        KX=I2NEXTX(JX,JY)
        KY=I2NEXTY(JX,JY)
        JX=KX
        JY=KY
        IF( I2UPGRID(JX,JY)==0 )THEN               !! GRIDS FIRSTLY CHECKED
          I2UPGRID(JX,JY)=UPG+1
          UPG=I2UPGRID(JX,JY)
        ELSE                                       !! GRIDS ALREADY CHECKED
          I2UPGRID(JX,JY)=I2UPGRID(JX,JY)+UPG
        ENDIF
      END DO
    ENDIF
  END DO
END DO

WRITE(LOGNAM,*)'  calculate basin'
I2BASIN(:,:)=0
IBASIN=0
DO ISEQ=RIVSEQMAX, 1, -1
  DO IY=1, NY
    DO IX=1, NX
      IF( I2RIVSEQ(IX,IY)==ISEQ .AND. I2NEXTX(IX,IY)<0 .AND. I2RIVSEQ(IX,IY)/=IMIS )THEN
        IBASIN=IBASIN+1
        I2BASIN(IX,IY)=IBASIN
      ENDIF
    END DO
  END DO
END DO
BASINMAX=IBASIN

DO IY=1, NY
  DO IX=1, NX
    JX=IX
    JY=IY
    DO WHILE( I2BASIN(JX,JY)==0 .and. I2NEXTX(JX,JY)>0 )
      KX=I2NEXTX(JX,JY)
      KY=I2NEXTY(JX,JY)
      JX=KX
      JY=KY
    END DO
    IBASIN=I2BASIN(JX,JY)
    JX=IX
    JY=IY
    DO WHILE( I2BASIN(JX,JY)==0 .and. I2NEXTX(JX,JY)>0  )
      I2BASIN(JX,JY)=IBASIN
      KX=I2NEXTX(JX,JY)
      KY=I2NEXTY(JX,JY)
      JX=KX
      JY=KY
    END DO
  END DO
END DO
ALLOCATE(BASINGRID(BASINMAX))
ALLOCATE(BASINREGION(BASINMAX))
ALLOCATE(REGIONGRID(REGIONALL))

WRITE(LOGNAM,*)'  allocate basin to cpu (MPI)'
DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY)<0 .and. I2NEXTX(IX,IY)/=IMIS )THEN
      IBASIN=I2BASIN(IX,IY)
      BASINGRID(IBASIN)=I2UPGRID(IX,IY)
    ENDIF
  END DO
END DO


REGIONGRID(:)=0
JREGION=1
DO IBASIN=1, BASINMAX
  GRIDMIN=NX*NY
  DO IREGION=1, REGIONALL
    IF( REGIONGRID(IREGION) < GRIDMIN )THEN
      GRIDMIN=REGIONGRID(IREGION)
      JREGION=IREGION
    ENDIF
  END DO
  BASINREGION(IBASIN)=JREGION
  REGIONGRID(JREGION)=REGIONGRID(JREGION)+BASINGRID(IBASIN)
!   write(*,*)IBASIN,BASINGRID(IBASIN),JREGION
END DO

NSEQMAX=0
DO IREGION=1, REGIONALL
  IF( REGIONTHIS==1 )THEN
    WRITE(LOGNAM,*) 'CALC_REGION: ', IREGION, REGIONGRID(IREGION)
  ENDIF
!   WRITE(*,*) 'CALC_REGION: ', IREGION, REGIONGRID(IREGION)
  NSEQMAX=MAX(NSEQMAX,REGIONGRID(IREGION))
END DO

DO IY=1, NY
  DO IX=1, NX
    IF( I2BASIN(IX,IY)>0 .and. I2NEXTX(IX,IY)/=-9999 )THEN
      IBASIN=I2BASIN(IX,IY)
      I2REGION(IX,IY)=BASINREGION(IBASIN)
    ENDIF
  END DO
END DO

DEALLOCATE(BASINGRID)
DEALLOCATE(BASINREGION)
DEALLOCATE(REGIONGRID)
DEALLOCATE(I2BASIN)
DEALLOCATE(I2UPGRID)

END SUBROUTINE CALC_REGION
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CALC_1D_SEQ
USE YOS_CMF_INPUT,           ONLY: IMIS
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)              :: IX,IY,JX,JY,ISEQ
!================================================
WRITE(LOGNAM,*) 'RIVMAP_INIT: convert 2D map to 1D sequence'

ALLOCATE( I1SEQX(NSEQMAX) )
ALLOCATE( I1SEQY(NSEQMAX) )
ALLOCATE( I1NEXT(NSEQMAX) )
ALLOCATE( I2VECTOR(NX,NY) )
I1SEQX(:)=0
I1SEQY(:)=0
I1NEXT(:)=0
I2VECTOR(:,:)=0

ISEQ=0
DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY).GT.0 .and. I2REGION(IX,IY)==REGIONTHIS )THEN
      ISEQ=ISEQ+1
      I1SEQX(ISEQ)=IX
      I1SEQY(ISEQ)=IY
      I2VECTOR(IX,IY)=ISEQ
    ENDIF
  END DO
END DO
NSEQRIV=ISEQ

DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY).LT.0 .AND. I2NEXTX(IX,IY).NE.IMIS .AND. I2REGION(IX,IY)==REGIONTHIS )THEN
      ISEQ=ISEQ+1
      I1SEQX(ISEQ)=IX
      I1SEQY(ISEQ)=IY
      I2VECTOR(IX,IY)=ISEQ
    ENDIF
  END DO
END DO
NSEQALL=ISEQ

DO ISEQ=1, NSEQALL
  IX=I1SEQX(ISEQ)
  IY=I1SEQY(ISEQ)
  IF( I2NEXTX(IX,IY)>0 )THEN
    JX=I2NEXTX(IX,IY)
    JY=I2NEXTY(IX,IY)
    I1NEXT(ISEQ)=I2VECTOR(JX,JY)
  ELSE
    I1NEXT(ISEQ)=-9
  ENDIF
END DO
      
END SUBROUTINE CALC_1D_SEQ
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE READ_BIFPARAM    !! evenly allocate pixels to mpi nodes (not used in v395)
USE YOS_CMF_INPUT,      ONLY: PMANRIV, PMANFLD
USE YOS_CMF_MAP,        ONLY: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN,&
                            & PTH_DST, PTH_ELV, PTH_WTH,  PTH_MAN
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)         :: IX,IY, JX,JY
INTEGER(KIND=JPIM)         :: IPTH,  ILEV,  NPTHOUT1
REAL(KIND=JPRB)            :: PELV,  PWTH,  PDPH
!================================================
WRITE(LOGNAM,*)"RIVMAP_INIT: Bifuraction channel:", TRIM(CPTHOUT)

TMPNAM=INQUIRE_FID()
OPEN(TMPNAM,FILE=CPTHOUT,FORM='FORMATTED')
READ(TMPNAM,*) NPTHOUT,NPTHLEV

WRITE(LOGNAM,*) "Bifurcation channel dimantion", NPTHOUT, NPTHLEV

ALLOCATE( PTH_UPST(NPTHOUT) )
ALLOCATE( PTH_DOWN(NPTHOUT) )
ALLOCATE( PTH_DST(NPTHOUT)  )
ALLOCATE( PTH_ELV(NPTHOUT,NPTHLEV) )
ALLOCATE( PTH_WTH(NPTHOUT,NPTHLEV) )
ALLOCATE( PTH_MAN(NPTHLEV)  )

NPTHOUT1=0
DO IPTH=1, NPTHOUT
  READ(TMPNAM,*) IX, IY, JX, JY, PTH_DST(IPTH), PELV, PDPH, (PTH_WTH(IPTH,ILEV),ILEV=1,NPTHLEV)
  PTH_UPST(IPTH)=I2VECTOR(IX,IY)
  PTH_DOWN(IPTH)=I2VECTOR(JX,JY)
  IF (PTH_UPST(IPTH) > 0 .AND. PTH_DOWN(IPTH) > 0) THEN
    NPTHOUT1=NPTHOUT1+1
  ENDIF
  DO ILEV=1, NPTHLEV
    IF( ILEV==1 )THEN            !!ILEV=1: water channel bifurcation. consider bifurcation channel depth
      PWTH=PTH_WTH(IPTH,ILEV)
      IF( PWTH>0 )then
        PTH_ELV(IPTH,ILEV)=PELV - PDPH
      ELSE
        PTH_ELV(IPTH,ILEV)=1.D20
      ENDIF
    ELSE
      PWTH=PTH_WTH(IPTH,ILEV)
      IF( PWTH>0 )then
        PTH_ELV(IPTH,ILEV)=PELV + ILEV - 2.0    !! ILEV=2: bank top level 
      ELSE
        PTH_ELV(IPTH,ILEV)=1.D20
      ENDIF
    ENDIF
  END DO
END DO
CLOSE(TMPNAM)

DO ILEV=1, NPTHLEV
  IF( ILEV==1 )THEN
    PTH_MAN(ILEV)=PMANRIV
  ELSE
    PTH_MAN(ILEV)=PMANFLD
  ENDIF
END DO

IF (NPTHOUT /= NPTHOUT1) THEN
  WRITE(LOGNAM,*)"Bifuraction channel outside of domain. Only valid:", NPTHOUT1
ENDIF

END SUBROUTINE READ_BIFPARAM
!==========================================================

END SUBROUTINE CMF_RIVMAP_INIT
!####################################################################





!####################################################################
SUBROUTINE CMF_TOPO_INIT
! read & set topography map 
! -- call from CMF_DRV_INIT
USE YOS_CMF_INPUT,  ONLY: TMPNAM,   NX,NY,NLFP, LMAPEND,  &
                        & LFPLAIN,  LMEANSL,  LGDWDLY,  LSLPMIX, &
                        & LSLOPEMOUTH
USE YOS_CMF_MAP,    ONLY: D2NXTDST, D2GRAREA, D2ELEVTN, D2RIVLEN, &
                        & D2RIVWTH, D2RIVHGT, D2FLDHGT, D2RIVELV, &
                        & D2FLDGRD, D2RIVMAN, D2RIVSTOMAX, D2FLDSTOMAX,  &
                        & DFRCINC,  NSEQALL,  NSEQMAX, D2MEANSL, D2DWNELV, &
                        & D2GDWDLY, I2MASK, D2ELEVSLOPE
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

WRITE(LOGNAM,*) 'CMF::TOPO_INIT: topography map initialization'

! *** 1. ALLOCATE ARRAYS
ALLOCATE( D2GRAREA(NSEQMAX,1) )
ALLOCATE( D2ELEVTN(NSEQMAX,1) )
ALLOCATE( D2NXTDST(NSEQMAX,1) )
ALLOCATE( D2RIVLEN(NSEQMAX,1) )
ALLOCATE( D2RIVWTH(NSEQMAX,1) )
ALLOCATE( D2RIVHGT(NSEQMAX,1) )
ALLOCATE( D2FLDHGT(NSEQMAX,1,NLFP) )
ALLOCATE( D2RIVMAN(NSEQMAX,1) )
ALLOCATE( D2MEANSL(NSEQMAX,1) )
ALLOCATE( D2DWNELV(NSEQMAX,1) )
ALLOCATE( D2GDWDLY(NSEQMAX,1) )
ALLOCATE( I2MASK(NSEQMAX,1) )

D2GRAREA(:,:)  =0._JPRB
D2ELEVTN(:,:)  =0._JPRB
D2NXTDST(:,:)  =0._JPRB
D2RIVLEN(:,:)  =0._JPRB
D2RIVWTH(:,:)  =0._JPRB
D2FLDHGT(:,:,:)=0._JPRB
D2RIVMAN(:,:)  =0._JPRB
D2MEANSL(:,:)  =0._JPRB
D2DWNELV(:,:)  =0._JPRB
D2GDWDLY(:,:)  =0._JPRB
I2MASK(:,:)    =0._JPIM     !! mask for calculation (only for IFS slopemix in v395)

!============================
! *** 2. Read topo map
WRITE(LOGNAM,*) 'CMF::TOPO_INIT: read topography maps'
IF ( .not. LMAPCDF ) THEN
  CALL READ_TOPO_BIN
ELSE
  CALL READ_TOPO_CDF
ENDIF

!============================
! *** 3a. Calc Channel Parameters
WRITE(LOGNAM,*) 'TOPO_INIT: calc river channel parameters'

ALLOCATE(D2RIVSTOMAX(NSEQMAX,1))
ALLOCATE(D2RIVELV(NSEQMAX,1))

IF ( LFPLAIN ) THEN
  D2RIVSTOMAX(:,:) = D2RIVLEN(:,:) * D2RIVWTH(:,:) * D2RIVHGT(:,:)
ELSE
  WRITE(LOGNAM,*) 'TOPO_INIT: no floodplain (rivstomax=1.D18)'
  D2RIVSTOMAX(:,:) = 1.D18
ENDIF
D2RIVELV(:,:) = D2ELEVTN(:,:) - D2RIVHGT(:,:)

!*** 3b. Calc Channel Parameters
WRITE(LOGNAM,*) 'TOPO_INIT: calc floodplain parameters'

ALLOCATE(D2FLDSTOMAX(NSEQMAX,1,NLFP))
ALLOCATE(D2FLDGRD(NSEQMAX,1,NLFP))
CALL SET_FLDSTG

!*** 3c. Calc downstream boundary
WRITE(LOGNAM,*) 'TOPO_INIT: calc downstream boundary elevation'
D2DWNELV(:,:)=D2ELEVTN(:,:)
IF( LMEANSL ) THEN
  D2DWNELV(:,:)=D2ELEVTN(:,:)+D2MEANSL(:,:)
ENDIF

CONTAINS
!==========================================================
!+ READ_TOPO_BIN
!+ READ_TOPO_CDF
!+ SET_FLDSTG
!+ SET_SLOPEMIX
!==========================================================
SUBROUTINE READ_TOPO_BIN
USE CMF_UTILS_MOD,       ONLY: MAP2VEC, CONV_END,  INQUIRE_FID
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)          :: ILFP
REAL(KIND=JPRM),ALLOCATABLE :: R2TEMP(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: D2TEMP(:,:)
!================================================
ALLOCATE(R2TEMP(NX,NY))
ALLOCATE(D2TEMP(NSEQMAX,1))

TMPNAM=INQUIRE_FID()

WRITE(LOGNAM,*)'TOPO_INIT: unit-catchment area : ',TRIM(CGRAREA) 
OPEN(TMPNAM,FILE=CGRAREA,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL MAP2VEC(R2TEMP,D2GRAREA)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: ground elevation : ',TRIM(CELEVTN)
OPEN(TMPNAM,FILE=CELEVTN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL MAP2VEC(R2TEMP,D2ELEVTN)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: downstream distance : ',TRIM(CNXTDST)
OPEN(TMPNAM,FILE=CNXTDST,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL MAP2VEC(R2TEMP,D2NXTDST)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: river channel length : ',TRIM(CRIVLEN)
OPEN(TMPNAM,FILE=CRIVLEN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL MAP2VEC(R2TEMP,D2RIVLEN)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: floodplain elevation profile : ',TRIM(CFLDHGT)
OPEN(TMPNAM,FILE=TRIM(CFLDHGT),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
DO ILFP=1,NLFP
  READ(TMPNAM,REC=ILFP) R2TEMP
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
  CALL MAP2VEC(R2TEMP,D2TEMP)
  D2FLDHGT(:,:,ILFP)= D2TEMP(:,:)
ENDDO
CLOSE(TMPNAM)

!*** river channel / groundwater parameters)

WRITE(LOGNAM,*)'TOPO_INIT: river channel depth : ',TRIM(CRIVHGT)
OPEN(TMPNAM,FILE=CRIVHGT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL MAP2VEC(R2TEMP,D2RIVHGT)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: river channel width : ',TRIM(CRIVWTH)
OPEN(TMPNAM,FILE=CRIVWTH,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL MAP2VEC(R2TEMP,D2RIVWTH)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: manning coefficient river: ',TRIM(CRIVMAN)
OPEN(TMPNAM,FILE=TRIM(CRIVMAN),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL MAP2VEC(R2TEMP,D2RIVMAN)
CLOSE(TMPNAM)

IF( LGDWDLY )THEN
  WRITE(LOGNAM,*)'TOPO_INIT: groundwater delay parameter: ',TRIM(CGDWDLY)
  OPEN(TMPNAM,FILE=TRIM(CGDWDLY),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
  READ(TMPNAM,REC=1) R2TEMP(:,:)
    IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
  CALL MAP2VEC(R2TEMP,D2GDWDLY)
  CLOSE(TMPNAM)
ENDIF

IF( LSLPMIX )THEN
  WRITE(LOGNAM,*)'TOPO_INIT: LSLPMIX only used in IFS, not availabke with binary map'
ENDIF

IF( LSLOPEMOUTH )THEN
  WRITE(LOGNAM,*)'TOPO_INIT: LSLOPEMOUTH only used in IFS, not availabke with binary map'
ENDIF


! ==========

IF( LMEANSL ) THEN
  WRITE(LOGNAM, *)'TOPO_INIT: mean sea level: ', TRIM(CMEANSL)
  OPEN(TMPNAM, FILE=CMEANSL, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*NX*NY)
  READ(TMPNAM, REC=1) R2TEMP(:,:)
    IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
  CALL MAP2VEC(R2TEMP, D2MEANSL)
  CLOSE(TMPNAM)
ENDIF

DEALLOCATE(R2TEMP)
DEALLOCATE(D2TEMP)

END SUBROUTINE READ_TOPO_BIN
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE READ_TOPO_CDF
#ifdef UseCDF_CMF
USE NETCDF 
USE CMF_UTILS_MOD,            ONLY: NCERROR,MAP2VEC
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)               :: NCID,VARID,STATUS
INTEGER(KIND=JPIM)               :: ILEV
REAL(KIND=JPRM),ALLOCATABLE      :: R2TEMP(:,:)
REAL(KIND=JPRB),ALLOCATABLE      :: D2TEMP(:,:)
!================================================
ALLOCATE(R2TEMP(NX,NY))
ALLOCATE(D2TEMP(NSEQMAX,1))

!! CLIM FILE
CALL NCERROR (NF90_OPEN(CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVCLINC) )

WRITE(LOGNAM,*)'TOPO_INIT: ctmare:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'ctmare',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' )
CALL MAP2VEC(R2TEMP,D2GRAREA)

WRITE(LOGNAM,*)'TOPO_INIT: elevtn:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'elevtn',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL MAP2VEC(R2TEMP,D2ELEVTN)

WRITE(LOGNAM,*)'TOPO_INIT: nxtdst:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'nxtdst',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL MAP2VEC(R2TEMP,D2NXTDST)

WRITE(LOGNAM,*)'TOPO_INIT: rivlen:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivlen',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL MAP2VEC(R2TEMP,D2RIVLEN)

WRITE(LOGNAM,*)'TOPO_INIT: fldhgt:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'fldhgt',VARID),'getting id' )
DO ILEV=1,NLFP
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP,(/1,1,ILEV/),(/NX,NY,1/)),'reading data' ) 
  CALL MAP2VEC(R2TEMP,D2TEMP)
  D2FLDHGT(:,:,ILEV)=D2TEMP(:,:)
ENDDO

IF ( LSLOPEMOUTH ) THEN
  ALLOCATE( D2ELEVSLOPE(NSEQMAX,1) )
  WRITE(LOGNAM,*)'TOPO_INIT: elevslope:',TRIM(CRIVPARNC)
  STATUS = NF90_INQ_VARID(NCID,'elevslope',VARID)
  IF (STATUS /= 0 ) THEN
    WRITE(LOGNAM,*)'TOPO_INIT: elevslope: not present, aborting'
    STOP 9 
  ELSE
    CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  ENDIF 
  CALL MAP2VEC(R2TEMP,D2ELEVSLOPE)
ENDIF

CALL NCERROR( NF90_CLOSE(NCID))

!!========== 
!! PAR FILE (river channel / groundwater parameters)
CALL NCERROR (NF90_OPEN(CRIVPARNC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVPARNC) )

WRITE(LOGNAM,*)'TOPO_INIT: rivwth:',TRIM(CRIVPARNC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivwth',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL MAP2VEC(R2TEMP,D2RIVWTH)

WRITE(LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CRIVPARNC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivhgt',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL MAP2VEC(R2TEMP,D2RIVHGT)

WRITE(LOGNAM,*)'TOPO_INIT: rivman:',TRIM(CRIVPARNC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivman',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL MAP2VEC(R2TEMP,D2RIVMAN)

IF ( LGDWDLY ) THEN
  WRITE(LOGNAM,*)'TOPO_INIT: GDWDLY:',TRIM(CRIVPARNC)
  STATUS = NF90_INQ_VARID(NCID,'gdwdly',VARID)
  IF (STATUS /= 0 ) THEN
    WRITE(LOGNAM,*)'TOPO_INIT: GDWDLY: not present, setting to zero'
    R2TEMP(:,:) = 0._JPRB
  ELSE
    CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  ENDIF 
  CALL MAP2VEC(R2TEMP,D2GDWDLY)
ENDIF

I2MASK(:,:)=0_JPIM
IF ( LSLPMIX ) THEN
  CALL SET_SLOPEMIX
ENDIF

CALL NCERROR( NF90_CLOSE(NCID))

!!========== 
!! MEAN SEA LEVEL FILE
IF( LMEANSL ) THEN
  CALL NCERROR (NF90_OPEN(CMEANSLNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMEANSLNC) )
  WRITE(LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CMEANSLNC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'meansl',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL MAP2VEC ( R2TEMP,D2MEANSL  )
  CALL NCERROR ( NF90_CLOSE(NCID) )
ENDIF 

DEALLOCATE(R2TEMP)
DEALLOCATE(D2TEMP)
#endif
END SUBROUTINE READ_TOPO_CDF
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE SET_FLDSTG
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)  ::  ISEQ, I
REAL(KIND=JPRB)     ::  DSTONOW
REAL(KIND=JPRB)     ::  DSTOPRE
REAL(KIND=JPRB)     ::  DHGTPRE
REAL(KIND=JPRB)     ::  DWTHINC
!================================================
D2FLDSTOMAX(:,:,:) = 0.D0
D2FLDGRD(:,:,:)    = 0.D0
DFRCINC=dble(NLFP)**(-1.)
!
DO ISEQ=1, NSEQALL
  DSTOPRE = D2RIVSTOMAX(ISEQ,1)
  DHGTPRE = 0.D0
  DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC
  DO I=1, NLFP
    DSTONOW = D2RIVLEN(ISEQ,1) * ( D2RIVWTH(ISEQ,1) + DWTHINC*(DBLE(I)-0.5) ) * (D2FLDHGT(ISEQ,1,I)-DHGTPRE)
    D2FLDSTOMAX(ISEQ,1,I) = DSTOPRE + DSTONOW
    D2FLDGRD(ISEQ,1,I) = (D2FLDHGT(ISEQ,1,I)-DHGTPRE) * DWTHINC**(-1.)
    DSTOPRE = D2FLDSTOMAX(ISEQ,1,I)
    DHGTPRE = D2FLDHGT(ISEQ,1,I)
  END DO
END DO
!
END SUBROUTINE SET_FLDSTG
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE SET_SLOPEMIX    !! only used in IFS0
#ifdef UseCDF_CMF
USE NETCDF 
USE CMF_UTILS_MOD,           ONLY: NCERROR,MAP2VECI

IMPLICIT NONE
INTEGER(KIND=JPIM),ALLOCATABLE  :: I2TEMP(:,:)
INTEGER(KIND=JPIM)              :: ISEQ, I0, I1
INTEGER(KIND=JPIM)              :: NCID,VARID,STATUS

  ALLOCATE(I2TEMP(NX,NY))
  WRITE(LOGNAM,*)'TOPO_INIT: mask_slope:',TRIM(CRIVPARNC)
  STATUS =  NF90_INQ_VARID(NCID,'mask_slope',VARID)
  IF (STATUS /= 0 ) THEN
    WRITE(LOGNAM,*)'TOPO_INIT: mask_slope: LSLPMIX should be set to FALSE: ABORTING!'
    STOP 9
  ENDIF 
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2TEMP),'reading data' ) 
  CALL MAP2VECI(I2TEMP,I2MASK)
  I0=0
  I1=0
  DO ISEQ=1,NSEQALL
    IF (I2MASK(ISEQ,1) == 1 ) THEN
      I1=I1+1
    ENDIF
    IF (I2MASK(ISEQ,1) == 0 ) THEN 
      I0=I0+1
    ENDIF
  ENDDO
  WRITE(LOGNAM,*)'TOPO_INIT: sum(mask==0), sum(mask==1)',I0,I1
  IF ( I0+I1 .NE. NSEQALL ) THEN 
     WRITE(LOGNAM,*)'TOPO_INIT: mask==0 + mask == 1 does not match NSEQALL.. something wrong, aborting'
     STOP 9
  ENDIF 

  DEALLOCATE(I2TEMP)
#endif
END SUBROUTINE SET_SLOPEMIX

END SUBROUTINE CMF_TOPO_INIT
!####################################################################


END MODULE CMF_CTRL_MAPS_MOD
