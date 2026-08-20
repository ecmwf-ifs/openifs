! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MPDECOMP(NPR, MAXLEN, LLIRANK, LLWVENVI)

!****  *MPDECOMP* - EVEN DECOMPOSITION OF THE GRID DOMAIN
!****               AMONG PROCESSES AFTER INPUT OF PREPROC GRID
!                   AND UBUF FILES.

!     J. BIDLOT    ECMWF   MARCH 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   JANUARY 1998 introduce NPR 
!     J. BIDLOT    ECMWF   OCTOBER 1998 COMPLETE READING OF
!                                       IU07
!     J. BIDLOT    ECMWF   FEBRUARY 1999 TAUT --> SQRT(TAUT)
!     J. BIDLOT    ECMWF   OCTOBER 2000 NOW READING SQRT(TAUT)

!     J. BIDLOT    ECMWF   FEBRUARY 2002 NEW 2-D DECOMPOSITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT    ECMWF       11/2003
!                  IF YOUR ARE RUNNING AT ECMWF:
!                  BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                  STRUCTURE OF THE UBUF FILE YOU WILL HAVE TO
!                  MAKE SURE THAT IT IS CREATED FOR YOUR RUN,
!                  OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                  THAT IS ALREADY ON DISK.
!                  YOU ALSO HAVE TO CHANGE MUBUF.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    P. de Rosnay  ECMWF   October 2009
!                  Re-initilize for SEKF surface analysis loops
!                  (in case of coupled Jacobians only) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     PURPOSE.
!     --------
!     IT WILL READ IU07 AND IU08 AND IF MESSAGE PASSING IT WILL
!     DETERMINE AN EVEN DECOMPOSITION OF THE GRID ARRAYS FOR USE ON
!     A DISTRIBUTED MEMORY COMPUTER USING A MESSAGE PASSING PROTOCOL
!     FOR THE EXCHANGE OF INFORMATION ACROSS THE DIFFERENT PE's
!     THE LENGTH OF THE MESSAGES BETWEEN DIFFERENT PE's IS ALSO
!     COMPUTED.

!*    INTERFACE.
!     ----------
!     CALL *MPDECOMP*(NPR,MAXLEN,LLIRANK,LLWVENVI)

!      *NPR*       NUMBER OF SUBDOMAINS (USUALLY THE NUMBER OF PE'S )
!      *MAXLEN*    MAXIMUM NUMBER OF POINTS IN ANY SUB DOMAIN
!      *LLIRANK*   IF TRUE THEN SOME OF THE ARRAYS THAT ARE ONLY USED
!      *LLWVENVI*  IF TRUE THEN WVENVI WILL BE ALLOCATED AND DEFINED

!     METHOD.
!     -------

!    THE DOMAIN WILL BE SUBDIVIDED INTO NPR SUBAREAS IN SUCH A WAY
!    THAT EACH ONE WILL CONTAIN THE SAME NUMBER OF SEA POINTS (+-1).
!    THE MODEL IS FIRST DIVIDED INTO NYDECOMP LATITUDINAL BANDS
!    THIS IS REFFERED TO AS THE FIRST 1D DECOMPOSITION AND IS
!    IDENTICAL TO THE DECOMPOSITION WHICH WAS USED BEFORE. 
!    THE SECOND 2D DECOMPOSITION USED A SIMILAR PROCEDURE TO SPLIT
!    EACH BANDS INTO SUBAREAS OF EQUAL NUMBER OF POINTS.
!    THERE WILL BE NXDECOMP SUBAREAS IN THE FIRST NYCUT LATITUDINAL
!    BANDS (starting from the southern boundary)
!    AND (NXDECOMP-1) SUBAREAS IN THE REMAINING (NYDECOMP-NYCUT) BANDS
!    IN SUCH A WAY THAT 
!    NPR=NXDECOMP*NYCUT+(NYDECOMP-NYCUT)*(NXDECOMP-1) 
 
!    NXDECOMP and NYDECOMP ARE DETERMINED IN SUCH A WAY THAT THE SUB
!    AREAS ARE AS SQUARE AS POSSIBLE BY ASSUMING THAT THE GLOBAL
!    EXTEND OF THE DOMAIN IS TWICE AS LONG IN THE LONGITUDINAL DIRECTION
!    THAN IN THE LATITUDINAL DIRECTION (as is the case for the globe).

!    BASED ON THE CURRENT ADVECTION SCHEME (SEE PROPAGS),
!    THE DECOMPOSITION THEN YIELDS FOR EACH PE
!    NGBTOPE  WHICH GIVES THE TOTAL NUMBER OF NEIGHBOURING PE'S TO
!    WHICH INFORMATION FROM THE LOCAL PE IS POTENTIALLY NEEDED.
!    NTOPELST(INGB) INGB=1,NGBTOPE  THE LIST OF NEIGHBOURING PE'S ,
!    NTOPE(IPROC) IPROC=1,NPR  THE NUMBER OF POINTS FOR WHICH INFORMATION
!    HAS TO BE SENT TO PE IPROC,
!    IJTOPE(IH,IPROC) IH=1,NTOPE(IPROC), IPROC=1,NPR THE IJ INDEX OF THOSE
!    POINTS FOR WHICH INFORMATION HAS TO BE SENT TO PE IPROC,
!
!    NGBFROMPE WHICH GIVES THE TOTAL NUMBER OF NEIGHBOURING PE'S FROM
!    WHICH INFORMATION IS POTENTIALLY NEEDED,
!    NFROMPELST(KNGB) KNGB=1,NGBFROMPE  THE LIST OF NEIGHBOURING PE'S ,
!    NFROMPE(IPROC) IPROC=1,NPR  THE NUMBER OF POINTS THE LOCAL PE HAS
!    TO RECEIVE INFORMATION FROM PE IPROC,
!    NIJSTART(IPROC) IPROC=1,NPR THE IJ INDEX OF THE FIRST HALO POINT
!    OBTAINED FROM PE IPROC IN THE BUFFERS THAT PADS BOTH SIDE OF THE
!    1-D SEA POINT ARRAY OF THE LOCAL PE. THE OTHER POINTS ARE STORED IN
!    SUCCESSIVE ORDER FROM THAT IJSTART. THE ARRAYS KLON, KLAT,
!    KRLAT, KRLON,
!    BY CONVENTION, CONTRIBUTIONS FROM PE WITH PE NUMBER LESS THAN
!    THE LOCAL ONE (IRANK) ARE PUT IN THE BOTTOM BUFFER
!    (I.E. THEIR IJ'S ARE LESS THAN NSTART) AND THOSE CONTRIBUTIONS WITH
!    PE NUMBER GREATER THAN THE LOCAL ONE (IRANK) ARE PUT IN TOP BUFFER 
!    (I.E. THEIR IJ'S ARE GREATER THAN NEND).

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCES.
!     -----------
!          NONE

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWFRED  , ONLY : FR       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : IJS, IJL, NTOTIJ,                           &
     &                      NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK,   &
     &                      IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET,        &
     &                      DELLAM   ,COSPH    ,DELPHI   ,              &
     &                      CDR      ,SDR      ,PRQRT
      USE YOWMAP   , ONLY : BLK2GLO  ,BLK2LOC  ,KXLTMIN  ,KXLTMAX  ,    &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,              &
     &            KMNOP    ,KMSOP    ,NIBLO    ,NGX      ,NGY
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP     ,KTAG     ,    &
     &                      NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE_RED ,LLUNSTR  ,LL1D
      USE YOWPCONS , ONLY : G        ,PI       ,ZPI
      USE YOWSHAL  , ONLY : BATHY    ,LLOCEANMASK, WVENVI
      USE YOWSTAT  , ONLY : IPROPAGS ,LSUBGRID
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,    &
     &            NFROMPE  ,NFROMPEMAX,NTOPE   ,NTOPEMAX ,NIJSTART ,    &
     &            IJTOPE   ,NGBTOPE  ,NTOPELST ,NGBFROMPE,NFROMPELST,   &
     &            IJ2NEWIJ ,NBLKS    ,NBLKE
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR      ,             &
     &            KRLAT    ,KRLON    ,                                  &
     &            WLAT     ,WCOR     ,WRLAT    ,WRLON    ,              &
     &            OBSLAT   ,OBSLON   ,OBSCOR   ,OBSRLAT  ,OBSRLON
      USE YOWUNIT  , ONLY : IREADG   ,IU07     ,IU08     ,LWVWAMINIT
      USE YOWWIND  , ONLY : NXFFS    ,NXFFE    ,NYFFS    ,NYFFE,        &
     &                      NXFFS_LOC,NXFFE_LOC,NYFFS_LOC,NYFFE_LOC

#ifdef WAM_HAVE_UNWAM
      USE YOWPD, ONLY : MNP => npa, RANK
#endif

      USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_ALLGATHERV, MPL_SCATTERV
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWABORT , ONLY : WAM_ABORT

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "jonswap.intfb.h"
#include "mchunk.intfb.h"
#include "propconnect.intfb.h"
#include "wvwaminit.intfb.h"
#include "wam_sortini.intfb.h"
#include "wam_sorti.intfb.h"
#include "wam_nproma.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NPR
      INTEGER(KIND=JWIM), INTENT(OUT) :: MAXLEN 
      LOGICAL, INTENT(IN) :: LLIRANK
      LOGICAL, INTENT(IN) :: LLWVENVI


      INTEGER(KIND=JWIM) :: IJ, M, K, I, J, IP, IPR, IAR, IX, JSN, IH
      INTEGER(KIND=JWIM) :: IC, ICC, JC, JCS, JCM, IIL, NGOU, NH, JH, INBNGH 
      INTEGER(KIND=JWIM) :: ITAG, IREAD 
      INTEGER(KIND=JWIM) :: NLAND 
      INTEGER(KIND=JWIM) :: ICHNK, IPRM, KIJS, IJSB, KIJL, IJLB, JKGLO 
      INTEGER(KIND=JWIM) :: NPROMA, MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS
      INTEGER(KIND=JWIM) :: NTEMP(1)
      INTEGER(KIND=JWIM) :: NLENHALO_MAX
      INTEGER(KIND=JWIM) :: ICOUNTS(NPR)
      INTEGER(KIND=JWIM) :: NPLEN(NPR)
      INTEGER(KIND=JWIM) :: NGAUSSW_sekf, NLON_sekf, NLAT_sekf
      INTEGER(KIND=JWIM) :: MPLENGTH, KCOUNT, ICL, ICR, ICOUNT, IPROC
      INTEGER(KIND=JWIM) :: NXDECOMP, NYDECOMP, NYCUT
      INTEGER(KIND=JWIM) :: ISTAGGER, NIJ, NTOT, NAREA
      INTEGER(KIND=JWIM) :: NMEAN, NREST, NPTS, IXLONMAX
      INTEGER(KIND=JWIM) :: KLATBOT, KLATTOP, KXLAT, NLONGMAX, KMIN, IXLONMIN
      INTEGER(KIND=JWIM) :: MAXPERMLEN, MXPRLEN 
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: NXS, NXE, NYS, NYE 
      INTEGER(KIND=JWIM), ALLOCATABLE :: ICOMBUF_S(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: ICOMBUF_R(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IDUM(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: KDUM(:)
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: INDEX, IJNDEX, NTOTSUB
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: NEWIJ2IJ
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: NSTART1D, NEND1D
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: KSTART1, KEND1, NLON, ILON
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: IXLON
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: ITEMP, NLENHALO, IJHALO 
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: KTEMP
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: IJFROMPE, IPROCFROM
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:)   :: IJFROMPEX
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSLON, KOBSLAT
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSCOR
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSRLON, KOBSRLAT

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: RS_sekf, RN_sekf
      REAL(KIND=JWRB) :: XDELLOINV, STAGGER, XLON, SQRT2O2, A, B
      REAL(KIND=JWRB) :: THETAMAX, SINTHMAX, DELTA
      REAL(KIND=JWRB) :: X4(2)
      REAL(KIND=JWRB),ALLOCATABLE :: RCOMBUF_S(:)
      REAL(KIND=JWRB),ALLOCATABLE :: RCOMBUF_R(:)
      REAL(KIND=JWRB),ALLOCATABLE :: RDUM(:)

      CHARACTER(LEN=80) :: LOGFILENAME

      LOGICAL :: LLEXIST
      LOGICAL :: LL_sekf
      LOGICAL :: LLRNL

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MPDECOMP',0,ZHOOK_HANDLE)

!     0. READ GRID INPUT FROM PREPROC 
!        ----------------------------

!AR: ITAG was not set
ITAG =0 
IREAD=IREADG

! CALL TO *READPRE* HAS BEEN MOVED TO WVWAMINIT
! Re-initilize for SEKF surface analysis loops
IF (LWVWAMINIT) THEN
  LL_sekf = .True.
  LLRNL=.TRUE.
  CALL WVWAMINIT(LL_sekf,IU06,LLRNL, NGAUSSW_sekf, NLON_sekf,NLAT_sekf,RS_sekf,RN_sekf)
ENDIF

WRITE(IU06,*) ' WAVE MODEL PREPROC GRID INFORMATION AVAILABLE'
CALL FLUSH (IU06)


IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM

      !! UNSTRUCTURED GRID : !! (This is still very experimental !!)

!AR: augmented domain ...
        NINF = 1
        NSUP = MNP
        MAXLEN = MAXVAL(RANK(:)%NP) 
! NEW MAPPINGS BASED ON THE DECOMPOSITION ...
        IJS = 1
        IJL = MNP
        NTOTIJ = MNP

        IF (ALLOCATED(NSTART)) DEALLOCATE(NSTART)
        ALLOCATE (NSTART(NPR))
        IF (ALLOCATED(NEND)) DEALLOCATE(NEND)
        ALLOCATE (NEND(NPR))
        NSTART(:) = 1
!!! this need to be sorted nstart and nend were designed to contain
!!!! the number of point on each PE excluding the ghost node but with unstructured 
!!! we compute on ghost node so it does work !!! 
!!! for now introduce NBLKS and NBLKE to do that job
        NEND(:) = NIBLO
        NEND(IRANK)=MNP

        IJSLOC=1
        IJLLOC=RANK(IRANK)%NP
        IJGLOBAL_OFFSET=RANK(IRANK)%ISTART-1

        IF (ALLOCATED(NBLKS)) DEALLOCATE(NBLKS)
        ALLOCATE (NBLKS(NPR))
        IF (ALLOCATED(NBLKE)) DEALLOCATE(NBLKE)
        ALLOCATE (NBLKE(NPR))
        DO IP=1,NPR
          NBLKS(IP)=RANK(IP)%ISTART
          NBLKE(IP)=NBLKS(IP)+RANK(IP)%NP-1
        ENDDO

!       the data structure of type FORCING_FIELDS are defacto uni-dimensional
!       when unstructured grid is used and only limited to the local+ghost points
        NXFFS=1
        NXFFE=MNP
        NYFFS=1
        NYFFE=1
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif

ELSE
      !! NON UNSTRUCTURED GRID : !!

      NXFFS=1
      NXFFE=NGX
      NYFFS=1
      NYFFE=NGY

!*    1. NEIGHBOURING GRID POINT INDICES
!        -------------------------------

!     DETERMINE NEIGHBOURING GRID POINT

!     FIND KNMOP AND KMSOP
      KMNOP=1
      KMSOP=NGY
      DO IJ = IJS, IJL
        KMNOP=MAX(BLK2GLO%KXLT(IJ),KMNOP)
        KMSOP=MIN(BLK2GLO%KXLT(IJ),KMSOP)
      ENDDO


!*    2. FIND THE NUMBER OF POINTS PER PE, THE START AND END INDEX
!        --------------------------------------------------------------

      IF (ALLOCATED(NSTART)) DEALLOCATE(NSTART)
      ALLOCATE (NSTART(NPR))
      IF (ALLOCATED(NEND)) DEALLOCATE(NEND)
      ALLOCATE (NEND(NPR))

      IF (ALLOCATED(KLENBOT)) DEALLOCATE(KLENBOT)
      ALLOCATE (KLENBOT(NPR))
      IF (ALLOCATED(KLENTOP)) DEALLOCATE(KLENTOP)
      ALLOCATE (KLENTOP(NPR))

      IF (ALLOCATED(NFROMPE)) DEALLOCATE(NFROMPE)
      ALLOCATE (NFROMPE(NPR))
      IF (ALLOCATED(NTOPE)) DEALLOCATE(NTOPE)
      ALLOCATE (NTOPE(NPR))
      IF (ALLOCATED(NIJSTART)) DEALLOCATE(NIJSTART)
      ALLOCATE (NIJSTART(NPR))

!     DETERMINE THE SRUCTURE OF THE MODEL DECOMPOSITION

      IF (LL1D) THEN
!     1D DECOMPOSITION ONLY (old)
        NXDECOMP=1
        NYDECOMP=NPR
        NYCUT=NYDECOMP
      ELSE
!     2D DECOMPOSITION (new)
        IF (NPR == 1) THEN
          NXDECOMP=1
          NYDECOMP=1
          NYCUT=1
        ELSEIF (NPR == 2) THEN
          NXDECOMP=2
          NYDECOMP=1
          NYCUT=1
        ELSE
          IPROC=0
          ICOUNT=0
!         find whether NPR can be expressed as 2*i**2 i=1,2,3,...
!         because in that case
!         NPR=NXDECOMP*NYCUT+(NYDECOMP-NYCUT)*(NXDECOMP-1)  (1) 
!         is satisfied with NXDECOMP=2*NYDECOMP and NYCUT=NYDECOMP
!         which is a perfect subdivision into identical squares 
!         of a rectangle of dimension 2L by L.
          DO WHILE (IPROC < NPR)
            ICOUNT=ICOUNT+1
            IPROC=2*ICOUNT**2
          ENDDO
          IF (IPROC == NPR) THEN
            NYDECOMP=INT(SQRT(REAL(NPR,JWRB)/2_JWRB))
            NXDECOMP=2*NYDECOMP
            NYCUT=NYDECOMP
          ELSE
!         if the even decomposition into squares is not possible
!         start with the following approximation for NYDECOMP
!         found by setting NXDECOMP=2*NYDECOMP into (1) and
!         play around NXDECOMP=2*NYDECOMP,NYDECOMP,-1 and
!                     NYCUT=NYDECOMP,1,-1 until a solution to (1)
!         is reached. 
            IPROC=0
            NYDECOMP=INT(SQRT(REAL(NPR,JWRB)/2_JWRB))+1
            DO NXDECOMP=2*NYDECOMP,NYDECOMP,-1
              DO NYCUT=NYDECOMP,1,-1
                IPROC=NYDECOMP*(NXDECOMP-1)+NYCUT
                IF (IPROC == NPR) EXIT
              ENDDO
              IF (IPROC == NPR) EXIT
            ENDDO
            IF (IPROC /= NPR) THEN
              WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!!'
              CALL ABORT1
            ENDIF
          ENDIF
        ENDIF
      ENDIF

!     FIRST 1-D DECOMPOSITION IN LATITUDINAL BANDS

      IF (ALLOCATED(NSTART1D)) DEALLOCATE(NSTART1D)
      ALLOCATE(NSTART1D(NYDECOMP))
      IF (ALLOCATED(NEND1D)) DEALLOCATE(NEND1D)
      ALLOCATE(NEND1D(NYDECOMP))

      IF (NYCUT == NYDECOMP) THEN
!       if the number of subareas per latitunal bands is the same
!       in all bands then the number of sea points in each band will
!       be determined to be as even as possible
        NMEAN=IJL/NYDECOMP
        NREST=IJL-NMEAN*NYDECOMP

     
        NSTART1D(1)=1
        IF (NREST > 0) THEN
          NPTS=NMEAN+1
          NREST=NREST-1 
        ELSE
          NPTS=NMEAN
        ENDIF
        NEND1D(1)=NSTART1D(1)+NPTS-1

        DO IP=2,NYDECOMP
          NSTART1D(IP)=NSTART1D(IP-1)+NPTS
          IF (NREST > 0) THEN
            NPTS=NMEAN+1
            NREST=NREST-1 
          ELSE
            NPTS=NMEAN
          ENDIF
          NEND1D(IP)=NSTART1D(IP)+NPTS-1
        ENDDO
      ELSE
!       if the number of subareas per latitunal bands is not the same
!       then the number of sea points per latitudinal bands will be
!       determined in such a way that number of points for the bands with
!       the least subareas (the top (nydecomp-nycut) bands) will
!       roughly scale like (nxdecomp-1)/nxdecomp the number of points
!       in the remaining bottom nycut bands

        IF (XDELLA < 0.05_JWRB) THEN
          NMEAN=INT( REAL(IJL,JWRU) * (REAL(NXDECOMP,JWRU)/REAL(((NXDECOMP-1)*NYDECOMP+NYCUT),JWRU)) )
        ELSE
          ! this is to insure that low resolution cases behaves exactly as before.
          ! But this should be removed as soon as possible !!!!
          NMEAN=NXDECOMP*IJL/((NXDECOMP-1)*NYDECOMP+NYCUT)
        ENDIF
      
        NSTART1D(1)=1
        NPTS=NMEAN
        NEND1D(1)=NSTART1D(1)+NPTS-1

        DO IP=2,NYCUT
          NSTART1D(IP)=NSTART1D(IP-1)+NPTS
          NPTS=NMEAN
          NEND1D(IP)=NSTART1D(IP)+NPTS-1
        ENDDO

        NMEAN=(IJL-NEND1D(NYCUT))/(NYDECOMP-NYCUT)
        NREST=(IJL-NEND1D(NYCUT))-NMEAN*(NYDECOMP-NYCUT)

        DO IP=NYCUT+1,NYDECOMP
          NSTART1D(IP)=NSTART1D(IP-1)+NPTS
          IF (NREST > 0) THEN
            NPTS=NMEAN+1
            NREST=NREST-1 
          ELSE
            NPTS=NMEAN
          ENDIF
          NEND1D(IP)=NSTART1D(IP)+NPTS-1
        ENDDO
      ENDIF 


!     SECOND 1-D DECOMPOSITION IN EACH LATITUDINAL BAND

      IF (LL1D .OR. NPR == 1) THEN
        DO IP=1,NYDECOMP
          NSTART(IP)=NSTART1D(IP)
          NEND(IP)=NEND1D(IP)
        ENDDO

        ALLOCATE(NEWIJ2IJ(0:NIBLO))
        DO IJ = 0, NIBLO
          NEWIJ2IJ(IJ)=IJ
        ENDDO

        IF (ALLOCATED(IJ2NEWIJ)) DEALLOCATE(IJ2NEWIJ)
        ALLOCATE(IJ2NEWIJ(0:NIBLO))
        DO IJ = 0, NIBLO
          IJ2NEWIJ(IJ)=IJ
        ENDDO

      ELSE

        ALLOCATE(NEWIJ2IJ(0:NIBLO))
        NEWIJ2IJ(0)=0

        IF (ALLOCATED(IJ2NEWIJ)) DEALLOCATE(IJ2NEWIJ)
        ALLOCATE(IJ2NEWIJ(0:NIBLO))
        IJ2NEWIJ(0)=0

        XDELLOINV=1.0_JWRB/XDELLO


        STAGGER=0.5_JWRB*(AMOEAP-AMOWEP+IPER*XDELLO)/NXDECOMP
        STAGGER=REAL(NINT(100*STAGGER),JWRB)/100.0_JWRB
        ISTAGGER=NINT(STAGGER*XDELLOINV)

        IPROC=0
        NIJ=0
        DO IPR=1,NYDECOMP
          IPROC=IPROC+1
          NSTART(IPROC)=NIJ+1
          NTOT=NEND1D(IPR)-NSTART1D(IPR)+1

!         find number of points per subarea
          IF (IPR <= NYCUT) THEN
            NAREA=NXDECOMP
          ELSE
            NAREA=NXDECOMP-1
          ENDIF 
          ALLOCATE(NTOTSUB(NAREA))

          NMEAN=NTOT/NAREA
          NREST=NTOT-NMEAN*NAREA
          DO IAR=1,NAREA
            IF (NREST > 0) THEN
              NTOTSUB(IAR)=NMEAN+1
              NREST=NREST-1
            ELSE
              NTOTSUB(IAR)=NMEAN
            ENDIF
          ENDDO


!         sort sea points in latitudinal bands by increasing longitudes
!         and increasing latitude. Note that we use the fact that the
!         IJ's are already ordered for each latitude 
!         by determining the array IJNDEX which contain the IJ's with
!         increasing longitude and latitude.

          KLATBOT=BLK2GLO%KXLT(NSTART1D(IPR))
          KLATTOP=BLK2GLO%KXLT(NEND1D(IPR))
          ALLOCATE(KSTART1(KLATBOT:KLATTOP))
          KSTART1(:)=0
          ALLOCATE(KEND1(KLATBOT:KLATTOP))
          KEND1(:)=0
          ALLOCATE(NLON(KLATBOT:KLATTOP))
          KXLAT=KLATBOT
          KSTART1(KXLAT) = NSTART1D(IPR)
          DO IJ=NSTART1D(IPR)+1,NEND1D(IPR)
            IF (KXLAT < BLK2GLO%KXLT(IJ)) THEN
              KXLAT = BLK2GLO%KXLT(IJ)
              KSTART1(KXLAT) = IJ 
              KEND1(KXLAT-1) = IJ-1
            ENDIF
          ENDDO
          KEND1(KLATTOP)=NEND1D(IPR)

          NLONGMAX=0
          DO KXLAT=KLATBOT,KLATTOP
            NLONGMAX=MAX(KEND1(KXLAT)-KSTART1(KXLAT)+1,NLONGMAX)
          ENDDO

          ALLOCATE(IXLON(NLONGMAX,KLATBOT:KLATTOP))

          IXLONMAX=INT(AMOWEP*XDELLOINV)-1
          DO KXLAT=KLATBOT,KLATTOP
            NLON(KXLAT)=0
          ENDDO
          KXLAT=KLATBOT
          DO IJ=NSTART1D(IPR),NEND1D(IPR)
            IF (KXLAT < BLK2GLO%KXLT(IJ)) THEN
              KXLAT = BLK2GLO%KXLT(IJ)
            ENDIF
            NLON(KXLAT)=NLON(KXLAT)+1
            IX = BLK2GLO%IXLG(IJ)
            JSN= BLK2GLO%KXLT(IJ)
            XLON=AMOWEP+(IX-1)*ZDELLO(JSN)
            XLON=REAL(NINT(100*XLON),JWRB)/100.0_JWRB
            IXLON(NLON(KXLAT),KXLAT)=NINT(XLON*XDELLOINV)
            IXLONMAX=MAX(IXLONMAX,IXLON(NLON(KXLAT),KXLAT))
          ENDDO

          ALLOCATE(IJNDEX(NTOT))

          ALLOCATE(ILON(KLATBOT:KLATTOP))
          DO KXLAT=KLATBOT,KLATTOP
            ILON(KXLAT)=1
          ENDDO
          JC=0
          KMIN=KLATBOT
          DO WHILE(KMIN > 0)
            IXLONMIN=IXLONMAX+1
            KMIN=0
            DO KXLAT=KLATBOT,KLATTOP
              IF (ILON(KXLAT) <= NLON(KXLAT)) THEN
                IF (IXLON(ILON(KXLAT),KXLAT) < IXLONMIN) THEN
                  KMIN=KXLAT
                  IXLONMIN=IXLON(ILON(KXLAT),KXLAT)
                ENDIF
              ENDIF
            ENDDO
            IF (KMIN > 0) THEN
              IJ=KSTART1(KMIN)+ILON(KMIN)-1 
              JC=JC+1
              IJNDEX(JC)=IJ
              ILON(KMIN)=ILON(KMIN)+1
            ENDIF
          ENDDO

!         find which points belong to a subarea

          JCS=1
          IF (MOD(IPR,2) == 0) THEN
!         staggering
            JCM=1
            DO KXLAT=KLATBOT,KLATTOP
              IIL=1
              DO WHILE (IXLON(MIN(IIL,NLON(KXLAT)),KXLAT) < ISTAGGER   &
     &                  .AND.   IIL <= NLON(KXLAT)                     &
     &                  .AND.   NLON(KXLAT) > 0 )
                IIL=IIL+1
                JCM=JCM+1 
              ENDDO
            ENDDO
          ELSE
            JCM=1
          ENDIF

          IAR=1
          IC=0
          DO JC=JCM,NTOT 
            NIJ=NIJ+1
            IC=IC+1 
            IF (IC == NTOTSUB(IAR)) THEN
              NEND(IPROC)=NIJ
            ELSEIF (IC > NTOTSUB(IAR)) THEN
              IC=1
              IAR=IAR+1
              IPROC=IPROC+1
              NSTART(IPROC)=NIJ
            ENDIF
            IJ=IJNDEX(JC)
            NEWIJ2IJ(NIJ)=IJ
            IJ2NEWIJ(IJ)=NIJ
          ENDDO

          DO JC=JCS,JCM-1
            NIJ=NIJ+1
            IC=IC+1 
            IF (IC == NTOTSUB(IAR)) THEN
              NEND(IPROC)=NIJ
            ELSEIF (IC > NTOTSUB(IAR)) THEN
              IC=1
              IAR=IAR+1
              IPROC=IPROC+1
              NSTART(IPROC)=NIJ
            ENDIF
            IJ=IJNDEX(JC)
            NEWIJ2IJ(NIJ)=IJ
            IJ2NEWIJ(IJ)=NIJ
          ENDDO


          DEALLOCATE(KSTART1)
          DEALLOCATE(KEND1)
          DEALLOCATE(NLON)
          DEALLOCATE(ILON)
          DEALLOCATE(IXLON)
          DEALLOCATE(IJNDEX)
          DEALLOCATE(NTOTSUB)

        ENDDO

!       RELABELLING OF THE ARRAYS

        ALLOCATE(KDUM(NIBLO))
        DO NIJ=NSTART(1),NEND(NPR)
          KDUM(NIJ)=BLK2GLO%IXLG(NEWIJ2IJ(NIJ))
        ENDDO
        DO NIJ=NSTART(1),NEND(NPR)
          BLK2GLO%IXLG(NIJ)=KDUM(NIJ)
        ENDDO
        DO NIJ=NSTART(1),NEND(NPR)
          KDUM(NIJ)=BLK2GLO%KXLT(NEWIJ2IJ(NIJ))
        ENDDO
        DO NIJ=NSTART(1),NEND(NPR)
          BLK2GLO%KXLT(NIJ)=KDUM(NIJ)
        ENDDO
        DEALLOCATE(KDUM)

      ENDIF ! END IF LL1D

      DEALLOCATE(NSTART1D)
      DEALLOCATE(NEND1D)


!     COMPUTE THE CONNECTION POINTERS FOR THE PROPAGATION SCHEME
!     ----------------------------------------------------------

      IF (ALLOCATED(KLAT)) DEALLOCATE(KLAT)
      ALLOCATE(KLAT(NSTART(IRANK):NEND(IRANK),2,2))
      IF (ALLOCATED(WLAT)) DEALLOCATE(WLAT)
      ALLOCATE(WLAT(NSTART(IRANK):NEND(IRANK),2))

      IF (ALLOCATED(KLON)) DEALLOCATE(KLON)
      ALLOCATE(KLON(NSTART(IRANK):NEND(IRANK),2))

      IF (IPROPAGS == 2) THEN
        IF (ALLOCATED(KCOR)) DEALLOCATE(KCOR)
        ALLOCATE(KCOR(NSTART(IRANK):NEND(IRANK),4,2))
        IF (ALLOCATED(WCOR)) DEALLOCATE(WCOR)
        ALLOCATE(WCOR(NSTART(IRANK):NEND(IRANK),4))

      ELSEIF (IPROPAGS == 1) THEN
        IF (ALLOCATED(KRLAT)) DEALLOCATE(KRLAT)
        ALLOCATE(KRLAT(NSTART(IRANK):NEND(IRANK),2,2))
        IF (ALLOCATED(KRLON)) DEALLOCATE(KRLON)
        ALLOCATE(KRLON(NSTART(IRANK):NEND(IRANK),2,2))
        IF (ALLOCATED(WRLAT)) DEALLOCATE(WRLAT)
        ALLOCATE(WRLAT(NSTART(IRANK):NEND(IRANK),2))
        IF (ALLOCATED(WRLON)) DEALLOCATE(WRLON)
        ALLOCATE(WRLON(NSTART(IRANK):NEND(IRANK),2))

      ENDIF

      NPROMA=NPROMA_WAM
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
      DO JKGLO = NSTART(IRANK), NEND(IRANK), NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1, NEND(IRANK))

        CALL PROPCONNECT(KIJS, KIJL, NEWIJ2IJ(KIJS))
      ENDDO
!$OMP END PARALLEL DO



!     3. DETERMINE THE LENGTH OF THE MESSAGE THAT WILL BE EXCHANGED 
!        BETWEEN NEIGHBORING SUB GRID DOMAINS
!        -----------------------------------------------------------

      MAXLEN=0
      DO IP=1,NPR
        NPLEN(IP)=NEND(IP)-NSTART(IP)+1
        MAXLEN=MAX(MAXLEN,NPLEN(IP))
      ENDDO

!     FIND INDEX AND PE OF THE POINTS IN THE HALO
      MAXPERMLEN=2*MAX(MAXLEN,NGX)+12
      IF (IPROPAGS == 0) THEN
        MXPRLEN=6*MAXPERMLEN
      ELSEIF (IPROPAGS == 1) THEN
        MXPRLEN=8*MAXPERMLEN
      ELSE
        MXPRLEN=12*MAXPERMLEN
      ENDIF

      ALLOCATE(IJFROMPE(MAXPERMLEN,NPR))
      ALLOCATE(IPROCFROM(MAXPERMLEN,NPR))
      ALLOCATE(NLENHALO(NPR))

      CALL GSTATS(1497,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IP,IH)
      DO IP=1,NPR
        DO IH=1,MAXPERMLEN
          IJFROMPE(IH,IP)=0
          IPROCFROM(IH,IP)=NPR+1
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1497,1)

!     DETERMINE IJFROMPE
      NLENHALO_MAX=0

      ALLOCATE(ITEMP(MXPRLEN))
      DO IH=1,MXPRLEN
        ITEMP(IH)=0
      ENDDO
      ALLOCATE(INDEX(MXPRLEN))


      IH=0

      IP=IRANK

!       CONTRIBUTION FROM NON-ROTATED GRID (NEEDED FOR ALL OPTIONS)
      DO IC=1,2
        DO IJ=NSTART(IP),NEND(IP)

          IF ( KLON(IJ,IC) > 0 .AND.                                  &
     &         KLON(IJ,IC) <= NIBLO .AND.                             &
     &        (KLON(IJ,IC) < NSTART(IP) .OR.                          &
     &         KLON(IJ,IC) > NEND(IP)       ) ) THEN
             IH=IH+1
             IF (IH > MXPRLEN) THEN
               WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
               WRITE(IU06,*) 'MXPRLEN TOO SMALL !!!'
               CALL ABORT1
             ENDIF
             ITEMP(IH)=KLON(IJ,IC)
          ENDIF

        ENDDO ! END DO ON IJ
      ENDDO ! END DO ON IC

      DO ICL=1,2
        DO IC=1,2
          DO IJ=NSTART(IP),NEND(IP)

          IF ( KLAT(IJ,IC,ICL) > 0 .AND.                              &
     &         KLAT(IJ,IC,ICL) <= NIBLO .AND.                         &
     &        (KLAT(IJ,IC,ICL) < NSTART(IP) .OR.                      &
     &         KLAT(IJ,IC,ICL) > NEND(IP)       ) ) THEN 
             IH=IH+1
             IF (IH > MXPRLEN) THEN
               WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
               WRITE(IU06,*) 'MXPRLEN TOO SMALL !!!'
               CALL ABORT1
             ENDIF
             ITEMP(IH)=KLAT(IJ,IC,ICL)
          ENDIF

          ENDDO ! END DO ON IJ
        ENDDO ! END DO ON IC
      ENDDO ! END DO ON ICL

      IF (IPROPAGS == 2) THEN
!       CONTRIBUTION FROM CORNER GRID POINT
!       (NEEDED FOR CTU SCHEME: IPROPAGS=2)
        DO ICL=1,2
          DO ICR=1,4
            DO IJ=NSTART(IP),NEND(IP)

              IF ( KCOR(IJ,ICR,ICL) > 0 .AND.                          &
     &             KCOR(IJ,ICR,ICL) <= NIBLO .AND.                     &
     &            (KCOR(IJ,ICR,ICL) < NSTART(IP) .OR.                  &
     &             KCOR(IJ,ICR,ICL) > NEND(IP)      ) ) THEN
                 IH=IH+1
                 IF (IH > MXPRLEN) THEN
                   WRITE(IU06,*) 'MPDECOMP : DECOMPOSITION PROBLEM !'
                   WRITE(IU06,*) 'MXPRLEN TOO SMALL !!!'
                   CALL ABORT1
                 ENDIF
                 ITEMP(IH)=KCOR(IJ,ICR,ICL)
              ENDIF

            ENDDO ! END DO ON IJ
          ENDDO ! END DO ON IC
        ENDDO ! END DO ON ICL
      ELSEIF (IPROPAGS == 1) THEN
!       CONTRIBUTION FROM ROTATED GRID
!       (NEEDED FOR DUAL ROTATED SCHEME SCHEME: IPROPAGS=1)
        DO ICL=1,2
          DO IC=1,2
            DO IJ=NSTART(IP),NEND(IP)

              IF ( KRLON(IJ,IC,ICL) > 0 .AND.                          &
     &             KRLON(IJ,IC,ICL) <= NIBLO .AND.                     &
     &            (KRLON(IJ,IC,ICL) < NSTART(IP) .OR.                  &
     &             KRLON(IJ,IC,ICL) > NEND(IP)       ) ) THEN
                 IH=IH+1
                 IF (IH > MXPRLEN) THEN
                   WRITE(IU06,*) 'MPDECOMP: decomposition problem !!'
                   WRITE(IU06,*) 'MXPRLEN TOO SMALL !!!'
                   CALL ABORT1
                 ENDIF
                 ITEMP(IH)=KRLON(IJ,IC,ICL)
              ENDIF

              IF ( KRLAT(IJ,IC,ICL) > 0 .AND.                          &
     &             KRLAT(IJ,IC,ICL) <= NIBLO .AND.                     &
     &            (KRLAT(IJ,IC,ICL) < NSTART(IP) .OR.                  &
     &             KRLAT(IJ,IC,ICL) > NEND(IP)       ) ) THEN 
                 IH=IH+1
                 IF (IH > MXPRLEN) THEN
                   WRITE(IU06,*) 'MPDECOMP: decomposition problem !!'
                   WRITE(IU06,*) 'MXPRLEN TOO SMALL !!!'
                   CALL ABORT1
                 ENDIF
                 ITEMP(IH)=KRLAT(IJ,IC,ICL)
              ENDIF

            ENDDO ! END DO ON IJ
          ENDDO ! END DO ON IC
        ENDDO ! END DO ON ICL
      ENDIF

      NH=IH

      IF (NH > 1) THEN
        CALL WAM_SORTINI(ITEMP,INDEX,NH)
        CALL WAM_SORTI(ITEMP,INDEX,NH)
      ENDIF

      JH=0
      IF (NH > 1) THEN
        JH=1
        IJFROMPE(JH,IP)=ITEMP(1)
        DO IH=2,NH
          IF (ITEMP(IH) > ITEMP(IH-1)) THEN
            JH=JH+1
            IF (JH > MAXPERMLEN) THEN
              WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
              WRITE(IU06,*) 'MAXPERMLEN TOO SMALL !!!'
              WRITE(IU06,*) 'IH = ',IH,' JH = ',JH
              CALL ABORT1
            ENDIF
            IJFROMPE(JH,IP)=ITEMP(IH)
          ENDIF
        ENDDO
      ENDIF
      NLENHALO(IP)=JH

      DEALLOCATE(ITEMP) 
      DEALLOCATE(INDEX)

!     UPDATE NLENHALO OVER ALL PROCS
      ICOUNTS(:)=1

      CALL GSTATS(694,0)
      NTEMP(1)=NLENHALO(IP)
      CALL MPL_ALLGATHERV(NTEMP,NLENHALO,ICOUNTS,CDSTRING='MPDECOMP:')
      CALL GSTATS(694,1)

      NLENHALO_MAX=MAXVAL(NLENHALO(:))

!     UPDATE IJFROMPE OVER ALL PROCS
      ICOUNTS(:)=MAXPERMLEN
      ALLOCATE(IJFROMPEX(MAXPERMLEN*NPR))

      CALL GSTATS(694,0)
      CALL MPL_ALLGATHERV(IJFROMPE(:,IP),IJFROMPEX,ICOUNTS,CDSTRING='MPDECOMP:')
      CALL GSTATS(694,1)

      DO J=1,NPR
        IJFROMPE(:,J)=IJFROMPEX((J-1)*MAXPERMLEN+1:J*MAXPERMLEN)
      ENDDO
      DEALLOCATE(IJFROMPEX)


!     RESIZE ARRAYS BASED ON MAXIMUM SIZE OF HALO

      NLENHALO_MAX=MAX(1,NLENHALO_MAX)
      ALLOCATE(KTEMP(NLENHALO_MAX,NPR))

      DO IP=1,NPR
        DO IH=1,NLENHALO_MAX
          KTEMP(IH,IP)=IJFROMPE(IH,IP)
        ENDDO
      ENDDO
      DEALLOCATE(IJFROMPE)
      ALLOCATE(IJFROMPE(NLENHALO_MAX,NPR))
      DO IP=1,NPR
        DO IH=1,NLENHALO_MAX
          IJFROMPE(IH,IP)=KTEMP(IH,IP)
        ENDDO
      ENDDO

      DO IP=1,NPR
        DO IH=1,NLENHALO_MAX
          KTEMP(IH,IP)=IPROCFROM(IH,IP)
        ENDDO
      ENDDO
      DEALLOCATE(IPROCFROM)
      ALLOCATE(IPROCFROM(NLENHALO_MAX,NPR))
      DO IP=1,NPR
        DO IH=1,NLENHALO_MAX
          IPROCFROM(IH,IP)=KTEMP(IH,IP)
        ENDDO
      ENDDO

      DEALLOCATE(KTEMP)

!     DETERMINE IPROCFROM, KLENBOT and KLENTOP
      DO IP=1,NPR
        KLENBOT(IP)=0
        KLENTOP(IP)=0
        DO IH=1,NLENHALO(IP)
          DO IPROC=1,NPR
            IF (IJFROMPE(IH,IP) >= NSTART(IPROC) .AND.                 &
     &          IJFROMPE(IH,IP) <=  NEND(IPROC) ) THEN
              IPROCFROM(IH,IP)=IPROC
              IF (IPROC < IP) THEN
                KLENBOT(IP)=KLENBOT(IP)+1
              ELSEIF (IPROC > IP) THEN
                KLENTOP(IP)=KLENTOP(IP)+1
              ENDIF
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      NINF = NSTART(IRANK)-KLENBOT(IRANK)
      NSUP = NEND(IRANK)+KLENTOP(IRANK)
      NLAND = NSUP+1

!     FIND THE LOCAL NTOPE
      DO IP=1,NPR
        NTOPE(IP)=0
      ENDDO

        DO IP=1,NPR
          DO IH=1,NLENHALO(IP)
            IF (IPROCFROM(IH,IP) == IRANK) THEN
              NTOPE(IP)=NTOPE(IP)+1
            ENDIF
          ENDDO
        ENDDO

        NTOPEMAX=0
        DO IP=1,NPR
          NTOPEMAX=MAX(NTOPEMAX,NTOPE(IP))
        ENDDO



!     FIND THE LOCAL NFROMPE
      DO IP=1,NPR
        NFROMPE(IP)=0
      ENDDO

        DO IH=1,NLENHALO(IRANK)
          NFROMPE(IPROCFROM(IH,IRANK))=NFROMPE(IPROCFROM(IH,IRANK))+1
        ENDDO

        NFROMPEMAX=0
        DO IP=1,NPR
          NFROMPEMAX=MAX(NFROMPEMAX,NFROMPE(IP))
        ENDDO


!     FIND NGBTOPE AND CREATE NTOPELST

      NGBTOPE=0
      DO IP=1,NPR
        IF (NTOPE(IP) > 0) NGBTOPE=NGBTOPE+1
      ENDDO

      IF (ALLOCATED(NTOPELST)) DEALLOCATE(NTOPELST)
      ALLOCATE(NTOPELST(NGBTOPE))
      INBNGH=0
      DO IP=1,NPR
        IF (NTOPE(IP) > 0) THEN 
          INBNGH=INBNGH+1
          NTOPELST(INBNGH)=IP
        ENDIF
      ENDDO

!     FIND NGBFROMPE AND CREATE NFROMPELST 

      NGBFROMPE=0 
      DO IP=1,NPR
        IF (NFROMPE(IP) > 0) NGBFROMPE=NGBFROMPE+1
      ENDDO

      IF (ALLOCATED(NFROMPELST)) DEALLOCATE(NFROMPELST)
      ALLOCATE(NFROMPELST(MAX(1,NGBFROMPE)))
      INBNGH=0
      DO IP=1,NPR
        IF (NFROMPE(IP) > 0) THEN 
          INBNGH=INBNGH+1
          NFROMPELST(INBNGH)=IP
        ENDIF
      ENDDO


!     DETERMINE WHICH IJ's NEED TO BE SEND TO THE OTHER PE'S
        IF (ALLOCATED(IJTOPE)) DEALLOCATE(IJTOPE)
        ALLOCATE(IJTOPE(NTOPEMAX,NPR))
        DO IP=1,NPR
          DO JH=1,NTOPEMAX
            IJTOPE(JH,IP)=NLAND
          ENDDO
        ENDDO
        DO IP=1,NPR
          JH=0
          DO IH=1,NLENHALO(IP)
            IF (IPROCFROM(IH,IP) == IRANK) THEN
              JH=JH+1
              IJTOPE(JH,IP)=IJFROMPE(IH,IP)
            ENDIF
          ENDDO
        ENDDO

      ALLOCATE(IJHALO(MAX(1,NLENHALO(IRANK))))

      DO IH=1,NLENHALO(IRANK)
        IF (IPROCFROM(IH,IRANK) < IRANK) THEN
          IJHALO(IH)=NINF+IH-1
        ELSEIF (IPROCFROM(IH,IRANK) > IRANK) THEN
          IJHALO(IH)=NEND(IRANK)+IH-KLENBOT(IRANK)
        ENDIF
      ENDDO


!     CHANGE THE LOCAL ADDRESSING OF KLAT KLON KRLAT
!     KRLON DEPTH FOR POINTS IN THE HALO
!     NOTE THAT THIS IMPLIES THAT THESE ARRAYS ARE LOCAL BECAUSE THEY
!     ARE DIFFERENT IN THE HALO REGIONS
      IF ((NPR > 1) .OR. LLIRANK) THEN

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            DO IH=1,NLENHALO(IRANK)
              IF (KLON(IJ,IC) == IJFROMPE(IH,IRANK)) THEN
                KLON(IJ,IC)=IJHALO(IH)
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        DO ICL=1,2
          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              DO IH=1,NLENHALO(IRANK)
                IF (KLAT(IJ,IC,ICL) == IJFROMPE(IH,IRANK)) THEN
                  KLAT(IJ,IC,ICL)=IJHALO(IH)
                  EXIT
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        IF (IPROPAGS == 2) THEN
!       (NEEDED FOR CTU SCHEME: IPROPAGS=2)
          DO ICL=1,2
            DO ICR=1,4
              DO IJ=NSTART(IRANK),NEND(IRANK)
                DO IH=1,NLENHALO(IRANK)
                  IF (KCOR(IJ,ICR,ICL) == IJFROMPE(IH,IRANK)) THEN
                    KCOR(IJ,ICR,ICL)=IJHALO(IH)
                    EXIT
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

        ELSEIF (IPROPAGS == 1) THEN
!       (NEEDED FOR DUAL ROTATED SCHEME SCHEME: IPROPAGS=1)
          DO ICL=1,2
            DO IC=1,2
              DO IJ=NSTART(IRANK),NEND(IRANK)
                DO IH=1,NLENHALO(IRANK)
                  IF (KRLON(IJ,IC,ICL) == IJFROMPE(IH,IRANK)) THEN
                    KRLON(IJ,IC,ICL)=IJHALO(IH)
                    EXIT
                  ENDIF
                ENDDO
                DO IH=1,NLENHALO(IRANK)
                  IF (KRLAT(IJ,IC,ICL) == IJFROMPE(IH,IRANK)) THEN
                    KRLAT(IJ,IC,ICL)=IJHALO(IH)
                    EXIT
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!     FIND NIJSTART

      DO IP=1,NPR
        NIJSTART(IP)=NLAND
      ENDDO
      IF ((NPR > 1) .OR. LLIRANK) THEN
        IF (IPROCFROM(1,IRANK) < IRANK) THEN
          NIJSTART(IPROCFROM(1,IRANK))=NINF
        ELSEIF (IPROCFROM(1,IRANK) > IRANK) THEN
          NIJSTART(IPROCFROM(1,IRANK))=NEND(IRANK)+1
        ENDIF
        DO IH=2,NLENHALO(IRANK)
          IF (IPROCFROM(IH,IRANK) /= IPROCFROM(IH-1,IRANK)) THEN
            IF (IPROCFROM(IH,IRANK) < IRANK) THEN
              NIJSTART(IPROCFROM(IH,IRANK))=NINF+IH-1
            ELSEIF (IPROCFROM(IH,IRANK) > IRANK) THEN
              NIJSTART(IPROCFROM(IH,IRANK))=NEND(IRANK)+IH-KLENBOT(IRANK) 
            ENDIF
          ENDIF
        ENDDO
      ENDIF

!     CDR AND SDR ARE ALSO USED IN PROPAGS1.
!     Note that they all will be defined locally for the halo points !!!

      IF (IPROPAGS == 1) THEN
        IF (ALLOCATED(CDR)) DEALLOCATE(CDR)
        ALLOCATE(CDR(NINF:NSUP+1,NANG))
        IF (ALLOCATED(SDR)) DEALLOCATE(SDR)
        ALLOCATE(SDR(NINF:NSUP+1,NANG))
!       PRQRT IS NOT DEFINED IN THE HALO !!!!
        IF (ALLOCATED(PRQRT)) DEALLOCATE(PRQRT)
        ALLOCATE(PRQRT(NSTART(IRANK):NEND(IRANK)))

        SQRT2O2=SIN(0.25_JWRB*PI)

        DO K=1,NANG
          CDR(NLAND,K) = 0.0_JWRB
          SDR(NLAND,K) = 0.0_JWRB
        ENDDO
        IF ( NPR > 1 ) THEN
          DO IJ=NSTART(IRANK),NEND(IRANK)
            JH = BLK2GLO%KXLT(IJ)
            DO K=1,NANG
              A=SQRT2O2*COSTH(K)*COSPH(JH)
              B=SQRT2O2*SINTH(K) 
              CDR(IJ,K)=A-B
              SDR(IJ,K)=A+B
            ENDDO
          ENDDO 
          IF (IRGG == 1) THEN
            DO IJ=NSTART(IRANK),NEND(IRANK)
              JH = BLK2GLO%KXLT(IJ)
              THETAMAX=ATAN2(1.0_JWRB,COSPH(JH))
              SINTHMAX=SIN(THETAMAX)
              DELTA=COSPH(JH)/(SINTHMAX*(1.0_JWRB+COSPH(JH)**2))
              PRQRT(IJ)=MIN(0.5_JWRB,DELTA)
            ENDDO
          ELSE
            DO IJ=NSTART(IRANK),NEND(IRANK)
              PRQRT(IJ)=0.5_JWRB
            ENDDO
          ENDIF
          DO IH=1,NLENHALO(IRANK)
            IJ=IJFROMPE(IH,IRANK)
            JH = BLK2GLO%KXLT(IJ)
            IJ=IJHALO(IH)
            DO K=1,NANG
              A=SQRT2O2*COSTH(K)*COSPH(JH)
              B=SQRT2O2*SINTH(K) 
              CDR(IJ,K)=A-B
              SDR(IJ,K)=A+B
            ENDDO
          ENDDO
        ELSE
          DO IJ=NINF,NSUP
            JH = BLK2GLO%KXLT(IJ)
            DO K=1,NANG
              A=SQRT2O2*COSTH(K)*COSPH(JH)
              B=SQRT2O2*SINTH(K) 
              CDR(IJ,K)=A-B
              SDR(IJ,K)=A+B
            ENDDO
          ENDDO 
          IF (IRGG == 1) THEN
            DO IJ=NSTART(IRANK),NEND(IRANK)
              JH = BLK2GLO%KXLT(IJ)
              THETAMAX=ATAN2(1.0_JWRB,COSPH(JH))
              SINTHMAX=SIN(THETAMAX)
              DELTA=COSPH(JH)/(SINTHMAX*(1.0_JWRB+COSPH(JH)**2))
              PRQRT(IJ)=MIN(0.5_JWRB,DELTA)
            ENDDO
          ELSE
            DO IJ=NSTART(IRANK),NEND(IRANK)
              PRQRT(IJ)=0.5_JWRB
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      DEALLOCATE(IJFROMPE)
      DEALLOCATE(IPROCFROM)
      DEALLOCATE(NLENHALO)
      DEALLOCATE(IJHALO)


!     5. MODIFY KLAT AND KLON SUCH THAT POINT INDICES FOR LAND IS NLAND
!        ---------------------------------------------------------------

        DO ICL=1,2
          DO IC=1,2
            DO IJ = NSTART(IRANK),NEND(IRANK)
              IF (KLAT(IJ,IC,ICL) == 0) KLAT(IJ,IC,ICL) = NLAND 
            ENDDO
          ENDDO
        ENDDO
        DO IC=1,2
          DO IJ = NSTART(IRANK),NEND(IRANK)
            IF (KLON(IJ,IC) == 0) KLON(IJ,IC) = NLAND 
          ENDDO
        ENDDO
        IF (IPROPAGS == 2) THEN
          DO ICL=1,2
            DO ICR=1,4
              DO IJ = NSTART(IRANK),NEND(IRANK)
                IF (KCOR(IJ,ICR,ICL) == 0) KCOR(IJ,ICR,ICL) = NLAND
              ENDDO
            ENDDO
          ENDDO
        ELSEIF (IPROPAGS == 1) THEN
          DO ICL=1,2
            DO IC=1,2
              DO IJ = NSTART(IRANK),NEND(IRANK)
                IF (KRLAT(IJ,IC,ICL) == 0) KRLAT(IJ,IC,ICL) = NLAND
                IF (KRLON(IJ,IC,ICL) == 0) KRLON(IJ,IC,ICL) = NLAND
              ENDDO
            ENDDO
          ENDDO
        ENDIF


!       READ IU08
!       =========

        IF (.NOT.LL1D .AND. NPR > 1) ALLOCATE(KDUM(NIBLO))

        MPLENGTH=MAXLEN
        ALLOCATE(ICOMBUF_S(MPLENGTH*NPR))
        ALLOCATE(ICOMBUF_R(MPLENGTH))

        ALLOCATE(IDUM(NIBLO))

        ALLOCATE(KOBSLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        ALLOCATE(KOBSLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        IF (IPROPAGS == 1) THEN
          ALLOCATE(KOBSRLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
          ALLOCATE(KOBSRLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        ENDIF
        IF (IPROPAGS == 2) THEN
          ALLOCATE(KOBSCOR(NSTART(IRANK):NEND(IRANK),NFRE_RED,4))
        ENDIF

        DO M=1,NFRE_RED ! loop over frequencies

!       READING KOBSLAT
!       ---------------
          DO IC=1,2
            ITAG=ITAG+1
            IF (IRANK == IREAD) THEN

              CALL GSTATS(1771,0)
              READ (IU08(IPROPAGS)) IDUM 
              CALL GSTATS(1771,1)

!             RELABELLING OF THE ARRAY
              IF (.NOT.LL1D .AND. NPR > 1) THEN
                DO NIJ=NSTART(1),NEND(NPR)
                  KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                  IDUM(NIJ)=KDUM(NIJ)
                ENDDO
              ENDIF

!             FILL THE SEND BUFFER
              DO IP=1,NPR
                KCOUNT=(IP-1)*MPLENGTH
                DO IJ=NSTART(IP),NEND(IP)
                  KCOUNT=KCOUNT+1
                  ICOMBUF_S(KCOUNT)=IDUM(IJ)
                ENDDO
              ENDDO
            ENDIF

            IF (NPR > 1) THEN
              CALL GSTATS(694,0)
              ICOUNTS(:)=MPLENGTH
              CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD,                  &
     &          KSENDBUF=ICOMBUF_S,                                     &
     &          KSENDCOUNTS=ICOUNTS,CDSTRING='MPDECOMP 1:')
              CALL GSTATS(694,1)
            ENDIF

!           KEEP THE RELEVANT PART OF KOBSLAT
            IF (IRANK == IREAD) THEN
              KCOUNT=(IRANK-1)*MPLENGTH
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLAT(IJ,M,IC)=ICOMBUF_S(KCOUNT)
              ENDDO
            ELSE
              KCOUNT=0
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLAT(IJ,M,IC)=ICOMBUF_R(KCOUNT)
              ENDDO
            ENDIF

          ENDDO

!       READING KOBSLON
!       ---------------

          DO IC=1,2

            ITAG=ITAG+1
            IF (IRANK == IREAD) THEN
              CALL GSTATS(1771,0)
              READ (IU08(IPROPAGS)) IDUM 
              CALL GSTATS(1771,1)

!             RELABELLING OF THE ARRAY
              IF (.NOT.LL1D .AND. NPR > 1) THEN
                DO NIJ=NSTART(1),NEND(NPR)
                   KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                  IDUM(NIJ)=KDUM(NIJ)
                ENDDO
              ENDIF
!             FILL THE SEND BUFFER
              DO IP=1,NPR
                KCOUNT=(IP-1)*MPLENGTH
                DO IJ=NSTART(IP),NEND(IP)
                  KCOUNT=KCOUNT+1
                  ICOMBUF_S(KCOUNT)=IDUM(IJ)
                ENDDO
              ENDDO
            ENDIF

            IF (NPR > 1) THEN
              CALL GSTATS(694,0)
              ICOUNTS(:)=MPLENGTH
              CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD,                  &
     &          KSENDBUF=ICOMBUF_S,                                     &
     &          KSENDCOUNTS=ICOUNTS,CDSTRING='MPDECOMP 2:')
              CALL GSTATS(694,1)
            ENDIF

!           KEEP THE RELEVANT PART OF KOBSLON
            IF (IRANK == IREAD) THEN
              KCOUNT=(IRANK-1)*MPLENGTH
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLON(IJ,M,IC)=ICOMBUF_S(KCOUNT)
              ENDDO
            ELSE
              KCOUNT=0
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLON(IJ,M,IC)=ICOMBUF_R(KCOUNT)
              ENDDO
            ENDIF

          ENDDO


!       READING KOBSRLAT
!       ----------------
          IF (IPROPAGS == 1) THEN

            DO IC=1,2
              ITAG=ITAG+1

              IF (IRANK == IREAD) THEN
                CALL GSTATS(1771,0)
                READ (IU08(IPROPAGS)) IDUM 
                CALL GSTATS(1771,1)

!               RELABELLING OF THE ARRAY
                IF (.NOT.LL1D .AND. NPR > 1) THEN
                  DO NIJ=NSTART(1),NEND(NPR)
                    KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    IDUM(NIJ)=KDUM(NIJ)
                  ENDDO
                ENDIF

!               FILL THE SEND BUFFER
                DO IP=1,NPR
                  KCOUNT=(IP-1)*MPLENGTH
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ICOMBUF_S(KCOUNT)=IDUM(IJ)
                  ENDDO
                ENDDO
              ENDIF

              IF (NPR > 1) THEN
                CALL GSTATS(694,0)
                ICOUNTS(:)=MPLENGTH
                CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD,                &
     &            KSENDBUF=ICOMBUF_S,                                   &
     &            KSENDCOUNTS=ICOUNTS,CDSTRING='MPDECOMP 3:')
                CALL GSTATS(694,1)
              ENDIF

!             KEEP THE RELEVANT PART OF KOBSRLAT
              IF (IRANK == IREAD) THEN
                KCOUNT=(IRANK-1)*MPLENGTH
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLAT(IJ,M,IC)=ICOMBUF_S(KCOUNT)
                ENDDO
              ELSE
                KCOUNT=0
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLAT(IJ,M,IC)=ICOMBUF_R(KCOUNT)
                ENDDO
              ENDIF

            ENDDO
          ENDIF

!       READING KOBSRLON
!       ----------------
          IF (IPROPAGS == 1) THEN

            DO IC=1,2
              ITAG=ITAG+1

              IF (IRANK == IREAD) THEN
                CALL GSTATS(1771,0)
                READ (IU08(IPROPAGS)) IDUM 
                CALL GSTATS(1771,1)

!               RELABELLING OF THE ARRAY
                IF (.NOT.LL1D .AND. NPR > 1) THEN
                  DO NIJ=NSTART(1),NEND(NPR)
                     KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    IDUM(NIJ)=KDUM(NIJ)
                  ENDDO
                ENDIF

!               FILL THE SEND BUFFER
                DO IP=1,NPR
                  KCOUNT=(IP-1)*MPLENGTH
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ICOMBUF_S(KCOUNT)=IDUM(IJ)
                  ENDDO
                ENDDO
              ENDIF

              IF (NPR > 1) THEN
                CALL GSTATS(694,0)
                ICOUNTS(:)=MPLENGTH
                CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD,                &
     &            KSENDBUF=ICOMBUF_S,                                   &
     &            KSENDCOUNTS=ICOUNTS,CDSTRING='MPDECOMP 4:')
                CALL GSTATS(694,1)
              ENDIF

!             KEEP THE RELEVANT PART OF KOBSRLON
              IF (IRANK == IREAD) THEN
                KCOUNT=(IRANK-1)*MPLENGTH
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLON(IJ,M,IC)=ICOMBUF_S(KCOUNT)
                ENDDO
              ELSE
                KCOUNT=0
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLON(IJ,M,IC)=ICOMBUF_R(KCOUNT)
                ENDDO
              ENDIF

            ENDDO
          ENDIF


!       READING KOBSCOR
!       ---------------
          IF (IPROPAGS == 2) THEN
            DO IC=1,4
              ITAG=ITAG+1

              IF (IRANK == IREAD) THEN
                CALL GSTATS(1771,0)
                READ (IU08(IPROPAGS)) IDUM 
                CALL GSTATS(1771,1)

!               RELABELLING OF THE ARRAY
                IF (.NOT.LL1D .AND. NPR > 1 ) THEN
                  DO NIJ=NSTART(1),NEND(NPR)
                   KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    IDUM(NIJ)=KDUM(NIJ)
                  ENDDO
                ENDIF

!               FILL THE SEND BUFFER
                DO IP=1,NPR
                  KCOUNT=(IP-1)*MPLENGTH
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ICOMBUF_S(KCOUNT)=IDUM(IJ)
                  ENDDO
                ENDDO
              ENDIF

              IF (NPR > 1) THEN
                CALL GSTATS(694,0)
                ICOUNTS(:)=MPLENGTH
                CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD,                &
     &            KSENDBUF=ICOMBUF_S,                                   &
     &            KSENDCOUNTS=ICOUNTS,CDSTRING='MPDECOMP 5:')
                CALL GSTATS(694,1)
              ENDIF           

!             KEEP THE RELEVANT PART OF KOBSCOR
              IF (IRANK == IREAD) THEN
                KCOUNT=(IRANK-1)*MPLENGTH
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSCOR(IJ,M,IC)=ICOMBUF_S(KCOUNT)
                ENDDO
              ELSE
                KCOUNT=0
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSCOR(IJ,M,IC)=ICOMBUF_R(KCOUNT)
                ENDDO
              ENDIF
            ENDDO

          ENDIF

        ENDDO ! end loop on frequencies

        IF (ALLOCATED(ICOMBUF_S)) DEALLOCATE(ICOMBUF_S)
        IF (ALLOCATED(ICOMBUF_R)) DEALLOCATE(ICOMBUF_R)

        DEALLOCATE(IDUM)

        IF (.NOT.LL1D .AND. NPR > 1 ) DEALLOCATE(KDUM)

      WRITE(IU06,*) ' WAVE MODEL PREPROC UBUF INFORMATION READ IN  (second part)'
      CALL FLUSH (IU06)

      DEALLOCATE(NEWIJ2IJ)


!     OBSTRUCTION COEFFICIENTS

!     NOTE: THE VALUE OF OBSLON WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      IF (ALLOCATED(OBSLON)) DEALLOCATE(OBSLON)
      ALLOCATE(OBSLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
!     NOTE: THE VALUE OF OBSLAT WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      IF (ALLOCATED(OBSLAT)) DEALLOCATE(OBSLAT)
      ALLOCATE(OBSLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))

      CALL GSTATS(1497,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICC,IC,M,IJ)
      DO ICC=1,2
        IC=ICC
        DO M=1,NFRE_RED
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF (.NOT. LSUBGRID) THEN
              OBSLON(IJ,M,IC)=1.0_JWRB
            ELSEIF (KOBSLON(IJ,M,IC) == 0) THEN
              OBSLON(IJ,M,IC)=0.0_JWRB
            ELSEIF (MOD(KOBSLON(IJ,M,IC),1000) == 0) THEN
              OBSLON(IJ,M,IC)=1.0_JWRB
            ELSE
              OBSLON(IJ,M,IC)=REAL(KOBSLON(IJ,M,IC),JWRB)*0.001_JWRB
            ENDIF
          ENDDO
        ENDDO

        DO M=1,NFRE_RED
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF (.NOT. LSUBGRID) THEN
              OBSLAT(IJ,M,IC)=1.0_JWRB
            ELSEIF (KOBSLAT(IJ,M,IC) == 0) THEN
              OBSLAT(IJ,M,IC)=0.0_JWRB
            ELSEIF (MOD(KOBSLAT(IJ,M,IC),1000) == 0) THEN
              OBSLAT(IJ,M,IC)=1.0_JWRB
            ELSE
              OBSLAT(IJ,M,IC)=REAL(KOBSLAT(IJ,M,IC),JWRB)*0.001_JWRB
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1497,1)

      DEALLOCATE(KOBSLON)
      DEALLOCATE(KOBSLAT)

      IF (IPROPAGS == 1) THEN
!       NOTE: THE VALUE OF OBSRLON WILL NOT BE RESET IN THE FIRST
        IF (ALLOCATED(OBSRLON)) DEALLOCATE(OBSRLON)
        ALLOCATE(OBSRLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
!       NOTE: THE VALUE OF OBSRLAT WILL NOT BE RESET IN THE FIRST
        IF (ALLOCATED(OBSRLAT)) DEALLOCATE(OBSRLAT)
        ALLOCATE(OBSRLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        CALL GSTATS(1497,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICC,IC,M,IJ)
        DO ICC=1,2
          IC=ICC
          DO M=1,NFRE_RED
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF (.NOT. LSUBGRID) THEN
                OBSRLON(IJ,M,IC)=1.0_JWRB
              ELSEIF (KOBSRLON(IJ,M,IC) == 0) THEN
                OBSRLON(IJ,M,IC)=0.0_JWRB
              ELSEIF (MOD(KOBSRLON(IJ,M,IC),1000) == 0) THEN
                OBSRLON(IJ,M,IC)=1.0_JWRB
              ELSE
                OBSRLON(IJ,M,IC)=REAL(KOBSRLON(IJ,M,IC),JWRB)*0.001_JWRB
              ENDIF
            ENDDO
          ENDDO

          DO M=1,NFRE_RED
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF (.NOT. LSUBGRID) THEN
                OBSRLAT(IJ,M,IC)=1.0_JWRB
              ELSEIF (KOBSRLAT(IJ,M,IC) == 0) THEN
                OBSRLAT(IJ,M,IC)=0.0_JWRB
               ELSEIF (MOD(KOBSRLAT(IJ,M,IC),1000) == 0) THEN
              OBSRLAT(IJ,M,IC)=1.0_JWRB
              ELSE
                OBSRLAT(IJ,M,IC)=REAL(KOBSRLAT(IJ,M,IC),JWRB)*0.001_JWRB
              ENDIF
            ENDDO
          ENDDO

        ENDDO
!$OMP END PARALLEL DO
        CALL GSTATS(1497,1)
        DEALLOCATE(KOBSRLON)
        DEALLOCATE(KOBSRLAT)
      ENDIF

      IF (IPROPAGS == 2) THEN
        IF (ALLOCATED(OBSCOR)) DEALLOCATE(OBSCOR)
        ALLOCATE(OBSCOR(NSTART(IRANK):NEND(IRANK),NFRE_RED,4))
        CALL GSTATS(1497,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICC,IC,M,IJ)
        DO ICC=1,4
          IC=ICC
          DO M=1,NFRE_RED
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF (.NOT. LSUBGRID) THEN
                OBSCOR(IJ,M,IC)=1.0_JWRB
              ELSEIF (KOBSCOR(IJ,M,IC) == 0) THEN
                OBSCOR(IJ,M,IC)=0.0_JWRB
              ELSEIF (MOD(KOBSCOR(IJ,M,IC),1000) == 0) THEN
                OBSCOR(IJ,M,IC)=1.0_JWRB
              ELSE
                OBSCOR(IJ,M,IC)=REAL(KOBSCOR(IJ,M,IC),JWRB)*0.001_JWRB
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        CALL GSTATS(1497,1)
        DEALLOCATE(KOBSCOR)
      ENDIF

!     7. DETERMINE KXLTMIN, KXLTMAX

      IF (ALLOCATED(KXLTMIN)) DEALLOCATE(KXLTMIN)
      ALLOCATE(KXLTMIN(NPR))
      IF (ALLOCATED(KXLTMAX)) DEALLOCATE(KXLTMAX)
      ALLOCATE(KXLTMAX(NPR))

      DO IP=1,NPR
        KXLTMIN(IP)=NGY
        KXLTMAX(IP)=1
        DO IJ=NSTART(IP),NEND(IP)
           KXLTMIN(IP)=MIN(KXLTMIN(IP),BLK2GLO%KXLT(IJ))
           KXLTMAX(IP)=MAX(KXLTMAX(IP),BLK2GLO%KXLT(IJ))
        ENDDO
      ENDDO

!     8. TEST AND RESET
!        -------------- 

! SPECIFY THE LIMITS OF THE SEGMENT OVER WHICH THE PE HAS DIRECT ACCESS
      IJS = NSTART(IRANK)
      IJL = NEND(IRANK)
      NTOTIJ = IJL-IJS+1

      IJSLOC = NSTART(IRANK)
      IJLLOC = NEND(IRANK)
      IJGLOBAL_OFFSET = 0

      IF (ALLOCATED(NBLKS)) DEALLOCATE(NBLKS)
      ALLOCATE (NBLKS(NPR))
      NBLKS(:)=NSTART(:)
      IF (ALLOCATED(NBLKE)) DEALLOCATE(NBLKE)
      ALLOCATE (NBLKE(NPR))
      NBLKE(:)=NEND(:)

      KTAG=KTAG+1

      IF (IRANK == IREAD) CLOSE (UNIT=IU08(IPROPAGS))
      ! For the SEKF surface analysis
      LWVWAMINIT=.TRUE.


ENDIF ! LLUNSTR


!! CREATE THE NPROMA CHUNKS STRUCTURES
!  -----------------------------------
! ADJUST NPROMA_WAM
MTHREADS=1
!$    MTHREADS=OMP_GET_MAX_THREADS()
NPROMA=NPROMA_WAM
CALL WAM_NPROMA(IJS, IJL, MTHREADS, NPROMA)
NPROMA_WAM=NPROMA

CALL MCHUNK



! CREATE IFROMIJ, JFROMIJ and WVENVI
! !!!! IT IS ONLY DEFINED FOR GRID POINTS ON A GIVEN PE  !!!!

IF (ALLOCATED(BLK2LOC%IFROMIJ)) THEN
    CALL BLK2LOC%DEALLOC()
ENDIF
CALL BLK2LOC%ALLOC(NPROMA_WAM, NCHNK)

ALLOCATE(NXS(NCHNK))
ALLOCATE(NXE(NCHNK))
ALLOCATE(NYS(NCHNK))
ALLOCATE(NYE(NCHNK))

IF (LLUNSTR) THEN

  ! the use of reduced values for the dimension of FIELDG not yet implemented
  NXFFS_LOC = NXFFS 
  NXFFE_LOC = NXFFE
  NYFFS_LOC = NYFFS
  NYFFE_LOC = NYFFE

  DO ICHNK = 1, NCHNK
    DO IPRM = 1, NPROMA_WAM 
      IJ = IJFROMCHNK(IPRM, ICHNK)
      IF (IJ > 0) THEN
        BLK2LOC%IFROMIJ(IPRM,ICHNK) = IJ
        BLK2LOC%KFROMIJ(IPRM,ICHNK) = 1
        BLK2LOC%JFROMIJ(IPRM,ICHNK) = 1
      ELSE
!!!     these are fictious points but will point to the first point in the chunk as it should always exist 
        BLK2LOC%IFROMIJ(IPRM,ICHNK) = BLK2LOC%IFROMIJ(1,ICHNK)
        BLK2LOC%KFROMIJ(IPRM,ICHNK) = BLK2LOC%KFROMIJ(1,ICHNK)
        BLK2LOC%JFROMIJ(IPRM,ICHNK) = BLK2LOC%JFROMIJ(1,ICHNK)
      ENDIF
    ENDDO
  ENDDO

ELSE

  IF (LLWVENVI) THEN
    IF (ALLOCATED(WVENVI%UCUR))THEN
        CALL WVENVI%DEALLOC()
    ENDIF
    CALL WVENVI%ALLOC(NPROMA_WAM, NCHNK)
  ENDIF

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(KIJS, IJSB, KIJL, IJLB, IJ, JH)
  DO ICHNK = 1, NCHNK
    KIJS = 1
    IJSB = IJFROMCHNK(KIJS, ICHNK)
    KIJL = KIJL4CHNK(ICHNK)
    IJLB = IJFROMCHNK(KIJL, ICHNK)

    BLK2LOC%IFROMIJ(KIJS:KIJL,ICHNK) = BLK2GLO%IXLG(IJSB:IJLB)
    BLK2LOC%KFROMIJ(KIJS:KIJL,ICHNK) = BLK2GLO%KXLT(IJSB:IJLB)
    BLK2LOC%JFROMIJ(KIJS:KIJL,ICHNK) = NGY - BLK2GLO%KXLT(IJSB:IJLB) + 1


    NXS(ICHNK) = NGX  ! initialisation (see below)
    NXE(ICHNK) = 1
    NYS(ICHNK) = NGY
    NYE(ICHNK) = 1
    DO IJ = KIJS, KIJL
      NXS(ICHNK) = MIN(NXS(ICHNK), BLK2LOC%IFROMIJ(IJ,ICHNK))
      NXE(ICHNK) = MAX(NXE(ICHNK), BLK2LOC%IFROMIJ(IJ,ICHNK))
      NYS(ICHNK) = MIN(NYS(ICHNK), BLK2LOC%JFROMIJ(IJ,ICHNK))
      NYE(ICHNK) = MAX(NYE(ICHNK), BLK2LOC%JFROMIJ(IJ,ICHNK))
    ENDDO

    IF (LLWVENVI) THEN
      DO IJ = KIJS, KIJL
        JH = BLK2LOC%KFROMIJ(IJ,ICHNK)
        WVENVI%COSPHM1(IJ,ICHNK) = 1.0_JWRB/COSPH(JH)
        WVENVI%DELLAM1(IJ,ICHNK) = 1.0_JWRB/DELLAM(JH)
      ENDDO
    ENDIF
 
    IF (LLWVENVI) THEN
      DO IJ = KIJS, KIJL
        WVENVI%DEPTH(IJ,ICHNK) = BATHY( BLK2LOC%IFROMIJ(IJ,ICHNK), BLK2LOC%KFROMIJ(IJ,ICHNK) )
      ENDDO

!!!!     when this is moved to reading depth as an input field, wvenvi should be allocated and intialised in wvalloc !!!
      WVENVI%UCUR(KIJS:KIJL,ICHNK) = 0.0_JWRB
      WVENVI%VCUR(KIJS:KIJL,ICHNK) = 0.0_JWRB
    ENDIF

!!! these are fictious points but will point to the first point in the chunk as it should always exist 
    IF (KIJL < NPROMA_WAM) THEN
      BLK2LOC%IFROMIJ(KIJL+1:NPROMA_WAM,ICHNK) = BLK2LOC%IFROMIJ(1,ICHNK)
      BLK2LOC%KFROMIJ(KIJL+1:NPROMA_WAM,ICHNK) = BLK2LOC%KFROMIJ(1,ICHNK)
      BLK2LOC%JFROMIJ(KIJL+1:NPROMA_WAM,ICHNK) = BLK2LOC%JFROMIJ(1,ICHNK)

      IF (LLWVENVI) THEN
        WVENVI%COSPHM1(KIJL+1:NPROMA_WAM,ICHNK) = WVENVI%COSPHM1(1,ICHNK)
        WVENVI%DELLAM1(KIJL+1:NPROMA_WAM,ICHNK) = WVENVI%DELLAM1(1,ICHNK)
 
        WVENVI%DEPTH(KIJL+1:NPROMA_WAM,ICHNK) = WVENVI%DEPTH(1,ICHNK)
        WVENVI%UCUR(KIJL+1:NPROMA_WAM,ICHNK) = WVENVI%UCUR(1,ICHNK)
        WVENVI%VCUR(KIJL+1:NPROMA_WAM,ICHNK) = WVENVI%VCUR(1,ICHNK)
      ENDIF
    ENDIF

  ENDDO
!$OMP   END PARALLEL DO

  NXFFS_LOC = MINVAL(NXS(:))
  NXFFE_LOC = MAXVAL(NXE(:))
  NYFFS_LOC = MINVAL(NYS(:))
  NYFFE_LOC = MAXVAL(NYE(:))


ENDIF ! LLUNSTR

IF(ALLOCATED(BATHY)) DEALLOCATE(BATHY)
IF(ALLOCATED(LLOCEANMASK)) DEALLOCATE(LLOCEANMASK)

DEALLOCATE(NXS)
DEALLOCATE(NXE)
DEALLOCATE(NYS)
DEALLOCATE(NYE)

WRITE(IU06,*) ' WAVE MODEL DECOMPOSITION FINISHED.'
WRITE(IU06,*) ''
CALL FLUSH(IU06)

IF (LHOOK) CALL DR_HOOK('MPDECOMP',1,ZHOOK_HANDLE)

END SUBROUTINE MPDECOMP
