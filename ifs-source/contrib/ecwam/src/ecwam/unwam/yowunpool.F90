! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!**********************************************************************
!*                                                                    *
!**********************************************************************
      MODULE YOWUNPOOL

         USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
         USE YOWPARAM , ONLY : MDC => NANG, MSC => NFRE

         IMPLICIT NONE

         LOGICAL :: LPREPROC
         LOGICAL :: LLUNBINOUT ! controls output into erghs.bin, ergtm.bin and ergtest.bin

         REAL(KIND=JWRU), PARAMETER   :: ZERO    = 0.0_JWRU
         REAL(KIND=JWRU), PARAMETER   :: ONE     = 1.0_JWRU
         REAL(KIND=JWRU), PARAMETER   :: TWO     = 2.0_JWRU
         REAL(KIND=JWRU), PARAMETER   :: ONEHALF = 0.5_JWRU
         REAL(KIND=JWRU), PARAMETER   :: PI      = 3.141592653589793_JWRU
         REAL(KIND=JWRU), PARAMETER   :: PI2     = 2.0_JWRU*PI
         REAL(KIND=JWRU), PARAMETER   :: INVPI   = 1.0_JWRU/PI
         REAL(KIND=JWRU), PARAMETER   :: INVPI2  = 1.0_JWRU/PI2

         REAL(KIND=JWRU), PARAMETER   :: REARTH  = 6378137_JWRU ! WGS84
         REAL(KIND=JWRU), PARAMETER   :: DEGRAD  = PI/180._JWRU
         REAL(KIND=JWRU), PARAMETER   :: RADDEG  = 180._JWRU/PI

         REAL(KIND=JWRU), PARAMETER   :: THR     = TINY(1.0_4)
         REAL(KIND=JWRU), PARAMETER   :: INVTHR  = 1.0_JWRU/TINY(1.0_4)

         REAL(KIND=JWRU), PARAMETER   :: THR8    = TINY(1.0_JWRU)
         REAL(KIND=JWRU), PARAMETER   :: INVTHR8 = 1.0_JWRU/TINY(1.0_JWRU)

         REAL(KIND=JWRU), PARAMETER   :: SMALL   = 1000._JWRU * THR8
         REAL(KIND=JWRU), PARAMETER   :: LARGE   = 1.0_JWRU/SMALL

         REAL(KIND=JWRU),  PARAMETER  :: ONESIXTH= 1.0_JWRU/6.0_JWRU
         REAL(KIND=JWRU),  PARAMETER  :: ONETHIRD= 1.0_JWRU/3.0_JWRU
         REAL(KIND=JWRU),  PARAMETER  :: TWOTHIRD= 2.0_JWRU/3.0_JWRU

         INTEGER(KIND=JWIM) :: OUT_METHOD = 2
         INTEGER(KIND=JWIM) :: MaxLen
         INTEGER(KIND=JWIM) :: NIBLO_FD
         
!         INTEGER(KIND=JWIM)                :: MNP ! Number of nodes ...
!         INTEGER(KIND=JWIM)                :: MNE ! Number of elements ...
         INTEGER(KIND=JWIM)                :: MAXMNECON ! max. number of connected elements on all nodes ...

         !INTEGER(KIND=JWIM), ALLOCATABLE   :: INE(:,:) ! element connection table ..
         INTEGER(KIND=JWIM), ALLOCATABLE   :: CCON(:)  ! number of connected elements per node ....
         INTEGER(KIND=JWIM), ALLOCATABLE   :: IE_CELL(:) ! element pointer in an ordered way. points do the needed element when looping over all nodes and all connectes elements 
         INTEGER(KIND=JWIM), ALLOCATABLE   :: POS_CELL(:) ! position of the above element in the connection table ...
         INTEGER(KIND=JWIM), ALLOCATABLE   :: POS_CELL2(:,:) ! position of the above element in the connection table ...
         INTEGER(KIND=JWIM), ALLOCATABLE   :: IE_CELL2(:,:) ! element pointer in an ordered way. points do the needed element when looping over all nodes and all connectes elements
         INTEGER(KIND=JWIM), ALLOCATABLE   :: CELLVERTEX(:,:,:) ! Temp array for pointer creation
         INTEGER(KIND=JWIM), ALLOCATABLE   :: ITER_EXP(:,:) ! needed amount of iterations to fullfill stability for each direction and each freq. 
         INTEGER(KIND=JWIM), ALLOCATABLE   :: ITER_EXPD(:) ! needed amount of iterations to fullfill stability for each freq.
         INTEGER(KIND=JWIM)                :: ITER_MAXA ! max. value of the above array
         INTEGER(KIND=JWIM), ALLOCATABLE   :: IE_OUTPTS(:,:) ! ELEMENT NUMBER CORRESPONDING TO ALL STRUCTURED GRID OUTPUT POINTS

!         REAL(KIND=JWRU), ALLOCATABLE    :: XP(:), YP(:) ! x.-y. coordinates 

         REAL(KIND=JWRU),  ALLOCATABLE   :: SI(:) ! median dual patch area for each node 
         REAL(KIND=JWRU),  ALLOCATABLE   :: TRIA(:) ! triangle area of each element 
         REAL(KIND=JWRU),  ALLOCATABLE   :: IEN(:,:) ! normal vectors for each element in x and y direction, inward pointing normalized by the edge length 

         REAL(KIND=JWRU), ALLOCATABLE    :: CFLCXY(:,:) ! cfl nunmber of each direction and freq. 

         INTEGER(KIND=JWIM)                :: IVECTOR 

         INTEGER(KIND=JWIM), SAVE          :: ITER_MAX

         LOGICAL                :: LCALC    = .TRUE.  ! if true do cfl number 
         LOGICAL                :: LCFL     = .FALSE. ! later use for writing cfl numbers for 
         LOGICAL                :: LCUR     = .FALSE. ! prepared for taking currnets into account 
         LOGICAL                :: LADVTEST = .FALSE. ! switch for testing advection ... 
         LOGICAL                :: LSPHE    = .TRUE.  ! spherical coordinates 
         LOGICAL                :: LVECTOR  = .FALSE. 

!         REAL(KIND=JWRU), ALLOCATABLE    :: WBACOLD(:,:,:)
!         REAL(KIND=JWRU), ALLOCATABLE    :: DSPEC(:,:,:)

         REAL(KIND=JWRU), ALLOCATABLE    :: CG(:,:), WK(:,:) ! group vel. and wave numbers ...
         
         TYPE FILEDEF
            CHARACTER(LEN=40)   :: FNAME ! unWAM files names ... 
            INTEGER(KIND=JWIM)  :: FHNDL ! unWAM file handles ...
         END TYPE

         TYPE (FILEDEF)         :: BND,WAV,STAT,DBG,WAMOUT,GRID,XFN_HS,XFN_TM, XFN_TEST
      END MODULE YOWUNPOOL
!**********************************************************************
!*                                                                    *
!**********************************************************************
