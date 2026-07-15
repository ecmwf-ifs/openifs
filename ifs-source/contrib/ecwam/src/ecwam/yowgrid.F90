! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWGRID

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *GRIDPAR*  GENERAL GRID INFORMATION.

      REAL(KIND=JWRB)              :: DELPHI 
      REAL(KIND=JWRB), ALLOCATABLE :: DELLAM(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: SINPH(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: COSPH(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: CDR(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: SDR(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: PRQRT(:)

      INTEGER(KIND=JWIM) :: IJS
      INTEGER(KIND=JWIM) :: IJL
      INTEGER(KIND=JWIM) :: NTOTIJ
      INTEGER(KIND=JWIM) :: IJSLOC, IJLLOC, IJGLOBAL_OFFSET
      INTEGER(KIND=JWIM) :: NPROMA_WAM
      INTEGER(KIND=JWIM) :: NCHNK
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: IJFROMCHNK
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: KIJL4CHNK
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: ICHNKFROMIJ 
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: IPRMFROMIJ 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *DELPHI*    REAL      GRID INCREMENT FOR LATITUDE (METRES).
!      *DELLAM*    REAL      GRID INCREMENT FOR LONGITUDE AT EQUATOR
!                            IN METRES.
!      *SINPH*     REAL      SIN OF LATITUDE.
!      *COSPH*     REAL      COS OF LATITUDE.
!      *CDR*       REAL      ADVECTION VELOCITY IN SE-NW DIRECTION WHEN
!                            GRID ROTATED BY 45 DEGREES.
!      *SDR*       REAL      ADVECTION VELOCITY IN SW-NE DIRECTION WHEN
!                            GRID ROTATED BY 45 DEGREES.
!      *PRQRT*     REAL      PERCENTAGE OF THE ROTATED (45 degree) GRID
!                            USED FOR THE PROPAGATION (see PROPAGS1).
!                            IT IS 50% EXCEPT IF FOR STABILITY REASON
!                            IT HAS TO BE DECREASED ABOVE ~75 DEGREE

!      *IJS*       INTEGER   INDEX OF FIRST POINT ON A GIVEN PROCESSOR
!      *IJL*       INTEGER   INDEX OF LAST POINT ON A GIVEN PROCESSOR
!      *NTOTIJ*    INTEGER   TOTAL NUMBER OF GRID POINTS PER PRECESSOR (IJL-IJS+1)

!!!!!! for the unstructured part of the code, it was coded so that halo points
!!!!!! are included as part of the points on a given processor (at the end of the vector of points).
!!!!!! For output, it was necessary to introduce IJSLOC and IJLLOC to point to the points that are purely local
!      *IJSLOC*    INTEGER   INDEX OF FIRST LOCAL POINT
!      *IJLLOC*    INTEGER   INDEX OF LAST LOCAL POINT
!      *IJGLOBAL_OFFSET INTEGER OFFSET TO PLACE FIRST LOCAL POINT IN GLOBAL ARRAY

!!!!! if not unstructured then IJSLOC = IJS and IJLLOC = IJL

!      FOR OPENMP:
!      *NPROMA_WAM* INTEGER  MAX NUMBERS OF GRID POINTS PER THREAD
!      *NCHNK*      INTEGER  NUMBER OF CHUNKS OF MAXIMUM NPROMA_WAM POINTS
!      *KIJL4CHNK*  INTEGER  GET INDEX OF LAST THREAD GRID POINT FOR A GIVEN CHUNK 
!      *IJFROMCHNK* INTEGER  GET IJ FROM THREAD GRID POINT INDEX AND CHUNK INDEX (=0 if not referencing any IJ)
!      *ICHNKFROMIJ*INTEGER  GET ICHNK FROM IJ 
!      *IPRMFROMIJ *INTEGER  GET IPRM FROM IJ
!
!       FOR EACH MPI PROCESSORm THE GRID POINTS INDEX GOES FROM IJS to IJL :
!       IJS ------------   IJ    -------------  IJL
!
!       WITH OPENMP, THE IJS to IJL POINTS ARE SPLIT INTO NCHNK CHUNKS
!       CHUNK_1  CHUNK_2 .....  CHUNK_ICHNK ...... CHUNK_NCHNK

!       FOR EACH CHUNK ICHNK,THE GRID POINT INDEX GOES FROM 1 to KIJL4CHNK :
!       1  ----IPRM ---- KIJL4CHNK
!       where KIJL4CHNK is at most NPROMA_WAM

!       IJFROMCHNK MAPS  (IPRM, ICHNK)  to IJ
!       ICHNKFROMIJ MAPS IJ to ICHNK
!       IPRMFROMIJ MAPS IJ to IPRM


! ----------------------------------------------------------------------
      END MODULE YOWGRID
