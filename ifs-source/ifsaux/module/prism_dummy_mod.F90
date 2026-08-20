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

!**** PRISM_DUMMY_MOD.F90
!
!     Purpose.
!     --------
!     Define interface and dummy routines needed to compile and link ifsMASTER 
!     with OASIS4 interfaces (ifs/prism/couplo4*) without linking 
!     to original prism libraries   
!
!     
!     Contains   
!     ----------
!     
!      MODULE PRISM_CONSTANTS 
!
!      MODULE PRISM_DUMMY 
!         INTERFACES for prism routines called in ifs/prism/couplo4* 
!
!      SUBROUTINES
!          prism_dummy_* referenced in PRISM_DUMMY    
!
!      Method
!      -------            
!      Provide interfaces with original prism rouinte names and refers to 
!        dummy subroutines named  prism_dummy_*  if original name is prism_*  
!
!      To run coupled runs with OASIS4 USE PRISM in ifs/prism/couplo4* instead !      of USE PRISM_DUMMY    
!          
!     Reference:
!     ---------
!       S. Valcke, 2006: OASIS4 User Guide  
!       PRISM Support Initiative Report No 3,
!       CERFACS, Toulouse, France, 64 pp.
!
!     Author:
!     -------
!       Johannes Flemming 
!
!     Modifications.
!     --------------
!      F. Vana  05-Mar-2015  Support for single precision




!**************************************************************************

MODULE PRISM_DUMMY_MOD

USE PARKIND1 , ONLY : JPRD, JPIM
USE YOMPRISM 
USE MPL_MODULE, ONLY : MPL_ABORT
IMPLICIT NONE
SAVE

INTERFACE PRISM_SET_POINTS
MODULE PROCEDURE PRISM_DUMMY_SET_POINTS_GRIDLESS,PRISM_DUMMY_SET_POINTS_3D_DBLE
END INTERFACE

INTERFACE PRISM_PUT
MODULE PROCEDURE PRISM_DUMMY_PUT
END INTERFACE 


INTERFACE PRISM_GET
MODULE PROCEDURE PRISM_DUMMY_GET
END INTERFACE 


INTERFACE PRISM_INITIALIZED
MODULE PROCEDURE PRISM_DUMMY_INITIALIZED
END INTERFACE

 
INTERFACE  PRISM_INIT 
MODULE PROCEDURE PRISM_DUMMY_INIT
END INTERFACE

INTERFACE PRISM_INIT_COMP 
MODULE PROCEDURE PRISM_DUMMY_INIT_COMP
END INTERFACE

INTERFACE   PRISM_ENDDEF
MODULE PROCEDURE PRISM_DUMMY_ENDDEF
END INTERFACE

INTERFACE PRISM_TERMINATE
MODULE PROCEDURE PRISM_DUMMY_TERMINATE
END INTERFACE

INTERFACE PRISM_GET_LOCALCOMM
MODULE PROCEDURE PRISM_DUMMY_GET_LOCALCOMM
END INTERFACE

INTERFACE  PRISM_DEF_VAR
MODULE PROCEDURE PRISM_DUMMY_DEF_VAR
END INTERFACE

INTERFACE PRISM_DEF_GRID
MODULE PROCEDURE PRISM_DUMMY_DEF_GRID
END INTERFACE


INTERFACE PRISM_DEF_PARTITION
MODULE PROCEDURE PRISM_DUMMY_DEF_PARTITION
END INTERFACE


INTERFACE PRISM_REDUCEDGRID_MAP
MODULE PROCEDURE PRISM_DUMMY_REDUCEDGRID_MAP
END INTERFACE

INTERFACE PRISM_SET_CORNERS
MODULE PROCEDURE PRISM_DUMMY_SET_CORNERS
END INTERFACE


INTERFACE PRISM_CALC_NEWDATE
MODULE PROCEDURE PRISM_DUMMY_CALC_NEWDATE
END INTERFACE




!**************************************************************************
CONTAINS

SUBROUTINE PRISM_DUMMY_SET_POINTS_GRIDLESS (METHOD_ID, POINT_NAME, GRID_ID, &
                                 &          NEW_POINTS, IERROR)
       CHARACTER (LEN=*), INTENT(IN)                   :: POINT_NAME
       INTEGER(KIND=JPIM),           INTENT(IN)                   :: GRID_ID
       LOGICAL,           INTENT(IN)                   :: NEW_POINTS
       INTEGER(KIND=JPIM),           INTENT(INOUT)                :: METHOD_ID
       INTEGER(KIND=JPIM),           INTENT(OUT)                  :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')


END SUBROUTINE PRISM_DUMMY_SET_POINTS_GRIDLESS

!**************************************************************************

SUBROUTINE  PRISM_DUMMY_SET_POINTS_3D_DBLE(METHOD_ID, POINT_NAME, GRID_ID, & 
                  &        POINTS_ACTUAL_SHAPE, POINTS_1ST_ARRAY, POINTS_2ND_ARRAY, & 
                  &        POINTS_3RD_ARRAY, NEW_POINTS, IERROR)
       CHARACTER (LEN=*), INTENT(IN)                   :: POINT_NAME
       INTEGER(KIND=JPIM),           INTENT(IN)                   :: GRID_ID
       INTEGER(KIND=JPIM),           INTENT(INOUT)                :: METHOD_ID
       REAL(KIND=JPRD),  INTENT(IN), DIMENSION (:)    :: POINTS_1ST_ARRAY, POINTS_2ND_ARRAY, POINTS_3RD_ARRAY
       INTEGER(KIND=JPIM),           INTENT(IN)                   :: POINTS_ACTUAL_SHAPE (2, *)
       LOGICAL,           INTENT(IN)                   :: NEW_POINTS
       INTEGER(KIND=JPIM),           INTENT(OUT)                  :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_SET_POINTS_3D_DBLE

!**************************************************************************

SUBROUTINE PRISM_DUMMY_INIT ( APPL_NAME, IERROR )
       CHARACTER (LEN=*), INTENT(IN) :: APPL_NAME
       INTEGER(KIND=JPIM), INTENT (OUT)         :: IERROR

        STOP 'ABORT PRISM_DUMMY'

END SUBROUTINE PRISM_DUMMY_INIT

!**************************************************************************

SUBROUTINE PRISM_DUMMY_INIT_COMP ( COMP_ID, COMP_NAME, IERROR )  
       CHARACTER (LEN=*), INTENT(IN) :: COMP_NAME
       INTEGER(KIND=JPIM), INTENT (OUT) :: COMP_ID, IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_INIT_COMP

!**************************************************************************

SUBROUTINE PRISM_DUMMY_ENDDEF ( IERROR )
       INTEGER(KIND=JPIM), INTENT (OUT) :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_ENDDEF

!**************************************************************************

SUBROUTINE PRISM_DUMMY_TERMINATE ( IERROR )
      INTEGER(KIND=JPIM), INTENT (OUT) :: IERROR

      CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_TERMINATE

!**************************************************************************

SUBROUTINE PRISM_DUMMY_GET_LOCALCOMM ( COMP_ID, LOCAL_COMM, IERROR )
       INTEGER(KIND=JPIM), INTENT (IN)  :: COMP_ID
       INTEGER(KIND=JPIM), INTENT (OUT) :: LOCAL_COMM, IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_GET_LOCALCOMM

!**************************************************************************

SUBROUTINE PRISM_DUMMY_DEF_VAR (VAR_ID, NAME, GRID_ID, METHOD_ID, MASK_ID, &
                     &      VAR_NODIMS, VAR_ACTUAL_SHAPE, VAR_TYPE, IERROR )
       INTEGER(KIND=JPIM), INTENT(OUT)         :: VAR_ID
       CHARACTER (LEN=*), INTENT(IN) :: NAME
       INTEGER(KIND=JPIM), INTENT(IN)          :: GRID_ID
       INTEGER(KIND=JPIM), INTENT(IN)          :: METHOD_ID
       INTEGER(KIND=JPIM), INTENT(IN)          :: MASK_ID
       INTEGER(KIND=JPIM), INTENT(IN)          :: VAR_NODIMS(2)
       INTEGER(KIND=JPIM), INTENT(IN)          :: VAR_ACTUAL_SHAPE(1:2, 1:VAR_NODIMS(1)+VAR_NODIMS(2))
       INTEGER(KIND=JPIM), INTENT(IN)          :: VAR_TYPE
       INTEGER(KIND=JPIM), INTENT(OUT)         :: IERROR 

       CALL MPL_ABORT('ABORT PRISM_DUMMY')
     END SUBROUTINE PRISM_DUMMY_DEF_VAR

!**************************************************************************  

SUBROUTINE PRISM_DUMMY_DEF_GRID ( GRID_ID, GRID_NAME, COMP_ID, GRID_VALID_SHAPE, &
                              &     GRID_TYPE, IERROR)
       CHARACTER(LEN=*), INTENT(IN) :: GRID_NAME
       INTEGER(KIND=JPIM), INTENT(IN)          :: COMP_ID
       INTEGER(KIND=JPIM), INTENT(IN)          :: GRID_VALID_SHAPE (2, *)
       INTEGER(KIND=JPIM), INTENT(IN)          :: GRID_TYPE
       INTEGER(KIND=JPIM), INTENT(OUT)         :: GRID_ID, IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_DEF_GRID

!**************************************************************************

SUBROUTINE PRISM_DUMMY_DEF_PARTITION ( GRID_ID, NBR_BLOCKS, &
                      &           PARTITION_ARRAY, EXTENT_ARRAY, IERROR )
       INTEGER(KIND=JPIM), INTENT (IN)                :: GRID_ID
       INTEGER(KIND=JPIM), INTENT (IN)                :: NBR_BLOCKS
       INTEGER(KIND=JPIM), INTENT (IN)                :: PARTITION_ARRAY(1:NBR_BLOCKS,*)
       INTEGER(KIND=JPIM), INTENT (IN)                :: EXTENT_ARRAY(1:NBR_BLOCKS,*)
       INTEGER(KIND=JPIM), INTENT (OUT)               :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_DEF_PARTITION

!**************************************************************************

SUBROUTINE PRISM_DUMMY_REDUCEDGRID_MAP( GRID_ID,  NBR_LATITUDES, &
             &    NBR_POINTS_PER_LAT, IERROR )
       INTEGER(KIND=JPIM), INTENT (IN)                :: GRID_ID
       INTEGER(KIND=JPIM), INTENT (IN)                :: NBR_LATITUDES
       INTEGER(KIND=JPIM), INTENT (IN)                :: NBR_POINTS_PER_LAT(NBR_LATITUDES)
       INTEGER(KIND=JPIM), INTENT (OUT)               :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')

END SUBROUTINE PRISM_DUMMY_REDUCEDGRID_MAP

!**************************************************************************

SUBROUTINE PRISM_DUMMY_SET_CORNERS(GRID_ID, NBR_CORNERS, CORNERS_ACTUAL_SHAPE, &
        &    CORNERS_1ST_ARRAY, CORNERS_2ND_ARRAY, CORNERS_3RD_ARRAY, IERROR)
       INTEGER(KIND=JPIM), INTENT (IN)                :: GRID_ID
       INTEGER(KIND=JPIM), INTENT (IN)                :: NBR_CORNERS
       REAL(KIND=JPRD), INTENT (IN)       :: CORNERS_1ST_ARRAY (:,:)
       REAL(KIND=JPRD), INTENT (IN)       :: CORNERS_2ND_ARRAY (:,:)
       REAL(KIND=JPRD), INTENT (IN)       :: CORNERS_3RD_ARRAY (:,:)
       INTEGER(KIND=JPIM), INTENT (IN)                :: CORNERS_ACTUAL_SHAPE (1:2, *)
       INTEGER(KIND=JPIM), INTENT(OUT)                :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')
END SUBROUTINE PRISM_DUMMY_SET_CORNERS

!**************************************************************************

SUBROUTINE PRISM_DUMMY_CALC_NEWDATE ( DATE, DATE_INCR, IERROR )
       REAL(KIND=JPRD), INTENT(IN)                     :: DATE_INCR
       TYPE (PRISM_TIME_STRUCT), INTENT(INOUT) :: DATE
       INTEGER(KIND=JPIM), INTENT(OUT)                    :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')
END SUBROUTINE PRISM_DUMMY_CALC_NEWDATE

!**************************************************************************

SUBROUTINE PRISM_DUMMY_PUT ( FIELD_ID, DATE, DATE_BOUNDS, DATA_ARRAY, INFO, IERROR )
      INTEGER(KIND=JPIM), INTENT (IN)                 :: FIELD_ID
      TYPE(PRISM_TIME_STRUCT), INTENT (IN) :: DATE
      TYPE(PRISM_TIME_STRUCT), INTENT (IN) :: DATE_BOUNDS(2)
      DOUBLE PRECISION, INTENT (IN)        :: DATA_ARRAY(:,:)
      INTEGER(KIND=JPIM), INTENT (OUT)               :: INFO
      INTEGER(KIND=JPIM), INTENT (OUT)               :: IERROR

      CALL MPL_ABORT('ABORT PRISM_DUMMY')

END   SUBROUTINE PRISM_DUMMY_PUT

!**************************************************************************

SUBROUTINE PRISM_DUMMY_GET ( FIELD_ID, DATE, DATE_BOUNDS, DATA_ARRAY, INFO, IERROR )
      INTEGER(KIND=JPIM), INTENT (IN)                 :: FIELD_ID
      TYPE(PRISM_TIME_STRUCT), INTENT (IN) :: DATE
      TYPE(PRISM_TIME_STRUCT), INTENT (IN) :: DATE_BOUNDS(2)
      DOUBLE PRECISION, INTENT (OUT)        :: DATA_ARRAY(:,:)
      INTEGER(KIND=JPIM), INTENT (OUT)               :: INFO
      INTEGER(KIND=JPIM), INTENT (OUT)               :: IERROR

       CALL MPL_ABORT('ABORT PRISM_DUMMY')
END   SUBROUTINE PRISM_DUMMY_GET

!**************************************************************************

SUBROUTINE PRISM_DUMMY_INITIALIZED ( FLAG, IERROR )
       LOGICAL, INTENT (OUT) :: FLAG
       INTEGER(KIND=JPIM), INTENT (OUT) :: IERROR
      
        CALL MPL_ABORT('ABORT PRISM_DUMMY')
END SUBROUTINE PRISM_DUMMY_INITIALIZED
END MODULE PRISM_DUMMY_MOD
