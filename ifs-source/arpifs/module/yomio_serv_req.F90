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

MODULE YOMIO_SERV_REQ

!**** *YOMIO_SERV_REQ*  - 

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO-FRANCE*
!      Original : 12-11-2012

!     Modifications :
!     P. Marguinaud : 10-10-2013 : Make it work

USE PARKIND1, ONLY : JPIM, JPRB

USE YOMIO_SERV_HDR, ONLY : IO_SERV_HDRG
USE IOFLDDESC_MOD,  ONLY : IOFLDDESC
USE IOCPTDESC_MOD,  ONLY : IOCPTDESC

IMPLICIT NONE

INTEGER (KIND=JPIM), PARAMETER :: NREQIDSIZE            = 4_JPIM,     &
                                & NIO_SERV_REQ_TYPE_IDX = 1_JPIM,     &
                                & NIO_SERV_REQ_STEP_IDX = 2_JPIM,     &
                                & NIO_SERV_REQ_TAG__IDX = 3_JPIM,     &
                                & NIO_SERV_REQ_SET__IDX = 4_JPIM
       

TYPE IO_SERV_REQ
  TYPE (IO_SERV_HDRG) :: YHDR                                                    ! Header
  INTEGER (KIND=JPIM) :: ID (NREQIDSIZE) = (/ 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM /)  ! Request id
  LOGICAL             :: LARRAY2D = .FALSE.
  REAL (KIND=JPRB),    POINTER :: ZFLBU1 (:)    => NULL ()     ! Field data
  REAL (KIND=JPRB),    POINTER :: ZFLBU2 (:, :) => NULL ()     ! Field data as 2D array (field size : nfields)
  REAL (KIND=JPRB),    POINTER :: ZFLBU3 (:, :, :) => NULL ()  ! Field data as 2D arrays ( field size x : field size y : nfields)
  TYPE (IOFLDDESC),    POINTER :: YFLDSC (:)    => NULL ()     ! Field descriptors
  TYPE (IOCPTDESC),    POINTER :: YCPDSC (:)    => NULL ()     ! Compressed field descriptors
  INTEGER (KIND=JPIM), POINTER :: IFLOFF (:)    => NULL ()     ! Field offsets
  INTEGER (KIND=JPIM), POINTER :: IFLDOM (:)    => NULL ()     ! Field domains
  INTEGER (KIND=JPIM) :: NSIZER = 0_JPIM                       ! Number of points received
  INTEGER (KIND=JPIM) :: NSIZED = 0_JPIM                       ! Number of 8-byte words received (including header data)
  LOGICAL             :: LHDROK = .FALSE.                      ! Field header set 
  LOGICAL             :: LLDONE = .FALSE.                      ! Request has completed
  REAL (KIND=JPRB)    :: ZTIME1 = 0._JPRB                      ! Time at which the request was created
  REAL (KIND=JPRB)    :: ZTIME2 = 0._JPRB                      ! Time at which the request has completed
  INTEGER (KIND=JPIM) :: IREQ_CREAT = 0_JPIM                   ! Rank of the request (creation order)
  INTEGER (KIND=JPIM) :: IREQ_COMPL = 0_JPIM                   ! Rank of the request (completion order)
END TYPE IO_SERV_REQ

TYPE IO_SERV_REQ_PTR
  TYPE (IO_SERV_REQ), POINTER :: P => NULL ()
END TYPE IO_SERV_REQ_PTR

SAVE

END MODULE YOMIO_SERV_REQ

