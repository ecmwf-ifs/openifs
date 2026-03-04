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

MODULE TOVSCV_MOD
!-----------------------------------------------------------------------
! Module containing the Type defenition and standard procedures for the
! so-called TOVS control variable.
! Written for OOPS but also used in pre-OOPS IFS
! In OOPS this is used for the representation of the state, increment and
! control-variable

!   Author.
!   -------
!   Mats Hamrud 01-Dec-2017 - From old code
!   Marcin Chrust  May-2020 - Parallel I/O
!-----------------------------------------------------------------------

USE PARKIND1,     ONLY : JPIM, JPRB, JPRD
USE YOMANCS,      ONLY : RMDI
USE TOVSCV_BASE_MOD, ONLY : TOVSCV_BASE, CREATE_TOVSCV_BASE, DEALLOCATE_TOVSCV_BASE, CHECK_DIMS
USE YOMVAR,       ONLY : LTOVSCV, LCLDSINK, LINC_TOVSCV
USE YOMLUN,       ONLY : NULOUT, NULERR
USE YOMMP0,       ONLY : MYPROC,NPROC
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMSATS,      ONLY : NVATOVINDX
USE MPL_MODULE,   ONLY : MPL_ALLREDUCE, MPL_ALLGATHERV, MPL_BROADCAST
USE IOSTREAM_MIX, ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,IO_GET, &
 &                       CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST,CLOSE_IOREQUEST
USE GENERIC_CTLVEC_MOD,  ONLY : GENERIC_CTLVEC, GATHER_GENERIC_CTLVEC, SCATTER_GENERIC_CTLVEC, &
 &                              MULTIGATHER_GENERIC_CTLVEC, MULTISCATTER_GENERIC_CTLVEC, &
 &                              GET_SIZE_GENERIC_CTLVEC, MULTI_IO
USE ORDER_INDEPENDENT_SUMMATION_MOD , ONLY : ORDER_INDEP_GLOBAL_SUM
USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA

IMPLICIT NONE
PRIVATE
  TYPE,EXTENDS(TOVSCV_BASE),PUBLIC :: TOVSCV
  REAL(KIND=JPRB),ALLOCATABLE :: TOVSCVX(:,:)
  CONTAINS
    PROCEDURE :: CREATE_TOVSCV
  END TYPE TOVSCV

  CONTAINS
  ! ------------------------------------------------------------------------------
  SUBROUTINE CREATE_TOVSCV(SELF)

  CLASS(TOVSCV), INTENT(INOUT)  :: SELF

  call abort1("Called dummy routine create_tovscv - exit")

  END SUBROUTINE CREATE_TOVSCV

END MODULE TOVSCV_MOD
