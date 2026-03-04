! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_MODULE
!
!    Message Passing Library (MPL)
!
!    Interface between parallel applications and the
!    Message Passing Interface (MPI standard) provided by the computer vendors
!    This version requires only MPI release 1.
!
!    designed and developed by
!      Mats Hamrud and David Dent,      ECMWF,      September 2000
!
!    all routines which wish to call MPL routines must contain:
!    USE MPL_MODULE

USE MPL_DATA_MODULE
USE MPL_INIT_MOD
USE MPL_BUFFER_METHOD_MOD
USE MPL_SEND_MOD
USE MPL_RECV_MOD
USE MPL_WAIT_MOD
USE MPL_BARRIER_MOD
USE MPL_BROADCAST_MOD
USE MPL_PROBE_MOD
USE MPL_END_MOD
USE MPL_MESSAGE_MOD
USE MPL_ABORT_MOD
USE MPL_COMM_CREATE_MOD
USE MPL_COMM_FREE_MOD
USE MPL_COMM_SPLIT_MOD
USE MPL_SETDFLT_COMM_MOD
USE MPL_ALLGATHER_MOD
USE MPL_MYRANK_MOD
USE MPL_NPROC_MOD
USE MPL_IOINIT_MOD
USE MPL_OPEN_MOD
USE MPL_CLOSE_MOD
USE MPL_READ_MOD
USE MPL_WRITE_MOD
USE MPL_ALLREDUCE_MOD
USE MPL_GATHERV_MOD
USE MPL_MYGATHERV_MOD
USE MPL_ALLGATHERV_MOD
USE MPL_ALLTOALLV_MOD
USE MPL_SCATTERV_MOD
USE MPL_GROUPS
USE MPL_ARG_MOD
USE MPL_LOCOMM_CREATE_MOD
USE MPL_TOUR_TABLE_MOD
USE MPL_TESTSOME_MOD
USE MPL_WAITANY_MOD
USE MPL_BYTES_MOD

END MODULE MPL_MODULE
