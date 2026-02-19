# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### mpi_serial library

ecbuild_add_library( TARGET mpi_serial
  TYPE STATIC
  SOURCES
            mpi_abort.F
            mpi_allgatherv.F
            mpi_allreduce.F
            mpi_alltoallv.F
            mpi_barrier.F
            mpi_bcast.F
            mpi_bsend.F
            mpi_buffer_attach.F
            mpi_buffer_detach.F
            mpi_cart_coords.F
            mpi_cart_create.F
            mpi_cart_rank.F
            mpi_cart_sub.F
            mpi_comm_create.F
            mpi_comm_free.F
            mpi_comm_group.F
            mpi_comm_rank.F
            mpi_comm_size.F
            mpi_comm_split.F
            mpi_end.F
            mpi_error_string.F
            mpi_file_close.F
            mpi_file_iread_shared.F
            mpi_file_iwrite_shared.F
            mpi_file_open.F
            mpi_file_read_ordered.F
            mpi_file_read_ordered_begin.F
            mpi_file_read_ordered_end.F
            mpi_file_read_shared.F
            mpi_file_write_ordered.F
            mpi_file_write_ordered_begin.F
            mpi_file_write_ordered_end.F
            mpi_file_write_shared.F
            mpi_finalize.F
            mpi_finalized.F
            mpi_gather.F
            mpi_gatherv.F
            mpi_get_count.F
            mpi_group_free.F
            mpi_group_incl.F
            mpi_iallgatherv.F
            mpi_ialltoallv.F
            mpi_ibcast.F
            mpi_ibsend.F
            mpi_igather.F
            mpi_igatherv.F
            mpi_init.F
            mpi_init_thread.F
            mpi_initialized.F
            mpi_iprobe.F
            mpi_irecv.F
            mpi_iscatter.F
            mpi_iscatterv.F
            mpi_isend.F
            mpi_probe.F
            mpi_recv.F
            mpi_scatterv.F
            mpi_send.F
            mpi_ssend.F
            mpi_testsome.F
            mpi_type_size.F
            mpi_wait.F
            mpi_waitall.F
            mpi_waitany.F
            mpi_waitsome.F
            vpp_abort.F
            vpp_barrier.F
    PUBLIC_INCLUDES 
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<INSTALL_INTERFACE:include/mpi_serial>
)

## Install and Export ##

install(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/mpif.h
    DESTINATION
        include/mpi_serial
)
