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


!     Initial determination of distributed versus shared memory configuration
!     Contains YOMMP0 and YOMGSTATS variables.
!     -----------------------------------------------------------------------
NAMELIST/NAMPAR0/ NPROC,NOUTPUT,&
                 &NPRGPNS,NPRGPEW,NPRTRW,NPRTRV,&
                 &NSPECRESMIN,&
                 &MP_TYPE,MBX_SIZE,&
                 &LMPOFF,LMPDIAG,&
                 &LOPT_SCALAR,LOPT_RS6K,NPRINTLEV,LSCMEC,&
                 &LSTATS,LSTATSCPU,LSYNCSTATS,LDETAILED_STATS,LXML_STATS,&
                 &LSTATS_MEM,NSTATS_MEM,LSTATS_ALLOC,LBARRIER_STATS,LBARRIER_STATS2,&
                 &NPRNT_STATS,NTRACE_STATS
!     ------------------------------------------------------------------

