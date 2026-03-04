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

module obs_store_options_mod

   use yommp0, only   : myproc
   use parkind1, only : jpim

   implicit none

   type obs_store_options_t
      logical            :: l_read_from_obs_store = .false.
      logical            :: l_write_to_obs_store  = .false.
      character(len=256) :: odb_obs_store_path    = ""
      integer(kind=jpim) :: npools_copedb         = 0
   contains
      procedure, public  :: setup => setup_obs_store_options
   end type obs_store_options_t

   type(obs_store_options_t)  :: ydobs_store_options

   
contains

   subroutine setup_obs_store_options(this)
      class(obs_store_options_t), intent(inout)  :: this
      character(len=256)                         :: clenv

      ! Pick up settings from environment variables
      call get_environment_variable("L_READ_FROM_OBS_STORE",clenv)
      if(trim(clenv) == "true") this%l_read_from_obs_store = .true.

      call get_environment_variable("L_WRITE_TO_OBS_STORE",clenv)
      if(trim(clenv) == "true") this%l_write_to_obs_store = .true.
  
      call get_environment_variable("ODB_OBS_STORE_PATH",clenv)
      this%odb_obs_store_path = trim(clenv)

      call get_environment_variable("NPOOLS_COPEDB",clenv)
      read(clenv,'(i6)') this%npools_copedb

      ! TODO: for OOPS configurations, pick up from JSON

      if(myproc == 1)then
        write(0,*) "ODB Observation Store options : L_READ_FROM_OBS_STORE = ",this%l_read_from_obs_store
        write(0,*) "ODB Observation Store options : L_WRITE_TO_OBS_STORE  = ",this%l_write_to_obs_store
        write(0,*) "ODB Observation Store options : ODB_OBS_STORE_PATH    = ",trim(this%odb_obs_store_path)
        write(0,*) "ODB Observation Store options : NPOOLS_COPEDB         = ",this%npools_copedb
      endif

   end subroutine setup_obs_store_options

end module obs_store_options_mod
