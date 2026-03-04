! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module dbase_mod
    !     Purpose.
    !     --------
    !     Provide DBASE_MOD from IFSOBS or a dummy if building without it.
    
    !     Author.
    !     -------
    !     Balthasar Reuter (ECMWF)
    
    !     Modifications.
    !     --------------
    !        Original : 2022-03-03
    
    !-----------------------------------------------------------------
    use parkind1, only: jpim
    
    use ifs_dbase_view_mod, only : abstract_dbase_view
    !-----------------------------------------------------------------
    
    implicit none
    
    type, public :: dbase
    contains
      procedure, public :: select => ifs_dummy_dbase_select
      procedure, public :: put    => ifs_dummy_dbase_put
    end type
    
    !-----------------------------------------------------------------
    contains
    !-----------------------------------------------------------------
    
    function ifs_dummy_dbase_select(this,viewname,view,poolno,sql) result(rc)
       !! Runs an SQL query on the database, returning the results as a dbase_view object.
       !! Example:
       !!   rc = mydb%select("robody.sql",myview)
       implicit none
       class(dbase),intent(inout) :: this
       character(len=*), intent(in)                :: viewname   !! Name of the SQL query that you want to run
       class(abstract_dbase_view), intent(out)     :: view       !! Output data object containing the results from the query
       integer(kind=jpim), intent(in),optional     :: poolno     !! Pool/parition id on which you want to run the query
       character(len=*), intent(in), optional      :: sql        !! SQL query string
       integer(kind=jpim)                          :: rc         !! Return code supplied from lower level library
       call abort1("Called dummy ifsobs function DBASE%SELECT")
    end function
    
    !-----------------------------------------------------------------
    
    function ifs_dummy_dbase_put(this,view,poolno) result(rc)
       !! Inserts data from the supplied dbase_view object into the database in the rows selected by the dbase_view%view_name SQL query.
       !! Following the select(),get(),put() access pattern used in ODB.
       !! Example:
       !!   rc = mydb%select("robody.sql",myview)
       !!   call myview%fill_mdi()
       !!   rc = mydb%put(myview)
       implicit none
       class(dbase), intent(inout)                   :: this
       class(abstract_dbase_view), intent(inout)     :: view     !! Data to be inserted into the database (the view still remembers which SQL view it represents)
       integer(kind=jpim), intent(in),optional       :: poolno   !! Pool/partition id into which you want to put the data
       integer(kind=jpim)                            :: rc       !! Return code supplied from lower level library
       call abort1("Called dummy ifsobs function DBASE%PUT")
    end function
    
    !-----------------------------------------------------------------
    
    end module dbase_mod