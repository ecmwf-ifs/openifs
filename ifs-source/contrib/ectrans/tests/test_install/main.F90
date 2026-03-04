! (C) Copyright 2020- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program main

! From library fiat
  use yomhook      ! assert found
  use mpl_module   ! assert found
! From library parkind_(dp|sp)
  use parkind1, only: JPRB ! assert found

implicit none

! assert includes are found
#include "setup_trans0.h"
#include "trans_end.h"

write(0,*) "JPRB =",JPRB ! depending on link with parkind_sp or parkind_dp this will print 4 or 8

end program