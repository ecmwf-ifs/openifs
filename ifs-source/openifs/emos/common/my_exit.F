! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

       SUBROUTINE MY_EXIT(K)
C
       IMPLICIT NONE
       INTEGER K
#ifdef __uxp__
C
C      VPP700 CALL TO SET EXIT NUMBER
       CALL SETRCD(K)
#endif
C
       CALL EXIT(K)
C
       RETURN
       END
