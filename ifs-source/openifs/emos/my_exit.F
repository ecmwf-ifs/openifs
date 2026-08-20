! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

       SUBROUTINE MY_EXIT(K)
C
       USE PARKIND1  ,ONLY : JPIM
       USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
       IMPLICIT NONE
       REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
       INTEGER, INTENT(IN) :: K
       INTEGER ISIGABRT
       
       IF (LHOOK) CALL DR_HOOK('MY_EXIT',0,ZHOOK_HANDLE)
       IF (K /= 0) THEN
         write(0,*)'MY_EXIT-code = ',K
         ISIGABRT = 6
         IF (LHOOK) CALL C_DRHOOK_RAISE(ISIGABRT)
       ENDIF
#ifdef __uxp__
C
C      VPP700 CALL TO SET EXIT NUMBER
       CALL SETRCD(K)
#endif
C
       CALL EXIT(K)
C
       IF (LHOOK) CALL DR_HOOK('MY_EXIT',1,ZHOOK_HANDLE)
       RETURN
       END SUBROUTINE MY_EXIT

       SUBROUTINE MY_EXIT2(K,CH)
C
       USE PARKIND1  ,ONLY : JPIM
       USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
       IMPLICIT NONE
       REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
       INTEGER, INTENT(IN) :: K
       CHARACTER(LEN=*), INTENT(IN) :: CH
       INTEGER ISIGABRT
       
       IF (LHOOK) CALL DR_HOOK('MY_EXIT2',0,ZHOOK_HANDLE)
       IF (K /= 0) THEN
         write(0,*)'MY_EXIT2-code = ',K
         write(0,'(1x,a)') CH
         ISIGABRT = 6
         IF (LHOOK) CALL C_DRHOOK_RAISE(ISIGABRT)
       ENDIF
#ifdef __uxp__
C
C      VPP700 CALL TO SET EXIT NUMBER
       CALL SETRCD(K)
#endif
C
       CALL EXIT(K)
C
       IF (LHOOK) CALL DR_HOOK('MY_EXIT2',1,ZHOOK_HANDLE)
       RETURN
       END SUBROUTINE MY_EXIT2
