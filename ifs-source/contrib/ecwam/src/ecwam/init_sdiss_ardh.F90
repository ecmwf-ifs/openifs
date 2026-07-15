! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INIT_SDISS_ARDH
! ----------------------------------------------------------------------

!**** *INIT_SDISS_ARDH* - INITIALISATION FOR SDISS_ARD

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *INIT_SDISS_ARDH

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR, TH, DELTH
      USE YOWPCONS , ONLY : RAD     ,G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPHYS  , ONLY : ISDSDTH, ISB, NSDSNTH, INDICESSAT, SATWEIGHTS

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: JD, K, M, I_INT, J_INT, NANGD

      REAL(KIND=JWRB) :: TPIINV, TMP01
      REAL(KIND=JWRB) :: DELTH_TRUNC, DELTH_LOC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INIT_SDISS_ARDH',0,ZHOOK_HANDLE)

      TPIINV = 1.0_JWRB/ZPI

      NANGD=NANG/2

!     COMPUTE SATWEIGHTS

      NSDSNTH  = MIN(NINT(ISDSDTH*RAD/(DELTH)),NANGD-1)
      DELTH_TRUNC=(TH(1)+ISDSDTH*RAD)-(TH(1+NSDSNTH)-0.5_JWRB*DELTH)
      DELTH_TRUNC=MAX(0.0_JWRB, MIN(DELTH_TRUNC,DELTH))

      IF (ALLOCATED(INDICESSAT)) DEALLOCATE(INDICESSAT)
      ALLOCATE(INDICESSAT(NANG,NSDSNTH*2+1))
      IF (ALLOCATED(SATWEIGHTS)) DEALLOCATE(SATWEIGHTS)
      ALLOCATE(SATWEIGHTS(NANG,NSDSNTH*2+1))

      DO K=1,NANG
        DO I_INT=K-NSDSNTH, K+NSDSNTH
          J_INT=I_INT
          IF (I_INT < 1)  J_INT=I_INT+NANG
          IF (I_INT > NANG) J_INT=I_INT-NANG
          INDICESSAT(K,I_INT-(K-NSDSNTH)+1)=J_INT

          IF (I_INT == K-NSDSNTH .OR. I_INT == K+NSDSNTH) THEN
            DELTH_LOC=DELTH_TRUNC
          ELSE
            DELTH_LOC=DELTH
          ENDIF
          SATWEIGHTS(K,I_INT-(K-NSDSNTH)+1)=DELTH_LOC*COS(TH(K)-TH(J_INT))**ISB
        END DO
      END DO

      IF (LHOOK) CALL DR_HOOK('INIT_SDISS_ARDH',1,ZHOOK_HANDLE)

      END SUBROUTINE INIT_SDISS_ARDH
