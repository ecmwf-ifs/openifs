! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE INITDPTHFLDS(WVENVI, WVPRPT, WVPRPT_LAND)
! ----------------------------------------------------------------------

!**** *INITDPTHFLDS* - 

!*    PURPOSE.
!     --------

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH:
!     EMAXDPT

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY:

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWFRED  , ONLY : ZPIFR
      USE YOWPARAM , ONLY : NFRE
      USE YOWSHAL  , ONLY : NDEPTH,                                &
     &                      TFAK, TCGOND, TSIHKD, TFAC_ST, GAM_B_J

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(ENVIRONMENT), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), DIMENSION(NPROMA_WAM, NFRE, NCHNK), INTENT(OUT) :: WVPRPT
      TYPE(FREQUENCY), DIMENSION(NFRE), INTENT(OUT) :: WVPRPT_LAND


      INTEGER(KIND=JWIM) :: M, IJ, ICHNK, KIJS, KIJL

      REAL(KIND=JWRB) :: GAM
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',0,ZHOOK_HANDLE)

ASSOCIATE(DEPTH => WVENVI%DEPTH, &
 &        INDEP => WVENVI%INDEP, &
 &        EMAXDPT => WVENVI%EMAXDPT, &
 &        WAVNUM => WVPRPT%WAVNUM, &
 &        CINV => WVPRPT%CINV, &
 &        CGROUP => WVPRPT%CGROUP, &
 &        XK2CG => WVPRPT%XK2CG, &
 &        OMOSNH2KD => WVPRPT%OMOSNH2KD, &
 &        STOKFAC => WVPRPT%STOKFAC, &
 &        CIWA => WVPRPT%CIWA )


      CALL GSTATS(1502,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, KIJL, IJ, M, GAM)
      DO ICHNK = 1, NCHNK
        KIJS = 1
        KIJL = NPROMA_WAM

!       COMPUTE THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
        DO IJ = KIJS, KIJL
!         REDUCE GAMMA FOR SMALL DEPTH ( < 4m)
!         (might need to be revisted when grid is fine resolution)
          IF (DEPTH(IJ, ICHNK) < 4.0_JWRB) THEN
            GAM=GAM_B_J*DEPTH(IJ, ICHNK)/4.0_JWRB
          ELSE
            GAM=GAM_B_J
          ENDIF
          EMAXDPT(IJ, ICHNK) = 0.0625_JWRB*(GAM*DEPTH(IJ, ICHNK))**2
        ENDDO

        DO M= 1, NFRE
          DO IJ = KIJS, KIJL
            WAVNUM(IJ, M, ICHNK) = TFAK(INDEP(IJ, ICHNK), M)
            CINV(IJ, M, ICHNK)   = WAVNUM(IJ, M, ICHNK)/ZPIFR(M)
            CGROUP(IJ, M, ICHNK) = TCGOND(INDEP(IJ, ICHNK), M)
            XK2CG(IJ, M, ICHNK) = WAVNUM(IJ, M, ICHNK)**2 * CGROUP(IJ, M, ICHNK)
            OMOSNH2KD(IJ, M, ICHNK) = TSIHKD(INDEP(IJ, ICHNK), M)
            STOKFAC(IJ, M, ICHNK) = TFAC_ST(INDEP(IJ, ICHNK), M)

            CIWA(IJ,M, ICHNK) = 1.0_JWRB
          ENDDO
        ENDDO

      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1502,1)


      ! Fictitious values for land point (NSUP+1)
      DO M=1,NFRE
        WVPRPT_LAND(M)%WAVNUM = TFAK(NDEPTH,M)
        WVPRPT_LAND(M)%CINV   = WVPRPT_LAND(M)%WAVNUM/ZPIFR(M)
        WVPRPT_LAND(M)%CGROUP = TCGOND(NDEPTH,M)
        WVPRPT_LAND(M)%XK2CG = WVPRPT_LAND(M)%WAVNUM**2 * WVPRPT_LAND(M)%CGROUP
        WVPRPT_LAND(M)%OMOSNH2KD = TSIHKD(NDEPTH,M)
        WVPRPT_LAND(M)%STOKFAC = TFAC_ST(NDEPTH,M)

        WVPRPT_LAND(M)%CIWA = 1.0_JWRB
      ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE INITDPTHFLDS
