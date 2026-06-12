MODULE KPP_TRIDRHS_MOD
CONTAINS
SUBROUTINE KPP_TRIDRHS &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
  &   LDKPPCAL ,PHO      ,PHO_INV  ,PTRI0    ,PTRI1    ,&
  &   PYO      ,PNTFLUX  ,PDIFF    ,PGHAT    ,PSTURFLUX,&
  &   PGHATFLUX,PRHS     ,PDTO     ,PADV     ,PRHO     ,&
  &   PCP      ,KSCLR    ,YDOCEAN_ML )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine computes right hand side of tridiagonal matrix 
!   for scalar fields:

! Interface :
! ---------
!   Call *KPP_TRIDRHS* from *KPP_OCNINT*

! Method :
! ------
!   implicit integration scheme (Backward Euler Method)

!     compute right hand side of tridiagonal matrix for scalar fields:
!     r.h.s. =  yo (old field) + flux-divergence of ghat
!               + flux-divergence of non-turbulant fluxes
!     note: surface layer needs +dto/h(1) * surfaceflux
!           bottom  ..... ..... +dto/h(jpnz)*diff(jpnz)/dzb(jpnz)*yo(jpnz+1)

! Externals :
! ---------

! Reference :
! ---------
! Large, W. G., J. C. McWilliams, S. C. Doney (1994), Rev. Geophys.

! Modifications :
! -------------
!     06-Jun-1994  Bill Large
!            2002  Steve Woolnough, Reading Univ.
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1,     ONLY : JPIM, JPRB
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML, NVELO, NSCLRO
USE ABORT_SURF_MOD
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVOP1
INTEGER(KIND=JPIM),INTENT(IN) :: KSCLR               ! 1: Temp. 2: Sal. 

REAL(KIND=JPRB),INTENT(IN)  :: PHO(KLON,KLEVOP1)     ! layer thickness
REAL(KIND=JPRB),INTENT(IN)  :: PHO_INV(KLON,KLEVOP1) ! 1/PHO
REAL(KIND=JPRB),INTENT(IN)  :: PTRI0(KLON,0:KLEVO)   ! array for diff. eq.
REAL(KIND=JPRB),INTENT(IN)  :: PTRI1(KLON,0:KLEVO)   ! array for diff. eq.
REAL(KIND=JPRB),INTENT(IN)  :: PYO(KIDIA:KFDIA,KLEVOP1)     ! old profile
REAL(KIND=JPRB),INTENT(IN)  :: PNTFLUX(KIDIA:KFDIA,0:KLEVO) 
                                                     ! non-turbulent flux 
                                                     ! = wxnt(0:jpnz,1:2)
REAL(KIND=JPRB),INTENT(IN)  :: PDIFF(KLON,0:KLEVO)   ! diff. profile on intfc
REAL(KIND=JPRB),INTENT(IN)  :: PGHAT(KLON,KLEVO)     ! ghat turbulent flux   
REAL(KIND=JPRB),INTENT(IN)  :: PSTURFLUX(KIDIA:KFDIA)! sfc turb (kinematic) 
                                                     ! flux = wx(0,n)
REAL(KIND=JPRB),INTENT(IN)  :: PGHATFLUX(KIDIA:KFDIA)! surface flux for ghat 
                                                     ! (includes solar flux)
REAL(KIND=JPRB),INTENT(IN)  :: PDTO                  ! time interval
REAL(KIND=JPRB),INTENT(IN)  :: PADV(KIDIA:KFDIA,KLEVO,NSCLRO) ! advection term
REAL(KIND=JPRB),INTENT(IN)  :: PRHO(KLON,0:KLEVOP1)
REAL(KIND=JPRB),INTENT(IN)  :: PCP(KLON,0:KLEVOP1)

REAL(KIND=JPRB),INTENT(OUT) :: PRHS(KIDIA:KFDIA,KLEVO) ! right hand side

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)

TYPE(TOCEAN_ML), INTENT(IN) :: YDOCEAN_ML

INTEGER(KIND=JPIM) :: JZ                       
INTEGER(KIND=JPIM) :: JL                       

REAL(KIND=JPRB)    :: ZDIVFLX                     
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_TRIDRHS_MOD:KPP_TRIDRHS',0,ZHOOK_HANDLE)
ASSOCIATE(NPD_KPP=>YDOCEAN_ML%NPD_KPP)

ZDIVFLX = 1.0_JPRB / REAL(NPD_KPP)

! 1. the surface layer (dto/h(1)=tri(0,1,ind)
!    ----------------------------------------

SELECT CASE ( KSCLR ) 
  CASE(1) ! Temperature
    DO JL = KIDIA, KFDIA
      IF ( LDKPPCAL(JL) ) THEN
        PRHS(JL,1) = PYO(JL,1) + PDTO * PHO_INV(JL,1)                    &
                   &      * ( PGHATFLUX(JL) * PDIFF(JL,1) * PGHAT(JL,1)  &
                   &          - PSTURFLUX(JL) * ZDIVFLX                  &
                   &          + PNTFLUX(JL,1) - PNTFLUX(JL,0) )          &
                   &  + PDTO * PADV(JL,1,KSCLR)                          &
                   &         / ( PRHO(JL,1)*PCP(JL,1) )
      ENDIF
    ENDDO
  CASE(2:NSCLRO) !Salinity
    DO JL = KIDIA, KFDIA
      IF ( LDKPPCAL(JL) ) THEN
        PRHS(JL,1) = PYO(JL,1) + PDTO * PHO_INV(JL,1)                    &
                   &      * ( PGHATFLUX(JL) * PDIFF(JL,1) * PGHAT(JL,1)  &
                   &          - PSTURFLUX(JL) * ZDIVFLX                  &
                   &          + PNTFLUX(JL,1) - PNTFLUX(JL,0) )          &
                   &  + PDTO * PADV(JL,1,KSCLR) 
      ENDIF
    ENDDO
  CASE DEFAULT
    CALL ABORT_SURF('KPP_TRIDRHS: KSCLR IS INVALID.')
END SELECT

! 2. inside the domain to npd
!    ------------------------

SELECT CASE ( KSCLR ) 
  CASE(1) ! Temperature
    IF(NPD_KPP >= 2) THEN
      DO JL = KIDIA, KFDIA
        IF ( LDKPPCAL(JL) ) THEN
          DO JZ = 2, NPD_KPP
            PRHS(JL,JZ) = PYO(JL,JZ) + PDTO * PHO_INV(JL,JZ)                 &
                        &   * ( PGHATFLUX(JL) * PDIFF(JL,JZ) * PGHAT(JL,JZ)  &
                        &      - PGHATFLUX(JL)*PDIFF(JL,JZ-1)*PGHAT(JL,JZ-1) &
                        &      - PSTURFLUX(JL) * ZDIVFLX                     &
                        &      + PNTFLUX(JL,JZ)- PNTFLUX(JL,JZ-1) )          &
                        &  + PDTO * PADV(JL,JZ,KSCLR)                        &
                        &         / ( PRHO(JL,JZ)*PCP(JL,JZ) )
          ENDDO
        ENDIF
      ENDDO 
    ENDIF
  CASE(2:NSCLRO) !Salinity and other scalers
    IF(NPD_KPP >= 2) THEN
      DO JL = KIDIA, KFDIA
        IF ( LDKPPCAL(JL) ) THEN
          DO JZ = 2, NPD_KPP
            PRHS(JL,JZ) = PYO(JL,JZ) + PDTO * PHO_INV(JL,JZ)               &
                      &   * ( PGHATFLUX(JL) * PDIFF(JL,JZ) * PGHAT(JL,JZ)  &
                      &      - PGHATFLUX(JL)*PDIFF(JL,JZ-1)*PGHAT(JL,JZ-1) &
                      &      - PSTURFLUX(JL) * ZDIVFLX                     &
                      &      + PNTFLUX(JL,JZ)- PNTFLUX(JL,JZ-1) )          &
                      &  + PDTO * PADV(JL,JZ,KSCLR) 
          ENDDO
        ENDIF
      ENDDO 
    ENDIF
END SELECT

! 3. inside the rest of the domain
!    -----------------------------

SELECT CASE ( KSCLR ) 
  CASE(1) ! Temperature
    DO JL = KIDIA, KFDIA
      IF ( LDKPPCAL(JL) ) THEN
        DO JZ = NPD_KPP+1, KLEVO-1
          PRHS(JL,JZ) = PYO(JL,JZ) + PDTO * PHO_INV(JL,JZ)               &
                      &  * ( PGHATFLUX(JL)*(PDIFF(JL,JZ)*PGHAT(JL,JZ)    & 
                      &      - PDIFF(JL,JZ-1) * PGHAT(JL,JZ-1) )         &
                      &      + PNTFLUX(JL,JZ) - PNTFLUX(JL,JZ-1) )       &
                      &  + PDTO * PADV(JL,JZ,KSCLR)                      &
                      &         / ( PRHO(JL,JZ)*PCP(JL,JZ) )
        ENDDO       
      ENDIF
    ENDDO 
  CASE(2:NSCLRO) !Salinity and other scalers
    DO JL = KIDIA, KFDIA
      IF ( LDKPPCAL(JL) ) THEN
        DO JZ = NPD_KPP+1, KLEVO-1
          PRHS(JL,JZ) = PYO(JL,JZ) + PDTO * PHO_INV(JL,JZ)               &
                      &  * ( PGHATFLUX(JL)*(PDIFF(JL,JZ)*PGHAT(JL,JZ)    & 
                      &      - PDIFF(JL,JZ-1) * PGHAT(JL,JZ-1) )         &
                      &      + PNTFLUX(JL,JZ) - PNTFLUX(JL,JZ-1) )       &
                      &  + PDTO * PADV(JL,JZ,KSCLR) 
        ENDDO       
      ENDIF
    ENDDO 
END SELECT

! 4. in the bottom layer     
!    -------------------

SELECT CASE ( KSCLR ) 
  CASE(1) ! Temperature
    IF(KLEVO > 1) THEN ! not for slab ocean
      DO JL = KIDIA, KFDIA
        IF ( LDKPPCAL(JL) ) THEN
          PRHS(JL,KLEVO)= PYO(JL,KLEVO) + PDTO * PHO_INV(JL,KLEVO)           &
                        & * ( PGHATFLUX(JL)*(PDIFF(JL,KLEVO)*PGHAT(JL,KLEVO) &
                        &            - PDIFF(JL,KLEVO-1)*PGHAT(JL,KLEVO-1) ) &
                        &        + PNTFLUX(JL,KLEVO) - PNTFLUX(JL,KLEVO-1) ) &
                        &   + PYO(JL,KLEVO+1)*PTRI1(JL,KLEVO)*PDIFF(JL,KLEVO)&
                        & + PDTO * PADV(JL,KLEVO,KSCLR)                      &
                        &        / ( PRHO(JL,KLEVO)*PCP(JL,KLEVO) )
        ENDIF
      ENDDO 
    ENDIF
  CASE(2:NSCLRO) !Salinity and other scalers
    IF(KLEVO > 1) THEN ! not for slab ocean
      DO JL = KIDIA, KFDIA
        IF ( LDKPPCAL(JL) ) THEN
          PRHS(JL,KLEVO)= PYO(JL,KLEVO) + PDTO * PHO_INV(JL,KLEVO)            &
                        &  * ( PGHATFLUX(JL)*(PDIFF(JL,KLEVO)*PGHAT(JL,KLEVO) &
                        &             - PDIFF(JL,KLEVO-1)*PGHAT(JL,KLEVO-1) ) &
                        &         + PNTFLUX(JL,KLEVO) - PNTFLUX(JL,KLEVO-1) ) &
                        &   + PYO(JL,KLEVO+1)*PTRI1(JL,KLEVO)*PDIFF(JL,KLEVO)&
                        & + PDTO * PADV(JL,KLEVO,KSCLR)
        ENDIF
      ENDDO 
    ENDIF
END SELECT

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_TRIDRHS_MOD:KPP_TRIDRHS',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_TRIDRHS
END MODULE KPP_TRIDRHS_MOD
