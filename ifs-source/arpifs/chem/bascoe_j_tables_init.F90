! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_J_TABLES_INIT

!**   DESCRIPTION 
!     ----------
!
!   initialize J tables and auxiliary variables (defined in module BASCOE_J_MODULE)
!       for C-IFS-BASCOE stratospheric chemistry, using TUV
!
! **** WARNING: THIS IS FOR THE OFFLINE MODE !
! **** WARNING: because this routine is called once (indirectly from CHEM_INIT) at
!               the initialization phase, the tables are computed for the start date of 
!               the run !
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_J_TABLES_INIT* IS CALLED FROM *BASCOE[TM5]_CHEM_INI*.

!
!     AUTHOR.
!     -------
!        Yves Christophe     *BIRA*


!-----------------------------------------------------------------------
!       Initialize photolysis tables for BASCOE chem
!-----------------------------------------------------------------------
USE YOMLUN             , ONLY : NULOUT
USE PARKIND1           , ONLY : JPIM, JPRB, JPRD              ! (JPRD needed for fcttim.func)
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMRIP0            , ONLY : NINDAT
USE YOMCST             , ONLY : RG
USE BASCOE_J_MODULE    , ONLY :  NDISS, JNAMES
USE BASCOE_J_EXT_MODULE , ONLY : NABSLAYER, LTHICK, WMOLE, TEMPER, DENS, HDENS, &
            & DENS_AIR, DENS_O2
USE BASCOE_J_TABLES_MODULE , ONLY : nsza, njo3, o3du_target, sza_tables, &
            & zjlev_tables, o3col_tables, alj_tables
USE BASCOE_TUV_MODULE  , ONLY :  nabspec, abs_molmass, abs_name, &
            & MXWVN, fbeamr, fbeamr2d, fbeam_dates, daily_solflux, &
            & air, O2abs, O3abs, NOabs, CO2abs, NO2abs

IMPLICIT NONE

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  CHARACTER(LEN=*), PARAMETER    :: CL_MY_NAME   = 'BASCOE_J_TABLES_INIT'

  INTEGER(KIND=JPIM), PARAMETER  :: IMXCLY     = NABSLAYER   ! number of absorbing layers= nb_levels-1
  INTEGER(KIND=JPIM), PARAMETER  :: jprocmax  = NDISS     ! number of photolysis reactions

  REAL(KIND=JPRB), PARAMETER     :: ZR0           = 6371.0E0    ! Effective earth radius (km)
  REAL(KIND=JPRB), PARAMETER     :: ZALBEDO       = 0.25E0      ! Surface albedo of earth
  REAL(KIND=JPRB), PARAMETER     :: ZRGAS          = 8.3144621   ! Perfect gas constant (J/K/mol)
  REAL(KIND=JPRB), PARAMETER     :: zkm_turbo     = 97.         ! altitude of turbopause
  REAL(KIND=JPRB), PARAMETER     :: ZVERY_SMALL    = 1.E-60

  INTEGER(KIND=JPIM)             :: IDXDATE(1)
  INTEGER(KIND=JPIM)             :: iproc, iabs, isza, JLEV, io3
  INTEGER(KIND=JPIM)             :: ILEVTURBO

  REAL(KIND=JPRB), DIMENSION(MXWVN)    :: ZFBEAMR
  REAL(KIND=JPRB)                :: ZDZKM
  REAL(KIND=JPRB), DIMENSION(0:IMXCLY) :: ZS, ZCP_LEV, ZG
  REAL(KIND=JPRB), DIMENSION(0:IMXCLY) :: ZO3REF, ZO3DU_REF
  REAL(KIND=JPRB), DIMENSION(0:IMXCLY,njo3) :: ZO3DU_TMP
  REAL(KIND=JPRB)                :: ZDRAT(nsza,jprocmax,0:IMXCLY)
  REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE

  CHARACTER(LEN=128)             :: CL_FILENM

#include "abor1.intfb.h"
#include "fcttim.func.h"
#include "bascoe_read_coldat.intfb.h"
#include "bascoe_tuv.intfb.h"
#include "bascoe_tuv_ini.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_J_TABLES_INIT',0,ZHOOK_HANDLE )

    !-----------------------------------------------------------------------
    !   ... Vertical profiles zs, g; dens(O2), dens(air) (from MSIS)
    !-----------------------------------------------------------------------

    ZDZKM = LTHICK
    ILEVTURBO = 0
    ZS(0) = NABSLAYER * ZDZKM        ! height of absorption top level
    DO JLEV = 1, IMXCLY
       ZS(JLEV) = ZS(JLEV-1) - ZDZKM
       IF( ZS(JLEV) > zkm_turbo ) ILEVTURBO = JLEV
    ENDDO
    IF( ZS(IMXCLY) < 0. .AND. ZS(IMXCLY) > -1.e-9 ) ZS(IMXCLY) = 0.
    WRITE(NULOUT,*) CL_MY_NAME // ': altitude grid: ZS(0)= ',ZS(0) &
   &   ,' to ZS(IMXCLY)= ',ZS(IMXCLY),' by ',ZDZKM,' km'
    WRITE(NULOUT,*) CL_MY_NAME // ': alt of turbopause(km) : ',ZS(ILEVTURBO)
    ZCP_LEV(0:IMXCLY) = 1.005
    ZG(0:IMXCLY) = 1.d2 * RG * ( ZR0 / (ZR0+ZS(0:IMXCLY)) )**2.                ! cm/s2
    DENS(0:NABSLAYER,O2abs) = DENS_O2
    DENS(0:NABSLAYER,air) = DENS_AIR

    !-----------------------------------------------------------------------
    !   ... Read nb density profiles of other absorbing species
    !-----------------------------------------------------------------------
    DO iabs = O3abs, nabspec
       CL_FILENM = TRIM(ADJUSTL( abs_name(iabs) )) // '_vertprof.dat'
       WRITE(NULOUT,*) CL_MY_NAME // ': reading >' // TRIM( CL_FILENM ) // '<'
       CALL BASCOE_READ_COLDAT( CL_FILENM, IMXCLY+1 &
   &                   , ZS(IMXCLY:0:-1) &
   &                   , DENS(IMXCLY:0:-1,iabs), 'asbndry' )
    ENDDO
    !-----------------------------------------------------------------------
    !       ... Find the total ozone column (Dobson Units) for ref profile
    !-----------------------------------------------------------------------
    ZO3REF(0:IMXCLY) = DENS(0:NABSLAYER,O3abs)
    CALL O3COL( IMXCLY, temper(0), ZG(0), zs, ZO3REF, ZO3DU_REF )
    WRITE(NULOUT,*) CL_MY_NAME // ': surf O3 col for ref O3 prof: ', ZO3DU_REF(IMXCLY)

    !-----------------------------------------------------------------------
    !       ... Scale heights
    !-----------------------------------------------------------------------
    HDENS(0:IMXCLY,air) = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( wmole(0:IMXCLY) * ZG(0:IMXCLY) )
    WRITE(NULOUT,*) CL_MY_NAME // ': air scale height at surface (cm): ', &
   &                                                  HDENS(IMXCLY,air)
    HDENS(0:IMXCLY,O2abs)  = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( abs_molmass(O2abs)  * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,O3abs)  = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( abs_molmass(O3abs)  * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,NOabs)  = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( abs_molmass(NOabs)  * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,CO2abs) = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( abs_molmass(CO2abs) * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,NO2abs) = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( abs_molmass(NO2abs) * ZG(0:IMXCLY) )

    !-----------------------------------------------------------------------
    !   ... Run BASCOE_TUV_INI i.e. initialize the photolysis code
    !-----------------------------------------------------------------------
    CALL BASCOE_TUV_INI( )

    !-----------------------------------------------------------------------
    !   obtain solflux for current date if daily solflux is requested...
    !-----------------------------------------------------------------------
    IF (daily_solflux) THEN
      ! check if date in range
      !
      IF ( NINDAT < fbeam_dates(1) .OR. NINDAT > fbeam_dates(SIZE(fbeam_dates)) ) THEN
        WRITE(NULOUT,*) CL_MY_NAME//' error: date ',NINDAT,'  outside range of solflux dates read from file', &
                    &   fbeam_dates(1), ' to ',  fbeam_dates(SIZE(fbeam_dates))
        CALL ABOR1(CL_MY_NAME//' error: date outside range of solflux dates read from file')
      ENDIF

      ! get closest date 
      IDXDATE = MINLOC( ABS(fbeam_dates - NINDAT ))
      IF ( fbeam_dates(IDXDATE(1)) /= NINDAT ) THEN
        WRITE(NULOUT,*) CL_MY_NAME//' warning solflux date ',fbeam_dates(IDXDATE(1)), &
        &                           ' differs from current date ',NINDAT
      ENDIF
      ZFBEAMR(:) = FBEAMR2D(IDXDATE(1),:)

    ELSE
      ZFBEAMR = FBEAMR

    ENDIF

    !-----------------------------------------------------------------------
    !  Transfer the vertical grid
    !-----------------------------------------------------------------------

    zjlev_tables(1:IMXCLY+1) = zs(IMXCLY:0:-1)

    !-----------------------------------------------------------------------
    !   ... Run BASCOE_TUV to calculate J and heating rates
    !-----------------------------------------------------------------------
    DO io3 = 1, njo3
       DENS(0:IMXCLY,O3abs) = ZO3REF(0:IMXCLY) * o3du_target(io3) / ZO3DU_REF(IMXCLY)
       CALL O3COL( IMXCLY, temper(0), ZG(0), zs(0:IMXCLY), DENS(0:IMXCLY,O3abs), ZO3DU_TMP(0:IMXCLY,io3) )
       !  Transfer the O3 overhead cols corresponding to O3 standard profiles
       o3col_tables(1:IMXCLY+1,1:NJO3) = ZO3DU_TMP(IMXCLY:0:-1,1:NJO3)  

       WRITE(NULOUT,'(2(a,f6.1))') CL_MY_NAME // ': surf O3du = ',ZO3DU_TMP(IMXCLY,io3),' ; target = ', o3du_target(io3)
       DO isza = 1, nsza
          CALL BASCOE_TUV( IMXCLY, ZALBEDO, ILEVTURBO, ZFBEAMR, &
   &                            HDENS, DENS, temper, sza_tables(isza), zs, ZCP_LEV, &
   &                            ZDRAT(isza,1:JPROCMAX,0:IMXCLY) )

          !-----------------------------------------------------------------------
          !   ... Basic error trapping; Set to zero all results < ZVERY_SMALL
          !-----------------------------------------------------------------------
          ! YC: TRICKY CHECK...
          !IF( ANY( ISNAN( ZDRAT(isza,:,:) ) ) .or. ANY( ZDRAT(isza,:,:) < 0.d0 ) ) THEN !  *WARNING* ISNAN is NOT standard !
          IF( ANY(ZDRAT(isza,1:JPROCMAX,0:IMXCLY) /= ZDRAT(isza,1:JPROCMAX,0:IMXCLY)) .OR. &
            & ANY(ZDRAT(isza,1:JPROCMAX,0:IMXCLY) < 0.E0_JPRB ) ) THEN
             WRITE(NULOUT,*) CL_MY_NAME // ', Fatal error: negative or NaN results for sza= ',sza_tables(isza)
             DO iproc = 1, jprocmax
                IF( ANY(ZDRAT(isza,iproc,:) /=  ZDRAT(isza,iproc,:)) ) THEN
                   WRITE(NULOUT,*) '   ISNAN(ZDRAT) for process ',jnames(iproc)
                ENDIF
                IF( ANY( ZDRAT(isza,iproc,:) < 0. ) ) THEN
                   WRITE(NULOUT,*) '   ZDRAT<0 for process ',jnames(iproc),MINVAL(ZDRAT(isza,iproc,:))
                ENDIF
             ENDDO
             CALL ABOR1( CL_MY_NAME//' ERROR: photorates NaN or negative' )
          ENDIF
          WHERE( ZDRAT(isza,:,:) < ZVERY_SMALL )
             ZDRAT(isza,:,:) = 0.d0
          END WHERE

          !-----------------------------------------------------------------------
          !   ... Transfer the results to global J-tables arrays
          !-----------------------------------------------------------------------
          DO iproc = 1, jprocmax
             alj_tables(isza,1:IMXCLY+1,io3,iproc) = ZDRAT(isza,iproc,IMXCLY:0:-1)  !  Transfer to J-table arrays
          ENDDO

       ENDDO    ! isza loop

    ENDDO   !   io3 loop

    !-----------------------------------------------------------------------
    !  Store the (natural) LOG of the J tables - set alj_global to default
    !-----------------------------------------------------------------------
    alj_tables = MAX( 1.E-30_JPRB, alj_tables )
    alj_tables = LOG( alj_tables )

    !-----------------------------------------------------------------------
    ! Diagnostics and error reporting
    !-----------------------------------------------------------------------
    DO isza = 1, nsza
       WRITE(NULOUT,'(a,i3,a,f9.3,a,f15.9)') CL_MY_NAME//': izen= ',isza &
  &          ,' ; sza_tables= ',sza_tables(isza),' ; alj_tables(jlev=31=30km,idiss=1=O2,' &
  &       //'izen,ijo3=2=260DU)= ',alj_tables(isza,31,2,1)
    ENDDO


IF (LHOOK) CALL DR_HOOK('BASCOE_J_TABLES_INIT',1,ZHOOK_HANDLE )


CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! from original BASCOE MODULE J_LIB
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE O3COL( KMXCLY, PTEMPER_TOP, PG_TOP, PZS, PO3DENS, PO3DU )
!-----------------------------------------------------------------------
!   Gives the O3 column (Dobson Units) above each altitude
!-----------------------------------------------------------------------
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN) :: KMXCLY
  REAL(KIND=JPRB), INTENT(IN)  :: PTEMPER_TOP, PG_TOP
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(IN)  :: PZS, PO3DENS
  REAL(KIND=JPRB), DIMENSION(0:KMXCLY), INTENT(OUT) :: PO3DU

  CHARACTER(LEN=*), PARAMETER   :: CL_MY_NAME   = 'O3COL'

  REAL(KIND=JPRB), PARAMETER    :: ZEPSILN = 1.E-5
  REAL(KIND=JPRB), PARAMETER    :: ZRGAS         = 8.3144621   ! Perfect gas constant (J/K/mol)
  REAL(KIND=JPRB), PARAMETER    :: ZO3MOLMASS     = 48.         ! O3 molar mass

  REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE
  REAL(KIND=JPRB) :: ZDELTAZ, ZO3LAYDENS
  INTEGER(KIND=JPIM) :: JLEV

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('BASCOE_J_TABLES_INIT:O3COL',0,ZHOOK_HANDLE )

  PO3DU(0) = PO3DENS(0) * 1.E7_JPRB * ZRGAS * PTEMPER_TOP / ( ZO3MOLMASS * PG_TOP )
  DO JLEV = 1,KMXCLY
    ZDELTAZ = (PZS(JLEV-1) - PZS(JLEV)) * 1.E5_JPRB
!-----------------------------------------------------------------------
!       ...Obtain the densities in each horizontal layer
!            (as in routine LAYDENS of module TUV_MIDATM)
!-----------------------------------------------------------------------
    IF( PO3DENS(JLEV-1) > 0._JPRB .AND. PO3DENS(JLEV) > 0._JPRB .AND. &
&             ABS(1.0 - PO3DENS(JLEV)/PO3DENS(JLEV-1)) > ZEPSILN ) THEN
       ZO3LAYDENS = 1. / (LOG(PO3DENS(JLEV) / PO3DENS(JLEV-1))) * &
&                         (PO3DENS(JLEV) - PO3DENS(JLEV-1)) * ZDELTAZ
    ELSE
       ZO3LAYDENS = &
&           .5_JPRB * ABS( (PO3DENS(JLEV-1) + PO3DENS(JLEV) )*ZDELTAZ )
    ENDIF
    PO3DU(JLEV) = PO3DU(JLEV-1) + ZO3LAYDENS
  ENDDO
  PO3DU(0:KMXCLY) = PO3DU(0:KMXCLY) * 1.E3_JPRB/2.687E19_JPRB        ! molec/cm2 -> Dobson Units

IF (LHOOK) CALL DR_HOOK('BASCOE_J_TABLES_INIT:O3COL',1,ZHOOK_HANDLE )

END SUBROUTINE O3COL

END SUBROUTINE BASCOE_J_TABLES_INIT

