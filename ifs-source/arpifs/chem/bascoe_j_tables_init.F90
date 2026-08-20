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
USE BASCOE_J_TABLES_MODULE , ONLY : NSZA, NJO3, O3DU_TARGET, SZA_TABLES, &
            & ZJLEV_TABLES, O3COL_TABLES, ALJ_TABLES
USE BASCOE_TUV_MODULE  , ONLY :  NABSPEC, ABS_MOLMASS, ABS_NAME, &
            & MXWVN, FBEAMR, FBEAMR2D, FBEAM_DATES, DAILY_SOLFLUX, &
            & AIR, O2ABS, O3ABS, NOABS, CO2ABS, NO2ABS

IMPLICIT NONE

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  CHARACTER(LEN=*), PARAMETER    :: CL_MY_NAME   = 'BASCOE_J_TABLES_INIT'

  INTEGER(KIND=JPIM), PARAMETER  :: IMXCLY     = NABSLAYER   ! number of absorbing layers= nb_levels-1
  INTEGER(KIND=JPIM), PARAMETER  :: JPROCMAX   = NDISS       ! number of photolysis reactions

  REAL(KIND=JPRB), PARAMETER     :: ZR0            = 6371.0E0    ! Effective earth radius (km)
  REAL(KIND=JPRB), PARAMETER     :: ZALBEDO        = 0.25E0      ! Surface albedo of earth
  REAL(KIND=JPRB), PARAMETER     :: ZRGAS          = 8.3144621   ! Perfect gas constant (J/K/mol)
  REAL(KIND=JPRB), PARAMETER     :: ZKM_TURBO      = 97.         ! altitude of turbopause
  REAL(KIND=JPRB), PARAMETER     :: ZVERY_SMALL    = 1.E-36

  INTEGER(KIND=JPIM)             :: IDXDATE(1)
  INTEGER(KIND=JPIM)             :: IPROC, IABS, ISZA, JLEV, IO3
  INTEGER(KIND=JPIM)             :: ILEVTURBO

  REAL(KIND=JPRB), DIMENSION(MXWVN)    :: ZFBEAMR
  REAL(KIND=JPRB)                :: ZDZKM
  REAL(KIND=JPRB), DIMENSION(0:IMXCLY) :: ZS, ZCP_LEV, ZG
  REAL(KIND=JPRB), DIMENSION(0:IMXCLY) :: ZO3REF, ZO3DU_REF
  REAL(KIND=JPRB), DIMENSION(0:IMXCLY,NJO3) :: ZO3DU_TMP
  REAL(KIND=JPRB)                :: ZDRAT(NSZA,JPROCMAX,0:IMXCLY)
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
       IF( ZS(JLEV) > ZKM_TURBO ) ILEVTURBO = JLEV
    ENDDO
    IF( ZS(IMXCLY) < 0. .AND. ZS(IMXCLY) > -1.e-9 ) ZS(IMXCLY) = 0.
    WRITE(NULOUT,*) CL_MY_NAME // ': altitude grid: ZS(0)= ',ZS(0) &
   &   ,' to ZS(IMXCLY)= ',ZS(IMXCLY),' by ',ZDZKM,' km'
    WRITE(NULOUT,*) CL_MY_NAME // ': alt of turbopause(km) : ',ZS(ILEVTURBO)
    ZCP_LEV(0:IMXCLY) = 1.005
    ZG(0:IMXCLY) = 1.d2 * RG * ( ZR0 / (ZR0+ZS(0:IMXCLY)) )**2.                ! cm/s2
    DENS(0:NABSLAYER,O2ABS) = DENS_O2
    DENS(0:NABSLAYER,AIR) = DENS_AIR

    !-----------------------------------------------------------------------
    !   ... Read nb density profiles of other absorbing species
    !-----------------------------------------------------------------------
    DO IABS = O3ABS, NABSPEC
       CL_FILENM = TRIM(ADJUSTL( ABS_NAME(IABS) )) // '_vertprof.dat'
       WRITE(NULOUT,*) CL_MY_NAME // ': reading >' // TRIM( CL_FILENM ) // '<'
       CALL BASCOE_READ_COLDAT( CL_FILENM, IMXCLY+1 &
   &                   , ZS(IMXCLY:0:-1) &
   &                   , DENS(IMXCLY:0:-1,IABS), 'asbndry' )
    ENDDO
    !-----------------------------------------------------------------------
    !       ... Find the total ozone column (Dobson Units) for ref profile
    !-----------------------------------------------------------------------
    ZO3REF(0:IMXCLY) = DENS(0:NABSLAYER,O3ABS)
    CALL O3COL( IMXCLY, temper(0), ZG(0), zs, ZO3REF, ZO3DU_REF )
    WRITE(NULOUT,*) CL_MY_NAME // ': surf O3 col for ref O3 prof: ', ZO3DU_REF(IMXCLY)

    !-----------------------------------------------------------------------
    !       ... Scale heights
    !-----------------------------------------------------------------------
    HDENS(0:IMXCLY,AIR) = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( wmole(0:IMXCLY) * ZG(0:IMXCLY) )
    WRITE(NULOUT,*) CL_MY_NAME // ': air scale height at surface (cm): ', &
   &                                                  HDENS(IMXCLY,AIR)
    HDENS(0:IMXCLY,O2ABS)  = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( ABS_MOLMASS(O2ABS)  * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,O3ABS)  = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( ABS_MOLMASS(O3ABS)  * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,NOABS)  = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( ABS_MOLMASS(NOABS)  * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,CO2ABS) = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( ABS_MOLMASS(CO2ABS) * ZG(0:IMXCLY) )
    HDENS(0:IMXCLY,NO2ABS) = 1.E7_JPRB * ZRGAS * temper(0:IMXCLY) / ( ABS_MOLMASS(NO2ABS) * ZG(0:IMXCLY) )

    !-----------------------------------------------------------------------
    !   ... Run BASCOE_TUV_INI i.e. initialize the photolysis code
    !-----------------------------------------------------------------------
    CALL BASCOE_TUV_INI( )

    !-----------------------------------------------------------------------
    !   obtain solflux for current date if daily solflux is requested...
    !-----------------------------------------------------------------------
    IF (DAILY_SOLFLUX) THEN
      ! check if date in range
      !
      IF ( NINDAT < FBEAM_DATES(1) .OR. NINDAT > FBEAM_DATES(SIZE(FBEAM_DATES)) ) THEN
        WRITE(NULOUT,*) CL_MY_NAME//' error: date ',NINDAT,'  outside range of solflux dates read from file', &
                    &   FBEAM_DATES(1), ' to ',  FBEAM_DATES(SIZE(FBEAM_DATES))
        CALL ABOR1(CL_MY_NAME//' error: date outside range of solflux dates read from file')
      ENDIF

      ! get closest date 
      IDXDATE = MINLOC( ABS(FBEAM_DATES - NINDAT ))
      IF ( FBEAM_DATES(IDXDATE(1)) /= NINDAT ) THEN
        WRITE(NULOUT,*) CL_MY_NAME//' warning solflux date ',FBEAM_DATES(IDXDATE(1)), &
        &                           ' differs from current date ',NINDAT
      ENDIF
      ZFBEAMR(:) = FBEAMR2D(IDXDATE(1),:)

    ELSE
      ZFBEAMR = FBEAMR

    ENDIF

    !-----------------------------------------------------------------------
    !  Transfer the vertical grid
    !-----------------------------------------------------------------------

    ZJLEV_TABLES(1:IMXCLY+1) = zs(IMXCLY:0:-1)

    !-----------------------------------------------------------------------
    !   ... Run BASCOE_TUV to calculate J and heating rates
    !-----------------------------------------------------------------------
    DO IO3 = 1, NJO3
       DENS(0:IMXCLY,O3ABS) = ZO3REF(0:IMXCLY) * O3DU_TARGET(IO3) / ZO3DU_REF(IMXCLY)
       CALL O3COL( IMXCLY, temper(0), ZG(0), zs(0:IMXCLY), DENS(0:IMXCLY,O3ABS), ZO3DU_TMP(0:IMXCLY,IO3) )
       !  Transfer the O3 overhead cols corresponding to O3 standard profiles
       O3COL_TABLES(1:IMXCLY+1,1:NJO3) = ZO3DU_TMP(IMXCLY:0:-1,1:NJO3)  

       WRITE(NULOUT,'(2(a,f6.1))') CL_MY_NAME // ': surf O3du = ',ZO3DU_TMP(IMXCLY,IO3),' ; target = ', O3DU_TARGET(IO3)
       DO ISZA = 1, NSZA
          CALL BASCOE_TUV( IMXCLY, ZALBEDO, ILEVTURBO, ZFBEAMR, &
   &                            HDENS, DENS, temper, SZA_TABLES(ISZA), zs, ZCP_LEV, &
   &                            ZDRAT(ISZA,1:JPROCMAX,0:IMXCLY) )

          !-----------------------------------------------------------------------
          !   ... Basic error trapping; Set to zero all results < ZVERY_SMALL
          !-----------------------------------------------------------------------
          ! YC: TRICKY CHECK...
          !IF( ANY( ISNAN( ZDRAT(ISZA,:,:) ) ) .or. ANY( ZDRAT(ISZA,:,:) < 0.0_JPRB ) ) THEN !  *WARNING* ISNAN is NOT standard !
          IF( ANY(ZDRAT(ISZA,1:JPROCMAX,0:IMXCLY) /= ZDRAT(ISZA,1:JPROCMAX,0:IMXCLY)) .OR. &
            & ANY(ZDRAT(ISZA,1:JPROCMAX,0:IMXCLY) < 0.E0_JPRB ) ) THEN
             WRITE(NULOUT,*) CL_MY_NAME // ', Fatal error: negative or NaN results for sza= ',SZA_TABLES(ISZA)
             DO IPROC = 1, JPROCMAX
                IF( ANY(ZDRAT(ISZA,IPROC,:) /=  ZDRAT(ISZA,IPROC,:)) ) THEN
                   WRITE(NULOUT,*) '   ISNAN(ZDRAT) for process ',jnames(IPROC)
                ENDIF
                IF( ANY( ZDRAT(ISZA,IPROC,:) < 0._JPRB ) ) THEN
                   WRITE(NULOUT,*) '   ZDRAT<0 for process ',jnames(IPROC),MINVAL(ZDRAT(ISZA,IPROC,:))
                ENDIF
             ENDDO
             CALL ABOR1( CL_MY_NAME//' ERROR: photorates NaN or negative' )
          ENDIF
          WHERE( ZDRAT(ISZA,:,:) < ZVERY_SMALL )
             ZDRAT(ISZA,:,:) = 0.0_JPRB
          END WHERE

          !-----------------------------------------------------------------------
          !   ... Transfer the results to global J-tables arrays
          !-----------------------------------------------------------------------
          DO IPROC = 1, JPROCMAX
             ALJ_TABLES(ISZA,1:IMXCLY+1,IO3,IPROC) = ZDRAT(ISZA,IPROC,IMXCLY:0:-1)  !  Transfer to J-table arrays
          ENDDO

       ENDDO    ! ISZA loop

    ENDDO   !   IO3 loop

    !-----------------------------------------------------------------------
    !  Store the (natural) LOG of the J tables - set alj_global to default
    !-----------------------------------------------------------------------
    ALJ_TABLES = MAX( 1.E-30_JPRB, ALJ_TABLES )
    ALJ_TABLES = LOG( ALJ_TABLES )

    !-----------------------------------------------------------------------
    ! Diagnostics and error reporting
    !-----------------------------------------------------------------------
    DO ISZA = 1, NSZA
       WRITE(NULOUT,'(a,i3,a,f9.3,a,f15.9)') CL_MY_NAME//': izen= ',ISZA &
  &          ,' ; SZA_TABLES= ',SZA_TABLES(ISZA),' ; ALJ_TABLES(jlev=31=30km,idiss=1=O2,' &
  &       //'izen,ijo3=2=260DU)= ',ALJ_TABLES(ISZA,31,2,1)
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

