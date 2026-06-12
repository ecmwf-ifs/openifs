SUBROUTINE SUDCDF
USE PARKIND1  ,ONLY : JPIM     ,JPRB,  JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

#ifdef DOC
! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SUDCDF * - Routine to initialize NetCDF output

!     Purpose.
!     --------
!           Initialize NetCDF output files:

!                o_efl.nc (energy fluxes)
!                o_wat.nc (water balance)
!                o_sus.nc (surface state)
!                o_sub.nc (subsurface state)
!                o_eva.nc (evaporation components)
!                o_cld.nc (cold-season processes)
!                o_fix.nc (fixed climate fields)
!                o_lke.nc (lake related fluxes)
!                o_ocd.nc (ocean mixed layer diagnostic variables)
!
!                o_co2.nc (CO2 fluxes)
!                o_bio.nc (biomass)
!		 o_veg.nc (vegetation)
!                o_ext.nc (extra)
!                o_til.nc (tiled output)
!                o_vty.nc (output per vegetation type)

!                o_d2m.nc (inst. 2 meters temp, dewpoint, rel. humi)


!**   Interface.
!     ----------
!        *CALL* *SUDCDF

!     Explicit arguments :
!     --------------------
!        none

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd HURK (KNMI)

!     Modifications.
!     --------------
!        Original     : 2000-7-14
!        E. Dutra  :  Add new output file for lakes data (o_lke.nc)  : 2008-07-04
!        E. Dutra    :11-2009  snow 2009
!        S. Boussetta/G.Balsamo May 2010 Add CTESSEL based on:  
!	             Marita Voogt (KNMI) new tiling (13), vegetation types (7)
!                            and climatology (dec 2004)
!                    Sebastien LAFONT (ECMWF) LAI,fract,Z0 are read for all the vegetation type
!                                 Then the value are computed only 
!                                 for the dominant low and dominant high
!        E. Dutra : 12/2013 :new o_d2m.nc output 
!        E. Dutra : June 2014 - netcf4 
!        A. Agusti-Panareda: 2020-11-17 Atmospheric CO2 forcing

#endif

USE YOMLUN1S , ONLY : NPOSEFL  ,NPOSSUS   ,NPOSWAT   ,NULOUT &
                    &,NPOSSUB  ,NPOSEVA   ,NPOSCLD   ,NPOSCLM &
                    &,NPOSLKE  ,NPOSOCD &                       !KPP
                    &,NPOSCO2  ,NPOSVEG  ,NPOSEXT    ,NPOSBIO &
		    &,NPOSTIL  ,NPOSVTY,NPOSD2M,NPOSGGD
    
USE YOMDPHY  , ONLY : NCSS, NLAT     ,NLON, NCOM &    !KPP
                     &,NVHILO, NTILES

USE YOMLOG1S , ONLY : LACCUMW  ,LRESET    ,CMODID &
                     &,CVERID   ,NACCUR    ,NDIMCDF &
                    &,LWREFL   ,LWRWAT   ,LWRSUB   ,LWRSUS   ,LWREVA &
                    &,LWRCLD   ,LWRCLM   ,LWRLKE   ,LWROCD &    !KPP
                    &,LWRCO2   ,LWRVEG   ,LWREXT &
	            &,LWRBIO   ,LWRTIL   ,LWRVTY,LWRD2M,LWRGGD

USE YOMLUN1S , ONLY : NCTYPE
USE NETCDF
USE NETCDF_UTILS, ONLY: NCERROR,TYPE_NCDF_VAR,INIT_NCDF_VAR,NC_DEF_VAR,&
                        DEF_NC_FILE,TYPE_NCDF_DIMSID
USE MPL_MODULE
IMPLICIT NONE


!Netcdf Interface
INTEGER(KIND=JPIM) :: IDIMID2(2),IDIMID3(3),IDIMID3LEV(3),IDIMID4VEG(4),IDIMID4TIL(4),IDIMID4LEV(4)
INTEGER(KIND=JPIM) :: IDIM3,IDIM2,IDIM4VEG,IDIM4TIL,IDIM3LEV,IDIM4LEV
INTEGER(KIND=JPIM) :: NPOS,VARID

     
!Local variables
INTEGER(KIND=JPIM), PARAMETER ::JPNCDF=16
INTEGER(KIND=JPIM) :: IPOS(JPNCDF)
INTEGER(KIND=JPIM) :: J
LOGICAL LPOS(JPNCDF)

TYPE(TYPE_NCDF_VAR)          :: YD_VARINFO
TYPE(TYPE_NCDF_DIMSID)       :: YD_DIMSID
INTEGER(KIND=JPIM)           :: NVARS2D,NVARS3D,IVAR,NVARS4D
CHARACTER(LEN=20)            :: CVARS2D(40),CVARS3D(40),CVARS4D(40)
INTEGER(KIND=JPIM)           :: MYPROC, NPROC
REAL(KIND=JPHOOK)              :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUDCDF',0,ZHOOK_HANDLE)

ASSOCIATE( NLATID=>YD_DIMSID%NLATID, NLONID=>YD_DIMSID%NLONID, NLEVSID=>YD_DIMSID%NLEVSID, &
         & NTILID=>YD_DIMSID%NTILID, NVTID=>YD_DIMSID%NVTID  , NMONID=>YD_DIMSID%NMONID  , &
         & NTIMID=>YD_DIMSID%NTIMID, NLEVSNID=>YD_DIMSID%NLEVSNID )

MYPROC = MPL_MYRANK()
NPROC  = MPL_NPROC()

IF( MYPROC == 1 ) THEN
  !* -- Initialize file arrays 
  IF(LWREFL)THEN
    CALL NCERROR( NF90_CREATE('o_efl.nc',NCTYPE,NPOS),'Creating file '//'o_efl.nc')
    NPOSEFL = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_efl.nc OPENED ON UNIT ',NPOSEFL
  ENDIF
  IF(LWRWAT)THEN
    CALL NCERROR( NF90_CREATE('o_wat.nc',NCTYPE,NPOS),'Creating file'//'o_wat.nc')
    NPOSWAT = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_wat.nc OPENED ON UNIT ',NPOSWAT
  ENDIF
  IF(LWRSUS)THEN
    CALL NCERROR( NF90_CREATE('o_sus.nc',NCTYPE,NPOS),'Creating file'//'o_sus')
    NPOSSUS = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_sus.nc OPENED ON UNIT ',NPOSSUS
  ENDIF
  IF(LWRSUB)THEN
    WRITE(NULOUT,*)'NETCDF-FILE o_sub.nc not Available ! check previous model versions'
    CALL ABOR1('SUDCDF:')
!   CALL NCERROR( NF90_CREATE('o_sub.nc',NCTYPE,NPOS),'Creating file'//'o_sub')
!   NPOSSUB = NPOS
!   WRITE(NULOUT,*)'NETCDF-FILE o_sub.nc OPENED ON UNIT ',NPOSSUB
  ENDIF
  IF(LWREVA)THEN
    CALL NCERROR( NF90_CREATE('o_eva.nc',NCTYPE,NPOS),'Creating file'//'o_eva.nc')
    NPOSEVA = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_eva.nc OPENED ON UNIT ',NPOSEVA
  ENDIF
  IF(LWRCLD)THEN
    CALL NCERROR( NF90_CREATE('o_cld.nc',NCTYPE,NPOS),'Creating file'//'o_cld.nc')
    NPOSCLD = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_cld.nc OPENED ON UNIT ',NPOSCLD
  ENDIF
  IF(LWRCLM)THEN
    CALL NCERROR( NF90_CREATE('o_fix.nc',NCTYPE,NPOS),'Creating file'//'o_fix.nc')
    NPOSCLM = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_fix.nc OPENED ON UNIT ',NPOSCLM
  ENDIF
  IF(LWRLKE)THEN
    WRITE(NULOUT,*)'NETCDF-FILE o_lke.nc not Available ! check previous model versions'
    CALL ABOR1('SUDCDF:')
!   CALL NCERROR( NF90_CREATE('o_lke.nc',NCTYPE,NPOS),'Creating file'//'o_lke.nc')
!   NPOSLKE = NPOS
!   WRITE(NULOUT,*)'NETCDF-FILE o_lke.nc OPENED ON UNIT ',NPOSLKE
  ENDIF
  IF(LWRCO2)THEN
    CALL NCERROR( NF90_CREATE('o_co2.nc',NCTYPE,NPOS),'Creating file'//'o_co2.nc')
    NPOSCO2 = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_co2.nc OPENED ON UNIT ',NPOSCO2
  ENDIF
  IF(LWRBIO)THEN
    CALL NCERROR( NF90_CREATE('o_bio.nc',NCTYPE,NPOS),'Creating file'//'o_bio.nc')
    NPOSBIO = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_bio.nc OPENED ON UNIT ',NPOSBIO
  ENDIF
  IF(LWRVEG)THEN
    CALL NCERROR( NF90_CREATE('o_veg.nc',NCTYPE,NPOS),'Creating file'//'o_veg.nc')
    NPOSVEG = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_veg.nc OPENED ON UNIT ',NPOSVEG
  ENDIF
  IF(LWREXT)THEN
!   WRITE(NULOUT,*)'NETCDF-FILE o_ext.nc not Available ! check previous model versions'
!   CALL ABOR1('SUDCDF:')
    CALL NCERROR( NF90_CREATE('o_ext.nc',NCTYPE,NPOS),'Creating file'//'o_ext.nc')
    NPOSEXT = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_ext.nc OPENED ON UNIT ',NPOSEXT
  ENDIF
  IF(LWRTIL)THEN
    CALL NCERROR( NF90_CREATE('o_til.nc',NCTYPE,NPOS),'Creating file'//'o_til.nc')
    NPOSTIL = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_til.nc OPENED ON UNIT ',NPOSTIL
  ENDIF
  IF(LWRVTY)THEN
    CALL NCERROR( NF90_CREATE('o_vty.nc',NCTYPE,NPOS),'Creating file'//'o_vty.nc')
    NPOSVTY = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_vty.nc OPENED ON UNIT ',NPOSVTY
  ENDIF
  IF(LWRD2M)THEN
    CALL NCERROR( NF90_CREATE('o_d2m.nc',NCTYPE,NPOS),'Creating file'//'o_d2m.nc')
    NPOSD2M = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_d2m.nc OPENED ON UNIT ',NPOSD2M
  ENDIF
  IF(LWRD2M)THEN
    CALL NCERROR( NF90_CREATE('o_ggd.nc',NCTYPE,NPOS),'Creating file'//'o_ggd.nc')
    NPOSGGD = NPOS
    WRITE(NULOUT,*)'NETCDF-FILE o_ggd.nc OPENED ON UNIT ',NPOSGGD
  ENDIF

  IPOS(1)=NPOSEFL
  IPOS(2)=NPOSWAT
  IPOS(3)=NPOSSUS
  IPOS(4)=NPOSSUB
  IPOS(5)=NPOSEVA
  IPOS(6)=NPOSCLD
  IPOS(7)=NPOSCLM
  IPOS(8)=NPOSLKE
  IPOS(9)=NPOSCO2
  IPOS(10)=NPOSBIO
  IPOS(11)=NPOSVEG
  IPOS(12)=NPOSEXT
  IPOS(13)=NPOSTIL
  IPOS(14)=NPOSVTY
  IPOS(15)=NPOSD2M
  IPOS(16)=NPOSGGD

  LPOS(1)=LWREFL
  LPOS(2)=LWRWAT
  LPOS(3)=LWRSUS
  LPOS(4)=LWRSUB
  LPOS(5)=LWREVA
  LPOS(6)=LWRCLD
  LPOS(7)=LWRCLM
  LPOS(8)=LWRLKE
  LPOS(9)=LWRCO2 !CTESSEL
  LPOS(10)=LWRBIO !CTESSEL
  LPOS(11)=LWRVEG !CTESSEL
  LPOS(12)=LWREXT !CTESSEL
  LPOS(13)=LWRTIL !CTESSEL
  LPOS(14)=LWRVTY !CTESSEL
  LPOS(15)=LWRD2M !CTESSEL
  LPOS(16)=LWRGGD !CTESSEL


!*  Main loop on files 
  DO J=1,JPNCDF
    IF (.NOT. LPOS(J)) CYCLE 
    NPOS = IPOS(J)

  !* -- 2. add header to the file and global attributes
    CALL DEF_NC_FILE(NPOS,YD_DIMSID)

    !! Dimensions definitions 
    IF(NDIMCDF == 2)THEN
      IDIMID3(1) = NLONID
      IDIMID3(2) = NLATID
      IDIMID3(3) = NTIMID
      IDIM3=3

      IDIMID3LEV(1) = NLONID
      IDIMID3LEV(2) = NLATID
      IDIMID3LEV(3) = NLEVSID
      IDIM3LEV=3
    
      IDIMID4VEG(1) = NLONID
      IDIMID4VEG(2) = NLATID
      IDIMID4VEG(3) = NVTID
      IDIMID4VEG(4) = NTIMID
      IDIM4VEG=4

      IDIMID4TIL(1) = NLONID
      IDIMID4TIL(2) = NLATID
      IDIMID4TIL(3) = NTILID
      IDIMID4TIL(4) = NTIMID
      IDIM4TIL=4

      IDIMID4LEV(1) = NLONID
      IDIMID4LEV(2) = NLATID
      IDIMID4LEV(3) = NLEVSID
      IDIMID4LEV(4) = NTIMID
      IDIM4LEV=4
    ELSE
      IDIMID3(1)=NLONID
      IDIMID3(2)=NTIMID
      IDIM3=2

      IDIMID3LEV(1) = NLONID
      IDIMID3LEV(2) = NLEVSID
      IDIM3LEV=2
    

      IDIMID4VEG(1) = NLONID
      IDIMID4VEG(2) = NVTID
      IDIMID4VEG(3) = NTIMID
      IDIM4VEG=3

      IDIMID4TIL(1) = NLONID
      IDIMID4TIL(2) = NTILID
      IDIMID4TIL(3) = NTIMID
      IDIM4TIL=3

      IDIMID4LEV(1) = NLONID
      IDIMID4LEV(2) = NLEVSID
      IDIMID4LEV(3) = NTIMID
      IDIM4LEV=3
    ENDIF

! * -- prognostics - mean values 
!******************************************************************
    IF( NPOS == NPOSGGD )THEN
      NVARS3D=2
      CVARS3D(1:NVARS3D)=(/'SoilMoist',&
                           'SoilTemp '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time nlevs lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID4LEV(1:IDIM4LEV),VARID)
      ENDDO

      NVARS3D=14
      CVARS3D(1:NVARS3D)=(/'AvgSurfT ','CanopInt ','SWE      ','SnowT    ',&
                           'SAlbedo  ','snowdens ','TLICE    ','TLMNW    ',&
                           'TLWML    ','TLBOT    ','TLSF     ','HLICE    ',&
                           'HLML     ','slw      '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! ggd

!* -- Surface energy balance
!******************************************************************
    IF(NPOS == NPOSEFL)THEN
      NVARS3D=13
      CVARS3D(1:NVARS3D)=(/'SWnet      ','LWnet      ','Qle        ',&
                           'Qh         ','Qg         ','Qf         ',&
                           'DelSoilHeat','DelColdCont','LWup       ',&
                           'SWdown     ','LWdown     ','Qgsn       ',&
                           'Qfsn       '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! efl 

! * -- Surface water balance 
!******************************************************************
  IF( NPOS == NPOSWAT )THEN
    NVARS3D=11
    CVARS3D(1:NVARS3D)=(/'Snowf       ','Rainf       ','Evap        ',&
                         'Qs          ','Qsb         ','Qsm         ',&
                         'DelSoilMoist','DelSWE      ','DelIntercept',&
                         'Intercept   ','fldfrc      '/)
    DO IVAR=1,NVARS3D
      CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
      CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
  ENDDO

  ENDIF ! wat

! * -- diagnostics at 2 m (different writing output timming !)
!******************************************************************
    IF( NPOS == NPOSD2M )THEN
      NVARS3D=3
      CVARS3D(1:NVARS3D)=(/'T2m   ','RH2m  ','D2m   '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! d2m

! * -- surface state variables
!******************************************************************
    IF( NPOS == NPOSSUS )THEN
      NVARS3D=7
      CVARS3D(1:NVARS3D)=(/'VegT        ','BaresoilT   ','RadT        ',&
                           'Albedo      ','lai         ','z0m         ',&
                           'lz0h        '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! sus

! * -- evaporation components
!******************************************************************
    IF( NPOS == NPOSEVA )THEN
      NVARS3D=10
      CVARS3D(1:NVARS3D)=(/'ECanop      ','TVeg        ','Conds       ',&
                           'ESoil       ','RootMoist   ','SubSnow     ',&
                           'PotEvapI    ','PotEvapWA   ','EWater      ',&
                           'PotEvapU    '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! eva

! * -- cold processes 
!******************************************************************
    IF( NPOS == NPOSCLD )THEN
      NVARS3D=5
      CVARS3D(1:NVARS3D)=(/'SnowFrac    ','IceFrac     ','Fdepth      ',&
                           'Tdepth      ','SnowDepth   '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! CLD

! * -- CO2 fluxes 
!******************************************************************
    IF( NPOS == NPOSCO2 )THEN
      NVARS3D=7
      CVARS3D(1:NVARS3D)=(/'Ag          ','Rd          ','An          ',&
                           'Rsoil_str   ','Reco        ','CO2flux     ','CH4flux     '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! co2

! * -- Biomass
!******************************************************************
    IF( NPOS == NPOSBIO )THEN
      NVARS3D=6
      CVARS3D(1:NVARS3D)=(/'lai         ','biomass     ','Bloss       ',&
                           'Bgain       ','Biomstr     ','Biomstr2    '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
      ENDDO

    ENDIF ! bio

! * -- vegetation variables 
!******************************************************************
    IF( NPOS == NPOSVEG )THEN
      NVARS4D=6
      CVARS4D(1:NVARS4D)=(/'rc          ','ra          ','LEvt        ',&
                         'f2vt        ','dsvt        ','dmaxvt      '/)
      DO IVAR=1,NVARS4D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS4D(IVAR)),CCOORD_IN="time vtype lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID4VEG(1:IDIM4VEG),VARID)
      ENDDO

    ENDIF ! veg

! * -- vegetation variables 2
!******************************************************************
    IF( NPOS == NPOSVTY )THEN
      NVARS4D=14
      CVARS4D(1:NVARS4D)=(/'vtfr        ','lai         ','biomass     ',&
                           'Bloss       ','Bgain       ','Biomstr     ',&
                           'Biomstr2    ','Ag          ','Rd          ',&
                           'An          ','Rsoil_str   ','Reco        ',&
                         'CO2flux     ','RnoQ10      '/)
      DO IVAR=1,NVARS4D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS4D(IVAR)),CCOORD_IN="time vtype lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID4VEG(1:IDIM4VEG),VARID)
      ENDDO

  ENDIF ! veg

! * -- tiled information
!******************************************************************
    IF( NPOS == NPOSTIL )THEN
      NVARS4D=9
      CVARS4D(1:NVARS4D)=(/'tifr        ','SWnet       ','LWnet       ',&
                           'Qle         ','Qh          ','LWup        ',&
                           'Evap        ','SkinT       ','SkinC       '/)
      DO IVAR=1,NVARS4D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS4D(IVAR)),CCOORD_IN="time tile lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID4TIL(1:IDIM4TIL),VARID)
      ENDDO

    ENDIF ! til

! * -- fixed climate fields 
!******************************************************************
    IF( NPOS == NPOSCLM )THEN
      NVARS3D=5
      CVARS3D(1:NVARS3D)=(/'SoilThick   ','SoilFC      ','SoilWilt    ',&
                           'SoilSat     ','RootDist    '/)
      DO IVAR=1,NVARS3D
        CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="nlevs lat lon")
        CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3LEV(1:IDIM3LEV),VARID)
      ENDDO

    ENDIF ! fix
  
  
! * -- Extra variables, forcing
!******************************************************************
  IF( NPOS == NPOSEXT )THEN
    NVARS3D=6
    CVARS3D(1:NVARS3D)=(/'Tair  ','Qair  ','Uwind ','Vwind ','psurf ','CO2air'/)
    DO IVAR=1,NVARS3D
      CALL INIT_NCDF_VAR(YD_VARINFO,TRIM(CVARS3D(IVAR)),CCOORD_IN="time lat lon")
      CALL NC_DEF_VAR(NPOS,YD_VARINFO,IDIMID3(1:IDIM3),VARID)
  ENDDO

     ENDIF ! sus


  !* -- end of definition mode
    CALL NCERROR( NF90_ENDDEF(NPOS) )
  ENDDO ! loop on files 

ENDIF

END ASSOCIATE

CALL MPL_BARRIER()

IF (LHOOK) CALL DR_HOOK('SUDCDF',1,ZHOOK_HANDLE)

END SUBROUTINE SUDCDF
