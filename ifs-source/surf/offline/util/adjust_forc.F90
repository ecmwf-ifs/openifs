! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!
program adjust_for
use grib_api
implicit none
integer                          :: igribQ
integer                          :: igribT
integer                          :: igribP
integer                          :: igribDZ
integer                          :: igribDTDZ
integer                          :: ifileQ
integer                          :: ifileT
integer                          :: ifileP
integer                          :: ifileDZ
integer                          :: ifileDTDZ
integer                          :: ofileT
integer                          :: ofileQ
integer                          :: ofileP
integer                          :: i
integer                          :: iret
integer                          :: iretQ
integer                          :: iretT
integer                          :: iretP
integer                          :: iretDZ
integer                          :: iretDTDZ
integer                          :: j
integer                          :: jmax=100
integer                          :: nb_pv
integer(4)                       :: nlev
integer(4)                       :: nlev1
character(len=256)               :: filename
character(len=256)               :: outfile
character(len=8)                 :: cdateQ
character(len=4)                 :: ctimeQ
character(len=3)                 :: censmemQ
integer(4)                       :: dataDateQ
integer(4)                       :: dataTimeQ
integer(4)                       :: perturbationNumberQ
integer(4)                       :: numberOfPointsQ
real(8)                          :: missingValueQ=9999
integer(4)                       :: dataDateT
integer(4)                       :: dataTimeT
integer(4)                       :: perturbationNumberT
integer(4)                       :: numberOfPointsT
real(8)                          :: missingValueT=9999
real(8)                          :: constA_lowml
real(8)                          :: constB_lowml
real(8),dimension(:),allocatable :: pv
real(8),dimension(:),allocatable :: constA
real(8),dimension(:),allocatable :: constB
real(8),dimension(:),allocatable :: lats
real(8),dimension(:),allocatable :: lons
real(8),dimension(:),allocatable :: Qair        ! Specific Humidity on Lowest Full Model Level
real(8),dimension(:),allocatable :: Tair        ! Temperature on Lowest Full Model Level
real(8),dimension(:),allocatable :: lnps        ! Log Surface Pressure [Pa]
real(8),dimension(:),allocatable :: DZ          ! Orography difference [m]
real(8),dimension(:),allocatable :: dtdz        ! Temperature Lapse Rate [K/m]
real(8)                          :: Pair        ! Pressure on Lowest Full Model Level
real(8)                          :: Tair_new    ! Temperature on Lowest Full Model Level after correction
real(8)                          :: Pair_new    ! Pressure on Lowest Full Model Level after correction
real(8)                          :: Qair_new    ! Specific Humidity on Lowest Full Model Level after correction 
real(8)                          :: Ps_new      ! Surface Pressure after correction
real(8)                          :: Tv          ! Virtual Temperature on Lowest Full Model Level
real(8)                          :: Qair_fg  
real(8)                          :: dQ
real(8)                          :: RG
real(8)                          :: PRHMAX
real(8)                          :: PRHMIN
real(8)                          :: RKBOL
real(8)                          :: RNAVO
real(8)                          :: R
real(8)                          :: RMD
real(8)                          :: RMV
real(8)                          :: RD
real(8)                          :: RV
real(8)                          :: RETV
real(8)                          :: PES     ! Saturation Pressure(es)
real(8)                          :: ZRH
real(8)                          :: PRH     ! Relative Humidity


!--- Prepare Constants ---
RG=9.80665d0          ! Gravitational Acceleration [m/s2]
PRHMAX=1.d0           ! Maximum Relative Humidity
PRHMIN=0.d0           ! Minimum Relative Humidity
RKBOL=1.380658E-23    ! Boltzmann constant (k) [J/K]
RNAVO=6.0221367E+23   ! Avogadro constant (N_A) [1/mol]
R=RNAVO*RKBOL         ! Perfect gas constant (= 8.314511211948600)
RMD=28.9644d0         ! Dry air mass
RMV=18.0153d0         ! Vapour  mass
RD=1000.d0*R/RMD      ! Dry air cst. (= 287.0596736665907 J/kg/K)
RV=1000.d0*R/RMV      ! Vapour  cst. (= 461.5249933083879 J/kg/K)
RETV=RV/RD-1.0d0      ! Rv/Rd-1 (= 0.608)


!--- Open Grib file for input ---
call grib_open_file(ifileQ, 'INPUTGRIB_Q','R')
call grib_open_file(ifileT, 'INPUTGRIB_T','R')
call grib_open_file(ifileP, 'INPUTGRIB_LNPS','R')
call grib_open_file(ifileDZ, 'INPUTGRIB_DZ','R')
call grib_open_file(ifileDTDZ, 'INPUTGRIB_DTDZ','R')
!--- Open Grib file for output ---
call grib_open_file(ofileT, 'OUTPUTGRIB_T','W')
call grib_open_file(ofileQ, 'OUTPUTGRIB_Q','W')
call grib_open_file(ofileP, 'OUTPUTGRIB_LNPS','W')


write(6,*) 
write(6,*) "adjustment of forcing data due to orography difference"

!--- Identify grib file ---
call grib_new_from_file(ifileQ,igribQ,iretQ)
call grib_new_from_file(ifileT,igribT,iretT)
call grib_new_from_file(ifileP,igribP,iretP)
call grib_new_from_file(ifileDZ,igribDZ,iretDZ)
call grib_new_from_file(ifileDTDZ,igribDTDZ,iretDTDZ)

!--- Get number of vertical layers and constant A and B ---
call grib_get_size(igribT,'pv',nb_pv)
allocate(pv(nb_pv))
call grib_get(igribT,'pv',pv)
nlev1=nb_pv/2
nlev=nlev1-1
allocate(constA(nlev1))
allocate(constB(nlev1))
constA(1:nlev1)=pv(1:nlev1)
constB(1:nlev1)=pv(nlev1+1:nb_pv)

!--- Constant A and B at lowest full model level ---
constA_lowml=constA(nlev)
constB_lowml=constB(nlev)

!--- Loop on all the messages in a file ---
do while (iretQ/=GRIB_END_OF_FILE .and. iretT/=GRIB_END_OF_FILE)
  !--- Get observation date and time ---
  call grib_get(igribQ,'dataDate',dataDateQ)
  call grib_get(igribQ,'dataTime',dataTimeQ)
  call grib_get(igribQ,'perturbationNumber',perturbationNumberQ)
  call grib_get(igribQ,'numberOfPoints',numberOfPointsQ)
  call grib_set(igribQ,'missingValue',missingValueQ)
  call grib_get(igribT,'dataDate',dataDateT)
  call grib_get(igribT,'dataTime',dataTimeT)
  call grib_get(igribT,'perturbationNumber',perturbationNumberT)
  call grib_get(igribT,'numberOfPoints',numberOfPointsT)
  call grib_set(igribT,'missingValue',missingValueT)


  if (dataDateQ==dataDateT .and. dataTimeQ==dataTimeT .and. &
     & perturbationNumberQ==perturbationNumberT .and. &
     & numberOfPointsQ==numberOfPointsT .and. missingValueQ==missingValueT) then

    write(cdateQ,'(i8.8)') dataDateQ
    write(ctimeQ,'(i4.4)') dataTimeQ
    write(censmemQ,'(i3.3)') perturbationNumberQ
    write(6,*) "processing... initialtime:",cdateQ,"_",ctimeQ,"UTC, member:",censmemQ

    allocate(lats(numberOfPointsQ))
    allocate(lons(numberOfPointsQ))
    allocate(Qair(numberOfPointsQ))
    allocate(Tair(numberOfPointsT))
    allocate(lnps(numberOfPointsT))
    allocate(DZ(numberOfPointsT))
    allocate(dtdz(numberOfPointsT))
 
    call grib_get_data(igribQ,lats,lons,Qair)
    call grib_get_data(igribT,lats,lons,Tair)
    call grib_get_data(igribP,lats,lons,lnps)
    call grib_get_data(igribDZ,lats,lons,DZ)
    call grib_get_data(igribDTDZ,lats,lons,dtdz)

    do i=1,numberOfPointsQ
      if (Qair(i) /= missingValueQ) then
        !===============================================================
        !  Calculation of relative humidity at lowest full model level
        !===============================================================
        !--- Pressure at lowest full model level ===
        Pair = (constA_lowml+constB_lowml*exp(lnps(i)) &
                &   + exp(lnps(i))) / 2.d0
        ! saturation vapor pressure with Teten's formula : eq(7.5) in IFS Document Cy38R1
        call cal_PES(Tair(i),PES)
        ! relative humidity : eq(7.72) in IFS Document Cy38R1
        ZRH=(Pair*Qair(i)*(RETV+1.0d0))&
          & /((1.0d0+RETV*Qair(i))*PES) 
        PRH=MAX(PRHMIN,MIN(ZRH,PRHMAX))

        !==============================================================
        !  Calculation of corrected variables for Tair
        !==============================================================
        call correct_Tair(Tair(i),dtdz(i),DZ(i),Tair_new)

        Qair_fg = Qair(i)
        do j = 1, jmax
          !==============================================================
          !  Calculation of corrected variables for lnps
          !==============================================================
          call correct_lnps(Tair(i),Tair_new,Qair(i),Qair_fg,lnps(i),DZ(i),RETV,RG,RD,Ps_new)

          !==============================================================
          !  Calculation of corrected variables for Qair
          !==============================================================
          call correct_Qair(Tair_new,Ps_new,constA_lowml,constB_lowml,RETV,PRH,Qair_new)

          dQ = abs(Qair_new - Qair_fg)
          if ( dQ > 0.0000001d0 ) then
            Qair_fg = Qair_fg + (Qair_new - Qair_fg)*0.5d0
          else
            exit
          endif 
        enddo

        ! for Output
        Tair(i) = Tair_new
        Qair(i) = Qair_new
        lnps(i) = log(Ps_new)
      end if
    enddo


    !==============================================================
    !  Outout of corrected variables
    !==============================================================
    ! write values to output grib file (T)
    call grib_set(igribT,"values",Tair)
    call grib_write(igribT,ofileT,iret)
    
    ! write values to output grib file (LNPS)
    call grib_set(igribP,"values",lnps)
    call grib_write(igribP,ofileP,iret)
    
    ! write values to output grib file (Q)
    call grib_set(igribQ,"values",Qair)
    call grib_write(igribQ,ofileQ,iret)
    

    deallocate(lats)
    deallocate(lons)
    deallocate(Qair)
    deallocate(Tair)
    deallocate(DZ)
    deallocate(dtdz)
    deallocate(lnps)

    call grib_release(igribQ)
    call grib_new_from_file(ifileQ,igribQ, iretQ)
    call grib_new_from_file(ifileT,igribT, iretT)
    call grib_new_from_file(ifileP,igribP, iretP)
    call grib_new_from_file(ifileDZ,igribDZ,iretDZ)
    call grib_new_from_file(ifileDTDZ,igribDTDZ,iretDTDZ)
  endif
end do 

call grib_close_file(ifileQ)
call grib_close_file(ifileT)
call grib_close_file(ifileP)
call grib_close_file(ifileDZ)
call grib_close_file(ifileDTDZ)
call grib_close_file(ofileT)
call grib_close_file(ofileQ)
call grib_close_file(ofileP)


contains
subroutine correct_Tair(Tair,dtdz,DZ,Tair_new)
  real(kind=8),intent(in   ) :: Tair
  real(kind=8),intent(in   ) :: dtdz
  real(kind=8),intent(in   ) :: DZ
  real(kind=8),intent(  out) :: Tair_new
 
  ! new Temperature
  Tair_new = Tair + DZ*dtdz
end subroutine
        
subroutine correct_lnps(Tair,Tair_new,Qair,Qair_fg,lnps,DZ,RETV,RG,RD,Ps_new)
  real(kind=8),intent(in   ) :: Tair
  real(kind=8),intent(in   ) :: Tair_new
  real(kind=8),intent(in   ) :: Qair
  real(kind=8),intent(in   ) :: Qair_fg
  real(kind=8),intent(in   ) :: lnps
  real(kind=8),intent(in   ) :: DZ
  real(kind=8),intent(  out) :: Ps_new
  real(kind=8)               :: Tv
  real(kind=8)               :: RETV
  real(kind=8)               :: RG
  real(kind=8)               :: RD
  real(kind=8)               :: Tave
  real(kind=8)               :: Qave

  ! new log surface pressure
  Tave = (Tair+Tair_new)*0.5d0
  Qave = (Qair+Qair_fg)*0.5d0
  Tv = Tave * (1.d0 + RETV * Qave)
  Ps_new = exp(lnps - RG*DZ/RD/Tv)
end subroutine

subroutine correct_Qair(Tair_new,Ps_new,constA_lowml,constB_lowml,RETV,PRH,Qair_new)
  real(kind=8),intent(in   ) :: Tair_new
  real(kind=8),intent(in   ) :: Ps_new
  real(kind=8),intent(in   ) :: constA_lowml
  real(kind=8),intent(in   ) :: constB_lowml
  real(kind=8),intent(in   ) :: RETV
  real(kind=8),intent(in   ) :: PRH
  real(kind=8),intent(  out) :: Qair_new
  real(kind=8)               :: PES

  ! Pressure at lowest full model level 
  Pair_new = (constA_lowml+constB_lowml*Ps_new + Ps_new) / 2.d0

  ! saturation vapor pressure with Teten's formula : eq(7.5) for IFS Document Cy38R1
  call cal_PES(Tair_new,PES)

  ! new specific humidity assuming constant RH
  Qair_new = PRH*PES / (Pair_new*(RETV+1.0d0)-PRH*PES*RETV)
end subroutine

subroutine cal_PES(Tair,PES)
  real(kind=8),intent(in   ) :: Tair
  real(kind=8),intent(  out) :: PES
  real(kind=8)               :: FOEALFA
  real(kind=8)               :: FOEEWM
  real(kind=8)               :: RTT
  real(kind=8)               :: RTWAT
  real(kind=8)               :: RTICE
  real(kind=8)               :: RKBOL
  real(kind=8)               :: RNAVO
  real(kind=8)               :: R
  real(kind=8)               :: RMD
  real(kind=8)               :: RMV
  real(kind=8)               :: RD
  real(kind=8)               :: RV
  real(kind=8)               :: R2ES
  real(kind=8)               :: R3LES
  real(kind=8)               :: R3IES
  real(kind=8)               :: R4LES
  real(kind=8)               :: R4IES

  RTT=273.16d0          ! Temperature of water fusion at normal pressure
  RTWAT=RTT             ! Criteria for temperature with respect to liquid water
  RTICE=RTT-23.d0       ! Criteria for temperature with respect to ice
  RKBOL=1.380658E-23    ! Boltzmann constant (k) [J/K]
  RNAVO=6.0221367E+23   ! Avogadro constant (N_A) [1/mol]
  R=RNAVO*RKBOL         ! Perfect gas constant (= 8.314511211948600)
  RMD=28.9644d0         ! Dry air mass
  RMV=18.0153d0         ! Vapour  mass
  RD=1000.d0*R/RMD      ! Dry air cst. (= 287.0596736665907 J/kg/K)
  RV=1000.d0*R/RMV      ! Vapour  cst. (= 461.5249933083879 J/kg/K)
  R2ES=611.21d0*RD/RV   ! Constant in Teten's formula: eq(7.5)
  R3LES=17.502d0        ! Mixing Ratio Over Liquid Water in Teten's formula: eq(7.5)
  R3IES=22.587d0        ! Mixing Ratio Over Ice in Teten's formula: eq(7.5)
  R4LES=32.19d0         ! Constant in Teten's formula for liquid water: eq(7.5)
  R4IES=-0.7d0          ! Constant in Teten's formula for ice: eq(7.5)
  RETV=RV/RD-1.0d0      ! Rv/Rd-1 (= 0.608)

  ! ratio between water and ice : Alpha in eq(7.73) for IFS Document Cy38R1
  FOEALFA = MIN(1.0d0,((MAX(RTICE,MIN(RTWAT,Tair))-RTICE)&
            &  /(RTWAT-RTICE))**2.d0)
  ! saturation vapor pressure with Teten's formula : eq(7.5) for IFS Document Cy38R1
  FOEEWM = R2ES *&
           &(FOEALFA*EXP(R3LES*(Tair-RTT)/(Tair-R4LES))+&
           &(1.0d0-FOEALFA)*EXP(R3IES*(Tair-RTT)/(Tair-R4IES)))
  PES = (RETV+1.0d0)*FOEEWM
end subroutine
 
end program 



