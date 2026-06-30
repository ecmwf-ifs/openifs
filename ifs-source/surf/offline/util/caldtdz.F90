! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!
program caldtdz
use grib_api
implicit none
integer                               :: igribQ
integer                               :: igribT
integer                               :: igribP
integer                               :: ifileQ
integer                               :: ifileT
integer                               :: ifileP
integer                               :: ofileT
integer                               :: ofileQ
integer                               :: ofileP
integer                               :: ofileDTDZ
integer                               :: iret
integer                               :: iretQ
integer                               :: iretT
integer                               :: iretP
integer                               :: countT
integer                               :: countQ
integer                               :: idxT
integer                               :: idxQ
integer                               :: idxP
integer                               :: i
integer                               :: j
integer                               :: k
integer                               :: idd
integer                               :: itt
integer                               :: iff
integer                               :: nb_pv
integer(4)                            :: numberOfPoints
integer(4)                            :: dataDateSize
integer(4)                            :: dataTimeSize
integer(4)                            :: stepSize
integer(4)                            :: levelSizeT
integer(4)                            :: levelSizeQ
integer(4)                            :: nlev
integer(4)                            :: nlev1
real(8)   ,dimension(:  ),allocatable :: pv
real(8)   ,dimension(:  ),allocatable :: constA
real(8)   ,dimension(:  ),allocatable :: constB
real(8)   ,dimension(:  ),allocatable :: ph          ! Half level pressure [Pa]
real(8)   ,dimension(:  ),allocatable :: dp          ! Thickness of half level pressure [Pa]
real(8)   ,dimension(:  ),allocatable :: zh          ! Half level height [m]
real(8)   ,dimension(:  ),allocatable :: zf          ! Full level height [m]
real(8)                               :: Ps          ! Surface Pressure before modification
real(8)                               :: RG
real(8)                               :: RKBOL
real(8)                               :: RNAVO
real(8)                               :: R
real(8)                               :: RMD
real(8)                               :: RMV
real(8)                               :: RD
real(8)                               :: RV
real(8)                               :: RETV
integer(4),dimension(:  ),allocatable :: level
integer(4),dimension(:  ),allocatable :: dataDates
integer(4),dimension(:  ),allocatable :: dataTimes
integer(4),dimension(:  ),allocatable :: step
real(8)   ,dimension(:  ),allocatable :: lnps        ! Log Surface Pressure [Pa]
real(8)   ,dimension(:  ),allocatable :: dtdz        ! Temperature Lapse Rate [K/m]
real(8)   ,dimension(:  ),allocatable :: DZ          ! Orography difference [m]
real(8)   ,dimension(:,:),allocatable :: Tairs       ! Temperature on Lowest Full Model Level
real(8)   ,dimension(:,:),allocatable :: Qairs       ! Specific Humidity on Lowest Full Model Level


!--- Prepare Constants ---
RG=9.80665d0          ! Gravitational Acceleration [m/s2]
RKBOL=1.380658E-23    ! Boltzmann constant (k) [J/K]
RNAVO=6.0221367E+23   ! Avogadro constant (N_A) [1/mol]
R=RNAVO*RKBOL         ! Perfect gas constant (= 8.314511211948600)
RMD=28.9644d0         ! Dry air mass
RMV=18.0153d0         ! Vapour  mass
RD=1000.d0*R/RMD      ! Dry air cst. (= 287.0596736665907 J/kg/K)
RV=1000.d0*R/RMV      ! Vapour  cst. (= 461.5249933083879 J/kg/K)
RETV=RV/RD-1.0d0      ! Rv/Rd-1 (= 0.608)



!===============================================================
!  Read variables from input grib file
!===============================================================
!------------------------------------------
! Get 'numberOfPoints' from input grib file
!------------------------------------------
call grib_open_file(ifileT, 'INPUTGRIB_T','R')
call grib_new_from_file(ifileT,igribT,iretT)
call grib_get(igribT,'numberOfPoints',numberOfPoints)

! Get number of vertical layers and constant A and B 
call grib_get_size(igribT,'pv',nb_pv)
allocate(pv(nb_pv))
call grib_get(igribT,'pv',pv)
nlev1=nb_pv/2
nlev=nlev1-1
write(6,*) "nb_pv, nlev1, nlev = ", nb_pv, nlev1,nlev
allocate(constA(nlev1))
allocate(constB(nlev1))
constA(1:nlev1)=pv(1:nlev1)
constB(1:nlev1)=pv(nlev1+1:nb_pv)
!write(6,*) "constA = ", constA
!write(6,*) "constB = ", constB

call grib_release(igribT)
call grib_close_file(ifileT)

!------------------------------------------
! Create grib keys and get basic information (Temperature)
!------------------------------------------
call grib_index_create(idxT,'INPUTGRIB_T','level,dataDate,dataTime,step')
call grib_index_get_size(idxT,'level',levelSizeT)
allocate(level(levelSizeT))
call grib_index_get(idxT,'level',level)
call grib_index_get_size(idxT,'dataDate',dataDateSize)
allocate(dataDates(dataDateSize))
call grib_index_get(idxT,'dataDate',dataDates)
call grib_index_get_size(idxT,'dataTime',dataTimeSize)
allocate(dataTimes(dataTimeSize))
call grib_index_get(idxT,'dataTime',dataTimes)
call grib_index_get_size(idxT,'step',stepSize)
allocate(step(stepSize))
call grib_index_get(idxT,'step',step)

!------------------------------------------
! Create grib keys (Specific Humidity)
!------------------------------------------
call grib_index_create(idxQ,'INPUTGRIB_Q','level,dataDate,dataTime,step')
call grib_index_get_size(idxQ,'level',levelSizeQ)
if (levelSizeT /= levelSizeQ) then
  write(6,*) "ERROR! : levelSizeT /= levelSizeQ"
  stop
endif

!------------------------------------------
! Create grib keys (Log Surface Pressure)
!------------------------------------------
call grib_index_create(idxP,'INPUTGRIB_LNPS','dataDate,dataTime,step')

!------------------------------------------
! Open grib file for output (Temperature Lapse Rate)
!------------------------------------------
call grib_open_file(ofileDTDZ, 'OUTPUTGRIB_DTDZ','W')


!------------------------------------------
! Allocate variables
!------------------------------------------
allocate(ph(nlev1))
allocate(dp(nlev1))
allocate(zh(nlev1))
allocate(zf(nlev1))

allocate(Tairs(numberOfPoints,levelSizeT))
allocate(Qairs(numberOfPoints,levelSizeQ))
allocate(lnps(numberOfPoints))
allocate(dtdz(numberOfPoints))

!------------------------------------------
! Read variables from input grib file 
!------------------------------------------
do idd=1,dataDateSize
  call grib_index_select(idxT,'dataDate',dataDates(idd))
  call grib_index_select(idxQ,'dataDate',dataDates(idd))
  call grib_index_select(idxP,'dataDate',dataDates(idd))
  do itt=1,dataTimeSize
    call grib_index_select(idxT,'dataTime',dataTimes(itt))
    call grib_index_select(idxQ,'dataTime',dataTimes(itt))
    call grib_index_select(idxP,'dataTime',dataTimes(itt))
    do iff=1,stepSize
      call grib_index_select(idxT,'step',step(iff))
      call grib_index_select(idxQ,'step',step(iff))
      call grib_index_select(idxP,'step',step(iff))

      Tairs(:,:) = 0.d0
      Qairs(:,:) = 0.d0

      do k=1,levelSizeT
        call grib_index_select(idxT,'level',level(k))
        call grib_index_select(idxQ,'level',level(k))

        call grib_new_from_index(idxT,igribT, iretT)
        call grib_new_from_index(idxQ,igribQ, iretQ)

        countT=0
        countQ=0
        do while (iretT /= GRIB_END_OF_INDEX)
          countT=countT+1
          if (iretQ /= GRIB_END_OF_INDEX) then
            countQ=countQ+1

            !------------------------------------------
            ! Read 3D Temperature
            !------------------------------------------
            call grib_get(igribT,'values',Tairs(:,k))
            call grib_release(igribT)
            call grib_new_from_index(idxT,igribT, iretT)
          
            !------------------------------------------
            ! Read 3D Specific Humidity
            !------------------------------------------
            call grib_get(igribQ,'values',Qairs(:,k))
            call grib_release(igribQ)
            call grib_new_from_index(idxQ,igribQ, iretQ)
          endif
        end do
        call grib_release(igribT)
        call grib_release(igribQ)
      end do 

      if (countT == 1 .and. countQ == 1) then
        write(6,'(A,i8,A,i4,A,i4,A,i4)') "Temperature Lapse Rate : date = ", dataDates(idd),",  time = ", dataTimes(itt), &
          & ",  step = ", step(iff)

        !------------------------------------------
        ! Read Log Surface Pressure
        !------------------------------------------
        call grib_new_from_index(idxP,igribP, iretP)
        call grib_get(igribP,'values',lnps(:))
        !call grib_release(igribP)

        !===============================================================
        !  Calculation of temperature lapse rate
        !===============================================================
        do i=1,numberOfPoints
     !if(i==10000) then
          Ps = exp(lnps(i))

          ph(:) = 0.d0
          zh(:) = 0.d0
          if (level(levelSizeT) == nlev) then
            ph(levelSizeT+1) = constA(level(levelSizeT)+1) + constB(level(levelSizeT)+1) * Ps
            zh(levelSizeT+1) = 0.0d0 
            do k = levelSizeT+1, 2, -1
              ph(k-1) = constA(level(k-1)) + constB(level(k-1)) * Ps
              zh(k-1) = zh(k) + RD*Tairs(i,k-1)*(1.d0+RETV*Qairs(i,k-1))*LOG(ph(k)/ph(k-1))/RG
            enddo
            do k = levelSizeT, 1, -1
              dp(k) = ph(k+1) - ph(k)
              zf(k) = zh(k+1) &
                    & + (1.d0 - ph(k)*LOG(ph(k+1)/ph(k))/dp(k)) &
                        *RD*Tairs(i,k)*(1.d0 + RETV*Qairs(i,k))/RG
            enddo
          else
            write(6,*) "ERROR! : level(levelSizeT) /= nlev"
            stop
          endif

          ! Temperature Lapse Rate [K/m]
          dtdz(i) = (Tairs(i,1)-Tairs(i,levelSizeT)) / (zf(1)-zf(levelSizeT))
          dtdz(i) = min(dtdz(i),0.d0)

     !  write(6,*)
     !  do k = 1,levelSizeT
     !    write(6,'(A,i4,3x,f10.5,3x,f15.5)') "k-1/2, zh, ph    = ",k  , zh(k)  , ph(k)
     !    write(6,'(A,i4,3x,f10.5,3x,f15.5)') "k    , zf, Tair  = ",k  , zf(k)  , Tairs(i,k)
     !    write(6,'(A,i4,3x,f10.5,3x,f15.5)') "k+1/2, zh, ph    = ",k+1, zh(k+1), ph(k+1)
     !    write(6,*)
     !  enddo
     !  write(6,*) "Temperature Lapse Rate =", dtdz(i)
     !endif
        end do ! loop for numberOfPoints


        !==============================================================
        !  Outout of modified variables using INPUTGRIB_DZ information
        !==============================================================
        call grib_set(igribP,"values",dtdz)
        call grib_write(igribP,ofileDTDZ,iret)
        call grib_release(igribP)
      endif
    end do ! loop for stepSize
  end do   ! loop for dataTimeSize
end do     ! loop for dataDateSize


!------------------------------------------
! Post process
!------------------------------------------
call grib_close_file(ofileDTDZ)
deallocate(level)
deallocate(Tairs)
deallocate(Qairs)
deallocate(lnps)
deallocate(dtdz)
deallocate(constA)
deallocate(constB)
deallocate(ph)
deallocate(dp)
deallocate(zh)
deallocate(zf)
call grib_index_release(idxT)
call grib_index_release(idxQ)
call grib_index_release(idxP)
    

end program 


