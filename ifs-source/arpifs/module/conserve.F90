! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

! module for mass and water conservation, mainly for long runs (climate)
! original by A Voldoire, routines conservini/conserv

module conserve
	use parkind1,only: jpim,jprb
	use yomlun,only: nulout,nulnam
	use yomcst,only: rg,rcpd,rcpv,rcw,rcs
	use yomct0,only: lcorwat
	use disgrid_mod,only: disgrid_send,disgrid_recv
	use diwrgrid_mod,only: diwrgrid_send,diwrgrid_recv
	use mpl_module,only: mpl_broadcast
	use yomgppb,only: gparbuf
	use yomcfu,only: tcfu
	use yommp0,only: myproc,my_region_ns,my_region_ew
	use type_geometry,only: geometry
	use yomdim,only: tdim
	use yomdimv,only: tdimv
	use yomgem,only: tgem
	use yomvert,only: tvab
	use yommp,only: tmp
	use yomrip,only: trip
	use yom_ygfl,only: type_gfld
	use yomgfl,only: tgfl
	use spectral_fields_data,only: spectral_field
	use surface_fields_mix,only: tsurf
	use yomphy,only: tphy
	use yomparar,only: tparar
	use yomhook,only: lhook,dr_hook,jphook
	use blocks

	implicit none

	private
	public conservini,conserv

	integer(kind=jpim),parameter :: proc0=1
	real(kind=jprb) :: vtotwat,vdrymass

#include "inv_trans.h"
#include "dir_trans.h"
contains
	subroutine conservini(geo,ygfl,phy,sp,gfl,ldcorm)
		type(geometry),intent(in) :: geo
		type(type_gfld),intent(in) :: ygfl
		type(tphy),intent(in) :: phy
		type(spectral_field),intent(in) :: sp
		real(kind=jprb),intent(in) :: gfl(:,:,:,:)
		logical,intent(in) :: ldcorm

		integer(kind=jpim) :: ifld,irw,igm,ipsf,ipdryf,ipvf,iener,liquid,ice,rain,snow
		real(kind=jprb),dimension(:),allocatable :: zweight,zphis
		real(kind=jprb),dimension(:,:),allocatable :: zsendg,zh,zrr,zrrcv
		real(kind=jprb),dimension(:,:,:),allocatable :: zpq
		real(kind=jprb),dimension(:,:),allocatable,target :: zsend
		real(kind=jprb),pointer :: zrw(:),zpsf(:),zpvf(:),zener(:)
		real(kind=jphook) :: zhook_handle

		if (lhook) call dr_hook('conservini',0,zhook_handle)
		associate(dimg=>geo%yrdim,dimv=>geo%yrdimv,gem=>geo%yrgem,vab=>geo%yrvab,&
			mp=>geo%yrmp,gsgeom=>geo%yrgsgeom_nb,csgleg=>geo%yrcsgleg,lcondwt=>phy%lcondwt,&
			lgpcmt=>phy%lgpcmt)

		associate(nproma=>dimg%nproma,ngpblks=>dimg%ngpblks,nflevg=>dimv%nflevg,&
			nspec2=>dimg%nspec2,ngptot=>gem%ngptot,ngptotg=>gem%ngptotg,rw=>csgleg%rw,&
			gm=>gsgeom%gm,nbsetsp=>mp%nbsetsp,nbsetlev=>mp%nbsetlev)

		irw = 1
		igm = 2
		ifld = 2

		if (ldcorm) then
			ipsf = getnaddindex(ifld)
			ipdryf = getnaddindex(ifld)
		end if

		if (lcorwat) ipvf = getnaddindex(ifld)

		allocate(zsend(ngptot,ifld))

		zrw => zsend(:,irw)
		call setprocfield(mp,ngptot,rw,zrw)
		zsend(:,igm) = gm(:)

		if (ldcorm) then
			zpsf => zsend(:,ipsf)
		else
			allocate(zpsf(ngptot))
		end if

		if (lcorwat) then
			zpvf => zsend(:,ipvf)
		else
			allocate(zpvf(ngptot))
		end if

		call sptrans_spsc2(nproma,ngpblks,nspec2,ngptot,nbsetsp,sp%sp,zpsf)
		zpsf(:) = exp(zpsf(:))

		allocate(zpq(nproma,nflevg,ngpblks),zh(ngptot,nflevg))
		call inv_trans(pspscalar=sp%q,kproma=nproma,pgp=zpq,kvsetsc=nbsetlev)
		call block2contiguous(nproma,ngpblks,ngptot,nflevg,zpq,zh)
		deallocate(zpq)

		zpvf(:) = 0
		call blockvsum(dimv,vab,zpsf,zh,zpvf)

		if (ldcorm) zsend(:,ipdryf) = zpsf(:)-zpvf(:)

		if (lcorwat) then
			allocate(zrr(ngptot,nflevg),zrrcv(ngptot,nflevg))

			if (lcondwt) then
				call condensatesind(ygfl,liquid,ice,rain,snow)
				call condpre(dimg,gem,gfl,liquid,ice,rain,snow,zrr)
				call blockvsum(dimv,vab,zpsf,zrr,zpvf)
			end if

			if (lgpcmt) then
				call convcondensatesind(ygfl,liquid,ice,rain,snow)
				call condpre(dimg,gem,gfl,liquid,ice,rain,snow,zrrcv)
				call blockvsum(dimv,vab,zpsf,zrrcv,zpvf)
			end if

			deallocate(zrr,zrrcv)
		end if

		deallocate(zh)

		if (myproc == proc0) then
			allocate(zsendg(ngptotg,ifld),zweight(ngptotg))

			call diwrgrid_recv(geo,ifld,zsend,1,zsendg)

			call gpweights(dimg,gem,zsendg(:,irw),zsendg(:,igm),zweight)

			if (ldcorm) then
				vdrymass = sum(zsendg(:,ipdryf)*zweight(:))
				write(nulout,"('conservini: total and dry masses',2(x,f10.3))")&
					sum(zsendg(:,ipsf)*zweight(:)),vdrymass
			end if

			if (lcorwat) then
				vtotwat = sum(zsendg(:,ipvf)*zweight(:))/rg
				write(nulout,"('conservini: total atm. water mass ',f8.3)") vtotwat
			end if

			deallocate(zweight,zsendg)
		else
			call diwrgrid_send(gem,1,ifld,zsend,1)
		end if

		deallocate(zsend)
		if (.not.ldcorm) deallocate(zpsf)
		if (.not.lcorwat) deallocate(zpvf)
		end associate
		end associate
		if (lhook) call dr_hook('conservini',1,zhook_handle)
	end subroutine

	subroutine conserv(geo,rip,ygfl,phy,parar,sp,surf,ydcfu,gfl,ldcorm)
		type(geometry),intent(in) :: geo
		type(trip),intent(in) :: rip
		type(type_gfld),intent(in) :: ygfl
		type(tphy),intent(in) :: phy
		type(tparar),intent(in) :: parar
		type(spectral_field),intent(in) :: sp
		type(tsurf),intent(in) :: surf
		type(tcfu),intent(inout) :: ydcfu
		real(kind=jprb),intent(in) :: gfl(:,:,:,:)
		logical,intent(in) :: ldcorm

		integer(kind=jpim) :: ifld,jfl,jp,jf,ibl,joff,k1,np,irw,igm,ipsf,ipdryf,ievasf,&
			iprpar,ipvf,ispvf,iprparl,iprparn,iphis,iflxt,ifts,ifsols,liquid,&
			ice,rain,snow
		real(kind=jprb)   :: zapcorm,zdrymass,ztotwat,zda,zdb,zesuf,zadd,zpr,zratio(1)
        real(kind=jphook) :: zhook_handle
		real(kind=jprb),dimension(:),allocatable :: zweight
		real(kind=jprb),dimension(:,:),allocatable :: zsendg,zh,zrr,zrrcv
		real(kind=jprb),dimension(:,:),allocatable,target :: zsend
		real(kind=jprb),dimension(:,:,:),allocatable :: zpq
		real(kind=jprb),pointer :: zrw(:),zpsf(:),zpvf(:),zspvf(:),zpdryf(:),zevasf(:),&
			zprpar(:),zprparl(:),zprparn(:)

		if (lhook) call dr_hook('conserv',0,zhook_handle)
		associate(dimg=>geo%yrdim,dimv=>geo%yrdimv,gem=>geo%yrgem,vab=>geo%yrvab,&
			mp=>geo%yrmp,gsgeom=>geo%yrgsgeom_nb,csgleg=>geo%yrcsgleg,lcondwt=>phy%lcondwt,&
			lgpcmt=>phy%lgpcmt)

		associate(nproma=>dimg%nproma,ngpblks=>dimg%ngpblks,nflevg=>dimv%nflevg,&
			nspec2=>dimg%nspec2,ngptot=>gem%ngptot,ngptotg=>gem%ngptotg,nbsetsp=>mp%nbsetsp,&
			nbsetlev=>mp%nbsetlev,tstep=>rip%tstep,rw=>csgleg%rw,gm=>gsgeom%gm,&
			sd_vh=>surf%sd_vh,sp_rr=>surf%sp_rr,ysd_vh=>surf%ysd_vh,ysp_rr=>surf%ysp_rr,&
			ycfupt=>ydcfu%ycfupt,gfubuf=>ydcfu%gfubuf)

		irw = 1
		igm = 2
		ipsf = 3
		ispvf = 4
		ifld = 4

		if (ldcorm) ipdryf = getnaddindex(ifld)

		if (lcorwat) then
			ievasf = getnaddindex(ifld)
			ipvf = getnaddindex(ifld)
			iprpar = getnaddindex(ifld)
			iprparl = getnaddindex(ifld)
			iprparn = getnaddindex(ifld)
		end if

		allocate(zsend(ngptot,ifld))

		zrw => zsend(:,irw)
		zsend(:,igm) = gm(:)
		zpsf => zsend(:,ipsf)
		zspvf => zsend(:,ispvf)

		if (ldcorm) zpdryf => zsend(:,ipdryf)

		if (lcorwat) then
			zpvf => zsend(:,ipvf)
			zevasf => zsend(:,ievasf)
			zprpar => zsend(:,iprpar)
			zprparl => zsend(:,iprparl)
			zprparn => zsend(:,iprparn)
		else
			allocate(zpvf(ngptot))
		end if

		call setprocfield(mp,ngptot,rw,zrw)

		call sptrans_spsc2(nproma,ngpblks,nspec2,ngptot,nbsetsp,sp%sp,zpsf)
		zpsf(:) = exp(zpsf(:))

		allocate(zpq(nproma,nflevg,ngpblks),zh(ngptot,nflevg))
		call inv_trans(pspscalar=sp%q,kproma=nproma,pgp=zpq,kvsetsc=nbsetlev)
		call block2contiguous(nproma,ngpblks,ngptot,nflevg,zpq,zh)
		deallocate(zpq)

		zpvf(:) = 0
		zspvf(:) = 0
		call blockvsum(dimv,vab,zpsf,zh,zpvf,zspvf)

		if (ldcorm) zpdryf(:) = zpsf(:)-zpvf(:)

		if (lcorwat) then
			allocate(zrr(ngptot,nflevg),zrrcv(ngptot,nflevg))

			if (lcondwt) then
				call condensatesind(ygfl,liquid,ice,rain,snow)
				call condpre(dimg,gem,gfl,liquid,ice,rain,snow,zrr)
			end if

			if (lgpcmt) then
				call convcondensatesind(ygfl,liquid,ice,rain,snow)
				call condpre(dimg,gem,gfl,liquid,ice,rain,snow,zrrcv)
			end if
		end if

		if (lcorwat) then
			if (lcondwt) call blockvsum(dimv,vab,zpsf,zrr,zpvf,zspvf)
			if (lgpcmt) call blockvsum(dimv,vab,zpsf,zrrcv,zpvf,zspvf)

			call block2contiguous1(nproma,ngpblks,ngptot,ysd_vh%yeva%mp,sd_vh,zevasf)
		end if

		deallocate(zh)

		if (lcorwat) then
			deallocate(zrr,zrrcv)

			if (parar%ngpar == 0) call abor1("Error: no gp parameter")

			call block2contiguous1(nproma,ngpblks,ngptot,parar%mrain,gparbuf,zprparl)
			call block2contiguous1(nproma,ngpblks,ngptot,parar%msnow,gparbuf,zprparn)
		end if

		if (lcorwat) zprpar(:) = zprparl(:) + zprparn(:)

		if (myproc == proc0) then
			allocate(zsendg(ngptotg,ifld),zweight(ngptotg))

			!$! call disgrid_recv(geo,proc0,ifld,zsendg,1)
			call diwrgrid_recv(geo,ifld,zsend,1,zsendg)

			call gpweights(dimg,gem,zsendg(:,irw),zsendg(:,igm),zweight)

			if (ldcorm.or.lcorwat) zadd = sum(zsendg(:,ispvf)*zweight(:))

			if (ldcorm) then
				zdrymass = sum(zsendg(:,ipdryf)*zweight(:))

				zapcorm = (vdrymass-zdrymass)/(1-zadd)

				write(nulout,"('conserv: current and former total dry masses',2(x,f10.3))")&
					zdrymass,vdrymass
				write(nulout,"('conserv: correction on total dry mass ',f8.3)") zapcorm

				zsendg(:,ipsf) = zsendg(:,ipsf)+zapcorm

				ifld = ifld+1
				call disgrid_send(geo,1,zsendg(:,ipsf),ifld,zpsf)
			else
				zapcorm = 0
			end if

			if (lcorwat) then
				zesuf = sum(zsendg(:,ievasf)*zweight(:))*tstep
				zpr = sum(zsendg(:,iprpar)*zweight(:))*tstep
				ztotwat = (sum(zsendg(:,ipvf)*zweight(:))+zapcorm*zadd)/rg

				if (zpr > epsilon(1._jprb)) then
					zratio = (vtotwat-ztotwat+zesuf)/zpr
				else
					zratio = 1
				end if

				write(nulout,"('conserv: correction ratio on water ',f7.4)") zratio

				vtotwat = ztotwat

				ifld = ifld+1
				call mpl_broadcast(zratio,ifld,proc0,cdstring="conserv-wat:")
			end if

			deallocate(zsendg)
		else
			call diwrgrid_send(gem,proc0,ifld,zsend,1)

			if (ldcorm) then
				ifld = ifld+1
				call disgrid_recv(geo,proc0,1,zpsf,ifld)
			end if

			if (lcorwat) then
				ifld = ifld+1
				call mpl_broadcast(zratio,ifld,proc0,cdstring="conserv-wat:")
			end if
		end if

		if (ldcorm) then
			zpsf(:) = log(zpsf(:))
			call gptrans_spsc2(nproma,ngpblks,nspec2,ngptot,nbsetsp,zpsf,sp%sp)
		end if

		deallocate(zsend)

		if (lcorwat) then
			gparbuf(:,parar%mrain,:) = gparbuf(:,parar%mrain,:)*zratio(1)
			gparbuf(:,parar%msnow,:) = gparbuf(:,parar%msnow,:)*zratio(1)

			do ibl=1,ngpblks
				np = blocksize(ibl,nproma,ngptot)

				if (ydcfu%nfdcfu > 0) then
					k1 = 1
					do jf=1,ydcfu%ntypcfu
						do jfl=k1,k1+ydcfu%type_cfu(jf)%ilev
							if (k1 == ycfupt%mfplcl) then
								gfubuf(1:np,jfl,ibl) = gfubuf(1:np,jfl,ibl)+&
									sd_vh(1:np,ysd_vh%ypcl%mp,ibl)*tstep*zratio(1)
							else if (k1 == ycfupt%mfplsl) then
								gfubuf(1:np,jfl,ibl) = gfubuf(1:np,jfl,ibl)+&
									sd_vh(1:np,ysd_vh%ypsl%mp,ibl)*tstep*zratio(1)
							else if (k1 == ycfupt%mfplcn) then
								gfubuf(1:np,jfl,ibl) = gfubuf(1:np,jfl,ibl)+&
									sd_vh(1:np,ysd_vh%ypcn%mp,ibl)*tstep*zratio(1)
							else if (k1 == ycfupt%mfplsn) then
								gfubuf(1:np,jfl,ibl) = gfubuf(1:np,jfl,ibl)+&
									sd_vh(1:np,ysd_vh%ypsn%mp,ibl)*tstep*zratio(1)
							end if
						end do

						k1 = k1+ydcfu%type_cfu(jf)%ilev
					end do
				end if
			end do
		else
			deallocate(zpvf)
		end if
		end associate
		end associate
		if (lhook) call dr_hook('conserv',1,zhook_handle)
	end subroutine

	integer function getnaddindex(index)
		integer,intent(inout) :: index

		index = index + 1
		getnaddindex = index
	end function

	subroutine setprocfield(mp,ngptot,xgeo,z)
		type(tmp),intent(in) :: mp
		integer(kind=jpim),intent(in) :: ngptot
		real(kind=jprb),intent(in) :: xgeo(:)
		real(kind=jprb),intent(out) :: z(:)

		integer(kind=jpim) :: k1,k2,jgl

		associate(nonl=>mp%nonl,nlstlat=>mp%nlstlat,nptrfloff=>mp%nptrfloff,&
			nfrstloff=>mp%nfrstloff,mylats=>mp%mylats)

		k1 = 1
		do jgl=1,nlstlat(my_region_ns)-nfrstloff
			k2 = k1 + nonl(nptrfloff+jgl,my_region_ew) - 1
			z(k1:k2) = xgeo(mylats(jgl))
			k1 = k2 + 1
		end do
		end associate
	end subroutine

	subroutine sptrans_spsc2(nproma,ngpblks,nspec2,ngptot,nbsetsp,zsps,zpsf)
		integer(kind=jpim),intent(in) :: nproma,ngpblks,nspec2,ngptot,nbsetsp
		real(kind=jprb),intent(in) :: zsps(1,nspec2)
		real(kind=jprb),intent(out) :: zpsf(:)

		integer(kind=jpim) :: ibl,joff,np,ivsetsc(1)
		real(kind=jprb) :: pgp2(nproma,1,ngpblks)

		ivsetsc(1) = nbsetsp
		call inv_trans(pspsc2=zsps,kproma=nproma,pgp2=pgp2,kvsetsc2=ivsetsc)

		call block2contiguous1(nproma,ngpblks,ngptot,1,pgp2,zpsf)
	end subroutine

	subroutine sptrans_spvordiv(dimg,dimv,sp,nbsetlev,zu,zv)
		type(tdim),intent(in) :: dimg
		type(tdimv),intent(in) :: dimv
		type(spectral_field),intent(in) :: sp
		integer(kind=jpim),intent(in) :: nbsetlev(:)
		real(kind=jprb),intent(out) :: zu(:,:,:),zv(:,:,:)

		real(kind=jprb),dimension(:,:,:),allocatable :: zuv

		allocate(zuv(dimg%nproma,2*dimv%nflevg,dimg%ngpblks))

		call inv_trans(pspvor=sp%vor,pspdiv=sp%div,kproma=dimg%nproma,pgp=zuv,&
			kvsetuv=nbsetlev)

		zu(:,:,:) = zuv(:,1:dimv%nflevg,:)
		zv(:,:,:) = zuv(:,dimv%nflevg+1:2*dimv%nflevg,:)

		deallocate(zuv)
	end subroutine

	subroutine gptrans_spsc2(nproma,ngpblks,nspec2,ngptot,nbsetsp,zpsf,zsps)
		integer(kind=jpim),intent(in) :: nproma,ngpblks,nspec2,ngptot,nbsetsp
		real(kind=jprb),intent(in) :: zpsf(:)
		real(kind=jprb),intent(out) :: zsps(1,nspec2)

		integer(kind=jpim) :: ibl,joff,np,ivsetsc(1)
		real(kind=jprb) :: pgp2(nproma,1,ngpblks)

		call contiguous2block1(nproma,ngpblks,ngptot,1,zpsf,pgp2)

		ivsetsc(1) = nbsetsp
		call dir_trans(pspsc2=zsps,kproma=nproma,pgp2=pgp2,kvsetsc2=ivsetsc)
	end subroutine

	subroutine energy(lcondwt,lgpcmt,dimg,dimv,gem,vab,sp,nbsetlev,zh,zrr,zrrcv,zpsf,&
		zener,zsener)
		logical,intent(in) :: lcondwt,lgpcmt
		type(tdim),intent(in) :: dimg
		type(tdimv),intent(in) :: dimv
		type(tgem),intent(in) :: gem
		type(tvab),intent(in) :: vab
		type(spectral_field),intent(in) :: sp
		integer(kind=jpim),intent(in) :: nbsetlev(:)
		real(kind=jprb),intent(in) :: zh(:,:),zrr(:,:),zrrcv(:,:),zpsf(:)
		real(kind=jprb),intent(out) :: zener(:)
		real(kind=jprb),optional,intent(out) :: zsener(:)

		integer(kind=jpim) :: ibl,joff,np
		real(kind=jprb),allocatable :: zt(:,:,:),zu(:,:,:),zv(:,:,:),ze(:,:)

		associate(nproma=>dimg%nproma,ngpblks=>dimg%ngpblks,nflevg=>dimv%nflevg,&
			ngptot=>gem%ngptot)

		allocate(zt(nproma,nflevg,ngpblks))
		allocate(zu(nproma,nflevg,ngpblks),zv(nproma,nflevg,ngpblks))

		call sptrans_spvordiv(dimg,dimv,sp,nbsetlev,zu,zv)
		call inv_trans(pspscalar=sp%t,kproma=nproma,pgp=zt,kvsetsc=nbsetlev)

		allocate(ze(ngptot,nflevg))

		!$! optimisation: minimisation of loops
		!$! ze is Cp in this loop
		if (lcondwt.and.lgpcmt) then
			ze(:,:) = rcpd+(rcpv-rcpd)*zh(:,:)+(rcw-rcpd)*zrr(:,:)+(rcs-rcpd)*zrrcv(:,:)
		else if (lcondwt) then
			ze(:,:) = rcpd+(rcpv-rcpd)*zh(:,:)+(rcw-rcpd)*zrr(:,:)
		else if (lgpcmt) then
			ze(:,:) = rcpd+(rcpv-rcpd)*zh(:,:)+(rcs-rcpd)*zrrcv(:,:)
		else
			ze(:,:) = rcpd+(rcpv-rcpd)*zh(:,:)
		end if

		joff = 0

		do ibl=1,ngpblks
			np = blocksize(ibl,nproma,ngptot)

			ze(joff+1:joff+np,:) = ze(joff+1:joff+np,:)*zt(1:np,:,ibl)+(zu(1:np,:,ibl)**2+&
				zv(1:np,:,ibl)**2)/2

			joff = joff + np
		end do

		zener(:) = 0
		if (present(zsener)) zsener(:) = 0
		call blockvsum(dimv,vab,zpsf,ze,zener,zsener)

		deallocate(ze,zv,zu,zt)
		end associate
	end subroutine

	subroutine gpweights(dimg,gem,zrw,zgm,zweight)
		type(tdim),intent(in) :: dimg
		type(tgem),intent(in) :: gem
		real(kind=jprb),intent(in) :: zrw(:),zgm(:)

		integer(kind=jpim) :: jgl,k1,k2
		real(kind=jprb),dimension(:),allocatable :: zweight

		associate(ndglg=>dimg%ndglg,nloeng=>gem%nloeng)

		k1 = 1
		do jgl=1,ndglg
			k2 = k1+nloeng(jgl)-1
			zweight(k1:k2) = zrw(k1:k2)/zgm(k1:k2)**2/nloeng(jgl)
			k1 = k2 + 1
		end do
		end associate
	end subroutine

	subroutine blockvsum(dimv,vab,zpsf,zmult,z,zs)
		type(tdimv),intent(in) :: dimv
		type(tvab),intent(in) :: vab
		real(kind=jprb),intent(in) :: zpsf(:),zmult(:,:)
		real(kind=jprb),intent(out) :: z(:)
		real(kind=jprb),optional,intent(out) :: zs(:)

		integer(kind=jpim) :: jl
		real(kind=jprb) :: zda,zdb

		associate(nflevg=>dimv%nflevg,vah=>vab%vah,vbh=>vab%vbh)

		if (present(zs)) then
			do jl=1,nflevg
				zda = vah(jl)-vah(jl-1)
				zdb = vbh(jl)-vbh(jl-1)

				z(:) = z(:)+zmult(:,jl)*(zda+zdb*zpsf(:))
				zs(:) = zs(:)+zmult(:,jl)*zdb
			end do
		else
			do jl=1,nflevg
				zda = vah(jl)-vah(jl-1)
				zdb = vbh(jl)-vbh(jl-1)

				z(:) = z(:)+zmult(:,jl)*(zda+zdb*zpsf(:))
			end do
		end if
		end associate
	end subroutine

	subroutine condensatesind(ygfl,liquid,ice,rain,snow)
		type(type_gfld),intent(in) :: ygfl
		integer(kind=jpim),intent(out) :: liquid,ice,rain,snow

		liquid = ygfl%yl%mp
		ice = ygfl%yi%mp
		rain = ygfl%yr%mp
		snow = ygfl%ys%mp
	end subroutine

	subroutine convcondensatesind(ygfl,liquid,ice,rain,snow)
		type(type_gfld),intent(in) :: ygfl
		integer(kind=jpim),intent(out) :: liquid,ice,rain,snow

		liquid = ygfl%ylconv%mp
		ice = ygfl%yiconv%mp
		rain = ygfl%yrconv%mp
		snow = ygfl%ysconv%mp
	end subroutine

	subroutine condpre(dimg,gem,gfl,liquid,ice,rain,snow,z)
		type(tdim),intent(in) :: dimg
		type(tgem),intent(in) :: gem
		integer(kind=jpim),intent(in) :: liquid,ice,rain,snow
		real(kind=jprb),intent(in) :: gfl(:,:,:,:)
		real(kind=jprb),intent(out) :: z(:,:)

		integer(kind=jpim) :: ibl,joff,np

		associate(nproma=>dimg%nproma,ngpblks=>dimg%ngpblks,ngptot=>gem%ngptot)

		joff = 0

		do ibl=1,ngpblks
			np = blocksize(ibl,nproma,ngptot)

			z(joff+1:joff+np,:) = gfl(1:np,:,liquid,ibl)+gfl(1:np,:,ice,ibl)+&
				gfl(1:np,:,snow,ibl)+gfl(1:np,:,rain,ibl)

			joff = joff + np
		end do
		end associate
	end subroutine
end module
