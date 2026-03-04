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

! module for 'nproma' blocks management: utility routines
! original by H Petithomme, 2019
module blocks
	use parkind1,only: jpim,jprb

	implicit none
contains
	integer function blockcount(nproma,ngptot)
		integer,intent(in) :: nproma,ngptot

		blockcount = (ngptot-1)/nproma+1
	end function

	integer function blocksize(ib,nproma,ngptot)
		integer,intent(in) :: ib,nproma,ngptot

		blocksize = min(ngptot-(ib-1)*nproma,nproma)
	end function

	subroutine block2contiguous1(nproma,ngpblks,ngptot,ip,zb,z)
		integer(kind=jpim),intent(in) :: nproma,ngpblks,ngptot,ip
		real(kind=jprb),intent(in) :: zb(:,:,:)
		real(kind=jprb),intent(out) :: z(:)

		integer(kind=jpim) :: ibl,joff,np

		joff = 0

		do ibl=1,ngpblks
			np = blocksize(ibl,nproma,ngptot)

			z(joff+1:joff+np) = zb(1:np,ip,ibl)

			joff = joff + np
		end do
	end subroutine

	subroutine block2contiguous1n(nproma,ngpblks,ngptot,ip,nlevel,zb,z)
		integer(kind=jpim),intent(in) :: nproma,ngpblks,ngptot,ip,nlevel
		real(kind=jprb),intent(in) :: zb(:,:,:,:)
		real(kind=jprb),intent(out) :: z(:,:)

		integer(kind=jpim) :: ibl,joff,np

		joff = 0

		do ibl=1,ngpblks
			np = blocksize(ibl,nproma,ngptot)

			z(joff+1:joff+np,:) = zb(1:np,:,ip,ibl)

			joff = joff + np
		end do
	end subroutine

	subroutine block2contiguous(nproma,ngpblks,ngptot,nfield,zb,z)
		integer(kind=jpim),intent(in) :: nproma,ngpblks,ngptot,nfield
		real(kind=jprb),intent(in) :: zb(nproma,nfield,ngpblks)
		real(kind=jprb),intent(out) :: z(ngptot,nfield)

		integer(kind=jpim) :: ibl,joff,np

		joff = 0

		do ibl=1,ngpblks
			np = blocksize(ibl,nproma,ngptot)

			z(joff+1:joff+np,:) = zb(1:np,:,ibl)

			joff = joff + np
		end do
	end subroutine

	subroutine contiguous2block1(nproma,ngpblks,ngptot,ip,z,zb)
		integer(kind=jpim),intent(in) :: nproma,ngpblks,ngptot,ip
		real(kind=jprb),intent(in) :: z(:)
		real(kind=jprb),intent(out) :: zb(:,:,:)

		integer(kind=jpim) :: ibl,joff,np

		joff = 0

		do ibl=1,ngpblks
			np = blocksize(ibl,nproma,ngptot)

			zb(1:np,ip,ibl) = z(joff+1:joff+np)

			joff = joff + np
		end do
	end subroutine
end module
