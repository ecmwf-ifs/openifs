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

module lececr_netcdf_mod
!---------------------------------------------------------------------------
!
! Procedures de lecture et ecriture de variables dans un fichier netcdf
!
! Fonction de lecture :
! ---------------------
!  lecnetcdf(ncid,cvarname,p_values) : procedure generique, lecture d'un champ netcdf.
!          ncid : identificateur du fichier netcdf
!          cvarname : nom du champ netcdf a recuperer
!          p_values : pointeur sur un tableau du meme type que le champ a recuperer 
!                     (real,double precision, integer), de dimension 1 a 4
!  lecatt(ncid,cvarname,cattname,value) :  procedure generique, lecture d'un attribut
!          cvarname : nom de la variable
!          cattname : nom de l'attribut
!          value : valeur de l'attribut (de type chaine, real, double precision  ou integer)
!
! Fonctions de definition :
! -------------------------
!  defdim(ncid,cdimname,idimsize,xtype) : definition d'une dimension
!          xtype : type netcdf de la dimension 
!          (NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT, NF90_FLOAT)
!
!  defvar(ncid,cvarname,ctabnomdim,xtype) : definition d'une variable netcdf
!          ctabmondim : vecteur contenant la liste des noms des dimensions de la variable
!
!  defatt(ncid,cvarname,cattname,value) : procedure generique, definition d'un attribut.
!          cvarname : nom de la variable
!          cattname : nom de l'attribut
!          value : valeur de l'attribut (de type chaine, real, double precision  ou integer)
!
! Fonction d'ecriture :
! ---------------------
!  putvar(ncid,cvarname,tabval) : procedure generique, ecriture d'un champ netcdf.
!      (les dimensions et le champ doivent avoir ete definis prealablement)
!          tabval : tableau de dimension 1 a 4, de type : real, double precision ou integer
!
! Florence Favot  08/2012 : original
!                 10/2012 : ajout lecatt
!---------------------------------------------------------------------------

interface lecnetcdf
    module procedure lecnetcdf_4d_f,lecnetcdf_3d_f,lecnetcdf_2d_f,lecnetcdf_1d_f, &
                     lecnetcdf_4d_i,lecnetcdf_3d_i,lecnetcdf_2d_i,lecnetcdf_1d_i
end interface

interface lecatt
    module procedure lecatt_c, lecatt_f, lecatt_i
end interface
 
interface defatt
    module procedure defatt_c, defatt_f, defatt_i
end interface

interface putvar
    module procedure putvar_1d_f, putvar_2d_f, putvar_3d_f, putvar_4d_f, &
                     putvar_1d_i, putvar_2d_i, putvar_3d_i, putvar_4d_i
end interface

contains
!===========================================================================
subroutine handle_err(status)

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim), intent ( in) :: status
     
if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
end if

end subroutine handle_err
!*******************************************************************

!===================================================================
! fonctions de l'interface generique lecnetcdf                     !
!===================================================================
! lecture des champs netcdf de dimension 1 a 4 :
!
!    lecnetcdf(ncid,cvarname,p_values)
!              ncid : identificateur du fichier netcdf
!              cvarname : nom du champ netcdf a recuperer
!              p_values : pointeur sur un tableau du meme type que le champ a recuperer
!                         (real,double precision, integer)
!*******************************************************************
! fonctions pour p_values de type : real
!*******************************************************************
!*******************************************************************
subroutine lecnetcdf_4d_f(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 4d de reels
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:,:,:,:),pointer       :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if (istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 4) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1),dimtab(2),dimtab(3),dimtab(4)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 4. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_4d_f
!*******************************************************************

subroutine lecnetcdf_3d_f(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 3d de reels
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:,:,:), pointer        :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if(istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 3) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1),dimtab(2),dimtab(3)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 3. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_3d_f
!*******************************************************************

subroutine lecnetcdf_2d_f(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 2d de reels
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:,:), pointer          :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if(istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 2) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1),dimtab(2)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 2. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_2d_f
!*******************************************************************

subroutine lecnetcdf_1d_f(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 1d de reels
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:), pointer            :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if(istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 1) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 1. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_1d_f

!*******************************************************************
! fonctions pour p_values de type : integer
!*******************************************************************
subroutine lecnetcdf_4d_i(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 4d d'entiers
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:,:,:,:),pointer    :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if(istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 4) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1),dimtab(2),dimtab(3),dimtab(4)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 4. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_4d_i
!*******************************************************************
subroutine lecnetcdf_3d_i(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 3d d'entier
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:,:,:),     pointer :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if(istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 3) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1),dimtab(2),dimtab(3)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 3. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_3d_i
!*******************************************************************
subroutine lecnetcdf_2d_i(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 2d d'entier
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:,:),       pointer :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if(istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 2) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1),dimtab(2)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 2. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_2d_i
!*******************************************************************
subroutine lecnetcdf_1d_i(ncid,cvarname,p_values)
!
! recupere le champ de nom "cvarname" dans un tableau 1d d'entier
!
!-----------------------------------------------------------------

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:),         pointer :: p_values

integer(kind=jpim)                               :: varid, istat
integer(kind=jpim)                               :: i, ndimvar
integer(kind=jpim), dimension(:),allocatable     :: dimtab
integer(kind=jpim), dimension(nf90_max_var_dims) :: dimIDs   

! recup de l'identificateur du champ
istat = nf90_inq_varid(ncid, cvarname, varid)
if (istat /= nf90_NoErr) then
    print *,'erreur variable : ',cvarname
    call handle_err(istat)
end if

! nombre de dimensions et identificateurs 
istat = nf90_inquire_variable(ncid, varid, ndims=ndimvar, dimids = dimIDs)
if(istat /= nf90_NoErr) call handle_err(istat)

! allocation du pointeur si la variable a le bon nombre de dimensions
if (ndimvar .eq. 1) then
    ! allocation du tableau contenant les dimensions
    allocate(dimtab(ndimvar))
    do i=1,ndimvar
        istat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimtab(i))
        if(istat /= nf90_NoErr) call handle_err(istat)
    enddo
    ! allocation du tableau contenant les donnees
    allocate(p_values(dimtab(1)))
    deallocate(dimtab)
    ! recup des donnees
    istat = nf90_get_var(ncid, varid, p_values)
    if (istat /= nf90_NoErr) then
        print *,'erreur lecture variable : ',cvarname
        call handle_err(istat)
    end if
else
    print *,'erreur :',cvarname,' n''est pas de dimension 1. Sa dimension est ',ndimvar
    stop
end if

end subroutine lecnetcdf_1d_i

!===================================================================
! fin fonctions de l'interface generique lecnetcdf                 !
!===================================================================

!===================================================================
! fonctions de l'interface generique lecatt                        !
!===================================================================
!*******************************************************************
subroutine lecatt_c(ncid,cvarname,cattname,value,istat)
!-------------------------------------------------------------------
! lecture d'un attribut d'une variable ayant une valeur de type chaine
use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
character(len=*),          intent(in) :: cattname
character(len=*),         intent(out) :: value
integer(kind=jpim),                  intent(out) :: istat

integer(kind=jpim)                               :: idvar

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'lecatt erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! lecture de l'attribut cattname de la variable 
istat=nf90_get_att(ncid,idvar,cattname,value)
!if (istat /= nf90_NoErr) then
!    print *,'lecatt erreur : pour attribut ',cattname,' de la variable ',cvarname
!end if

end subroutine lecatt_c
!*******************************************************************
subroutine lecatt_f(ncid,cvarname,cattname,value,istat)
!-------------------------------------------------------------------
! lecture d'un attribut d'une variable ayant une valeur de type chaine
use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
character(len=*),          intent(in) :: cattname
real(kind=jprb),                     intent(out) :: value
integer(kind=jpim),                  intent(out) :: istat

integer(kind=jpim)                               :: idvar

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'lecatt erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! lecture de l'attribut cattname de la variable 
istat=nf90_get_att(ncid,idvar,cattname,value)
!if (istat /= nf90_NoErr) then
!    print *,'lecatt erreur : pour attibut ',cattname,' de la variable ',cvarname
!end if

end subroutine lecatt_f
!*******************************************************************
subroutine lecatt_i(ncid,cvarname,cattname,value,istat)
!-------------------------------------------------------------------
! lecture d'un attribut d'une variable ayant une valeur de type chaine
use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
character(len=*),          intent(in) :: cattname
integer(kind=jpim),                  intent(out) :: value
integer(kind=jpim),                  intent(out) :: istat

integer(kind=jpim)                               :: idvar

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'lecatt erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! lecture de l'attribut cattname de la variable 
istat=nf90_get_att(ncid,idvar,cattname,value)
!if (istat /= nf90_NoErr) then
!    print *,'lecatt erreur : pour attibut ',cattname,' de la variable ',cvarname
!end if

end subroutine lecatt_i
!*******************************************************************
!===================================================================
! fin fonctions de l'interface generique lecatt                    !
!===================================================================

!===================================================================
! fonctions pour ecriture dans un fichier netcdf                   !
!===================================================================
!*******************************************************************
subroutine defdim(ncid,cdimname,idimsize,xtype)
! definition d'une dimension
! xtype : type netcdf de la dimension 
!         (NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT, NF90_FLOAT)

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cdimname
integer(kind=jpim),                   intent(in) :: idimsize
integer(kind=jpim),                   intent(in) :: xtype   

integer(kind=jpim)                               :: iddim, idvar, istat

! definition de la dimension
istat=nf90_def_dim(ncid,trim(cdimname),idimsize,iddim)
if (istat /= nf90_NoErr) then
    print *,'defdim erreur : dimension ',trim(cdimname)
    call handle_err(istat)
end if

! definition de la variable dimension
istat=nf90_def_var(ncid,trim(cdimname),xtype,(/ iddim /),idvar)
if (istat /= nf90_NoErr) then
    print *,'defdim erreur : variable ',trim(cdimname)
    call handle_err(istat)
end if

end subroutine defdim
!*******************************************************************

!*******************************************************************
subroutine defvar(ncid,cvarname,ctabnomdim,xtype)
!-------------------------------------------------------------------
! definition d'une variable netcdf
!
! xtype : type de la variable 
!         (NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT, NF90_FLOAT)

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
character(len=*),dimension(:),intent(in) :: ctabnomdim
integer(kind=jpim),                   intent(in) :: xtype   

integer(kind=jpim)                               :: idvar, istat
integer(kind=jpim)                               :: i, ndimvar

integer(kind=jpim), dimension(:),allocatable     :: dimidtab

! recup des identificateurs des dimensions
ndimvar=size(ctabnomdim)
!print *,'def_var : ',cvarname,' nb dim =',ndimvar

allocate(dimidtab(size(ctabnomdim)))
do i=1,ndimvar
    istat = nf90_inq_dimid(ncid, trim(ctabnomdim(i)), dimidtab(i))
    if (istat /= nf90_NoErr) then 
        print *,'erreur dimension ',trim(ctabnomdim(i))
        call handle_err(istat)
    end if
enddo

!------ définition de la variable
istat = nf90_def_var(ncid,trim(cvarname),xtype,dimidtab,idvar)
deallocate(dimidtab)
if (istat /= nf90_NoErr) call handle_err(istat)

end subroutine defvar

!===================================================================
! fonctions de l'interface generique defatt                        !
!===================================================================
subroutine defatt_c(ncid,cvarname,cattname,value)
!-------------------------------------------------------------------
! definition d'un attribut d'une variable ayant une valeur de type chaine
use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
character(len=*),          intent(in) :: cattname
character(len=*),          intent(in) :: value

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'deffatt erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! definition de l'attribut associe a la variable 
istat=nf90_put_att(ncid,idvar,cattname,value)

end subroutine defatt_c
!*******************************************************************
subroutine defatt_f(ncid,cvarname,cattname,value)
!-------------------------------------------------------------------
! definition d'un attribut d'une variable ayant une valeur de type reel
use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
character(len=*),          intent(in) :: cattname
real(kind=jprb),                      intent(in) :: value

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'deffatt erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! definition de l'attribut associe a la variable 
istat=nf90_put_att(ncid,idvar,cattname,value)

end subroutine defatt_f
!*******************************************************************
subroutine defatt_i(ncid,cvarname,cattname,value)
!-------------------------------------------------------------------
! definition d'un attribut d'une variable ayant une valeur de type entier
use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
character(len=*),          intent(in) :: cattname
integer(kind=jpim),                   intent(in) :: value

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'deffatt erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! definition de l'attribut associe a la variable 
istat=nf90_put_att(ncid,idvar,cattname,value)

end subroutine defatt_i
!===================================================================
! fin fonctions de l'interface generique defatt                    !
!===================================================================

!===================================================================
! fonctions de l'interface generique putvar                        !
!===================================================================
subroutine putvar_1d_f(ncid,cvarname,tabval)
! ecrit le champ 1d de reels dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:),         intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_1d_f
!*******************************************************************
subroutine putvar_2d_f(ncid,cvarname,tabval)
! ecrit le champ 2d de reels dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:,:),       intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_2d_f
!*******************************************************************
subroutine putvar_3d_f(ncid,cvarname,tabval)
! ecrit le champ 3d de reels dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:,:,:),     intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_3d_f
!*******************************************************************
subroutine putvar_4d_f(ncid,cvarname,tabval)
! ecrit le champ 4d de reels dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
real(kind=jprb),dimension(:,:,:,:),   intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_4d_f
!*******************************************************************
subroutine putvar_1d_i(ncid,cvarname,tabval)
! ecrit le champ 1d d'entier dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:),      intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_1d_i
!*******************************************************************
subroutine putvar_2d_i(ncid,cvarname,tabval)
! ecrit le champ 2d d'entier dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:,:),    intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_2d_i
!*******************************************************************
subroutine putvar_3d_i(ncid,cvarname,tabval)
! ecrit le champ 3d d'entier dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:,:,:),  intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_3d_i
!*******************************************************************
subroutine putvar_4d_i(ncid,cvarname,tabval)
! ecrit le champ 4d d'entier dans un fichier netcdf,
! le nom de la variable dans le fichier netcdf est cvarname

use netcdf
use parkind1, only : jpim,jprb
implicit none
include 'netcdf.inc'

integer(kind=jpim),                   intent(in) :: ncid
character(len=*),          intent(in) :: cvarname
integer(kind=jpim),dimension(:,:,:,:),intent(in) :: tabval

integer(kind=jpim)                               :: idvar, istat

! recup de l'identificateur de la variable
istat = nf90_inq_varid(ncid, cvarname, idvar)
if (istat /= nf90_NoErr) then
    print *,'putvar erreur : identificateur variable ',cvarname
    call handle_err(istat)
end if

! ecriture des valeurs 
istat = nf90_put_var(ncid,idvar,tabval)
if (istat /= nf90_noerr) then
    print *,'putvar erreur : ecriture variable ',cvarname
    call handle_err(istat)
end if

end subroutine putvar_4d_i
!*******************************************************************
!===================================================================
! fin fonctions de l'interface generique putvar                    !
!===================================================================

end module lececr_netcdf_mod
