program drhook_ex5
use mpl_module
use yomhook, only : jphook, dr_hook
use sdl_mod, only : sdl_traceback
implicit none
integer j, numargs
integer jpe, npes, mype
character(len=256) arg, env
real(jphook) :: zhook_handle
call mpl_init(ldinfo=.false.)
call dr_hook('drhook_ex5',0,zhook_handle)
npes = mpl_nproc()
mype = mpl_myrank()
do jpe=1,npes
   if (jpe == mype) then
      write(0,1000) mype,&
           & ': Basic MPL/MPI implementation works : # of MPI-tasks = ',npes
      numargs = mpl_iargc()
      write(0,1000) mype,&
           & ': Number of args = ',numargs
      do j=0,numargs
         call mpl_getarg(j,arg)
         write(0,1001) mype, ': arg#', j, ' "'//trim(arg)//'"'
      enddo
      call get_environment_variable('MPICH_ROOT',env)
      write(0,1002) mype, ': env MPICH_ROOT="'//trim(env)//'"'
      call ec_flush(0)
   endif
   call mpl_barrier()
enddo
call mpl_barrier()
call sdl_traceback() ! Testing traceback, too
call mpl_end()
1000  format(i5,a,i5)
1001  format(i5,a,i2,a)
1002  format(i5,a)
call dr_hook('drhook_ex5',1,zhook_handle)
end program drhook_ex5
