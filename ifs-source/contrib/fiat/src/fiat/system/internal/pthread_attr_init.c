/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#if defined(__GNUC__) || defined(__PGI)

/* pthread_attr_init() interception to reset guard region size 
   between thread stacks, by S.Saarinen, 30-Sep-2016 */

/* 
   See : man pthread_attr_init

   A custom pthread_attr_init() mainly to control the memory "gap" or protected "guard size" between slave threads 
   Slave threads allocate memory from the heap -- only master thread allocates from (genuine) stack
   Guard region is usually very small and it is possible that a slave thread may (accidentally) overwrite
   to the adjacent slave threads memory arena. By setting the guard size big enough, usually a few MiB suffice,
   then an overwrite could hit this protected memory area easier triggering usually a SIGSEGV.
   The traceback that follows usually shows the code location and allows to fix the issue for good.

   To set the guard size use : export EC_THREAD_GUARDSIZE=<value>[G|M|K]

   Caveat: Some MPI drivers may not like this -- check especially the Mellanox HPC-X/OpenMPI
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <pthread.h>
#include <dlfcn.h>
#include <sys/types.h>

#undef RNDUP_DIV
#define RNDUP_DIV(i,n) (( (i) + (n) - 1 ) / (n))

#undef RNDUP
#define RNDUP(i,n) ( RNDUP_DIV(i,n) * (n))

#if defined(RTLD_NEXT)
#define PTR_LIBC RTLD_NEXT
#else
#define PTR_LIBC ((void*) -1L)
#endif

static int (*ptr_pthread_attr_init)(pthread_attr_t *attr) = NULL;
int pthread_attr_init(pthread_attr_t *attr)
{
  int rc;
  static int done = 0;
  FILE *fp = NULL;
  extern int ec_mpirank(), ec_gettid(); // ec_env.c
  int me = ec_mpirank(); // Could be called before the MPI_Init*() or w/o presence of MPI (i.e. rank 0)
  pid_t pid = getpid();
  pid_t tid = ec_gettid();
  int master = (pid == tid) ? 1 : 0;
  if (!ptr_pthread_attr_init) {
    ptr_pthread_attr_init = (int (*)(pthread_attr_t *a))dlsym(PTR_LIBC, "pthread_attr_init");
    if (!ptr_pthread_attr_init) {
      fprintf(stderr,"***Error: Dynamic linking to pthread_attr_init() failed : errno = %d\n",errno);
      abort();
    }
    /* We intend to output only from MPI-task 0, master thread */
    if (!done && me == 0 && master) fp = stderr;
    done = 1;
  }
  rc = ptr_pthread_attr_init(attr);
  {
    char *env_gs = getenv("EC_THREAD_GUARDSIZE");
    if (!env_gs) env_gs = getenv("THREAD_GUARDSIZE"); // Compatibility
    if (env_gs) {
      size_t pgsize = getpagesize();
      size_t guardsize = atoll(env_gs);
      if (strchr(env_gs,'G')) guardsize *= 1073741824; /* hence, in GiB */
      else if (strchr(env_gs,'M')) guardsize *= 1048576; /* hence, in MiB */
      else if (strchr(env_gs,'K')) guardsize *= 1024; /* hence, in KiB */
      guardsize = RNDUP(guardsize,pgsize);
      if (guardsize > pgsize) { /* Now we *DO* bother */
        char *env_omp = getenv("OMP_STACKSIZE");
        size_t omp_stacksize = env_omp ? atoll(env_omp) : 0;
        size_t stacksize = 0;
        int iret = pthread_attr_getstacksize(attr,&stacksize);
#if 1
        if (fp) fprintf(fp,
                        "[%s@%s:%d] [pid=%ld:tid=%ld]: Requesting guard region size "
                        "between thread stacks : %lld bytes (%s PAGESIZE = %ld)\n",
                        __FUNCTION__,__FILE__,__LINE__,
                        (long int)pid,(long int)tid,
                        (long long int)guardsize,
                        (guardsize > pgsize) ? ">" : "<=",
                        (long)pgsize);
#endif
        if (env_omp) {
          if (strchr(env_omp,'G')) omp_stacksize *= 1073741824; /* hence, in GiB */
          else if (strchr(env_omp,'M')) omp_stacksize *= 1048576; /* hence, in MiB */
          else if (strchr(env_omp,'K')) omp_stacksize *= 1024; /* hence, in KiB */
        }
        if (fp) fprintf(fp,
                        "[%s@%s:%d] [pid=%ld:tid=%ld]: Stack size(s) : %lld bytes (def), %lld bytes (OMP) : [iret=%d]\n",
                        __FUNCTION__,__FILE__,__LINE__,
                        (long int)pid,(long int)tid,
                        (long long int)stacksize,
                        (long long int)omp_stacksize,
                        iret);
        if (iret == 0 && omp_stacksize > guardsize) {
          iret = pthread_attr_setguardsize(attr,guardsize);
          (void) pthread_attr_getguardsize(attr,&guardsize);
          if (fp) fprintf(fp,
                          "[%s@%s:%d] [pid=%ld:tid=%ld]: Guard region size now : %lld bytes : [iret=%d]\n",
                          __FUNCTION__,__FILE__,__LINE__,
                          (long int)pid,(long int)tid,
                          (long long int)guardsize,iret);
        }
      }
    }
  }
  if (fp) fflush(fp);
  return rc;
}

#endif /* defined(__GNUC__) || defined(__PGI) */
