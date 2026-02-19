/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


/* env.c */

/* Implement Fortran-callable ec_getenv and ec_putenv,
   since not all environments have getenv & putenv,
   but Unix/C library always have them */

/* Author: Sami Saarinen, ECMWF, 15-Mar-2006 */

//if (defined(__GNUC__) || defined(__PGI))
#define _GNU_SOURCE
//endif

#if defined(__APPLE__)
#include <pthread.h>
#endif

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
#include <sys/syscall.h>

#include "raise.h"

#if !defined(HOST_NAME_MAX) && defined(_POSIX_HOST_NAME_MAX)
#define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
#endif

#if !defined(HOST_NAME_MAX) && defined(_SC_HOST_NAME_MAX)
#define HOST_NAME_MAX _SC_HOST_NAME_MAX
#endif

#if defined(HOST_NAME_MAX)
#define EC_HOST_NAME_MAX HOST_NAME_MAX
#else
#define EC_HOST_NAME_MAX 512
#endif

extern char **environ; /* Global Unix var */
static int numenv = 0;

void
ec_numenv_bind_c(int *n)
{ /* Returns the number of environment variables currently active */
  int j=0;
  if (environ) {
    for (; environ[j]; j++) { }
  }
  if (n) *n = j;
  numenv = j; /* Not thread-safe */
}

void
ec_environ_bind_c(const int *i,
	   char *value,
	   /* Hidden arguments */
	   const int valuelen)
{ /* Returns (*i)'th environment number; 
     Note: "Fortran", not "C" range between [1..numenv] */
  int j = (i && environ) ? (*i) : 0;
  memset(value, ' ', valuelen);
  if (j >= 1 && j <= numenv) {
    char *p = environ[--j];
    if (p) {
      int len = strlen(p);
      if (valuelen < len) len = valuelen;
      memcpy(value,p,len);
    }
  }
}

void
ec_getenv_bind_c(const char *s,
	   char *value,
	   /* Hidden arguments */
	   int slen,
	   const int valuelen)
{
  char *env = NULL;
  char *p = malloc(slen+1);
  if (!p) {
    fprintf(stderr,"ec_getenv_(): Unable to allocate %d bytes of memory\n", slen+1);
    ABOR1("ec_getenv_(): Unable to allocate memory");
  }
  memcpy(p,s,slen);
  p[slen]='\0';
  memset(value, ' ', valuelen);
  env = getenv(p);
  if (env) {
    int len = strlen(env);
    if (valuelen < len) len = valuelen;
    memcpy(value,env,len);
  }
  free(p);
}

void
ec_putenv_overwrite_bind_c(const char *s,
	   /* Hidden argument */
	   int slen)
{
  const char *x = &s[slen-1];
  /* strip trailing blanks first */
  while (slen > 0 && *x == ' ') { --slen; --x; }
  /* now go ahead */
  if (slen > 0) {
    char *p = malloc(slen+1);
    if (!p) {
      fprintf(stderr,"ec_putenv_(): Unable to allocate %d bytes of memory\n", slen+1);
      ABOR1("ec_putenv_(): Unable to allocate memory");
    }
    memcpy(p,s,slen);
    p[slen]='\0';
    putenv(p);
    /* Cannot free(p); , since putenv() uses this memory area for good ;-( */
  }
}

void
ec_putenv_nooverwrite_bind_c(const char *s,
		       /* Hidden argument */
		       int slen)
{
  const char *x = &s[slen-1];
  /* strip trailing blanks first */
  while (slen > 0 && *x == ' ') { --slen; --x; }
  /* now go ahead */
  if (slen > 0) {
    char *eq = NULL;
    char *p = malloc(slen+1);
    if (!p) {
      fprintf(stderr,"ec_putenv_nooverwrite_(): Unable to allocate %d bytes of memory\n", slen+1);
      ABOR1("ec_putenv_nooverwrite_(): Unable to allocate memory");
    }
    memcpy(p,s,slen);
    p[slen]='\0';
    eq = strchr(p,'=');
    if (eq) {
      char *env = NULL;
      *eq = '\0';
      env = getenv(p);
      if (env) {
	/* Already found ==> do not overwrite */
	free(p);
	return;
      }
      else {
	/* Reset '=' back and continue with putenv() */
	*eq = '=';
      }
    }
    putenv(p);
    /* Cannot free(p); , since putenv() uses this memory area for good ;-( */
  }
}

/*--- sleep_by_spinning ---*/

static int sleep_by_spinning(long secs, long nanosecs) { /* see also drhook.c */
  /* This does not call sleep() at all i.e. is not SIGALRM driven */
  int rc;
  struct timespec req, rem;
  req.tv_sec = secs;
  req.tv_nsec = nanosecs;
  rc = nanosleep(&req, &rem);
  if (rc == -1) {
    if (errno == EINTR) {
      rc = rem.tv_sec;
    }
    else
      rc = 0; /* Can't do much more about this */
  }
  return rc;
}

unsigned int
ec_sleep_(const int *nsec)
{
  //return sleep((nsec && *nsec > 0) ? *nsec : 0);
  return sleep_by_spinning((nsec && *nsec > 0) ? *nsec : 0, 0);
}

unsigned int
ec_sleep(const int nsec)
{
  return ec_sleep_(&nsec);
}

/* Microsecond-sleep, by S.Saarinen, 25-Jan-2008 */

void  /* Global, C-callable, too */
ec_microsleep(int usecs) {
  if (usecs > 0) {
    struct timeval t;
    t.tv_sec =  usecs/1000000;
    t.tv_usec = usecs%1000000;
    // (void) select(0, NULL, NULL, NULL, &t);
    (void) sleep_by_spinning(t.tv_sec, (long)1000*t.tv_usec);
  }
}

void
ec_usleep_(const int *usecs)
{
  if (usecs && *usecs > 0) ec_microsleep(*usecs);
}

void
ec_usleep(const int *usecs)
{
  ec_usleep_(usecs);
}

/* ec_gethostname   , by S.Saarinen, 30-Sep-2016 */
/* ec_getpaddedhost ,       -"-    , 13-Jul-2021 */ 

void ec_getpaddedhost_(char a[],
                       const int *padding, /* a char */
                       /* Hidden argument */
                       int alen)
{
  char s[EC_HOST_NAME_MAX];
  int c = padding ? *padding : (int) ' ';
  memset(a,c,alen);
  if (gethostname(s,sizeof(s)) == 0) {
    int len;
    char *pdot = strchr(s,'.');
    if (pdot) *pdot = '\0'; // cut short from "." char e.g. hostname.fmi.fi becomes just "hostname"
    len = strlen(s);
    if (len > alen) len = alen;
    memcpy(a,s,len);
  }
}

void ec_getpaddedhost(char a[],
                      const int *padding, /* a char */
                      /* Hidden argument */
                      int alen)
{
  ec_getpaddedhost_(a,padding,alen);
}

void ec_gethostname_(char a[], 
		     /* Hidden argument */
		     int alen)
{
  char s[EC_HOST_NAME_MAX];
  memset(a,' ',alen);
  if (gethostname(s,sizeof(s)) == 0) {
    int len;
    char *pdot = strchr(s,'.');
    if (pdot) *pdot = '\0'; // cut short from "." char e.g. hostname.fmi.fi becomes just "hostname"
    len = strlen(s);
    if (len > alen) len = alen;
    memcpy(a,s,len);
  }
}

void ec_gethostname(char a[], 
		     /* Hidden argument */
		     int alen)
{
  ec_gethostname_(a,alen);
}

/* ec_coreid(): For checking runtime affinities (not setting them, though) */

#if defined(LINUX) && !defined(__NEC__)
#include <sched.h>
int sched_getcpu(void);
#define getcpu() sched_getcpu()
#else
#define getcpu() -1
#endif

int ec_coreid_()
{
  return getcpu();
}

int ec_coreid()
{
  return ec_coreid_();
}


#ifndef SYS_gettid
#define SYS_gettid __NR_gettid
#endif

static pid_t GETtid() {
#if defined(__APPLE__)
  uint64_t tid64;
  pthread_threadid_np(NULL, &tid64);
  pid_t tid = (pid_t)tid64;
#else
  pid_t tid = syscall(SYS_gettid);
#endif
  return tid;
}

int ec_getpid_() { // Your Fortran did not recognize GETPID() ? Use then: mypid = ec_getpid()
  return (int)getpid();
}

int ec_getpid() { // C-callable
  return ec_getpid_();
}

int ec_gettid_() { // Your Fortran did not recognize GETTID() ? Use then: mypid = ec_gettid()
  return (int)GETtid();
}

int ec_gettid() { // C-callable
  return ec_gettid_();
}

#if defined(__APPLE__)
  // These variables and functions are Linux specific
  static int CPU_SETSIZE = 0;
  typedef struct cpu_set_t cpu_set_t;
  struct cpu_set_t {};
  static int CPU_ISSET(int index, const cpu_set_t* cpu_set) { return 0; }
  static void CPU_ZERO(cpu_set_t* cpu_set) {}
  static int sched_getaffinity(pid_t pid, size_t size, cpu_set_t* cpu_set) { return 0; }
#endif

static char *cpuset_to_cstr(const cpu_set_t *cpuset, char str[], const int lenstr)
{
  char *ptr = str;
  int len = lenstr;
  int off = 0;
  int i, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, cpuset)) {
      int j, run = 0;
      if (!entry_made) entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, cpuset)) run++;
        else break;
      }
      if (!run) {
        snprintf(ptr+off, len-off, "%d,", i);
      }
      else if (run == 1) {
        snprintf(ptr+off, len-off, "%d,%d,", i, i + 1);
        i++;
      }
      else {
        snprintf(ptr+off, len-off, "%d-%d,", i, i + run);
        i += run;
      }
      off = strlen(ptr);
    }
  }
  *(ptr+off-entry_made) = '\0'; // remove char (usually ",")
  return ptr;
}

/* ec_affinity(): Returns compact thread affinity string to Fortran */

void ec_affinity_(char s[],
                  // hidden
                  const int slen)
{
  if (slen > 0) {
    cpu_set_t coremask;
    char *p, str[CPU_SETSIZE * 8];
    pid_t tid = GETtid();
    CPU_ZERO(&coremask);
    (void) sched_getaffinity(tid,sizeof(coremask),&coremask);
    p = cpuset_to_cstr(&coremask, str, sizeof(str));
    memset(s,' ',slen);
    if (p) {
      int plen = strlen(p);
      if (plen > slen) plen = slen;
      memcpy(s,p,plen);
    }
  }
}

void ec_affinity(char s[],
                 // hidden
                 const int slen)
{
  ec_affinity_(s,slen);
}

#ifdef DARSHAN
/* Some issues with Darshan -- better to use our own version of MPI_Wtime (mpi_wtime_ in Fortran) */
double mpi_wtime_()
{
  extern double util_walltime_(); /* from drhook.c */
  return util_walltime_();
}
#endif

/* ec_mpi_epoch(): Is used (indirectly) to calculate pretty accurate MPI_Init*() overhead */

double ec_mpi_epoch_()
{
  double epoch = 0;
#if 1
  // Avoids fork()
  struct timeval tbuf;
  if (gettimeofday(&tbuf,NULL) == -1) perror("ec_mpi_epoch_");
  epoch = (double) tbuf.tv_sec + (tbuf.tv_usec / 1000000.0); // should be aligned with /bin/date
#else
  // With (v)fork() -- do NOT use this
  FILE *fp = popen("/bin/date +%s.%N","r");
  if (fp) {
    fscanf(fp,"%lf",&epoch);
    pclose(fp);
  }
#endif
  return epoch;
}

double ec_mpi_epoch() { return ec_mpi_epoch_(); }

/* ec_cpumodel(): Returns CPU model of this platform */

void ec_cpumodel_(char *s,
                  /* hidden arg */
                  const int slen)
{
  FILE *fp = fopen("/proc/cpuinfo","r");
  memset(s,' ',slen);
  if (fp) {
    const char target[] = "model name\t: ";
    int tlen = strlen(target);
    char p[4096];
    while (!feof(fp) && fgets(p,sizeof(p),fp)) {
      if (strncmp(p,target,tlen) == 0) { // egrep ^<target>
        int len;
        char *nl = strchr(p,'\n');
        if (nl) *nl = '\0';
        len = strlen(p+tlen);
        if (slen < len) len = slen;
        memcpy(s,p+tlen,len);
        break;
      }
    }
    fclose(fp);
  }
}

void ec_cpumodel(char *s,
                 /* hidden arg */
                 const int slen)
{
  ec_cpumodel_(s, slen);
}

/* ec_mpirank(): Sometimes we need to know MPI_COMM_WORLD's "myrank" aka "me" before MPI_Init*() was even been called */
/* ec_mpisize(): Sometimes we need to know MPI_COMM_WORLD's number of ranks before MPI_Init*() was even been called */

int ec_mpirank_()
{
  static int me = -1; /* MPI task id >= 0 && <= max tasks - 1 */
  if (me < 0) {
    /* Trying to figure out MPI task id *before* MPI_Init*() was called */
    char *env_mpirank = NULL;
    if (!env_mpirank) env_mpirank = getenv("PMI_FORK_RANK"); // Cray MPICH : when invoked with PMI_NO_FORK=1
    if (!env_mpirank) env_mpirank = getenv("ALPS_APP_PE"); // Cray ALPS
    if (!env_mpirank) env_mpirank = getenv("PMIX_RANK"); // OpenMPI when using "srun" and SLURM_MPI_TYPE=pmix
    if (!env_mpirank) env_mpirank = getenv("PMI_RANK"); // MPICH (except Cray MPICH) -- also SLURM "srun" regardless of MPI implementation
    if (!env_mpirank) env_mpirank = getenv("OMPI_COMM_WORLD_RANK"); // Genuine OpenMPI
    if (!env_mpirank) env_mpirank = getenv("EC_FARM_ID"); // ECMWF extension
    if (env_mpirank) me = atoi(env_mpirank);
    if (me < 0) me = 0; // Bailing out
  }
  return me;
}

int ec_mpirank() { return ec_mpirank_(); }

int ec_mpisize_()
{
  static int numranks = 0;
  if (numranks < 1) {
    char *env_mpisize = NULL;
    // Please note : Cray ALPS does NOT have any specific "ALPS_APP_SIZE" or similar !!!
    if (!env_mpisize) env_mpisize = getenv("PMIX_SIZE"); // OpenMPI when using srun and SLURM_MPI_TYPE=pmix : actually this ALSO does NOT exist !!!
    if (!env_mpisize) env_mpisize = getenv("PMI_SIZE"); // MPICH (except Cray MPICH) -- also SLURM "srun" regardless of MPI implementation (also for SLURM_MPI_TYPE=pmix)
    if (!env_mpisize) env_mpisize = getenv("OMPI_COMM_WORLD_SIZE"); // Genuine OpenMPI
    if (!env_mpisize) env_mpisize = getenv("SLURM_NTASKS"); // If SLURM (--ntasks=value or -n value)
    if (!env_mpisize) env_mpisize = getenv("SLURM_NPROCS"); // If SLURM (--ntasks=value or -n value) -- backward compatible with SLURM_NTASKS
    if (!env_mpisize) env_mpisize = getenv("EC_FARM_SIZE"); // ECMWF extension
    if (env_mpisize) numranks = atoi(env_mpisize);
    if (numranks < 1) numranks = 1;
  }
  return numranks;
}

int ec_mpisize() { return ec_mpisize_(); }
