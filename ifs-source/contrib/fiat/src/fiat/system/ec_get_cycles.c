/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ec_get_cycles.h"
#include "cas.h"

// Borrowed from PAPI -- currently (27-Sep-2019/SS) just x86, Power & ARM versions
// Fortran callable

#if (defined(__i386__)||defined(__x86_64__))

// X86

long long int ec_get_cycles_()
{
  long long int ret = 0;
#ifdef __x86_64__
  do {
    unsigned int a, d;
    asm volatile ( "rdtsc":"=a" ( a ), "=d"( d ) );
    ( ret ) = ( ( long long int ) a ) | ( ( ( long long int ) d ) << 32 );
  } while ( 0 );
#else
  __asm__ __volatile__( "rdtsc":"=A"( ret ): );
#endif
  return ret;
}

#elif defined(__powerpc__)

// POWER PC

#include <stdint.h>
#include <sys/platform/ppc.h>
#include <sched.h>
#include <stdio.h>
#include <linux/limits.h>

long long int ec_get_cycles_()
{
  uint64_t result;
  int64_t retval;
  unsigned long int dummy;
  static volatile uint64_t multiplier = 1;

  if (multiplier == 1) { // once only
    static volatile sig_atomic_t mylock = 0;
    cas_lock(&mylock);
    if (multiplier == 1) {
      FILE *fp;
      extern int ec_coreid(); // from ec_env.c
      int cpuid = ec_coreid();
      const char max_freq_file_fmt[] = "/sys/devices/system/cpu/cpu%d/cpufreq/cpuinfo_max_freq"; // in kHZ
      char file[PATH_MAX];
      snprintf(file,sizeof(file),max_freq_file_fmt,cpuid);
      fp = fopen(file,"r");
      if (fp) {
        int max_khz = 0;
        if (fscanf(fp,"%d",&max_khz) == 1) {
          multiplier = ((uint64_t)max_khz * 1000000ULL) / __ppc_get_timebase_freq();
        }
        else {
          multiplier = 0; // failed
        }
        fclose(fp);
      }
      else {
        perror(file);
        multiplier = 0; // failed
      }
    }
    cas_unlock(&mylock);
  }

#ifdef __powerpc64__
  /*
    This reads timebase in one 64bit go.  Does *not* include a workaround for the cell (see 
    http://ozlabs.org/pipermail/linuxppc-dev/2006-October/027052.html)
  */
  __asm__ volatile(
      "mftb    %0"
      : "=r" (result));

#else
  /*
    Read the high 32bits of the timer, then the lower, and repeat if high order has changed in the meantime.
    See http://ozlabs.org/pipermail/linuxppc-dev/1999-October/003889.html
  */
  __asm__ volatile(
      "mfspr   %1,269\n\t"  /* mftbu */
      "mfspr   %L0,268\n\t" /* mftb */
      "mfspr   %0,269\n\t"  /* mftbu */
      "cmpw    %0,%1\n\t"   /* check if the high order word has chanegd */
      "bne     $-16"
      : "=r" (result), "=r" (dummy));
#endif
  retval = (result*multiplier)/1000ULL;
  return retval;
}

#elif defined(__aarch64__)

// ARM 64

long long ec_get_cycles_()
{
  register unsigned long ret;
  __asm__ __volatile__ ("isb; mrs %0, cntvct_el0" : "=r" (ret));
  return ret;
}


#else
long long int ec_get_cycles_()
{
  return 0;
}
#endif

long long int ec_get_cycles() { return ec_get_cycles_(); }


