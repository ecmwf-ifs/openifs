/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

//
// Explicit hugepages (user application program needs to be linked or LD_PRELOADed by /usr/lib64/libhugetlbfs.so)
//
// To build:
//
// gcc -DMAIN_HUGEPAGES reserve_hugepages.c -o reserve_hugepages.x
//
// Need to be a root to run this -- or as a "sudo" do the following:
//
// (1) sudo chown root reserve_hugepages.x
// (2) sudo chmod u+s reserve_hugepages.x
// (3) ./reserve_hugepages.x pernode_MB verbosity  # NB: no need for sudo anymore (verbosity=0 turns messages off)
// or with mpirun/mpiexec: one task per reserved node
//    mpirun -np number_of_nodes -ppn 1 ./reserve_hugepages.x pernode_MB verbosity
// (4) Release : run the same as in (3) but put pernode_MB=0
//
// See also from /proc/meminfo how much MemTotal you have per node and take max 80% of that -- rest for small pages
//
// Adapted from Peter Boyle (see https://arxiv.org/pdf/1711.04883.pdf)
// by Sami Saarinen, ECMWF, 08-Mar-2018
//

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

static int reserve_2M_hugepages(long long int bytes, int verbose)
// reserve_2M_hugepages : reserves 2M hugepages
{
  int rc = 0;
  if (bytes >= 0) {
    rc = setuid(0);
    if (rc == 0) {
#if 1
      // Mikko's method
      long long int npages = bytes / (long long int)2097152;
      char cmd[4096];
      snprintf(cmd,sizeof(cmd),"/bin/echo %lld > /sys/kernel/mm/hugepages/hugepages-2048kB/nr_hugepages",npages);
      system(cmd);
      snprintf(cmd,sizeof(cmd),"/bin/echo %lld > /proc/sys/vm/nr_hugepages",npages);
      system(cmd);
      if (verbose) {
	//snprintf(cmd,sizeof(cmd),"/usr/bin/more /sys/kernel/mm/hugepages/hugepages-2048kB/nr_hugepages /proc/sys/vm/nr_hugepages < /dev/null | /usr/bin/cat");
	//system(cmd);
	snprintf(cmd,sizeof(cmd),"h=$(hostname -s); egrep -i ^hugepage /proc/meminfo | perl -pe 's/^/'${h}': /'");
	system(cmd);
      }
#else
      rc = access("/usr/bin/hugeadm",X_OK);
      if (rc == 0) {
	char host[512];
	(void) gethostname(host,sizeof(host));
	snprintf(cmd,sizeof(cmd),"/usr/bin/hugeadm --pool-pages-min=2M:%lld",npages);
	if (verbose) {
	  printf("%s: %s\n",host,cmd);
	  fflush(stdout);
	}
	system(cmd);
	if (verbose) {
	  system("/usr/bin/hugeadm --pool-list");
	  system("/usr/bin/cat /sys/kernel/mm/hugepages/hugepages-2048kB/nr_hugepages");
	  system("/usr/bin/cat /proc/sys/vm/nr_hugepages");
	}
      }
      else {
	fprintf(stderr,"Error: Unable to execute /usr/bin/hugeadm\n");
      }
#endif
    }
    else {
      fprintf(stderr,"Error: Unable to run as root\n");
    }
  }
  return rc;
}

#ifdef MAIN_HUGEPAGES
int main(int argc, char **argv)
{
  // Usage: reserve_hugepages.x pernode_MB verbosity
  long long int pernode_MB = 0;
  long long int bytes = 0;
  int verbose = 1;
  if (argc >= 2) {
    pernode_MB = atoll(argv[1]);
    if (pernode_MB < 0) pernode_MB = 0;
    if (argc > 2) verbose = atoi(argv[2]);
  }
  bytes = pernode_MB * (long long int)1048576;
  return reserve_2M_hugepages(bytes,verbose);
}
#endif
