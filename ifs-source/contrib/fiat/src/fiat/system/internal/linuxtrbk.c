/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


/* linuxtrbk.c : Print traceback on linux */

/* 
   Author: Sami Saarinen, ECMWF, 28-Apr-2006
   The code "nicked" from ifsaux/support/drhook.c

*/

//#if (defined(__GNUC__) || defined(__PGI))
#define _GNU_SOURCE
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <alloca.h>
#include <math.h>
#include "ec_args.h"
#include "drhook.h"

#define PRETOSTR(x) #x
#define TOSTR(x) PRETOSTR(x)

#define strequ(s1,s2)     ((void *)s1 && (void *)s2 && strcmp(s1,s2) == 0)
#define strnequ(s1,s2,n)  ((void *)s1 && (void *)s2 && memcmp(s1,s2,n) == 0)

#ifdef WITHOUT_CXXDEMANGLE
static char *cxxdemangle(const char *mangled_name, int *status) {
  if( status ) *status = 1;
  return NULL;
}
#endif

typedef struct {
  const char *func;
  const char *file;
  unsigned int lineno;
} BFD_t;

static int ResolveViaBFD(void *address, BFD_t *b, const char *str);
static void InitBFD();

#ifndef LINELEN
#define LINELEN 1024

#define FFL __FUNCTION__,__FILE__,__LINE__

#if (defined(__GNUC__) || defined(__PGI))

#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>

#if defined(__APPLE__)
// define _XOPEN_SOURCE to enable deprecated use of ucontext.h
#define _XOPEN_SOURCE
#endif

#if !defined(CYGWIN) && !defined(__NEC__)
#include <execinfo.h>
#include <ucontext.h>
#endif

#endif /* defined(__GNUC__) */

#if !defined(ADDR2LINE)
#define ADDR2LINE /usr/bin/addr2line
#endif

#define len_addr2linecmd (sizeof(TOSTR(ADDR2LINE)) + 80 + 4096 + GNUC_BTRACE * 30)
static char prealloc_addr2linecmd[len_addr2linecmd] = "";

#ifdef BFDLIB
#include <bfd.h>
static const int has_bfd = 1;
#else
static const int has_bfd = 0;
#endif

void gdb_trbk(); // defined below
void dbx_trbk(); // defined below
void LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr); // defined below

#ifdef __NEC__
void
LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr)
{
#if 1
  // When VE_TRACEBACK=ALL the following should have the same impact as with implicit traceback (when compiled & linked with -traceback)
  // NB: No control on output channel
  if (!sigcontextptr) sigcontextptr = __builtin_frame_address(0);
  __builtin_traceback(sigcontextptr);
#else
  gdb_trbk();
#endif
}
#elif defined(CYGWIN)
void
LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr)
{
  const char *pfx = prefix  ? prefix  : drhook_PREFIX(0);
  const char *ts  = timestr ? timestr : drhook_TIMESTR(0);
  int sigcontextptr_given = sigcontextptr ? 1 : 0;
  static int recur = 0;
  const char *a_out = ec_argv()[0];
  fprintf(stderr,"%s %s [LinuxTraceBack] Backtrace(s) for program '%s' : sigcontextptr=%p\n", 
      pfx,ts,a_out ? a_out : "/dev/null", sigcontextptr);

  if (++recur > 1) {
    fprintf(stderr,
      "%s %s [LinuxTraceBack] I don't handle recursive calls very well (recursion level = %d)\n",
      pfx,ts,recur);
    if (recur > 10) {
      fprintf(stderr,"%s %s [LinuxTraceBack] Recursion too deep. Exiting immediately with _exit(%d)\n",
        pfx,ts,recur);
      fflush(NULL);
      _exit(recur); /* Exit immediately */
    }
  }
  if (!sigcontextptr_given) goto finish;

  gdb_trbk();
  dbx_trbk();

 finish:
  fprintf(stderr,"%s %s [%s@%s:%d] End of backtrace(s)\n",pfx,ts,FFL);
  recur--;
}
#else
void
LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr)
{
  int sigcontextptr_given = sigcontextptr ? 1 : 0;
  static int recur = 0;
#if (defined(__GNUC__) || defined(__PGI))
  ucontext_t ctx;
  if (!sigcontextptr) {
      sigcontextptr = (getcontext(&ctx) == 0) ? &ctx : NULL;
  }
#endif
  const char *pfx   = prefix  ? prefix  : drhook_PREFIX(0);
  const char *ts    = timestr ? timestr : drhook_TIMESTR(0);
  const char *a_out = ec_argv()[0];
  const char *s1    = strlen(pfx) ? " " : "";
  const char *s2    = strlen(ts)  ? " " : ""; 
 
  fprintf(stderr,"%s%s%s%s[LinuxTraceBack] Backtrace(s) for program '%s' : sigcontextptr=%p\n", 
    pfx,s1,ts,s2,a_out ? a_out : "/dev/null", sigcontextptr);

  if (++recur > 1) {
    fprintf(stderr,
      "%s%s%s%s[LinuxTraceBack] I don't handle recursive calls very well (recursion level = %d)\n",
      pfx,s1,ts,s2,recur);
    if (recur > 10) {
      fprintf(stderr,"%s%s%s%s[LinuxTraceBack] Recursion too deep. Exiting immediately with _exit(%d)\n",
        pfx,s1,ts,s2,recur);
      _exit(recur); /* Exit immediately */
    }
  }

#if (defined(__GNUC__) || defined(__PGI))
  if (sigcontextptr) {
    /* To have a desired effect, 
       compile with -g (and maybe -O1 or greater to get some optimization)
       and link with -g -Wl,-export-dynamic */
    char *linuxtrbk_fullpath = getenv("LINUXTRBK_FULLPATH");
    int linuxtrbk_fullpath_on = (linuxtrbk_fullpath && *linuxtrbk_fullpath == '1') ? 1 : 0;
    void *trace[GNUC_BTRACE];
    ucontext_t *uc = (ucontext_t *)sigcontextptr;
    int fd = fileno(stderr);
    int trace_size = backtrace(trace, GNUC_BTRACE);
    char *addr2linecmd = (has_bfd || (access(TOSTR(ADDR2LINE),X_OK) != 0)) ? NULL : prealloc_addr2linecmd;
    char **strings = NULL;
    if (trace_size > 1) {
      /* overwrite sigaction with caller's address */
#ifdef __powerpc64__
      trace[1] = uc ? (void *) uc->uc_mcontext.regs->nip : NULL;   // Trick from PAPI_overflow() 
#elif defined(__x86_64__) && defined(REG_RIP) // gcc specific
      trace[1] = uc ? (void *) uc->uc_mcontext.gregs[REG_RIP] : NULL; // RIP: x86_64 specific ; only available in 64-bit mode */
#elif defined(__i386__) && defined(REG_EIP) // gcc specific
      trace[1] = uc ? (void *) uc->uc_mcontext.gregs[REG_EIP] : NULL; // EIP: x86 specific ; only available in 32-bit mode */
#endif
    }
    strings = backtrace_symbols(trace, trace_size);
    fprintf(stderr,"%s%s%s%s[LinuxTraceBack] Backtrace (size = %d) with %s\n",
        pfx,s1,ts,s2,trace_size,
        addr2linecmd ? "addr2line-cmd" : has_bfd ? "BFD-method" : "plain hex-dump");
    if (strings && trace_size > 0) {
      int i;
      FILE *fp = NULL;
      if (addr2linecmd) {
        /* Use ADDR2LINE to obtain source file & line numbers for each trace-address */
        //snprintf(addr2linecmd, len_addr2linecmd, "%s -fs -e '%s'", TOSTR(ADDR2LINE), a_out);
        strcpy(addr2linecmd,TOSTR(ADDR2LINE));
        strcat(addr2linecmd," -fs -e '");
        strcat(addr2linecmd,a_out);
        strcat(addr2linecmd,"'");
        for (i = 0; i < trace_size; i++) {
          char s[30];
          if (trace[i])
            snprintf(s,sizeof(s)," %p",trace[i]);
          else
            snprintf(s,sizeof(s)," 0x0");
          strcat(addr2linecmd,s);
        }
        if (getenv("LD_PRELOAD")) {
          static char ld_preload[] = "LD_PRELOAD=";
          putenv(ld_preload);
        }
        fprintf(stderr,"%s%s%s%s[LinuxTraceBack] %s\n", pfx,s1,ts,s2,addr2linecmd);
        fp = popen(addr2linecmd,"r");
        /* free(addr2linecmd); */
      }
      if (fp) {
        int ndigits = (trace_size > 0) ? 1 + (int)log10(trace_size) : 0;
        extern char *cxxdemangle(const char *mangled_name, int *status); // cxxdemangle.cc (C++ code) : returned string must be free'd
        for (i = 0; i < trace_size; i++) {
          int ok = 0;
          char func[LINELEN];
          if (!feof(fp) && fgets(func, LINELEN, fp)) {
            char line[LINELEN];
            if (!feof(fp) && fgets(line, LINELEN, fp)) {
              char *cxxfunc = NULL;
              char *nl = strchr(func,'\n');
              char *leftB, *plus;
              const char *last_slash = linuxtrbk_fullpath_on ? NULL : strrchr(strings[i],'/');
              if (last_slash) last_slash++; else last_slash = strings[i];
              if (nl) *nl = '\0';
              cxxfunc = cxxdemangle(func,NULL);
              nl = strchr(line,'\n');
              if (nl) *nl = '\0';
              leftB = strchr(last_slash,'(');
              plus = strrchr(last_slash,'+');
              if (leftB && plus && (int)(plus-leftB) > 1) {
                int istat = 0;
                char *cxx = NULL;
                char *therest = plus + 1;
                *plus = '\0';
                cxx = cxxdemangle(leftB + 1,&istat);
                if (cxx) *leftB = '\0';
                fprintf(stderr, "%s%s%s%s[LinuxTraceBack] [%*.*d]: %s%s%s+%s : %s%s at %s\n", 
                  pfx,s1,ts,s2, ndigits, ndigits, i,
                  last_slash,
                  cxx ? "(" : "",
                  cxx ? cxx : "",
                  therest,
                  cxxfunc ? cxxfunc : func,
                  cxxfunc ? "" : "()",
                  line);
                if (cxx) free(cxx);
              }
              else {
                fprintf(stderr, "%s%s%s%s[LinuxTraceBack] [%*.*d]: %s : %s%s at %s\n", 
                  pfx,s1,ts,s2, ndigits, ndigits, i,
                  last_slash,
                  cxxfunc ? cxxfunc : func,
                  cxxfunc ? "" : "()",
                  line);
              }
              if (cxxfunc) free(cxxfunc);
              ok = 1;
            }
          }
          if (!ok) {
              char *cxx = cxxdemangle(strings[i],NULL);
            fprintf(stderr, "%s%s%s%s[LinuxTraceBack] [%*.*d]: %s\n", 
              pfx,s1,ts,s2, ndigits, ndigits, i, 
              cxx ? cxx : strings[i]);
            if (cxx) free(cxx);
          }
        } /* for (i = 0; i < trace_size; i++) */
        fflush(stderr);
        pclose(fp);
      } /* if (fp) */
      else {
        int ndigits = (trace_size > 0) ? 1 + (int)log10(trace_size) : 0;
        InitBFD();
        for (i = 0 ; i < trace_size; ++i) {
          BFD_t b;
          const char *last_slash = linuxtrbk_fullpath_on ? NULL : strrchr(strings[i],'/');
          if (last_slash) last_slash++; else last_slash = strings[i];
          if (ResolveViaBFD(trace[i], &b, last_slash) == 0 /*success*/) {
            fprintf(stderr,"%s%s%s%s[LinuxTraceBack] [%*.*d]: %s : %s() at %s:%u\n",
              pfx,s1,ts,s2,ndigits,ndigits,i,
              last_slash,b.func,b.file,b.lineno);
          }
          else {
            fprintf(stderr,"%s%s%s%s[LinuxTraceBack] [%*.*d]: %s : %p\n",
              pfx,s1,ts,s2,ndigits,ndigits,i,
              last_slash,trace[i]);
          }
        }
      }
    }
    else {
      /* Print traceback directly to fd=2 (stderr) */
      backtrace_symbols_fd(trace, trace_size, fd);
    } /* if (addr2linecmd) else ... */
    if (strings) free(strings); /* Could we live without this free() ? */
  }
#endif /* __GNUC__ */

  if (!sigcontextptr_given) goto finish;

  gdb_trbk();
  dbx_trbk();

finish:
  fprintf(stderr,"%s%s%s%s[LinuxTraceBack] End of backtrace(s)\n",pfx,s1,ts,s2);
  recur--;
}
#endif
 
void linux_trbk_(void)
{
  LinuxTraceBack(NULL,NULL,NULL);
}

#else

/* Non-Linux: A dummy call which does nothing */

void LinuxTraceBack(const char *prefix, const char *timestr, void *sigcontextptr) { }

void linux_trbk_(void) { }

#endif

void linux_trbk(void)
{
  linux_trbk_();
}

/* GNU-debugger traceback */

#if !defined(GNUDEBUGGER)
#define GNUDEBUGGER /usr/bin/gdb
#endif

void gdb_trbk_()
{
  char *gdb = getenv("GNUDEBUGGER");
  if (gdb && 
      (access(TOSTR(GNUDEBUGGER),X_OK) == 0) && /* GNUDEBUGGER was set */
      (strequ(gdb,"1")    || 
       strequ(gdb,"true") || 
       strequ(gdb,"TRUE"))) {
    char gdbcmd[65536];
    pid_t pid = getpid();
    const char *a_out = ec_argv()[0];
    fprintf(stderr,
      "[gdb_trbk] : Invoking %s ...\n",
      TOSTR(GNUDEBUGGER));
    snprintf(gdbcmd,sizeof(gdbcmd),
       "set +eux; %s -batch -n -q -ex 'thread apply all bt' %s %ld < /dev/null",
       TOSTR(GNUDEBUGGER), a_out, (long int)pid);
    
    /* fprintf(stderr,"%s\n",gdbcmd); */
    { int idummy = system(gdbcmd); }
  }
}

void gdb_trbk() { gdb_trbk_(); }


/* DBX-debugger traceback */

#if !defined(DBXDEBUGGER)
#define DBXDEBUGGER /usr/bin/dbx
#endif

void dbx_trbk_()
{
  char *dbx = getenv("DBXDEBUGGER");
  if (dbx && 
      (access(TOSTR(DBXDEBUGGER),X_OK) == 0) && /* DBXDEBUGGER was set */
      (strequ(dbx,"1")    || 
       strequ(dbx,"true") || 
       strequ(dbx,"TRUE"))) {
    pid_t pid = getpid();
    const char *a_out = ec_argv()[0];
    char dbxcmd[65536];
    const char *qopt = "";
    fprintf(stderr,
        "[dbx_trbk] : Invoking %s ...\n",
        TOSTR(DBXDEBUGGER));
    if (a_out && (access(a_out,X_OK|R_OK) == 0)) {
      snprintf(dbxcmd,sizeof(dbxcmd),
           "set +e; /bin/echo 'where; quit; '"
           " | %s%s %s %d ",
           TOSTR(DBXDEBUGGER), qopt, a_out, pid);
    }
    else {
      snprintf(dbxcmd,sizeof(dbxcmd),
           "set +e; /bin/echo 'where; quit; '"
           " | %s%s - %d ",
           TOSTR(DBXDEBUGGER), qopt, pid);
    }
    
    /* fprintf(stderr,"%s\n",dbxcmd); */
    { int idummy = system(dbxcmd); }
  }
}

void dbx_trbk() { dbx_trbk_(); }

#ifdef BFDLIB
static bfd *abfd = 0;
static asymbol **syms = 0;
static asection *text = 0;
#endif

static void InitBFD()
{
#ifdef BFDLIB
  if (!abfd) {
    const char *a_out = ec_argv()[0];
 
    bfd_init();
 
    abfd = bfd_openr(a_out, 0);
    if (!abfd) {
      perror("bfd_openr failed: ");
      return;
    }
 
    bfd_check_format(abfd,bfd_object);
 
    unsigned int storage_needed = bfd_get_symtab_upper_bound(abfd);
    syms = (asymbol **) malloc(storage_needed);
    unsigned int cSymbols = bfd_canonicalize_symtab(abfd, syms);
 
    text = bfd_get_section_by_name(abfd, ".text");
  }
#endif
}
 
static int ResolveViaBFD(void *address, BFD_t *b, const char *str)
{
  int rc = -1;
#ifdef BFDLIB
  if (!abfd) InitBFD();
  if (b) {
    long offset = (long)(address - text->vma);
    if (offset > 0) {
      memset(b,0x0,sizeof(*b));
      if (bfd_find_nearest_line(abfd, text, syms, offset, &b->file, &b->func, &b->lineno) && b->file && b->func) {
        const char *last_slash = strrchr(b->file,'/');
        if (last_slash) b->file = last_slash + 1;
        rc = 0;
      }
      else if (str) {
        static const char qmarks[] = "??";
        if (!b->func) {
          static char *s = NULL; // Not thread safe
          char *loc_lbr;
          if (s) free(s);
          s = strdup(str);
          loc_lbr = strchr(s,'(');
          if (loc_lbr) {
            char *loc_add = NULL;
            *loc_lbr++ = 0;
            loc_add = strchr(loc_lbr,'+');
            if (loc_add) *loc_add = 0;
            b->func = loc_lbr;
          }
          else {
            b->func = qmarks;
          }
        }
        if (!b->file) b->file = qmarks;
        rc = 0;
      } /* if (bfd_find_nearest_line(...)) else if (str) ... */
    }
  }
#endif
  return rc;
}

