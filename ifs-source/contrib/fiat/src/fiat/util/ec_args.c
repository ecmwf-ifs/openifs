/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#if defined(LINUX)
#include <linux/limits.h>
#endif
#ifdef __APPLE__
#include <libproc.h>
#include <errno.h>
#endif

#include "ec_args.h"

#define EC_MAX_ARGS 512

typedef struct {
  char *name;
  int len;
} arg_t;

static arg_t *args = NULL;
static const char* * c_argv = NULL;
static int numargs = -1;
static char *cl_terminate = NULL;
static char *a_out = NULL;

#if !defined(PSCMD)
#define PSCMD "/bin/ps"
#endif /* !defined(PSCMD) */

#if !defined(TAILCMD)
#define TAILCMD "/usr/bin/tail"
#endif /* !defined(TAILCMD) */

static const char *get_a_out()
{
  /* progname is a blank string;
     this is most likely due to a Fortran-call to getarg
     from program that has a C-main program, thus Fortran getarg
     may return a blank string */
  
  /* Using an alternative method of getting a.out :
     
  ps -p<pid> | tail -1 | awk '{print $NF}' 

  This can give false positives, as we search the path for the binary: we may have
  executed ./MASTER rather than the first MASTER in the path.

  Naturally this cannot be the nicest method around ;-(
  So, a competition is launched: make this better and 
  you may win a week in Bahamas!!

  */

#if defined(LINUX)
   /* On Linux, the following is more reliable, if /proc is available. 
    */
  if (!a_out) {
    int len;
    char symlink[NAME_MAX], real_exe[NAME_MAX];
    snprintf(symlink, NAME_MAX, "/proc/%d/exe", getpid());
    if ((len = readlink(symlink, real_exe, NAME_MAX)) > 0)  {
      a_out = malloc((len+1)*sizeof(*a_out));
      strncpy(a_out, real_exe, len);
	    a_out[len] = '\0';
    }
  }
#elif defined(__APPLE__)
   /* On macOS, the following is more reliable
    */
  if (!a_out) {
    char path[PROC_PIDPATHINFO_MAXSIZE];
    pid_t pid = getpid();
    int ret = proc_pidpath( pid, path, sizeof(path)); // For name only, use 'proc_name' instead
    if ( ret <= 0 ) {
        fprintf(stderr, "PID %d: proc_pidpath ();\n", pid);
        fprintf(stderr, "    %s\n", strerror(errno));
    } else {
        a_out = strdup(path);
    }
  }
#endif

  if (!a_out && (access(PSCMD,X_OK) == 0)) {
    char cmd[sizeof(PSCMD) + sizeof(TAILCMD) + 100];
    FILE *fp = NULL;
    pid_t pid = getpid();
    sprintf(cmd,"%s -p%d | %s -1 | awk '{if (NF>4) {print $4} else {print $NF}}'", PSCMD, (int)pid, TAILCMD);
    fp = popen(cmd, "r");
    if (fp) {
      char c[65536];
      if (fscanf(fp,"%s",c) == 1) {
        if (!strchr(c,'/')) { 
          /* The file path was NOT embedded in the name 
             ==> Must search from $PATH f.ex. /bin:/usr/bin:/some/thing/else:/etc/bin */
          char *path = getenv("PATH");
          if (path) {
            int lenc = strlen(c);
            char *saved = strdup(path);
            char *start = saved;
            char *token = strtok(saved,":");
            do {
              int lenf = strlen(start) + 1 + lenc + 1;
              char *fullpath = malloc(lenf * sizeof(*fullpath));
              snprintf(fullpath,lenf,"%s/%s",start,c);

              if (access(fullpath,X_OK) == 0) { /* It's this one!! */
                a_out = fullpath;
                break; /* do { ... } while (token) */
              }
              free(fullpath);
              start = token;
              token = strtok(NULL,":");
            } while (token);
            free(saved);
          }
        } /* if (!strchr(c,'/')) */
        if (!a_out) a_out = strdup(c);
      }
      pclose(fp);
    }
  }
  if (!a_out) a_out = strdup("/unknown/executable");

  return a_out;
}

static void reset_argv() {
  int i;
  if( c_argv == NULL ) {
    c_argv = calloc(EC_MAX_ARGS, sizeof(char*));
  }
  for( i=0; i<EC_MAX_ARGS; ++i ) {
    c_argv[i] = NULL;
  }
}

/*
 * Set 0-terminated arguments, including program name, as in C
 */
void ec_args(int argc, char* argv[]) {
  if (numargs == -1 && !args) {
    if (argc > 0) {
      int j;
      args = calloc(argc, sizeof(arg_t));
      reset_argv();
      /* cl_terminate: see ifsaux/module/mpl_arg_mod.F90 */
      if (!cl_terminate) {
        char *env = getenv("MPL_CL_TERMINATE");
        cl_terminate = env ? strdup(env) : strdup("-^");
      }
      numargs = 0;
      for (j=0; j<argc; j++) {
        if (!argv[j] || strcmp(argv[j],cl_terminate) == 0) break;
        args[j].name = strdup(argv[j]);
        args[j].len = strlen(argv[j]);
        c_argv[j] = args[j].name;
        numargs++;
        if ( j==0 ) {
          if (a_out) free(a_out);
          a_out = strdup(args[0].name);
        }
      }
      if (numargs == 0) {
        const char *arg0 = get_a_out();
        args[0].name = strdup(arg0);
        args[0].len = strlen(arg0);
        c_argv[0] = args[j].name;
      }
      else {
        if (a_out) free(a_out);
        a_out = strdup(args[0].name);
        numargs--; /* Fortran # of args == C # of args - 1 */
      }
    } /* if (argc > 0) */
  } /* if (numargs == -1 && !args) */
}

/*
 * Return number of arguments + program name, as in C
 */
int ec_argc(void) { return 1 + numargs; }

// Return array of 0-terminated arguments as in C
const char* const* ec_argv(void)
{
  if( c_argv == NULL ) {
    reset_argv(); // not thread safe
    c_argv[0] = a_out ? a_out : get_a_out();
  }
  return c_argv;
}


/* 
 * Legacy Fortran functions, currently only used within mpl_arg_mod
 * Please do not use! Use "ec_args_mod" instead
 */

int iargc_c_(void) { return numargs; }

void getarg_c_(const int *argno, char *arg
	       /* Hidden argument */
	       , const int arg_len)
{
  int Argno = argno ? *argno : -1;
  int len = 0;
  const char *s = NULL;
  /* Special case : Argno := 0 i.e. the executable name */
  if (arg && arg_len > 0 && Argno == 0) {
    if( a_out ) {
      s = a_out;
    }
    else {
      s = get_a_out();
    }
    len = strlen(s);
  }
  else if ( arg && arg_len > 0 && 
            Argno > 0 && Argno <= numargs &&
            args && args[Argno].name) {
    s = args[Argno].name;
    len = args[Argno].len;
  }
  if (arg && arg_len > 0) memset(arg,' ',arg_len);
  if (s && len > 0) {
    if (arg_len < len) len = arg_len;
    strncpy(arg,s,len);
    if (arg_len > len) memset(&arg[len],' ',arg_len-len);
  }
}

void putarg_c_(const int *argno, const char *arg
	       /* Hidden argument */
	       , int arg_len)
{

  int Argno = argno ? *argno : -1;
  if (arg && arg_len >= 0 && 
      Argno >= 0 && Argno <= numargs && args) {
    char *s = calloc(arg_len+1,sizeof(*s));
    strncpy(s,arg,arg_len);
    s[arg_len] = '\0';
    if (args[Argno].name) free(args[Argno].name);
    args[Argno].name = s;
    args[Argno].len = arg_len;
    c_argv[Argno] = args[Argno].name;
  }
}


void putarg_info_(const int *argc, const char *cterm
		  /* Hidden argument */
		  , int cterm_len)
{
  int Argc = argc ? *argc : 0;
  if (cterm && cterm_len >= 0) {
    if (cl_terminate) free(cl_terminate);
    cl_terminate = calloc(cterm_len+1,sizeof(*cl_terminate));
    strncpy(cl_terminate,cterm,cterm_len);
    cl_terminate[cterm_len] = '\0';
  }
  if (numargs >= 0 || args) {
    if (args) {
      int j;
      for (j=0; j<=numargs; j++) { /* note:  "j<=", not "j<" */
        if (args[j].name) {
          free(args[j].name);
        }
      }
      free(args);
      args = NULL;
    } /* if (args) */
    numargs = -1;
  }
  /* Re-initialize args & numargs */
  if (Argc < 0) Argc = 0;
  numargs = Argc;
  args = calloc(1 + numargs, sizeof(arg_t));
  reset_argv();
}
