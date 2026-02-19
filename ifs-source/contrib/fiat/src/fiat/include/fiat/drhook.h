/*
 * (C) Copyright 2005- ECMWF.
 * (C) Copyright 2013- Meteo-France.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef _DRHOOK_H_
#define _DRHOOK_H_

#include "fiat_pp.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef GNUC_BTRACE
#define GNUC_BTRACE 128
#endif

extern int drhook_lhook;

/* drhook.c external interfaces */

extern void
c_drhook_getenv_(const char *s,
                 char *value,
                 /* Hidden arguments */
                 int slen,
                 const int valuelen);

extern void
c_drhook_memcounter_(const int *thread_id,
                     const long long int *size,
                     long long int *keyptr_addr);

extern void
c_drhook_raise_(const int *sig);

extern void
c_drhook_print_(const int *ftnunitno,
                const int *thread_id,
                const int *print_option, /* 
                                            1=raw call counts 
                                            2=calling tree
                                            3=profiling info
                                         */
                int *level);

extern void
c_drhook_init_signals_(const int *enforce);

extern void
c_drhook_set_lhook_(const int *lhook);

extern void 
c_drhook_init_(const char *progname,
               const int *num_threads
               /* Hidden length */
               ,int progname_len);

extern void
c_drhook_start_(const char *name, 
                const int *thread_id, 
                double *key,
                const char *filename,
                const int *sizeinfo
                /* Hidden length */
                ,int name_len, int filename_len);

extern void
c_drhook_end_(const char *name,
              const int *thread_id,
              const double *key,
              const char *filename,
              const int *sizeinfo
              /* Hidden length */
              ,int name_len, int filename_len);

extern void
c_drhook_watch_(const int *onoff,
		const char *array_name,
		const void *array_ptr,
		const int *nbytes,
		const int *abort_if_changed,
		const int *printkey,
		const int *nvals,
		const int *print_traceback_when_set
		/* Hidden length */
		,int array_name_len);

extern void
c_drhook_check_watch_(const char *where,
		      const int *allow_abort
		      /* Hidden length */
		      , int where_len);

/* see dr_hook_prt.F90 for below */
extern void dr_hook_prt_logical_(const int *ftnunitno, const void *ptr, const int *nmax);
extern void dr_hook_prt_char_(const int *ftnunitno, const void *ptr, const int *nmax);
extern void dr_hook_prt_i4_(const int *ftnunitno, const void *ptr, const int *nmax);
extern void dr_hook_prt_i8_(const int *ftnunitno, const void *ptr, const int *nmax);
extern void dr_hook_prt_r4_(const int *ftnunitno, const void *ptr, const int *nmax);
extern void dr_hook_prt_r8_(const int *ftnunitno, const void *ptr, const int *nmax);

extern void ec_meminfo_(const int *KU, const char *CDSTRING, 
			const int *KCOMM, const int *KBARR, 
			const int *KIOTASK, const int *KCALL,
			int len_CDSTRING); /* see ec_meminfo.F90 */

/* see drhook.c */
extern const char *drhook_TIMESTR(int tid);
extern const char *drhook_PREFIX(int tid);

/**** C-interface to Dr.Hook ****/

extern void drhook_init(int argc, char *argv[]);

extern int drhook_active( void );

extern void
Dr_Hook(const char *name, int option, double *handle, 
        const char *filename, int sizeinfo,
        int name_len, int filename_len);

#define DRHOOK_START_RECUR(name,recur) \
  static const char *drhook_name = #name; \
  static const int drhook_name_len = sizeof(#name) - 1; /* Compile time eval */ \
  static const char *drhook_filename = __FILE__; \
  static const int drhook_filename_len = sizeof(__FILE__) - 1; /* Compile time eval */ \
  double zhook_handle; \
  if (!recur && drhook_lhook) Dr_Hook(drhook_name, 0, &zhook_handle, \
                                      drhook_filename, 0, \
                                      drhook_name_len, drhook_filename_len); {


#define DRHOOK_START_BY_STRING_RECUR(name, recur) \
  static const char *drhook_name = name; \
  static const int drhook_name_len = sizeof(name) - 1; /* Compile time eval */ \
  static const char *drhook_filename = __FILE__; \
  static const int drhook_filename_len = sizeof(__FILE__) - 1; /* Compile time eval */ \
  double zhook_handle; \
  if (!recur && drhook_lhook) Dr_Hook(drhook_name, 0, &zhook_handle, \
                                      drhook_filename, 0, \
                                      drhook_name_len, drhook_filename_len); {

#define DRHOOK_START_BY_STRING(name) DRHOOK_START_BY_STRING_RECUR(name,0)

#define DRHOOK_RETURN_RECUR(sizeinfo,recur) \
  if (!recur && drhook_lhook) Dr_Hook(drhook_name, 1, &zhook_handle, \
                                      drhook_filename, sizeinfo, \
                                      drhook_name_len, drhook_filename_len)

#define DRHOOK_RETURN(sizeinfo) DRHOOK_RETURN_RECUR(sizeinfo,0)

#define DRHOOK_END_RECUR(sizeinfo,recur) ; } DRHOOK_RETURN_RECUR(sizeinfo,recur)

#define DRHOOK_END_DEFAULT()          DRHOOK_END_RECUR(0,0)
#define DRHOOK_END_SIZEINFO(sizeinfo) DRHOOK_END_RECUR(sizeinfo,0)

#define DRHOOK_END( ... ) FIAT_PP_CAT( _DRHOOK_END_, FIAT_PP_VARIADIC_SIZE(__VA_ARGS__) )( __VA_ARGS__ )
// With 0 args --> DRHOOK_END_DEFAULT()
// With 1 arg  --> DRHOOK_END_SIZEINFO(sizeinfo)
#define _DRHOOK_END_0    DRHOOK_END_DEFAULT
#define _DRHOOK_END_1    DRHOOK_END_SIZEINFO

#define DRHOOK_START(name) DRHOOK_START_RECUR(name,0)


typedef void (*drhook_abort_t)(const char *file, int line, const char *text);
extern void drhook_set_abort( drhook_abort_t );
extern void drhook_abort( const char *file, int line, const char *txt );

extern void drhook_calltree( void );


/* Fortran routines */

extern void
dr_hook_prt_(const int *ftnunitno,
             const char *s
             /* Hidden arguments */
             , int s_len);

extern void
dr_hook_procinfo_(int *myproc, int *nproc);

#ifdef __cplusplus
} // extern "C"
#endif

#endif  /* _DRHOOK_H_ */

