/*
 * (C) Copyright 1989- ECMWF.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * 
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction
 * 
 * (C) Copyright 1989- Meteo-France.
 * 
*/

#include "fi_libc.h"

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>


void fi_gettimeofday_ (double * t)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  *t = (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1e6;
}

void fi_fopen_ (FILE ** fpp, const char * path, const char * mode, int path_len, int mode_len)
{
  char path1[path_len+1];
  char mode1[mode_len+1];
  FILE * fp;
  int i;

  for (i = path_len; (path[i-1] == ' ') && (i > 0); i--);
  path_len = i;

  for (i = 0; i < path_len; i++)
    path1[i] = path[i];
  path1[path_len] = '\0';

  for (i = mode_len; (mode[i-1] == ' ') && (i > 0); i--);
  mode_len = i;

  for (i = 0; i < mode_len; i++)
    mode1[i] = mode[i];
  mode1[mode_len] = '\0';

  fp = fopen (path1, mode1);

  *fpp = fp;
}

void fi_fclose_ (fi_integer4 * err, FILE ** fpp)
{
  *err = fclose (*fpp);
}

void fi_fread_ (fi_integer8 * err, void * ptr, fi_integer8 * size, fi_integer8 * nmemb, FILE ** fpp)
{
  *err = fread (ptr, *size, *nmemb, *fpp);
}

void fi_fwrite_ (fi_integer8 * err, const void * ptr, fi_integer8 * size, fi_integer8 * nmemb, FILE ** fpp)
{
  *err = fwrite (ptr, *size, *nmemb, *fpp);
}

void fi_fseek_ (fi_integer4 * err, FILE ** fpp, fi_integer8 * offset, fi_integer4 * whence)
{
  *err = fseek (*fpp, *offset, *whence);
}

void fi_fileno_ (fi_integer4 * err, FILE ** fpp)
{
  *err = fileno (*fpp);
}

void fi_ftell_ (fi_integer8 * err, FILE ** fpp)
{
  *err = ftell (*fpp);
}

void fi_fstat_ (fi_integer4 *err, fi_integer4 *fd, fi_integer8 buf[13])
{
  struct stat buf1;
  *err = fstat (*fd, &buf1);

  buf[0]  = buf1.st_dev;      
  buf[1]  = buf1.st_ino;      
  buf[2]  = buf1.st_mode;     
  buf[3]  = buf1.st_nlink;    
  buf[4]  = buf1.st_uid;      
  buf[5]  = buf1.st_gid;      
  buf[6]  = buf1.st_rdev;     
  buf[7]  = buf1.st_size;     
  buf[8]  = buf1.st_blksize;  
  buf[9]  = buf1.st_blocks;   
  buf[10] = buf1.st_atime;    
  buf[11] = buf1.st_mtime;    
  buf[12] = buf1.st_ctime;    
}

void fi_ftruncate_ (fi_integer4* err, fi_integer4 *fd, fi_integer8 *length)
{
  *err = ftruncate (*fd, *length);
}

void fi_mkdir_ (const char * path, int path_len)
{
  char p[path_len+1];
  memcpy (p, path, path_len);
  p[path_len] = '\0';
  mkdir (p, 0755);
}

void fi_chdir_ (fi_integer4 * ierr, const char * path, int path_len)
{
  char p[path_len+1];
  memcpy (p, path, path_len);
  p[path_len] = '\0';
  *ierr = chdir (p);
}

void fi_rename_ (fi_integer4 * ierr, const char * oldpath, const char * newpath, int oldpath_len, int newpath_len)
{
  char _oldpath[oldpath_len+1];
  char _newpath[newpath_len+1];

  strncpy (_oldpath, oldpath, oldpath_len); _oldpath[oldpath_len] = '\0';
  strncpy (_newpath, newpath, newpath_len); _newpath[newpath_len] = '\0';

  *ierr = rename (_oldpath, _newpath);
}

void fi_fflush_ (fi_integer4* err, FILE **fpp)
{
  *err = fflush (*fpp);
}

void fi_errno_ (fi_integer4* err)
{
  *err = errno;
}

void fi_unlink_ (fi_integer4* err, const char * path, int path_len)
{
  char p[path_len+1];
  memcpy (p, path, path_len);
  p[path_len] = '\0';
  *err = unlink (p);
}

