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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>

void spawn_ (const char * cmd, int * kret, const unsigned long cmd_len)
{
  char _cmd[cmd_len+1];
  pid_t pid;

  *kret = 0;


  memset (_cmd, 0, sizeof (_cmd));
  strncpy (_cmd, cmd, cmd_len);

  fprintf (stderr, "SPAWN: %s\n", _cmd);

  /* We use vfork because it requires 
   * less memory than plain fork */

  errno = 0;

#if defined(LINUX) || defined(__APPLE__)
  pid = vfork ();

  if (pid < 0)
    {
      perror ("SPAWN");
      goto error;
    }
  else if (pid > 0)
    {
      int status;
      waitpid (pid, &status, 0);
      if (! WIFEXITED (status))
        goto error;
      if (WEXITSTATUS (status))
        goto error;
    }
  else
    {
      execl (_cmd, _cmd, (char *)NULL);
      _exit (1);
    }
#else

  if (system (_cmd) != 0)
    goto error;

#endif

  return;

error:

  *kret = 1;

  return;
}




