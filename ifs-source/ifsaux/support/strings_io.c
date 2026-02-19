/**
* Copyright 1981-2022 ECMWF.
*
* This software is licensed under the terms of the Apache Licence
* Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation
* nor does it submit to any jurisdiction.
*/


/**

This file describes a number of subroutines which can be called
from FORTRAN to handle a text file:

STRINGS_IO_OPEN          to open the file
STRINGS_IO_CLOSE         to close the file
STRINGS_IO_READ          to read from the file
STRINGS_IO_WRITE         to write to the file
STRINGS_IO_FLUSH         to flush the file

The subroutines are written in C and use standard library functions
for file handling (fopen, fclose, fseek, ftell, fread and fwrite).

STRINGS_IO_OPEN
===============

A subroutine which can be called from FORTRAN to open a text file and return 
a suitable unit number for use in calls to STRINGS_IO_READ, STRINGS_IO_WRITE
and STRINGS_IO_CLOSE.

where:

Input parameters are CHARACTERs:
--------------------------------

FILENAME =      a character string describing the file

MODE     =      a character string describing the mode of
                access to the file:
                        r       for read
                        w       for write
                        a       for append



Output parameters are INTEGERs:
-------------------------------

KUNIT =         unit number for the file - it is a C FILE pointer
                and not a FORTRAN unit number.

KRET     =      -1 = Could not open file.
                -2 = Invalid file name.
                -3 = Invalid open mode specified
                 0 = OK.

STRINGS_IO_CLOSE
==============

A subroutine which can be called from FORTRAN to close a text
file previously opened with STRINGS_IO_OPEN.

The format and arguments for the subroutine are as follows:

SUBROUTINE STRINGS_IO_CLOSE(KUNIT,KRET)

where:


Input parameter is an INTEGER:
------------------------------

KUNIT =         unit number for the file; this must have been
                obtained by calling STRINGS_IO_OPEN (see below) - it is
                a C FILE pointer and not a FORTRAN unit number.



Output parameter is an INTEGER:
-------------------------------

KRET     =      -1 error in handling the file.
                 0 = OK.

STRINGS_IO_READ
=============

A subroutine which can be called from FORTRAN to read a line from a 
text file.

The format and arguments for the subroutine are as follows:

SUBROUTINE STRINGS_IO_READ(KUNIT,CDSTRING,KOUNT,KRET)

where:


Input parameters are INTEGERs:
------------------------------

KUNIT =         unit number for the file; this must have been
                obtained by calling STRINGS_IO_OPEN (see below) - it
                is a C FILE pointer and not a FORTRAN unit number.

KOUNT   =       number of CHARACTERS to read from the file.



Output parameters are a CHARACTER and an INTEGER:
------------------------------

CDSTRING =        a CHARACTER array to accept the text from the read.

KRET    =        -2      if there is an error in handling the file

                 -1      if end-of-file is encountered
                         (Note that EOF does not cause a program fail, 
                         so this value must be explicitly caught by 
                         the caller to avoid looping at the EOF)

                 >= 0    number of CHARACTERs read from the file.



STRINGS_IO_WRITE
==============

A subroutine which can be called from FORTRAN to write a line 
from a text file.

The format and arguments for the subroutine are as follows:

SUBROUTINE STRINGS_IO_WRITE(KUNIT,KARRAY,KOUNT,KRET)

where:


Input parameters are INTEGERs:
------------------------------

KUNIT =         unit number for the file; this must have been
                obtained by calling STRINGS_IO_OPEN (see below) - it is
                a C FILE pointer and not a FORTRAN unit number.

CDARRAY =       a CHARACTER array holding the text to write.

KOUNT   =       number of STRINGS to write to the file.



Output parameter is an INTEGER:
-------------------------------

KRET    =       -1      if there is an error in writing to the file

                >= 0    number of CHARACTERs written to the file.
*/

/*
// strings_io.c
*/
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <errno.h>

#ifdef PTHREADS
#include <pthread.h>
#endif

static FILE** fptable = NULL;
static int fptableSize = 0;

#ifdef PTHREADS
static pthread_mutex_t fpTableBusy = PTHREAD_MUTEX_INITIALIZER;
#endif
#define BUFFLEN 4096

/*
// Default buffer size for I/O operations (set via setvbuf)
*/
#define SIZE BUFSIZ
static long size = SIZE;
static int sizeSet = 0;
static char * envSize;
static char** fileBuffer = NULL;

/*
// Debug flags.
*/
#define DEBUGOFF 1
#define DEBUG (debugSet > DEBUGOFF )
static char * debugLevel;
static int debugSet = 0;

#define NAMEBUFFLEN 256
#define MODEBUFFLEN 10

#define CURRENT_FILE (fptable[*unit])

/*
//------------------------------------------------------------------------
// STRINGS_IO_OPEN - Open file (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_strings_io_open_( int* unit, char* name, char* mode, int* iret,
                       int l1, int l2 ) {

/*
// Purpose:
//  Opens file, returns the index of a UNIX FILE pointer held in
//  an internal table (fptable).
//
// First time through, reads value in environment variable STRINGS_IO_BUFSIZE
// (if it is set) and uses it as the size to be used for internal file
// buffers; the value is passed to setvbuf. If STRINGS_IO_BUFSIZE is not set,
// a default value is used.
//
// Function  accepts:
//    name = filename
//    mode = r, r+, w
//
//    Note: l1 and l2 are the lengths of the FORTRAN character strings
//          in name and mode, and should not be passed.
//
// Function returns:
//   INTEGER iret:
//     -1 = Could not open file.
//     -2 = Invalid file name.
//     -3 = Invalid open mode specified
//      0 = OK.
*/
int n;
char *p;
char  flags[4];

char namebuff[NAMEBUFFLEN+1], modebuff[MODEBUFFLEN+1];

/*
// See if DEBUG switched on.
*/
    if( ! debugSet ) {
      debugLevel = getenv("STRINGS_IO_DEBUG");
      if( debugLevel == NULL )
        debugSet = DEBUGOFF;              /* off */
      else {
        int loop;
        for( loop = 0; loop < strlen(debugLevel) ; loop++ ) {
          if( ! isdigit(debugLevel[loop]) ) {
            printf("Invalid number string in STRINGS_IO_DEBUG: %s\n", debugLevel);
            printf("STRINGS_IO_DEBUG must comprise only digits [0-9].\n");
            debugSet = DEBUGOFF;
          }
        }
        debugSet = DEBUGOFF + atol( debugLevel );

      }
      if( DEBUG ) printf("STRINGS_IO_OPEN: debug switched on\n");
    }

/*
// Put the character strings into buffers and ensure that there is a
// null terminator (for SGI case when FORTRAN CHARACTER variable is full
// right to end with characters
*/
    {
      int n1, n2;
      n1 = (l1>NAMEBUFFLEN) ? NAMEBUFFLEN : l1;
      n2 = (l2>MODEBUFFLEN) ? MODEBUFFLEN : l2;

      strncpy( namebuff, name, n1);
      strncpy( modebuff, mode, n2);
      namebuff[n1] = '\0';
      modebuff[n2] = '\0';
    }

    strcpy(flags,"");

    /* *unit = (int) NULL;  sami bug fix */
    *unit = 0;
    *iret = 0;

/*
// Strip trailing blanks
*/
    p  = namebuff + strlen(namebuff) - 1 ;
    while(*p == ' ') {
      *p = 0;
      p--;
    }

    if( DEBUG ) printf("STRINGS_IO_OPEN: filename = [%s]\n", namebuff);
/*
// Build open flags from "modes"
*/
    p = modebuff;

    switch(*p) {

      case 'a':
      case 'A': strcat(flags, "a");
                      break;

      case 'c':
      case 'C':
      case 'w':
      case 'W': strcat(flags, "w");
                break;

      case 'r':
      case 'R':
                if( *(p+1) == '+' )
                  strcat(flags, "r+");
                else
                  strcat(flags, "r");
                break;

      default:  *iret = -3;
                return;

    }
    if( DEBUG ) printf("STRINGS_IO_OPEN: file open mode = %s\n", flags);

/*
// Look for a free slot in fptable.
// (Create the table the first time through).
*/
#ifdef PTHREADS
/*
// Wait if another thread opening a file
*/
    pthread_mutex_lock(&fpTableBusy);
#endif

    n = 0;
    if( fptableSize == 0 ) {
      int i;
      fptableSize = 2;
      fptable = (FILE **) malloc(fptableSize*sizeof(FILE *));
      if( fptable == NULL ) {
        perror("Unable to allocate space for table of FILE pointers");
        exit(1);
      }

      fileBuffer = (char **) malloc(fptableSize*sizeof(char *));
      if( fileBuffer == NULL ) {
        perror("Unable to allocate space for FILE buffers");
        exit(1);
      }

      for( i = 0; i < fptableSize; i++ ) {
        fptable[i] = 0;
        fileBuffer[i] = NULL;
      }
    }
    else {
      while( n < fptableSize ) {
        if(fptable[n]==0) {
          *unit = n;
          break;
        }
        n++;
      }
    }
/*
// If the table overflows, double its size.
*/
    if( n == fptableSize) {
      int i;
      fptableSize = 2*fptableSize;
      fptable = (FILE **) realloc(fptable, fptableSize*sizeof(FILE *));
      if( fptable == NULL ) {
        perror("Unable to reallocate space for table of FILE pointers");
        exit(1);
      }
      n = fptableSize/2;

      fileBuffer = (char **) realloc(fileBuffer, fptableSize*sizeof(char *));
      if( fileBuffer == NULL ) {
        perror("Unable to allocate space for FILE buffers");
        exit(1);
      }

      n = fptableSize/2;
      for( i = n; i < fptableSize; i++ ) {
        fptable[i] = 0;
        fileBuffer[i] = NULL;
      }

      *unit = n;
    }

    if( DEBUG ) printf("STRINGS_IO_OPEN: fptable slot = %d\n", *unit);

    if( DEBUG ) printf("STRINGS_IO_OPEN: using fopen\n");
    fptable[n] = fopen(namebuff, flags );

    if(fptable[n] == NULL) {
      perror(namebuff);
      *iret = -1;
#ifdef PTHREADS
      pthread_mutex_unlock(&fpTableBusy);
#endif
      return;
    }

/*
// Now allocate a buffer for the file, if necessary.
*/
    if( ! sizeSet ) {
      envSize = getenv("STRINGS_IO_BUFSIZE");
      if( envSize == NULL )
        size = SIZE;             /* default */
      else {
        int loop;
        for( loop = 0; loop < strlen(envSize) ; loop++ ) {
          if( ! isdigit(envSize[loop]) ) {
            printf("Invalid number string in STRINGS_IO_BUFSIZE: %s\n", envSize);
            printf("STRINGS_IO_BUFSIZE must comprise only digits [0-9].\n");
            exit(1);
          }
        }
        size = atol( envSize );
      }
      if( size <= 0 ) {
        printf("Invalid buffer size in STRINGS_IO_BUFSIZE: %s\n", envSize);
        printf("Buffer size defined by STRINGS_IO_BUFSIZE must be positive.\n");
        exit(1);
      }
      sizeSet = 1;
    }

    if( DEBUG ) printf("STRINGS_IO_OPEN: file buffer size = %ld\n", size);

    if( fileBuffer[n] == NULL ) {
      fileBuffer[n] = (char *) malloc(size);
    }

    if( setvbuf(CURRENT_FILE, fileBuffer[*unit], _IOFBF, size) ) {
        perror("setvbuf failed");
        *iret = -1;
    }

#ifdef PTHREADS
    pthread_mutex_unlock(&fpTableBusy);
#endif
}

void c_strings_io_open__( int* unit, char* name, char* mode, int* iret, int l1, int l2 ) {
  c_strings_io_open_(unit,name,mode,iret,l1,l2);
}
void c_strings_io_open( int* unit, char* name, char* mode, int* iret, int l1, int l2 ) {
  c_strings_io_open_(unit,name,mode,iret,l1,l2);
}

/*
//------------------------------------------------------------------------
//  STRINGS_IO_READ - Read (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_strings_io_read_(int* unit,char* buffer,int* nchars,int* iret) {
/*
// Purpose:
//  Reads a line of text from a file..
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    nchars = number of chars to read.
//
//  Returns:
//    iret:
//      -2 = error in reading file,
//      -1 = end-of-file,
//      otherwise, = number of chars read.
*/

    if( DEBUG ) {
      printf("STRINGS_IO_READ: fptable slot = %d. ", *unit);
      printf("Number of chars to read = %d\n", *nchars);
    }
    if( fgets(buffer, *nchars, CURRENT_FILE) == NULL ) {
      if( ! feof(CURRENT_FILE) ) {
        *iret = -2;             /*  error in file-handling  */
        perror("strings_io_read");
        clearerr(CURRENT_FILE);
        return;
      }
      else {
        *iret = -1;             /*  end-of-file */
        clearerr(CURRENT_FILE);
      }
    }
    else {
      *iret = strlen(buffer);
    }

    if( DEBUG ) {
      printf("STRINGS_IO_READ: fptable slot = %d. ", *unit);
      printf("Number of chars read = %d\n", *iret);
      printf("STRINGS READ: %s\n",buffer);
    }

    return;
}

void c_strings_io_read__(int* unit,char* buffer,int* nchars,int* iret) {
  c_strings_io_read_(unit,buffer,nchars,iret);
}

void c_strings_io_read(int* unit,char* buffer,int* nchars,int* iret) {
  c_strings_io_read_(unit,buffer,nchars,iret);
}

/*
//------------------------------------------------------------------------
//  STRINGS_IO_WRITE - Write (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_strings_io_write_(int* unit,char* buffer,int* nchars,int* iret) {
/*
// Purpose:
//  Writes a line of text to a file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    nchars = number of chars to write.
//
//  Returns:
//    iret:
//      -1 = Could not write to file.
//     >=0 = Number of chars written.
*/
    if( DEBUG ) {
      printf("STRINGS_IO_WRITE: fptable slot = %d. ", *unit);
    }

    if( fprintf(CURRENT_FILE, "%.*s\n", *nchars, buffer) < 0 ) {
      perror("strings_io_write");
      *iret = -1;
    }
    else {
      *iret = 0;
    }

    if( DEBUG ) {
      printf("STRINGS_IO_WRITE: fptable slot = %d. ", *unit);
      printf("STRINGS_IO_WRITE: number of chars written = %d\n", *iret+1);
    }

    return;
}

void c_strings_io_write__(int* unit,char* buffer,int* nchars,int* iret) {
  c_strings_io_write_(unit,buffer,nchars,iret);
}

void c_strings_io_write(int* unit,char* buffer,int* nchars,int* iret) {
  c_strings_io_write_(unit,buffer,nchars,iret);
}

/*
//------------------------------------------------------------------------
//  STRINGS_IO_FSYNC - fsync
//------------------------------------------------------------------------
*/
void c_strings_io_fsync(int fd, int* iret) {
    // Same implementation of see eckit::fsync
    *iret = fsync(fd);
    while (*iret < 0 && errno == EINTR) {
        *iret = fsync(fd);
    }
    return;
}


/*
//------------------------------------------------------------------------
//  STRINGS_IO_FLUSH - flush + fsync
//------------------------------------------------------------------------
*/
void c_strings_io_flush_(int* unit, int* iret) {
/*
// Purpose:     Flushes file.
*/
    if( DEBUG )
      printf("STRINGS_IO_FLUSH: fptable slot = %d\n", *unit);

    // Implementation matching eckit FileHandle

    // flush
    if( (*iret = fflush(CURRENT_FILE)) != 0) {
      perror("strings_io_flush: fflush failed");
      return;
    }

    // fsync
    c_strings_io_fsync(fileno(CURRENT_FILE),iret);
    while (*iret < 0 && errno == EINTR) {
      c_strings_io_fsync(fileno(CURRENT_FILE),iret);
    }
    if (*iret < 0) {
      perror("strings_io_flush: Cannot fsync");
      return;
    }

    *iret = 0;
    return;
}

void c_strings_io_flush__(int* unit, int* iret) {
  c_strings_io_flush_(unit,iret);
}

void c_strings_io_flush(int* unit, int* iret) {
  c_strings_io_flush_(unit,iret);
}

/*
//------------------------------------------------------------------------
//   STRINGS_IO_CLOSE - close (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_strings_io_close_(int* unit,int* iret) {
/*
// Purpose:
//  Closes file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
////  Returns:
//    iret:
//      0 = OK.
//      otherwise = error in handling file.
*/
    *iret = 0;
    if( DEBUG ) {
      printf("STRINGS_IO_CLOSE: fptable slot = %d\n", *unit);
    }

    if( !CURRENT_FILE ) {
      printf("WARNING: strings_io_close: File (fptable slot = %d) was already closed.\n", *unit);
      return;
    }

    // Flush before closing
    c_strings_io_flush(unit,iret);
    if( *iret != 0 ) {
      perror("strings_io_close: could not flush");
      return;
    }

    // Close
    if( ( *iret = fclose(CURRENT_FILE) ) != 0 ) {
      perror("strings_io_close");
      return;
    }

    CURRENT_FILE = 0;
    return;
}

void c_strings_io_close__(int* unit,int* iret) {
  c_strings_io_close_(unit,iret);
}

void c_strings_io_close(int* unit,int* iret) {
  c_strings_io_close_(unit,iret);
}


