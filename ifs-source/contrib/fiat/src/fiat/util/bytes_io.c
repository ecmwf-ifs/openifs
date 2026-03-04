/*
 * (C) Copyright 2005- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/**

This file describes a number of subroutines which can be called
from FORTRAN to handle an unblocked binary file:

BYTES_IO_OPEN          to open the file
BYTES_IO_CLOSE         to close the file
BYTES_IO_TELL          to give position in the file (optional)
BYTES_IO_SEEK          to position the file (optional)
BYTES_IO_READ          to read from the file
BYTES_IO_WRITE         to write to the file
BYTES_IO_FLUSH         to flush the file

The subroutines are written in C and use standard library functions
for file handling (fopen, fclose, fseek, ftell, fread and fwrite).

These routines are a direct replacement for:
  - PBOPEN
  - PBCLOSE
  - PBTELL
  - PBSEEK
  - PBREAD
  - PBWRITE
  - PBFLUSH


BYTES_IO_OPEN
=============

A subroutine which can be called from FORTRAN to open an unblocked 
binary file and return a suitable unit number for use in calls to 
BYTES_IO_READ, and BYTES_IO_WRITE, BYTES_IO_SEEK, BYTES_IO_TELL.

The format and arguments for the subroutine are as follows:

SUBROUTINE BYTES_IO_OPEN(KUNIT,FILENAME,MODE,KRET)

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

BYTES_IO_CLOSE
==============

A subroutine which can be called from FORTRAN to close an unblocked 
binary file previously opened with BYTES_IO_OPEN.

The format and arguments for the subroutine are as follows:

SUBROUTINE BYTES_IO_CLOSE(KUNIT,KRET)

where:


Input parameter is an INTEGER:
------------------------------

KUNIT =         unit number for the file; this must have been
                obtained by calling BYTES_IO_OPEN (see below) - it is
                a C FILE pointer and not a FORTRAN unit number.



Output parameter is an INTEGER:
-------------------------------

KRET     =      -1 error in handling the file.
                 0 = OK.



BYTES_IO_SEEK
=============

A subroutine which can be called from FORTRAN to position an 
unblocked binary file at any desired byte position.

The format and arguments for the subroutine are as follows:

SUBROUTINE BYTES_IO_SEEK(KUNIT,KOFFSET,KSTART,KRET)

where:


Input parameters are INTEGERs:
------------------------------

KUNIT =         unit number for the file; this must have been
                obtained by calling BYTES_IO_OPEN (see below) - it is
                a C FILE pointer and not a FORTRAN unit number.

KOFFSET =       number of bytes to offset the file; this is used
                as either an absolute or relative offset depending
                on the value of KSTART.

KSTART  =       0       if KOFFSET is an absolute count from the
                        beginning of the file,

                1       if KOFFSET is a relative offset from the
                        current byte position in the file,

                2       if KOFFSET is an absolute offset from the
                        end of file.


Output parameter is an INTEGER:
-------------------------------

KRET    =       -2      if there is an error in handling the file

                -1      if end-of-file is encountered
                        (Note that EOF does not cause a program fail, 
                        so this value must be explicitly caught by 
                        the caller to avoid looping at the EOF)

                >= 0    byte offset from the start of file after
                        positioning.




BYTES_IO_READ
=============

A subroutine which can be called from FORTRAN to read a block of 
bytes from an unblocked binary file.  ( Note that this routine
behaves differently from BYTES_IO_READ2 on hitting end-of-file. )

The format and arguments for the subroutine are as follows:

SUBROUTINE BYTES_IO_READ(KUNIT,KARRAY,KOUNT,KRET)

where:


Input parameters are INTEGERs:
------------------------------

KUNIT =         unit number for the file; this must have been
                obtained by calling BYTES_IO_OPEN (see below) - it is
                a C FILE pointer and not a FORTRAN unit number.

KOUNT   =       number of BYTES to read from the file.



Output parameters are INTEGERs:
------------------------------

KARRAY =        an INTEGER array to accept the bytes from the read.

KRET    =       -2      if there is an error in handling the file

                -1      if end-of-file is encountered
                        (Note that EOF does not cause a program fail, 
                        so this value must be explicitly caught by 
                        the caller to avoid looping at the EOF)

                >= 0    number of BYTES read from the file.



BYTES_IO_WRITE
==============

A subroutine which can be called from FORTRAN to write a block of 
bytes to an unblocked binary file.

The format and arguments for the subroutine are as follows:

SUBROUTINE BYTES_IO_WRITE(KUNIT,KARRAY,KOUNT,KRET)

where:


Input parameters are INTEGERs:
------------------------------

KUNIT =         unit number for the file; this must have been
                obtained by calling BYTES_IO_OPEN (see below) - it is
                a C FILE pointer and not a FORTRAN unit number.

KARRAY =        an INTEGER array holding the bytes for the write.

KOUNT   =       number of BYTES to write to the file.



Output parameter is an INTEGER:
-------------------------------

KRET    =       -1      if there is an error in writing to the file

                >= 0    number of BYTES written to the file.
*/



/*
// bytes_io.c
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
// BYTES_IO_OPEN - Open file (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_bytes_io_open_( int* unit, char* name, char* mode, int* iret,
                       int l1, int l2 ) {

/*
// Purpose:
//  Opens file, returns the index of a UNIX FILE pointer held in
//  an internal table (fptable).
//
// First time through, reads value in environment variable BYTES_IO_BUFSIZE
// (if it is set) and uses it as the size to be used for internal file
// buffers; the value is passed to setvbuf. If BYTES_IO_BUFSIZE is not set,
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
      debugLevel = getenv("BYTES_IO_DEBUG");
      if( debugLevel == NULL )
        debugSet = DEBUGOFF;              /* off */
      else {
        size_t loop;
        for( loop = 0; loop < strlen(debugLevel) ; loop++ ) {
          if( ! isdigit(debugLevel[loop]) ) {
            printf("Invalid number string in BYTES_IO_DEBUG: %s\n", debugLevel);
            printf("BYTES_IO_DEBUG must comprise only digits [0-9].\n");
            debugSet = DEBUGOFF;
          }
        }
        debugSet = DEBUGOFF + atol( debugLevel );
      }
      if( DEBUG ) printf("BYTES_IO_OPEN: debug switched on\n");
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

    if( DEBUG ) printf("BYTES_IO_OPEN: filename = [%s]\n", namebuff);
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
    if( DEBUG ) printf("BYTES_IO_OPEN: file open mode = %s\n", flags);

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

    if( DEBUG ) printf("BYTES_IO_OPEN: fptable slot = %d\n", *unit);

    if( DEBUG ) printf("BYTES_IO_OPEN: using fopen\n");
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
      envSize = getenv("BYTES_IO_BUFSIZE");
      if( envSize == NULL )
        size = SIZE;             /* default */
      else {
        size_t loop;
        for( loop = 0; loop < strlen(envSize) ; loop++ ) {
          if( ! isdigit(envSize[loop]) ) {
            printf("Invalid number string in BYTES_IO_BUFSIZE: %s\n", envSize);
            printf("BYTES_IO_BUFSIZE must comprise only digits [0-9].\n");
            exit(1);
          }
        }
        size = atol( envSize );
      }
      if( size <= 0 ) {
        printf("Invalid buffer size in BYTES_IO_BUFSIZE: %s\n", envSize);
        printf("Buffer size defined by BYTES_IO_BUFSIZE must be positive.\n");
        exit(1);
      }
      sizeSet = 1;
    }

    if( DEBUG ) printf("BYTES_IO_OPEN: file buffer size = %ld\n", size);

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

void c_bytes_io_open__( int* unit, char* name, char* mode, int* iret, int l1, int l2 ) {
  c_bytes_io_open_(unit,name,mode,iret,l1,l2);
}
void c_bytes_io_open( int* unit, char* name, char* mode, int* iret, int l1, int l2 ) {
  c_bytes_io_open_(unit,name,mode,iret,l1,l2);
}

/*
//------------------------------------------------------------------------
// BYTES_IO_SEEK - Seek (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_bytes_io_seek_(int* unit,int* offset,int* whence,int* iret) {
/*
//
// Purpose:
//   Seeks to a specified location in file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    offset = byte count
//
//    whence  = 0, from start of file
//            = 1, from current position
//            = 2, from end of file.
//
//  Returns:
//    iret:
//      -2 = error in handling file,
//      -1 = end-of-file
//      otherwise,  = byte offset from start of file.
*/
int my_offset = (int) *offset;
int my_whence = (int) *whence;

/*
// Must use negative offset if working from end-of-file
*/
    if( DEBUG ) {
      printf("BYTES_IO_SEEK: fptable slot = %d\n", *unit);
      printf("BYTES_IO_SEEK: Offset = %d\n", my_offset);
      printf("BYTES_IO_SEEK: Type of offset = %d\n", my_whence);
    }

    if( my_whence == 2) my_offset = - abs(my_offset);

    *iret = ftell(CURRENT_FILE);
    if( DEBUG ) printf("BYTES_IO_SEEK: current position = %d\n", *iret);
    if( *iret == my_offset && my_whence == 0)
      *iret = 0;
    else
      *iret = fseek(CURRENT_FILE, my_offset, my_whence);

    if( DEBUG ) printf("BYTES_IO_SEEK: fseek return code = %d\n",*iret);

    if( *iret != 0 ) {
      if( ! feof(CURRENT_FILE) ) {
        *iret = -2;             /* error in file-handling */
        perror("bytes_io_seek");
      }
      else
        *iret = -1;             /* end-of-file  */

      clearerr(CURRENT_FILE);
      return;
    }

/*
// Return the byte offset from start of file
*/
    *iret = ftell(CURRENT_FILE);

    if( DEBUG )
      printf("BYTES_IO_SEEK: byte offset from start of file = %d\n",*iret);

    return;

}

void c_bytes_io_seek(int* unit,int* offset,int* whence,int* iret) {
  c_bytes_io_seek_(unit,offset,whence,iret);
}


/*
//------------------------------------------------------------------------
// BYTES_IO_TELL - Tells current file position (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_bytes_io_tell_(int* unit,int* iret) {
/*
//
// Purpose:
//   Tells current byte offset in file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//  Returns:
//    iret:
//      -2 = error in handling file,
//      otherwise,  = byte offset from start of file.
*/


/*
// Return the byte offset from start of file
*/
    *iret = ftell(CURRENT_FILE);

    if( *iret < 0 ) {
      if( DEBUG ) {           /* error in file-handling */
        printf("BYTES_IO_TELL: fptable slot = %d. ", *unit);
        printf("Error status = %d\n", *iret);
      }
      perror("bytes_io_tell");
      *iret = -2;
    }

    if( DEBUG ) {
      printf("BYTES_IO_TELL: fptable slot = %d. ", *unit);
      printf("Byte offset from start of file = %d\n",*iret);
    }

    return;
}

void c_bytes_io_tell__(int* unit,int* iret) {
  c_bytes_io_tell_(unit,iret);
}

void c_bytes_io_tell(int* unit,int* iret) {
  c_bytes_io_tell_(unit,iret);
}


/*
//------------------------------------------------------------------------
//  BYTES_IO_READ - Read (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_bytes_io_read_(int* unit,char* buffer,int* nbytes,int* iret) {
/*
// Purpose:
//  Reads a block of bytes from a file..
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    nbytes = number of bytes to read.
//
//  Returns:
//    iret:
//      -2 = error in reading file,
//      -1 = end-of-file,
//      otherwise, = number of bytes read.
*/
    if( DEBUG ) {
      printf("BYTES_IO_READ: fptable slot = %d. ", *unit);
      printf("Number of bytes to read = %d\n", *nbytes);
    }

    if( (*iret = fread(buffer, 1, *nbytes, CURRENT_FILE) ) != *nbytes) {
      if( ! feof(CURRENT_FILE) ) {
        *iret = -2;             /*  error in file-handling  */
        perror("bytes_io_read");
        clearerr(CURRENT_FILE);
        return;
      }
      else {
        *iret = -1;             /*  end-of-file */
        clearerr(CURRENT_FILE);
      }
    }

    if( DEBUG ) {
      printf("BYTES_IO_READ: fptable slot = %d. ", *unit);
      printf("Number of bytes read = %d\n", *nbytes);
    }

    return;
}

void c_bytes_io_read__(int* unit,char* buffer,int* nbytes,int* iret) {
  c_bytes_io_read_(unit,buffer,nbytes,iret);
}

void c_bytes_io_read(int* unit,char* buffer,int* nbytes,int* iret) {
  c_bytes_io_read_(unit,buffer,nbytes,iret);
}


/*
//------------------------------------------------------------------------
//  BYTES_IO_WRITE - Write (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_bytes_io_write_(int* unit,char* buffer,int* nbytes,int* iret) {
/*
// Purpose:
//  Writes a block of bytes to a file.
//
//  Function  accepts:
//    unit = the index of a UNIX FILE pointer held in
//           an internal table (fptable).
//
//    nbytes = number of bytes to write.
//
//  Returns:
//    iret:
//      -1 = Could not write to file.
//     >=0 = Number of bytes written.
*/
    if( DEBUG ) {
      printf("BYTES_IO_WRITE: fptable slot = %d. ", *unit);
      printf("Number of bytes to write = %d\n", *nbytes);
    }

    if( (*iret = fwrite(buffer, 1, *nbytes, CURRENT_FILE) ) != *nbytes) {
      perror("bytes_io_write");
      *iret = -1;
    }

    if( DEBUG ) {
      printf("BYTES_IO_WRITE: fptable slot = %d. ", *unit);
      printf("BYTES_IO_WRITE: number of bytes written = %d\n", *iret);
    }

    return;
}

void c_bytes_io_write__(int* unit,char* buffer,int* nbytes,int* iret) {
  c_bytes_io_write_(unit,buffer,nbytes,iret);
}

void c_bytes_io_write(int* unit,char* buffer,int* nbytes,int* iret) {
  c_bytes_io_write_(unit,buffer,nbytes,iret);
}


/*
//------------------------------------------------------------------------
//  BYTES_IO_FSYNC - fsync
//------------------------------------------------------------------------
*/
void c_bytes_io_fsync(int fd, int* iret) {
    // Same implementation of see eckit::fsync
    *iret = fsync(fd);
    while (*iret < 0 && errno == EINTR) {
        *iret = fsync(fd);
    }
    return;
}


/*
//------------------------------------------------------------------------
//  BYTES_IO_FLUSH - flush + fsync
//------------------------------------------------------------------------
*/
void c_bytes_io_flush_(int* unit, int* iret) {
/*
// Purpose:     Flushes file.
*/
    if( DEBUG )
      printf("BYTES_IO_FLUSH: fptable slot = %d\n", *unit);

    // Implementation matching eckit FileHandle

    // flush
    if( (*iret = fflush(CURRENT_FILE)) != 0) {
      perror("bytes_io_flush: fflush failed");
      return;
    }

    // fsync
    c_bytes_io_fsync(fileno(CURRENT_FILE),iret);
    while (*iret < 0 && errno == EINTR) {
      c_bytes_io_fsync(fileno(CURRENT_FILE),iret);
    }
    if (*iret < 0) {
      perror("bytes_io_flush: Cannot fsync");
      return;
    }

    *iret = 0;
    return;
}

void c_bytes_io_flush__(int* unit, int* iret) {
  c_bytes_io_flush_(unit,iret);
}

void c_bytes_io_flush(int* unit, int* iret) {
  c_bytes_io_flush_(unit,iret);
}

/*
//------------------------------------------------------------------------
//   BYTES_IO_CLOSE - close (from FORTRAN)
//------------------------------------------------------------------------
*/
void c_bytes_io_close_(int* unit,int* iret) {
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
      printf("BYTES_IO_CLOSE: fptable slot = %d\n", *unit);
    }

    if( !CURRENT_FILE ) {
      printf("WARNING: bytes_io_close: File (fptable slot = %d) was already closed.\n", *unit);
      return;
    }

    // Flush before closing
    c_bytes_io_flush(unit,iret);
    if( *iret != 0 ) {
      perror("bytes_io_close: could not flush");
      return;
    }

    // Close
    if( ( *iret = fclose(CURRENT_FILE) ) != 0 ) {
      perror("bytes_io_close");
      return;
    }

    CURRENT_FILE = 0;
    return;
}

void c_bytes_io_close__(int* unit,int* iret) {
  c_bytes_io_close_(unit,iret);
}

void c_bytes_io_close(int* unit,int* iret) {
  c_bytes_io_close_(unit,iret);
}
