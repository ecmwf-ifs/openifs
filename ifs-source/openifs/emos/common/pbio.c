/**
* Copyright 1981- ECMWF.
*
* This software is licensed under the terms of the Apache Licence
* Version 2.0 which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation
* nor does it submit to any jurisdiction.
*/

/*
// pbio.c
*/
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>

#ifdef PTHREADS
#include <pthread.h>
#endif

#ifdef FOPEN64
#define OFF_T off64_t
extern FILE *fopen64(const char *, const char *);
#else
#define OFF_T off_t
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

#include "bufrgrib.h"
#include "fort2c.h"
#include "common/fortint.h"
#include "fileRead.h"

#define NAMEBUFFLEN 256
#define MODEBUFFLEN 10

#define CURRENT_FILE (fptable[*unit])

FILE* pbfp(long index) {
  if( (fptable == NULL) || ((int)index < 0) || ((int)index >= fptableSize) )
    return (FILE *) NULL;
  else
    return fptable[index];
}

/*
//------------------------------------------------------------------------
// PBOPEN - Open file (from FORTRAN)
//------------------------------------------------------------------------
*/
#if defined hpR64 || defined hpiaR64
void pbopen_(fortint* unit,_fcd name,_fcd mode,fortint* iret,long l1,long l2) {
#else
void pbopen_(fortint* unit,_fcd name,_fcd mode,fortint* iret,fortint l1,fortint l2) {
#endif
/*
// Purpose:
//  Opens file, returns the index of a UNIX FILE pointer held in
//  an internal table (fptable).
//
// First time through, reads value in environment variable PBIO_BUFSIZE
// (if it is set) and uses it as the size to be used for internal file
// buffers; the value is passed to setvbuf. If PBIO_BUFSIZE is not set,
// a default value is used.
//
// Function  accepts:
//    name = filename
//    mode = r, r+, w
//
//    Note: l1 and l2 are the lengths of the FORTRAN character strings
//          in name and mode.
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

#if (!defined CRAY) && (!defined VAX)
char namebuff[NAMEBUFFLEN+1], modebuff[MODEBUFFLEN+1];
#else
char * namebuff, * modebuff;
#endif

/*
// See if DEBUG switched on.
*/
    if( ! debugSet ) {
      debugLevel = getenv("PBIO_DEBUG");
      if( debugLevel == NULL )
        debugSet = DEBUGOFF;              /* off */
      else {
        int loop;
        for( loop = 0; loop < strlen(debugLevel) ; loop++ ) {
          if( ! isdigit(debugLevel[loop]) ) {
            printf("Invalid number string in PBIO_DEBUG: %s\n", debugLevel);
            printf("PBIO_DEBUG must comprise only digits [0-9].\n");
            debugSet = DEBUGOFF;
          }
        }
        debugSet = DEBUGOFF + atol( debugLevel );
      }
      if( DEBUG ) printf("PBIO_PBOPEN: debug switched on\n");
    }

#if (!defined CRAY) && (!defined VAX)
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
#else
    if(!(namebuff = fcd2char(name))) {
      *iret = -2;
      return;
    }
    if(!(modebuff = fcd2char(mode))) {
      free(namebuff);
      *iret = -2;
      return;
    }
#endif

    strcpy(flags,"");

    /* *unit = (fortint) NULL;  sami bug fix */
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
    if( DEBUG ) printf("PBIO_PBOPEN: filename = %s\n", namebuff);
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
    if( DEBUG ) printf("PBIO_PBOPEN: file open mode = %s\n", flags);

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

    if( DEBUG ) printf("PBIO_PBOPEN: fptable slot = %d\n", *unit);

#ifdef FOPEN64
    if( DEBUG ) printf("PBIO_PBOPEN: using fopen64\n");
    fptable[n] = fopen64(namebuff, flags );
#else
    if( DEBUG ) printf("PBIO_PBOPEN: using fopen\n");
    fptable[n] = fopen(namebuff, flags );
#endif

    if(fptable[n] == NULL) {
      perror(namebuff);
      *iret = -1;
#if (defined CRAY) || (defined VAX)
      free(namebuff);
      free(modebuff);
#endif
#ifdef PTHREADS
      pthread_mutex_unlock(&fpTableBusy);
#endif
      return;
    }

/*
// Now allocate a buffer for the file, if necessary.
*/
    if( ! sizeSet ) {
      envSize = getenv("PBIO_BUFSIZE");
      if( envSize == NULL )
        size = SIZE;             /* default */
      else {
        int loop;
        for( loop = 0; loop < strlen(envSize) ; loop++ ) {
          if( ! isdigit(envSize[loop]) ) {
            printf("Invalid number string in PBIO_BUFSIZE: %s\n", envSize);
            printf("PBIO_BUFSIZE must comprise only digits [0-9].\n");
            exit(1);
          }
        }
        size = atol( envSize );
      }
      if( size <= 0 ) {
        printf("Invalid buffer size in PBIO_BUFSIZE: %s\n", envSize);
        printf("Buffer size defined by PBIO_BUFSIZE must be positive.\n");
        exit(1);
      }
      sizeSet = 1;
    }

    if( DEBUG ) printf("PBIO_PBOPEN: file buffer size = %d\n", size);

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

#if (defined CRAY) || (defined VAX)
    free(namebuff);
    free(modebuff);
#endif

}

#if (defined hpR64) || (defined hpiaR64)
void pbopen(fortint* unit,_fcd name,_fcd mode,fortint* iret,long l1,long l2) {
#else
void pbopen(fortint* unit,_fcd name,_fcd mode,fortint* iret,fortint l1,fortint l2) {
#endif

  pbopen_(unit,name,mode,iret,l1,l2);
}

/*
//------------------------------------------------------------------------
// PBSEEK - Seek (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbseek_(fortint* unit,fortint* offset,fortint* whence,fortint* iret) {
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
      printf("PBIO_PBSEEK: fptable slot = %d\n", *unit);
      printf("PBIO_PBSEEK: Offset = %d\n", my_offset);
      printf("PBIO_PBSEEK: Type of offset = %d\n", my_whence);
    }

    if( my_whence == 2) my_offset = - abs(my_offset);

    *iret = fileTell(CURRENT_FILE);
    if( DEBUG ) printf("PBIO_PBSEEK: current position = %d\n", *iret);
    if( *iret == my_offset && my_whence == 0)
      *iret = 0;
    else
      *iret = fileSeek(CURRENT_FILE, my_offset, my_whence);

    if( DEBUG ) printf("PBIO_PBSEEK: fileSeek return code = %d\n",*iret);

    if( *iret != 0 ) {
      if( ! feof(CURRENT_FILE) ) {
        *iret = -2;             /* error in file-handling */
        perror("pbseek");
      }
      else
        *iret = -1;             /* end-of-file  */

      clearerr(CURRENT_FILE);
      return;
    }

/*
// Return the byte offset from start of file
*/
    *iret = fileTell(CURRENT_FILE);

    if( DEBUG )
      printf("PBIO_PBSEEK: byte offset from start of file = %d\n",*iret);

    return;

}

void pbseek(fortint* unit,fortint* offset,fortint* whence,fortint* iret) {

  pbseek_(unit,offset,whence,iret);
}


void pbseek64_(fortint* unit,long long* offset,fortint* whence,long long* iret) {
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
long long my_offset = (long long) *offset;
int my_whence = (int) *whence;

/*
// Must use negative offset if working from end-of-file
*/
    if( DEBUG ) {
      printf("PBIO_PBSEEK: fptable slot = %d\n", *unit);
      printf("PBIO_PBSEEK: Offset = %d\n", my_offset);
      printf("PBIO_PBSEEK: Type of offset = %d\n", my_whence);
    }

    if( my_whence == 2) my_offset = - abs(my_offset);

    *iret = fileTell(CURRENT_FILE);
    if( DEBUG ) printf("PBIO_PBSEEK: current position = %d\n", *iret);
    if( *iret == my_offset && my_whence == 0)
      *iret = 0;
    else
      *iret = fileSeek(CURRENT_FILE, my_offset, my_whence);

    if( DEBUG ) printf("PBIO_PBSEEK: fileSeek return code = %d\n",*iret);

    if( *iret != 0 ) {
      if( ! feof(CURRENT_FILE) ) {
        *iret = -2;             /* error in file-handling */
        perror("pbseek");
      }
      else
        *iret = -1;             /* end-of-file  */

      clearerr(CURRENT_FILE);
      return;
    }

/*
// Return the byte offset from start of file
*/
    *iret = fileTell(CURRENT_FILE);

    if( DEBUG )
      printf("PBIO_PBSEEK: byte offset from start of file = %d\n",*iret);

    return;

}

void pbseek64(fortint* unit,long long* offset,fortint* whence,long long* iret) {

  pbseek64_(unit,offset,whence,iret);
}

/*
//------------------------------------------------------------------------
// PBTELL - Tells current file position (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbtell_(fortint* unit,fortint* iret) {
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
    *iret = fileTell(CURRENT_FILE);

    if( *iret < 0 ) {
      if( DEBUG ) {           /* error in file-handling */
        printf("PBIO_PBTELL: fptable slot = %d. ", *unit);
        printf("Error status = %d\n", *iret);
      }
      perror("pbtell");
      *iret = -2;
    }

    if( DEBUG ) {
      printf("PBIO_PBTELL: fptable slot = %d. ", *unit);
      printf("Byte offset from start of file = %d\n",*iret);
    }

    return;
}

void pbtell(fortint* unit,fortint* iret) {

  pbtell_(unit,iret);
}

/*
//------------------------------------------------------------------------
// PBTELL64 - Tells current file position (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbtell64_(fortint* unit,long long* iret) {
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
    *iret = fileTell(CURRENT_FILE);

    if( *iret < 0 ) {
      if( DEBUG ) {           /* error in file-handling */
        printf("PBIO_PBTELL64: fptable slot = %d. ", *unit);
        printf("Error status = %d\n", *iret);
      }
      perror("pbtell64");
      *iret = -2;
    }

    if( DEBUG ) {
      printf("PBIO_PBTELL: fptable slot = %d. ", *unit);
      printf("Byte offset from start of file = %d\n",*iret);
    }

    return;
}

void pbtell64(fortint* unit,long long* iret) {

  pbtell64_(unit,iret);
}

/*
//------------------------------------------------------------------------
//  PBREAD - Read (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbread_(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {
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
      printf("PBIO_READ: fptable slot = %d. ", *unit);
      printf("Number of bytes to read = %d\n", *nbytes);
    }

    if( (*iret = fread(buffer, 1, *nbytes, CURRENT_FILE) ) != *nbytes) {
      if( ! feof(CURRENT_FILE) ) {
        *iret = -2;             /*  error in file-handling  */
        perror("pbread");
        clearerr(CURRENT_FILE);
        return;
      }
      else {
        *iret = -1;             /*  end-of-file */
        clearerr(CURRENT_FILE);
      }
    }

    if( DEBUG ) {
      printf("PBIO_READ: fptable slot = %d. ", *unit);
      printf("Number of bytes read = %d\n", *nbytes);
    }

    return;
}

void pbread(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {

  pbread_(unit,buffer,nbytes,iret);
}

/*
//------------------------------------------------------------------------
//  PBREAD2 - Read (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbread2_(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {
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
      printf("PBIO_READ2: fptable slot = %d. ", *unit);
      printf("Number of bytes to read = %d\n", *nbytes);
    }

   if( (*iret = fread(buffer, 1, *nbytes, CURRENT_FILE) ) != *nbytes) {
     if( ! feof(CURRENT_FILE) ) {
       *iret = -2;             /*  error in file-handling  */
       perror("pbread2");
       clearerr(CURRENT_FILE);
     }
   }

  if( DEBUG )
    printf("PBIO_READ2: Number of bytes read = %d\n", *iret);

   return;
}

void pbread2(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {

  pbread2_(unit,buffer,nbytes,iret);
}

/*
//------------------------------------------------------------------------
//  PBWRITE - Write (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbwrite_(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {
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
      printf("PBIO_WRITE: fptable slot = %d. ", *unit);
      printf("Number of bytes to write = %d\n", *nbytes);
    }

    if( (*iret = fwrite(buffer, 1, *nbytes, CURRENT_FILE) ) != *nbytes) {
      perror("pbwrite");
      *iret = -1;
    }

    if( DEBUG ) {
      printf("PBIO_WRITE: fptable slot = %d. ", *unit);
      printf("PBIO_WRITE: number of bytes written = %d\n", *iret);
    }

    return;
}

void pbwrite(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {

  pbwrite_(unit,buffer,nbytes,iret);
}

/*
//------------------------------------------------------------------------
//   PBCLOSE - close (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbclose_(fortint* unit,fortint* iret) {
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
    if( DEBUG )
      printf("PBIO_CLOSE: fptable slot = %d\n", *unit);

    if( ( *iret = fclose(CURRENT_FILE) ) != 0 ) perror("pbclose");
    CURRENT_FILE = 0;

    return;
}

void pbclose(fortint* unit,fortint* iret) {

  pbclose_(unit,iret);
}

/*
//------------------------------------------------------------------------
//  PBFLUSH - flush (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbflush_(fortint * unit) {
/*
// Purpose:     Flushes file.
*/
    if( DEBUG )
      printf("PBIO_FLUSH: fptable slot = %d\n", *unit);

    fflush(CURRENT_FILE);
}

void pbflush(fortint * unit) {

  pbflush_(unit);
}

/*
//------------------------------------------------------------------------
//  PBREAD3 - read (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbread3_(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {
/*
// Purpose:
//  Reads a block of bytes from a file..
//
//  Function returns:
//    status :   -2 = error in reading file,
//               -1 = end-of-file,
//               otherwise, = number of bytes read.
*/
    if( DEBUG )
      printf("PBIO_READ3: number of bytes to read = %d\n", *nbytes);

    *iret = read(*unit, buffer, *nbytes);

    if( DEBUG )
      printf("PBIO_READ3: number of bytes read = %d\n", *iret);
/*
// Error in file-handling
*/
    if(*iret == -1) {
      *iret = -2;
      perror("pbread3");
      return;
    }
/*
// Read problem
*/
    else if(*iret != *nbytes) {
      printf("EOF; pbread3; bytes requested %d; read in: %d\n",
             *nbytes,*iret);
      *iret = -1;
      return;
    }

}

void pbread3(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {

  pbread3_(unit,buffer,nbytes,iret);
}

void pbread4_(fortint* unit,char* buffer,size_t* nbytes,size_t* iret) {
/*
// Purpose:
//  Reads a block of bytes from a file..
//
//  Function returns:
//    status :   -2 = error in reading file,
//               -1 = end-of-file,
//               otherwise, = number of bytes read.
      printf("PBIO_READ4: sizeof = %ld\n", sizeof(buffer));
*/
/*    if( DEBUG ) */
     printf("PBIO_READ4: number of bytes to read = %ld\n", *nbytes);

/* Is 4MB ok ? */
size_t left   = *nbytes;
size_t len    = 0;
size_t chunk = 4*1024*1024; /* Is 4MB ok ? */
size_t pos    = 0;

   while((left > 0) && (len = read(*unit,(buffer+pos),chunk)) > 0)
   {
       left -= len;
       pos += len;
       if(left < chunk)
           chunk = left;
   }
   if(len>0)
       *iret = pos;
   else
       *iret = len;

    if( DEBUG )
         printf("PBIO_READ4: number of bytes READ = %ld\n", pos);

    if(*iret == -1) {
      *iret = -2;
      perror("pbread4");
      return;
    }
/*
// Read problem
*/
    else if(pos != *nbytes) {
      printf("EOF; pbread4; bytes requested %ld; read in: %ld\n",
             *nbytes,pos);
      *iret = -1;
      return;
    }
      printf("PBIO_READ4: number of bytes read = %ld\n", pos);
      printf("PBIO_READ4: return code *iret = %ld\n", *iret);

}

void pbread4(fortint* unit,char* buffer,size_t* nbytes,size_t* iret) {

  pbread4_(unit,buffer,nbytes,iret);
}
void pbread5_(fortint* unit,char* buffer,size_t* nbytes,size_t* iret) {
/*
// Purpose:
//  Reads a block of bytes from a file..
//
//  Function returns:
//    status :   -2 = error in reading file,
//               -1 = end-of-file,
//               otherwise, = number of bytes read.
*/
    if( DEBUG )
      printf("PBIO_READ5: number of bytes to read = %lu\n", *nbytes);

    *iret = read(*unit, buffer, *nbytes);

    if( DEBUG )
      printf("PBIO_READ5: number of bytes read = %lun", *iret);
/*
// Error in file-handling
*/
    if(*iret == -1) {
      *iret = -2;
      perror("pbread5");
      return;
    }
/*
// Read problem
*/
    else if(*iret != *nbytes) {
      printf("EOF; pbread5; bytes requested %lu; read in: %lu\n",
             *nbytes,*iret);
      *iret = -1;
      return;
    }

}

void pbread5(fortint* unit,char* buffer,size_t* nbytes,size_t* iret) {
  pbread5_(unit,buffer,nbytes,iret);
}

/*
//------------------------------------------------------------------------
//  oct_bin3
//------------------------------------------------------------------------
*/
static int oct_bin3(int onum) {
/*
// Converts an integer to octal digits
*/
   char tmp[20];
   int  rc;

   sprintf(tmp,"%d",onum);
   sscanf(tmp,"%o",&rc);
   return rc;
}

/*
//------------------------------------------------------------------------
//  PBOPEN3 - Open file (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbopen3_(fortint* unit,_fcd name,_fcd mode,fortint* iret,fortint l1,fortint l2) {
/*
// Purpose:
//   Opens file, return UNIX FILE pointer.
//
// Function returns:
//   iret:  -1 = Could not open file.
//          -2 = Invalid file name.
//          -3 = Invalid open mode specified
//           0 = OK.
//
// Note: l1 and l2 are the lengths of the character strings in
//       name and mode on SGI.
*/
char *p;
int oflag;
int dmas;
int filemode;
char  flags[4];

#if (!defined CRAY) && (!defined VAX)
char namebuff[NAMEBUFFLEN], modebuff[MODEBUFFLEN];
#else
char * namebuff, * modebuff;
#endif

/*
// See if DEBUG switched on.
*/
    if( ! debugSet ) {
      debugLevel = getenv("PBIO_DEBUG");
      if( debugLevel == NULL )
        debugSet = DEBUGOFF;              /* off */
      else {
        int loop;
        for( loop = 0; loop < strlen(debugLevel) ; loop++ ) {
          if( ! isdigit(debugLevel[loop]) ) {
            printf("Invalid number string in PBIO_DEBUG: %s\n", debugLevel);
            printf("PBIO_DEBUG must comprise only digits [0-9].\n");
            debugSet = DEBUGOFF;
          }
        }
        debugSet = DEBUGOFF + atol( debugLevel );
      }
      if( DEBUG ) printf("PBIO_PBOPEN3: debug switched on\n");
    }

#if (!defined CRAY) && (!defined VAX)
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
#else
    if(!(namebuff = fcd2char(name))) {
      *iret = -2;
      return;
    }
    if(!(modebuff = fcd2char(mode))) {
      free(namebuff);
      *iret = -2;
      return;
    }
#endif

    strcpy(flags,"");

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
    if( DEBUG ) printf("PBIO_PBOPEN: filename = %s\n", namebuff);

/*
// Build open flags from "modes"
*/
    p = modebuff;

    switch(*p) {

      case 'a':
      case 'A': oflag = 0x100 | O_CREAT | 2 | 0x08;
                filemode = 766;
                break;

      case 'c':
      case 'C':
      case 'w':
      case 'W': oflag = 0x100 | O_CREAT | 1;
                filemode = 766;
                break;

      case 'r':
      case 'R': oflag = 0;
                filemode = 444;
                break;

      default:  *iret = -3;
                return;

    }

    if( DEBUG ) printf("PBIO_PBOPEN: file open mode = %s\n", modebuff);

    dmas = umask(000);
    *unit = open(namebuff, oflag, oct_bin3(filemode));
    umask(dmas);

    if(*unit == -1) {
      perror(namebuff);
      perror("pbopen3");
      *iret = -2;
    }

    if( DEBUG ) printf("PBIO_PBOPEN3: file pointer = %0x\n", *unit);

#if (defined CRAY) || (defined VAX)
    free(namebuff);
    free(modebuff);
#endif

}

void pbopen3(fortint* unit,_fcd name,_fcd mode,fortint* iret,fortint l1,fortint l2) {

  pbopen3_(unit,name,mode,iret,l1,l2);
}

/*
//------------------------------------------------------------------------
//  PBCLOSE3 - Close file (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbclose3_(fortint* unit,fortint* iret) {
/*
//
// Purpose:  Closes file.
//
// Function returns:
//   status : non-0 = error in handling file.
//            0 = OK.
*/
    if( DEBUG ) printf("PBIO_PBCLOSE3: file pointer = %0x\n", *unit);

    *iret = close(*unit);

    if(*iret != 0) perror("pbclose3");
}

void pbclose3(fortint* unit,fortint* iret) {

  pbclose3_(unit,iret);
}

/*
//------------------------------------------------------------------------
//  PBSEEK3 - Seek (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbseek3_(fortint* unit,fortint* offset,fortint* whence,fortint* iret) {
/*
//
// Purpose:  Seeks to specified location in file.
//
// Function returns:
//   status: -2 = error in handling file,
//           -1 = end-of-file
//           otherwise,         = byte offset from start of file.
//
//   whence  = 0, from start of file
//           = 1, from current position
//           = 2, from end of file.
*/
fortint my_offset = *offset;
int my_whence;

    if( DEBUG ) {
      printf("PBIO_PBSEEK3: file pointer = %0x\n", *unit);
      printf("PBIO_PBSEEK3: offset = %d\n", my_offset);
      printf("PBIO_PBSEEK3: type of offset = %d\n", *whence);
    }

/*
// Must use negative offset if working from end-of-file
*/
    if( *whence == 2) {
      my_offset = - abs(my_offset);
      my_whence = 2;
    }
    else if(*whence == 0)
    {
      my_whence = 0;
    }
    else
    {
      my_whence = 1;
    }

    if((*iret=lseek(*unit, my_offset, my_whence)) < 0) {
      perror("pbseek3;");
      *iret = -1;           /* end-of-file  */
    }

    if( DEBUG )
      printf("PBIO_PBSEEK3: byte offset from start of file = %d\n",*iret);

}

void pbseek3(fortint* unit,fortint* offset,fortint* whence,fortint* iret) {

  pbseek3_(unit,offset,whence,iret);
}

/*
//------------------------------------------------------------------------
//  PBWRITE3 - Write (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbwrite3_(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {
/*
// Purpose:  Writes a block of bytes to a file.
//
// Function returns:
//   status: -1 = Could not write to file.
//           >=0 = Number of bytes written.
*/
    if( DEBUG ) {
      printf("PBIO_PBWRITE3: file pointer = %0x\n", *unit);
      printf("PBIO_WRITE#: number of bytes to write = %d\n", *nbytes);
    }

    *iret = write(*unit, buffer, *nbytes);
    if( DEBUG )
      printf("PBIO_WRITE3: number of bytes written = %d\n", *iret);

    if( *iret != *nbytes ) {
      perror("pbwrite3: ");
      *iret = -1;
    }

}

void pbwrite3(fortint* unit,char* buffer,fortint* nbytes,fortint* iret) {

  pbwrite3_(unit,buffer,nbytes,iret);
}

/*
//------------------------------------------------------------------------
//  GRIBREAD
//------------------------------------------------------------------------
*/
void gribread_(
char* buffer,fortint* buffsize,fortint* readsize,fortint* status,fortint* unit){
/*
//  Called as a FORTRAN subroutine:
//
//    CALL GRIBREAD( KARRAY, KINLEN, KOUTLEN, IRET, KUNIT )
//
*/
fortint holdsize = *buffsize;

/*
// Read GRIB product
*/

    *status = readprod("GRIB",buffer,&holdsize,fileRead,fileSeek,fileTell,
                       CURRENT_FILE);
    *readsize = abs(holdsize );

    if( DEBUG ) {
      printf("PBIO_GRIBREAD: fptable slot = %d. ", *unit);
      printf("Number of bytes read = %d\n", *readsize);
    }

    return;
}

void gribread(
char* buffer,fortint* buffsize,fortint* readsize,fortint* status,fortint* unit){

  gribread_(buffer,buffsize,readsize,status,unit);
}

/*
//------------------------------------------------------------------------
//  BUFRREAD
//------------------------------------------------------------------------
*/
void bufrread_(
char* buffer,fortint* buffsize,fortint* readsize,fortint* status,fortint* unit){
/*
//  Called as a FORTRAN subroutine:
//
//    CALL BUFRREAD( KARRAY, KINLEN, KOUTLEN, IRET, KUNIT )
//
*/
fortint holdsize = *buffsize;

/*
// Read BUFR product
*/
    *status = readprod("BUFR",buffer,&holdsize,fileRead,fileSeek,fileTell,
                       CURRENT_FILE);
    *readsize =  abs(holdsize );

    if( DEBUG ) {
      printf("PBIO_BUFRREAD: fptable slot = %d. ", *unit);
      printf("Number of bytes read = %d\n", *readsize);
    }

    return;
}

void bufrread(
char* buffer,fortint* buffsize,fortint* readsize,fortint* status,fortint* unit){

  bufrread_(buffer,buffsize,readsize,status,unit);
}

/*
//------------------------------------------------------------------------
//  PSEUREAD
//------------------------------------------------------------------------
*/
void pseuread_(
char* buffer,fortint* buffsize,fortint* readsize,fortint* status,fortint* unit){
/*
//  Called as a FORTRAN subroutine:
//
//    CALL PSEUREAD( KARRAY, KINLEN, KOUTLEN, IRET, KUNIT )
//
*/
fortint holdsize = *buffsize;

/*
// Read GRIB product
*/
    *status = readprod(NULL,buffer,&holdsize,fileRead,fileSeek,fileTell,
                       CURRENT_FILE);
    *readsize = abs(holdsize );

    if( DEBUG ) {
      printf("PBIO_PSEUREAD: fptable slot = %d. ", *unit);
      printf("Number of bytes read = %d\n", *readsize);
    }

    return;
}

void pseuread(
char* buffer,fortint* buffsize,fortint* readsize,fortint* status,fortint* unit){

  pseuread_(buffer,buffsize,readsize,status,unit);
}

/*
//------------------------------------------------------------------------
//  PBSIZE
//------------------------------------------------------------------------
*/
void pbsize_(fortint* unit,fortint* plen) {
/*
//  Returns the size in bytes of the next GRIB, BUFR, TIDE, BUDG, DIAG
//  product.
//
//  Called from FORTRAN:
//      CALL PBSIZE( KUNIT, LENGTH)
//
//  unit  = file id returned from PBOPEN.
//  plen  = size in bytes of the next product.
//        = -2 if error allocating memory for internal buffer.
//
//  The input file is left positioned where it started.
*/
fortint iret;
char statbuff[BUFFLEN];
char * buff;
long offset, loop = 1;

/*
//  Use a smallish buffer for processing; this should suffice for all cases
//  except versions -1 and 0 of GRIB and large BUFR products
*/
    offset = (fortint) fileTell( CURRENT_FILE);
    if( DEBUG ) {
      printf("PBIO_SIZE: fptable slot = %d. ", *unit);
      printf("Current file position = %lu\n", offset);
    }

    *plen = BUFFLEN;
    if( DEBUG )
      printf("PBIO_SIZE: current buffer size = %d\n", *plen);

    iret = readprod(NULL,statbuff,plen,fileRead,fileSeek,fileTell,CURRENT_FILE);
    if( iret == -2 ) {
      printf("readprod error %d\n", iret);
      *plen = -2;
      return;
    }
/*
//  If the smallish buffer is too small, progressively increase it until
//  big enough
*/
    while ( iret == -4 ) {
      loop++;
      buff = (char *) malloc( BUFFLEN*loop);
      if( buff == NULL) {
        perror("malloc failed in PBSIZE");
        *plen = -2;
        return;
      }
      *plen = BUFFLEN*loop;
      if( DEBUG )
        printf("PBIO_SIZE: buffer size increased to: %d\n", *plen);

      offset = (fortint) fileSeek( CURRENT_FILE, offset, SEEK_SET);
      offset = (fortint) fileTell( CURRENT_FILE);
      iret = readprod(NULL,buff,plen,fileRead,fileSeek,fileTell,CURRENT_FILE);
      free(buff);
    }

    if( iret == -2 ) {
      printf("readprod error %d\n", iret);
      *plen = -2;
    }
/*
//  Put the file pointer back where it started
*/
    if( DEBUG ) {
      printf("PBIO_SIZE: file pointer set back to: %lu\n", offset);
      printf("PBIO_SIZE: Product size = %d\n", *plen);
    }
    offset = (fortint) fileSeek( CURRENT_FILE, offset, SEEK_SET);

    return ;
}

void pbsize(fortint* unit,fortint* plen) {

  pbsize_(unit,plen);
}

/*
//------------------------------------------------------------------------
//  CREXRD
//------------------------------------------------------------------------
*/
#define END_OF_FILE -1
#define FILE_READ_ERROR -2
#define USER_BUFFER_TOO_SMALL -3
#define FILE_TOO_SMALL -5
#define MINIMUM_CREX_SIZE 13
#define CREX 0x43524558

typedef char * String;

void crexrd_(String buffer,int* bufflen,int* size,int* status,fortint * unit) {
/*
//  Called from FORTRAN:
//    CALL CREXRD( KARRAY, KINLEN, NREAD, IRET, KUNIT )
*/
int loop;
OFF_T foundPosition;
int number, crexFound = 0, endFound = 0;
String endBuffer;
String next;
char plplcrcrlf7777[10] = {0,0,0,0,0,0,0,0,0,0};
char PlPlCrCrLf7777[10] = {0x2b,0x2b,0x0d,0x0d,0x0a,0x37,0x37,0x37,0x37,0x00};

/*
// Check buffer big enough for CREX search
*/
  if( *bufflen < MINIMUM_CREX_SIZE ) {
    *status = USER_BUFFER_TOO_SMALL;
    return;
  }

/*
// Look for CREX
*/
  for( loop = 0; loop <= 4; loop++) buffer[loop] = '\0';

  while( !crexFound ) {
    buffer[0] = buffer[1];
    buffer[1] = buffer[2];
    buffer[2] = buffer[3];
    number = fread((buffer+3), 1, 1, CURRENT_FILE);
    if( feof(CURRENT_FILE) ) {
      *status = END_OF_FILE;
      return;
    }
    if( (number != 1) || ferror(CURRENT_FILE) ) {
      perror("crexrd file read error");
      *status = FILE_READ_ERROR;
      return;
    }
    if( strcmp(buffer,"CREX") == 0 ) {
      crexFound = 1;
#ifdef FOPEN64
      foundPosition = ftello64(CURRENT_FILE) - 4;
#else
      foundPosition = ftell(CURRENT_FILE) - 4;
#endif
    }
  }

/*
// Read some more characters into the buffer
*/
  number = fread((buffer+4), 1, ((*bufflen)-4), CURRENT_FILE);
  if( ferror(CURRENT_FILE) ) {
    perror("crexrd file read error");
    *status = FILE_READ_ERROR;
    return;
  }
  endBuffer = buffer + number + 3;

/*
// Look for ++CrCrLf7777 at end of product
*/
  next = buffer+4;
  endFound = 0;
  for( loop = 0; loop < 8; loop++ )
    plplcrcrlf7777[loop] = *(next++);
  plplcrcrlf7777[9] = '\0';

  while( (!endFound) && (next<=endBuffer) ) {
    plplcrcrlf7777[8] = *(next++);
    if( strcmp(plplcrcrlf7777,PlPlCrCrLf7777) == 0 ) {
      endFound = 1;
      *size = (int) (next - buffer);
    }

    for( loop = 0; loop < 8; loop++ )
      plplcrcrlf7777[loop] = plplcrcrlf7777[loop+1];
  }

  if( !endFound ) {
    if( feof(CURRENT_FILE) )
      *status = END_OF_FILE;
    else
      *status = USER_BUFFER_TOO_SMALL;
    return;
  }

/*
// Position file at end of CREX product
*/
#ifdef FOPEN64
  *status = fseeko64(CURRENT_FILE, (foundPosition+(*size)), 0);
#else
  *status = fseek(CURRENT_FILE, (foundPosition+(*size)), 0);
#endif

}

void crexrd(String buffer,int* bufflen,int* size,int* status,fortint * unit) {

  crexrd_(buffer,bufflen,size,status,unit);
}
