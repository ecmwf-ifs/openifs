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
//    readprod.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#ifdef FOPEN64
#define OFF_T off64_t
#else
#define OFF_T off_t
#endif

#ifdef sun4
#include <memory.h>
#endif

#include "bufrgrib.h"
#include "common/fortint.h"
#include "fileRead.h"
#include "sizeRoutines.h"

#define THREE_BYTE_LONG(p) ((((*(p))<<16) & 0xff0000) | (((*(p+1))<<8) & 0xff00) | (*(p+2)) & 0xff)
#define CHECK(stat,message,code) if((stat)) {perror(message);return(code);}

#define WMOBIT1 0x80
#define WMOBIT2 0x40
#define WMOBIT7 0x02
#define WMOBIT8 0x01

#define END_OF_FILE      -1
#define INTERNAL_ERROR   -2
#define BUFFER_TOO_SMALL -3
#define USER_BUFFER_TINY -4
#define MISPLACED_7777   -5

#define SECT_0_LEN 4
#define LEN_7777 4
#define SMALL 40
#define LARGEBUF 200000

#define GRIB 0x47524942
#define BUFR 0x42554652
#define BUDG 0x42554447
#define TIDE 0x54494445
#define DIAG 0x44494147
#define CREX 0x43524558
#define CODE_7777 0x37373737

static int grab(char * , char * , long ,long ,long * );
static fortint waveLength(char * );

#ifdef FOPEN64
fortint readprod( char * prod_id, char * buffer, fortint * size,
              fortint (*fileRead)(char *, fortint , void *),
              OFF_T (*fileSeek)(void *, OFF_T, fortint),
              OFF_T (*fileTell)(void *),
              void * stream)
#else
fortint readprod( char * prod_id, char * buffer, fortint * size,
              fortint (*fileRead)(char *, fortint , void *),
              fortint (*fileSeek)(void *, fortint, fortint),
              fortint (*fileTell)(void *),
              void * stream)
#endif
/*
//  Reads a BUFR, GRIB, BUDG, TIDE, DIAG product.
//
//  prod_id = "BUFR", "GRIB", "BUDG", "TIDE", "DIAG" ...
//  buffer = buffer to hold product,
//  size = on input, size of buffer in bytes,
//         on output, number of bytes read.
//  fileRead = function used to read the product,
//  fileSeek = function used to reposition within the product byte stream,
//  fileTell = function used to report the current byte position within
//             the product byte stream,
//  stream = data describing the input stream (eg. FILE *).
//
//  Returns
//    -1 if can't find product  (eg. end-of-file)
//    -2 if internal processing error (malloc fail, file seek error)
//    -3 if buffer too small for whole product.
//    -4 if user buffer too small to even start product processing.
//    -5 if the 7777 group is in the wrong place
//
*/
{
static char * shold = NULL;           /* buffer used to read product */
char * hold = NULL;
fortint holdsize = SMALL, numread;
fortint found = 0, length, prodlen;
fortint given_buflen = *size;
char p;
unsigned long code=0, wanted = 0;
OFF_T status = 0, i, present_position = 0;
fortint one = 1;

/* See if user gave a buffer for holding the product
*/
  if ( buffer == NULL ) {         /* No buffer given, get some */
    if ( shold == NULL ) {        /* but only first time round */
      shold = (char *) malloc( LARGEBUF);
      CHECK((shold == NULL),"malloc failed in readnext",INTERNAL_ERROR);
    }
    given_buflen = LARGEBUF;
    hold = shold;
  }
  else                             /* User buffer given */
    hold = buffer;

/* Make sure the user gave some buffer
*/
  if ( given_buflen < holdsize ) return USER_BUFFER_TINY;

/* Look for product identifier
*/
  if ( prod_id != NULL )
    for ( i = 0; i < 4; i++ )
      wanted = ( (wanted << 8) + *(prod_id+i) ) & 0xFFFFFFFF;

/* Read some bytes from the product
*/
  do {
    numread = fileRead( &p, one, stream );
    if ( numread <= 0 ) {
      *size = 0;
      return END_OF_FILE;
    }

    code = ( (code << 8) + p ) & 0xFFFFFFFF;

    if ( prod_id == NULL ) {
      switch(code) {
        case BUDG:
        case BUFR:
        case GRIB:
        case TIDE:
        case DIAG:
        case CREX: found = 1;
      }
    }
    else
      if ( code == wanted ) found = 1;

  } while ( ! found );

/* Find the product length
*/
  prodlen = prodsize( code, hold, given_buflen, &holdsize, fileRead, stream);

  if( prodlen == -4 ) return USER_BUFFER_TINY;

  if( prodlen < -4 ) {           /* special case for very large wave products */
    fortint estimate;

    estimate = -prodlen;
    prodlen = estimate - 128;
    *size = prodlen;

    present_position = fileTell((FILE *)stream);
    status = fileSeek( (FILE *)stream, prodlen - holdsize , SEEK_CUR);
    CHECK(status,"fileSeek error in readprod",INTERNAL_ERROR);

/*  Read bytes from the product upto and including 7777
*/
    code = 0;
    do {
      numread = fileRead( &p, one, stream );
      if ( numread <= 0 ) return MISPLACED_7777;

      *size = *size + 1;
      if( *size > estimate ) return MISPLACED_7777;

      code = ( (code << 8) + p ) & 0xFFFFFFFF;

    } while ( code != CODE_7777 );

    if( given_buflen < *size )
      return BUFFER_TOO_SMALL;
    else {
      status = fileSeek( (FILE *)stream, present_position, SEEK_SET);
      CHECK(status,"fileSeek error in readprod",INTERNAL_ERROR);

      numread = fileRead(hold+present_position,(*size-present_position),stream);
      if ( numread <= 0 ) return INTERNAL_ERROR;
      return 0;
    }
  }

/* Move to end of product and check position of 7777 group
   in the case that the user gave a NULL buffer and zero length
*/
  if ( *size == 0 ) {
    char grp_7777[5];

    *size = prodlen;             /* report the actual product length */

    status = fileSeek( (FILE *)stream, prodlen - holdsize - LEN_7777, SEEK_CUR);
    CHECK(status,"fileSeek error in readprod",INTERNAL_ERROR);

    numread = fileRead( grp_7777, LEN_7777 ,stream );
    grp_7777[4] = '\0';
    if( strcmp(grp_7777,"7777") != 0 ) return MISPLACED_7777;

    return BUFFER_TOO_SMALL;
  }

/* Otherwise read as much of the product as possible into the buffer, then
   move to end of product and check position of 7777 group
*/
  *size = prodlen;                /* report the actual product length */

  length = (given_buflen > prodlen) ? prodlen : given_buflen ;

  if ( length > holdsize ) {

/*
//  Before reading rest of product, put zero at end of expected 7777
//  group to ensure that it is not left over from an earlier product
//  of the same length
*/
    *(hold+length-1) = '\0';
    numread = fileRead( hold+holdsize, length-holdsize ,stream );
    if( numread <= 0 ) {
      *size = holdsize-numread;
      return END_OF_FILE;
    }

    if ( given_buflen < prodlen ) {
      char grp_7777[5];

      status = fileSeek((FILE*)stream,prodlen-given_buflen - LEN_7777,SEEK_CUR);
      CHECK(status,"fileSeek error in readprod",INTERNAL_ERROR);

      numread = fileRead( grp_7777, LEN_7777 ,stream );
      if( numread <= 0 ) return END_OF_FILE;

      grp_7777[4] = '\0';
      if( strcmp(grp_7777,"7777") != 0 ) return MISPLACED_7777;

      return BUFFER_TOO_SMALL;
    }
  }

  if(strncmp(hold+prodlen-LEN_7777,"7777",(size_t)4)!=0) return MISPLACED_7777;

  return *size;

}

fortint readprod_decode_unsigned_byte_long(const unsigned char* p, long o, int l)
{
    fortint accum = 0;
    int i = 0;
    unsigned char b = p[o++];
    accum  <<= 8;
    accum |= b ;

    for ( i=1; i<l; i++ )
    {
        b = p[o++];
        accum <<= 8;
        accum |= b ;
    }
    return accum;
}

static fortint gribsize(char * hold, fortint leng, fortint * holdsize,
                     fortint (*fileRead)(), void * stream)
/*
//  Calculates the size in bytes of a GRIB product.
//
//  hold = buffer to hold product,
//  leng = length of buffer,
//  holdsize = number of bytes of the product in the hold buffer.
//             Note that this increases if necessary as more bytes are read.
//  fileRead = function used to read the product,
//  stream = data describing the input stream (eg. FILE *).
*/
{
fortint length, numread, num, hsize = *holdsize;
fortint section2, section3 ;
unsigned char * phold = (unsigned char *) hold;
const char* big_prod_flag = NULL;

/* Need more bytes to decide on which version of GRIB
*/
  if ( leng < 24 ) return USER_BUFFER_TINY;

/* Put first 24 bytes in buffer
*/
  num = 24;
  if ( hsize < num) {
    numread = fileRead( hold + hsize, num - hsize, stream);
    if ( numread <= 0 ) {          /* eg. on END_OF_FILE */
      *holdsize -= numread;
      return *holdsize;
    }
    hsize = num;
  }
  *holdsize = hsize;

/* See if the GRIB version is 0 ...
*/
  if ( THREE_BYTE_LONG(phold+4) == 24 ) {
    length = 28;                   /* GRIB + section 1 */

/* Check for presence of sections 2 and 3
*/
    section2 = ( hold[11] & WMOBIT1 ) ;
    section3 = ( hold[11] & WMOBIT2 ) ;

/* Add up all lengths of sections
*/
    return lentotal(hold, holdsize, leng, length, section2, section3,
                    fileRead, stream);
  }


/* ... or version 2 ...
*/
if ( hold[7] == 2 ){
     return readprod_decode_unsigned_byte_long((const unsigned char*)hold,8,8);
    }
/* ... or version 1 ...
*/

  if ( ( hold[21] != 0 ) || ( hold[22] != 0 ) ) {

/* Nightmare fixup for very large GRIB products (eg 2D wave spectra).

   If the most-significant of the 24 bits is set, this indicates a
   very large product; the size has to be rescaled by a factor of 120.
   This is a fixup to get round the GRIB code practice of representing
   a product length by 24 bits. It is only possible because the
   (default) rounding for GRIB products is 120 bytes.

*/
    fortint fixlen;

    fixlen = THREE_BYTE_LONG(phold+4);
    if( fixlen <= 0x7fffff )
      return fixlen;
    else
    {
      fortint largeLen;
      unsigned long code = 0;

      big_prod_flag = getenv("ECMWF_CODED_BIG_PRODUCT");
      if(big_prod_flag)
         if(strcmp(big_prod_flag,"OFF") == 0) return fixlen;

      fixlen = (fixlen & 0x7fffff)*120;  /* first guess at length */

      num = fixlen;

      if( leng < num ) { /* -unsure - apparently test is given buffer is too small */
        num = fixlen - 128;
        numread = fileRead( hold + hsize, leng - hsize, stream);
        if ( numread <= 0 ) {          /* eg. on END_OF_FILE */
          *holdsize -= numread;
          return *holdsize;
        }
        *holdsize = leng;
        return (-fixlen);
      }
      else {
        fortint fileReadResult =0;
        if ( hsize < num) { /* buffer is big enough */

          /* go to the estimate minus 128 bytes to look after the 7777*/
          numread = fileRead( hold + hsize, num - hsize - 128, stream);

          if ( numread <= 0 ) {          /* eg. on END_OF_FILE */
            *holdsize -= numread;
            return *holdsize;
          }

          /* search byte after byte for the 7777 in the stream and stop when found*/
          do {
            if (num<numread) return MISPLACED_7777;
            fileReadResult = fileRead( hold + numread + hsize, 1, stream);

            if ( fileReadResult <= 0 ) {          /* eg. on END_OF_FILE ---should not happen*/
                 *holdsize += numread;
                 *holdsize -= fileReadResult;
                 return *holdsize; /*returns the number of bytes actually read */
            }
            numread += fileReadResult;
            code = ( (code << 8) + hold[hsize + numread-1] ) & 0xFFFFFFFF;
          } while ( code != CODE_7777 );
          hsize = numread;
        }
        *holdsize = hsize;
        largeLen = waveLength(hold);
      }

      num = largeLen;
      *holdsize = num;
      return num;
    }

  }

  length = 24;                     /* GRIB + section 1 */

/* Check for presence of sections 2 and 3
*/
  section2 = ( hold[7] & WMOBIT8 ) ;
  section3 = ( hold[7] & WMOBIT7 ) ;

/* Add up all lengths of sections
*/
  return lentotal(hold, holdsize, leng, length, section2, section3,
                    fileRead, stream);

}


static fortint lentotal(char *hold, fortint *holdsize, fortint leng, fortint length,
                     fortint section2, fortint section3,
                     fortint (*fileRead)(), void *stream)
/*
//  Returns the total length in bytes of all sections of the GRIB product.
//
//  hold = buffer to hold product,
//  leng = length of buffer,
//  holdsize = number of bytes of the product in the hold buffer.
//             Note that this increases if necessary as more bytes
//             are read.
//  length = length of (section 0 + section 1).
//  section2 is TRUE if section 2 is present in the product.
//  section3 is TRUE if section 3 is present in the product.
//  fileRead = function used to read the product,
//  stream = data describing the input stream (eg. FILE *).
*/
{
fortint numread, hsize = *holdsize;
unsigned char * phold = (unsigned char *) hold;
fortint next, next_sec = 4;

/* Adjust count of sections to check
*/
  if ( section2 ) next_sec--;
  if ( section3 ) next_sec--;

/* Get the size of the next section
*/
  if ( leng < length ) return USER_BUFFER_TINY;
  if ( hsize < length) {
    numread = fileRead( hold + hsize, length - hsize, stream);
    if ( numread <= 0 ) {        /* eg. on END_OF_FILE */
      *holdsize -= numread;
      return *holdsize;
    }
    hsize = length;
  }
  *holdsize = hsize;

/* Get the size of remaining sections
*/
  for ( next = next_sec; next < 5 ; next++ ) {
    if ( leng < (length+4) ) return USER_BUFFER_TINY;
    if ( hsize < (length+4)) {
      numread = fileRead( hold+hsize, (length+4)-hsize, stream);
      if ( numread <= 0 ) {    /* eg. on END_OF_FILE */
        *holdsize -= numread;
        return *holdsize;
      }
      hsize = length + 4;
    }
    *holdsize = hsize;
    length += THREE_BYTE_LONG(phold+length);
  }

/* Add on the size of section 5
*/
  length += LEN_7777;

  return length;
}


static fortint bufrsize(char * hold, fortint leng, fortint * holdsize,
                     fortint (*fileRead)(), void * stream)
/*
//  Returns the size in bytes of the BUFR code product.
//
//  hold = buffer to hold product,
//  leng = length of buffer,
//  holdsize = number of bytes of the product in the hold buffer.
//         Note that this increases if necessary as more bytes are read.
//  fileRead = function used to read the product,
//  stream = data describing the input stream (eg. FILE *).
*/
{
unsigned char * phold = (unsigned char *) hold;
fortint numread, hsize = *holdsize;
fortint num, length;
fortint next, next_sec = 3;


/* Need more bytes to decide on which version of BUFR
*/
  if ( leng < 24 ) return USER_BUFFER_TINY;
  num = 24;                        /* put first 24 bytes in buffer*/
  if ( hsize < num) {
    numread = fileRead( hold + hsize, num - hsize, stream);
    if ( numread <= 0 ) {        /* eg. on END_OF_FILE */
      *holdsize -= numread;
      return *holdsize;
    }
    hsize = num;
  }
  *holdsize = hsize;

/* If it's Edition 2, or later, octets 5-7 give full product size
*/
  if ( hold[7] > 1 ) return THREE_BYTE_LONG(phold+4);

/* Otherwise, we have to step through the individual sections
   adding up the lengths
*/

/* Add on the length of section 1 and ensure enough of product is in
   memory to continue
*/
  length = SECT_0_LEN + THREE_BYTE_LONG(phold+4);
  if ( leng < (length+4) ) return USER_BUFFER_TINY;
  if ( hsize < (length+4)) {
    numread = fileRead( hold+hsize, (length+4)-hsize, stream);
    if ( numread <= 0 ) {        /* eg. on END_OF_FILE */
      *holdsize -= numread;
      return *holdsize;
    }
    hsize = length + 4;
    *holdsize = hsize;
  }

/* Check for presence of section 2
*/
  if ( hold[11] & WMOBIT1 ) next_sec = 2;

/* Get the size of remaining sections
*/
  for ( next = next_sec; next < 5 ; next++ ) {
    length += THREE_BYTE_LONG(phold+length);
    if ( leng < (length+4) ) return USER_BUFFER_TINY;
    if ( hsize < (length+4)) {
      numread = fileRead( hold+hsize,(length+4)-hsize, stream);
      if ( numread <= 0 ) {    /* eg. on END_OF_FILE */
        *holdsize -= numread;
        return *holdsize;
      }
      hsize = length + 4;
      *holdsize = hsize;
    }
  }

/* Add on the size of section 5
*/
  length += LEN_7777;
  if ( leng < length ) return USER_BUFFER_TINY;

  return length;

}

static fortint tide_budg_size(char * hold, fortint leng, fortint * holdsize,
                           fortint (*fileRead)(), void * stream)
/*
//  Returns the size in bytes of the TIDE/BUDG/DIAG code product.
//
//  hold = buffer to hold product,
//  leng = length of buffer,
//  holdsize = number of bytes of the product in the hold buffer.
//             Note that this increases if necessary as more bytes are read.
//  fileRead = function used to read the product,
//  stream = data describing the input stream (eg. FILE *).
*/
{
unsigned char * phold = (unsigned char *) hold;
fortint numread, hsize = *holdsize;
fortint num, length;

/* Need more bytes to get length of section 1
*/
  num = 8;                        /* put first 8 bytes in buffer */
  if ( leng < num ) return USER_BUFFER_TINY;
  if ( hsize < num) {
    numread = fileRead( hold + hsize, num - hsize, stream);
    if ( numread <= 0 ) {       /* eg. on END_OF_FILE */
      *holdsize -= numread;
      return *holdsize;
    }
    hsize = num;
  }
  *holdsize = hsize;

/* Have to step through individual sections adding up the lengths
*/

/* Add on the length of section 1 and ensure enough of product is in
   memory to continue
*/
  length = SECT_0_LEN + THREE_BYTE_LONG(phold+4);
  if ( leng < (length+4) ) return USER_BUFFER_TINY;
  if ( hsize < (length+4)) {
    numread = fileRead( hold+hsize, (length+4)-hsize, stream);
    if ( numread <= 0 ) {        /* eg. on END_OF_FILE */
      *holdsize -= numread;
      return *holdsize;
    }
    hsize = length + 4;
    *holdsize = hsize;
  }

/* Get the size of remaining section
*/
  length += THREE_BYTE_LONG(phold+length);

/* Add on the size of section 5
*/
  length += LEN_7777;

  if ( leng < length ) return USER_BUFFER_TINY;
  return length;
}

static fortint prodsize(fortint code, char * hold, fortint leng, fortint * holdsize,
                     fortint (*fileRead)(), void * stream)
/*
//  Returns size of BUFR, GRIB, BUDG, TIDE, DIAG product in bytes.
//
//  hold = buffer holding product,
//  leng = size of buffer in bytes,
//  holdsize = number of bytes of product already read into hold.
//  fileRead = function used to read the product,
//  stream = data describing the input stream (eg. FILE *).
*/
{
long lcode = code;

  *holdsize = 4;
  switch( lcode )
  {
     case BUFR: memcpy(hold,"BUFR",4);
                return bufrsize( hold, leng, holdsize, fileRead, stream);

     case BUDG: memcpy(hold,"BUDG",4);
                return tide_budg_size( hold, leng, holdsize, fileRead, stream);

     case GRIB: memcpy(hold,"GRIB",4);
                return gribsize( hold, leng, holdsize, fileRead, stream);

     case TIDE: memcpy(hold,"TIDE",4);
                return tide_budg_size( hold, leng, holdsize, fileRead, stream);

     case DIAG: memcpy(hold,"DIAG",4);
                return tide_budg_size( hold, leng, holdsize, fileRead, stream);

     case CREX: memcpy(hold,"CREX",4);
                return crex_size( stream);

     default: return 0;
    }
}


#define BIT1 0x80
#define BIT2 0x40

static fortint waveLength(char * buffer)
{
/*
//  On entry, buffer contains a GRIB edition 1 product (somewhere).
//
//  The value returned is the length of the GRIB product calculated
//  from the individual lengths of sections 0, 1, 2, 3, 4 and 5;
//  sections 2 and 3 are optional and may or may not be present.
//
//  If there is a problem processing the product, or if GRIB edition
//  -1 or 0 is encountered, the return value is -1.
*/
int found = 0;
int code = 0;
long bytes_read = 0, advance;
char p, edit_num, flag23;
char size[3];
int section0 = 8, section1, section2, section3, section4, section5 = 4;
long total = 0;
size_t rtotal;
long s0;

/*  Read bytes until "GRIB" found */

    do
    {
        if( grab(buffer, &p, 1, 1, &bytes_read) != 0) return (-1);
        code = ( (code << 8) + p ) & 0xFFFFFFFF;
        if (code == GRIB ) found = 1;
    } while ( ! found );
    s0 = bytes_read - 4;
    bytes_read = 4;

    if( grab(buffer, size, 3, 1, &bytes_read) != 0) return (-1);
    rtotal = THREE_BYTE_LONG(size);

    if( rtotal > 0x800000 ) {
        total = (rtotal&0x7fffff) * 120;
    }

/*  Check the edition number */

    if( grab(buffer, &edit_num, 1, 1, &bytes_read) != 0) return (-1);
    if( edit_num != 1 ) {
      printf("Cannot handle GRIB edition 0\n");
      return (-1);
    }                             /* reject edition 0 */

    if( (*(buffer+21-s0) == '\0') && (*(buffer+22-s0) == '\0') ) {
      printf("Cannot handle GRIB edition -1\n");
      return (-1);
    }                             /* reject edition -1 */

/*  Read length of section 1 */
    if( grab(buffer, size, 3, 1, &bytes_read) != 0) return (-1);
    section1 = THREE_BYTE_LONG(size);

/*  Now figure out if sections 2/3 are present */

    advance = 4;
    bytes_read += advance;
    if( grab(buffer, &flag23, 1, 1, &bytes_read) != 0) return (-1);
    section2 = flag23 & BIT1;
    section3 = flag23 & BIT2;

/*  Advance to end of section 1 */

    advance = section1 - (bytes_read - section0);
    bytes_read += advance;

/*  Read section 2 length if it is given*/

    if( section2 )
    {
        if( grab(buffer, size, 3, 1, &bytes_read) != 0) return (-1);
        section2 = THREE_BYTE_LONG(size);
        advance = section2 - (bytes_read - section0 - section1);
        bytes_read += advance;
    }

/*  Read section 3 length if it is given*/

    if( section3 )
    {
        if( grab(buffer, size, 3, 1, &bytes_read) != 0) return (-1);
        section3 = THREE_BYTE_LONG(size);
        advance = section3 - (bytes_read - section0 - section1 - section2);
        bytes_read += advance;
    }

/*  Read section 4 length */

    if( grab(buffer, size, 3, 1, &bytes_read) != 0) return (-1);
    section4 = THREE_BYTE_LONG(size);
    /*  If it happens that the section 2 is not the difference of bytes
        between the real message length and the coded message length, that mean that
        the 120 trick is not used, if so return the real length   */
     if(section4 > 120 )
       return rtotal;
    section4 = total + 3 - bytes_read - section4;

    return (section0+section1+section2+section3+section4+section5);
}

static int grab(char * buffer, char * where, long size,long cnt,long * num_bytes_read)
{
long number = size*cnt;

    memcpy(where, (buffer+(*num_bytes_read)), number);
    *num_bytes_read += number;

    return 0;
}


#define BUFFLEN 1000

static fortint crex_size( void * stream) {
/*
//  Returns the size in bytes of the TIDE/BUDG/DIAG code product.
//
//  stream = data describing the input stream (eg. FILE *).
*/
char buffer[BUFFLEN];
int startPosition, inBuffer;
char * next, * endBuffer;
char plplcrcrlf7777[10] = {0,0,0,0,0,0,0,0,0,0};
char PlPlCrCrLf7777[10] = {0x2b,0x2b,0x0d,0x0d,0x0a,0x37,0x37,0x37,0x37,0x00};
int endFound, loop, status, size = 0;
int size7777 = sizeof(PlPlCrCrLf7777) - 1;

/*
// Record current file position
*/
  startPosition = fileTell((FILE *) stream);
  if( startPosition < 0 ) {
    perror("crex_size: error recording current file position.");
    exit(1);
  }
/*
// Read some bytes to start with
*/
  inBuffer = fileRead(buffer, BUFFLEN, (FILE *) stream);
  if( ferror((FILE *) stream) ) {
    perror("crex_size: file read error");
    exit(1);
  }
  next = buffer;
  endBuffer = next + abs(inBuffer);
/*
// Look for ++CrCrLf7777 at end of product
*/
  endFound = 0;
  while( !endFound ) {
    for( loop = 0; loop < 8; loop++ )
      plplcrcrlf7777[loop] = *(next++);
      plplcrcrlf7777[9] = '\0';

    while( next<=endBuffer ) {
      plplcrcrlf7777[8] = *(next++);
      if( strcmp(plplcrcrlf7777,PlPlCrCrLf7777) == 0 ) {
        endFound = 1;
/*
//      Position file where it started
*/
        status = fileSeek(stream, startPosition, SEEK_SET);
        if( status != 0 ) {
          perror("crex_size: file repositioning error");
          exit(1);
        }
        return ( (int) (next - buffer + size + 4) );
      }

      for( loop = 0; loop < 8; loop++ )
        plplcrcrlf7777[loop] = plplcrcrlf7777[loop+1];
    }

    if( !endFound ) {
      if( feof((FILE *) stream) ) {
        printf("crex_size: end-of-file hit before end of CREX found\n");
        exit(1);
      }
      else {
        size += (BUFFLEN - size7777);
        memmove(buffer, (buffer+(BUFFLEN-size7777)), size7777);
        inBuffer = fileRead((buffer+size7777),
                            (BUFFLEN-size7777),
                            (FILE *) stream);
        if( ferror((FILE *) stream) ) {
          perror("crex_size: file read error");
          exit(1);
        }
        if( inBuffer == 0 ) return 0;
        next = buffer;
        endBuffer = next + abs(inBuffer);
      }
    }
  }

/*
// Position file where it started
*/
  status = fileSeek(stream, startPosition, SEEK_SET);
  if( status != 0 ) {
    perror("crex_size: file repositioning error");
    exit(1);
  }
  return 0;
}
