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

#include <stdio.h>
#include <assert.h>
#define EFDEBUG 0
#define fortint int
#define GBYTES gbytes_
#define GBYTE  gbyte_
#define SBYTES sbytes_
#define SBYTE  sbyte_

#define SWORD 32
#define MASKSWORD 0x3f


static int onbit[32] = {
  0x00000001,0x00000002,0x00000004,0x00000008,
  0x00000010,0x00000020,0x00000040,0x00000080,
  0x00000100,0x00000200,0x00000400,0x00000800,
  0x00001000,0x00002000,0x00004000,0x00008000,
  0x00010000,0x00020000,0x00040000,0x00080000,
  0x00100000,0x00200000,0x00400000,0x00800000,
  0x01000000,0x02000000,0x04000000,0x08000000,
  0x10000000,0x20000000,0x40000000,0x80000000};


static unsigned long masks[] = {0,1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535,
                                131071, 262143, 524287, 1048575, 2097151, 4194303, 8388607, 16777215, 33554431,
                                67108863, 134217727, 268435455, 536870911, 1073741823, 2147483647, -1};


void GBYTE(
    void* Source,
    void* Destination,
    int* bitToSkip,
    int* bitsPerValue) {
  /*
  //  GBYTE:
  //
  //   Unpacks one value from source to destination.
  //
  //   Makes an initial skip over bits in source.
  //
  //   startSkip >= 0            Number of bits to be skipped preceding first
  //                             value in source
  //
  //   0 < bitsPerValue < SWORD <= 32  Size in bits of each value
  //
  //
  //  Rewritten to get better preformance on little endian machines.
  //  This code is portable to big endian machines.
  //                   November 2006  Enrico Fucile
  //
  */
  fortint mask,nextWord,nextValueFirstBit,moveBitsLeft,moveBitsRight,nextByte,nextBitOfByte;
  int bytesPerValue,extraBitsPerValue;
  unsigned long i0,i1,i2,i3,temp=0,temp1;
  unsigned char* source = (unsigned char*) Source;

  fortint* destination = (fortint*) Destination;
  nextBitOfByte = *bitToSkip % 8;
  extraBitsPerValue = *bitsPerValue % 8;

#if EFDEBUG
  printf("---- nextBitOfByte = %d\n",nextBitOfByte);
  if (nextBitOfByte != 0) printf("---- (not zero) nextBitOfByte = %d\n",nextBitOfByte);
  printf("---- bitsPerValue= %d ---\n",*bitsPerValue);
#endif

  if (*bitsPerValue > 32) {
    int byte, bitShift, bit;
    fortint t = 0;
    int bpv = *bitsPerValue;
    int nextBit = *bitToSkip;

    for ( bit = 0; bit < bpv; bit++ ) {
      byte     = ( nextBit ) >> 3;
      bitShift = 7 - ( nextBit++ ) & 0x7;

      temp <<= 1 ;
      if( (source[byte] & onbit[bitShift]) ) temp |= 1;
    }

    *destination = temp;
    return;
  }

  if ( extraBitsPerValue==0 && nextBitOfByte == 0 ) {
    nextByte = *bitToSkip / 8;
#if EFDEBUG
    printf("---- nextByte = %d\n",nextByte);
    printf("---- bytesPerValue = %d\n",bytesPerValue);
    fflush(stdout);
#endif
    switch (*bitsPerValue) {
    case 8:
      *destination=(fortint) source[nextByte];
      break;
    case 16:
      *destination=(fortint) ((source[nextByte]<<8) + source[nextByte+1]);
      break;
    case 24:
      *destination=(fortint) ((source[nextByte]<<16) + (source[nextByte+1]<<8) +source[nextByte+2]);
      break;
    case 32:
      *destination=(fortint) ((source[nextByte]<<24) + (source[nextByte+1]<<16) +(source[nextByte+2]<<8)+source[nextByte+3]);
      break;
    default:
      destination=0;
    }
#if EFDEBUG
    printf("-------------------- destination=%ld\n",(long)*destination);
#endif
  } else {
    mask=masks[*bitsPerValue];
    nextWord = *bitToSkip / SWORD;
    nextByte = nextWord*4;
    nextValueFirstBit = *bitToSkip % SWORD;
    moveBitsLeft = 0;
    moveBitsRight = SWORD - *bitsPerValue - nextValueFirstBit;
#if EFDEBUG
    printf("---- nextWord = %d\n",nextWord);
    printf("---- nextValueFirstBit = %d\n",nextValueFirstBit);
    printf("---- moveBitsRight= %d\n",moveBitsRight);
#endif
    if (moveBitsRight > 0) {

      i0=(unsigned long)(source[nextByte]) << 24;
      i1=(unsigned long)(source[nextByte+1]) << 16;
      i2=(unsigned long)(source[nextByte+2]) << 8;
      i3=(unsigned long)(source[nextByte+3]) ;
      *destination = (fortint)( ((i0 + i1 + i2 + i3) >> moveBitsRight) & mask );
#if EFDEBUG
      printf("-------------------- destination=%ld\n",(char)*destination,(long)*destination);
#endif

    } else if ( moveBitsRight < 0 ) {

      moveBitsLeft=-moveBitsRight;
      moveBitsRight=SWORD-moveBitsLeft;
      i0=(unsigned long)(source[nextByte]) << 24;
      i1=(unsigned long)(source[nextByte+1]) << 16;
      i2=(unsigned long)(source[nextByte+2]) << 8;
      i3=(unsigned long)(source[nextByte+3]) ;
      temp1 = (unsigned long)(i0 + i1 + i2 + i3) << moveBitsLeft ;

      i0=(unsigned long)(source[nextByte+4]) << 24;
      i1=(unsigned long)(source[nextByte+5]) << 16;
      i2=(unsigned long)(source[nextByte+6]) << 8;
      i3=(unsigned long)(source[nextByte+7]) ;
      temp= (unsigned long)(i0 + i1 + i2 + i3) >> moveBitsRight;

      *destination = (fortint) (temp1 | temp);
      *destination &= mask;
#if EFDEBUG
      printf("+++++++++ destination=%ld temp1=%ld temp=%ld mask=%ld bitsPerValue=%d moveBitsRight=%d moveBitsLeft=%d\n",*destination,temp1,temp,mask,*bitsPerValue,moveBitsRight,moveBitsLeft);
#endif

    } else {

      i0=(fortint)(source[nextByte]) << 24;
      i1=(fortint)(source[nextByte+1]) << 16;
      i2=(fortint)(source[nextByte+2]) << 8;
      i3=(fortint)(source[nextByte+3]) ;
      *destination = ( i0 + i1 + i2 + i3 ) & mask ;
#if EFDEBUG
      printf("-------------------- destination=%ld\n",(char)*destination,(long)*destination);
#endif
    }
  }
}


#if 0
void GBYTE(
    void* Source,
    void* Destination,
    int* startSkip,
    int* bitsPerValue) {
  /*
  //  GBYTE:
  //
  //   Unpacks one value from source to destination.
  //
  //   Makes an initial skip over bits in source.
  //
  //   startSkip >= 0            Number of bits to be skipped preceding first
  //                             value in source
  //
  //   0 < bitsPerValue < SWORD  Size in bits of each value
  //
  // Move one bit at a time (!!)
  // Uses shifts and masks instead of divides.
  // Sets destination bit to 0/1 as source bit is 0/1.
  */
  int bpv, temp, nextBit, byte, bitShift, bit;
  unsigned char* source = (unsigned char*) Source;
  int* destination = (int*) Destination;

  temp = 0;
  bpv = *bitsPerValue;
  nextBit = *startSkip;

  for ( bit = 0; bit < bpv; bit++ ) {
    byte     = ( nextBit ) >> 3;
    bitShift = 7 - ( nextBit++ ) & 0x7;

    temp <<= 1 ;
    if( (source[byte] & onbit[bitShift]) ) temp |= 1;
  }

  *destination = temp;
}
#endif


/*
//  Retrieve or store arbitrary bit-size values from or to SWORD-bit words
//
//  Rewritten April 2000, J.D.Chambers, ECMWF.
//
*/


void GBYTES(
    void* Source,
    void* Destination,
    int* startSkip,
    int* bitsPerValue,
    int* skipBetweenValues,
    int* numberOfValues) {
  /*
  //  GBYTES:
  //
  //   Unpacks values from source to destination.
  //
  //   Makes an initial skip over bits in source, then skips bits between values.
  //
  //   startSkip >= 0            Number of bits to be skipped preceding first
  //                             value in source
  //
  //   0 < bitsPerValue < SWORD  Size in bits of each value
  //
  //   skipBetweenValues >= 0    Number of bits to be skipped between values
  //
  //   numberOfValues >= 0       Number of values to be packed/unpacked
  //
  // Move one bit at a time (!!)
  // Uses shifts and masks instead of divides.
  // Sets destination bit to 0/1 as source bit is 0/1.
  */
  unsigned char* source = (unsigned char*) Source;
  int* destination = (int*) Destination;
  int bpv, bit, byte, bitShift, nextBit, nextValueFirstBit, next, step;

  bpv=*bitsPerValue;
  step = *bitsPerValue + *skipBetweenValues;
  nextValueFirstBit = *startSkip;

  /* assert(*bitsPerValue <= 8*sizeof(*destination)); */
  for (next = 0; next < *numberOfValues; ++next) {
    destination[next] = 0;
    nextBit=nextValueFirstBit;

    /* printf("--- nextValueFirstBit=%d\n",nextValueFirstBit); */
    GBYTE(source,&destination[next],&nextValueFirstBit,&bpv);

    /*
      for ( bit = 0; bit < bpv; bit++ ) {
        byte     = ( nextBit ) >> 3;
        bitShift = 7 - ( nextBit++ ) & 0x7;

        destination[next] <<= 1 ;
        if( (source[byte] & onbit[bitShift]) ) destination[next] |= 1;
      }
      */

    nextValueFirstBit += step;
  }
}


void SBYTES(
    void* Destination,
    void* Source,
    int* startSkip,
    int* bitsPerValue,
    int* skipBetweenValues,
    int* numberOfValues) {
  /*
  //  SBYTES:
  //
  //   Packs values from source to destination.
  //
  //   Makes an initial skip over bits in source, then skips bits between values.
  //
  //   startSkip >= 0            Number of bits to be skipped preceding first
  //                             value in destination
  //
  //   0 < bitsPerValue < SWORD  Size in bits of each value
  //
  //   skipBetweenValues >= 0    Number of bits to be skipped between values
  //
  //   numberOfValues >= 0       Number of values to be packed/unpacked
  //
  // Move one bit at a time (!!)
  // Uses shifts and masks instead of divides.
  // Sets destination bit to 0/1 as source bit is 0/1.
  */
  unsigned char* destination = (unsigned char*) Destination;
  int* source = (int*) Source;
  int bitShift, bitPosition;
  int nextValueFirstBit, bpvm1;
  int nextDestination, nextBit, bit;
  int step, nextValue;
  unsigned char setBit0[8] = {0xfe,0xfd,0xfb,0xf7,0xef,0xdf,0xbf,0x7f};
  unsigned char setBit1[8] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};

  step = *bitsPerValue + *skipBetweenValues;
  nextValueFirstBit = *startSkip;
  bpvm1 = (*bitsPerValue - 1);

  for( nextValue = 0; nextValue < *numberOfValues; ++nextValue) {
    nextBit = nextValueFirstBit;

    for ( bit = 0; bit < *bitsPerValue; bit++ ) {
      nextDestination = ( nextBit ) >> 3;
      bitShift    = (bpvm1 - bit) & MASKSWORD;
      bitPosition = 7 - (( nextBit++ ) & 0x7);

      if( (source[nextValue] & onbit[bitShift]) )
        destination[nextDestination] |= setBit1[bitPosition];
      else
        destination[nextDestination] &= setBit0[bitPosition];
    }

    nextValueFirstBit += step;
  }
}


void SBYTE(
    void* destination,
    void* source,
    int* startSkip,
    int* bitsPerValue) {
  /*
  //  SBYTE:
  //
  //   Packs one value from source to destination.
  //
  //   Makes an initial skip over bits in source.
  //
  //   startSkip >= 0            Number of bits to be skipped preceding first
  //                             value in source
  //
  //   0 < bitsPerValue < SWORD  Size in bits of each value
  //
  */
  int skipBetweenValues=0, numberOfValues=1;

  SBYTES(destination,
         source,
         startSkip,
         bitsPerValue,
         &skipBetweenValues,
         &numberOfValues);
}

