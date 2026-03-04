/* 
 * SHA-256 hash in C and x86 assembly
 * 
 * Copyright (c) 2014 Project Nayuki
 * http://www.nayuki.io/page/fast-sha2-hashes-in-x86-assembly
 * 
 * (MIT License)
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */
 
// 2-Sep-2015   Alan Geer (ECMWF)   Allow discontinuous streams of data 
// 5-Aug-2016   Alan Geer (ECMWF)   Fix undersized integers that allowed infinite loop

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

extern void sha256_compress(uint32_t state[8], const uint8_t block[64]);

/* Full message hasher */

void sha256_hash(const uint8_t *message, int64_t len, uint32_t hash[8], 
  uint8_t block[64], int64_t *len_block, uint32_t action) {
  
/* Actions (enumeration interoperabiltiy with fortran is iffy, so just use numbers):
 *   0 - Complete message to be processed in one call
 *   1 - Start processing incomplete message (state and len_state are kept by user between calls)
 *   2 - Continue processing incomplete message
 *   3 - Complete processing 
 */
  int64_t istart;
  int64_t rem;
  int64_t i;

  if (action <= 1) {

    // Initialise the hash
    hash[0] = UINT32_C(0x6A09E667);
    hash[1] = UINT32_C(0xBB67AE85);
    hash[2] = UINT32_C(0x3C6EF372);
    hash[3] = UINT32_C(0xA54FF53A);
    hash[4] = UINT32_C(0x510E527F);
    hash[5] = UINT32_C(0x9B05688C);
    hash[6] = UINT32_C(0x1F83D9AB);
    hash[7] = UINT32_C(0x5BE0CD19);
    istart = 0;
    *len_block = 0;
    rem = len;

  } else {

    istart = 64 - *len_block;
    rem = len - istart;
    if( rem >= 0) {

      // We now have a full block to be hashed
      memcpy(block+*len_block, message, istart);
      sha256_compress(hash, block);
      *len_block = 0;

    } else {

      // Only a partial block is available as yet
      memcpy(block+*len_block, message, len);
      *len_block += len;

    }

  }

  if(rem > 0) {

    // All complete blocks in the message
    for (i = istart; len - i >= 64; i += 64)
      sha256_compress(hash, message + i);

    // Remaining bit of message
    rem = len - i;
    memcpy(block, message + i, rem);
    *len_block = rem;

  } else {

    rem = *len_block;

  }

  if (action==0 || action==3) {

    // Hash padding and completion
    block[rem] = 0x80;
    rem++;
    if (64 - rem >= 8)
      memset(block + rem, 0, 56 - rem);
    else {
      memset(block + rem, 0, 64 - rem);
      sha256_compress(hash, block);
      memset(block, 0, 56);
    }

    uint64_t longLen = ((uint64_t)len) << 3;
    for (i = 0; i < 8; i++)
      block[64 - 1 - i] = (uint8_t)(longLen >> (i * 8));
    sha256_compress(hash, block);

  }
}

