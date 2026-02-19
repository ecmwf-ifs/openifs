/*
 * SHA-256 hash in C
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

#include <stdint.h>


void sha256_compress(uint32_t state[8], const uint8_t block[64]) {
	// 32-bit right rotation
	#define ROR(x, i)  \
		(((x) << (32 - (i))) | ((x) >> (i)))
	
	#define LOADSCHEDULE(i)  \
		schedule[i] =                           \
			  (uint32_t)block[i * 4 + 0] << 24  \
			| (uint32_t)block[i * 4 + 1] << 16  \
			| (uint32_t)block[i * 4 + 2] <<  8  \
			| (uint32_t)block[i * 4 + 3];
	
	#define SCHEDULE(i)  \
		schedule[i] = schedule[i - 16] + schedule[i - 7]  \
			+ (ROR(schedule[i - 15], 7) ^ ROR(schedule[i - 15], 18) ^ (schedule[i - 15] >> 3))  \
			+ (ROR(schedule[i - 2], 17) ^ ROR(schedule[i - 2], 19) ^ (schedule[i - 2] >> 10));
	
	#define ROUND(a, b, c, d, e, f, g, h, i, k) \
		h += (ROR(e, 6) ^ ROR(e, 11) ^ ROR(e, 25)) + (g ^ (e & (f ^ g))) + UINT32_C(k) + schedule[i];  \
		d += h;  \
		h += (ROR(a, 2) ^ ROR(a, 13) ^ ROR(a, 22)) + ((a & (b | c)) | (b & c));
	
	uint32_t schedule[64];
	LOADSCHEDULE( 0)
	LOADSCHEDULE( 1)
	LOADSCHEDULE( 2)
	LOADSCHEDULE( 3)
	LOADSCHEDULE( 4)
	LOADSCHEDULE( 5)
	LOADSCHEDULE( 6)
	LOADSCHEDULE( 7)
	LOADSCHEDULE( 8)
	LOADSCHEDULE( 9)
	LOADSCHEDULE(10)
	LOADSCHEDULE(11)
	LOADSCHEDULE(12)
	LOADSCHEDULE(13)
	LOADSCHEDULE(14)
	LOADSCHEDULE(15)
	SCHEDULE(16)
	SCHEDULE(17)
	SCHEDULE(18)
	SCHEDULE(19)
	SCHEDULE(20)
	SCHEDULE(21)
	SCHEDULE(22)
	SCHEDULE(23)
	SCHEDULE(24)
	SCHEDULE(25)
	SCHEDULE(26)
	SCHEDULE(27)
	SCHEDULE(28)
	SCHEDULE(29)
	SCHEDULE(30)
	SCHEDULE(31)
	SCHEDULE(32)
	SCHEDULE(33)
	SCHEDULE(34)
	SCHEDULE(35)
	SCHEDULE(36)
	SCHEDULE(37)
	SCHEDULE(38)
	SCHEDULE(39)
	SCHEDULE(40)
	SCHEDULE(41)
	SCHEDULE(42)
	SCHEDULE(43)
	SCHEDULE(44)
	SCHEDULE(45)
	SCHEDULE(46)
	SCHEDULE(47)
	SCHEDULE(48)
	SCHEDULE(49)
	SCHEDULE(50)
	SCHEDULE(51)
	SCHEDULE(52)
	SCHEDULE(53)
	SCHEDULE(54)
	SCHEDULE(55)
	SCHEDULE(56)
	SCHEDULE(57)
	SCHEDULE(58)
	SCHEDULE(59)
	SCHEDULE(60)
	SCHEDULE(61)
	SCHEDULE(62)
	SCHEDULE(63)
	
	uint32_t a = state[0];
	uint32_t b = state[1];
	uint32_t c = state[2];
	uint32_t d = state[3];
	uint32_t e = state[4];
	uint32_t f = state[5];
	uint32_t g = state[6];
	uint32_t h = state[7];
	ROUND(a, b, c, d, e, f, g, h,  0, 0x428A2F98)
	ROUND(h, a, b, c, d, e, f, g,  1, 0x71374491)
	ROUND(g, h, a, b, c, d, e, f,  2, 0xB5C0FBCF)
	ROUND(f, g, h, a, b, c, d, e,  3, 0xE9B5DBA5)
	ROUND(e, f, g, h, a, b, c, d,  4, 0x3956C25B)
	ROUND(d, e, f, g, h, a, b, c,  5, 0x59F111F1)
	ROUND(c, d, e, f, g, h, a, b,  6, 0x923F82A4)
	ROUND(b, c, d, e, f, g, h, a,  7, 0xAB1C5ED5)
	ROUND(a, b, c, d, e, f, g, h,  8, 0xD807AA98)
	ROUND(h, a, b, c, d, e, f, g,  9, 0x12835B01)
	ROUND(g, h, a, b, c, d, e, f, 10, 0x243185BE)
	ROUND(f, g, h, a, b, c, d, e, 11, 0x550C7DC3)
	ROUND(e, f, g, h, a, b, c, d, 12, 0x72BE5D74)
	ROUND(d, e, f, g, h, a, b, c, 13, 0x80DEB1FE)
	ROUND(c, d, e, f, g, h, a, b, 14, 0x9BDC06A7)
	ROUND(b, c, d, e, f, g, h, a, 15, 0xC19BF174)
	ROUND(a, b, c, d, e, f, g, h, 16, 0xE49B69C1)
	ROUND(h, a, b, c, d, e, f, g, 17, 0xEFBE4786)
	ROUND(g, h, a, b, c, d, e, f, 18, 0x0FC19DC6)
	ROUND(f, g, h, a, b, c, d, e, 19, 0x240CA1CC)
	ROUND(e, f, g, h, a, b, c, d, 20, 0x2DE92C6F)
	ROUND(d, e, f, g, h, a, b, c, 21, 0x4A7484AA)
	ROUND(c, d, e, f, g, h, a, b, 22, 0x5CB0A9DC)
	ROUND(b, c, d, e, f, g, h, a, 23, 0x76F988DA)
	ROUND(a, b, c, d, e, f, g, h, 24, 0x983E5152)
	ROUND(h, a, b, c, d, e, f, g, 25, 0xA831C66D)
	ROUND(g, h, a, b, c, d, e, f, 26, 0xB00327C8)
	ROUND(f, g, h, a, b, c, d, e, 27, 0xBF597FC7)
	ROUND(e, f, g, h, a, b, c, d, 28, 0xC6E00BF3)
	ROUND(d, e, f, g, h, a, b, c, 29, 0xD5A79147)
	ROUND(c, d, e, f, g, h, a, b, 30, 0x06CA6351)
	ROUND(b, c, d, e, f, g, h, a, 31, 0x14292967)
	ROUND(a, b, c, d, e, f, g, h, 32, 0x27B70A85)
	ROUND(h, a, b, c, d, e, f, g, 33, 0x2E1B2138)
	ROUND(g, h, a, b, c, d, e, f, 34, 0x4D2C6DFC)
	ROUND(f, g, h, a, b, c, d, e, 35, 0x53380D13)
	ROUND(e, f, g, h, a, b, c, d, 36, 0x650A7354)
	ROUND(d, e, f, g, h, a, b, c, 37, 0x766A0ABB)
	ROUND(c, d, e, f, g, h, a, b, 38, 0x81C2C92E)
	ROUND(b, c, d, e, f, g, h, a, 39, 0x92722C85)
	ROUND(a, b, c, d, e, f, g, h, 40, 0xA2BFE8A1)
	ROUND(h, a, b, c, d, e, f, g, 41, 0xA81A664B)
	ROUND(g, h, a, b, c, d, e, f, 42, 0xC24B8B70)
	ROUND(f, g, h, a, b, c, d, e, 43, 0xC76C51A3)
	ROUND(e, f, g, h, a, b, c, d, 44, 0xD192E819)
	ROUND(d, e, f, g, h, a, b, c, 45, 0xD6990624)
	ROUND(c, d, e, f, g, h, a, b, 46, 0xF40E3585)
	ROUND(b, c, d, e, f, g, h, a, 47, 0x106AA070)
	ROUND(a, b, c, d, e, f, g, h, 48, 0x19A4C116)
	ROUND(h, a, b, c, d, e, f, g, 49, 0x1E376C08)
	ROUND(g, h, a, b, c, d, e, f, 50, 0x2748774C)
	ROUND(f, g, h, a, b, c, d, e, 51, 0x34B0BCB5)
	ROUND(e, f, g, h, a, b, c, d, 52, 0x391C0CB3)
	ROUND(d, e, f, g, h, a, b, c, 53, 0x4ED8AA4A)
	ROUND(c, d, e, f, g, h, a, b, 54, 0x5B9CCA4F)
	ROUND(b, c, d, e, f, g, h, a, 55, 0x682E6FF3)
	ROUND(a, b, c, d, e, f, g, h, 56, 0x748F82EE)
	ROUND(h, a, b, c, d, e, f, g, 57, 0x78A5636F)
	ROUND(g, h, a, b, c, d, e, f, 58, 0x84C87814)
	ROUND(f, g, h, a, b, c, d, e, 59, 0x8CC70208)
	ROUND(e, f, g, h, a, b, c, d, 60, 0x90BEFFFA)
	ROUND(d, e, f, g, h, a, b, c, 61, 0xA4506CEB)
	ROUND(c, d, e, f, g, h, a, b, 62, 0xBEF9A3F7)
	ROUND(b, c, d, e, f, g, h, a, 63, 0xC67178F2)
	state[0] += a;
	state[1] += b;
	state[2] += c;
	state[3] += d;
	state[4] += e;
	state[5] += f;
	state[6] += g;
	state[7] += h;
}

