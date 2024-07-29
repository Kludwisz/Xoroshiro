#pragma once
#include "xoroshiro.h"
#include "inttypes.h"

#define ONLY_TOP_BIT (0x8000000000000000ULL)

// hacky way to turn off warnings and still be able to compile
// the code just fine even without 128-bit int support
#ifndef __SIZEOF_INT128__
#define __int128 int
#endif
typedef unsigned __int128 uint128_t;

// structures for the fast xoroshiro transform
// matrix multiplication procedure

typedef struct FXTMatrix FXTMatrix;
struct FXTMatrix {
    uint64_t M[2][2][64];
};

typedef struct XMatrix XMatrix;
struct XMatrix {
    uint128_t M[128];
};

// ------------------------------------------------------------------------------------