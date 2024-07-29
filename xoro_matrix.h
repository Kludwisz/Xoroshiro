#pragma once
#include "inttypes.h"

// structures for the fast xoroshiro transform
// matrix multiplication procedure

typedef struct FXTMatrix FXTMatrix;
struct FXTMatrix {
    uint64_t M[128][2];
};

typedef struct FXTDiagMatrix FXTDiagMatrix;
struct FXTDiagMatrix {
    uint64_t M[128 * 2 - 1][2];
};

