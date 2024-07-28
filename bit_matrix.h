#pragma once
#include "inttypes.h"
#include "stdlib.h"

// -------------------------------------------------------------------------------

enum ErrType {
    SUCCESS = 0,
    BAD_MATRIX_SIZE,
    BAD_MATRIX_SHAPE,
    OUT_OF_BOUNDS,
    MATRIX_UNINITIALIZED,
    UNEXPECTED_ERROR
};
typedef enum ErrType ErrType;

const char* getErrorMessage(ErrType err);

// -------------------------------------------------------------------------------

typedef struct BitMatrix BitMatrix;
struct BitMatrix 
{
    int sizeI, sizeJ;
    int **M;
};

void allocMatrix(BitMatrix *mx, int sizeI, int sizeJ);
void deallocMatrix(BitMatrix *mx);

ErrType makeZeroMatrix(BitMatrix *mx);
ErrType makeUnitMatrix(BitMatrix *mx);
ErrType makeLeftShiftMatrix(BitMatrix *mx, int shift);
ErrType makeLeftRotationMatrix(BitMatrix *mx, int rot);

void printMatrix(BitMatrix *mx);
ErrType matrixInject(BitMatrix *injected, BitMatrix *target, int iPos, int jPos);
ErrType matrixSum(BitMatrix *a, BitMatrix *b, BitMatrix *sum);
ErrType matrixTimesMatrix(BitMatrix *a, BitMatrix *b, BitMatrix *product);