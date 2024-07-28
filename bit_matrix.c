#include "bit_matrix.h"
#include "stdio.h"

const char* getErrorMessage(ErrType err) {
    static const char* MSG[] = {
        "Operation successful",
        "ERR: Bad matrix size",
        "ERR: Bad matrix shape",
        "ERR: Out of matrix bounds",
        "ERR: Matrix uninitialized",
        "ERR: An unexpected error occurred"
    };

    return MSG[(int)err];
}

// ---------------------------------------------------------
// Bit matrix allocation
// ---------------------------------------------------------

void allocMatrix(BitMatrix *mx, int sizeI, int sizeJ)
{
    mx->sizeI = sizeI;
    mx->sizeJ = sizeJ;

    mx->M = (int**)malloc(sizeof(int*) * sizeI);
    for (int i = 0; i < sizeI; i++)
    {
        (mx->M)[i] = (int*)malloc(sizeof(int) * sizeJ);
    }
}

void deallocMatrix(BitMatrix *mx)
{
    for (int i = 0; i < mx->sizeI; i++)
        free((mx->M)[i]);
    free(mx->M);
    mx->sizeI = 0;
    mx->sizeJ = 0;
}

// ---------------------------------------------------------
// Initializers for all the different types of matrices
// ---------------------------------------------------------

// 0
ErrType makeZeroMatrix(BitMatrix* mx) 
{
    if (mx->sizeI <= 0 || mx->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (mx->sizeI != mx->sizeJ)
        return BAD_MATRIX_SHAPE;

    for (int i = 0; i < mx->sizeI; i++)
    {
        for (int j = 0; j < mx->sizeJ; j++)
        {
            (mx->M)[i][j] = 0;
        }
    }

    return SUCCESS;
}

// I
ErrType makeUnitMatrix(BitMatrix* mx) 
{
    if (mx->sizeI <= 0 || mx->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (mx->sizeI != mx->sizeJ)
        return BAD_MATRIX_SHAPE;

    for (int i = 0; i < mx->sizeI; i++)
    {
        for (int j = 0; j < mx->sizeJ; j++)
        {
            (mx->M)[i][j] = (i == j ? 1 : 0);
        }
    }

    return SUCCESS;
}

// S(a) - left-shift by a bits
ErrType makeLeftShiftMatrix(BitMatrix* mx, int shift) 
{
    if (mx->sizeI <= 0 || mx->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (mx->sizeI != mx->sizeJ)
        return BAD_MATRIX_SHAPE;

    (void)makeZeroMatrix(mx);

    for (int i = shift, j = 0; i < mx->sizeI; i++, j++)
    {
        (mx->M)[i][j] = 1;
    }

    return SUCCESS;
}

// R(a) - left-rotate by a bits
ErrType makeLeftRotationMatrix(BitMatrix* mx, int rot) 
{
    if (mx->sizeI <= 0 || mx->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (mx->sizeI != mx->sizeJ)
        return BAD_MATRIX_SHAPE;

    rot %= mx->sizeI;
    (void)makeZeroMatrix(mx);

    for (int i = rot, j = 0; j < mx->sizeJ; i++, j++)
    {
        if (i == mx->sizeI)
            i -= mx->sizeI;
        (mx->M)[i][j] = 1;
    }

    return SUCCESS;
}

// ---------------------------------------------------------
// Matrix operations and utilities
// ---------------------------------------------------------

void printMatrix(BitMatrix* mx)
{
    printf("Size = (%d, %d)\n", mx->sizeI, mx->sizeJ);
    for (int i = 0; i < mx->sizeI; i++)
    {
        for (int j = 0; j < mx->sizeJ; j++)
        {
            printf("%d ", (mx->M)[i][j]);
        }
        printf("\n");
    }
}

ErrType matrixInject(BitMatrix* injected, BitMatrix* target, int iPos, int jPos) 
{
    if (injected->sizeI <= 0 || injected->M == NULL || target->sizeI <= 0 || target->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (target->sizeI < injected->sizeI || target->sizeJ < injected->sizeJ)
        return BAD_MATRIX_SIZE;

    for (int i = 0; i < injected->sizeI; i++)
    {
        for (int j = 0; j < injected->sizeJ; j++)
        {
            const int iTarget = i + iPos;
            const int jTarget = j + jPos;
            if (iTarget >= target->sizeI || jTarget >= target->sizeJ)
                return OUT_OF_BOUNDS;

            (target->M)[iTarget][jTarget] = (injected->M)[i][j];
        }
    }

    return SUCCESS;
}

// A + B
ErrType matrixSum(BitMatrix* a, BitMatrix* b, BitMatrix* sum) 
{
    if (a->sizeI <= 0 || a->M == NULL || b->sizeI <= 0 || b->M == NULL || sum->sizeI <= 0 || sum->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (a->sizeI != b->sizeI || b->sizeI != sum->sizeI || a->sizeJ != b->sizeJ || b->sizeJ != sum->sizeJ)
        return BAD_MATRIX_SIZE;

    for (int i = 0; i < a->sizeI; i++)
    {
        for (int j = 0; j < a->sizeJ; j++)
        {
            const int newBit = (a->M)[i][j] ^ (b->M)[i][j];
            (sum->M)[i][j] = newBit;
        }
    }

    return SUCCESS;
}

// A x B
ErrType matrixTimesMatrix(BitMatrix* a, BitMatrix* b, BitMatrix* product) 
{
    if (a->sizeI <= 0 || a->M == NULL || b->sizeI <= 0 || b->M == NULL || product->sizeI <= 0 || product->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (a->sizeJ != b->sizeI || a->sizeI != product->sizeI || b->sizeJ != product->sizeJ)
        return BAD_MATRIX_SIZE;

    int** aM = a->M;
    int** bM = b->M;

    for (int i = 0; i < product->sizeI; i++)
    {
        for (int j = 0; j < product->sizeJ; j++)
        {
            int newBit = 0;
            for (int k = 0; k < a->sizeJ; k++)
                newBit ^= (aM[i][k] & bM[k][j]);
            
            (product->M)[i][j] = newBit;
        }
    }

    return SUCCESS;
}
