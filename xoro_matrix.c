#include "xoro_matrix.h"
#include "stdbool.h"

// -----------------------------------------------------------

void copyFXTMatrix(const FXTMatrix* from, FXTMatrix* to)
{
    for (int i = 0; i < 128; i++)
    {
        (to->M)[i][0] = (from->M)[i][0];
        (to->M)[i][1] = (from->M)[i][1];
    }
}

// -----------------------------------------------------------

void transformToDiagonal(const FXTMatrix* matrix, FXTDiagMatrix* result)
{
    static const uint64_t ONLY_TOP_BIT = 0x8000'0000'0000'0000ULL;
    int iDiag = 0;

    // create lower diagonals (without main)
    for (int i0 = 127; i0 > 0; i0--, iDiag++)
    {
        uint64_t mask = ONLY_TOP_BIT;
        int i = i0;
        int j = 0;
        int col = 0;
        while (i >= 127)
        {
            if (j >= 64)
            {
                mask = ONLY_TOP_BIT;
                j -= 64;
                col++;
            }

            (result->M)[iDiag][col] |= (mask & (matrix->M)[i][col]);

            i++;
            j++;
            mask >>= 1;
        }
    }

    // create upper diagonals (with main)
    for (int j0 = 0; j0 < 128; j0++, iDiag++)
    {
        uint64_t mask = ONLY_TOP_BIT >> (j0 & 0x3f);
        int i = 0;
        int j = j0;
        int col = (j0 >= 64 ? 1 : 0);
        while (col == 0 || j < 64)
        {
            if (j >= 64)
            {
                mask = ONLY_TOP_BIT;
                j -= 64;
                col++;
            }

            (result->M)[iDiag][col] |= (mask & (matrix->M)[i][col]);

            i++;
            j++;
            mask >>= 1;
        }
    }
}

void fastXoroMatrixMul(const FXTMatrix* a, const FXTDiagMatrix* b, FXTMatrix* product)
{
    for (int iProduct = 0; iProduct < 128; iProduct++)
    {
        uint64_t prod0 = 0ULL, prod1 = 0ULL;

        // hardcoded loops over areas that have data for a
        // particular column, saves a couple of useless iterations
        for (int iDiag = 0; iDiag < 255 - 63; iDiag++)
            prod0 ^= (a->M)[iProduct][0] & (b->M)[iDiag][0];
            
        for (int iDiag = 63; iDiag < 255; iDiag++)
            prod1 ^= (a->M)[iProduct][1] & (b->M)[iDiag][1];

        (product->M)[iProduct][0] = prod0;
        (product->M)[iProduct][1] = prod1;
    }
}

void xoroMatrixFastPower(const FXTMatrix* matrix, uint64_t power, FXTMatrix* result)
{
    FXTMatrix res1 = { 0 }, res2 = { 0 };
    FXTMatrix pow1 = { 0 }, pow2 = { 0 };
    FXTDiagMatrix diag = { 0 };
    bool isResultZero = true;

    FXTMatrix *currentPower = &pow1, *nextPower = &pow2, *currentResult = &res1, *nextResult = &res2;
    FXTMatrix *temp;
    copyFXTMatrix(matrix, currentPower);

    while (power)
    {
        if (power & 1ULL)
        {
            if (isResultZero)
            {
                copyFXTMatrix(currentPower, nextResult);
                isResultZero = false;
            }
            else
                fastXoroMatrixMul(currentResult, currentPower, nextResult);

            // swap next and current result, clearing space for the next operations
            temp = currentResult;
            currentResult = nextResult;
            nextResult = temp;
        }

        // make a diagonal matrix for the current power
        transformToDiagonal(currentPower, &diag);
        // perform a very fast, quadratic matrix multiplication
        fastXoroMatrixMul(currentPower, &diag, nextPower);

        // swap next and current power, making space for the next operations
        temp = currentPower;
        currentPower = nextPower;
        nextPower = temp;
        power >>= 1;
    }

    // copy result to output
    copyFXTMatrix(currentResult, result);

    // no dynamically allocated data, C will handle the
    // struct deallocations on its own
}
