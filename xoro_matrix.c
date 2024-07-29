#include "xoro_matrix.h"
#include "stdbool.h"

#include "stdio.h"
#define DEBUG_ON true
#define DEBUG(...) { if (DEBUG_ON) { printf(__VA_ARGS__); } }

// -----------------------------------------------------------

void copyFXTMatrix(const FXTMatrix* from, FXTMatrix* to)
{
    for (int qi = 0; qi < 2; qi++)
    {
        for (int qj = 0; qj < 2; qj++)
        {
            for (int k = 0; k < 64; k++)
                (to->M)[qi][qj][k] = (from->M)[qi][qj][k];
        }
    }
}

// -----------------------------------------------------------

/*
void transformToDiagonal(const FXTMatrix* matrix, FXTDiagMatrix* result)
{
    int iDiag = 0;

    // create lower diagonals (without main)
    for (int i0 = 127; i0 > 0; i0--, iDiag++)
    {
        uint64_t mask = ONLY_TOP_BIT;
        int i = i0;
        int j = 0;
        int col = 0;
        while (i <= 127)
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
*/

void transposeFXTM(const FXTMatrix* matrix, FXTMatrix* transposed)
{
    for (int j = 0; j < 64; j++)
    {
        const int sh = 63-j;
        uint64_t val = 0ULL;

        for (int i = 0; i < 64; i++)
        {
            val <<= 1;
            val |= ((matrix->M)[i][0] >> sh) & 1ULL;
        }
        (transposed->M)[j][0] = val;

        val = 0ULL;
        for (int i = 64; i < 128; i++)
        {
            val <<= 1;
            val |= ((matrix->M)[i][0] >> sh) & 1ULL;
        }
        (transposed->M)[j][1] = val;
    }

    for (int j = 0; j < 64; j++)
    {
        const int sh = 63-j;
        uint64_t val = 0ULL;

        for (int i = 0; i < 64; i++)
        {
            val <<= 1;
            val |= ((matrix->M)[i][1] >> sh) & 1ULL;
        }
        (transposed->M)[j+64][0] = val;

        val = 0ULL;
        for (int i = 64; i < 128; i++)
        {
            val <<= 1;
            val |= ((matrix->M)[i][1] >> sh) & 1ULL;
        }
        (transposed->M)[j+64][1] = val;
    }
}

// aT should be the FXTM trasnposition of b
void multiplyFTXM(const FXTMatrix* aT, const FXTMatrix* b, FXTMatrix* c)
{
    // main diagonal
}

void fastXoroMatrixPower(const FXTMatrix* matrix, uint64_t power, FXTMatrix* result)
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
        bool isDiagonalReady = false;

        if (power & 1ULL)
        {
            DEBUG("power = %llu\n", power);
            if (isResultZero)
            {
                copyFXTMatrix(currentPower, nextResult);
                DEBUG("copied successfully\n");
                isResultZero = false;
            }
            else
            {
                transformToDiagonal(currentPower, &diag);
                DEBUG("else: diagonalized successfully\n");
                isDiagonalReady = true;
                fastXoroMatrixMul(currentResult, &diag, nextResult);
                DEBUG("else: performed fast mul successfully\n");
            }
                
            // swap next and current result, clearing space for the next operations
            temp = currentResult;
            currentResult = nextResult;
            nextResult = temp;
        }

        // only calculate next power when needed
        if (power >>= 1)
        {
            // make a diagonal matrix for the current power (if needed)
            if (!isDiagonalReady) {
                DEBUG("power = %llu", power)
                transformToDiagonal(currentPower, &diag);
                for (int i = 0; i < 255; i++)
                {
                    printf("%-4llx %-4llx\n", (diag.M)[i][0], (diag.M)[i][1]);
                }
            }

            // perform a very fast, quadratic matrix multiplication
            fastXoroMatrixMul(currentPower, &diag, nextPower);

            // swap next and current power, making space for the next operations
            temp = currentPower;
            currentPower = nextPower;
            nextPower = temp;
        }
    }

    // copy result to output
    copyFXTMatrix(currentResult, result);

    // no dynamically allocated data, C will handle the
    // struct deallocations on its own
}

// ---------------------------------------
#include "stdio.h"

// TODO optimize the hell out of this
void advanceXoroshiroFXTM(Xoroshiro *state, const FXTMatrix* fxtm)
{
    uint64_t newLo = 0ULL, newHi = 0ULL;

    for (int j = 0; j < 64; j++)
    {
        for (int i = 0; i < 64; i++)
        {
            uint64_t bit1 = ((fxtm->M)[i][0] >> j) & 1ULL;
            uint64_t bit2 = (state->lo >> (63 - i)) & 1ULL;
            newLo ^= (bit1 & bit2) << j;
        }
        for (int i = 0; i < 64; i++)
        {
            uint64_t bit1 = ((fxtm->M)[i+64][0] >> j) & 1ULL;
            uint64_t bit2 = (state->hi >> (63 - i)) & 1ULL;
            newLo ^= (bit1 & bit2) << j;
        }
    }

    for (int j = 0; j < 64; j++)
    {
        for (int i = 0; i < 64; i++)
        {
            uint64_t bit1 = ((fxtm->M)[i][1] >> j) & 1ULL;
            uint64_t bit2 = (state->lo >> (63 - i)) & 1ULL;
            newHi ^= (bit1 & bit2) << j;
        }
        for (int i = 0; i < 64; i++)
        {
            uint64_t bit1 = ((fxtm->M)[i+64][1] >> j) & 1ULL;
            uint64_t bit2 = (state->hi >> (63 - i)) & 1ULL;
            newHi ^= (bit1 & bit2) << j;
        }
    }

    state->lo = newLo;
    state->hi = newHi;
}