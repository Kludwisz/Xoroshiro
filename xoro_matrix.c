#include "xoro_matrix.h"
#include "stdbool.h"

#include "stdio.h"
#define DEBUG_ON false
#define DEBUG(...) { if (DEBUG_ON) { printf(__VA_ARGS__); } }

// -----------------------------------------------------------

void copyFXTM(const FXTMatrix* from, FXTMatrix* to)
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

void transposeFXTM(const FXTMatrix* matrix, FXTMatrix* transposed)
{
    for (int qi = 0; qi < 2; qi++) for (int qj = 0; qj < 2; qj++)
    {
        for (int j = 0; j < 64; j++)
        {
            const int sh = 63-j;
            uint64_t val = 0ULL;

            for (int i = 0; i < 64; i++)
            {
                val <<= 1;
                val |= ((matrix->M)[qi][qj][i] >> sh) & 1ULL;
            }
            (transposed->M)[qi][qj][j] = val;
        }
    }
}

// -----------------------------------------------------------

// void clearQuadrant(uint64_t* quad)
// {
//     for (int i = 0; i < 64; i++)
//         quad[i] = 0ULL;
// }

// void multiplyQuadrant(const uint64_t firstQuad[], const uint64_t secondQuad[], uint64_t resultQuad[])
// {
//     // TODO
// }

// bT is the FXTM "transposition" of a
void multiplyFXTM(const FXTMatrix* a, const FXTMatrix* bT, FXTMatrix* c)
{
    // 0,0
    // clearQuadrant(c->M[0][0]);
    // multiplyQuadrant(aT->M[0][0], b->M[0][0], c->M[0][0]);
    // multiplyQuadrant(aT->M[0][1], b->M[0][1], c->M[0][0]);
    // // 0,1
    // clearQuadrant(c->M[0][1]);
    // multiplyQuadrant(aT->M[1][0], b->M[0][0], c->M[0][1]);
    // multiplyQuadrant(aT->M[1][1], b->M[0][1], c->M[0][1]);
    // // 1,0
    // clearQuadrant(c->M[1][0]);
    // multiplyQuadrant(aT->M[0][0], b->M[1][0], c->M[1][0]);
    // multiplyQuadrant(aT->M[0][1], b->M[1][1], c->M[1][0]);
    // // 1,1
    // clearQuadrant(c->M[1][1]);
    // multiplyQuadrant(aT->M[1][0], b->M[1][0], c->M[1][1]);
    // multiplyQuadrant(aT->M[1][1], b->M[1][1], c->M[1][1]);

    // calculate each bit in the result matrix separately
    for (int qi = 0; qi < 2; qi++) for (int qj = 0; qj < 2; qj++)
    {
        for (int i = 0; i < 64; i++)
        {
            uint64_t res = 0ULL;
            const uint64_t aLo = a->M[qi][0][i], aHi = a->M[qi][1][i];
            for (int j = 0; j < 64; j++)
            {
                res <<= 1;
                const uint64_t mulXor = (aLo & bT->M[0][qj][j]) ^ (aHi & bT->M[1][qj][j]);
                const int newBit = __builtin_parityll(mulXor); // aka: GCC, do your thing!
                res |= newBit; 
            }

            (c->M)[qi][qj][i] = res;
            DEBUG("%llx\n", res);
        }
    }
}

void fastXoroMatrixPower(const FXTMatrix* matrix, uint64_t power, FXTMatrix* result)
{
    FXTMatrix res1 = { 0 }, res2 = { 0 };
    FXTMatrix pow1 = { 0 }, pow2 = { 0 };
    FXTMatrix transposed = { 0 };
    bool isResultZero = true;

    FXTMatrix *currentPower = &pow1, *nextPower = &pow2, *currentResult = &res1, *nextResult = &res2;
    FXTMatrix *temp;
    copyFXTM(matrix, currentPower);

    while (power)
    {
        bool transpositionDone = false;

        if (power & 1ULL)
        {
            DEBUG("if power & 1:  power = %llu\n", power);
            if (isResultZero)
            {
                copyFXTM(currentPower, nextResult);
                DEBUG("copied successfully\n");
                isResultZero = false;
            }
            else
            {
                transposeFXTM(currentPower, &transposed);
                transpositionDone = true;
                DEBUG("else: transposed successfully\n");
                multiplyFXTM(currentPower, &transposed, nextResult);
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
            // create the transposition of currentPower
            DEBUG("calc next:  power = %llu\n", power)
            if (!transpositionDone)
                transposeFXTM(currentPower, &transposed);

            //for (int i = 0; i < 64; i++)
            //    DEBUG("%llu\n", transposed.M[0][0][i])

            // perform a very fast, quasi-quadratic matrix multiplication
            multiplyFXTM(currentPower, &transposed, nextPower);
            for (int i = 0; i < 64; i++)
                DEBUG("%llu\n", nextPower->M[0][0][i])

            // swap next and current power, making space for the next operations
            temp = currentPower;
            currentPower = nextPower;
            nextPower = temp;
        }
    }

    // copy result to output
    copyFXTM(currentResult, result);

    // no dynamically allocated data, C will handle the
    // struct deallocations on its own
}

// ---------------------------------------

void advanceXoroshiroFXTM(Xoroshiro *state, const FXTMatrix* fxtm)
{
    uint64_t oldState[2] = { state->lo, state->hi };
    uint64_t newState[2] = { 0ULL, 0ULL };

    for (int qi = 0; qi < 2; qi++) for (int qj = 0; qj < 2; qj++)
    {
        for (int i = 0; i < 64; i++)
        {
            const int sh = 63 - i;

            for (int j = 0; j < 64; j++)
            {
                const uint64_t bit1 = ((fxtm->M)[qi][qj][i] >> j) & 1ULL;
                const uint64_t bit2 = (oldState[qi] >> sh) & 1ULL;
                newState[qj] ^= (bit1 & bit2) << j;
            }
        }
    }

    state->lo = newState[0];
    state->hi = newState[1];
}

// --------------------------------------------------------------
// For 128-bit integers
// --------------------------------------------------------------

// void clearXM(XMatrix* xm)
// {
//     for (int i = 0; i < 128; i++)
//         (xm->M)[i] = 0;
// }

// void copyXM(const XMatrix* from, XMatrix* to)
// {
//     for (int k = 0; k < 128; k++)
//         (to->M)[k] = (from->M)[k];
// }

// void transposeXM(const XMatrix* matrix, XMatrix* transposed)
// {
//     for (int j = 0; j < 128; j++)
//     {
//         const int sh = 128-j;
//         uint128_t val = 0;

//         for (int i = 0; i < 128; i++)
//         {
//             val <<= 1;
//             val |= ((matrix->M)[i] >> sh) & 1;
//         }

//         (transposed->M)[j] = val;
//     }
// }

// void multiplyXM(const XMatrix* aT, const XMatrix* b, XMatrix* res)
// {
//     // TODO
// }