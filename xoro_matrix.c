#include "xoro_matrix.h"
#include "stdbool.h"

#include "stdio.h"
#define DEBUG_ON true
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

// -----------------------------------------------------------


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

void clearQuadrant(uint64_t* quad)
{
    for (int i = 0; i < 64; i++)
        quad[i] = 0ULL;
}

void multiplyQuadrant(const uint64_t firstQuad[], const uint64_t secondQuad[], uint64_t resultQuad[])
{
    // TODO
}

// aT should be the FXTM trasnposition of b
void multiplyFXTM(const FXTMatrix* aT, const FXTMatrix* b, FXTMatrix* c)
{
    // In the future, if 128-bit integers are supported by more CPUs
    // natively, this entire function could be reduced to just a single 
    // quadrant multiplication, speeding it up significantly.

    // 0,0
    clearQuadrant(c->M[0][0]);
    multiplyQuadrant(aT->M[0][0], b->M[0][0], c->M[0][0]);
    multiplyQuadrant(aT->M[0][1], b->M[0][1], c->M[0][0]);
    // 0,1
    clearQuadrant(c->M[0][1]);
    multiplyQuadrant(aT->M[1][0], b->M[0][0], c->M[0][1]);
    multiplyQuadrant(aT->M[1][1], b->M[0][1], c->M[0][1]);
    // 1,0
    clearQuadrant(c->M[1][0]);
    multiplyQuadrant(aT->M[0][0], b->M[1][0], c->M[1][0]);
    multiplyQuadrant(aT->M[0][1], b->M[1][1], c->M[1][0]);
    // 1,1
    clearQuadrant(c->M[1][1]);
    multiplyQuadrant(aT->M[1][0], b->M[1][0], c->M[1][1]);
    multiplyQuadrant(aT->M[1][1], b->M[1][1], c->M[1][1]);
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
        if (power & 1ULL)
        {
            DEBUG("power = %llu\n", power);
            if (isResultZero)
            {
                copyFXTM(currentPower, nextResult);
                DEBUG("copied successfully\n");
                isResultZero = false;
            }
            else
            {
                transposeFXTM(currentResult, &transposed);
                DEBUG("else: transposed successfully\n");
                multiplyFXTM(&transposed, currentPower, nextResult);
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
            DEBUG("power = %llu", power)
            transposeFXTM(currentPower, &transposed);
            // perform a very fast, quasi-quadratic matrix multiplication
            multiplyFXTM(&transposed, currentPower, nextPower);

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