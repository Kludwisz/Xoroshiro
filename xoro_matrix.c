#include "xoro_matrix.h"
#include "utils.h"
static const bool DEBUG_MODE = false;

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

// bT is the FXTM "transposition" of a
void multiplyFXTM(const FXTMatrix* a, const FXTMatrix* bT, FXTMatrix* c)
{
    // calculate each bit in the result matrix separately
    for (int qi = 0; qi < 2; qi++) for (int qj = 0; qj < 2; qj++)
    {
        for (int i = 0; i < 64; i++)
        {
            uint64_t res = 0ULL;
            // these will be constant throughout the j-loop
            const uint64_t aLo = a->M[qi][0][i], aHi = a->M[qi][1][i];

            for (int j = 0; j < 64; j++)
            {
                res <<= 1;
                const uint64_t mulXor = (aLo & bT->M[0][qj][j]) ^ (aHi & bT->M[1][qj][j]);
                const int newBit = __builtin_parityll(mulXor); // aka: GCC, do your thing!
                res |= newBit; 
            }

            (c->M)[qi][qj][i] = res;
            // DEBUG("%llx\n", res);
        }
    }
}

// void xoroMatrixPowerHelper(
//     uint64_t partialPower, bool shouldFillPowers, bool* isResultZero, 
//     FXTMatrix **currentPower, FXTMatrix **nextPower, FXTMatrix **currentResult, FXTMatrix **nextResult
// )
// {
//     FXTMatrix transposed = { 0 };
//     FXTMatrix *temp;

//     int i = 1;
//     while (partialPower != 0ULL || (shouldFillPowers && i <= 64))
//     {
//         bool transpositionDone = false;

//         if (partialPower & 1ULL)
//         {
//             DEBUG("if power & 1:  power = %llu\n", partialPower);
//             if (*isResultZero)
//             {
//                 copyFXTM(*currentPower, *nextResult);
//                 DEBUG("copied successfully\n");
//                 *isResultZero = false;
//             }
//             else
//             {
//                 transposeFXTM(*currentPower, &transposed);
//                 transpositionDone = true;
//                 DEBUG("else: transposed successfully\n");
//                 multiplyFXTM(*currentResult, &transposed, *nextResult);
//                 DEBUG("else: performed fast mul successfully\n");
//             }
                
//             // swap next and current result, clearing space for the next operations
//             temp = *currentResult;
//             *currentResult = *nextResult;
//             *nextResult = temp;
//         }

//         // only calculate next power when needed
//         partialPower >>= 1;
//         if (partialPower != 0ULL || shouldFillPowers)
//         {
//             // create the transposition of currentPower
//             DEBUG("calc next:  power = %llu\n", partialPower)
//             if (!transpositionDone)
//                 transposeFXTM(*currentPower, &transposed);

//             // perform a very fast, quasi-quadratic matrix multiplication
//             multiplyFXTM(*currentPower, &transposed, *nextPower);

//             // swap next and current power, making space for the next operations
//             temp = *currentPower;
//             *currentPower = *nextPower;
//             *nextPower = temp;
//         }

//         printf("%d\n", i);
//         i++;
//     }
// }

void fastXoroMatrixPower(const FXTMatrix *matrix, uint64_t power, FXTMatrix *result)
{
    fastXoroMatrixPower128(matrix, 0ULL, power, result);
}

// the two uint64_t fields give us the ability to advance by an arbitrary number
// of Xoroshiro states. TODO: consider adding a struct for the uint64_t's
void fastXoroMatrixPower128(const FXTMatrix* matrix, uint64_t powerUpper64, uint64_t powerLower64,  FXTMatrix* result)
{
    FXTMatrix res1 = { 0 }, res2 = { 0 };
    FXTMatrix pow1 = { 0 }, pow2 = { 0 };
    FXTMatrix transposed = { 0 };
    FXTMatrix *temp;
    
    bool isResultZero = true;

    FXTMatrix *currentPower = &pow1, *nextPower = &pow2, *currentResult = &res1, *nextResult = &res2;
    copyFXTM(matrix, currentPower);

    uint64_t power = powerLower64;
    for (int bitID = 0; bitID < 128; bitID++)
    {
        if (power == 0ULL)
        {
            if (bitID > 63)
                break;  // already did upper powers
            else if (powerUpper64 == 0ULL)
                break;  // there is no upper power, just break
            // else continue until condition 3. reached, filling all the lower powers
        }

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
                multiplyFXTM(currentResult, &transposed, nextResult);
                DEBUG("else: performed fast mul successfully\n");
            }
                
            // swap next and current result, clearing space for the next operations
            temp = currentResult;
            currentResult = nextResult;
            nextResult = temp;
        }

        // create the transposition of currentPower
        DEBUG("calc next:  power = %llu\n", power)
        if (!transpositionDone)
            transposeFXTM(currentPower, &transposed);

        // perform a very fast, quasi-quadratic matrix multiplication
        multiplyFXTM(currentPower, &transposed, nextPower);
        power >>= 1;
        DEBUG("Multiplied, current bitID = %d\n", bitID)

        // swap next and current power, making space for the next operations
        temp = currentPower;
        currentPower = nextPower;
        nextPower = temp;

        if (bitID == 63)
        {
            // at key point, switch out bits to upper power
            DEBUG("Switch to upper power\n")
            power = powerUpper64;
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