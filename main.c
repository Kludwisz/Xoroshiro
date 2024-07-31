#include "bit_matrix.h"
#include "xoro_matrix.h"
#include "xoroshiro.h"
#include "gaussian_elim.h"
#include "utils.h"
static const bool DEBUG_MODE = false;


#define CHECKED_OPERATION(code) \
{ \
    ErrType e = code; \
    if (e != SUCCESS) { \
        printf("Line %d: %s\n", __LINE__, getErrorMessage(e)); \
        goto ERROR; \
    } \
}

// ------------------------------------------------------------

ErrType createXoroshiroTransformationMatrix(BitMatrix* mx) 
{
    if (mx->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (mx->sizeI != 128 || mx->sizeJ != 128)
        return BAD_MATRIX_SIZE;

    makeZeroMatrix(mx);

    BitMatrix I64_const;
    allocMatrix(&I64_const, 64, 64);
    CHECKED_OPERATION( makeUnitMatrix(&I64_const) );

    BitMatrix Rot64;
    BitMatrix Shift64;
    BitMatrix ShiftPlusUnit64;
    BitMatrix TopLeftTransform64;
    allocMatrix(&Rot64, 64, 64);
    allocMatrix(&Shift64, 64, 64);
    allocMatrix(&ShiftPlusUnit64, 64, 64);
    allocMatrix(&TopLeftTransform64, 64, 64);
    
    // create all submatrices and inject into transformation matrix
    const int a = 49;
    const int b = 21;
    const int c = 28;

    CHECKED_OPERATION( makeLeftRotationMatrix(&Rot64, a) );
    CHECKED_OPERATION( makeLeftShiftMatrix(&Shift64, b) );

    CHECKED_OPERATION( matrixSum(&Shift64, &I64_const, &ShiftPlusUnit64) );
    CHECKED_OPERATION( matrixInject(&ShiftPlusUnit64, mx, 64, 0) ); // bottom left

    CHECKED_OPERATION( matrixSum(&Rot64, &ShiftPlusUnit64, &TopLeftTransform64) );
    CHECKED_OPERATION( matrixInject(&TopLeftTransform64, mx, 0, 0) ); // top left

    CHECKED_OPERATION( makeLeftRotationMatrix(&Rot64, c) );
    CHECKED_OPERATION( matrixInject(&Rot64, mx, 64, 64) ); // bottom right
    CHECKED_OPERATION( matrixInject(&Rot64, mx, 0, 64) );  // top right

    // -----------------------
    deallocMatrix(&I64_const);
    deallocMatrix(&Rot64);
    deallocMatrix(&Shift64);
    deallocMatrix(&ShiftPlusUnit64);
    deallocMatrix(&TopLeftTransform64);
    return SUCCESS;

    ERROR:
    deallocMatrix(&I64_const);
    deallocMatrix(&Rot64);
    deallocMatrix(&Shift64);
    deallocMatrix(&ShiftPlusUnit64);
    deallocMatrix(&TopLeftTransform64);
    return UNEXPECTED_ERROR;
}

ErrType createXoroshiroState(BitMatrix* stateVector, const Xoroshiro* const xr)
{
    if (stateVector->M == NULL)
        return MATRIX_UNINITIALIZED;
    if (stateVector->sizeI != 1 || stateVector->sizeJ != 128)
        return BAD_MATRIX_SIZE;

    makeZeroMatrix(stateVector);

    int i = 127;
    for (uint64_t hi = xr->hi; hi != 0ULL; hi >>= 1)
    {
        int bit = (int)(hi & 1ULL);
        (stateVector->M)[0][i] = bit;
        i--;
    }

    i = 63;
    for (uint64_t lo = xr->lo; lo != 0ULL; lo >>= 1)
    {
        int bit = (int)(lo & 1ULL);
        (stateVector->M)[0][i] = bit;
        i--;
    }

    return SUCCESS;
}

void printStateAsInt(const BitMatrix* const state)
{
    uint64_t lo = 0, hi = 0;

    for (int i = 0; i < 64; i++) {
        lo <<= 1;
        lo |= (state->M)[0][i];
    }
    for (int i = 64; i < 128; i++) {
        hi <<= 1;
        hi |= (state->M)[0][i];
    }

    printf("State after matrix multiplication(s): ( %llu , %llu )\n", lo, hi);
}

void printXoro(const Xoroshiro* const xr)
{
    printf("Current Xoroshiro state:    ( %llu , %llu )\n", xr->lo, xr->hi);
}

// -----------------------------------------------------------------------------------

ErrType testMatrixAdvancement(const int advanceCount)
{
    BitMatrix xoroshiroTransform;
    BitMatrix stateVector;
    BitMatrix newStateVector;
    allocMatrix(&xoroshiroTransform, 128, 128);
    allocMatrix(&stateVector, 1, 128);
    allocMatrix(&newStateVector, 1, 128);

    CHECKED_OPERATION( createXoroshiroTransformationMatrix(&xoroshiroTransform) );

    Xoroshiro xrsr = { 6128965821029732131ULL, 1902848264284609573ULL };
    CHECKED_OPERATION( createXoroshiroState(&stateVector, &xrsr) );
    
    BitMatrix *a, *product, *temp;
    a = &stateVector;
    product = &newStateVector;

    for (int i = 0; i < advanceCount; i++)
    {
        // advance the standard way
        xNextLong(&xrsr);

        // advance the matrix multiplication way
        CHECKED_OPERATION( matrixTimesMatrix(a, &xoroshiroTransform, product) );
        temp = a;
        a = product;
        product = temp;

        // check the states
        printf("Advanced by %d:\n", i+1);
        printXoro(&xrsr);
        printStateAsInt(a);
        printf("---------------------------------\n");
    }

    // ---------------------------------
    deallocMatrix(&xoroshiroTransform);
    deallocMatrix(&stateVector);
    deallocMatrix(&newStateVector);
    return SUCCESS;

    ERROR:
    deallocMatrix(&xoroshiroTransform);
    deallocMatrix(&stateVector);
    deallocMatrix(&newStateVector);
    return UNEXPECTED_ERROR;
}

// -----------------------------------------------------------------------------------

void getStandardFXTM(FXTMatrix *fxtm)
{
    BitMatrix mx;
    allocMatrix(&mx, 128, 128);
    (void)createXoroshiroTransformationMatrix(&mx);

    // copy bit data to FXTM
    for (int i = 0; i < 128; i++)
    {
        uint64_t low = 0ULL, high = 0ULL;

        for (int j = 0; j < 64; j++)
        {
            low <<= 1;
            low |= (mx.M)[i][j];
        }
        for (int j = 64; j < 128; j++)
        {
            high <<= 1;
            high |= (mx.M)[i][j];
        }

        (fxtm->M)[i / 64][0][i & 63] = low;
        (fxtm->M)[i / 64][1][i & 63] = high;
    }

    deallocMatrix(&mx);
}

void printStandardFXTM()
{
    FXTMatrix stdXrsr = { 0 };
    getStandardFXTM(&stdXrsr);

    for (int qi = 0; qi < 2; qi++) 
    {
        printf("{\n");
        for (int qj = 0; qj < 2; qj++)
        {
            printf("{\n");

            for (int i = 0; i < 64; i++)
            {
                printf("%lluULL", stdXrsr.M[qi][qj][i]);
                printf(i != 63 ? ",\n" : "\n");
            }

            printf(qj == 1 ? "}\n" : "},\n");
        }
        printf(qi == 1 ? "};\n" : "},\n");
    }
}

int testFXTM() 
{
    Xoroshiro x = { 145982789338521ULL, 70932749124324ULL };
    Xoroshiro x2 = { 145982789338521ULL, 70932749124324ULL };

    FXTMatrix fxtm = { 0 };

    fastXoroMatrixPower128(&XOROSHIRO_STANDARD_FXTM, ONLY_TOP_BIT, 0ULL, &fxtm);
    advanceXoroshiroFXTM(&x2, &fxtm); // 2^127
    advanceXoroshiroFXTM(&x2, &fxtm); // 2^127
    // advanced by 2^128, should be equivalent to 1
    // as the period is 2^128 - 1
    printXoro(&x2);
    
    xNextLong(&x);
    printXoro(&x);

    Xoroshiro x3 = { 123456789ULL, 987654321ULL };
    fastXoroMatrixPower128(&XOROSHIRO_STANDARD_FXTM, FULL_64, FULL_64, &fxtm);

    printXoro(&x3);
    advanceXoroshiroFXTM(&x3, &fxtm); // 2^128 - 1
    printXoro(&x3);

    return 0;
}

int main()
{
    Equation eqnsTest[4];
    initEquation(eqnsTest, 0xfULL, 11ULL, 1);
    initEquation(eqnsTest+1, 0x7ULL, 22ULL, 0);
    initEquation(eqnsTest+2, 0xaULL, 15ULL, 1);
    initEquation(eqnsTest+3, 0x3ULL, 9ULL, 0);

    Solution result;
    gaussianElimGF2(eqnsTest, 4, &result);

    if (result.isContradictory)
    {
        printf("No solutions exist.\n");
        return 0;
    }

    printf("Reduced equation count: %d\n", result.equationCount);
    printf("Known variables: %d, parameters: %d\n", result.knownVariableCount, result.parameterCount);

    return 0;
}