#include "xoro_matrix.h"
#include "xoroshiro.h"
#include "gaussian_elim.h"
#include "state_solver.h"
#include "utils.h"
static const bool DEBUG_MODE = false;

// -----------------------------------------------------------------------------------

void printXoro(const Xoroshiro* state)
{
    printf("Current Xoroshiro state:   ( %llu, %llu )\n", state->lo, state->hi);
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

int testEquationSysSolver() {
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

int main()
{
    // testFXTM();
    // testEquationSysSolver();
    // return 0;

    return batchTestSolverCorrectness(100);
}