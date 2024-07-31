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

int testSolverWithFileInput()
{
    ParameterSet pSet;
    initParamsFromFile(&pSet, "input\\example.txt");
    Solution sol;
    solveForStartingState(&pSet, &sol);

    uint64_t maxParamValue = 1ULL << sol.knownVariableCount;
    Xoroshiro* possibleStates = (Xoroshiro*)malloc(maxParamValue * sizeof(Xoroshiro));
    getAllXoroshiroStates(&sol, possibleStates);

    for (int i = 0; i < maxParamValue; i++)
    {
        printf("Found possible Xoroshiro state: ( 0x%llx , 0x%llx )\n", possibleStates[i].lo, possibleStates[i].hi);
    }

    freeParamSet(&pSet);
}

int main()
{
    testSolverWithFileInput();
    // testFXTM();
    // testEquationSysSolver();
    // return 0;

    // return batchTestSolverCorrectness(100);

    // Xoroshiro state = { 0xdeadbeef, 0xdeadbeef };

    // for (int i = 0; i < 10; i++) {
    //     xNextLong(&state);
    //     printf("%d: %d %d\n", i+1, (state.lo >> i) & 1, 63-i);
    //     printf("%d: %d %d\n", i+1, (state.hi >> (63-i)) & 1, i+64);
    // }

    // state.lo = 0xdeadbeef;
    // state.hi = 0xdeadbeef;

    // for (int i = 0; i < 271; i++)
    //     xNextLong(&state);

    // printf("271: hi >> 32: %08llx\n", state.hi >> 32);
    // uint64_t v = state.hi >> 32;
    // for (int i = 95; i >= 64; i--, v >>= 1)
    //     printf("%d %d %d\n", 271, i, v&1);

    // state.lo = 0xdeadbeef;
    // state.hi = 0xdeadbeef;

    // for (int i = 0; i < 314; i++)
    //     xNextLong(&state);

    // printf("314: lo >> 32: %08llx\n", state.lo >> 32);
    // v = state.lo >> 32;
    // for (int i = 31; i >= 0; i--, v >>= 1)
    //     printf("%d %d %d\n", 314, i, v&1);

    // state.lo = 0xdeadbeef;
    // state.hi = 0xdeadbeef;

    // for (int i = 0; i < 27182; i++)
    //     xNextLong(&state);

    // printf("27182: bits 0-11: %x\n", (int)(state.hi >> (63 - 11)));
    // v = state.hi >> (63 - 11);
    // for (int i = 11+64; i >= 64; i--, v >>= 1)
    //     printf("%d %d %d\n", 27182, i, v&1);

    // state.lo = 0xdeadbeef;
    // state.hi = 0xdeadbeef;

    // for (int i = 0; i < 31415; i++)
    //     xNextLong(&state);

    // printf("31415: hi >> 32: %08x\n", (int)(state.hi >> 32));
    // v = state.hi >> 32;
    // for (int i = 95; i >= 64; i--, v >>= 1)
    //     printf("%d %d %d\n", 31415, i, v&1);
}