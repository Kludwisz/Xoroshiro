#include "state_solver.h"
#include "gaussian_elim.h"
#include "xoroshiro.h"
#include "xoro_matrix.h"
#include "utils.h"
const bool DEBUG_MODE = true;

#include "time.h"


void initParamSet(ParameterSet *obj, int size)
{
    obj->paramCount = size;
    obj->stateIDs = (int*)malloc(size * sizeof(int));
    obj->bitIDs = (int*)malloc(size * sizeof(int));
    obj->bitValues = (int*)malloc(size * sizeof(int));
}

void freeParamSet(ParameterSet *obj)
{
    free(obj->stateIDs);
    free(obj->bitIDs);
    free(obj->bitValues);
}

// ------------------------------------------------------------

ParameterSet* genRandomParamSet(uint64_t seed)
{
    ParameterSet ps;

}

int testSolverCorrectness(uint64_t seed)
{
    srand((int)seed);
    return 0;
}

int batchTestSolverCorrectness(int testCount)
{
    uint64_t startSeed = (uint64_t)time(NULL);

    for (int i = 0; i < testCount; i++)
    {
        if (testSolverCorrectness(startSeed + i) != 0)
        {
            DEBUG("Failed for seed %llu\n", startSeed + i)
            return 1;
        }
        DEBUG("Test %d OK!\n", i)
    }
    return 0;
}
