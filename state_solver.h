#pragma once
#include "gaussian_elim.h"
#include "xoroshiro.h"

typedef struct ParameterSet ParameterSet;
struct ParameterSet {
    int paramCount;

    int* stateIDs;
    int* bitIDs;
    int* bitValues;
};

void initParamSet(ParameterSet *obj, int size);
void freeParamSet(ParameterSet *obj);

void initParamsFromFile(ParameterSet *obj, const char* filename);

int solveForStartingState(const ParameterSet* paramSet, Solution* sol);

ParameterSet* genRandomParamSet(uint64_t seed, const Xoroshiro state);
int testSolverCorrectness(uint64_t seed);
int batchTestSolverCorrectness(int testCount);