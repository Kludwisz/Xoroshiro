#pragma once
#include "inttypes.h"
#include "stdbool.h"

typedef struct Equation Equation;
struct Equation {
    uint64_t lhs[2];
    int rhs;
};

typedef struct Solution Solution;
struct Solution {
    bool isContradictory;

    int knownVariableCount;
    int parameterCount;
    bool isParameter[128];

    int equationCount;
    Equation equations[128];
};

// ----------------------------------

void initEquation(Equation *eq, uint64_t lhsLo, uint64_t lhsHi, int rhs);
void xorEquationBy(Equation *a, const Equation *b);

int gaussianElimGF2(Equation equations[], int equationCount, Solution* sol);