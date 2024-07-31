#pragma once
#include "inttypes.h"

typedef struct Equation Equation;
struct Equation {
    uint64_t lhs[2];
    int rhs;
};

// ----------------------------------

void initEquation(Equation *eq, uint64_t lhsLo, uint64_t lhsHi, int rhs);
void xorEquationBy(Equation *a, const Equation *b);

int gaussianElimGF2(Equation equations[], int equationCount);