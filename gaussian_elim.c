#include "gaussian_elim.h"
#include "utils.h"
static const bool DEBUG_MODE = false;


void initEquation(Equation *eq, uint64_t lhsLo, uint64_t lhsHi, int rhs)
{
    eq->lhs[0] = lhsLo;
    eq->lhs[1] = lhsHi;
    eq->rhs = rhs;
}

void xorEquationBy(Equation *a, const Equation *b)
{
    a->lhs[0] ^=  b->lhs[0];
    a->lhs[1] ^=  b->lhs[1];
    a->rhs ^= b->rhs;
}

int eqnCompare(const void* a_void, const void* b_void)
{
    const Equation *a = (Equation*)a_void, *b = (Equation*)b_void;

    if (a->lhs[0] > b->lhs[0])
        return 1;
    else if (a->lhs[0] < b->lhs[0])
        return -1;
    
    if (a->lhs[1] > b->lhs[1])
        return 1;
    else if (a->lhs[1] < b->lhs[1])
        return -1;
    
    return 0;
}

int gaussianElimGF2(Equation equations[], int equationCount, Solution* sol)
{
    bool* usedEqns = (bool*)malloc(sizeof(bool) * equationCount);
    for (int i = 0; i < equationCount; i++)
        usedEqns[i] = false;

    sol->parameterCount = 128;
    sol->knownVariableCount = 0;
    for (int i = 0; i < 128; i++)
        sol->isParameter[i] = true;

    for (int seg = 0; seg <= 1; seg++)
    {
        uint64_t mask = ONLY_TOP_BIT;
        for (int i = 0; i < 64; i++, mask >>= 1)
        {
            for (int eID = 0; eID < equationCount; eID++)
            {
                if (usedEqns[eID] || (equations[eID].lhs[seg] & mask) == 0ULL)
                    continue;
                usedEqns[eID] = true;
                sol->isParameter[i + seg*64] = false;
                (sol->parameterCount)--;
                (sol->knownVariableCount)++;

                for (int eID2 = 0; eID2 < equationCount; eID2++)
                {
                    if (eID == eID2 || (equations[eID2].lhs[seg] & mask) == 0ULL)
                        continue;
                    xorEquationBy(&(equations[eID2]), &(equations[eID]));
                }
                break;
            }
        }
    }
    free(usedEqns);

    // sort the equations in ascending order for readability
    qsort(equations, equationCount, sizeof(equations[0]), eqnCompare);

    DEBUG("Sorted!\n");
    for (int i = 0; i < equationCount; i++) {
        DEBUG("Equation %d:  %016llx %016llx = %d\n", i+1, equations[i].lhs[0], equations[i].lhs[1], equations[i].rhs);
    }

    // check if there are any redundant equations (0 = 0) and/or falsehoods (0 = 1)
    int firstNonRedundantEq = 0;

    for (int eID = 0; eID < equationCount; eID++)
    {
        if (equations[eID].lhs[0] != 0ULL || equations[eID].lhs[1] != 0ULL)
            break;
        firstNonRedundantEq++;

        if (equations[eID].rhs == 1)
        {
            sol->isContradictory = true;
            return -1;
        }
    }

    // fill the Solution struct
    sol->isContradictory = false;
    sol->equationCount = equationCount - firstNonRedundantEq;

    for (int eID = 0; eID < sol->equationCount; eID++)
        sol->equations[eID] = equations[eID + firstNonRedundantEq];
    
    return 0;
}
