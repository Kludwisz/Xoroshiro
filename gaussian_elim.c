#include "gaussian_elim.h"
#include "utils.h"
const bool DEBUG_MODE = false;


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

int gaussianElimGF2(Equation equations[], int equationCount)
{
    bool* usedEqns = (bool*)malloc(sizeof(bool) * equationCount);
    for (int i = 0; i < equationCount; i++)
        usedEqns[i] = false;

    for (int seg = 0; seg <= 1; seg++)
    {
        uint64_t mask = 0x8000'0000'0000'0000ULL;
        for (int i = 0; i < 64; i++, mask >>= 1)
        {
            for (int eID = 0; eID < equationCount; eID++)
            {
                if ((equations[eID].lhs[seg] & mask) == 0ULL || usedEqns[eID])
                    continue;
                usedEqns[eID] = true;

                for (int eID2 = 0; eID2 < equationCount; eID2++)
                {
                    if (eID == eID2)
                        continue;
                    xorEquationBy(&(equations[eID2]), &(equations[eID]));
                }
                break;
            }
        }
    }

    free(usedEqns);
    return 0;
}
