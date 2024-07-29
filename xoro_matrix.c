#include "xoro_matrix.h"

void transformToDiagonal(const FXTMatrix* matrix, FXTDiagMatrix* result)
{
    static const uint64_t ONLY_TOP_BIT = 0x8000'0000'0000'0000ULL;
    int iDiag = 0;

    // create lower diagonals (without main)
    for (int i0 = 127; i0 > 0; i0--, iDiag++)
    {
        uint64_t mask = ONLY_TOP_BIT;
        int i = i0;
        int j = 0;
        int col = 0;
        while (i >= 127)
        {
            if (j >= 64)
            {
                mask = ONLY_TOP_BIT;
                j -= 64;
                col++;
            }

            (result->M)[iDiag][col] |= (mask & (matrix->M)[i][col])

            i++;
            j++;
            mask >>= 1;
        }
    }
}

void fastXoroMatrixMul(const FXTMatrix* a, const FXTDiagMatrix* b, FXTMatrix* product)
{

}

void xoroMatrixFastPower(const FXTMatrix* matrix, const uint64_t power, FXTMatrix* result)
{

}
