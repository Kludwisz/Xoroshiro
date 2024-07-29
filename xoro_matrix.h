#pragma once
#include "xoroshiro.h"
#include "inttypes.h"

#define ONLY_TOP_BIT (0x8000000000000000ULL)

// structures for the fast xoroshiro transform
// matrix multiplication procedure

typedef struct FXTMatrix FXTMatrix;
struct FXTMatrix {
    uint64_t M[2][2][64];
};

// typedef struct FXTDiagMatrix FXTDiagMatrix;
// struct FXTDiagMatrix {
//     uint64_t M[128 * 2 - 1][2];
// };

void copyFXTMatrix(const FXTMatrix *from, FXTMatrix *to);
//void transformToDiagonal(const FXTMatrix *matrix, FXTDiagMatrix *result);
//void fastXoroMatrixMul(const FXTMatrix *a, const FXTDiagMatrix *b, FXTMatrix *product);
void fastXoroMatrixPower(const FXTMatrix *matrix, uint64_t power, FXTMatrix *result);

void advanceXoroshiroFXTM(Xoroshiro *state, const FXTMatrix* fxtm);


// ------------------------------------------------------------------------------------

/*
static const FXTMatrix XOROSHIRO_STANDARD_MATRIX = {{
    { 9223653511831486464ULL, 134217728ULL },
    { 4611826755915743232ULL, 67108864ULL },
    { 2305913377957871616ULL, 33554432ULL },
    { 1152956688978935808ULL, 16777216ULL },
    { 576478344489467904ULL, 8388608ULL },
    { 288239172244733952ULL, 4194304ULL },
    { 144119586122366976ULL, 2097152ULL },
    { 72059793061183488ULL, 1048576ULL },
    { 36029896530591744ULL, 524288ULL },
    { 18014948265295872ULL, 262144ULL },
    { 9007474132647936ULL, 131072ULL },
    { 4503737066323968ULL, 65536ULL },
    { 2251868533161984ULL, 32768ULL },
    { 1125934266580992ULL, 16384ULL },
    { 562967133290496ULL, 8192ULL },
    { 281483566645248ULL, 4096ULL },
    { 140741783322624ULL, 2048ULL },
    { 70370891661312ULL, 1024ULL },
    { 35185445830656ULL, 512ULL },
    { 17592722915328ULL, 256ULL },
    { 8796361457664ULL, 128ULL },
    { 9223376435035504640ULL, 64ULL },
    { 4611688217517752320ULL, 32ULL },
    { 2305844108758876160ULL, 16ULL },
    { 1152922054379438080ULL, 8ULL },
    { 576461027189719040ULL, 4ULL },
    { 288230513594859520ULL, 2ULL },
    { 144115256797429760ULL, 1ULL },
    { 72057628398714880ULL, 9223372036854775808ULL },
    { 36028814199357440ULL, 4611686018427387904ULL },
    { 18014407099678720ULL, 2305843009213693952ULL },
    { 9007203549839360ULL, 1152921504606846976ULL },
    { 4503601774919680ULL, 576460752303423488ULL },
    { 2251800887459840ULL, 288230376151711744ULL },
    { 1125900443729920ULL, 144115188075855872ULL },
    { 562950221864960ULL, 72057594037927936ULL },
    { 281475110932480ULL, 36028797018963968ULL },
    { 140737555466240ULL, 18014398509481984ULL },
    { 70368777733120ULL, 9007199254740992ULL },
    { 35184388866560ULL, 4503599627370496ULL },
    { 17592194433280ULL, 2251799813685248ULL },
    { 8796097216640ULL, 1125899906842624ULL },
    { 4398048608320ULL, 562949953421312ULL },
    { 2199024304160ULL, 281474976710656ULL },
    { 1099512152080ULL, 140737488355328ULL },
    { 549756076040ULL, 70368744177664ULL },
    { 274878038020ULL, 35184372088832ULL },
    { 137439019010ULL, 17592186044416ULL },
    { 68719509505ULL, 8796093022208ULL },
    { 9223372071214530560ULL, 4398046511104ULL },
    { 4611686035607265280ULL, 2199023255552ULL },
    { 2305843017803632640ULL, 1099511627776ULL },
    { 1152921508901816320ULL, 549755813888ULL },
    { 576460754450908160ULL, 274877906944ULL },
    { 288230377225454080ULL, 137438953472ULL },
    { 144115188612727040ULL, 68719476736ULL },
    { 72057594306363520ULL, 34359738368ULL },
    { 36028797153181760ULL, 17179869184ULL },
    { 18014398576590880ULL, 8589934592ULL },
    { 9007199288295440ULL, 4294967296ULL },
    { 4503599644147720ULL, 2147483648ULL },
    { 2251799822073860ULL, 1073741824ULL },
    { 1125899911036930ULL, 536870912ULL },
    { 562949955518465ULL, 268435456ULL },
    { 9223372036854775808ULL, 134217728ULL },
    { 4611686018427387904ULL, 67108864ULL },
    { 2305843009213693952ULL, 33554432ULL },
    { 1152921504606846976ULL, 16777216ULL },
    { 576460752303423488ULL, 8388608ULL },
    { 288230376151711744ULL, 4194304ULL },
    { 144115188075855872ULL, 2097152ULL },
    { 72057594037927936ULL, 1048576ULL },
    { 36028797018963968ULL, 524288ULL },
    { 18014398509481984ULL, 262144ULL },
    { 9007199254740992ULL, 131072ULL },
    { 4503599627370496ULL, 65536ULL },
    { 2251799813685248ULL, 32768ULL },
    { 1125899906842624ULL, 16384ULL },
    { 562949953421312ULL, 8192ULL },
    { 281474976710656ULL, 4096ULL },
    { 140737488355328ULL, 2048ULL },
    { 70368744177664ULL, 1024ULL },
    { 35184372088832ULL, 512ULL },
    { 17592186044416ULL, 256ULL },
    { 8796093022208ULL, 128ULL },
    { 9223376434901286912ULL, 64ULL },
    { 4611688217450643456ULL, 32ULL },
    { 2305844108725321728ULL, 16ULL },
    { 1152922054362660864ULL, 8ULL },
    { 576461027181330432ULL, 4ULL },
    { 288230513590665216ULL, 2ULL },
    { 144115256795332608ULL, 1ULL },
    { 72057628397666304ULL, 9223372036854775808ULL },
    { 36028814198833152ULL, 4611686018427387904ULL },
    { 18014407099416576ULL, 2305843009213693952ULL },
    { 9007203549708288ULL, 1152921504606846976ULL },
    { 4503601774854144ULL, 576460752303423488ULL },
    { 2251800887427072ULL, 288230376151711744ULL },
    { 1125900443713536ULL, 144115188075855872ULL },
    { 562950221856768ULL, 72057594037927936ULL },
    { 281475110928384ULL, 36028797018963968ULL },
    { 140737555464192ULL, 18014398509481984ULL },
    { 70368777732096ULL, 9007199254740992ULL },
    { 35184388866048ULL, 4503599627370496ULL },
    { 17592194433024ULL, 2251799813685248ULL },
    { 8796097216512ULL, 1125899906842624ULL },
    { 4398048608256ULL, 562949953421312ULL },
    { 2199024304128ULL, 281474976710656ULL },
    { 1099512152064ULL, 140737488355328ULL },
    { 549756076032ULL, 70368744177664ULL },
    { 274878038016ULL, 35184372088832ULL },
    { 137439019008ULL, 17592186044416ULL },
    { 68719509504ULL, 8796093022208ULL },
    { 34359754752ULL, 4398046511104ULL },
    { 17179877376ULL, 2199023255552ULL },
    { 8589938688ULL, 1099511627776ULL },
    { 4294969344ULL, 549755813888ULL },
    { 2147484672ULL, 274877906944ULL },
    { 1073742336ULL, 137438953472ULL },
    { 536871168ULL, 68719476736ULL },
    { 268435584ULL, 34359738368ULL },
    { 134217792ULL, 17179869184ULL },
    { 67108896ULL, 8589934592ULL },
    { 33554448ULL, 4294967296ULL },
    { 16777224ULL, 2147483648ULL },
    { 8388612ULL, 1073741824ULL },
    { 4194306ULL, 536870912ULL },
    { 2097153ULL, 268435456ULL }
}};
/**/