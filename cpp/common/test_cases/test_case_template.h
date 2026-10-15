///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

enum class TestCase {
    KLEIN_5TO3,
    HYPERSHPERE_5TO1,
    CIRCLE_9TO8,
    HYPERSHPERE_N_TO_K,
    S1_K_TIMES,
    BOOK_SPHERE3,
    BOOK_SPHERE4,
    BOOK_TORUS4,
    BOOK_ZCOSW,
    BOOK_Z2COSW,
    BOOK_Z2W5,
    BOOK_S2XS1,
    BOOK_S2XS1_COARSE,
    BOOK_AMBIG_SPLIT,
    BOOK_AMBIG_TUBE,
    BOOK_J1A_ELLIPSE,
    BOOK_J1A_TORUS,
    BOOK_J1A_S3,
    BOOK_J1A_ZCOSW,
    MOD_TWO_SPHERES,
    MOD_TWO_SPHERES_P8,
    MOD_EIGHT_TORI,
    MOD_KNOT,
    MOD_INTERSECT3,
    MOD_R4_TWO_BALLS,
    MOD_R4_DIFF,
    MOD_R4_TUBE,
    MOD_PIPEBOX,
    MOD_PIPEBOX_P8,
    BOOK_HYPERBOLA4,
    BOOK_T3_R6_B,
    BOOK_CONV_S2,
    BOOK_CONV_T2,
    BOOK_BRATU3,
    BOOK_BRATU5,
    BOOK_BLOBS3,
    BOOK_BLOBS4,
};

template <TestCase CASE>
struct TestCaseImpl;
