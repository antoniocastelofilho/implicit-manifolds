// Marching Simplex: visits every cube of the grid and every simplex of its K1 triangulation.
#include <chrono>
#include <iostream>
#include "simplicial.h"

int main(int argc, char **argv) {
    std::string out = (argc > 1) ? argv[1] : "";
    PolWriter W(out);
    const bool build = !out.empty();
    auto t0 = std::chrono::high_resolution_clock::now();
    uint64_t ncubes = 1; for (size_t i = 0; i < NDIM; ++i) ncubes *= DOMAIN_DIV[i];
    uint64_t n_eval = 0, n_simplices = 0, n_crossed = 0;
    std::vector<std::array<double, NDIM>> cx(1ULL << NDIM);
    std::vector<std::array<double, KDIM>> cf(1ULL << NDIM);
    Simplex s; SimplexData d; Polytope P;
    for (uint64_t g = 0; g < ncubes; ++g) {
        uint64_t c = g;
        for (int64_t i = NDIM - 1; i >= 0; --i) { s.y[i] = (int64_t)(c % DOMAIN_DIV[i]); c /= DOMAIN_DIV[i]; }
        // F at the 2^n vertices of the cube
        for (size_t b = 0; b < (1ULL << NDIM); ++b) {
            vertexCoords(s.y, b, cx[b]);
            func(cx[b].data(), cf[b].data());
            ++n_eval;
        }
        std::array<uint8_t, NDIM> perm;
        for (size_t i = 0; i < NDIM; ++i) perm[i] = (uint8_t)i;
        do {
            s.p = perm;
            for (size_t j = 0; j <= NDIM; ++j) {
                const size_t b = localLabel(s, j);
                d.lab[j] = b; d.x[j] = cx[b]; d.f[j] = cf[b];
            }
            ++n_simplices;
            if (processSimplex(d, P, build)) { ++n_crossed; W.write(s, d, P); }
        } while (std::next_permutation(perm.begin(), perm.end()));
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time (MS): " << std::chrono::duration<double>(t1 - t0).count() << " s" << std::endl;
    std::cout << "cubes: " << ncubes << "  simplices: " << n_simplices << "  crossed: " << n_crossed << std::endl;
    std::cout << "F evaluations: " << n_eval << std::endl;
    return 0;
}
