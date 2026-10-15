// Marching Simplex on the adaptive triangulation J1a (Chap. 4), in the framework of the K1/J1 codes:
// same test cases (DOMAIN_MIN, DOMAIN_MAX, DOMAIN_DIV, func), same output format (Appendix A).
//
// usage: ms_j1a <levels> [out.pol|-] [check]
//   The DOMAIN_DIV level-0 blocks are refined <levels> times around the submanifold N = M ∩ G^{-1}(0)
//   when the test case defines `sub` (with SUBDIM components), or around M itself otherwise.
//   With "check", the program verifies that the simplices form a triangulation (facets shared by two).
#include <chrono>
#include <iostream>
#include <map>
#include <type_traits>
#include "../common/definitions.h"
#include "../common/utility.h"
#include "j1a.h"

using TC = TestCaseImpl<TestCase::SELECTED_CASE>;

template <typename T, typename = void> struct HasSub : std::false_type {};
template <typename T> struct HasSub<T, std::void_t<decltype(T::SUBDIM)>> : std::true_type {};

template <typename T> static void subMap(const double *x, double *f) {
    double fm[KDIM]; func(x, fm);
    for (size_t i = 0; i < KDIM; ++i) f[i] = fm[i];
    T::sub(x, f + KDIM);
}
template <typename T> static int subDim() { if constexpr (HasSub<T>::value) return (int)T::SUBDIM; else return 0; }
template <typename T> static j1a::Map refinementMap() {
    if constexpr (HasSub<T>::value) return subMap<T>;
    else return [](const double *x, double *f) { func(x, f); };
}

int main(int argc, char **argv) {
    const int levels = argc > 1 ? std::atoi(argv[1]) : 2;
    const std::string out = (argc > 2 && std::string(argv[2]) != "-") ? argv[2] : "";
    const bool check = argc > 3 && std::string(argv[3]) == "check";
    const int n = NDIM, k = KDIM;
    std::vector<int> N(DOMAIN_DIV.begin(), DOMAIN_DIV.end());
    j1a::Grid G(n, N, j1a::Vec(DOMAIN_MIN.begin(), DOMAIN_MIN.end()), j1a::Vec(DOMAIN_MAX.begin(), DOMAIN_MAX.end()), levels);

    auto t0 = std::chrono::high_resolution_clock::now();
    // ---- adaptive refinement
    const int kr = k + subDim<TC>(); j1a::Map R = refinementMap<TC>();
    j1a::Approximator AR(G, kr, R, PERTURBATION_AMPLITUDE);
    std::vector<j1a::LocalSimplex> L;
    for (int lev = 0; lev < levels; ++lev) {
        std::vector<j1a::Block> mark;
        for (auto &b : G.blocks()) {
            if (b.level != lev) continue;
            j1a::enumerateSimplices(n, G.refinedFaces(b), L);
            j1a::Cell C; bool hit = false;
            for (auto &S : L) {
                std::vector<const std::vector<double> *> V;
                for (auto &w : S.v) V.push_back(&AR.vertex(j1a::latticeCoords(G, b, w)));
                if (j1a::buildCell(n, kr, V, C, false)) { hit = true; break; }
            }
            if (hit) mark.push_back(b);
        }
        for (auto &b : mark) G.refine(b);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    // ---- approximation of M on the final triangulation
    j1a::Approximator A(G, k, [](const double *x, double *f) { func(x, f); }, PERTURBATION_AMPLITUDE);
    FILE *fp = nullptr;
    if (!out.empty()) {
        fp = std::fopen(out.c_str(), "w");
        std::fprintf(fp, "%d %d\n", n, k);
        for (int i = 0; i < n; ++i) std::fprintf(fp, "%d ", N[i]);
        std::fprintf(fp, "\n\n");
    }
    std::map<int, uint64_t> perLevel, crossedPerLevel;
    uint64_t nSimp = 0, nCross = 0, nTransition = 0;
    std::unordered_map<uint64_t, int> facets;                       // for the conformity check
    double volume = 0;
    j1a::Cell C;
    for (auto &b : G.blocks()) {
        auto Phi = G.refinedFaces(b);
        if (!Phi.empty()) ++nTransition;
        j1a::enumerateSimplices(n, Phi, L);
        perLevel[b.level] += L.size(); nSimp += L.size();
        for (size_t q = 0; q < L.size(); ++q) {
            auto &S = L[q];
            std::vector<std::vector<int64_t>> Z;
            std::vector<const std::vector<double> *> V;
            for (auto &w : S.v) { Z.push_back(j1a::latticeCoords(G, b, w)); V.push_back(&A.vertex(Z.back())); }
            if (check) {
                std::vector<uint64_t> keys; for (auto &z : Z) keys.push_back(j1a::latticeKey(G, z));
                for (int i = 0; i <= n; ++i) {
                    std::vector<uint64_t> f; for (int j = 0; j <= n; ++j) if (j != i) f.push_back(keys[j]);
                    std::sort(f.begin(), f.end());
                    uint64_t hsh = 1469598103934665603ULL; for (auto x : f) hsh = (hsh ^ x) * 1099511628211ULL;
                    facets[hsh]++;
                }
                // volume in lattice units
                std::vector<double> M(n * n);
                for (int r = 1; r <= n; ++r) for (int c = 0; c < n; ++c) M[(r - 1) * n + c] = double(Z[r][c] - Z[0][c]);
                double det = 1;
                for (int c = 0; c < n; ++c) {
                    int p = c; for (int r = c + 1; r < n; ++r) if (std::fabs(M[r * n + c]) > std::fabs(M[p * n + c])) p = r;
                    if (M[p * n + c] == 0) { det = 0; break; }
                    if (p != c) { for (int q2 = 0; q2 < n; ++q2) std::swap(M[c * n + q2], M[p * n + q2]); det = -det; }
                    det *= M[c * n + c];
                    for (int r = c + 1; r < n; ++r) { double f = M[r * n + c] / M[c * n + c]; for (int q2 = c; q2 < n; ++q2) M[r * n + q2] -= f * M[c * n + q2]; }
                }
                volume += std::fabs(det);
            }
            if (!j1a::buildCell(n, k, V, C, fp != nullptr)) continue;
            ++nCross; crossedPerLevel[b.level]++;
            if (fp) {
                std::fprintf(fp, "%llu", (unsigned long long)(G.key(b) * 1000000ULL + q));
                for (int i = 0; i < n; ++i) std::fprintf(fp, " %lld", (long long)b.g[i]);
                std::fprintf(fp, "\n1\n%zu\n", C.vmask.size());
                for (size_t a = 0; a < C.vmask.size(); ++a) {
                    std::fprintf(fp, "%d ", b.level);
                    for (int i = 0; i < n; ++i) std::fprintf(fp, " %.12g", C.pts[a][i]);
                    std::fprintf(fp, "\n");
                }
                for (int j = 1; j <= n - k; ++j) {
                    std::fprintf(fp, "%zu\n", C.cmask[j].size());
                    for (auto &bb : C.bnd[j]) { for (int x : bb) std::fprintf(fp, "%d ", x + 1); std::fprintf(fp, "\n"); }
                }
            }
        }
    }
    if (fp) { std::fprintf(fp, "-1\n"); std::fclose(fp); }
    auto t2 = std::chrono::high_resolution_clock::now();

    std::map<int, int> blocksPerLevel; for (auto &b : G.blocks()) blocksPerLevel[b.level]++;
    std::cout << "n=" << n << " k=" << k << " levels=" << levels << (HasSub<TC>::value ? " (refined around N)" : " (refined around M)") << "\n";
    std::cout << "blocks per level:"; for (auto &p : blocksPerLevel) std::cout << " " << p.first << ":" << p.second; std::cout << "\n";
    std::cout << "blocks: " << G.leaves.size() << "  transition blocks: " << nTransition << "\n";
    std::cout << "simplices: " << nSimp << "  crossed: " << nCross << "\n";
    std::cout << "crossed per level:"; for (auto &p : crossedPerLevel) std::cout << " " << p.first << ":" << p.second; std::cout << "\n";
    std::cout << "F evaluations: " << A.nEval << " (refinement: " << AR.nEval << ")\n";
    std::cout << "time refinement: " << std::chrono::duration<double>(t1 - t0).count() << " s, approximation: "
              << std::chrono::duration<double>(t2 - t1).count() << " s\n";
    uint64_t finest = 1; for (int i = 0; i < n; ++i) finest *= (uint64_t)N[i] << levels;
    uint64_t fact = 1; for (int i = 2; i <= n; ++i) fact *= i;
    std::cout << "uniform grid at the finest level: " << finest << " blocks, " << finest * (fact << n) << " simplices\n";
    if (check) {
        std::map<int, uint64_t> cnt; for (auto &p : facets) cnt[p.second]++;
        std::cout << "facet multiplicities:"; for (auto &p : cnt) std::cout << " " << p.first << "x" << p.second; std::cout << "\n";
        double unit = 1.0; uint64_t f = 1; for (int i = 2; i <= n; ++i) f *= i;
        double boxv = 1; for (int i = 0; i < n; ++i) boxv *= double(N[i]) * 4.0 * double(1 << levels);
        std::cout << "volume of simplices / volume of box: " << volume / double(f) / boxv << "\n"; (void)unit;
    }
    return 0;
}
