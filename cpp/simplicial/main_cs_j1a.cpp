// Continuation Simplex on the adaptive triangulation J1a (Chaps. 4 and 6).
// The refinement itself is driven by continuation: at each level the component of M through
// FIRST_POINT is traced, and the blocks of the traced simplices that meet N = M ∩ G^{-1}(0)
// (or M itself, if the test case has no `sub`) are refined.  The neighbour of a simplex across a
// facet is found by locating the block that contains a point just beyond the facet and searching
// its simplices for the one that shares the facet.
//
// usage: cs_j1a <levels> [out.pol|-]
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

struct BlockSimplices {                       // simplices of a block, with lattice keys of their vertices
    std::vector<j1a::LocalSimplex> L;
    std::vector<std::vector<uint64_t>> keys;
    std::vector<std::vector<std::vector<int64_t>>> Z;
};

class Tracer {
public:
    const j1a::Grid &G; j1a::Approximator &A; int n, k;
    std::unordered_map<uint64_t, BlockSimplices> cache;
    uint64_t lookups = 0;
    Tracer(const j1a::Grid &G_, j1a::Approximator &A_) : G(G_), A(A_), n(G_.n), k(A_.k) {}
    BlockSimplices &simplices(const j1a::Block &b) {
        uint64_t bk = G.key(b);
        auto it = cache.find(bk);
        if (it != cache.end()) return it->second;
        BlockSimplices S;
        j1a::enumerateSimplices(n, G.refinedFaces(b), S.L);
        for (auto &s : S.L) {
            std::vector<uint64_t> ks; std::vector<std::vector<int64_t>> zs;
            for (auto &w : s.v) { zs.push_back(j1a::latticeCoords(G, b, w)); ks.push_back(j1a::latticeKey(G, zs.back())); }
            S.keys.push_back(ks); S.Z.push_back(zs);
        }
        return cache.emplace(bk, std::move(S)).first->second;
    }
    std::vector<const std::vector<double> *> values(BlockSimplices &S, size_t q) {
        std::vector<const std::vector<double> *> V;
        for (auto &z : S.Z[q]) V.push_back(&A.vertex(z));
        return V;
    }
    // neighbour across the facet opposite vertex i; returns false at the boundary of the box
    bool neighbour(const j1a::Block &b, size_t q, int i, j1a::Block &nb, size_t &nq) {
        BlockSimplices &S = simplices(b);
        auto V = values(S, q);
        std::vector<double> fc(n, 0.0);
        for (int j = 0; j <= n; ++j) if (j != i) for (int c = 0; c < n; ++c) fc[c] += (*V[j])[c] / n;
        std::vector<double> p(n);
        for (int c = 0; c < n; ++c) p[c] = fc[c] + 1e-4 * (fc[c] - (*V[i])[c]);
        if (!G.leafAt(p, nb)) return false;
        std::vector<uint64_t> f; for (int j = 0; j <= n; ++j) if (j != i) f.push_back(S.keys[q][j]);
        std::sort(f.begin(), f.end());
        BlockSimplices &T = simplices(nb);
        for (size_t r = 0; r < T.L.size(); ++r) {
            if (nb.level == b.level && G.key(nb) == G.key(b) && r == q) continue;
            ++lookups;
            std::vector<uint64_t> kk = T.keys[r]; std::sort(kk.begin(), kk.end());
            if (std::includes(kk.begin(), kk.end(), f.begin(), f.end())) { nq = r; return true; }
        }
        return false;
    }
};

int main(int argc, char **argv) {
    const int levels = argc > 1 ? std::atoi(argv[1]) : 2;
    const std::string out = (argc > 2 && std::string(argv[2]) != "-") ? argv[2] : "";
    const int n = NDIM, k = KDIM;
    std::vector<int> N(DOMAIN_DIV.begin(), DOMAIN_DIV.end());
    j1a::Grid G(n, N, j1a::Vec(DOMAIN_MIN.begin(), DOMAIN_MIN.end()), j1a::Vec(DOMAIN_MAX.begin(), DOMAIN_MAX.end()), levels);
    const int kr = k + subDim<TC>();
    j1a::Approximator AF(G, k, [](const double *x, double *f) { func(x, f); }, PERTURBATION_AMPLITUDE);
    j1a::Approximator AR(G, kr, refinementMap<TC>(), PERTURBATION_AMPLITUDE);
    uint64_t totalTraced = 0;
    auto t0 = std::chrono::high_resolution_clock::now();

    // trace the component of M through FIRST_POINT; calls visit(block, index, cell) for each crossed simplex
    auto trace = [&](bool build, auto &&visit) {
        Tracer T(G, AF);
        j1a::Block b0; j1a::Vec x0(FIRST_POINT.begin(), FIRST_POINT.end());
        if (!G.leafAt(x0, b0)) { std::cerr << "first point outside\n"; std::exit(1); }
        std::vector<std::pair<uint64_t, uint32_t>> front;
        std::unordered_set<uint64_t> done;
        auto id = [&](const j1a::Block &b, size_t q) { return G.key(b) * 1000003ULL + q; };
        j1a::Cell C;
        { // a crossed simplex of the first block (closest to the first point)
            BlockSimplices &S = T.simplices(b0); double best = 1e300; size_t bq = S.L.size();
            for (size_t q = 0; q < S.L.size(); ++q) {
                auto V = T.values(S, q);
                if (!j1a::buildCell(n, k, V, C, false)) continue;
                double d = 0; for (auto *v : V) for (int c = 0; c < n; ++c) d += ((*v)[c] - x0[c]) * ((*v)[c] - x0[c]);
                if (d < best) { best = d; bq = q; }
            }
            if (bq == S.L.size()) { std::cerr << "no crossed simplex in the first block\n"; std::exit(1); }
            front.push_back({G.key(b0), (uint32_t)bq}); done.insert(id(b0, bq));
        }
        uint64_t count = 0;
        while (!front.empty()) {
            auto [bk, q] = front.back(); front.pop_back();
            j1a::Block b = G.decode(bk);
            BlockSimplices &S = T.simplices(b);
            auto V = T.values(S, q);
            if (!j1a::buildCell(n, k, V, C, build)) continue;
            ++count; visit(b, q, S, C);
            uint64_t facets = 0, all = (1ULL << (n + 1)) - 1;
            for (auto m : C.vmask) facets |= all & ~m;
            for (int i = 0; i <= n; ++i) {
                if (!((facets >> i) & 1)) continue;
                j1a::Block nb; size_t nq;
                if (!T.neighbour(b, q, i, nb, nq)) continue;
                if (done.insert(id(nb, nq)).second) front.push_back({G.key(nb), (uint32_t)nq});
            }
        }
        return count;
    };

    // ---- refinement driven by continuation
    for (int lev = 0; lev < levels; ++lev) {
        std::vector<j1a::Block> mark;
        totalTraced += trace(false, [&](const j1a::Block &b, size_t q, BlockSimplices &S, j1a::Cell &) {
            if (b.level != lev) return;
            std::vector<const std::vector<double> *> V;
            for (auto &z : S.Z[q]) V.push_back(&AR.vertex(z));
            j1a::Cell D;
            if (j1a::buildCell(n, kr, V, D, false)) mark.push_back(b);
        });
        for (auto &b : mark) G.refine(b);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    // ---- final trace with output
    FILE *fp = nullptr;
    if (!out.empty()) {
        fp = std::fopen(out.c_str(), "w");
        std::fprintf(fp, "%d %d\n", n, k);
        for (int i = 0; i < n; ++i) std::fprintf(fp, "%d ", N[i]);
        std::fprintf(fp, "\n\n");
    }
    std::map<int, uint64_t> perLevel;
    uint64_t nCross = trace(fp != nullptr, [&](const j1a::Block &b, size_t q, BlockSimplices &, j1a::Cell &C) {
        perLevel[b.level]++;
        if (!fp) return;
        std::fprintf(fp, "%llu", (unsigned long long)(G.key(b) * 1000003ULL + q));
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
    });
    if (fp) { std::fprintf(fp, "-1\n"); std::fclose(fp); }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::map<int, int> bl; for (auto &b : G.blocks()) bl[b.level]++;
    std::cout << "n=" << n << " k=" << k << " levels=" << levels << "\n";
    std::cout << "blocks:"; for (auto &p : bl) std::cout << " " << p.first << ":" << p.second; std::cout << " total " << G.leaves.size() << "\n";
    std::cout << "crossed simplices: " << nCross << "  per level:"; for (auto &p : perLevel) std::cout << " " << p.first << ":" << p.second; std::cout << "\n";
    std::cout << "simplices traced during refinement: " << totalTraced << "\n";
    std::cout << "F evaluations: " << AF.nEval << " (refinement map: " << AR.nEval << ")\n";
    std::cout << "time refinement: " << std::chrono::duration<double>(t1 - t0).count() << " s, final trace: "
              << std::chrono::duration<double>(t2 - t1).count() << " s\n";
    return 0;
}
