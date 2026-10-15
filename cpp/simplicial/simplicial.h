// Marching Simplex and Continuation Simplex on the K1 (Freudenthal-Kuhn) triangulation
// of the perturbed grid used by the pta/gcch/fch programs.  Shared definitions.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include "../common/definitions.h"
#include "../common/utility.h"


// open-addressing hash set of 64-bit keys
struct HashSet {
    std::vector<uint64_t> t; size_t cnt = 0, mask;
    explicit HashSet(int lg = 10) { t.assign(1ULL << lg, 0); mask = t.size() - 1; }
    static uint64_t mix(uint64_t x) { x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31); }
    bool insert(uint64_t key) {
        key += 1;
        if (10 * (cnt + 1) > 9 * t.size()) grow();
        size_t h = mix(key) & mask;
        while (t[h]) { if (t[h] == key) return false; h = (h + 1) & mask; }
        t[h] = key; ++cnt; return true;
    }
    void grow() {
        std::vector<uint64_t> o; o.swap(t); t.assign(o.size() * 2, 0); mask = t.size() - 1; cnt = 0;
        for (uint64_t k : o) if (k) { size_t h = mix(k) & mask; while (t[h]) h = (h + 1) & mask; t[h] = k; ++cnt; }
    }
};

// n-simplex of K1: cube y and permutation p; vertices v_0=y, v_j=v_{j-1}+e_{p_{j-1}}
struct Simplex { std::array<int64_t, NDIM> y; std::array<uint8_t, NDIM> p; };

// Triangulation: K1 (default) or J1 (compile with -DJ1_TRIANGULATION).
// In K1 every cube is a translate of CFK; in J1 the cube g is reflected so that v0 is its even corner
// and the sign vector is s_j = (-1)^{g_j} (Chap. 4).
#ifdef J1_TRIANGULATION
constexpr bool IS_J1 = true;
#else
constexpr bool IS_J1 = false;
#endif

inline size_t parityMask(const Simplex &s) {
    size_t b = 0;
    for (size_t i = 0; i < NDIM; ++i) if (s.y[i] & 1) b |= (1ULL << i);
    return b;
}
inline int signOf(const Simplex &s, size_t dir) { return (s.y[dir] & 1) ? -1 : 1; }

inline uint64_t factorial(size_t n) { uint64_t f = 1; for (size_t i = 2; i <= n; ++i) f *= i; return f; }

inline uint64_t cubeLabel(const Simplex &s) {
    uint64_t c = 0;
    for (size_t i = 0; i < NDIM; ++i) c = c * DOMAIN_DIV[i] + (uint64_t)s.y[i];
    return c;
}
inline uint64_t permLabel(const Simplex &s) {
    uint64_t l = 0;
    for (size_t i = 0; i < NDIM; ++i) {
        uint64_t smaller = 0;
        for (size_t j = i + 1; j < NDIM; ++j) smaller += (s.p[j] < s.p[i]);
        l = l * (NDIM - i) + smaller;
    }
    return l;
}
inline uint64_t key(const Simplex &s) { return cubeLabel(s) * factorial(NDIM) + permLabel(s); }

inline Simplex decode(uint64_t k) {
    const uint64_t fact = factorial(NDIM);
    uint64_t l = k % fact, c = k / fact;
    Simplex s;
    for (int64_t i = NDIM - 1; i >= 0; --i) { s.y[i] = (int64_t)(c % DOMAIN_DIV[i]); c /= DOMAIN_DIV[i]; }
    std::array<uint64_t, NDIM> code{};
    for (int64_t i = NDIM - 1; i >= 0; --i) { code[i] = l % (NDIM - i); l /= (NDIM - i); }
    std::vector<uint8_t> avail; for (size_t i = 0; i < NDIM; ++i) avail.push_back((uint8_t)i);
    for (size_t i = 0; i < NDIM; ++i) { s.p[i] = avail[code[i]]; avail.erase(avail.begin() + code[i]); }
    return s;
}

inline bool inside(const Simplex &s) {
    for (size_t i = 0; i < NDIM; ++i) if (s.y[i] < 0 || s.y[i] >= (int64_t)DOMAIN_DIV[i]) return false;
    return true;
}

// pivoting rules of K1 and J1 across the facet opposite vertex i (Props. 4.2 and 4.3)
inline Simplex pivot(const Simplex &s, size_t i) {
    Simplex r = s;
    if (IS_J1) {
        if (i == 0) r.y[s.p[0]] += signOf(s, s.p[0]);
        else if (i == NDIM) r.y[s.p[NDIM - 1]] -= signOf(s, s.p[NDIM - 1]);
        else std::swap(r.p[i - 1], r.p[i]);
        return r;
    }
    if (i == 0) {
        r.y[s.p[0]] += 1;
        for (size_t j = 0; j + 1 < NDIM; ++j) r.p[j] = s.p[j + 1];
        r.p[NDIM - 1] = s.p[0];
    } else if (i == NDIM) {
        r.y[s.p[NDIM - 1]] -= 1;
        r.p[0] = s.p[NDIM - 1];
        for (size_t j = 1; j < NDIM; ++j) r.p[j] = s.p[j - 1];
    } else {
        std::swap(r.p[i - 1], r.p[i]);
    }
    return r;
}

// local (cube) binary label of vertex j of the simplex
inline size_t localLabel(const Simplex &s, size_t j) {
    size_t b = IS_J1 ? parityMask(s) : 0;
    for (size_t i = 0; i < j; ++i) b ^= (1ULL << s.p[i]);   // K1: sets the bit; J1: moves towards the odd side
    return b;
}

// perturbed coordinates of a vertex given by cube y and local label b
inline void vertexCoords(const std::array<int64_t, NDIM> &y, size_t b, std::array<double, NDIM> &x) {
    std::array<UINT_COORD, NDIM> g{};
    for (size_t i = 0; i < NDIM; ++i) {
        g[i] = (UINT_COORD)y[i];
        x[i] = DOMAIN_MIN[i] + (double)(y[i] + (int64_t)((b >> i) & 1)) * DOMAIN_STEP[i];
    }
    addPerturbationVertex(g, b, x);
}

struct SimplexData {
    std::array<std::array<double, NDIM>, NDIM + 1> x;   // vertex coordinates
    std::array<std::array<double, KDIM>, NDIM + 1> f;   // F at the vertices
    std::array<size_t, NDIM + 1> lab;                   // local labels
};

// zero of the affine interpolant on the k-face idx; returns false if not transverse
inline bool kFaceZero(const SimplexData &d, const std::array<size_t, KDIM + 1> &idx, std::array<double, NDIM> *pt) {
    constexpr size_t M = KDIM + 1;
    double A[M][M + 1];
    for (size_t c = 0; c < M; ++c) {
        A[0][c] = 1.0;
        for (size_t r = 0; r < KDIM; ++r) A[r + 1][c] = d.f[idx[c]][r];
    }
    A[0][M] = 1.0;
    for (size_t r = 1; r < M; ++r) A[r][M] = 0.0;
    for (size_t c = 0; c < M; ++c) {
        size_t piv = c; double best = std::abs(A[c][c]);
        for (size_t r = c + 1; r < M; ++r) if (std::abs(A[r][c]) > best) { best = std::abs(A[r][c]); piv = r; }
        if (best < 1e-300) return false;
        if (piv != c) for (size_t q = 0; q <= M; ++q) std::swap(A[c][q], A[piv][q]);
        for (size_t r = c + 1; r < M; ++r) {
            const double fct = A[r][c] / A[c][c];
            for (size_t q = c; q <= M; ++q) A[r][q] -= fct * A[c][q];
        }
    }
    double lam[M];
    for (int64_t r = M - 1; r >= 0; --r) {
        double s = A[r][M];
        for (size_t q = r + 1; q < M; ++q) s -= A[r][q] * lam[q];
        lam[r] = s / A[r][r];
        if (!(lam[r] >= 0)) return false;
    }
    if (pt) {
        pt->fill(0.0);
        for (size_t c = 0; c < M; ++c) for (size_t i = 0; i < NDIM; ++i) (*pt)[i] += lam[c] * d.x[idx[c]][i];
    }
    return true;
}

// result of processing one simplex
struct Polytope {
    std::vector<uint64_t> vmask;                     // transverse k-faces (vertex masks)
    std::vector<std::array<double, NDIM>> pts;       // their zeros
    std::array<std::vector<uint64_t>, MDIM + 1> cmask;            // j-cells: masks of (k+j)-faces
    std::array<std::vector<std::vector<size_t>>, MDIM + 1> bnd;   // boundaries (0-based)
    uint64_t facets = 0;                             // bit i: facet opposite v_i is crossed
};

// test the k-faces of the simplex; if build, also build the combinatorial skeleton
inline bool processSimplex(const SimplexData &d, Polytope &P, bool build) {
    P.vmask.clear(); P.pts.clear(); P.facets = 0;
    const uint64_t all = (1ULL << (NDIM + 1)) - 1;
    std::array<size_t, KDIM + 1> idx{};
    for (size_t i = 0; i <= KDIM; ++i) idx[i] = i;
    std::array<double, NDIM> pt;
    while (true) {
        if (kFaceZero(d, idx, build ? &pt : nullptr)) {
            uint64_t m = 0;
            for (size_t i = 0; i <= KDIM; ++i) m |= 1ULL << idx[i];
            P.vmask.push_back(m);
            if (build) P.pts.push_back(pt);
            P.facets |= (all & ~m);
        }
        int64_t t = KDIM;
        while (t >= 0 && idx[t] == NDIM - KDIM + (size_t)t) --t;
        if (t < 0) break;
        ++idx[t];
        for (size_t j = t + 1; j <= KDIM; ++j) idx[j] = idx[j - 1] + 1;
    }
    if (P.vmask.empty() || !build) return !P.vmask.empty();
    // j-cells lie in the (k+j)-faces that contain a crossed k-face (Sect. "combinatorial skeleton")
    P.cmask[0] = P.vmask;
    for (size_t j = 1; j <= MDIM; ++j) {
        P.cmask[j].clear(); P.bnd[j].clear();
        for (uint64_t M = 1; M <= all; ++M) {
            if ((size_t)__builtin_popcountll(M) != KDIM + 1 + j) continue;
            std::vector<size_t> b;
            for (size_t q = 0; q < P.cmask[j - 1].size(); ++q)
                if ((P.cmask[j - 1][q] & ~M) == 0) b.push_back(q);
            if (b.empty()) continue;
            P.cmask[j].push_back(M);
            P.bnd[j].push_back(b);
        }
    }
    return true;
}

// writer of the .pol format (Appendix A)
struct PolWriter {
    FILE *fp = nullptr;
    explicit PolWriter(const std::string &fn) {
        if (fn.empty()) return;
        fp = std::fopen(fn.c_str(), "w");
        std::fprintf(fp, "%zu %zu\n", NDIM, KDIM);
        for (size_t i = 0; i < NDIM; ++i) std::fprintf(fp, "%zu ", DOMAIN_DIV[i]);
        std::fprintf(fp, "\n\n");
    }
    void write(const Simplex &s, const SimplexData &d, const Polytope &P) {
        if (!fp) return;
        std::fprintf(fp, "%llu", (unsigned long long)key(s));
        for (size_t i = 0; i < NDIM; ++i) std::fprintf(fp, " %lld", (long long)s.y[i]);
        std::fprintf(fp, "\n1\n%zu\n", P.vmask.size());
        for (size_t q = 0; q < P.vmask.size(); ++q) {
            for (size_t i = 0; i <= NDIM; ++i) if (P.vmask[q] >> i & 1) std::fprintf(fp, "%zu ", d.lab[i]);
            for (size_t i = 0; i < NDIM; ++i) std::fprintf(fp, " %.12g", P.pts[q][i]);
            std::fprintf(fp, "\n");
        }
        for (size_t j = 1; j <= MDIM; ++j) {
            std::fprintf(fp, "%zu\n", P.cmask[j].size());
            for (auto &b : P.bnd[j]) {
                for (size_t q : b) std::fprintf(fp, "%zu ", q + 1);
                std::fprintf(fp, "\n");
            }
        }
    }
    ~PolWriter() { if (fp) { std::fprintf(fp, "-1\n"); std::fclose(fp); } }
};
