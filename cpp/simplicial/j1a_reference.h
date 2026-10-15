// J1a: adaptive triangulation of a box in any dimension
// (A. Castelo, L.G. Nonato, M. Siqueira, R. Minghim, G. Tavares, Computers & Graphics 30, 2006).
// Blocks are the leaves of a 2^n-tree over a grid of level-0 blocks, with the balance restriction.
// The simplices of a block are generated from (s, pi) and, in transition blocks, (h, t).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace j1a {

using Vec = std::vector<double>;
using Map = std::function<void(const double *x, double *f)>;

struct Block { int level; std::vector<int64_t> g; };

class Grid {
public:
    int n, maxLevel;
    std::vector<int> N;          // level-0 blocks per direction
    Vec lo, hi, h0;              // box and half-size of level-0 blocks
    std::unordered_set<uint64_t> leaves;

    Grid(int n_, std::vector<int> N_, Vec lo_, Vec hi_, int maxLevel_)
        : n(n_), maxLevel(maxLevel_), N(N_), lo(lo_), hi(hi_), h0(n_) {
        for (int i = 0; i < n; ++i) h0[i] = (hi[i] - lo[i]) / (2.0 * N[i]);
        std::vector<int64_t> g(n, 0);
        while (true) {
            leaves.insert(key({0, g}));
            int i = 0;
            while (i < n && ++g[i] == N[i]) g[i++] = 0;
            if (i == n) break;
        }
    }
    // ---- block keys
    uint64_t key(const Block &b) const {
        uint64_t k = 0;
        for (int i = n - 1; i >= 0; --i) k = k * ((uint64_t)N[i] << maxLevel) + (uint64_t)b.g[i];
        return k * (maxLevel + 1) + b.level;
    }
    Block decode(uint64_t k) const {
        Block b; b.level = (int)(k % (maxLevel + 1)); k /= (maxLevel + 1); b.g.resize(n);
        for (int i = 0; i < n; ++i) { uint64_t r = (uint64_t)N[i] << maxLevel; b.g[i] = (int64_t)(k % r); k /= r; }
        return b;
    }
    double hs(int level, int i) const { return h0[i] / double(1 << level); }
    Vec center(const Block &b) const {
        Vec c(n);
        for (int i = 0; i < n; ++i) c[i] = lo[i] + (2 * b.g[i] + 1) * hs(b.level, i);
        return c;
    }
    bool inBox(const Vec &x) const {
        for (int i = 0; i < n; ++i) if (x[i] < lo[i] || x[i] > hi[i]) return false;
        return true;
    }
    // leaf containing x (x strictly inside a block)
    bool leafAt(const Vec &x, Block &out) const {
        if (!inBox(x)) return false;
        for (int r = 0; r <= maxLevel; ++r) {
            Block b; b.level = r; b.g.resize(n);
            for (int i = 0; i < n; ++i) {
                int64_t gi = (int64_t)std::floor((x[i] - lo[i]) / (2 * hs(r, i)));
                b.g[i] = std::min<int64_t>(gi, ((int64_t)N[i] << r) - 1);
            }
            if (leaves.count(key(b))) { out = b; return true; }
        }
        return false;
    }
    void split(const Block &b) {
        leaves.erase(key(b));
        for (int o = 0; o < (1 << n); ++o) {
            Block c; c.level = b.level + 1; c.g.resize(n);
            for (int i = 0; i < n; ++i) c.g[i] = 2 * b.g[i] + ((o >> i) & 1);
            leaves.insert(key(c));
        }
    }
    // refine b, enforcing the balance restriction (blocks that meet differ by at most one level)
    void refine(const Block &b) {
        if (b.level >= maxLevel || !leaves.count(key(b))) return;
        Vec c = center(b);
        int nd = 1; for (int i = 0; i < n; ++i) nd *= 3;
        for (int d = 0; d < nd; ++d) {
            int dd = d; Vec x(n); bool zero = true;
            for (int i = 0; i < n; ++i) {
                int di = dd % 3 - 1; dd /= 3; if (di) zero = false;
                x[i] = c[i] + di * hs(b.level, i) * 1.125;
            }
            if (zero) continue;
            Block nb;
            if (leafAt(x, nb) && nb.level < b.level) refine(nb);
        }
        if (leaves.count(key(b))) split(b);
    }
    std::vector<Block> blocks() const {
        std::vector<Block> v; v.reserve(leaves.size());
        for (auto k : leaves) v.push_back(decode(k));
        return v;
    }
    // ---- refined faces: masks over 2n signed axes (bit 2i: x_i=-1, bit 2i+1: x_i=+1)
    std::vector<uint32_t> refinedFaces(const Block &b) const {
        std::vector<uint32_t> out;
        Vec c = center(b);
        int nf = 1; for (int i = 0; i < n; ++i) nf *= 3;
        for (int f = 0; f < nf; ++f) {
            // digit 0: free, 1: side -1, 2: side +1
            int ff = f, fixed = 0; uint32_t mask = 0; std::vector<int> side(n);
            for (int i = 0; i < n; ++i) { side[i] = ff % 3; ff /= 3; if (side[i]) { ++fixed; mask |= 1u << (2 * i + (side[i] == 2)); } }
            if (fixed == 0 || fixed == n) continue;      // the block itself, or a vertex
            std::vector<int> fx; for (int i = 0; i < n; ++i) if (side[i]) fx.push_back(i);
            bool ref = false;
            for (int sub = 0; sub < (1 << fixed) && !ref; ++sub) {   // blocks around the face
                Vec x(n);
                for (int i = 0; i < n; ++i) x[i] = c[i] + hs(b.level, i) * (side[i] ? 0.0 : 0.4375);
                for (int j = 0; j < fixed; ++j) {
                    int i = fx[j]; double sg = side[i] == 2 ? 1.0 : -1.0;
                    x[i] = c[i] + sg * hs(b.level, i) * (((sub >> j) & 1) ? 1.0625 : 0.9375);
                }
                Block nb;
                if (leafAt(x, nb) && nb.level == b.level + 1) ref = true;
            }
            if (ref) out.push_back(mask);
        }
        return out;
    }
};

// A simplex of J1a in the frame of its block, with integer coordinates in units of hs/2 (values in -2..2)
struct LocalSimplex {
    std::vector<int> s, pi, t; int h;          // h = n for simplices that are not subdivided
    std::vector<std::vector<int>> v;           // n+1 vertices
};

inline void enumerateSimplices(int n, const std::vector<uint32_t> &Phi, std::vector<LocalSimplex> &out) {
    out.clear();
    std::vector<int> pi(n); std::iota(pi.begin(), pi.end(), 0);
    do {
        for (int sm = 0; sm < (1 << n); ++sm) {
            std::vector<int> s(n);
            for (int i = 0; i < n; ++i) s[i] = ((sm >> i) & 1) ? 1 : -1;
            int h = n; uint32_t prefix = 0;
            for (int j = 0; j < n && h == n; ++j) {
                prefix |= 1u << (2 * pi[j] + (s[pi[j]] > 0));
                for (uint32_t phi : Phi) if ((phi & ~prefix) == 0) { h = j + 1; break; }
            }
            if (h == n) {
                LocalSimplex L; L.s = s; L.pi = pi; L.h = n;
                std::vector<int> v(n, 0); L.v.push_back(v);
                for (int j = 0; j < n; ++j) { v[pi[j]] += 2 * s[pi[j]]; L.v.push_back(v); }
                out.push_back(L);
                continue;
            }
            for (int tm = 0; tm < (1 << (n - h)); ++tm) {
                LocalSimplex L; L.s = s; L.pi = pi; L.h = h; L.t.resize(n - h);
                for (int j = 0; j < n - h; ++j) L.t[j] = ((tm >> j) & 1) ? 1 : -1;
                std::vector<int> v(n, 0); L.v.push_back(v);
                for (int j = 0; j < h - 1; ++j) { v[pi[j]] += 2 * s[pi[j]]; L.v.push_back(v); }
                v[pi[h - 1]] += 2 * s[pi[h - 1]];
                for (int j = h; j < n; ++j) v[pi[j]] += s[pi[j]];       // centre of the sub-face
                L.v.push_back(v);
                for (int j = h; j < n; ++j) { v[pi[j]] += L.t[j - h]; L.v.push_back(v); }
                out.push_back(L);
            }
        }
    } while (std::next_permutation(pi.begin(), pi.end()));
}

// global integer lattice: unit = (level-0 half size) / 2^(maxLevel+1)
inline std::vector<int64_t> latticeCoords(const Grid &G, const Block &b, const std::vector<int> &w) {
    std::vector<int64_t> z(G.n);
    const int64_t scale = int64_t(1) << (G.maxLevel - b.level);
    for (int i = 0; i < G.n; ++i) z[i] = (2 * b.g[i] + 1) * 2 * scale + w[i] * scale;
    return z;
}
inline uint64_t latticeKey(const Grid &G, const std::vector<int64_t> &z) {
    uint64_t k = 0;
    for (int i = G.n - 1; i >= 0; --i) k = k * (((uint64_t)G.N[i] << (G.maxLevel + 2)) + 1) + (uint64_t)z[i];
    return k;
}

// Piecewise linear approximation of F^{-1}(0) on J1a, with vertex values cached and perturbed
class Approximator {
public:
    const Grid &G; int k; Map F; double eps;
    std::unordered_map<uint64_t, std::vector<double>> cache;   // x (n) followed by F(x) (k)
    uint64_t nEval = 0;
    Approximator(const Grid &G_, int k_, Map F_, double eps_ = 1e-9) : G(G_), k(k_), F(F_), eps(eps_) {}
    const std::vector<double> &vertex(const std::vector<int64_t> &z) {
        uint64_t key = latticeKey(G, z);
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
        std::vector<double> d(G.n + k);
        const double unit = 1.0 / double(int64_t(1) << (G.maxLevel + 1));
        std::mt19937_64 rng(key * 0x9E3779B97F4A7C15ULL + 12345);
        std::uniform_real_distribution<double> U(-1, 1);
        for (int i = 0; i < G.n; ++i) d[i] = G.lo[i] + z[i] * unit * G.h0[i] + eps * G.h0[i] * U(rng);
        F(d.data(), d.data() + G.n); ++nEval;
        return cache.emplace(key, std::move(d)).first->second;
    }
};

// zero of the affine interpolant on a k-face; returns false if not crossed
inline bool kFaceZero(int n, int k, const std::vector<const std::vector<double> *> &V, std::vector<double> *pt) {
    const int M = k + 1;
    std::vector<double> A(M * (M + 1));
    auto a = [&](int r, int c) -> double & { return A[r * (M + 1) + c]; };
    for (int c = 0; c < M; ++c) { a(0, c) = 1; for (int r = 0; r < k; ++r) a(r + 1, c) = (*V[c])[n + r]; }
    a(0, M) = 1; for (int r = 1; r < M; ++r) a(r, M) = 0;
    for (int c = 0; c < M; ++c) {
        int p = c; double best = std::fabs(a(c, c));
        for (int r = c + 1; r < M; ++r) if (std::fabs(a(r, c)) > best) { best = std::fabs(a(r, c)); p = r; }
        if (best < 1e-300) return false;
        if (p != c) for (int q = 0; q <= M; ++q) std::swap(a(c, q), a(p, q));
        for (int r = c + 1; r < M; ++r) { double f = a(r, c) / a(c, c); for (int q = c; q <= M; ++q) a(r, q) -= f * a(c, q); }
    }
    std::vector<double> lam(M);
    for (int r = M - 1; r >= 0; --r) {
        double s = a(r, M); for (int q = r + 1; q < M; ++q) s -= a(r, q) * lam[q];
        lam[r] = s / a(r, r); if (!(lam[r] >= 0)) return false;
    }
    if (pt) { pt->assign(n, 0.0); for (int c = 0; c < M; ++c) for (int i = 0; i < n; ++i) (*pt)[i] += lam[c] * (*V[c])[i]; }
    return true;
}

// combinatorial skeleton of the cell of the approximation inside one n-simplex
struct Cell {
    std::vector<uint64_t> vmask; std::vector<std::vector<double>> pts;
    std::vector<std::vector<uint64_t>> cmask; std::vector<std::vector<std::vector<int>>> bnd;
};
inline bool buildCell(int n, int k, const std::vector<const std::vector<double> *> &V, Cell &C, bool build) {
    C.vmask.clear(); C.pts.clear();
    std::vector<int> idx(k + 1); std::iota(idx.begin(), idx.end(), 0);
    std::vector<const std::vector<double> *> W(k + 1); std::vector<double> pt;
    while (true) {
        for (int i = 0; i <= k; ++i) W[i] = V[idx[i]];
        if (kFaceZero(n, k, W, build ? &pt : nullptr)) {
            uint64_t m = 0; for (int i = 0; i <= k; ++i) m |= 1ULL << idx[i];
            C.vmask.push_back(m); if (build) C.pts.push_back(pt);
        }
        int t = k; while (t >= 0 && idx[t] == n - k + t) --t;
        if (t < 0) break;
        ++idx[t]; for (int j = t + 1; j <= k; ++j) idx[j] = idx[j - 1] + 1;
    }
    if (C.vmask.empty() || !build) return !C.vmask.empty();
    const int m = n - k; const uint64_t all = (1ULL << (n + 1)) - 1;
    C.cmask.assign(m + 1, {}); C.bnd.assign(m + 1, {});
    C.cmask[0] = C.vmask;
    for (int j = 1; j <= m; ++j)
        for (uint64_t Mk = 1; Mk <= all; ++Mk) {
            if (__builtin_popcountll(Mk) != k + 1 + j) continue;
            std::vector<int> b;
            for (size_t q = 0; q < C.cmask[j - 1].size(); ++q) if ((C.cmask[j - 1][q] & ~Mk) == 0) b.push_back((int)q);
            if (b.empty()) continue;
            C.cmask[j].push_back(Mk); C.bnd[j].push_back(b);
        }
    return true;
}

} // namespace j1a
