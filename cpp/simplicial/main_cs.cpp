// Continuation Simplex: visits only the n-simplices crossed by the component of the
// approximation that contains FIRST_POINT, using the pivoting rules of K1.
#include <chrono>
#include <cstdlib>
#include <iostream>
#include "simplicial.h"

static uint64_t n_eval = 0;

static void evalSimplex(const Simplex &s, SimplexData &d) {
    for (size_t j = 0; j <= NDIM; ++j) {
        d.lab[j] = localLabel(s, j);
        vertexCoords(s.y, d.lab[j], d.x[j]);
        func(d.x[j].data(), d.f[j].data());
        ++n_eval;
    }
}

int main(int argc, char **argv) {
    // usage: cs [out.pol|-] [log2 initial hash size]
    std::string out = (argc > 1 && std::string(argv[1]) != "-") ? argv[1] : "";
    const int lg = (argc > 2) ? std::atoi(argv[2]) : 10;
    PolWriter W(out);
    const bool build = !out.empty();

    // simplex of K1 containing FIRST_POINT (sorting the fractional coordinates)
    Simplex s0;
    std::array<double, NDIM> frac{};
    std::array<size_t, NDIM> ord{};
    for (size_t i = 0; i < NDIM; ++i) {
        double t = (FIRST_POINT[i] - DOMAIN_MIN[i]) / DOMAIN_STEP[i];
        s0.y[i] = std::min<int64_t>((int64_t)std::floor(t), DOMAIN_DIV[i] - 1);
        frac[i] = t - s0.y[i];
        ord[i] = i;
    }
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return frac[a] > frac[b]; });
    for (size_t i = 0; i < NDIM; ++i) s0.p[i] = (uint8_t)ord[i];

    auto t0 = std::chrono::high_resolution_clock::now();
    SimplexData d; Polytope P;
    std::vector<uint64_t> front;
    {   // if necessary, walk to a nearby crossed simplex
        std::vector<Simplex> q{s0}; HashSet seen; seen.insert(key(s0));
        for (size_t qi = 0; qi < q.size() && qi < 100000; ++qi) {
            evalSimplex(q[qi], d);
            if (processSimplex(d, P, false)) { front.push_back(key(q[qi])); break; }
            for (size_t i = 0; i <= NDIM; ++i) {
                Simplex r = pivot(q[qi], i);
                if (inside(r) && seen.insert(key(r))) q.push_back(r);
            }
        }
        if (front.empty()) { std::cout << "no crossed simplex found" << std::endl; return 1; }
    }
    HashSet done(lg);
    done.insert(front[0]);
    uint64_t n_simplices = 0;
    while (!front.empty()) {
        Simplex s = decode(front.back()); front.pop_back();
        evalSimplex(s, d);
        if (!processSimplex(d, P, build)) continue;
        ++n_simplices;
        W.write(s, d, P);
        for (size_t i = 0; i <= NDIM; ++i) {
            if (!(P.facets >> i & 1)) continue;
            Simplex r = pivot(s, i);
            if (!inside(r)) continue;
            const uint64_t kr = key(r);
            if (done.insert(kr)) front.push_back(kr);
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time (CS): " << std::chrono::duration<double>(t1 - t0).count() << " s" << std::endl;
    std::cout << "crossed n-simplices: " << n_simplices << std::endl;
    std::cout << "F evaluations: " << n_eval << std::endl;
    return 0;
}
