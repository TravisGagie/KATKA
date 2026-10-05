// Run-length tag array (document = genus of each BWT row) with document listing by Muthukrishnan's
// method without the C array (Sadakane's trick, as in Olbrich & Ohlebusch, CPM 2026): an RMQ over
// C[i] = previous run with the same tag, plus marks for the tags already reported.
// Sampling (s > 1): following LF from a run's head stays in the same document except at the first
// character of a document's forward string, so a run's tag equals the tag of the run containing LF of
// its head; the tags of only some runs are stored (roots = runs whose LF leaves the document, at least
// one run on every cycle, and greedily enough that every run reaches a sampled run in fewer than s steps).
#pragma once
#include "rz_index.hpp"
#include <sdsl/rmq_support.hpp>

namespace rz {

struct tag_index {
    u64 n = 0, rho = 0, s = 1;
    sdsl::sd_vector<> RB; sdsl::rank_support_sd<> RBr; sdsl::select_support_sd<> RBs;   // run starts (rows)
    sdsl::bit_vector SB; sdsl::rank_support_v5<> SBr;                                  // sampled runs (s > 1)
    sdsl::int_vector<> L;                                                              // tags of sampled runs
    sdsl::rmq_succinct_sct<> rmq;                                                      // over C
    const rlbwt *bwt = nullptr;                                                        // for LF (s > 1)
    std::vector<char> mark; std::vector<u64> marked;
    u64 lf_steps = 0;

    inline u64 run_of(u64 row) const { return RBr(row + 1) - 1; }
    inline u64 tag(u64 i) {
        if (s == 1) return L[i];
        while (!SB[i]) {
            u64 b = RBs(i + 1); unsigned char c = bwt->at(b);
            i = run_of(bwt->C[c] + bwt->rank(b, c)); ++lf_steps;
        }
        return L[SBr(i)];
    }
    // distinct tags of rows [sp, ep), appended to out (unsorted)
    void list(u64 sp, u64 ep, std::vector<u64> &out) {
        if (mark.empty()) mark.assign(1 << 16, 0);
        u64 a = run_of(sp), b = run_of(ep - 1);
        auto rep = [&](u64 t) { if (t >= mark.size()) mark.resize(t + 1, 0); mark[t] = 1; marked.push_back(t); out.push_back(t); };
        rep(tag(a));
        if (b > a) rec(a + 1, b, out, rep);
        for (u64 t : marked) mark[t] = 0;
        marked.clear();
    }
    template <class F> void rec(u64 l, u64 r, std::vector<u64> &out, F &rep) {
        u64 k = rmq(l, r), t = tag(k);
        if (t < mark.size() && mark[t]) return;      // every run in [l, r] has its tag earlier in the interval
        rep(t);
        if (k > l) rec(l, k - 1, out, rep);
        if (k < r) rec(k + 1, r, out, rep);
    }
    u64 bytes() const {
        return sdsl::size_in_bytes(RB) + sdsl::size_in_bytes(RBr) + sdsl::size_in_bytes(RBs) + (s > 1 ? sdsl::size_in_bytes(SB) + sdsl::size_in_bytes(SBr) : 0)
               + sdsl::size_in_bytes(L) + sdsl::size_in_bytes(rmq);
    }
    void save(const std::string &f) const {
        std::ofstream o(f, std::ios::binary); o.write("RZTAG001", 8); o.write((char *)&n, 8); o.write((char *)&rho, 8); o.write((char *)&s, 8);
        RB.serialize(o); RBr.serialize(o); RBs.serialize(o); if (s > 1) { SB.serialize(o); SBr.serialize(o); } L.serialize(o); rmq.serialize(o);
    }
    bool load(const std::string &f) {
        std::ifstream i(f, std::ios::binary); char mg[8];
        if (!i.read(mg, 8) || std::string(mg, 8) != "RZTAG001") return false;
        i.read((char *)&n, 8); i.read((char *)&rho, 8); i.read((char *)&s, 8);
        RB.load(i); RBr.load(i, &RB); RBs.load(i, &RB); if (s > 1) { SB.load(i); SBr.load(i, &SB); } L.load(i); rmq.load(i);
        return (bool)i;
    }
};

}  // namespace rz
