// Run-length tag array (document = genus of each BWT row) with document listing by Muthukrishnan's
// method without the C array (Sadakane's trick, as in Olbrich & Ohlebusch, CPM 2026): an RMQ over
// C[i] = previous run with the same tag, plus marks for the tags already reported.
// Sampling (s > 1): following LF from a run's head stays in the same document except at the first
// character of a document's forward string, so a run's tag equals the tag of the run containing LF of
// its head; the tags of only some runs are stored (roots = runs whose LF leaves the document, at least
// one run on every cycle, and greedily enough that every run reaches a sampled run in fewer than s steps).
// Counting needs no RMQ, so it is stored in its own file, out.tag.rmq, built only for listing (rz-tagbuild
// without -n) and loaded only when listing; files of the old format (RZTAG001) carry it at the end.
#pragma once
#include "rz_index.hpp"
#include <sdsl/rmq_support.hpp>
#include "gtag.hpp"
#include <memory>

namespace rz {

struct tag_index {
    u64 n = 0, rho = 0, s = 1;
    sdsl::sd_vector<> RB; sdsl::rank_support_sd<> RBr; sdsl::select_support_sd<> RBs;   // run starts (rows)
    sdsl::bit_vector SB; sdsl::rank_support_v5<> SBr;                                  // sampled runs (s > 1)
    sdsl::int_vector<> L;                                                              // tags of sampled runs
    sdsl::rmq_succinct_sct<> rmq; bool has_rmq = false;                                // over C (listing only)
    std::shared_ptr<gram_tags> gram;                                                   // optional: genera of the runs from a grammar (gtag.hpp; counting only)
    const rlbwt *bwt = nullptr;                                                        // for LF (s > 1)
    static inline thread_local std::vector<char> mark; static inline thread_local std::vector<u64> marked;   // per thread (rz-classify -j)
    static inline thread_local u64 lf_steps = 0;

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
    // number of rows of [sp, ep) with each tag, as (tag, count) pairs sorted by tag
    static inline thread_local std::vector<uint32_t> cnt;   // per thread (rz-classify -j)
    void count(u64 sp, u64 ep, std::vector<std::pair<u64, u64>> &out) {
        if (cnt.empty()) cnt.assign(1 << 16, 0);
        out.clear();
        u64 a = run_of(sp), b = run_of(ep - 1);
        if (gram) {   // genera from the grammar, decoded sequentially from run a
            u64 i = a, en = RBs(a + 1);
            gram->decode(a, b - a + 1, [&](uint32_t t) {
                u64 st = std::max(sp, en); en = i + 1 < rho ? (u64)RBs(i + 2) : n; ++i;
                if (t >= cnt.size()) cnt.resize(t + 1, 0);
                if (cnt[t] == 0) marked.push_back(t);
                cnt[t] += (uint32_t)(std::min(ep, en) - st);
            });
        } else
        for (u64 i = a; i <= b; ++i) {
            u64 st = std::max(sp, (u64)RBs(i + 1)), en = std::min(ep, i + 1 < rho ? (u64)RBs(i + 2) : n);
            u64 t = tag(i);
            if (t >= cnt.size()) cnt.resize(t + 1, 0);
            if (cnt[t] == 0) marked.push_back(t);
            cnt[t] += (uint32_t)(en - st);
        }
        for (u64 t : marked) { out.push_back({t, cnt[t]}); cnt[t] = 0; }
        marked.clear();
        std::sort(out.begin(), out.end());
    }
    u64 bytes() const {
        return sdsl::size_in_bytes(RB) + sdsl::size_in_bytes(RBr) + sdsl::size_in_bytes(RBs) + (s > 1 ? sdsl::size_in_bytes(SB) + sdsl::size_in_bytes(SBr) : 0)
               + sdsl::size_in_bytes(L) + (has_rmq ? sdsl::size_in_bytes(rmq) : 0) + (gram ? gram->bytes() : 0);
    }
    void save(const std::string &f) const {
        std::ofstream o(f, std::ios::binary); o.write("RZTAG002", 8); o.write((char *)&n, 8); o.write((char *)&rho, 8); o.write((char *)&s, 8);
        RB.serialize(o); RBr.serialize(o); RBs.serialize(o); if (s > 1) { SB.serialize(o); SBr.serialize(o); } L.serialize(o);
        if (has_rmq) { std::ofstream r(f + ".rmq", std::ios::binary); r.write("RZRMQ001", 8); rmq.serialize(r); }
    }
    // with the genera of the runs from a grammar instead of L (counting only; s must be 1)
    void save_gram(const std::string &f) const {
        std::ofstream o(f, std::ios::binary); o.write("RZGTAG01", 8); o.write((char *)&n, 8); o.write((char *)&rho, 8); o.write((char *)&s, 8);
        RB.serialize(o); RBr.serialize(o); RBs.serialize(o); gram->serialize(o);
    }
    // need_rmq: load the RMQ too (for listing); fails if it was not built
    bool load(const std::string &f, bool need_rmq = true) {
        std::ifstream i(f, std::ios::binary); char mg[8];
        if (!i.read(mg, 8)) return false;
        std::string m(mg, 8);
        if (m == "RZGTAG01") {
            if (need_rmq) { fprintf(stderr, "%s is a grammar-compressed tag array, which supports counting only (RZ_COUNTS=1)\n", f.c_str()); return false; }
            i.read((char *)&n, 8); i.read((char *)&rho, 8); i.read((char *)&s, 8);
            RB.load(i); RBr.load(i, &RB); RBs.load(i, &RB);
            gram = std::make_shared<gram_tags>(); gram->load(i);
            return (bool)i && gram->rho == rho;
        }
        if (m != "RZTAG001" && m != "RZTAG002") return false;
        i.read((char *)&n, 8); i.read((char *)&rho, 8); i.read((char *)&s, 8);
        RB.load(i); RBr.load(i, &RB); RBs.load(i, &RB); if (s > 1) { SB.load(i); SBr.load(i, &SB); } L.load(i);
        if (!i) return false;
        if (!need_rmq) return true;
        if (m == "RZTAG001") { rmq.load(i); has_rmq = (bool)i; return has_rmq; }
        std::ifstream r(f + ".rmq", std::ios::binary);
        if (!r.read(mg, 8) || std::string(mg, 8) != "RZRMQ001") { fprintf(stderr, "%s.rmq missing: listing needs the RMQ (rz-tagbuild without -n)\n", f.c_str()); return false; }
        rmq.load(r); has_rmq = (bool)r; return has_rmq;
    }
};

}  // namespace rz
