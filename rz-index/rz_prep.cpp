// Build S = G_1 X rc(G_1) X ... G_k X rc(G_k) X from a list of genome FASTA files.
//
// usage: rz-prep [-d] <list.tsv> <out-prefix>
//   -d : also write one dataset file per species (maximal run of consecutive lines with the
//        same species) as <out>.ds<j>.S, listed in order in <out>.datasets, for building a
//        multi-string BWT with pfp-merge (one string per dataset).
//   -s : each FASTA file may hold many genomes; split it into genomes by sample ID (header text
//        before the first '.'), e.g. AllTheBacteria's SAMEA1410869.contig00001 -> SAMEA1410869.
//   -U : homopolymer compression: replace every run of equal characters of each genome by one character
//        (before taking the reverse complement, which commutes with it); classify with rz-classify -U.
//   list.tsv : one genome per line (one file per line; with -s, one or more genomes per file), in tree (left-to-right) order:  <species>\t<fasta[.gz]>
//              (a line with a single field is taken as the FASTA path, species "-")
// writes:
//   <out>.S     the text (bytes A,C,G,T,X; no terminator)
//   <out>.tbl   one line per genome: index, species, file, length, start-in-S
//
// Each FASTA file is one genome; its records (contigs) are joined with X.
// Characters other than ACGT (case-insensitive) become X.

#include <zlib.h>
#include <algorithm>
#include "digest.hpp"
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

static std::string read_genome(const std::string &path) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) { std::cerr << "cannot open " << path << "\n"; exit(1); }
    std::string g;
    std::vector<char> buf(1 << 20);
    bool header = false, bol = true, any = false;
    int got;
    while ((got = gzread(f, buf.data(), buf.size())) > 0) {
        for (int i = 0; i < got; ++i) {
            char c = buf[i];
            if (bol && c == '>') { header = true; if (any) g.push_back('X'); }
            if (c == '\n') { header = false; bol = true; continue; }
            bol = false;
            if (header || c == '\r' || c == ' ' || c == '\t') continue;
            switch (c) {
                case 'A': case 'a': g.push_back('A'); break;
                case 'C': case 'c': g.push_back('C'); break;
                case 'G': case 'g': g.push_back('G'); break;
                case 'T': case 't': g.push_back('T'); break;
                default: g.push_back('X');
            }
            any = true;
        }
    }
    gzclose(f);
    return g;
}

// Read a multi-genome FASTA file; a genome is a maximal run of consecutive records whose
// headers share the same sample ID (the header text up to the first '.' or whitespace).
static void read_genomes(const std::string &path, std::vector<std::pair<std::string, std::string>> &out) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) { std::cerr << "cannot open " << path << "\n"; exit(1); }
    std::vector<char> buf(1 << 20);
    bool header = false, bol = true;
    std::string hdr;
    int got;
    auto start_record = [&](const std::string &h) {
        std::string id = h.substr(0, h.find_first_of(". \t"));
        if (out.empty() || out.back().first != id) out.push_back({id, std::string()});
        else out.back().second.push_back('X');
    };
    while ((got = gzread(f, buf.data(), buf.size())) > 0) {
        for (int i = 0; i < got; ++i) {
            char c = buf[i];
            if (bol && c == '>') { header = true; hdr.clear(); bol = false; continue; }
            if (c == '\n') { if (header) start_record(hdr); header = false; bol = true; continue; }
            bol = false;
            if (header) { hdr.push_back(c); continue; }
            if (c == '\r' || c == ' ' || c == '\t') continue;
            if (out.empty()) out.push_back({"-", std::string()});
            std::string &g = out.back().second;
            switch (c) {
                case 'A': case 'a': g.push_back('A'); break;
                case 'C': case 'c': g.push_back('C'); break;
                case 'G': case 'g': g.push_back('G'); break;
                case 'T': case 't': g.push_back('T'); break;
                default: g.push_back('X');
            }
        }
    }
    gzclose(f);
}

static char comp(char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'X'; }
}

int main(int argc, char **argv) {
    bool ds = false, split = false, verify = false, hpc = false;
    int mk = 0, mw = 0, ik = 0, iw = 0, bk = 0, bw = 0;
    int a = 1;
    for (; a < argc && argv[a][0] == '-'; ++a) {
        std::string o = argv[a];
        if (o == "-d") ds = true; else if (o == "-U") hpc = true; else if (o == "-s") split = true; else if (o == "-V") verify = true;
        else if (o == "-M" && a + 1 < argc) {
            std::string v = argv[++a]; size_t c = v.find(',');
            if (c == std::string::npos) { std::cerr << "-M needs k,w\n"; return 1; }
            mk = std::stoi(v.substr(0, c)); mw = std::stoi(v.substr(c + 1));
            if (mk < 1 || mk > 3 || mw < mk) { std::cerr << "-M: need 1 <= k <= 3 and w >= k\n"; return 1; }
        }
        else if (o == "-B" && a + 1 < argc) {   // byte-mapped digest for any k <= 4 (map in <out>.map)
            std::string v = argv[++a]; size_t c = v.find(',');
            if (c == std::string::npos) { std::cerr << "-B needs k,w\n"; return 1; }
            bk = std::stoi(v.substr(0, c)); bw = std::stoi(v.substr(c + 1));
            if (bk < 1 || bk > 4 || bw < bk) { std::cerr << "-B: need 1 <= k <= 4 and w >= k\n"; return 1; }
        }
        else if (o == "-I" && a + 1 < argc) {   // integer digest: out.S holds uint32 codes, no datasets
            std::string v = argv[++a]; size_t c = v.find(',');
            if (c == std::string::npos) { std::cerr << "-I needs k,w\n"; return 1; }
            ik = std::stoi(v.substr(0, c)); iw = std::stoi(v.substr(c + 1));
            if (ik < 1 || ik > 12 || iw < ik) { std::cerr << "-I: need 1 <= k <= 12 and w >= k\n"; return 1; }
        }
        else { std::cerr << "unknown option " << o << "\n"; return 1; }
    }
    if (argc - a != 2) { std::cerr << "usage: rz-prep [-d] [-s] [-U] [-M k,w [-V]] [-I k,w] [-B k,w] <list.tsv> <out-prefix>\n"; return 1; }
    if (hpc && (bk || ik || mk)) { std::cerr << "-U cannot be combined with digests\n"; return 1; }
    uint64_t raw = 0, asym = 0;
    if (bk) {   // pass 1: which minimizer codes occur (on either strand); give each a byte
        std::vector<uint64_t> seen(257, 0);
        std::ifstream l1(argv[a]); std::string ln;
        while (std::getline(l1, ln)) {
            if (ln.empty() || ln[0] == '#') continue;
            size_t tab = ln.find('\t'); std::string path = tab == std::string::npos ? ln : ln.substr(tab + 1);
            std::vector<std::pair<std::string, std::string>> gs;
            if (split) read_genomes(path, gs); else gs.push_back({"", read_genome(path)});
            std::vector<uint32_t> c;
            for (auto &gp : gs) {
                rz::digest_codes(gp.second, bk, bw, 1u << (2 * bk), c);
                for (uint32_t x : c) if (x < (1u << (2 * bk))) { seen[x]++; seen[rz::dg_rc_code(x, bk)]++; }
            }
        }
        rz::dg_bytemap &M = rz::bytemap();
        M.k = bk; M.w = bw;
        std::fill(M.code2byte, M.code2byte + 257, 0);
        for (int i = 0; i < 256; ++i) M.comp[i] = (unsigned char)i;
        // Only 252 bytes are free (0..2 are terminators, X separates).  If more minimizers occur, the rarest
        // ones (in reverse-complement pairs, so the digest stays strand-symmetric) share one byte.
        uint32_t K = 1u << (2 * bk);
        std::vector<uint32_t> occ; uint64_t total = 0;
        for (uint32_t x = 0; x < K; ++x) if (seen[x]) { occ.push_back(x); total += seen[x]; }
        std::vector<char> merged(K, 0);
        uint64_t nmerged = 0, merged_occ = 0;
        if (occ.size() > 252) {
            std::vector<uint32_t> byf(occ);
            std::sort(byf.begin(), byf.end(), [&](uint32_t a, uint32_t b) { return seen[a] < seen[b]; });
            for (uint32_t x : byf) {
                if (occ.size() - nmerged + 1 <= 252) break;
                for (uint32_t y : {x, rz::dg_rc_code(x, bk)}) if (!merged[y]) { merged[y] = 1; ++nmerged; merged_occ += seen[y]; }
            }
        }
        int next = 3, used = 0, shared = 0;
        for (uint32_t x : occ) {
            if (merged[x]) continue;
            if (next == 'X') ++next;
            M.code2byte[x] = (unsigned char)next++; ++used;
        }
        if (nmerged) { if (next == 'X') ++next; shared = next; for (uint32_t x : occ) if (merged[x]) M.code2byte[x] = (unsigned char)shared; }
        if (next > 255) { std::cerr << "-B: too many distinct minimizers\n"; return 1; }
        M.code2byte[K] = 'X';
        for (uint32_t x = 0; x < K; ++x) if (M.code2byte[x]) M.comp[M.code2byte[x]] = M.code2byte[rz::dg_rc_code(x, bk)];
        if (nmerged) std::cerr << "-B: " << occ.size() << " distinct minimizers; the " << nmerged << " rarest share byte " << shared
                               << " (" << merged_occ << " of " << total << " occurrences, " << 100.0 * merged_occ / total << "%)\n";
        M.on = true; M.save(std::string(argv[a + 1]) + ".map");
        std::cerr << "byte-mapped digest k=" << bk << " w=" << bw << ": " << used << " distinct minimizers\n";
        mk = bk; mw = bw;      // the -M path below now digests with the map (rz::digest_mapped)
    }
    std::ifstream list(argv[a]);
    std::string out = argv[a + 1];
    FILE *fd = nullptr;
    std::string last_species;
    uint64_t nds = 0;
    std::ofstream dsl;
    if (ds) dsl.open(out + ".datasets");
    FILE *fs = fopen((out + ".S").c_str(), "wb");
    std::ofstream tbl(out + ".tbl");
    std::string line;
    uint64_t pos = 0, idx = 0;
    while (std::getline(list, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::string species = "-", path = line;
        size_t tab = line.find('\t');
        if (tab != std::string::npos) { species = line.substr(0, tab); path = line.substr(tab + 1); }
        std::vector<std::pair<std::string, std::string>> gs;
        if (split) read_genomes(path, gs); else gs.push_back({"", read_genome(path)});
        for (auto &gp : gs) {
        std::string gname = split ? path + "#" + gp.first : path;
        if (gp.second.empty()) { std::cerr << "warning: empty genome " << gname << "\n"; continue; }
        if (hpc) { raw += gp.second.size(); gp.second.erase(std::unique(gp.second.begin(), gp.second.end()), gp.second.end()); }
        std::string rc;
        if (ik) {                               // G's codes, then the reverse complement's, each followed by sep
            uint32_t sep = (uint32_t)(1u << (2 * ik));
            std::vector<uint32_t> cg;
            raw += gp.second.size();
            rz::digest_codes(gp.second, ik, iw, sep, cg);
            if (cg.empty()) continue;
            std::vector<uint32_t> cr(cg.rbegin(), cg.rend());
            for (auto &x : cr) if (x != sep) x = rz::dg_rc_code(x, ik);
            cg.push_back(sep); cr.push_back(sep);
            tbl << idx << "\t" << species << "\t" << gname << "\t" << cg.size() - 1 << "\t" << pos << "\n";
            fwrite(cg.data(), 4, cg.size(), fs); fwrite(cr.data(), 4, cr.size(), fs);
            pos += cg.size() + cr.size(); ++idx;
            continue;
        }
        if (mk) {
            raw += gp.second.size();
            std::string dgst = bk ? rz::digest_mapped(gp.second) : rz::digest(gp.second, mk, mw);
            if (verify) {                       // digest(rc(G)) must equal rc of digest(G)
                std::string grc(gp.second.rbegin(), gp.second.rend());
                for (auto &c : grc) c = comp(c);
                if ((bk ? rz::digest_mapped(grc) : rz::digest(grc, mk, mw)) != rz::rc_string(dgst)) ++asym;
            }
            gp.second.swap(dgst);
            if (gp.second.empty()) { std::cerr << "warning: empty digest " << gname << "\n"; continue; }
            rc = rz::rc_string(gp.second);
        } else {
            rc.assign(gp.second.rbegin(), gp.second.rend());
            for (auto &c : rc) c = comp(c);
        }
        std::string &g = gp.second;
        tbl << idx << "\t" << species << "\t" << gname << "\t" << g.size() << "\t" << pos << "\n";
        fwrite(g.data(), 1, g.size(), fs); fputc('X', fs);
        fwrite(rc.data(), 1, rc.size(), fs); fputc('X', fs);
        if (ds) {
            if (!fd || species != last_species) {
                if (fd) fclose(fd);
                std::string name = out + ".ds" + std::to_string(nds++) + ".S";
                fd = fopen(name.c_str(), "wb");
                dsl << name << "\n";
                last_species = species;
            }
            fwrite(g.data(), 1, g.size(), fd); fputc('X', fd);
            fwrite(rc.data(), 1, rc.size(), fd); fputc('X', fd);
        }
        pos += 2 * (g.size() + 1);
        ++idx;
        }
    }
    fclose(fs);
    if (fd) fclose(fd);
    std::cerr << "genomes=" << idx << " |S|=" << pos;
    if (ds) std::cerr << " datasets=" << nds;
    if (ik) std::cerr << " integer digest k=" << ik << " w=" << iw << ": " << raw << " bases -> " << pos << " codes (uint32, both strands)";
    if (mk) std::cerr << " digest k=" << mk << " w=" << mw << ": " << raw << " bases -> " << pos / 2 << " symbols per strand ("
                      << (double)pos / 2 / (raw ? raw : 1) << ")";
    if (hpc) std::cerr << " homopolymer-compressed: " << raw << " bases -> " << pos / 2 << " per strand (" << (double)pos / 2 / (raw ? raw : 1) << ")";
    if (verify) std::cerr << " strand-asymmetric genomes=" << asym;
    std::cerr << "\n";
    return verify && asym ? 3 : 0;
}
