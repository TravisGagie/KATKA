# Grammar-compressing the run-length tag array

**Status (October 2026): implemented.** `rz-index/gtag.hpp` (encoding), `rz-index/rz_gtagbuild.cpp` (builder from
RePair's .R/.C output, checking every run), `.gtag` files load in `rz-classify -T` for counting (`RZ_COUNTS=1`), and
`silva/run_gtag.sh` builds, checks and times it.  On SILVA, RePair on the whole genus sequence (Navarro's `irepair`
from BigRePair's `largeb_repair`, 4.5 min, 8.9 GB): 19.0 M rules, final sequence 38.4 M symbols, height 46; the
`.gtag` file is 0.459 GB (run starts 0.228 + rules 0.089 + final sequence 0.136 + tables) instead of 0.857 GB; the
default index 1.04 GB instead of 1.44 GB, 75.0 us/read instead of 66.0 (L = 30), identical answers; on the Pareto
surface.  BigRePair's grammar (57.1 M rules) gives 0.530 GB and 76.7 us.  The notes below are the original plan.


Goal: an optional second level of compression for KATKA's tag array, to add to the paper if the result is on the
time-space-accuracy Pareto surface (paper: "Taxonomic Classification with Complete Tag Arrays", Section 3.3 and
Appendix B.3; `pareto_ft.py --counts`).  Framing for the paper: the tag array is **run-length compressible** (that is
the main idea and stays the default); because LF maps a BWT run's rows to consecutive rows whose suffixes are one
character longer, the run-length encoded tag array is also made of copies of its own substrings (like the gap-coded
suffix array), so it can **optionally** be grammar-compressed on top.

## What the index stores now (rz-index/tag.hpp)
* run starts of the tag array in an Elias-Fano bitvector `RB` (rank `RBr`, select `RBs`): 0.228 GB on SILVA;
* the genus of every run, `L` (int_vector, 14 bits): 0.629 GB;
* the RMQ for listing only, in `bac.s1.tag.rmq` (not loaded when counting).
Counting a BWT interval [sp, ep): runs a = run_of(sp) .. b = run_of(ep-1); for each run, its genus `tag(i)` and its
overlap with [sp, ep) from `RBs`.  Run lengths therefore never need to be in the grammar: **only the genus sequence of
the runs is grammar-compressed**, and expansion lengths are counted in runs.

## Measurements so far (SILVA, undigested, 359,351,512 runs, 9,118 genera; desktop i7-8700)
* LZ77 (tagz.cpp): genus sequence z = 42.3 M phrases; run lengths 13.0 M; both interleaved 52.2 M (BWT r = 68.3 M).
* BigRePair -i on the genus sequence (tokwrite.cpp genus; run_bigrepair.sh): 57.1 M rules, |C| = 2.31 M, height 40,
  3,654 distinct expansion lengths; C symbols expand to 156 runs on average (max 8,716) (glen.c).
* Size estimates (gtrick.c): plain rules 2 x 26 bits + C = 0.379 GB, +29-bit lengths = 0.586 GB; SPIRE'20 encoding
  (Gagie, I, Manzini, Navarro, Sakamoto, Seelbach Benkner, Takabatake, "Practical random access to SLP-compressed
  texts", SPIRE 2020: non-terminals grouped by expansion length, rule = (left child's length, offsets of both children
  in their groups), MPHF from length to group) = 0.380 GB with all lengths; storing the left child's length in
  ceil(log2(l-1)) bits within the group of length l = **0.319 GB**.  Total tag array: 0.228 + 0.319 = 0.55 GB vs
  0.857 GB now; default index 1.44 GB -> about 1.13 GB.
* Cutting the grammar at expansion length w and storing the cut symbols' expansions (gcut.c): 0.60-0.70 GB, not
  worth it (the dictionary loses the sharing).
* Speed (rz-gbench, plain arrays, all lengths stored; desktop's sandboxed VM, 40,000 reads, 213 k intervals, 16.3
  runs per interval on average): counting 12.7 us per read with the tag array, 20.9 us with the grammar, identical
  answers.  Whole classification there: 83 us per read (66 us on the desktop itself).

## Next steps
1. Implement the compact grammar as a counting backend (`tag.hpp` or a new `gtag.hpp`): run starts as now; rules in
   the SPIRE'20 grouped encoding with group-relative left lengths (fixed width per group, bit-packed; a simple
   length->group map is fine instead of an MPHF since there are only ~3.7 k distinct lengths); C with cumulative run
   counts (Elias-Fano).  Decode a..b: predecessor in C, descend (height <= 40) using lengths, then sequential decode
   with a stack.  Optional per-group speedups only if needed.
2. `rz-tagbuild` (or a new tool) builds it from bigrepair's `.R`/`.C`; `rz-classify -T` loads either kind (magic).
3. Check: identical outputs to the run-length tag array (rz-gbench does this for plain arrays; extend to the compact
   encoding), on the mini reference first, then SILVA.
4. Time on the desktop with the app closed (`run_pareto_ft.sh` / `run_katka.sh` style, COUNTS=1), sizes from the
   files; rerun `pareto_ft.py --counts` with the new configuration.
5. If on the Pareto surface: paper Section 3.3 (one paragraph after the run-length compression), tables
   tab:choices / tab:pareto, abstract/conclusions sentence on grammar compression (now future work) becomes a result.

## Testing without the full SILVA build
demo_mini.sh builds a 304-genus subset (every 30th genus) in silva/work/mini/ref.  To test the grammar code, build a
KATKA index on that subset (rz-prep / Big-BWT / rz-build -G / rz-csabuild / rz-tagbuild, as in silva/build_katka.sh,
pointed at the subset), write its tokens with tokwrite, run bigrepair -i (https://gitlab.com/manzai/bigrepair; `-i`
reads 32-bit integers < 2^30), and compare outputs.  Full-size timing has to be on the desktop.
