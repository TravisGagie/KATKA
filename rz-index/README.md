# rz-index: leftmost and rightmost occurrences in genome collections

The rz-index answers one query. Given a pattern P over {A,C,G,T}, it reports the starting positions of the **leftmost and rightmost** occurrences of P in

    S = G_1 X rc(G_1) X G_2 X rc(G_2) X ... G_k X rc(G_k) X

Here the genomes G_j are listed in the left-to-right order of the leaves of a phylogenetic tree. Mapping the two positions to genomes and taking their LCA gives the smallest subtree containing every occurrence of P or rc(P).

## Structure

* **RLBWT(S$).** This is the same RLBWT an r-index for S uses. Since S is closed under reverse complement, one RLBWT serves both searches: backward search for P, and for rc(P) in place of P^rev.
* **Left side.** It is built from a left-referencing parse of S, where every phrase copies from an earlier source; greedy LZ77 is one such parse, and a windowed version is another. The leftmost occurrence of P always contains a phrase end e. For each phrase end there is one grid point:
  * y is the rank of the BWT row of S[e+1..];
  * x is the rank of the row of mirror(e), the position in the partner strand at which rc(S[..e]) begins;
  * the value is the phrase's rank in text order.

  The side also has three sparse bitvectors (Elias–Fano) and a wavelet-matrix 2D range-min with a succinct RMQ per level.
* **Right side.** It is the mirror image, built from a left-referencing parse of S^rev. The rightmost occurrence contains a phrase start b of that parse (read in S), and the 2D range query is a max.
* **Genome starts.** A sparse bitvector marks where each strand starts, so a position maps to its genome.

## Query

1. Backward search for P gives the intervals of every P[i..m).
2. Backward search for rc(P) gives the intervals of every rc(P[0..i)).
3. For each split, map the two intervals through the sparse bitvectors and run the 2D range-min (left side) or range-max (right side).

`rz-bench` reports the total time and a breakdown into backward searches, left grid, and right grid. For the r-index it reports the total and a breakdown into backward search with toehold, and the φ loop that tracks the minimum and maximum.

## Building

```
mkdir deps && cd deps
git clone https://github.com/simongog/sdsl-lite && (cd sdsl-lite && ./install.sh "$PWD/../sdsl")
git clone https://github.com/alshai/Big-BWT && (cd Big-BWT && make)     # or any BWT(S$) builder
cd ../rz-index && make        # needs zlib; builds rz-prep rz-lz77 rz-build rz-genpat rz-bench
# make test-tools builds mbwt (a small multi-string BWT builder for tests) and ri-build-from (needs alshai/r-index)
```

## Quick start

```
BIGBWT=path/to/Big-BWT/bigbwt PFPMERGE=path/to/pfp-merge THREADS=8 ./rz-pipeline.sh list.tsv out/zymo
```

This runs every step below: `rz-prep -d`, Big-BWT per species dataset, `PFPmerge.py`, both parses, `rz-build -a`, and a benchmark. pfp-merge needs `git submodule update --init --recursive` and its `build/` directory made. Its output is exactly the multi-string BWT that `rz-build` expects, with one 0x00 per dataset. On both test collections it was byte-identical to the output of the test helper `mbwt`.

## Scalable construction (no suffix array of S anywhere)

| step | tool | memory | notes |
|---|---|---|---|
| 1. S from FASTA files | `rz-prep [-d] [-s] list.tsv out` | one genome | `list.tsv`: `species<TAB>fasta[.gz]` per line, in tree order. Contigs are joined with X, and non-ACGT characters become X. Writes `out.S` and `out.tbl`. With `-s`, a FASTA file may hold many genomes, split by sample ID (header text before the first `.`, e.g. AllTheBacteria's `SAMEA1410869.contig00001`). With `-d` it also writes one dataset file per species (`out.ds<j>.S`, listed in `out.datasets`) for pfp-merge. |
| 2. BWT | Big-BWT on each file in `out.datasets`, then `PFPmerge.py -o out.bwt <datasets>` (run from the pfp-merge directory); or Big-BWT on `out.S`; or any tool | tool-dependent | Either a **multi-string BWT/eBWT with one terminator per species dataset** (pfp-merge), or **BWT(S$) of S as one string** (Big-BWT). One byte per row; terminator bytes 0x00–0x02. |
| 3. parses | `rz-lz77 [-t T] [-w W] [-b B] out.S out.left` and `rz-lz77 -r ... out.S out.right` | ≈ 13(W+B) bytes per thread | Windowed greedy LZ77 (KKP), in independent fixed blocks, parallel. S^rev is read backwards and never written. |
| 4. indexes | `rz-build [-a] out.S out.bwt out.left out.right out.tbl out`; `rz-build -A out.S out.bwt - - out.tbl out` builds only the r-index | ≈ 21r + 25z bytes plus the RLBWT (rz-index); ≈ 30r plus the RLBWT (`-A`) | Streams the BWT twice to build the RLBWT. It then counts the terminators to decide between the two BWT cases. In the multi-string case it matches each terminator row to its dataset by spelling the dataset's tail backwards with LF. Next it does one LF pass over each string, using move-style run arrays, to find the BWT rows of the O(z) positions needed. Writes `out.rz`. With `-a`, the same pass collects the SA samples and it also writes the baseline r-index `out.rix` on the same RLBWT. |
| 5. benchmark | `rz-genpat out.S m N pats` then `rz-bench [-t out.S] out.rz out.rix pats` | index sizes | `-t` adds a brute-force check (find/rfind) for small S. |

Correctness only needs every phrase to have an earlier source. A windowed parse is therefore valid, and it is barely larger when the genomes of a species are consecutive.

**Baseline r-index (`rz::rindex`).** It is the standard Gagie–Navarro–Prezza algorithm:
- a backward search with the toehold (SA at run ends);
- then φ from predecessor samples at run starts, tracking the minimum and maximum.

It is re-implemented here on the same RLBWT class so that it also works on multi-string BWTs. There, the rows whose BWT character is a terminator are sampled too, so φ never runs across a string boundary. On single-string BWTs it gave the same answers as Prezza's r-index and within 5% of its size (90.4 MB vs 86.3 MB on the 300 MB test). `ri-build-from` builds Prezza's r-index from Big-BWT-style samples, for cross-checking only.

## Test results (2-core cloud VM, 7 GB RAM)

The test data are synthetic species-structured collections from `test/gen.py`: species with independent random base genomes, and genomes derived along a random within-species tree (0.05% SNPs and a few indels per generation). The reverse complements are interleaved by `rz-prep`.

* **Correctness.** Leftmost and rightmost positions matched brute force (find/rfind) and the r-index in every run, with zero disagreements. This covered a 3.8 MB collection with blocks as small as 50 K (m = 1 to 300) and a 300 MB collection (m = 10, 32, 100). On the 2.2 GB collection, the answers matched the r-index for all 5000 patterns.
* **Multi-string BWTs.** pfp-merge's actual output (from the zip you sent) on the 3.8 MB and 300 MB collections was byte-identical to that of `mbwt`, which builds a multi-string BWT with one terminator per species dataset. The rz-index built from the pfp-merge output of the 3.8 MB collection matched brute force for m = 1 to 300. On the 3.8 MB and 300 MB collections, the rz-index built from it gave the same answers as the one built from Big-BWT's single-string BWT, and both matched brute force. A small collection whose BWT has two adjacent terminators, the case where a terminator row doesn't start a run, was also correct at m = 1 to 1000.
* **Windowed parse vs exact greedy LZ77** (300 MB collection; one genome pair ≈ 2 MB):

  | window W = block B | exact | 64 M | 16 M | 4 M |
  |---|---|---|---|---|
  | z | 1,152,592 | 1,163,561 (+1%) | 1,219,014 (+6%) | 1,498,044 (+30%) |
  | peak memory (1 thread) | 3.6 GB | 2.1 GB | 0.54 GB | 0.14 GB |

* **Scaling run: N = 2.2 G** (5 species × 110 genomes × 2 Mbp, both strands; r = 25.5 M, zL = 4.19 M, zR = 4.68 M):

  | step | time | peak memory |
  |---|---|---|
  | rz-prep | 25 s | 0.01 GB |
  | Big-BWT | 262 s | 1.75 GB |
  | rz-lz77, left (2 threads, W = B = 64 M) | 466 s | 3.3 GB |
  | rz-lz77, right | 447 s | 3.1 GB |
  | rz-build -a (LF pass at 3.5 M steps/s) | 667 s | 1.65 GB |
  | ri-build-from | 64 s | 1.2 GB |

  Sizes: the rz-index is 255 MB (RLBWT 74 MB, left side 84 MB, right side 97 MB). The r-index is 272 MB.

  Queries, µs per pattern (1000 random patterns each):

  | m | rz total (bs / left / right) | r-index total (bs / φ) | r-index avg occ |
  |---|---|---|---|
  | 10 | 174 (29 / 72 / 72) | 1423 (47 / 1376) | 2224 |
  | 20 | 120 (53 / 31 / 36) | 155 (63 / 92) | 105 |
  | 50 | 227 (117 / 56 / 54) | 179 (101 / 77) | 97 |
  | 100 | 427 (231 / 92 / 103) | 262 (186 / 76) | 88 |
  | 200 | 719 (410 / 152 / 156) | 356 (293 / 63) | 71 |

  With 110 genomes per species, a species-specific pattern occurs about 110 times. That is where the two indexes cross over. The r-index's φ cost grows linearly with occ, while the rz-index's cost doesn't depend on occ. With hundreds to thousands of genomes per species, as in AllTheBacteria, the balance moves toward the rz-index.
