# KATKA

KATKA classifies DNA reads taxonomically using long maximal exact matches (MEMs) and compressed
indexes. This repository has the code and experiment scripts for

> Travis Gagie and Gonzalo Navarro. *Taxonomic Classification with Long Maximal Exact Matches,
> Run-Length Compressed Indexes and Tag Arrays.* 2026.

The name comes from two earlier projects: KATKA (SPIRE 2022), a Kraken-like index with *k* given at query
time, and its MEM-based successor (SEA 2024). Almost everything has changed since then. The present
version works in three steps:

* It finds the MEMs of length at least *L* in a read with **Boyer–Moore–Li**, over the sequences or
  over strand-symmetric **minimizer digests** of them.
* It answers each MEM either by **listing** the genera that contain it, or by their **LCA**.
* It answers with a subsampled r-index, a run-length compressed suffix array (RLCSA), or a
  **run-length compressed tag array** with Muthukrishnan/Sadakane document listing. The *rz-index*
  finds only the leftmost and rightmost occurrences.

The experiments reproduce the SILVA 16S rRNA benchmark of Cliffy (Ahmed, Boucher and Langmead, Genome
Research 2025), using its scripts and read simulator.

## Layout

    rz-index/   C++17 sources: index construction (rz-prep, rz-lz77, rz-build, rz-aux, sr-build,
                rz-csabuild, rz-tagbuild, rz-vfybuild) and the classifier (rz-classify)
    runrz.sh    builds the text, BWT, parses, rz-index, r-index and sr-indexes from a list of FASTA files
    silva/      SILVA experiment scripts, scoring (rz_score.py, rz_score2.py) and the Pareto table
    setup_deps.sh  fetches the dependencies (pinned commits) and SILVA 138.1

## Building

    ./setup_deps.sh                 # sdsl-lite, Big-BWT, pfp-merge, Cliffy's scripts, MicrobeMixer, SILVA
    make -C rz-index DEPS=$PWD/deps

You also need g++ with C++17 support, cmake, zlib, seqtk, ART (`art_illumina`) and GNU time.
For Python 3 you need `regex`, `requests`, `aiohttp`, `pandas` and `multiprocess`.
If a SILVA download in `setup_deps.sh` fails, get the three files from the SILVA 138.1 release on
https://www.arb-silva.de and put them, uncompressed, in `silva/exp1_data/`.

## Reproducing the SILVA experiments

Every script can be resumed: it skips any step whose output already exists. Results go under `silva/work/`.

    silva/run_silva.sh        # reference (one document per genus, in tree order), 10M simulated read pairs
                              # per region, undigested index, and Cliffy-style classification with
                              # Ziv-Merhav phrases
    silva/run_silva_dg.sh     # digested reference (k = 4, w = 11) and its indexes; ZMPs vs MEMs
    silva/run_silva_mem.sh    # undigested MEMs (forward-backward)
    silva/build_indexes.sh    # RLCSAs, sr-indexes (s = 0, 4, 16), tag arrays (s = 1, 4, 16), etc.
    silva/run_bml.sh          # BML over digests: accuracy for L = 15..100 (LCA and listing), and speed
    silva/run_bml_ud.sh       # the same without digestion
    silva/run_tag.sh          # speed of tag-array listing and LCA
    python3 silva/pareto_bml.py dg ud   # time/space/accuracy table and Pareto set

Negative results from the paper:

    silva/run_hyb.sh          # answering MEMs with many occurrences by their LCA (rz_classify -H / -A)
    silva/run_vfy_acc.sh      # verifying digested MEMs against the DNA (rz-classify -V)
    silva/run_kebab.sh        # ideal KeBaB pseudo-MEMs before undigested BML (RZ_KEBAB=k)
    silva/run_trim.sh         # what a one-level phrase index would lose (RZ_TRIM, RZ_TRIMFIX)

`silva/test_bml.sh` and `silva/test_mem.sh` are quick consistency checks. They confirm that all backends give
identical answers, and `RZ_BML_CHECK=1` compares BML against brute force.

The paper's appendix gives every parameter of the experiments. In particular, timings are
single-threaded and exclude index loading, and Cliffy's own indexes are not built by default
because their construction needs over a terabyte of temporary disk space.

## Classifying your own reads

    ./rz-index/rz-classify -L 30 -l -C idx/bac.csa -T idx/bac.s1.tag idx/bac.rz reads.fq out.listings

This lists the genera of every MEM of at least 30 bases, using an explicit RLCSA and an unsampled tag array.
For a digested index, add `-B idx/bac.map`. Run `rz-classify` without arguments to see all the options.
The output uses Cliffy's listing format, one line per read: `[start,end] {doc,doc,...}` for each MEM.

## License

KATKA is free software, distributed under the GNU General Public License, version 3 (see `LICENSE`).
