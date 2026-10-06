# KATKA

KATKA classifies DNA reads taxonomically with long maximal exact matches (MEMs) and compressed indexes.
This repository has the code and the experiments for

> Travis Gagie and Gonzalo Navarro. *Taxonomic Classification with Tag Arrays.* 2026.

By default, KATKA works in three steps:

* It finds the MEMs of at least *L* = 30 bases between a read and the reference with **Boyer–Moore–Li**, on a
  run-length compressed suffix array (RLCSA) with a lookup table of all 10-mers.
* It lists the genera containing each MEM with a **run-length compressed tag array** and Muthukrishnan/Sadakane
  document listing, counting the occurrences in each genus.
* It gives each genus **credit in proportion to its occurrences**, summed over the MEMs of both mates.

Two options trade accuracy or memory for speed:

* **Minimizer digests:** matching over strand-symmetric digests is about three times smaller and faster, but
  slightly less accurate.
* **Two-level parse index:** with closed syncmers, slightly faster but uses more memory.

The name comes from two earlier projects, KATKA (SPIRE 2022) and its MEM-based successor (SEA 2024). Almost
everything has changed since then.

The experiments use the SILVA 16S rRNA benchmark of Cliffy (Ahmed, Boucher and Langmead, Genome Research 2025),
with its scripts and read simulator. SILVA (https://www.arb-silva.de) is licensed under CC BY 4.0. The simulated
reads in `silva/data/` and in the `silva-sample` release were made with MicrobeMixer from SILVA 138.1.

## Requirements

You need Linux with g++ (C++17), cmake, git, wget, zlib, seqtk and GNU time. For Python 3 you need `regex`,
`requests`, `aiohttp`, `pandas` and `multiprocess`. On Ubuntu:

    sudo apt install build-essential cmake git wget zlib1g-dev seqtk time python3-pip
    pip install regex requests aiohttp pandas multiprocess      # add --break-system-packages if pip asks

## Quick demo

    git clone https://github.com/TravisGagie/KATKA.git && cd KATKA
    ./demo_quick.sh

The script works in three steps:

1. It fetches and builds the dependencies and downloads SILVA 138.1.
2. It builds the reference, with one document per genus, and the default index.
3. It classifies 5,000 simulated read pairs for each 16S region (V1–V2, V3–V4, V4 and V4–V5), first with one
   thread and then with all cores.

It prints the index size, the time per read, the throughput and the genus-level accuracy, and writes them to
`silva/work/compare/katka-<host>/summary.txt`. It needs about 20 GB of disk and 16 GB of memory, mostly for
building the index. It can be resumed.

On a laptop CPU with performance and efficiency cores, run `PIN=0 ./demo_quick.sh` so the single-thread timing
stays on one core.

`silva/build_katka.sh` builds only what you ask for: by default, just the default configuration's index. Its
first lines list the options, such as `INDEXES="rz rzdg"` for digests too, `TAGS`, `CSA`, `SR` and `GRIDS` for the
other structures in the paper, and `ALL=1` for everything.

## Full demo: the paper's experiments

    ./demo_full.sh

The script works in three steps:

1. It builds every index in the paper.
2. It downloads the paper's reads from the `silva-sample` release (258 MB). These are every 50th simulated
   pair, 200,000 per region.
3. It runs every experiment on them, in the order listed at the top of the script, logging to
   `silva/work/logs/`.

It needs about 60 GB of disk and 16 GB of memory and takes a day or more, so run it on an otherwise idle
machine. It can be resumed, and `STEPS="..."` runs only some of the experiments. The time/space/accuracy table
and its Pareto set end up in `silva/work/pareto_counts.txt`.

## Comparing with Cliffy on the same machine

`silva/run_cliffy.sh` and `silva/run_katka.sh` run Cliffy and KATKA on the paper's reads, on one machine:

    ./setup_deps.sh && make -C rz-index DEPS=$PWD/deps
    INDEXES="rz rzdg" silva/build_katka.sh   # reference and our default indexes, with and without digests
    silva/get_reads.sh                    # the paper's reads
    git clone https://github.com/oma219/cliffy.git silva/cliffy && git -C silva/cliffy checkout 3763529
    silva/build_cliffy.sh                 # check silva/build_cliffy.log
    silva/run_cliffy.sh                   # builds Cliffy's indexes if they are missing, then runs Cliffy
    silva/run_katka.sh                    # then KATKA, with the same number of threads

Cliffy's index construction writes over a terabyte of temporary files for the undigested index (its paper's
Table 1). Set `TMPDIR_CLIFFY` to a directory on a large disk.

Each script does three things:

1. **Single thread.** Each configuration classifies the reads once, then classifies a single read. The
   difference gives the time per read without loading the index. `PIN=c` runs these on logical CPU *c*.
2. **Parallel.** Each configuration classifies the same reads with P threads.
   * KATKA uses `rz-classify -j P`, which shares one copy of the index.
   * Cliffy queries with one thread, so it gets P processes, each on one slice of the reads.
   * For Cliffy, P is the number of logical cores, reduced if P copies of its largest index would not fit in
     memory. `run_katka.sh` then uses the same P. `P=...` overrides it.
   * Every parallel output must be identical to the single-thread one.
3. **Accuracy** at the genus level.
   * KATKA is scored with proportional credit.
   * Cliffy is scored with its LCA query and, unless `THEIRS=0`, with its own `classify_silva_readset.py`. That
     covers all its methods, including approximate listing, but is slow.

By default, Cliffy's `full_text` and `minimizers` indexes are compared with KATKA at L = 15 and 30, with and
without digests. Change them with `CLIFFY=...` and `OURS=...`. The results go to
`silva/work/compare/cliffy-<host>/` and `silva/work/compare/katka-<host>/`.

## Classifying your own reads

    cd rz-index
    RZ_COUNTS=1 ./rz-classify -j 8 -L 30 -l -C idx/bac.csa -T idx/bac.s1.tag -F 10 idx/bac.rz reads.fq out.listings

For every MEM of at least 30 bases, this lists the genera containing it and its number of occurrences in each,
using 8 threads. The output is the same for any number of threads.

The output has one line per read in Cliffy's listing format: `[start,end] {doc:count,doc:count,...}` for each
MEM. Without `RZ_COUNTS=1`, it lists just the documents.

* `silva/build_katka.sh` shows how to build the index from FASTA files.
* `silva/rz_score2.py` shows how to turn listings into classifications.
* For a digested index, add `-B idx/bac.map` and use `-F 2`.
* Run `rz-classify` without arguments to see all the options.

## Layout

    demo_quick.sh, demo_full.sh   the two demos
    setup_deps.sh       fetches the dependencies (pinned commits) and SILVA 138.1
    runrz.sh            builds the text, BWT, parses, rz-index, r-index and sr-indexes from a list of FASTA files
    rz-index/           C++17 sources: index construction (rz-prep, rz-lz77, rz-build, rz-aux, sr-build,
                        rz-csabuild, rz-tagbuild, rz-parsebuild, rz-vfybuild) and the classifier (rz-classify)
    silva/build_katka.sh, get_reads.sh, reads.sh   reference, indexes and reads for the experiments
    silva/run_*.sh      the experiments (each says at the top what it measures and where its results go)
    silva/rz_score2.py, pareto_ft.py   scoring and the time/space/accuracy table
    silva/run_silva.sh  Cliffy's complete pipeline, including simulating all 10 million read pairs per region
    silva/make_samples.sh   how the read samples were taken from the full simulation
    silva/data/         the quick demo's reads

## License

GPL-3.0; see `LICENSE`.
