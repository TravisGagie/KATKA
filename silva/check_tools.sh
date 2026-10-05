#!/bin/bash
# Records which tools are installed on this machine, for the Cliffy/SILVA experiment.
cd "$(dirname "$0")"
{
for t in seqtk art_illumina kraken2 kraken2-build snakemake python3 g++ cmake git wget; do
  printf "%-14s " "$t:"; command -v $t || echo MISSING
done
echo "--- versions"
kraken2 --version 2>&1 | head -1
snakemake --version 2>&1 | head -1
art_illumina 2>&1 | grep -i -m1 "version"
seqtk 2>&1 | grep -i -m1 version
python3 -c "import regex, requests; print('python regex', regex.__version__, 'requests', requests.__version__)" 2>&1
g++ --version | head -1
nproc; free -g | head -2; df -h ~ | tail -1
} > tools.txt 2>&1
cat tools.txt
