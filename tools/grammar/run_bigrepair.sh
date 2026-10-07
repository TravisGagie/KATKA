#!/bin/bash
# BigRePair (-i, 32-bit symbols) on the run-length encoded tag array: genus IDs only, then lengths and genera
# interleaved (tags_*.int32, written by tokwrite).  Results and timings in bigrepair.out.
#   cd tools/grammar && ./run_bigrepair.sh
cd "$(dirname "$0")"; BR=${BIGREPAIR:-../../deps/bigrepair/bigrepair}
for f in tags_genus.int32 tags_both.int32; do
  echo "=== $f"
  /usr/bin/time -f "  (%e s, peak %M KB)" $BR -i -r $f 2>&1 | grep -v "^$"
  echo "  .R $(stat -c %s $f.R) bytes, .C $(stat -c %s $f.C) bytes"
done 2>&1 | tee bigrepair.out
