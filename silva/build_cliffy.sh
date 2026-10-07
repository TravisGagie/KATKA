#!/bin/bash
# Builds Cliffy (pfp_doc / pfp_doc64) in silva/cliffy (clone it first; see the README).  Log: silva/build_cliffy.log
# Compiles through silva/gxx-compat, which adds fixes needed for Cliffy's dependencies under recent GCC when the
# compiler supports them.  JOBS (default 8) parallel compile jobs.  On failure, prints the first errors.
S=$(cd "$(dirname "$0")" && pwd); LOG=$S/build_cliffy.log
[ -d $S/cliffy ] || { echo "FAILED: no silva/cliffy (git clone https://github.com/oma219/cliffy.git silva/cliffy)"; exit 1; }
for t in g++ cmake make; do command -v $t > /dev/null || { echo "FAILED: missing $t"; exit 1; }; done
export CMAKE_POLICY_VERSION_MINIMUM=3.5      # its dependencies ask for CMake versions older than 3.5
F="-include cstdint"                          # headers that forget <cstdint>
# sdsl-lite template bodies that no longer compile under GCC 15 but are never used (only if g++ knows the flag)
echo 'int main() { return 0; }' | g++ -Werror -Wno-template-body -x c++ - -o /dev/null 2> /dev/null && F="$F -Wno-template-body"
export GXX_COMPAT_FLAGS="$F"
{ echo "g++: $(g++ --version | head -1)"; echo "cmake: $(cmake --version | head -1)"; echo "extra flags: $F"
  echo "cliffy: $(git -C $S/cliffy rev-parse --short HEAD 2> /dev/null)"; } | tee $LOG
cd $S/cliffy && rm -rf build && mkdir -p build && cd build
{ cmake -DCMAKE_CXX_COMPILER="$S/gxx-compat" .. && make -j${JOBS:-8} install; } >> $LOG 2>&1
st=$?
if [ $st = 0 ] && [ -x $S/cliffy/build/cliffy ]; then echo "Cliffy built: silva/cliffy/build/cliffy"; exit 0; fi
echo "FAILED to build Cliffy (exit status $st).  First errors in silva/build_cliffy.log:"
grep -n -m 12 -E "error:|Error [0-9]|CMake Error|fatal|No such file|not found" $LOG
echo "--- last lines of the log:"; tail -15 $LOG
exit 1
