#!/bin/bash
# Builds Cliffy (pfp_doc / pfp_doc64) in silva/cliffy (clone it first; see the README).  Log: silva/build_cliffy.log
# Compiles through silva/gxx-compat, which adds fixes needed for Cliffy's dependencies under recent GCC when the
# compiler supports them.  JOBS (default 8) parallel compile jobs.  On failure, prints the first errors.
S=$(cd "$(dirname "$0")" && pwd); LOG=$S/build_cliffy.log
[ -d $S/cliffy ] || { echo "FAILED: no silva/cliffy (git clone https://github.com/oma219/cliffy.git silva/cliffy)"; exit 1; }
for t in g++ cmake make; do command -v $t > /dev/null || { echo "FAILED: missing $t"; exit 1; }; done
# Cliffy needs CMake 3.13 or later; a newer one can be installed without root: python3 -m pip install --user cmake
cv=$(cmake --version | head -1 | grep -o "[0-9][0-9.]*" | head -1)
if [ "$(printf '%s\n3.13\n' "$cv" | sort -V | head -1)" != 3.13 ]; then
  echo "FAILED: Cliffy needs CMake 3.13 or later, and this is CMake $cv."
  echo "  Install a newer one without root:  python3 -m pip install --user cmake  (then put ~/.local/bin first in PATH)"
  exit 1
fi
# silva/cliffy_rz.patch fixes Cliffy for recent CMake (thirdparty/CMakeLists.txt defines libdivsufsort's targets
# again after sdsl-lite has) and older glibc (WEXITSTATUS of an rvalue in src/pfp_doc.cpp).  The patched files
# are restored to Cliffy's version and the whole patch reapplied, so an older version of the patch is replaced.
git -C $S/cliffy apply --reverse --check $S/cliffy_rz.patch 2> /dev/null || {
  git -C $S/cliffy checkout -- thirdparty/CMakeLists.txt src/pfp_doc.cpp && git -C $S/cliffy apply $S/cliffy_rz.patch
} || { echo "FAILED: could not apply silva/cliffy_rz.patch to silva/cliffy"; exit 1; }
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
