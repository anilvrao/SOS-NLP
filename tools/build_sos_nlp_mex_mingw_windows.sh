#!/usr/bin/env bash
set -euo pipefail

export PATH=/mingw64/bin:/usr/bin:$PATH

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATLAB_ROOT="${MATLAB_ROOT:-/c/Program Files/MATLAB/R2026a}"
MATLAB_RELEASE="${MATLAB_RELEASE:-R2026a}"

SRC="$ROOT/src/soscode"
BRIDGE="$ROOT/src/highs_bridge"
HIGHS_SRC="$ROOT/third_party/HiGHS"
HIGHS_BUILD="$ROOT/build/highs-build-win64"
HIGHS_PREFIX="$ROOT/third_party/highs-install-win64"
HIGHS_LIB="$HIGHS_PREFIX/lib/libhighs.a"
IOSRC="$ROOT/src/fortran/sos_mex_io.f"
BASE_MEXSRC="$ROOT/src/mex/sos_nlp_mex_clean.cpp"
WIN_MEXSRC="$ROOT/build/generated_mex_src_win64/sos_nlp_mex_mingw_windows.cpp"
OBJ="$ROOT/build/obj_mingw_windows"
SOSLIB="$ROOT/build/libsos_nlp_mingw_windows.a"
OUTDIR="$ROOT/build/mex/mingw_windows/mexw64"

mkdir -p "$HIGHS_BUILD" "$OBJ" "$OUTDIR" "$(dirname "$WIN_MEXSRC")"

for f in "$SRC" "$BRIDGE" "$HIGHS_SRC/CMakeLists.txt" "$IOSRC" "$BASE_MEXSRC" \
         "$MATLAB_ROOT/extern/include/mex.h" \
         "$MATLAB_ROOT/extern/lib/win64/mingw64/mexFunction.def" \
         "$MATLAB_ROOT/extern/lib/win64/mingw64/libmex.lib" \
         "$MATLAB_ROOT/extern/lib/win64/mingw64/libmx.lib" \
         "$MATLAB_ROOT/extern/lib/win64/mingw64/libmat.lib"; do
  if [[ ! -e "$f" ]]; then
    echo "Missing required path: $f" >&2
    exit 2
  fi
done

cmake -S "$HIGHS_SRC" -B "$HIGHS_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HIGHS_PREFIX" \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DBUILD_CXX=ON \
  -DBUILD_CXX_EXE=OFF \
  -DBUILD_TESTING=OFF \
  -DBUILD_EXAMPLES=OFF \
  -DZLIB=OFF \
  -DFORTRAN=OFF \
  -DCSHARP=OFF \
  -DPYTHON_BUILD_SETUP=OFF

cmake --build "$HIGHS_BUILD" --target install --config Release --parallel "${BUILD_JOBS:-8}"

if [[ ! -f "$HIGHS_LIB" ]]; then
  echo "Missing source-built HiGHS archive: $HIGHS_LIB" >&2
  exit 2
fi

rm -f "$OBJ"/*.o "$SOSLIB" "$OUTDIR/sos_nlp_mex.mexw64"

sed \
  -e 's|mxCreateString("maca64")|mxCreateString("win64")|g' \
  -e "s|mxCreateString(\"R2025b\")|mxCreateString(\"$MATLAB_RELEASE\")|g" \
  -e 's|gfortran + HiGHS|gfortran MinGW Windows + source-built HiGHS|g' \
  -e 's|sos_nlp_highs_experiment/build/libsos_nlp.a + HiGHS|build/libsos_nlp_mingw_windows.a + source-built HiGHS|g' \
  "$BASE_MEXSRC" > "$WIN_MEXSRC"

count=0
while IFS= read -r -d '' f; do
  b="$(basename "$f")"
  case "$b" in
    qpopt.f|qpcore.f|lpcore.f) continue ;;
  esac
  o="$OBJ/${b%.f}.o"
  gfortran -O2 -std=legacy -fallow-argument-mismatch \
    -ffunction-sections -fdata-sections \
    -c "$f" -o "$o"
  ar rcs "$SOSLIB" "$o"
  count=$((count+1))
done < <(find "$SRC" -maxdepth 1 -name '*.f' -print0 | sort -z)

for f in "$BRIDGE"/*.c; do
  b="$(basename "$f")"
  o="$OBJ/${b%.c}.o"
  gcc -O2 -ffunction-sections -fdata-sections \
    -I"$HIGHS_PREFIX/include" \
    -I"$HIGHS_PREFIX/include/highs" \
    -c "$f" -o "$o"
  ar rcs "$SOSLIB" "$o"
  count=$((count+1))
done

gfortran -O2 -std=legacy -fallow-argument-mismatch \
  -ffunction-sections -fdata-sections \
  -c "$IOSRC" -o "$OBJ/sos_mex_io.o"

echo "Built SOS archive from $count source objects: $SOSLIB"

if ar t "$SOSLIB" | grep -E '^(qpopt|qpcore|lpcore)\.o$'; then
  echo "ERROR: excluded QPOPT/QPCORE/LPCORE object present in archive" >&2
  exit 3
fi

g++ -O2 -std=gnu++17 -shared \
  -DMATLAB_MEX_FILE -DMX_COMPAT_64 \
  -I"$MATLAB_ROOT/extern/include" \
  -o "$OUTDIR/sos_nlp_mex.mexw64" \
  "$WIN_MEXSRC" \
  "$OBJ/sos_mex_io.o" \
  "$MATLAB_ROOT/extern/lib/win64/mingw64/mexFunction.def" \
  "$SOSLIB" \
  -Wl,--whole-archive "$HIGHS_LIB" -Wl,--no-whole-archive \
  -Wl,--gc-sections \
  -Wl,-Bstatic \
  -lgfortran -lquadmath -lstdc++ -lwinpthread -lgcc -lgcc_eh \
  -Wl,-Bdynamic \
  "$MATLAB_ROOT/extern/lib/win64/mingw64/libmex.lib" \
  "$MATLAB_ROOT/extern/lib/win64/mingw64/libmx.lib" \
  "$MATLAB_ROOT/extern/lib/win64/mingw64/libmat.lib" \
  -lkernel32 -luser32 -ladvapi32 -lshell32 -lws2_32

ls -l "$HIGHS_LIB" "$SOSLIB" "$OUTDIR/sos_nlp_mex.mexw64"
