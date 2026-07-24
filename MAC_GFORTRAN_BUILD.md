# macOS gfortran Source Build

This file describes the open-source macOS build of the SOS NLP solvers using
Homebrew `gfortran` and source-built HiGHS.

The build uses:

```text
src/soscode/                 no-STOP SOS NLP Fortran source
src/commons/                 SOS common-block include files
src/highs_bridge/            SOS-to-HiGHS QP bridge
src/fortran/sos_mex_io.f     Fortran MEX I/O helper
src/mex/sos_nlp_mex_clean.cpp MATLAB MEX gateway source
third_party/HiGHS/           HiGHS source
matlab/build_sos_nlp_mex_gfortran_mac.m
```

The optimal-control Fortran source files are not part of this repository and
are not required for this NLP solver build.

## Expected Tools

Install the Xcode command line tools and Homebrew.  Then install the build
tools:

```sh
brew install gcc cmake
```

Expected tool locations on Apple Silicon are typically:

```text
/opt/homebrew/bin/gfortran
/opt/homebrew/bin/cmake
```

## Build HiGHS

From the repository root:

```sh
/opt/homebrew/bin/cmake \
  -S third_party/HiGHS \
  -B build/highs-build-gfortran-mac \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=third_party/highs-install-gfortran-mac \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DBUILD_CXX=ON \
  -DBUILD_CXX_EXE=OFF \
  -DBUILD_TESTING=OFF \
  -DBUILD_EXAMPLES=OFF \
  -DZLIB=OFF \
  -DFORTRAN=OFF \
  -DCSHARP=OFF \
  -DPYTHON_BUILD_SETUP=OFF \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=12.0
```

Then build and install HiGHS:

```sh
/opt/homebrew/bin/cmake \
  --build build/highs-build-gfortran-mac \
  --target install \
  --config Release \
  --parallel 8
```

## Build The SOS Static Archive

```sh
make \
  FC=/opt/homebrew/bin/gfortran \
  HIGHS_PREFIX=third_party/highs-install-gfortran-mac \
  OBJ_DIR=build/obj_gfortran_mac \
  LIB=build/libsos_nlp_gfortran_mac.a \
  FFLAGS='-O2 -std=legacy -fallow-argument-mismatch -mmacosx-version-min=12.0 -fPIC' \
  lib
```

## Build The MATLAB MEX File

Start MATLAB from the repository root and run:

```matlab
addpath('matlab')
mexPath = build_sos_nlp_mex_gfortran_mac()
```

The expected output is:

```text
build/mex/gfortran_mac/<mexext>/sos_nlp_mex.<mexext>
```

On Apple Silicon MATLAB, `<mexext>` is usually:

```text
mexmaca64
```

## Test In MATLAB

Add the generated MEX directory to the front of the MATLAB path:

```matlab
addpath(fullfile(pwd,'build','mex','gfortran_mac',mexext),'-begin')
which sos_nlp_mex
info = sos_nlp_mex('info')
```

The `info.fortranCompiler` field should identify a gfortran macOS build with
source-built HiGHS.

## Verify The no-STOP Source

Run:

```sh
rg -n '^[[:space:]]*[sS][tT][oO][pP]([[:space:]]|$|[0-9])' src/soscode
```

The command should produce no matches.

## Notes

- `qpopt.f`, `qpcore.f`, and `lpcore.f` are excluded from the SOS archive.
- `qpopt_highs.f` and `src/highs_bridge/sos_highs_qp.c` provide the HiGHS QP
  path.
- Generated files under `build/` and `third_party/highs-install-*` should not
  be committed.
- The SOS source clears saved solver work arrays at the start of a new run so
  an abandoned reverse-communication sequence does not reuse stale allocated
  state.  A hard MATLAB CTRL-C interrupt can still stop native code at an
  arbitrary point; if MATLAB behaves unpredictably after such an interrupt,
  restart MATLAB before another SOS MEX solve.
