# macOS Source Build

This build starts from source only:

- `src/soscode`: no-STOP SOS Fortran source;
- `src/highs_bridge`: SOS-to-HiGHS bridge source;
- `src/fortran/sos_mex_io.f`: Fortran MEX I/O helper;
- `src/mex/sos_nlp_mex_clean.cpp`: MEX gateway source;
- `third_party/HiGHS`: HiGHS source.

Generated build products are under `build/` and
`third_party/highs-install-source-mac/`.

## Configure HiGHS

```sh
/opt/homebrew/bin/cmake \
  -S third_party/HiGHS \
  -B build/highs-build-source-mac \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=third_party/highs-install-source-mac \
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

## Build And Install HiGHS

```sh
/opt/homebrew/bin/cmake \
  --build build/highs-build-source-mac \
  --target install \
  --config Release \
  --parallel 8
```

## Build The No-STOP SOS Archive

```sh
make \
  FC=/opt/homebrew/bin/gfortran \
  HIGHS_PREFIX=third_party/highs-install-source-mac \
  OBJ_DIR=build/obj_gfortran_mac \
  LIB=build/libsos_nlp_gfortran_mac.a \
  FFLAGS='-O2 -std=legacy -fallow-argument-mismatch -mmacosx-version-min=12.0 -fPIC' \
  lib
```

## Build The MATLAB MEX

```sh
/usr/local/bin/matlab -batch "cd(pwd); addpath('matlab'); build_sos_nlp_mex_gfortran_mac"
```

The output is:

```text
build/mex/gfortran_mac/<mexext>/sos_nlp_mex.<mexext>
```

## MATLAB Addpath

```matlab
addpath(fullfile(pwd,'build','mex','gfortran_mac',mexext),'-begin')
which sos_nlp_mex
info = sos_nlp_mex('info')
```

## Optional NAG Build

Build the archive with NAG:

```sh
make \
  FC=/usr/local/bin/unagfor \
  HIGHS_PREFIX=third_party/highs-install-source-mac \
  OBJ_DIR=build/obj_nag_mac \
  LIB=build/libsos_nlp_nag_mac.a \
  FFLAGS='-O2 -dusty -PIC' \
  lib
```

Then run:

```sh
/usr/local/bin/matlab -batch "cd(pwd); addpath('matlab'); build_sos_nlp_mex_nag_mac"
```

If NAG licensing requires an environment variable, set it before invoking
MATLAB or in MATLAB before building:

```matlab
setenv('NAG_KUSARI_FILE','/path/to/nag.key')
```

## Verification

```sh
rg -n '^[[:space:]]*[sS][tT][oO][pP]([[:space:]]|$|[0-9])' src/soscode
```

The command should produce no matches.
