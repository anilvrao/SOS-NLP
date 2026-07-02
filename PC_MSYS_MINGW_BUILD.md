# Windows MSYS/MINGW Source Build

This build starts from source only:

- `src/soscode`: SOS Fortran source with executable `STOP` statements removed;
- `src/highs_bridge`: SOS-to-HiGHS bridge source;
- `src/fortran/sos_mex_io.f`: Fortran MEX I/O helper;
- `src/mex/sos_nlp_mex_clean.cpp`: MEX gateway source;
- `third_party/HiGHS`: HiGHS source;
- `matlab/build_sos_nlp_mex_mingw_windows.m`: MATLAB wrapper;
- `tools/build_sos_nlp_mex_mingw_windows.sh`: MSYS/MINGW build script.

Generated build products are under `build/` and
`third_party/highs-install-win64/`.

## Expected Tools

Install MSYS2 with the MinGW-w64 64-bit toolchain:

```sh
pacman -S --needed mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-cmake mingw-w64-x86_64-make mingw-w64-x86_64-python
```

The script assumes:

```text
C:\msys64\usr\bin\bash.exe
C:\Program Files\MATLAB\R2026a
```

If MATLAB is installed elsewhere, set `MATLAB_ROOT` before running the script
in MSYS:

```sh
export MATLAB_ROOT="/c/Program Files/MATLAB/R2026a"
export MATLAB_RELEASE="R2026a"
```

## Build From MATLAB

```matlab
cd('C:\path\to\sos_nlp_highs_no_stop')
addpath('matlab')
mexPath = build_sos_nlp_mex_mingw_windows()
```

The expected output MEX is:

```text
build\mex\mingw_windows\mexw64\sos_nlp_mex.mexw64
```

## Build Directly From MSYS

```sh
cd /c/path/to/sos_nlp_highs_no_stop
bash tools/build_sos_nlp_mex_mingw_windows.sh
```

## Test In MATLAB

```matlab
addpath('C:\path\to\sos_nlp_highs_no_stop\build\mex\mingw_windows\mexw64','-begin')
which sos_nlp_mex
info = sos_nlp_mex('info')
```

The `info.fortranCompiler` field should contain:

```text
gfortran MinGW Windows + source-built HiGHS
```

## Verification

```sh
rg -n '^[[:space:]]*[sS][tT][oO][pP]([[:space:]]|$|[0-9])' src/soscode
```

The command should produce no matches.
