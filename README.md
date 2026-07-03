# SOS NLP HiGHS MEX Build

This repository contains the source needed to build the SOS nonlinear
programming MEX gateway used by GPOPS-II with HiGHS as the QP backend.

The SOS Fortran source in `src/soscode` has been modified so executable
Fortran `STOP` statements are removed and solver exits can return through the
MEX gateway instead of terminating MATLAB.

## Attribution

The SOS solver source code was written by John T. Betts.  The MATLAB MEX
interface, MATLAB build wrappers, platform build scripts, and SOS/HiGHS
integration support were written by Anil V. Rao.

See `ATTRIBUTION.md`, `PERMISSION`, `LICENSE`, and
`THIRD_PARTY_NOTICES.txt` for additional attribution and licensing details.

## Contents

```text
ATTRIBUTION.md                   Authorship and attribution notes
Makefile                         SOS static archive build
src/commons/                     SOS common-block include files
src/soscode/                     no-STOP SOS Fortran source
src/highs_bridge/                SOS-to-HiGHS bridge
src/fortran/sos_mex_io.f         Fortran MEX I/O helper
src/mex/sos_nlp_mex_clean.cpp    MATLAB MEX gateway
matlab/                          MATLAB build wrappers
tools/                           platform build scripts
third_party/HiGHS/               HiGHS source
docs/sosdoc.2025.02.pdf          SOS user's guide
```

Generated files are intentionally not tracked.  Builds create output under
`build/` and install HiGHS under `third_party/highs-install-*`.

## Verify No STOP Statements

Run this from the repository root:

```sh
rg -n '^[[:space:]]*[sS][tT][oO][pP]([[:space:]]|$|[0-9])' src/soscode
```

The command should produce no matches.

## macOS Build

See `SOURCE_BUILD.md`.

Build identifiers:

```text
gfortran_mac wrapper: matlab/build_sos_nlp_mex_gfortran_mac.m
nag_mac wrapper:      matlab/build_sos_nlp_mex_nag_mac.m
gfortran_mac output:  build/mex/gfortran_mac/<mexext>/sos_nlp_mex.<mexext>
nag_mac output:       build/mex/nag_mac/<mexext>/sos_nlp_mex.<mexext>
```

The main macOS MATLAB wrapper is:

```matlab
addpath('matlab')
mexPath = build_sos_nlp_mex_gfortran_mac()
```

There is also a NAG wrapper:

```matlab
addpath('matlab')
mexPath = build_sos_nlp_mex_nag_mac()
```

## Windows MSYS/MINGW Build

See `PC_MSYS_MINGW_BUILD.md`.

Build identifiers:

```text
mingw_windows MATLAB wrapper: matlab/build_sos_nlp_mex_mingw_windows.m
mingw_windows MSYS script:    tools/build_sos_nlp_mex_mingw_windows.sh
mingw_windows output:         build/mex/mingw_windows/mexw64/sos_nlp_mex.mexw64
```

From MATLAB on Windows:

```matlab
cd('C:\path\to\sos_nlp_highs_no_stop')
addpath('matlab')
mexPath = build_sos_nlp_mex_mingw_windows()
```

From an MSYS2 shell:

```sh
cd /c/path/to/sos_nlp_highs_no_stop
bash tools/build_sos_nlp_mex_mingw_windows.sh
```

## MATLAB Usage After Build

Add the generated MEX directory to the MATLAB path:

```matlab
addpath(fullfile(pwd,'build','mex','gfortran_mac',mexext),'-begin')
which sos_nlp_mex
info = sos_nlp_mex('info')
```

Expected `info.fortranCompiler` text identifies the compiler/toolchain and
source-built HiGHS.

## Notes

- `qpopt.f`, `qpcore.f`, and `lpcore.f` are excluded from the SOS archive.
- `qpopt_highs.f` and `src/highs_bridge/sos_highs_qp.c` provide the HiGHS QP
  path.
- On macOS, interrupting a running SOS MEX with CTRL-C may leave MATLAB's MEX
  runtime in a non-reusable state.  Restart MATLAB before another SOS run if
  this occurs.
