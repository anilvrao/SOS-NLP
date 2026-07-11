# SOS-NLP Repository Manifest

This repository contains the source code and build support needed to build the
Sparse Optimization Suite nonlinear programming solvers `sprnlp` and `barnlp`
with the open-source HiGHS QP backend.

This repository is NLP-solver-only.  It does not include the SOS optimal
control package or the SOS optimal-control Fortran routines.

## Included Top-Level Files

```text
README.md                    Repository overview
LICENSE                      MIT License
PERMISSION                   Permission to open-source from John T. Betts
ATTRIBUTION.md               Authorship and attribution notes
THIRD_PARTY_NOTICES.txt      Third-party license notices
SOURCE_BUILD.md              Source build notes
PC_MSYS_MINGW_BUILD.md       Windows MSYS2/MinGW build notes
MAC_GFORTRAN_BUILD.md        macOS gfortran build notes
Makefile                     SOS static archive build
```

## Included Source Directories

```text
src/soscode/                 SOS NLP Fortran source
src/commons/                 SOS common-block include files
src/highs_bridge/            SOS-to-HiGHS QP bridge
src/fortran/                 Fortran helper source for the MATLAB MEX gateway
src/mex/                     MATLAB MEX gateway source
third_party/HiGHS/           HiGHS source code
matlab/                      MATLAB build wrappers
tools/                       Platform build scripts
docs/                        SOS documentation
```

## Excluded Material

The following are intentionally not part of this repository:

```text
SOS optimal control Fortran source files
generated object files
generated static libraries
generated MEX files
installed HiGHS build trees
MATLAB temporary files
```

The historical dense QP path based on QPOPT/LPCORE/QPCORE is not used in this
open-source package.  The HiGHS bridge supplies the QP backend used by the SOS
NLP solvers.

## Build Outputs

Build products are generated locally under:

```text
build/
third_party/highs-install-*
```

These generated outputs are ignored by Git.

## No Executable STOP Statements

The SOS Fortran source in this package is the no-STOP variant intended for
library and MATLAB MEX use.  To verify that executable `STOP` statements are
not present in the SOS source, run:

```sh
rg -n '^[[:space:]]*[sS][tT][oO][pP]([[:space:]]|$|[0-9])' src/soscode
```

The command should produce no matches.
