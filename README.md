# SOS NLP Solvers

This repository contains the Fortran source code for the SOS nonlinear
programming solvers:

- `SPRNLP`: active-set SQP NLP solver
- `BARNLP`: barrier NLP solver

The package is intended to be the standalone NLP-solver source package.  It
includes only the NLP solver source code and the dependencies required to
compile it.  It does not include MATLAB MEX files, GPOPS-II files, MATLAB
wrappers, documentation, examples, drivers, or unrelated SOS distribution
material.

## What Is Included

```text
src/soscode/          Fortran source closure for SPRNLP and BARNLP
src/commons/          Common-block include files required by the solvers
src/highs_bridge/     C bridge from the legacy QPOPT-compatible call to HiGHS
third_party/HiGHS/    HiGHS source code
Makefile              Convenience build for a static SOS NLP library
LICENSE               MIT license for the SOS NLP package
PERMISSION            Permission statement for release of the SOS NLP source
THIRD_PARTY_NOTICES.txt
                      HiGHS MIT license notice
```

## What Is Excluded

This package deliberately excludes:

- MATLAB interface files and MEX build files
- GPOPS-II files
- SOS delay-equation source files
- examples and drivers
- SOS documentation and non-build distribution material
- compiled objects, libraries, and generated build artifacts
- original QPOPT/Optimates core files `lpcore.f`, `qpcore.f`, and `qpopt.f`

The QP solve used by the SOS NLP routines is provided through the
HiGHS-backed compatibility layer in `src/soscode/qpopt_highs.f` and
`src/highs_bridge/sos_highs_qp.c`.

## Build

The default build uses `gfortran`, `cc`, `cmake`, and a C++ compiler.  It first
builds and installs HiGHS locally under `third_party/highs-install`, then builds
the SOS NLP static archive:

```sh
make
```

The output archive is:

```text
build/libsos_nlp.a
```

To use a different Fortran compiler:

```sh
make clean
make FC=/path/to/compiler
```

For the NAG Fortran compiler:

```sh
make clean
make FC=nagfor FFLAGS="-O2 -dusty"
```

For a custom wrapper such as `unagfor`:

```sh
make clean
make FC=/usr/local/bin/unagfor FFLAGS="-O2 -dusty"
```

To remove both SOS and HiGHS build products:

```sh
make distclean
```

## Linking

Applications should link their Fortran driver against `build/libsos_nlp.a` and
the HiGHS library built under `third_party/highs-install`.  Because HiGHS is a
C++ library, final executables may also need the C++ runtime library for the
platform/compiler being used.

## License

The SOS NLP package is distributed under the MIT License.  HiGHS is also
distributed under the MIT License; see `THIRD_PARTY_NOTICES.txt` and
`third_party/HiGHS/LICENSE.txt`.
