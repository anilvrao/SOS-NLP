# Package Manifest

## Included

- `src/soscode`: Fortran source closure for the NLP solvers `sprnlp` and
  `barnlp`.
- `src/commons`: Common-block include files required by the copied Fortran
  source files.
- `src/highs_bridge`: C bridge that maps the legacy QPOPT-compatible call
  signature to HiGHS.
- `third_party/HiGHS`: HiGHS source code and its license files.
- `Makefile`: Convenience build for `build/libsos_nlp.a`.
- `LICENSE`: SOS NLP package MIT license.
- `PERMISSION`: Permission statement for release of the SOS NLP solver source.
- `THIRD_PARTY_NOTICES.txt`: HiGHS license notice.

## Excluded

- MATLAB code, MEX gateways, and GPOPS-II integration files.
- Source code unrelated to the standalone NLP solvers.
- Documentation, examples, drivers, and non-build distribution material.
- Original QPOPT/Optimates core files:
  - `lpcore.f`
  - `qpcore.f`
  - `qpopt.f`
- Compiled objects, static libraries, MEX files, and generated build products.
