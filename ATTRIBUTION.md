# Attribution

This repository combines the SOS nonlinear programming solver source, a MATLAB
MEX interface, and build support for compiling SOS with HiGHS.

## SOS Solver Source

The SOS solver source code contained primarily in `src/soscode` and
`src/commons` was written by John T. Betts.

The SOS source is distributed in this repository with explicit written
permission from John T. Betts.  See `PERMISSION` and `LICENSE`.

## MATLAB Interface And Build Support

The MATLAB MEX interface, MATLAB build wrappers, platform build scripts, and
SOS/HiGHS integration support were written by Anil V. Rao.
This includes:

```text
src/mex/
src/fortran/
src/highs_bridge/
matlab/
tools/
Makefile
```

## HiGHS

The HiGHS source code is included under `third_party/HiGHS` and is distributed
under its own MIT License.  See `third_party/HiGHS/LICENSE.txt` and
`THIRD_PARTY_NOTICES.txt`.
