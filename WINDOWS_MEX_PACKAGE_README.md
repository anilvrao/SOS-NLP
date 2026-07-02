# Windows MEX Notes

This repository does not track the built Windows MEX file.  Build it from
source using the instructions in `PC_MSYS_MINGW_BUILD.md`.

The expected build output is:

```text
build/mex/mingw_windows/mexw64/sos_nlp_mex.mexw64
```

Build notes:

- HiGHS is built from source using MSYS2/MINGW.
- The SOS Fortran source has no executable `STOP` statements.
- `qpopt.f`, `qpcore.f`, and `lpcore.f` are excluded from the archive.
- `qpopt_highs.f` and `src/highs_bridge/sos_highs_qp.c` are included.

For inclusion in GPOPS-II on Windows, copy the generated MEX to:

```text
nlp/sos/sos_nlp_mex.mexw64
```
