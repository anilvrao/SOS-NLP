# GitHub Package Manifest

This package was assembled from the local no-STOP SOS/HiGHS source build tree
on 2026-07-02.

It includes source and build scripts only.  Generated object files, archives,
installed HiGHS libraries, and MEX binaries were intentionally excluded.

Important included source:

```text
ATTRIBUTION.md
src/soscode
src/commons
src/highs_bridge
src/fortran/sos_mex_io.f
src/mex/sos_nlp_mex_clean.cpp
third_party/HiGHS
```

Included documentation:

```text
docs/sosdoc.2025.02.pdf
```

Important included build entry points:

```text
Makefile
matlab/build_sos_nlp_mex_gfortran_mac.m
matlab/build_sos_nlp_mex_nag_mac.m
matlab/build_sos_nlp_mex_mingw_windows.m
tools/build_sos_nlp_mex_mingw_windows.sh
```

No executable `STOP` statements were found in `src/soscode` by:

```sh
rg -n '^[[:space:]]*[sS][tT][oO][pP]([[:space:]]|$|[0-9])' src/soscode
```
