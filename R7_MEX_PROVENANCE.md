# R7 SOS MEX Provenance

This note records the Apple Silicon macOS SOS MEX used in the GPOPS-II R7
release.

## R7 MEX Fingerprint

```text
File:   nlp/sos/sos_nlp_mex.mexmaca64
SHA256: 2f0a365b8ea5ef4586a977e1fe88f30241e3233fb7d1de4108fbbcbd21536b4c
Size:   4872840 bytes
```

The matching local build artifact was verified at:

```text
/Users/anilvrao/Documents/MATLAB/SOS-NLP-old-20260711/build/mex/gfortran_mac/mexmaca64/sos_nlp_mex.mexmaca64
```

That artifact is byte-for-byte identical to the R7 MEX file located at:

```text
/Users/anilvrao/Documents/MATLAB/GPOPS-II-Assembla/R7-June-2026/nlp/sos/sos_nlp_mex.mexmaca64
```

## Build Lineage

The R7 MEX reports the following embedded build identifiers:

```text
gfortran macOS + source-built HiGHS
build/libsos_nlp_gfortran_mac.a + source-built HiGHS
```

The source, MATLAB wrapper, and tool files in this repository match the
R7-producing local tree. Generated files such as `build/`, static archives,
installed HiGHS libraries, and platform-specific MEX binaries are intentionally
not tracked in Git.

To rebuild on Apple Silicon macOS, follow `MAC_GFORTRAN_BUILD.md` and run:

```matlab
addpath('matlab')
mexPath = build_sos_nlp_mex_gfortran_mac()
```

A newly rebuilt MEX may not be byte-for-byte identical to the R7 binary because
static archives and linker output can contain build-time ordering or toolchain
metadata. The source and build recipe are the relevant reproducibility record.
