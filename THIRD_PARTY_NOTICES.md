# Third-Party Notices

This document identifies third-party software and components used by or associated with wFAST.

The wFAST source code itself is licensed under the MIT License. See the
[LICENSE](LICENSE) file for details.

The MIT License for wFAST does not replace, modify, or override the license
terms of any third-party software described below.

---

## Intel oneAPI Math Kernel Library (oneMKL)

wFAST uses Intel oneAPI Math Kernel Library (oneMKL) included with Intel oneAPI 2025.3.0 for numerical
linear algebra.

In particular, wFAST uses the PARDISO sparse direct solver provided through
Intel oneMKL.

Intel oneMKL is third-party software and is not licensed under the wFAST MIT
License.

The wFAST source repository does not include Intel oneMKL library files.

Users who build wFAST must obtain Intel oneAPI 2025.3.0 (including Intel oneMKL) separately and comply with the
license terms applicable to the Intel software they use.

Intel's software-development-tools licensing overview lists Intel oneMKL
under the Intel Simplified Software License. Users and distributors should
also consult the license files supplied with their installed Intel oneAPI
software because the terms applicable to a particular release may change.

Official Intel licensing overview (August 2024):

https://cdrdv2-public.intel.com/777665/intel-software-development-tools-license%20-overview-august-2024.pdf

The official wFAST build environment uses Intel oneAPI 2025.3.0.
The corresponding Intel license and redistribution notices should be retained
with any binary distribution, as applicable.

---

## Intel OpenMP Runtime

Parallel builds of wFAST use the parallel threading layer of Intel oneMKL and
may depend on the Intel OpenMP runtime.

On 64-bit Windows, this dependency normally includes:

```text
libiomp5md.dll
```

`libiomp5md.dll` is not part of wFAST and is not licensed under the wFAST MIT
License.

The DLL is intentionally not included in the wFAST source repository.

Users of a parallel wFAST build must ensure that the appropriate Intel OpenMP
runtime is installed or otherwise available in the Windows DLL search path.

Any redistribution of Intel runtime components must comply with the Intel
license terms and redistribution conditions applicable to the Intel software
version used to build wFAST.

Official wFAST release packages do not include Intel runtime files.

Users of a parallel wFAST build must install the required Intel OpenMP
runtime separately. In particular, `libiomp5md.dll` must be available through
the user's Intel runtime installation or another valid Windows DLL search
location.

---

## Sequential and Parallel Builds

wFAST can be built using either the sequential or parallel Intel oneMKL
threading layer.

### Parallel build

The parallel build is intended for improved performance on larger analyses.

Typical dependency:

```text
Intel oneMKL
Intel OpenMP runtime (libiomp5md.dll on Windows)
```

### Sequential build

The sequential build uses the sequential Intel oneMKL threading layer.

Depending on the compiler and linker configuration, a sequential build can
avoid the `libiomp5md.dll` runtime dependency.

Even when Intel libraries are statically linked into a wFAST executable,
the incorporated Intel components remain subject to the applicable Intel
license terms.

---

## GiD

wFAST includes a GiD Problem Type that provides preprocessing, analysis
launching, and post-processing integration.

GiD itself is third-party software and is not included in this repository.
Users must obtain GiD separately and use it under the license terms applicable
to their GiD installation.

The files distributed in:

```text
gid/wFAST.gid/
```

are part of the wFAST project unless otherwise identified.

---

## Image and Other Project Assets

The image assets included in the wFAST repository are cleared for
redistribution and use with the project.

These assets are distributed as part of wFAST unless otherwise identified.

---

## No Third-Party License Grant by wFAST

Nothing in the wFAST MIT License grants rights to Intel oneMKL, Intel OpenMP,
GiD, or any other third-party software.

Users are responsible for obtaining any required third-party software and for
complying with its applicable license terms.

---