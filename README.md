# wFAST

wFAST is a geomechanical finite-element solver with a GiD Problem Type interface for model preparation, analysis execution, and post-processing.

This repository contains:

- the wFAST Fortran solver source code,
- the GiD Problem Type (`wFAST.gid`),
- Windows build files for Intel Fortran / Visual Studio,
- configuration files, and
- example models and supporting documentation.

> **Status:** Initial open-source release preparation.

---

## Features

wFAST currently supports, through the GiD Problem Type:

- solid materials,
- joint materials,
- initial stress conditions,
- pore-pressure conditions,
- temperature conditions,
- displacement boundary conditions,
- traction/loading conditions,
- density changes,
- GiD mesh generation and model-data export,
- execution of the wFAST solver from GiD, and
- GiD post-processing output.

The sparse linear systems are solved using Intel oneMKL PARDISO.

---

## Repository Structure

```text
wFAST/
|
|-- README.md
|-- LICENSE
|-- .gitignore
|-- THIRD_PARTY_NOTICES.md
|
|-- src/
|   `-- wFAST.f90
|
|-- windows/
|   |-- wFAST.sln
|   `-- wFAST.vfproj
|
|-- gid/
|   `-- wFAST.gid/
|       |-- wFAST.xml
|       |-- wFAST.spd
|       |-- wFAST.tcl
|       |-- wFAST.cnd
|       |-- wFAST.bat
|       |
|       |-- exec/
|       |   `-- wFAST.cfg
|       |
|       |-- scripts/
|       |   |-- writing.tcl
|       |   `-- geometry_window.tcl
|       |
|       |-- xml/
|       |   |-- materials.xml
|       |   `-- materials_joint.xml
|       |
|       |-- images/
|       |   |-- wFAST.png
|       |   |-- general.png
|       |   |-- group.png
|       |   |-- initial.png
|       |   |-- loading.png
|       |   |-- material.png
|       |   |-- material1.png
|       |   `-- material2.png
|       |
|       |-- symbols/
|       |
|       `-- resources/
|           `-- msgs/
|
`-- examples/
    `-- ...
```

The exact contents of `src/`, `windows/`, and `examples/` may evolve as the repository is cleaned up for release.

---

## Requirements

### GiD

- GiD 16 or later.

The current Problem Type metadata specifies:

```text
MinimumGiDVersion = 16
```

### Windows build environment

To build the Windows executable, the following are expected:

- Microsoft Visual Studio 2022,
- Intel Fortran Compiler (`ifx` or compatible Intel Fortran compiler),
- Intel oneAPI 2025.3.0 (Intel oneMKL),
- 64-bit Windows build environment.

---

## Development Toolchain

The current wFAST build uses Intel oneAPI 2025.3.0 (Intel oneMKL).


## Building wFAST

wFAST uses Intel oneMKL PARDISO. Two build configurations are possible.

### Parallel oneMKL build

The parallel build is recommended for larger analyses where solver performance is important.

Recommended settings:

```text
Platform:                 x64
Configuration:            Release
Intel oneMKL:             Parallel
PARDISO:                  Parallel oneMKL threading layer
```

A typical Intel Fortran compiler configuration corresponds to the parallel oneMKL option.

The resulting executable depends on the Intel OpenMP runtime on Windows:

```text
libiomp5md.dll
```

`libiomp5md.dll` is not included in the source repository or official
wFAST release packages.

Users of a parallel build must install the required Intel OpenMP runtime
separately and ensure that it is available in the Windows DLL search path.

### Sequential oneMKL build

A sequential build can be used when avoiding the Intel OpenMP runtime dependency is more important than multi-core PARDISO performance.

Recommended settings:

```text
Platform:                 x64
Configuration:            Release
Intel oneMKL:             Sequential
PARDISO:                  Sequential oneMKL threading layer
```

With an appropriate static runtime configuration, this build can avoid the `libiomp5md.dll` dependency.

### Output executable

The solver executable used by the GiD Problem Type must be named:

```text
wFAST.exe
```

Place the executable in:

```text
gid/wFAST.gid/exec/
```

so that the directory contains, at minimum:

```text
exec/
|-- wFAST.exe
`-- wFAST.cfg
```

For source-only use, `wFAST.exe` is intentionally not tracked by Git.

---

## Checking Runtime Dependencies

After building `wFAST.exe`, Windows users can inspect its dynamic dependencies from a Visual Studio Developer Command Prompt:

```bat
dumpbin /DEPENDENTS wFAST.exe
```

For the parallel build, `libiomp5md.dll` is expected.

For the sequential build, `libiomp5md.dll` should normally not appear if the intended sequential/static configuration has been used.

---

## Installing the GiD Problem Type

Copy the complete:

```text
gid/wFAST.gid
```

directory into a GiD Problem Types directory.

After installation, start GiD and select `wFAST` from the Problem Type menu.


---

## Running an Analysis in GiD

A typical workflow is:

1. Start GiD.
2. Select the `wFAST` Problem Type.
3. Create or import the model geometry.
4. Define solid and, if required, joint materials.
5. Assign material properties.
6. Define initial conditions.
7. Define boundary and loading conditions.
8. Generate the mesh.
9. Run `Calculate`.
10. Open the generated results in GiD Postprocess.

The GiD Problem Type creates the wFAST input data, launches `wFAST.exe`, and converts/renames the generated post-processing files so they can be read by GiD.

---

## Solver Configuration

The default solver configuration file is:

```text
gid/wFAST.gid/exec/wFAST.cfg
```

The GiD execution script launches the solver using this configuration file.

Do not rename or relocate this file unless the corresponding path in the GiD execution script is updated.

---

## GiD Execution Files

The GiD Problem Type uses the following main files:

```text
wFAST.xml     Problem Type metadata
wFAST.spd     GiD data-tree and analysis-condition definitions
wFAST.tcl     Main Tcl initialization and callbacks
wFAST.cnd     GiD condition definitions
wFAST.bat     Windows solver-launch script
```

Additional Tcl scripts are stored in:

```text
scripts/
```

Material databases are stored in:

```text
xml/
```

---

## Error and Log Files

The solver and GiD launcher use separate files for diagnostic output and execution errors.

Typical files include:

```text
wfastinp.log
wfast_windows.err
wfast_gid.err
```

Their intended roles are:

- `wfastinp.log`: normal solver log output,
- `wfast_windows.err`: wFAST diagnostic output,
- `wfast_gid.err`: GiD execution error file, created only when the launcher detects a failure.

A successful GiD calculation should not leave a `wfast_gid.err` file containing an error message.

---

## Generated Files

wFAST and GiD generate intermediate and result files during an analysis. These files are not intended to be tracked in the source repository.

Examples include:

```text
gid.dat
wfastinp.dat
wfastinp.res
wfastinp.log
wfast_windows.err
wfast_gid.err
*.post.res
*.post.msh
```

Additional intermediate files may be generated depending on the analysis configuration.

---

## Examples

Example GiD models should be placed under:

```text
examples/
```

A small verification model is recommended so that users can confirm that:

- the GiD Problem Type is installed correctly,
- `wFAST.exe` can be launched,
- the solver completes normally, and
- GiD can read the generated post-processing results.

---

## Binary Releases

The Git repository is intended to contain source files rather than compiled executables.

Precompiled Windows packages may be distributed separately through GitHub Releases.

Possible release packages include:

```text
wFAST-vX.Y.Z-windows-x64-parallel.zip
wFAST-vX.Y.Z-windows-x64-sequential.zip
```

The parallel package requires the Intel OpenMP runtime.

Official wFAST releases do not bundle Intel runtime files. Users of the
parallel build must install the required Intel runtime separately.

---

## Image Assets

All image assets included in the wFAST repository are cleared for
redistribution and use.

## Third-Party Software

wFAST uses Intel oneMKL, including PARDISO, for sparse linear-system solution.

Intel oneMKL and the Intel OpenMP runtime are third-party software and are not part of the wFAST source code.

See:

```text
THIRD_PARTY_NOTICES.md
```

for additional information.

---

## License

wFAST is licensed under the MIT License.

See the [LICENSE](LICENSE) file for details.

Third-party software used by wFAST, including Intel oneMKL and PARDISO,
is subject to its respective license terms. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for details.

The wFAST license does not replace or modify the licenses of third-party software such as Intel oneMKL.

---

## Citation

To be updated.

---

## Authors and Contact

Maintainer: Kenji Furui
Affiliation: Geomechanics and Petroleum Production Engineering Laboratory, Waseda University
Contact: wfast-support@list.waseda.jp


---

## Version

Initial open-source release: v1.0.0
