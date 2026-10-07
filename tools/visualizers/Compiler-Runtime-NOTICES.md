# Compiler Runtime Notices

MML Visualizers archives may rely on compiler and C/C++ runtime components for the platform used to
build each artifact.

## Current Packaging Notes

- Windows Qt deployment uses `windeployqt --no-compiler-runtime`; MSVC runtime components are not
  intentionally copied by the Qt deployment step.
- Windows WPF artifacts are framework-dependent .NET applications.
- Linux artifacts may rely on system glibc, libstdc++, libgcc, and related runtime libraries unless
  explicitly bundled.
- macOS artifacts may rely on system libc++, system frameworks, and deployment-tool output.

## License

Compiler runtimes and system libraries remain under their own licenses and are not relicensed by MML
Visualizers.

## Required Release Check

For every artifact, inspect the packaged `bin/`, `lib/`, framework, and plugin directories. If MSVC
redistributables, libstdc++, libgcc, libc++, libc, universal CRT, or other compiler/runtime files are
bundled, include the exact redistribution license and notice files required by that toolchain.
