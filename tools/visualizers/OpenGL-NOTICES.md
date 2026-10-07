# OpenGL And Platform Graphics Notices

MML Visualizers Qt releases use OpenGL through Qt OpenGL Widgets and platform OpenGL/GLU support.

## Components

Current platform handling links or relies on:

- Windows: `opengl32` and `glu32` system libraries.
- macOS: the platform OpenGL framework.
- Linux: OpenGL/GLU libraries discovered by CMake and Qt OpenGL Widgets.

## License

OpenGL implementations, GLU libraries, GPU drivers, and platform graphics frameworks are supplied by
operating-system vendors, GPU vendors, package managers, or system distributions. Those components
remain under their own licenses and are not relicensed by MML Visualizers.

## Required Release Check

Most OpenGL/GLU dependencies are system-provided and are not redistributed in MML Visualizers
archives. If an artifact bundles Mesa, GLU, ANGLE, GPU driver components, or other graphics support
libraries, include the exact upstream license and notice files for those redistributed components.
