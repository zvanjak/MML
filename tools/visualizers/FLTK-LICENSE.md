# FLTK Notices

MML Visualizers FLTK release archives may redistribute FLTK runtime libraries or binaries that
statically link FLTK code for the FLTK-backed 2D visualizer family.

## Components

The FLTK-backed visualizers use FLTK through CMake `find_package(FLTK)`.

## License

FLTK is distributed under the GNU Library General Public License version 2 with the FLTK static
linking exception. FLTK components are not relicensed by MML Visualizers. The MML Visualizers
application license remains in `LICENSE.md`.

## Required Release Check

For every FLTK archive, verify whether FLTK is redistributed as shared libraries or statically linked
into the visualizer executables. Include the FLTK upstream license text and any notice files provided
by the package manager or source distribution used for that artifact.

Official FLTK licensing information:

- https://www.fltk.org/COPYING.php
- https://www.fltk.org/doc-1.4/license.html
