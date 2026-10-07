# Qt Notices

MML Visualizers Qt release archives may redistribute Qt 6 runtime libraries, Qt plugins, and Qt
deployment-tool output for the Qt-backed visualizer family.

## Components

The Qt-backed visualizers link against these Qt modules in the current build configuration:

- Qt Core
- Qt Gui
- Qt Widgets
- Qt OpenGL Widgets
- Qt platform plugins and support plugins copied by Qt deployment tools when applicable

## License

Qt is available under multiple licensing options from The Qt Company, including commercial Qt
licenses and open-source licenses. MML Visualizers release artifacts that are not built under a
commercial Qt license must comply with the applicable open-source Qt license for the redistributed
Qt components, typically GNU Lesser General Public License version 3 (LGPLv3) for the dynamically
linked Qt runtime libraries and plugins used here.

Qt components are not relicensed by MML Visualizers. Qt license terms apply only to Qt components;
the MML Visualizers application license remains in `LICENSE.md`.

## Required Release Check

For every Qt archive, verify the exact Qt DLLs, shared libraries, frameworks, and plugins shipped by
`windeployqt`, `macdeployqt`, or the Linux install rules. Include Qt's upstream license files and
any module-specific notices required by the Qt distribution used for that artifact.

Official Qt licensing information:

- https://www.qt.io/licensing/
- https://doc.qt.io/qt-6/lgpl.html
- https://doc.qt.io/qt-6/licenses-used-in-qt.html
