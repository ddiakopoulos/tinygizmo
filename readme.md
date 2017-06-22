# immediate-mode gizmo

This project is a lightweight, self-contained library for gizmo editing. It includes gizmos for translation, rotation, and scale. Implemented in C++11, the library does not perform rendering directly and instead provides a per-frame buffer of world-space triangles. 

An included example illustrates a rendering implementation using GLFW with an OpenGL 3.3 context. Known limitations include hardcoded assumptions about a Y-up coordinate system. 

# Motivation

This project was born out of frustrations with other immediate-mode gizmo implementations. Instead of transform matrices and euler angles, the library exposes a `rigid_transform` struct consisting of a position, rotation quaternion, and scale. 

# Attribution

This project would not have been possible without the gizmo-math reference implementations in @sgorsten's public-domain [workbench](https://github.com/sgorsten/workbench) project. 

# License 

This is free and unencumbered software released into the public domain. For more information, please refer to <http://unlicense.org>