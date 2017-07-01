# immediate-mode gizmo

<p align="center">
  <img src="https://raw.githubusercontent.com/ddiakopoulos/tinygizmo/master/preview.png"/>
</p>

This project is a lightweight, self-contained library for gizmo editing. It includes gizmos for translation, rotation, and scale. Implemented in C++11, the library does not perform rendering directly and instead provides a per-frame buffer of world-space triangles.

An included example illustrates a rendering and key/mouse implementation using GLFW with an OpenGL 3.3 context. Known limitations include hardcoded assumptions about a right-handed, Y-up coordinate system. 

# Motivation

This project was born out of frustrations with other immediate-mode gizmo implementations. Instead of 4x4 matrices and euler angles, the library exposes a `rigid_transform` struct consisting of a position, rotation quaternion, and scale.

# Features
* Both global and local transform modes for the translational and rotational gizmos
* Snap-to-unit (both linear and angular)
  * Set any of the `snap_` values in the `interaction_state` struct. 
* VR ready (the user must call `update(...)` and `draw()` for each eye)
* Hotkeys for transitioning between translation, rotation, and scaling
  * `ctrl-t` to activate the translation gizmo
  * `ctrl-r` to activate the rotation gizmo
  * `ctrl-s` to activate the scale gizmo
  * `ctrl-l` to toggle between global and local transform modes

# Attribution

This project would not have been possible without the gizmo-math reference implementations in the public-domain [workbench](https://github.com/sgorsten/workbench) project. 

# License 

This is free and unencumbered software released into the public domain. For more information, please refer to <http://unlicense.org>