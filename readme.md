# tinygizmo

<p align="center">
  <img src="https://raw.githubusercontent.com/ddiakopoulos/tinygizmo/master/preview.png"/>
</p>

This project is a lightweight, self-contained library for gizmo editing commonly found in many game engines. It includes mechanisms for manipulating 3d position, rotation, and scale. Implemented in C++11, the library does not perform rendering directly and instead provides a per-frame buffer of world-space triangles. 

An included example is built on top of GLFW (with an OpenGL 3.3 context). Known limitations include hardcoded assumptions about a right-handed, Y-up coordinate system. While the gizmos are provided with vertex normals, the example does not perform any fancy shading. Furthermore, mouse-drag input with certain gizmos at extreme interaction grazing angles is known to produce anomalous output. 

# Motivation

This project was born out of mild frustration with other immediate-mode gizmo libraries. Instead of 4x4 matrices and euler angles, the library exposes a `rigid_transform` consisting of a 3d position, rotation quaternion, and scale. The library is implemented in around 1100 lines of code, which also includes a complete 3d math library in ~400 LoC - [linalg.h](https://github.com/sgorsten/linalg). Alternatives include [ImGuizmo](https://github.com/CedricGuillemet/ImGuizmo) and [Im3D](https://github.com/john-chapman/im3d). Tinygizmo fits in-between these projects by being fully self-contained (no dependency on Dear ImGui), and by being provided in the public domain. 

# Features
* Both global and local transform modes for translational and rotational gizmos
* Snap-to-unit (both linear and angular)
  * Set any of the `snap_` values in the `interaction_state` struct. 
* VR ready (the user must call `update(...)` and `draw()` for each eye)
* Hotkeys for transitioning between translation, rotation, and scaling:
  * `ctrl-t` to activate the translation gizmo
  * `ctrl-r` to activate the rotation gizmo
  * `ctrl-s` to activate the scale gizmo
  * `ctrl-l` to toggle between global and local transform modes

# Attribution

This project would not have been possible without reference implementations in the public-domain [workbench](https://github.com/sgorsten/workbench) project. 

# License 

This is free and unencumbered software released into the public domain. For more information, please refer to <http://unlicense.org>