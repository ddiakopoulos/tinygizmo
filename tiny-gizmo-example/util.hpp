// This is free and unencumbered software released into the public domain.
// For more information, please refer to <http://unlicense.org>

#pragma once

#ifndef tinygizmo_example_util_hpp
#define tinygizmo_example_util_hpp

#include <functional>
#include <vector>

#define GLEW_STATIC
#define GL_GLEXT_PROTOTYPES
#include "glew.h"

#define GLFW_INCLUDE_GLU
#include "GLFW\glfw3.h"

#include "../tiny-gizmo.hpp"
#include "linalg.h"

///////////////////////////////////
//   Windowing & App Lifecycle   //
///////////////////////////////////

struct camera
{
    float yfov, near_clip, far_clip;
    linalg::aliases::float3 position;
    float pitch, yaw;
    linalg::aliases::float4 get_orientation() const { return qmul(rotation_quat(linalg::aliases::float3(0, 1, 0), yaw), rotation_quat(linalg::aliases::float3(1, 0, 0), pitch)); }
    linalg::aliases::float4x4 get_view_matrix() const { return mul(rotation_matrix(qconj(get_orientation())), translation_matrix(-position)); }
    linalg::aliases::float4x4 get_projection_matrix(const float aspectRatio) const { return linalg::perspective_matrix(yfov, aspectRatio, near_clip, far_clip); }
    linalg::aliases::float4x4 get_viewproj_matrix(const float aspectRatio) const { return mul(get_projection_matrix(aspectRatio), get_view_matrix()); }
};

class Window
{
    GLFWwindow * window;
public:
    std::function<void(unsigned int codepoint)> on_char;
    std::function<void(int key, int action, int mods)> on_key;
    std::function<void(int button, int action, int mods)> on_mouse_button;
    std::function<void(linalg::aliases::float2 pos)> on_cursor_pos;
    std::function<void(int numFiles, const char ** paths)> on_drop;

    Window(int width, int height, const char * title)
    {
        if (glfwInit() == GL_FALSE) throw std::runtime_error("glfwInit() failed");
        window = glfwCreateWindow(width, height, title, nullptr, nullptr);
        if (!window) throw std::runtime_error("glfwCreateWindow() failed");

        glfwMakeContextCurrent(window);

        if (GLenum err = glewInit())
        {
            throw std::runtime_error(std::string("glewInit() failed - ") + (const char *)glewGetErrorString(err));
        }

        glfwSetCharCallback(window, [](GLFWwindow * window, unsigned int codepoint) {
            auto w = (Window *)glfwGetWindowUserPointer(window); if (w->on_char) w->on_char(codepoint);
        });

        glfwSetKeyCallback(window, [](GLFWwindow * window, int key, int, int action, int mods) {
            auto w = (Window *)glfwGetWindowUserPointer(window); if (w->on_key) w->on_key(key, action, mods);
        });

        glfwSetMouseButtonCallback(window, [](GLFWwindow * window, int button, int action, int mods) {
            auto w = (Window *)glfwGetWindowUserPointer(window); if (w->on_mouse_button) w->on_mouse_button(button, action, mods);
        });

        glfwSetCursorPosCallback(window, [](GLFWwindow * window, double xpos, double ypos) {
            auto w = (Window *)glfwGetWindowUserPointer(window); if (w->on_cursor_pos) w->on_cursor_pos(linalg::aliases::float2(linalg::aliases::double2(xpos, ypos)));
        });

        glfwSetDropCallback(window, [](GLFWwindow * window, int numFiles, const char ** paths) {
            auto w = (Window *)glfwGetWindowUserPointer(window); if (w->on_drop) w->on_drop(numFiles, paths);
        });

        glfwSetWindowUserPointer(window, this);
    }

    ~Window()
    {
        glfwMakeContextCurrent(window);
        glfwDestroyWindow(window);
        glfwTerminate();
    }

    Window(const Window &) = delete;
    Window(Window &&) = delete;
    Window & operator = (const Window &) = delete;
    Window & operator = (Window &&) = delete;

    GLFWwindow * get_glfw_window_handle() { return window; };
    bool should_close() const { return !!glfwWindowShouldClose(window); }
    int get_window_attrib(int attrib) const { return glfwGetWindowAttrib(window, attrib); }
    linalg::aliases::int2 get_window_size() const { linalg::aliases::int2 size; glfwGetWindowSize(window, &size.x, &size.y); return size; }
    void set_window_size(linalg::aliases::int2 newSize) { glfwSetWindowSize(window, newSize.x, newSize.y); }
    linalg::aliases::int2 get_framebuffer_size() const { linalg::aliases::int2 size; glfwGetFramebufferSize(window, &size.x, &size.y); return size; }
    linalg::aliases::float2 get_cursor_pos() const { linalg::aliases::double2 pos; glfwGetCursorPos(window, &pos.x, &pos.y); return linalg::aliases::float2(pos); }

    void swap_buffers() { glfwSwapBuffers(window); }
    void close() { glfwSetWindowShouldClose(window, 1); }
};

#endif // end tinygizmo_example_util_hpp
