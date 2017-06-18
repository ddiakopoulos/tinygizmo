#pragma once

#ifndef gizmin_example_util_hpp
#define gizmin_example_util_hpp

#include <functional>
#include <vector>

#define GLEW_STATIC
#define GL_GLEXT_PROTOTYPES
#include "glew.h"

#define GLFW_INCLUDE_GLU
#include "GLFW\glfw3.h"

#include "gizmo.hpp"
#include "linalg.h"

///////////////////////////////////
//   Windowing & App Lifecycle   //
///////////////////////////////////

struct InputEvent
{
    enum Type { CURSOR, MOUSE, KEY, CHAR, SCROLL };

    GLFWwindow * window;
    linalg::aliases::int2 windowSize;

    Type type;
    int action;
    int mods;

    linalg::aliases::float2 cursor;
    bool drag = false;

    linalg::aliases::uint2 value; // button, key, codepoint, scrollX, scrollY

    bool is_down() const { return action != GLFW_RELEASE; }
    bool is_up() const { return action == GLFW_RELEASE; }

    bool using_shift_key() const { return mods & GLFW_MOD_SHIFT; };
    bool using_control_key() const { return mods & GLFW_MOD_CONTROL; };
    bool using_alt_key() const { return mods & GLFW_MOD_ALT; };
    bool using_super_key() const { return mods & GLFW_MOD_SUPER; };
};

static InputEvent make_input_event(GLFWwindow * window, InputEvent::Type type, const linalg::aliases::float2 cursor, int action)
{
    static bool isDragging = false;

    InputEvent e;
    e.window = window;
    e.type = type;
    e.cursor = cursor;
    e.action = action;
    e.mods = 0;

    glfwGetWindowSize(window, &e.windowSize.x, &e.windowSize.y);

    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) | glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT)) e.mods |= GLFW_MOD_SHIFT;
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) | glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL)) e.mods |= GLFW_MOD_CONTROL;
    if (glfwGetKey(window, GLFW_KEY_LEFT_ALT) | glfwGetKey(window, GLFW_KEY_RIGHT_ALT)) e.mods |= GLFW_MOD_ALT;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SUPER) | glfwGetKey(window, GLFW_KEY_RIGHT_SUPER)) e.mods |= GLFW_MOD_SUPER;

    if (type == InputEvent::MOUSE)
    {
        if (e.is_down()) isDragging = true;
        else if (e.is_up()) isDragging = false;
    }
    e.drag = isDragging;

    return e;
}

class Window
{
    GLFWwindow * window;
public:
    std::function<void(unsigned int codepoint)> on_char;
    std::function<void(int key, int action, int mods)> on_key;
    std::function<void(int button, int action, int mods)> on_mouse_button;
    std::function<void(linalg::aliases::float2 pos)> on_cursor_pos;
    std::function<void(int numFiles, const char ** paths)> on_drop;
    std::function<void(InputEvent e)> on_input;

    Window(int width, int height, const char * title)
    {
        if (glfwInit() == GL_FALSE)
        {
            throw std::runtime_error("glfwInit() failed");
        }

        window = glfwCreateWindow(width, height, title, nullptr, nullptr);

        if (window == nullptr)
        {
            throw std::runtime_error("glfwCreateWindow() failed");
        }

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

#endif // end gizmin_example_util_hpp
