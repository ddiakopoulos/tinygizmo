#include <iostream>
#include <functional>
#include <map>
#include <memory>
#include <chrono>
#include <vector>
#include <stdint.h>
#include <complex>
#include <type_traits>
#include <fstream>

#include "util.hpp"
#include "gl-api.hpp"

constexpr const char basic_wireframe_vert[] = R"(#version 330
    layout(location = 0) in vec3 vertex;
    uniform mat4 u_mvp;
    uniform vec3 u_color;
    out vec3 color;
    void main()
    {
        gl_Position = u_mvp * vec4(vertex.xyz, 1);
        color = u_color;
    }
)";

constexpr const char basic_wireframe_frag[] = R"(#version 330
    in vec3 color;
    out vec4 f_color;
    void main()
    {
        f_color = vec4(color.rgb, 1);
    }
)";

//////////////////////////
//   Main Application   //
//////////////////////////

void draw_mesh(GlShader * shader, GlMesh * mesh, const float4x4 & viewProj, const float4x4 & model, const float3 & color) 
{
    auto modelViewProjectionMatrix = mul(viewProj, model);
    shader->bind();
    shader->uniform("u_color", color);
    shader->uniform("u_mvp", modelViewProjectionMatrix);
    mesh->draw_elements();
    shader->unbind();
}

/* 
void upload_mesh(const runtime_mesh & cpu, GlMesh * gpu)
{
    gpu->set_vertices(cpu.vertices, GL_STATIC_DRAW);
    gpu->set_attribute(0, 3, GL_FLOAT, GL_FALSE, sizeof(float3), ((float*)0) + 0);
    gpu->set_non_indexed(GL_LINES);
}
*/

std::unique_ptr<Window> win;
std::unique_ptr<GlShader> wireframeShader;

int main(int argc, char * argv[])
{
    try
    {
        win.reset(new Window(1280, 800, "gizmin-example"));
    }
    catch (const std::exception & e)
    {
        std::cout << "Caught GLFW window exception: " << e.what() << std::endl;
    }

    gizmo_editor::interaction_state gis;
    gizmo_editor gizmoEditor;

    gizmoEditor.render = []()
    {

    };

    int2 windowSize = win->get_window_size();

    wireframeShader.reset(new GlShader(basic_wireframe_vert, basic_wireframe_frag));

    win->on_char = [&](int codepoint)
    {
        auto e = make_input_event(win->get_glfw_window_handle(), InputEvent::CHAR, win->get_cursor_pos(), 0);
        e.value[0] = codepoint;
        if (win->on_input) win->on_input(e);
    };

    win->on_key = [&](int key, int action, int mods)
    {
        auto e = make_input_event(win->get_glfw_window_handle(), InputEvent::KEY, win->get_cursor_pos(), action);
        e.value[0] = key;
        if (win->on_input) win->on_input(e);
    };

    win->on_mouse_button = [&](int button, int action, int mods)
    {
        auto e = make_input_event(win->get_glfw_window_handle(), InputEvent::MOUSE, win->get_cursor_pos(), action);
        e.value[0] = button;
        if (win->on_input) win->on_input(e);
    };

    win->on_cursor_pos = [&](float2 pos)
    {
        auto e = make_input_event(win->get_glfw_window_handle(), InputEvent::CURSOR, win->get_cursor_pos(), 0);
        if (win->on_input) win->on_input(e);
    };

    win->on_input = [&](const InputEvent & event)
    {
        if (event.type == InputEvent::KEY)
        {
            if (event.value[0] == GLFW_KEY_L) gis.hotkey_local = event.is_down();
            if (event.value[0] == GLFW_KEY_T) gis.hotkey_translate = event.is_down();
            if (event.value[0] == GLFW_KEY_R) gis.hotkey_rotate = event.is_down();
            if (event.value[0] == GLFW_KEY_S) gis.hotkey_scale = event.is_down();
        }

        if (event.type == InputEvent::MOUSE)
        {
            gis.mouse_left = event.is_down();
        }

        if (event.type == InputEvent::CURSOR)
        {
            gis.cursor = event.cursor;
        }
    };

    auto t0 = std::chrono::high_resolution_clock::now();
    while (!win->should_close())
    {
        glfwPollEvents();

        auto t1 = std::chrono::high_resolution_clock::now();
        float timestep = std::chrono::duration<float>(t1 - t0).count();
        t0 = t1;

        glViewport(0, 0, windowSize.x, windowSize.y);

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        
        // gizmo input interaction state populated from win->on_input(...) callback above
        // now populate other app parameters: 
        gis.viewport = rect{ 0, 0, windowSize.x, windowSize.y };
        gis.timestep = timestep;
        gis.cam.near_clip = 0.01f;
        gis.cam.near_clip = 64.0f;
        gis.cam.yfov = tau / 2.f;
        gis.cam.position = float3(0, 0, 0);
        gis.cam.orientation = float4(0, 0, 0, 1);

        gizmoEditor.update(gis);

        //glPushMatrix();
        //glOrtho(0, windowSize.x, windowSize.y, 0, -1, +1);
        //glPopMatrix();

        gizmoEditor.draw();

        gl_check_error(__FILE__, __LINE__);

        win->swap_buffers();
    }
    return EXIT_SUCCESS;
}
