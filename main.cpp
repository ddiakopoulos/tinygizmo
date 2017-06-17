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

void draw_mesh(GlShader * shader, GlMesh * mesh, const float4x4 & viewProj, const float4x4 & model) 
{
    auto modelViewProjectionMatrix = mul(viewProj, model);
    shader->bind();
    shader->uniform("u_color", float3(1, 1, 1));
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

    int2 windowSize = win->get_window_size();

    wireframeShader.reset(new GlShader(basic_wireframe_vert, basic_wireframe_frag));

    win->on_key = [&](int key, int action, int mods)
    {
    };

    win->on_drop = [&](int numFiles, const char ** paths)
    {
    };

    win->on_mouse_button = [&](int button, int action, int mods)
    {
    };

    win->on_cursor_pos = [&](float2 pos)
    {
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

        //glPushMatrix();
        //glOrtho(0, windowSize.x, windowSize.y, 0, -1, +1);
        //glPopMatrix();

        gl_check_error(__FILE__, __LINE__);

        win->swap_buffers();
    }
    return EXIT_SUCCESS;
}
