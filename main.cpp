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
    layout(location = 1) in vec3 normal;
    layout(location = 2) in vec3 color;
    out vec3 v_color;
    uniform mat4 u_mvp;
    void main()
    {
        gl_Position = u_mvp * vec4(vertex.xyz, 1);
        v_color = color;
    }
)";

constexpr const char basic_wireframe_frag[] = R"(#version 330
    in vec3 v_color;
    out vec4 f_color;

    void main()
    {
        f_color = vec4(v_color, 1);
    }
)";

//////////////////////////
//   Main Application   //
//////////////////////////

geometry_mesh make_cube(const float3 & min_bounds, const float3 & max_bounds)
{
    const auto a = min_bounds, b = max_bounds;
    geometry_mesh mesh;
    mesh.vertices = {
        { { a.x, a.y, a.z },{ -1,0,0 } },{ { a.x, a.y, b.z },{ -1,0,0 } },
        { { a.x, b.y, b.z },{ -1,0,0 } },{ { a.x, b.y, a.z },{ -1,0,0 } },
        { { b.x, a.y, a.z },{ +1,0,0 } },{ { b.x, b.y, a.z },{ +1,0,0 } },
        { { b.x, b.y, b.z },{ +1,0,0 } },{ { b.x, a.y, b.z },{ +1,0,0 } },
        { { a.x, a.y, a.z },{ 0,-1,0 } },{ { b.x, a.y, a.z },{ 0,-1,0 } },
        { { b.x, a.y, b.z },{ 0,-1,0 } },{ { a.x, a.y, b.z },{ 0,-1,0 } },
        { { a.x, b.y, a.z },{ 0,+1,0 } },{ { a.x, b.y, b.z },{ 0,+1,0 } },
        { { b.x, b.y, b.z },{ 0,+1,0 } },{ { b.x, b.y, a.z },{ 0,+1,0 } },
        { { a.x, a.y, a.z },{ 0,0,-1 } },{ { a.x, b.y, a.z },{ 0,0,-1 } },
        { { b.x, b.y, a.z },{ 0,0,-1 } },{ { b.x, a.y, a.z },{ 0,0,-1 } },
        { { a.x, a.y, b.z },{ 0,0,+1 } },{ { b.x, a.y, b.z },{ 0,0,+1 } },
        { { b.x, b.y, b.z },{ 0,0,+1 } },{ { a.x, b.y, b.z },{ 0,0,+1 } },
    };
    mesh.triangles = { { 0,1,2 },{ 0,2,3 },{ 4,5,6 },{ 4,6,7 },{ 8,9,10 },
    { 8,10,11 },{ 12,13,14 },{ 12,14,15 },{ 16,17,18 },
    { 16,18,19 },{ 20,21,22 },{ 20,22,23 } };
    return mesh;
}

void draw_mesh(GlShader * shader, GlMesh * mesh, const linalg::aliases::float4x4 & viewProj, const linalg::aliases::float4x4 & model)
{
    linalg::aliases::float4x4 modelViewProjectionMatrix = mul(viewProj, model);
    shader->bind();
    shader->uniform("u_mvp", modelViewProjectionMatrix);
    mesh->draw_elements();
    shader->unbind();
}

void upload_mesh(const geometry_mesh & cpu, GlMesh * gpu)
{
    const std::vector<linalg::aliases::float3> & verts = reinterpret_cast<const std::vector<linalg::aliases::float3> &>(cpu.vertices);
    const std::vector<linalg::aliases::uint3> & tris = reinterpret_cast<const std::vector<linalg::aliases::uint3> &>(cpu.triangles);
    gpu->set_vertices(verts, GL_DYNAMIC_DRAW);
    gpu->set_attribute(0, 3, GL_FLOAT, GL_FALSE, sizeof(geometry_vertex), (GLvoid*) offsetof(geometry_vertex, position));
    gpu->set_attribute(1, 3, GL_FLOAT, GL_FALSE, sizeof(geometry_vertex), (GLvoid*) offsetof(geometry_vertex, normal));
    gpu->set_attribute(2, 3, GL_FLOAT, GL_FALSE, sizeof(geometry_vertex), (GLvoid*) offsetof(geometry_vertex, color));
    gpu->set_elements(tris, GL_DYNAMIC_DRAW);
}

const linalg::aliases::float4x4 identity4x4 = { { 1, 0, 0, 0 },{ 0, 1, 0, 0 },{ 0, 0, 1, 0 },{ 0, 0, 0, 1 } };

std::unique_ptr<Window> win;
std::unique_ptr<GlShader> wireframeShader;
std::unique_ptr<GlMesh> gizmoEditorMesh, cubeMesh;

int main(int argc, char * argv[])
{
    bool ml = 0, mr = 0, bf = 0, bl = 0, bb = 0, br = 0;

    camera cam = {};
    cam.yfov = 1.0f;
    cam.near_clip = 0.01f;
    cam.far_clip = 32.0f;
    cam.position = { 0,1.5f,4 };

    try
    {
        win.reset(new Window(1280, 800, "gizmin-example"));
        glfwSwapInterval(1);
    }
    catch (const std::exception & e)
    {
        std::cout << "Caught GLFW window exception: " << e.what() << std::endl;
    }

    linalg::aliases::int2 windowSize = win->get_window_size();
    wireframeShader.reset(new GlShader(basic_wireframe_vert, basic_wireframe_frag));

    gizmoEditorMesh.reset(new GlMesh());
    cubeMesh.reset(new GlMesh());

    gizmo_application_state gizmo_state;
    gizmo_context gizmo_ctx;

    geometry_mesh cube = make_cube(float3(-0.5), float3(0.5));
    upload_mesh(cube, cubeMesh.get());

    gizmo_ctx.render = [&](const geometry_mesh & r)
    {
        upload_mesh(r, gizmoEditorMesh.get());
        draw_mesh(wireframeShader.get(), gizmoEditorMesh.get(), cam.get_viewproj_matrix((float) windowSize.x / (float) windowSize.y), identity4x4);
    };

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

    win->on_cursor_pos = [&](linalg::aliases::float2 position)
    {
        auto e = make_input_event(win->get_glfw_window_handle(), InputEvent::CURSOR, position, 0);
        if (win->on_input) win->on_input(e);
    };

    minalg::float2 lastCursor;

    win->on_input = [&](const InputEvent & event)
    {
        if (event.type == InputEvent::KEY)
        {
            if (event.value[0] == GLFW_KEY_LEFT_CONTROL) gizmo_state.hotkey_ctrl = event.is_down();
            if (event.value[0] == GLFW_KEY_L) gizmo_state.hotkey_local = event.is_down();
            if (event.value[0] == GLFW_KEY_T) gizmo_state.hotkey_translate = event.is_down();
            if (event.value[0] == GLFW_KEY_R) gizmo_state.hotkey_rotate = event.is_down();
            if (event.value[0] == GLFW_KEY_S) gizmo_state.hotkey_scale = event.is_down();
            if (event.value[0] == GLFW_KEY_W) bf = event.is_down();
            if (event.value[0] == GLFW_KEY_A) bl = event.is_down();
            if (event.value[0] == GLFW_KEY_S) bb = event.is_down();
            if (event.value[0] == GLFW_KEY_D) br = event.is_down();
        }
        else if (event.type == InputEvent::MOUSE)
        {
            gizmo_state.mouse_left = event.is_down();
            if (event.value[0] == GLFW_MOUSE_BUTTON_LEFT) ml = event.is_down();
            if (event.value[0] == GLFW_MOUSE_BUTTON_RIGHT) mr = event.is_down();
        }
        else if (event.type == InputEvent::CURSOR)
        {
            gizmo_state.cursor = minalg::float2(event.cursor.x, event.cursor.y);
            auto deltaCursorMotion = gizmo_state.cursor - lastCursor;
            if (mr)
            {
                cam.yaw -= deltaCursorMotion.x * 0.01f;
                cam.pitch -= deltaCursorMotion.y * 0.01f;
            }
            lastCursor = minalg::float2(event.cursor.x, event.cursor.y);
        }
    };

    rigid_transform transform;

    auto t0 = std::chrono::high_resolution_clock::now();
    while (!win->should_close())
    {
        glfwPollEvents();

        auto t1 = std::chrono::high_resolution_clock::now();
        float timestep = std::chrono::duration<float>(t1 - t0).count();
        t0 = t1;

        if (mr)
        {
            const linalg::aliases::float4 orientation = cam.get_orientation();
            linalg::aliases::float3 move;
            if (bf) move -= qzdir(orientation);
            if (bl) move -= qxdir(orientation);
            if (bb) move += qzdir(orientation);
            if (br) move += qxdir(orientation);
            if (length2(move) > 0) cam.position += normalize(move) * (timestep * 10);
        }

        glViewport(0, 0, windowSize.x, windowSize.y);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);  

        linalg::aliases::float4 cameraOrientation = cam.get_orientation();

        // Gizmo input interaction state populated via win->on_input(...) callback above.
        // Now populate other app parameters: 
        gizmo_state.viewport_size = minalg::float2(windowSize.x, windowSize.y);
        gizmo_state.cam.near_clip = cam.near_clip;
        gizmo_state.cam.far_clip = cam.far_clip;
        gizmo_state.cam.yfov = cam.yfov;
        gizmo_state.cam.position = minalg::float3(cam.position.x, cam.position.y, cam.position.z);
        gizmo_state.cam.orientation = minalg::float4(cameraOrientation.x, cameraOrientation.y, cameraOrientation.z, cameraOrientation.w);

        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);

        auto cubeModelMatrix = reinterpret_cast<const linalg::aliases::float4x4 &>(transform.matrix());
        draw_mesh(wireframeShader.get(), cubeMesh.get(), cam.get_viewproj_matrix((float)windowSize.x / (float)windowSize.y), cubeModelMatrix);

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);

        gizmo_ctx.update(gizmo_state);
        transform_gizmo("xform-example-gizmo", gizmo_ctx, transform);
        gizmo_ctx.draw();

        gl_check_error(__FILE__, __LINE__);

        win->swap_buffers();
    }
    return EXIT_SUCCESS;
}
