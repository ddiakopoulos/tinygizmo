// This is free and unencumbered software released into the public domain.
// For more information, please refer to <http://unlicense.org>

#include <iostream>
#include <chrono>

#include "util.hpp"
#include "gl-api.hpp"
#include "teapot.h"

using namespace tinygizmo;
using namespace minalg;

const linalg::aliases::float4x4 identity4x4 = { { 1, 0, 0, 0 },{ 0, 1, 0, 0 },{ 0, 0, 1, 0 },{ 0, 0, 0, 1 } };

constexpr const char gizmo_vert[] = R"(#version 330
    layout(location = 0) in vec3 vertex;
    layout(location = 1) in vec3 normal;
    layout(location = 2) in vec4 color;
    out vec4 v_color;
    out vec3 v_world, v_normal;
    uniform mat4 u_mvp;
    void main()
    {
        gl_Position = u_mvp * vec4(vertex.xyz, 1);
        v_color = color;
        v_world = vertex;
        v_normal = normal;
    }
)";

constexpr const char gizmo_frag[] = R"(#version 330
    in vec4 v_color;
    in vec3 v_world, v_normal;
    out vec4 f_color;
    uniform vec3 u_eye;
    void main()
    {
        vec3 light = vec3(1) * max(dot(v_normal, normalize(u_eye - v_world)), 0.50) + 0.25;
        f_color = v_color * vec4(light, 1);
    }
)";

constexpr const char lit_vert[] = R"(#version 330
    uniform mat4 u_modelMatrix;
    uniform mat4 u_viewProj;

    layout(location = 0) in vec3 inPosition;
    layout(location = 1) in vec3 inNormal;

    out vec3 v_position, v_normal;

    void main()
    {
        vec4 worldPos = u_modelMatrix * vec4(inPosition, 1);
        v_position = worldPos.xyz;
        v_normal = normalize((u_modelMatrix * vec4(inNormal,0)).xyz);
        gl_Position = u_viewProj * worldPos;
    }
)";

constexpr const char lit_frag[] = R"(#version 330
    uniform vec3 u_diffuse = vec3(1, 1, 1);
    uniform vec3 u_eye;

    in vec3 v_position;
    in vec3 v_normal;

    out vec4 f_color;
    
    vec3 compute_lighting(vec3 eyeDir, vec3 position, vec3 color)
    {
        vec3 light = vec3(0, 0, 0);
        vec3 lightDir = normalize(position - v_position);
        light += color * u_diffuse * max(dot(v_normal, lightDir), 0);
        vec3 halfDir = normalize(lightDir + eyeDir);
        light += color * u_diffuse * pow(max(dot(v_normal, halfDir), 0), 128);
        return light;
    }

    void main()
    {
        vec3 eyeDir = vec3(0, 1, -2);
        vec3 light = vec3(0, 0, 0);
        light += compute_lighting(eyeDir, vec3(+3, 1, 0), vec3(235.0/255.0, 43.0/255.0, 211.0/255.0));
        light += compute_lighting(eyeDir, vec3(-3, 1, 0), vec3(43.0/255.0, 236.0/255.0, 234.0/255.0));
        f_color = vec4(light + vec3(0.5, 0.5, 0.5), 1.0);
    }
)";

//////////////////////////
//   Main Application   //
//////////////////////////

geometry_mesh make_teapot()
{
    geometry_mesh mesh;
    for (int i = 0; i < 4974; i+=6)
    {
        geometry_vertex v;
        v.position = float3(teapot_vertices[i + 0], teapot_vertices[i + 1], teapot_vertices[i + 2]);
        v.normal = float3(teapot_vertices[i + 3], teapot_vertices[i + 4], teapot_vertices[i + 5]);
        mesh.vertices.push_back(v);
    }
    for (int i = 0; i < 4680; i+=3) mesh.triangles.push_back(uint3(teapot_triangles[i + 0], teapot_triangles[i + 1], teapot_triangles[i + 2]));
    return mesh;
}

void draw_mesh(GlShader & shader, GlMesh & mesh, const linalg::aliases::float3 eye, const linalg::aliases::float4x4 & viewProj, const linalg::aliases::float4x4 & model)
{
    linalg::aliases::float4x4 modelViewProjectionMatrix = mul(viewProj, model);
    shader.bind();
    shader.uniform("u_mvp", modelViewProjectionMatrix);
    shader.uniform("u_eye", eye);
    mesh.draw_elements();
    shader.unbind();
}

void draw_lit_mesh(GlShader & shader, GlMesh & mesh, const linalg::aliases::float3 eye, const linalg::aliases::float4x4 & viewProj, const linalg::aliases::float4x4 & model)
{
    shader.bind();
    shader.uniform("u_viewProj", viewProj);
    shader.uniform("u_modelMatrix", model);
    shader.uniform("u_eye", eye);
    mesh.draw_elements();
    shader.unbind();
}

void upload_mesh(const geometry_mesh & cpu, GlMesh & gpu)
{
    const auto & verts = reinterpret_cast<const std::vector<linalg::aliases::float3> &>(cpu.vertices);
    const auto & tris = reinterpret_cast<const std::vector<linalg::aliases::uint3> &>(cpu.triangles);
    gpu.set_vertices(verts, GL_DYNAMIC_DRAW);
    gpu.set_attribute(0, 3, GL_FLOAT, GL_FALSE, sizeof(geometry_vertex), (GLvoid*) offsetof(geometry_vertex, position));
    gpu.set_attribute(1, 3, GL_FLOAT, GL_FALSE, sizeof(geometry_vertex), (GLvoid*) offsetof(geometry_vertex, normal));
    gpu.set_attribute(2, 4, GL_FLOAT, GL_FALSE, sizeof(geometry_vertex), (GLvoid*) offsetof(geometry_vertex, color));
    gpu.set_elements(tris, GL_DYNAMIC_DRAW);
}

std::unique_ptr<Window> win;

int main(int argc, char * argv[])
{
    bool ml = 0, mr = 0, bf = 0, bl = 0, bb = 0, br = 0;

    camera cam = {};
    cam.yfov = 1.0f;
    cam.near_clip = 0.01f;
    cam.far_clip = 32.0f;
    cam.position = { 0,1.5f,4 };

    gizmo_application_state gizmo_state;
    gizmo_context gizmo_ctx;

    try
    {
        win.reset(new Window(1280, 800, "tiny-gizmo-example-app"));
        glfwSwapInterval(1);
    }
    catch (const std::exception & e)
    {
        std::cout << "Caught GLFW window exception: " << e.what() << std::endl;
    }

    auto windowSize = win->get_window_size();

    GlShader wireframeShader, litShader;
    GlMesh gizmoEditorMesh, teapotMesh;

    wireframeShader = GlShader(gizmo_vert, gizmo_frag);
    litShader = GlShader(lit_vert, lit_frag);

    geometry_mesh teapot = make_teapot();
    upload_mesh(teapot, teapotMesh);

    gizmo_ctx.render = [&](const geometry_mesh & r)
    {
        upload_mesh(r, gizmoEditorMesh);
        draw_mesh(wireframeShader, gizmoEditorMesh, cam.position, cam.get_viewproj_matrix((float) windowSize.x / (float) windowSize.y), identity4x4);
    };

    win->on_key = [&](int key, int action, int mods)
    {
        if (key == GLFW_KEY_LEFT_CONTROL) gizmo_state.hotkey_ctrl = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_L) gizmo_state.hotkey_local = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_T) gizmo_state.hotkey_translate = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_R) gizmo_state.hotkey_rotate = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_S) gizmo_state.hotkey_scale = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_W) bf = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_A) bl = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_S) bb = (action != GLFW_RELEASE);
        if (key == GLFW_KEY_D) br = (action != GLFW_RELEASE);
    };

    win->on_mouse_button = [&](int button, int action, int mods)
    {
        gizmo_state.mouse_left = (action != GLFW_RELEASE);
        if (button == GLFW_MOUSE_BUTTON_LEFT) ml = (action != GLFW_RELEASE);
        if (button == GLFW_MOUSE_BUTTON_RIGHT) mr = (action != GLFW_RELEASE);
    };

    minalg::float2 lastCursor;
    win->on_cursor_pos = [&](linalg::aliases::float2 position)
    {
        gizmo_state.cursor = minalg::float2(position.x, position.y);
        auto deltaCursorMotion = gizmo_state.cursor - lastCursor;
        if (mr)
        {
            cam.yaw -= deltaCursorMotion.x * 0.01f;
            cam.pitch -= deltaCursorMotion.y * 0.01f;
        }
        lastCursor = gizmo_state.cursor;
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

        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glClearColor(0.725f, 0.725f, 0.725f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        auto cameraOrientation = cam.get_orientation();

        // Gizmo input interaction state populated via win->on_input(...) callback above. Update app parameters: 
        gizmo_state.viewport_size = minalg::float2(windowSize.x, windowSize.y);
        gizmo_state.cam.near_clip = cam.near_clip;
        gizmo_state.cam.far_clip = cam.far_clip;
        gizmo_state.cam.yfov = cam.yfov;
        gizmo_state.cam.position = minalg::float3(cam.position.x, cam.position.y, cam.position.z);
        gizmo_state.cam.orientation = minalg::float4(cameraOrientation.x, cameraOrientation.y, cameraOrientation.z, cameraOrientation.w);
  
        glDisable(GL_CULL_FACE);
        auto teapotModelMatrix = reinterpret_cast<const linalg::aliases::float4x4 &>(transform.matrix());
        draw_lit_mesh(litShader, teapotMesh, cam.position, cam.get_viewproj_matrix((float)windowSize.x / (float)windowSize.y), teapotModelMatrix);

        glClear(GL_DEPTH_BUFFER_BIT);

        gizmo_ctx.update(gizmo_state);
        transform_gizmo("xform-example-gizmo", gizmo_ctx, transform);
        gizmo_ctx.draw();

        //std::cout << "Position:    " << transform.position << std::endl;
        //std::cout << "Orientation: " << transform.orientation << std::endl;
        //std::cout << "Scale:       " << transform.scale << std::endl;

        gl_check_error(__FILE__, __LINE__);

        win->swap_buffers();
    }
    return EXIT_SUCCESS;
}
