// This is free and unencumbered software released into the public domain.
// For more information, please refer to <http://unlicense.org>

#include "tiny-gizmo.hpp"

#include <assert.h>
#include <memory>       
#include <vector>
#include <iostream>
#include <functional>
#include <map>

using namespace minalg;
using namespace tinygizmo;

///////////////////////
//   Utility Math    //
///////////////////////

static const float4x4 Identity4x4 = { { 1, 0, 0, 0 },{ 0, 1, 0, 0 },{ 0, 0, 1, 0 },{ 0, 0, 0, 1 } };
static const float3x3 Identity3x3 = { { 1, 0, 0 },{ 0, 1, 0 },{ 0, 0, 1 } };
static const float tau = 6.28318530718f;

void flush_to_zero(float3 & f)
{
    if (std::abs(f.x) < 0.02f) f.x = 0.f;
    if (std::abs(f.y) < 0.02f) f.y = 0.f;
    if (std::abs(f.z) < 0.02f) f.z = 0.f;
}

// 32 bit Fowler–Noll–Vo Hash
uint32_t hash_fnv1a(const std::string & str)
{
    static const uint32_t fnv1aBase32 = 0x811C9DC5u;
    static const uint32_t fnv1aPrime32 = 0x01000193u;

    uint32_t result = fnv1aBase32;

    for (auto & c : str)
    {
        result ^= static_cast<uint32_t>(c);
        result *= fnv1aPrime32;
    }
    return result;
}

float3 snap(const float3 & value, const float snap)
{
    if (snap > 0.0f) return float3(floor(value / snap) * snap);
    return value;
}

float4 make_rotation_quat_axis_angle(const float3 & axis, float angle)
{
    return{ axis * std::sin(angle / 2), std::cos(angle / 2) };
}

float4 make_rotation_quat_between_vectors_snapped(const float3 & from, const float3 & to, const float angle)
{
    auto a = normalize(from);
    auto b = normalize(to);
    auto snappedAcos = std::floor(std::acos(dot(a, b)) / angle) * angle;
    return make_rotation_quat_axis_angle(normalize(cross(a, b)), snappedAcos);
}

template<typename T> T clamp(const T & val, const T & min, const T & max) { return std::min(std::max(val, min), max); }

struct gizmo_mesh_component { geometry_mesh mesh; float3 base_color, highlight_color; };
struct gizmo_renderable { geometry_mesh mesh; float3 color; };

struct ray { float3 origin, direction; };
inline ray transform(const rigid_transform & p, const ray & r) { return{ p.transform_point(r.origin), p.transform_vector(r.direction) }; }
inline ray detransform(const rigid_transform & p, const ray & r) { return{ p.detransform_point(r.origin), p.detransform_vector(r.direction) }; }
inline float3 transform_coord(const float4x4 & transform, const float3 & coord) { auto r = mul(transform, float4(coord, 1)); return (r.xyz() / r.w); }

struct rect
{
    int x0, y0, x1, y1;
    int width() const { return x1 - x0; }
    int height() const { return y1 - y0; }
    int2 dims() const { return{ width(), height() }; }
    float aspect_ratio() const { return (float)width() / height(); }
};

/////////////////////////////////////////
// Ray-Geometry Intersection Functions //
/////////////////////////////////////////

bool intersect_ray_plane(const ray & ray, const float4 & plane, float * hit_t)
{
    float denom = dot(plane.xyz(), ray.direction);
    if (std::abs(denom) == 0) return false;
    if (hit_t) *hit_t = -dot(plane, float4(ray.origin, 1)) / denom;
    return true;
}

bool intersect_ray_triangle(const ray & ray, const float3 & v0, const float3 & v1, const float3 & v2, float * hit_t)
{
    auto e1 = v1 - v0, e2 = v2 - v0, h = cross(ray.direction, e2);
    auto a = dot(e1, h);
    if (std::abs(a) == 0) return false;

    float f = 1 / a;
    auto s = ray.origin - v0;
    auto u = f * dot(s, h);
    if (u < 0 || u > 1) return false;

    auto q = cross(s, e1);
    auto v = f * dot(ray.direction, q);
    if (v < 0 || u + v > 1) return false;

    auto t = f * dot(e2, q);
    if (t < 0) return false;

    if (hit_t) *hit_t = t;
    return true;
}

bool intersect_ray_mesh(const ray & ray, const geometry_mesh & mesh, float * hit_t)
{
    float best_t = std::numeric_limits<float>::infinity(), t;
    int best_tri = -1;
    for (auto & tri : mesh.triangles)
    {
        if (intersect_ray_triangle(ray, mesh.vertices[tri[0]].position, mesh.vertices[tri[1]].position, mesh.vertices[tri[2]].position, &t) && t < best_t)
        {
            best_t = t;
            best_tri = &tri - mesh.triangles.data();
        }
    }
    if (best_tri == -1) return false;
    if (hit_t) *hit_t = best_t;
    return true;
}

/////////////////////////////
// Camera Helper Functions //
/////////////////////////////

float4x4 get_view_matrix(const camera_parameters & cam) { return mul(rotation_matrix(qconj(cam.orientation)), translation_matrix(-cam.position)); }
float4x4 get_projection_matrix(const rect & viewport, const camera_parameters & cam) { return perspective_matrix(cam.yfov, viewport.aspect_ratio(), cam.near_clip, cam.far_clip); }
float4x4 get_viewproj_matrix(const rect & viewport, const camera_parameters & cam) { return mul(get_projection_matrix(viewport, cam), get_view_matrix(cam)); }

ray get_ray_from_pixel(const float2 & pixel, const rect & viewport, const camera_parameters & cam)
{
    const float x = 2 * (pixel.x - viewport.x0) / viewport.width() - 1, y = 1 - 2 * (pixel.y - viewport.y0) / viewport.height();
    const float4x4 inv_view_proj = inverse(get_viewproj_matrix(viewport, cam));
    const float4 p0 = mul(inv_view_proj, float4(x, y, -1, 1)), p1 = mul(inv_view_proj, float4(x, y, +1, 1));
    return{ cam.position, p1.xyz()*p0.w - p0.xyz()*p1.w };
}

///////////////////////////////
// Geometry + Mesh Utilities //
///////////////////////////////

void compute_normals(geometry_mesh & mesh)
{
    for (geometry_vertex & v : mesh.vertices) v.normal = float3();
    for (auto & t : mesh.triangles)
    {
        geometry_vertex & v0 = mesh.vertices[t.x], &v1 = mesh.vertices[t.y], &v2 = mesh.vertices[t.z];
        const float3 n = cross(v1.position - v0.position, v2.position - v0.position);
        v0.normal += n; v1.normal += n; v2.normal += n;
    }
    for (geometry_vertex & v : mesh.vertices)
    {
        v.normal = normalize(v.normal);
    }
}

geometry_mesh make_box_geometry(const float3 & min_bounds, const float3 & max_bounds)
{
    const auto a = min_bounds, b = max_bounds;
    geometry_mesh mesh;
    mesh.vertices = {
        { { a.x, a.y, a.z },{ -1,0,0 } }, { { a.x, a.y, b.z },{ -1,0,0 } },
        { { a.x, b.y, b.z },{ -1,0,0 } }, { { a.x, b.y, a.z },{ -1,0,0 } },
        { { b.x, a.y, a.z },{ +1,0,0 } }, { { b.x, b.y, a.z },{ +1,0,0 } },
        { { b.x, b.y, b.z },{ +1,0,0 } }, { { b.x, a.y, b.z },{ +1,0,0 } },
        { { a.x, a.y, a.z },{ 0,-1,0 } }, { { b.x, a.y, a.z },{ 0,-1,0 } },
        { { b.x, a.y, b.z },{ 0,-1,0 } }, { { a.x, a.y, b.z },{ 0,-1,0 } },
        { { a.x, b.y, a.z },{ 0,+1,0 } }, { { a.x, b.y, b.z },{ 0,+1,0 } },
        { { b.x, b.y, b.z },{ 0,+1,0 } }, { { b.x, b.y, a.z },{ 0,+1,0 } },
        { { a.x, a.y, a.z },{ 0,0,-1 } }, { { a.x, b.y, a.z },{ 0,0,-1 } },
        { { b.x, b.y, a.z },{ 0,0,-1 } }, { { b.x, a.y, a.z },{ 0,0,-1 } },
        { { a.x, a.y, b.z },{ 0,0,+1 } }, { { b.x, a.y, b.z },{ 0,0,+1 } },
        { { b.x, b.y, b.z },{ 0,0,+1 } }, { { a.x, b.y, b.z },{ 0,0,+1 } },
    };
    mesh.triangles = { { 0,1,2 },{ 0,2,3 },{ 4,5,6 },{ 4,6,7 },{ 8,9,10 },
                       { 8,10,11 },{ 12,13,14 },{ 12,14,15 },{ 16,17,18 },
                       { 16,18,19 },{ 20,21,22 },{ 20,22,23 } };
    return mesh;
}

geometry_mesh make_cylinder_geometry(const float3 & axis, const float3 & arm1, const float3 & arm2, int slices)
{
    // Generated curved surface
    geometry_mesh mesh;

    for (uint32_t i = 0; i <= slices; ++i)
    {
        const float tex_s = static_cast<float>(i) / slices, angle = (float)(i%slices) * tau / slices;
        const float3 arm = arm1 * std::cos(angle) + arm2 * std::sin(angle);
        mesh.vertices.push_back({ arm, normalize(arm) });
        mesh.vertices.push_back({ arm + axis, normalize(arm) });
    }
    for (uint32_t i = 0; i<slices; ++i)
    {
        mesh.triangles.push_back({ i * 2, i * 2 + 2, i * 2 + 3 });
        mesh.triangles.push_back({ i * 2, i * 2 + 3, i * 2 + 1 });
    }

    // Generate caps
    uint32_t base = (uint32_t) mesh.vertices.size();
    for (uint32_t i = 0; i<slices; ++i)
    {
        const float angle = static_cast<float>(i%slices) * tau / slices, c = std::cos(angle), s = std::sin(angle);
        const float3 arm = arm1 * c + arm2 * s;
        mesh.vertices.push_back({ arm + axis, normalize(axis) });
        mesh.vertices.push_back({ arm, -normalize(axis) });
    }
    for (uint32_t i = 2; i<slices; ++i)
    {
        mesh.triangles.push_back({ base, base + i * 2 - 2, base + i * 2 });
        mesh.triangles.push_back({ base + 1, base + i * 2 + 1, base + i * 2 - 1 });
    }
    return mesh;
}

geometry_mesh make_lathed_geometry(const float3 & axis, const float3 & arm1, const float3 & arm2, int slices, std::initializer_list<float2> points)
{
    geometry_mesh mesh;
    for (int i = 0; i <= slices; ++i)
    {
        const float angle = (static_cast<float>(i % slices) * tau / slices) + (tau/8.f), c = std::cos(angle), s = std::sin(angle);
        const float3x2 mat = { axis, arm1 * c + arm2 * s };
        float3 n = normalize(mat.y); // TODO: Proper normals for each segment
        for (auto & p : points) mesh.vertices.push_back({ mul(mat, p), n });

        if (i > 0)
        {
            for (size_t j = 1; j < points.size(); ++j)
            {
                uint32_t i0 = (i - 1)*points.size() + (j - 1);
                uint32_t i1 = (i - 0)*points.size() + (j - 1);
                uint32_t i2 = (i - 0)*points.size() + (j - 0);
                uint32_t i3 = (i - 1)*points.size() + (j - 0);
                mesh.triangles.push_back({ i0,i1,i2 });
                mesh.triangles.push_back({ i0,i2,i3 });
            }
        }
    }
    compute_normals(mesh);
    return mesh;
}

//////////////////////////////////
// Gizmo Context Implementation //
//////////////////////////////////

enum class interact
{
    none,
    translate_x, translate_y, translate_z,
    translate_yz, translate_zx, translate_xy,
    translate_xyz,
    rotate_x, rotate_y, rotate_z,
    scale_x, scale_y, scale_z,
    scale_xyz
};

struct gizmo_context::gizmo_context_impl
{
    gizmo_context * ctx;

    gizmo_context_impl(gizmo_context * ctx);

    std::map<interact, gizmo_mesh_component> mesh_components;
    std::vector<gizmo_renderable> drawlist;

    interact interaction_mode;
    transform_mode mode{ transform_mode::translate };

    std::map<uint32_t, bool> active;

    gizmo_application_state active_state, last_state;
    float3 original_position;               // Original position of an object being manipulated with a gizmo
    float4 original_orientation;            // Original orientation of an object being manipulated with a gizmo
    float3 original_scale;                  // Original scale of an object being manipulated with a gizmo
    float3 click_offset;                    // Offset from position of grabbed object to coordinates of clicked point
    bool local_toggle{ false };             // State to describe if the gizmo should use transform-local math
    bool has_clicked{ false };              // State to describe if the user has pressed the left mouse button during the last frame
    bool has_released{ false };             // State to describe if the user has released the left mouse button during the last frame

    ray get_ray_from_cursor(const camera_parameters & cam) const;

    // Public methods
    void update(const gizmo_application_state & state);
    void draw();
};

gizmo_context::gizmo_context_impl::gizmo_context_impl(gizmo_context * ctx) : ctx(ctx)
{
    std::initializer_list<float2> arrow_points  = { { 0.25f, 0 }, { 0.25f, 0.05f },{ 1, 0.05f },{ 1, 0.10f },{ 1.2f, 0 } };
    std::initializer_list<float2> mace_points   = { { 0.25f, 0 }, { 0.25f, 0.05f },{ 1, 0.05f },{ 1, 0.1f },{ 1.25f, 0.1f }, { 1.25f, 0.0f } };
    std::initializer_list<float2> ring_points   = { { +0.025f, 1 },{ -0.025f, 1 },{ -0.025f, 1 },{ -0.025f, 1.1f },{ -0.025f, 1.1f },{ +0.025f, 1.1f },{ +0.025f, 1.1f },{ +0.025f, 1 } };
    mesh_components[interact::translate_x]      = { make_lathed_geometry({ 1,0,0 },{ 0,1,0 },{ 0,0,1 }, 16, arrow_points), { 1,0.5f,0.5f }, { 1,0,0 } };
    mesh_components[interact::translate_y]      = { make_lathed_geometry({ 0,1,0 },{ 0,0,1 },{ 1,0,0 }, 16, arrow_points), { 0.5f,1,0.5f }, { 0,1,0 } };
    mesh_components[interact::translate_z]      = { make_lathed_geometry({ 0,0,1 },{ 1,0,0 },{ 0,1,0 }, 16, arrow_points), { 0.5f,0.5f,1 }, { 0,0,1 } };
    mesh_components[interact::translate_yz]     = { make_box_geometry({ -0.01f,0.25,0.25 },{ 0.01f,0.75f,0.75f }), { 0.5f,1,1 }, { 0,1,1 } };
    mesh_components[interact::translate_zx]     = { make_box_geometry({ 0.25,-0.01f,0.25 },{ 0.75f,0.01f,0.75f }), { 1,0.5f,1 }, { 1,0,1 } };
    mesh_components[interact::translate_xy]     = { make_box_geometry({ 0.25,0.25,-0.01f },{ 0.75f,0.75f,0.01f }), { 1,1,0.5f }, { 1,1,0} };
    mesh_components[interact::translate_xyz]    = { make_box_geometry({ -0.05f,-0.05f,-0.05f },{ 0.05f,0.05f,0.05f }),{ 0.9f, 0.9f, 0.9f },{ 1,1,1 } };
    mesh_components[interact::rotate_x]         = { make_lathed_geometry({ 1,0,0 },{ 0,1,0 },{ 0,0,1 }, 32, ring_points), { 1, 0.5f, 0.5f }, { 1, 0, 0} };
    mesh_components[interact::rotate_y]         = { make_lathed_geometry({ 0,1,0 },{ 0,0,1 },{ 1,0,0 }, 32, ring_points), { 0.5f,1,0.5f }, { 0,1,0 } };
    mesh_components[interact::rotate_z]         = { make_lathed_geometry({ 0,0,1 },{ 1,0,0 },{ 0,1,0 }, 32, ring_points), { 0.5f,0.5f,1}, { 0,0,1 } };
    mesh_components[interact::scale_x]          = { make_lathed_geometry({ 1,0,0 },{ 0,1,0 },{ 0,0,1 }, 4, mace_points),{ 1,0.5f,0.5f },{ 1,0,0 } };
    mesh_components[interact::scale_y]          = { make_lathed_geometry({ 0,1,0 },{ 0,0,1 },{ 1,0,0 }, 4, mace_points),{ 0.5f,1,0.5f },{ 0,1,0 } };
    mesh_components[interact::scale_z]          = { make_lathed_geometry({ 0,0,1 },{ 1,0,0 },{ 0,1,0 }, 4, mace_points),{ 0.5f,0.5f,1 },{ 0,0,1 } };
}

ray gizmo_context::gizmo_context_impl::get_ray_from_cursor(const camera_parameters & cam) const
{
    return get_ray_from_pixel(active_state.cursor, { 0, 0, (int) active_state.viewport_size.x, (int) active_state.viewport_size.y }, cam);
}

void gizmo_context::gizmo_context_impl::update(const gizmo_application_state & state)
{
    active_state = state;
    local_toggle = (!last_state.hotkey_local && active_state.hotkey_local && active_state.hotkey_ctrl) ? !local_toggle : local_toggle;
    has_clicked = (!last_state.mouse_left && active_state.mouse_left) ? true : false;
    has_released = (last_state.mouse_left && !active_state.mouse_left) ? true : false;
    drawlist.clear();
}

void gizmo_context::gizmo_context_impl::draw()
{
    if (ctx->render)
    {
        geometry_mesh r; // Combine all gizmo sub-meshes into one super-mesh
        for (auto & m : drawlist)
        {
            uint32_t numVerts = r.vertices.size();
            auto it = r.vertices.insert(r.vertices.end(), m.mesh.vertices.begin(), m.mesh.vertices.end());
            for (auto & f : m.mesh.triangles) r.triangles.push_back({numVerts + f.x, numVerts + f.y, numVerts + f.z });
            for (; it != r.vertices.end(); ++it) it->color = m.color; // Take the color and shove it into a per-vertex attribute
        }
        ctx->render(r);
    }
    last_state = active_state;
}

///////////////////////////////////
// Private Gizmo Implementations //
///////////////////////////////////

void axis_rotation_dragger(gizmo_context::gizmo_context_impl & g, const float3 & axis, const float3 & center, const float4 & start_orientation, float4 & orientation)
{
    if (g.active_state.mouse_left)
    {
        rigid_transform original_pose = { start_orientation, g.original_position };
        float3 the_axis = original_pose.transform_vector(axis);
        float4 the_plane = { the_axis, -dot(the_axis, g.click_offset) };
        const ray r = g.get_ray_from_cursor(g.active_state.cam);

        float t;
        if (intersect_ray_plane(r, the_plane, &t))
        {
            float3 center_of_rotation = g.original_position + the_axis * dot(the_axis, g.click_offset - g.original_position);
            float3 arm1 = normalize(g.click_offset - center_of_rotation);
            float3 arm2 = normalize(r.origin + r.direction * t - center_of_rotation);

            float d = dot(arm1, arm2);
            if (d > 0.999f) { orientation = start_orientation; return; }

            float angle = std::acos(d);
            if (angle < 0.001f) { orientation = start_orientation; return; }

            if (g.active_state.snap_rotation)
            {
                auto snapped = make_rotation_quat_between_vectors_snapped(arm1, arm2, g.active_state.snap_rotation);
                orientation = qmul(snapped, start_orientation);
            }
            else
            {
                auto a = normalize(cross(arm1, arm2));
                orientation = qmul(rotation_quat(a, angle), start_orientation);
            }
        }
    }
}

void plane_translation_dragger(gizmo_context::gizmo_context_impl & g, const float3 & plane_normal, float3 & point)
{
    // Mouse clicked
    if (g.has_clicked) g.original_position = point; 

    if (g.active_state.mouse_left)
    {
        // Define the plane to contain the original position of the object
        const float3 plane_point = g.original_position;

        // Define a ray emitting from the camera underneath the cursor
        const ray ray = g.get_ray_from_cursor(g.active_state.cam);

        // If an intersection exists between the ray and the plane, place the object at that point
        const float denom = dot(ray.direction, plane_normal);
        if (std::abs(denom) == 0) return;

        const float t = dot(plane_point - ray.origin, plane_normal) / denom;
        if (t < 0) return;

        point = ray.origin + ray.direction * t;

        if (g.active_state.snap_translation) point = snap(point, g.active_state.snap_translation);
    }
}

void axis_translation_dragger(gizmo_context::gizmo_context_impl & g, const float3 & axis, float3 & point)
{
    if (g.active_state.mouse_left)
    {
        // First apply a plane translation dragger with a plane that contains the desired axis and is oriented to face the camera
        const float3 plane_tangent = cross(axis, point - g.active_state.cam.position);
        const float3 plane_normal = cross(axis, plane_tangent);
        plane_translation_dragger(g, plane_normal, point);

        // Constrain object motion to be along the desired axis
        point = g.original_position + axis * dot(point - g.original_position, axis);
    }
}

///////////////////////////////
//   Gizmo Implementations   //
///////////////////////////////

bool intersect(gizmo_context::gizmo_context_impl & g, const ray & r, interact i, float & t, const float best_t)
{
    if (intersect_ray_mesh(r, g.mesh_components[i].mesh, &t) && t < best_t) return true;
    return false;
}

void position_gizmo(const std::string & name, gizmo_context::gizmo_context_impl & g, const float4 & orientation, float3 & position)
{
    auto p = rigid_transform(g.local_toggle ? orientation : float4(0, 0, 0, 1), position);

    auto h = hash_fnv1a(name);

    if (g.has_clicked)
    {
        g.interaction_mode = interact::none;
        auto ray = detransform(p, g.get_ray_from_cursor(g.active_state.cam));

        float best_t = std::numeric_limits<float>::infinity(), t;
        if (intersect(g, ray, interact::translate_x, t, best_t))   { g.interaction_mode = interact::translate_x;   best_t = t; }
        if (intersect(g, ray, interact::translate_y, t, best_t))   { g.interaction_mode = interact::translate_y;   best_t = t; }
        if (intersect(g, ray, interact::translate_z, t, best_t))   { g.interaction_mode = interact::translate_z;   best_t = t; }
        if (intersect(g, ray, interact::translate_yz, t, best_t))  { g.interaction_mode = interact::translate_yz;  best_t = t; }
        if (intersect(g, ray, interact::translate_zx, t, best_t))  { g.interaction_mode = interact::translate_zx;  best_t = t; }
        if (intersect(g, ray, interact::translate_xy, t, best_t))  { g.interaction_mode = interact::translate_xy;  best_t = t; }
        if (intersect(g, ray, interact::translate_xyz, t, best_t)) { g.interaction_mode = interact::translate_xyz; best_t = t; }

        if (g.interaction_mode != interact::none)
        {
            g.click_offset = g.local_toggle ? p.transform_vector(ray.origin + ray.direction*t) : ray.origin + ray.direction*t;
            g.active[h] = true;
        }
        else g.active[h] = false;
    }

    std::vector<float3> axes;
    if (g.local_toggle) axes = { qxdir(p.orientation), qydir(p.orientation), qzdir(p.orientation) };
    else axes = { { 1, 0, 0 },{ 0, 1, 0 },{ 0, 0, 1 } };

    if (g.active[h])
    {
        position += g.click_offset;
        switch (g.interaction_mode)
        {
        case interact::translate_x: axis_translation_dragger(g, axes[0], position); break;
        case interact::translate_y: axis_translation_dragger(g, axes[1], position); break;
        case interact::translate_z: axis_translation_dragger(g, axes[2], position); break;
        case interact::translate_yz: plane_translation_dragger(g, axes[0], position); break;
        case interact::translate_zx: plane_translation_dragger(g, axes[1], position); break;
        case interact::translate_xy: plane_translation_dragger(g, axes[2], position); break;
        case interact::translate_xyz: plane_translation_dragger(g, -minalg::qzdir(g.active_state.cam.orientation), position); break;
        }
        position -= g.click_offset;
    }

    if (g.has_released) g.interaction_mode = interact::none;

    std::vector<interact> draw_interactions
    {
        interact::translate_x, interact::translate_y, interact::translate_z,
        interact::translate_yz, interact::translate_zx, interact::translate_xy,
        interact::translate_xyz
    };

    float4x4 model = p.matrix();

    for (auto c : draw_interactions)
    {
        gizmo_renderable r;
        r.mesh = g.mesh_components[c].mesh;
        r.color = (c == g.interaction_mode) ? g.mesh_components[c].base_color : g.mesh_components[c].highlight_color;
        for (auto & v : r.mesh.vertices) v.position = transform_coord(model, v.position); // transform local coordinates into worldspace
        g.drawlist.push_back(r);
    }
}

void orientation_gizmo(const std::string & name, gizmo_context::gizmo_context_impl & g, const float3 & center, float4 & orientation)
{
    assert(length2(orientation) > float(1e-6));

    auto h = hash_fnv1a(name);

    auto p = rigid_transform(g.local_toggle ? orientation : float4(0, 0, 0, 1), center); // Orientation is local by default

    if (g.has_clicked)
    {
        g.interaction_mode = interact::none;
        auto ray = detransform(p, g.get_ray_from_cursor(g.active_state.cam));
        float best_t = std::numeric_limits<float>::infinity(), t;

        if (intersect(g, ray, interact::rotate_x, t, best_t)) { g.interaction_mode = interact::rotate_x; best_t = t; }
        if (intersect(g, ray, interact::rotate_y, t, best_t)) { g.interaction_mode = interact::rotate_y; best_t = t; }
        if (intersect(g, ray, interact::rotate_z, t, best_t)) { g.interaction_mode = interact::rotate_z; best_t = t; }

        if (g.interaction_mode != interact::none)
        {
            g.original_position = center;
            g.original_orientation = orientation;
            g.click_offset = p.transform_point(ray.origin + ray.direction * t);
            g.active[h] = true;
        }
        else g.active[h] = false;
    }

    float3 activeAxis;
    if (g.active[h])
    {
        const float4 starting_orientation = g.local_toggle ? g.original_orientation : float4(0, 0, 0, 1);
        switch (g.interaction_mode)
        {
        case interact::rotate_x: axis_rotation_dragger(g, { 1, 0, 0 }, center, starting_orientation, p.orientation); activeAxis = { 1, 0, 0 }; break;
        case interact::rotate_y: axis_rotation_dragger(g, { 0, 1, 0 }, center, starting_orientation, p.orientation); activeAxis = { 0, 1, 0 }; break;
        case interact::rotate_z: axis_rotation_dragger(g, { 0, 0, 1 }, center, starting_orientation, p.orientation); activeAxis = { 0, 0, 1 }; break;
        }
    }

    if (g.has_released) g.interaction_mode = interact::none;

    const auto model = p.matrix();

    std::vector<interact> draw_interactions;
    if (!g.local_toggle && g.interaction_mode != interact::none) draw_interactions = { g.interaction_mode };
    else draw_interactions = { interact::rotate_x, interact::rotate_y, interact::rotate_z };

    for (auto c : draw_interactions)
    {
        gizmo_renderable r;
        r.mesh = g.mesh_components[c].mesh;
        r.color = (c == g.interaction_mode) ? g.mesh_components[c].base_color : g.mesh_components[c].highlight_color;
        for (auto & v : r.mesh.vertices) v.position = transform_coord(model, v.position); // transform local coordinates into worldspace
        g.drawlist.push_back(r);
    }

    // For non-local transformations, we only present one rotation ring 
    // and draw an arrow from the center of the gizmo to indicate the degree of rotation
    if (g.local_toggle == false && g.interaction_mode != interact::none)
    {
        // Create orthonormal basis for drawing the arrow
        float3 a = qrot(p.orientation, g.click_offset - g.original_position);
        float3 zDir = normalize(activeAxis), xDir = normalize(cross(a, zDir)), yDir = cross(zDir, xDir);

        // Ad-hoc geometry
        std::initializer_list<float2> arrow_points = { { 0.0f, 0.f },{ 0.0f, 0.05f },{ 0.8f, 0.05f },{ 0.9f, 0.10f },{ 1.0f, 0 } };
        auto geo = make_lathed_geometry(yDir, xDir, zDir, 16, arrow_points);

        gizmo_renderable r;
        r.mesh = geo;
        r.color = float3(1, 1, 1);
        for (auto & v : r.mesh.vertices) v.position = transform_coord(model, v.position);
        g.drawlist.push_back(r);

        orientation = qmul(p.orientation, g.original_orientation);
    }
    else if (g.local_toggle == true && g.interaction_mode != interact::none) orientation = p.orientation;
}

void axis_scale_dragger(gizmo_context::gizmo_context_impl & g, const float3 & axis, const float3 & center, float3 & scale, const bool uniform)
{
    if (g.active_state.mouse_left)
    {
        const float3 plane_tangent = cross(axis, center - g.active_state.cam.position);
        const float3 plane_normal = cross(axis, plane_tangent);

        float3 distance;
        if (g.active_state.mouse_left)
        {
            // Define the plane to contain the original position of the object
            const float3 plane_point = center;

            // Define a ray emitting from the camera underneath the cursor
            const ray ray = g.get_ray_from_cursor(g.active_state.cam);

            // If an intersection exists between the ray and the plane, place the object at that point
            const float denom = dot(ray.direction, plane_normal);
            if (std::abs(denom) == 0) return;

            const float t = dot(plane_point - ray.origin, plane_normal) / denom;
            if (t < 0) return;

            distance = ray.origin + ray.direction * t;
        }

        float3 offset_on_axis = (distance - g.click_offset) * axis;
        flush_to_zero(offset_on_axis);
        float3 new_scale = scale = g.original_scale + offset_on_axis;

        if (uniform) scale = float3(clamp(dot(distance, new_scale), 0.01f, 1000.f));
        else scale = float3(clamp(new_scale.x, 0.01f, 1000.f), clamp(new_scale.y, 0.01f, 1000.f), clamp(new_scale.z, 0.01f, 1000.f));
        if (g.active_state.snap_scale) scale = snap(scale, g.active_state.snap_scale);
    }
}

void scale_gizmo(const std::string & name, gizmo_context::gizmo_context_impl & g, const float4 & orientation, const float3 & center, float3 & scale)
{
    auto h = hash_fnv1a(name);

    auto p = rigid_transform(orientation, center);

    if (g.has_clicked)
    {
        g.interaction_mode = interact::none;
        auto ray = detransform(p, g.get_ray_from_cursor(g.active_state.cam));
        float best_t = std::numeric_limits<float>::infinity(), t;

        if (intersect(g, ray, interact::scale_x, t, best_t)) { g.interaction_mode = interact::scale_x; best_t = t; }
        if (intersect(g, ray, interact::scale_y, t, best_t)) { g.interaction_mode = interact::scale_y; best_t = t; }
        if (intersect(g, ray, interact::scale_z, t, best_t)) { g.interaction_mode = interact::scale_z; best_t = t; }

        if (g.interaction_mode != interact::none)
        {
            g.original_scale = scale;
            g.click_offset = p.transform_point(ray.origin + ray.direction*t);
            g.active[h] = true;
        }
        else g.active[h] = false;
    }

    if (g.has_released) g.interaction_mode = interact::none;

    if (g.active[h])
    {
        switch (g.interaction_mode)
        {
        case interact::scale_x: axis_scale_dragger(g, { 1,0,0 }, center, scale, g.active_state.hotkey_ctrl); break;
        case interact::scale_y: axis_scale_dragger(g, { 0,1,0 }, center, scale, g.active_state.hotkey_ctrl); break;
        case interact::scale_z: axis_scale_dragger(g, { 0,0,1 }, center, scale, g.active_state.hotkey_ctrl); break;
        }
    }

    auto model = p.matrix();

    std::vector<interact> draw_components { interact::scale_x, interact::scale_y, interact::scale_z };

    for (auto c : draw_components)
    {
        gizmo_renderable r;
        r.mesh = g.mesh_components[c].mesh;
        r.color = (c == g.interaction_mode) ? g.mesh_components[c].base_color : g.mesh_components[c].highlight_color;
        for (auto & v : r.mesh.vertices) v.position = transform_coord(model, v.position); // transform local coordinates into worldspace
        g.drawlist.push_back(r);
    }
}

//////////////////////////////////
// Public Gizmo Implementations //
//////////////////////////////////

gizmo_context::gizmo_context() { impl.reset(new gizmo_context_impl(this)); };
gizmo_context::~gizmo_context() { }
void gizmo_context::update(const gizmo_application_state & state) { impl->update(state); }
void gizmo_context::draw() { impl->draw(); }
transform_mode gizmo_context::get_mode() const { return impl->mode; }

void tinygizmo::transform_gizmo(const std::string & name, gizmo_context & g, rigid_transform & t)
{
    if (g.impl->active_state.hotkey_ctrl == true)
    {
        if (g.impl->last_state.hotkey_translate == false && g.impl->active_state.hotkey_translate == true) g.impl->mode = transform_mode::translate;
        else if (g.impl->last_state.hotkey_rotate == false && g.impl->active_state.hotkey_rotate == true) g.impl->mode = transform_mode::rotate;
        else if (g.impl->last_state.hotkey_scale == false && g.impl->active_state.hotkey_scale == true) g.impl->mode = transform_mode::scale;
    }

    if (g.impl->mode == transform_mode::translate) position_gizmo(name, *g.impl, t.orientation, t.position);
    else if (g.impl->mode == transform_mode::rotate) orientation_gizmo(name, *g.impl, t.position, t.orientation);
    else if (g.impl->mode == transform_mode::scale) scale_gizmo(name, *g.impl, t.orientation, t.position, t.scale);
}