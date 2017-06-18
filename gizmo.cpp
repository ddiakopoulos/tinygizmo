#include "gizmo.hpp"

inline bool has_clicked(const gizmo_editor::interaction_state & last, const gizmo_editor::interaction_state & active)
{
    if (last.mouse_left == false && active.mouse_left == true) return true;
    return false;
}

inline bool has_released(const gizmo_editor::interaction_state & last, const gizmo_editor::interaction_state & active)
{
    if (last.mouse_left == true && active.mouse_left == false) return true;
    return false;
}

inline bool is_dragging(const gizmo_editor::interaction_state & last, const gizmo_editor::interaction_state & active)
{
    if (last.mouse_left == true && active.mouse_left == true) return true;
    return false;
}

///////////////////////
// Translation Gizmo //
///////////////////////

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

bool intersect_ray_mesh(const ray & ray, const geometry_mesh & mesh, float * hit_t, int * hit_tri)
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
    if (hit_tri) *hit_tri = best_tri;
    return true;
}

///////////////////////
// Translation Gizmo //
///////////////////////

float4x4 get_view_matrix(const camera_parameters & cam) { return mul(rotation_matrix(qconj(cam.orientation)), translation_matrix(-cam.position)); }
float4x4 get_projection_matrix(const rect & viewport, const camera_parameters & cam) { return perspective_matrix(cam.yfov, viewport.aspect_ratio(), cam.near_clip, cam.far_clip); }
float4x4 get_viewproj_matrix(const rect & viewport, const camera_parameters & cam) { return mul(get_projection_matrix(viewport, cam), get_view_matrix(cam)); }

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

ray get_ray_from_pixel(const float2 & pixel, const rect & viewport, const camera_parameters & cam)
{
    const float x = 2 * (pixel.x - viewport.x0) / viewport.width() - 1, y = 1 - 2 * (pixel.y - viewport.y0) / viewport.height();
    const float4x4 inv_view_proj = inverse(get_viewproj_matrix(viewport, cam));
    const float4 p0 = mul(inv_view_proj, float4(x, y, -1, 1)), p1 = mul(inv_view_proj, float4(x, y, +1, 1));
    return{ cam.position, p1.xyz()*p0.w - p0.xyz()*p1.w };
}

geometry_mesh make_box_geometry(const float3 & min_bounds, const float3 & max_bounds)
{
    const auto a = min_bounds, b = max_bounds;
    geometry_mesh mesh;
    mesh.vertices = {
        { { a.x, a.y, a.z },{ -1,0,0 } },
        { { a.x, a.y, b.z },{ -1,0,0 } },
        { { a.x, b.y, b.z },{ -1,0,0 } },
        { { a.x, b.y, a.z },{ -1,0,0 } },
        { { b.x, a.y, a.z },{ +1,0,0 } },
        { { b.x, b.y, a.z },{ +1,0,0 } },
        { { b.x, b.y, b.z },{ +1,0,0 } },
        { { b.x, a.y, b.z },{ +1,0,0 } },
        { { a.x, a.y, a.z },{ 0,-1,0 } },
        { { b.x, a.y, a.z },{ 0,-1,0 } },
        { { b.x, a.y, b.z },{ 0,-1,0 } },
        { { a.x, a.y, b.z },{ 0,-1,0 } },
        { { a.x, b.y, a.z },{ 0,+1,0 } },
        { { a.x, b.y, b.z },{ 0,+1,0 } },
        { { b.x, b.y, b.z },{ 0,+1,0 } },
        { { b.x, b.y, a.z },{ 0,+1,0 } },
        { { a.x, a.y, a.z },{ 0,0,-1 } },
        { { a.x, b.y, a.z },{ 0,0,-1 } },
        { { b.x, b.y, a.z },{ 0,0,-1 } },
        { { b.x, a.y, a.z },{ 0,0,-1 } },
        { { a.x, a.y, b.z },{ 0,0,+1 } },
        { { b.x, a.y, b.z },{ 0,0,+1 } },
        { { b.x, b.y, b.z },{ 0,0,+1 } },
        { { a.x, b.y, b.z },{ 0,0,+1 } },
    };
    mesh.triangles = { { 0,1,2 },{ 0,2,3 },{ 4,5,6 },{ 4,6,7 },{ 8,9,10 },{ 8,10,11 },{ 12,13,14 },{ 12,14,15 },{ 16,17,18 },{ 16,18,19 },{ 20,21,22 },{ 20,22,23 } };
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
        const float angle = static_cast<float>(i%slices) * tau / slices, c = std::cos(angle), s = std::sin(angle);
        const float3x2 mat = { axis, arm1 * c + arm2 * s };
        float3 n = normalize(mat.y); // TODO: Proper normals for each segment
        for (auto & p : points) mesh.vertices.push_back({ mul(mat,p), n });

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

///////////////////////
// Translation Gizmo //
///////////////////////

gizmo_editor::gizmo_editor()
{
    std::initializer_list<float2> arrow_points = { { 0, 0.05f },{ 1, 0.05f },{ 1, 0.10f },{ 1.2f, 0 } };
    std::initializer_list<float2> ring_points = { { +0.05f, 1 },{ -0.05f, 1 },{ -0.05f, 1 },{ -0.05f, 1.2f },{ -0.05f, 1.2f },{ +0.05f, 1.2f },{ +0.05f, 1.2f },{ +0.05f, 1 } };
    geomeshes[0] = make_lathed_geometry({ 1,0,0 }, { 0,1,0 }, { 0,0,1 }, 12, arrow_points);
    geomeshes[1] = make_lathed_geometry({ 0,1,0 }, { 0,0,1 }, { 1,0,0 }, 12, arrow_points);
    geomeshes[2] = make_lathed_geometry({ 0,0,1 }, { 1,0,0 }, { 0,1,0 }, 12, arrow_points);
    geomeshes[3] = make_box_geometry({ -0.01f,0,0 }, { 0.01f,0.4f,0.4f });
    geomeshes[4] = make_box_geometry({ 0,-0.01f,0 }, { 0.4f,0.01f,0.4f });
    geomeshes[5] = make_box_geometry({ 0,0,-0.01f }, { 0.4f,0.4f,0.01f });
    geomeshes[6] = make_lathed_geometry({ 1,0,0 }, { 0,1,0 }, { 0,0,1 }, 24, ring_points);
    geomeshes[7] = make_lathed_geometry({ 0,1,0 }, { 0,0,1 }, { 1,0,0 }, 24, ring_points);
    geomeshes[8] = make_lathed_geometry({ 0,0,1 }, { 1,0,0 }, { 0,1,0 }, 24, ring_points);
}

void gizmo_editor::update(interaction_state & state)
{
    active_state = state;
    drawlist.clear();
}

void gizmo_editor::draw()
{
    if (render)
    {
        geometry_mesh r;
        
        // Combine all gizmo sub-meshes into one super-mesh
        for (auto & m : drawlist)
        {
            uint32_t numVerts = m.mesh.vertices.size();

            r.vertices.insert(r.vertices.end(), m.mesh.vertices.begin(), m.mesh.vertices.end());

            for (auto & f : r.triangles)
            {
                r.triangles.push_back({numVerts + f.x, numVerts + f.y, numVerts + f.z });
            }

            for (auto & v : r.vertices) v.color = m.color; // Take the color and shove it into a per-vertex attribute
        }

        render(r);
    }

    last_state = active_state;
}

void plane_translation_dragger(gizmo_editor & g, const float3 & plane_normal, float3 & point)
{
    // Mouse clicked
    if (has_clicked(g.last_state, g.active_state)) 
    { 
        g.original_position = point; 
    }

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
    }
}

void axis_translation_dragger(gizmo_editor & g, const float3 & axis, float3 & point)
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

void position_gizmo(gizmo_editor & g, int id, float3 & position)
{
    // On click, set the gizmo mode based on which component the user clicked on
    if (has_clicked(g.last_state, g.active_state))
    {
        g.gizmode = gizmo_mode::none;
        auto ray = g.get_ray_from_cursor(g.active_state.cam);
        ray.origin -= position;

        float best_t = std::numeric_limits<float>::infinity(), t;

        if (intersect_ray_mesh(ray, g.geomeshes[0], &t) && t < best_t) { g.gizmode = gizmo_mode::translate_x; best_t = t; }
        if (intersect_ray_mesh(ray, g.geomeshes[1], &t) && t < best_t) { g.gizmode = gizmo_mode::translate_y; best_t = t; }
        if (intersect_ray_mesh(ray, g.geomeshes[2], &t) && t < best_t) { g.gizmode = gizmo_mode::translate_z; best_t = t; }
        if (intersect_ray_mesh(ray, g.geomeshes[3], &t) && t < best_t) { g.gizmode = gizmo_mode::translate_yz; best_t = t; }
        if (intersect_ray_mesh(ray, g.geomeshes[4], &t) && t < best_t) { g.gizmode = gizmo_mode::translate_zx; best_t = t; }
        if (intersect_ray_mesh(ray, g.geomeshes[5], &t) && t < best_t) { g.gizmode = gizmo_mode::translate_xy; best_t = t; }

        if (g.gizmode != gizmo_mode::none)
        {
            g.click_offset = ray.origin + ray.direction*t;
            // g.g.set_pressed(id);
        }
    }

    // TODO
    // If the user has previously clicked on a gizmo component, allow the user to interact with that gizmo
    if (is_dragging(g.last_state, g.active_state))
    {
        position += g.click_offset;
        switch (g.gizmode)
        {
            case gizmo_mode::translate_x: axis_translation_dragger(g, { 1,0,0 }, position); break;
            case gizmo_mode::translate_y: axis_translation_dragger(g, { 0,1,0 }, position); break;
            case gizmo_mode::translate_z: axis_translation_dragger(g, { 0,0,1 }, position); break;
            case gizmo_mode::translate_yz: plane_translation_dragger(g, { 1,0,0 }, position); break;
            case gizmo_mode::translate_zx: plane_translation_dragger(g, { 0,1,0 }, position); break;
            case gizmo_mode::translate_xy: plane_translation_dragger(g, { 0,0,1 }, position); break;
        }
        position -= g.click_offset;
    }

    // On release, deactivate the current gizmo mode
    if (has_released(g.last_state, g.active_state))
    {
        g.gizmode = gizmo_mode::none;
    }

    // Add the gizmo to our 3D draw list
    const float3 colors[] = 
    {
        g.gizmode == gizmo_mode::translate_x ? float3(1,0.5f,0.5f) : float3(1,0,0),
        g.gizmode == gizmo_mode::translate_y ? float3(0.5f,1,0.5f) : float3(0,1,0),
        g.gizmode == gizmo_mode::translate_z ? float3(0.5f,0.5f,1) : float3(0,0,1),
        g.gizmode == gizmo_mode::translate_yz ? float3(0.5f,1,1) : float3(0,1,1),
        g.gizmode == gizmo_mode::translate_zx ? float3(1,0.5f,1) : float3(1,0,1),
        g.gizmode == gizmo_mode::translate_xy ? float3(1,1,0.5f) : float3(1,1,0),
    };

    auto model = translation_matrix(position);

    for (int i = 0; i < 6; ++i)
    {
        gizmo_renderable r;
        r.mesh = g.geomeshes[i];
        r.color = colors[i];
        for (auto & v : r.mesh.vertices) v.position = transform_coord(model, v.position); // transform local coordinates into worldspace
        g.drawlist.push_back(r);
    }
}

void axis_rotation_dragger(gizmo_editor & g, const float3 & axis, const float3 & center, float4 & orientation)
{
    if (g.active_state.mouse_left)
    {
        pose original_pose = { g.original_orientation, g.original_position };
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
            if (d > 0.999f) { orientation = g.original_orientation; return; }

            float angle = std::acos(d);
            if (angle < 0.001f) { orientation = g.original_orientation; return; }

            auto a = normalize(cross(arm1, arm2));
            orientation = qmul(rotation_quat(a, angle), g.original_orientation);
        }
    }
}

void orientation_gizmo(gizmo_editor & g, int id, const float3 & center, float4 & orientation)
{
    auto p = pose(orientation, center);

    // On click, set the gizmo mode based on which component the user clicked on
    if (has_clicked(g.last_state, g.active_state))
    {
        g.gizmode = gizmo_mode::none;
        auto ray = detransform(p, g.get_ray_from_cursor(g.active_state.cam));
        float best_t = std::numeric_limits<float>::infinity(), t;

        if (intersect_ray_mesh(ray, g.geomeshes[6], &t) && t < best_t) { g.gizmode = gizmo_mode::rotate_yz; best_t = t; }
        if (intersect_ray_mesh(ray, g.geomeshes[7], &t) && t < best_t) { g.gizmode = gizmo_mode::rotate_zx; best_t = t; }
        if (intersect_ray_mesh(ray, g.geomeshes[8], &t) && t < best_t) { g.gizmode = gizmo_mode::rotate_xy; best_t = t; }

        if (g.gizmode != gizmo_mode::none)
        {
            g.original_position = center;
            g.original_orientation = orientation;
            g.click_offset = p.transform_point(ray.origin + ray.direction*t);
           // g.g.set_pressed(id);
        }
    }

    // TODO
    /*
    // If the user has previously clicked on a gizmo component, allow the user to interact with that gizmo
    if (g.g.is_pressed(id))
    {
        switch (g.gizmode)
        {
        case gizmo_mode::rotate_yz: axis_rotation_dragger(g, { 1,0,0 }, center, orientation); break;
        case gizmo_mode::rotate_zx: axis_rotation_dragger(g, { 0,1,0 }, center, orientation); break;
        case gizmo_mode::rotate_xy: axis_rotation_dragger(g, { 0,0,1 }, center, orientation); break;
        }
    }
    */

    // On release, deactivate the current gizmo mode
    if (has_released(g.last_state, g.active_state))
    {
        g.gizmode = gizmo_mode::none;
    }

    // Add the gizmo to our 3D draw list
    const float3 colors[] = {
        g.gizmode == gizmo_mode::rotate_yz ? float3(0.5f,1,1) : float3(0,1,1),
        g.gizmode == gizmo_mode::rotate_zx ? float3(1,0.5f,1) : float3(1,0,1),
        g.gizmode == gizmo_mode::rotate_xy ? float3(1,1,0.5f) : float3(1,1,0),
    };

    const auto model = p.matrix();

    for (int i = 6; i < 9; ++i)
    {
        gizmo_renderable r;
        r.mesh = g.geomeshes[i];
        r.color = colors[i - 6];
        for (auto & v : r.mesh.vertices) v.position = transform_coord(model, v.position); // transform local coordinates into worldspace
        g.drawlist.push_back(r);
    }
}
