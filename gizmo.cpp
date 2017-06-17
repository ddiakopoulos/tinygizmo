#include "gizmo.hpp"
#include "procedural_mesh.hpp"

using namespace avl;

///////////////////////
// Translation Gizmo //
///////////////////////

/*
struct TranslationGizmo : public IGizmo
{
    Renderable & object;
    ViewportRaycast rc;
    
    float3 axis;
    float3 initialPosition;
    float initialOffset;
    
    float compute_offset(float2 cursor)
    {
        Ray first = {initialPosition, axis}; // From the object origin along the axis of travel
        Ray second = rc.from(cursor); // world space click
        float3 lineDir = second.origin - first.origin;
        const float e1e2 = dot(first.direction, second.direction);
        const float denom = 1 - e1e2 * e1e2;
        const float distance = (dot(lineDir, first.direction) - dot(lineDir, second.direction) * e1e2) / denom;
        return distance;
    }
    
    TranslationGizmo(ViewportRaycast rc, Renderable & object, const float3 & axis, const float2 & cursor) : rc(rc), object(object), axis(qrot(object.pose.orientation, axis)), initialPosition(object.pose.position)
    {
        initialOffset = compute_offset(cursor);
    }
    
    void on_drag(float2 cursor, bool extra) final
    {
        const float offset = compute_offset(cursor);
        object.pose.position = initialPosition + (axis * (offset - initialOffset));
    }
    
    void on_release() final {}
    void on_cancel() final { object.pose.position = initialPosition; }
};

////////////////////
// Rotation Gizmo //
////////////////////

struct RotationGizmo : public IGizmo
{
    Renderable & object;
    ViewportRaycast rc;
    
    float3 axis, edge1;
    float4 initialOrientation;
    
    float3 compute_edge(const float2 & cursor)
    {
        auto ray = rc.from(cursor);
        float3 intersectionPt;
        float t = 0.0;
        intersect_ray_plane(ray, Plane(axis, object.pose.position), &intersectionPt, &t);
        return ray.calculate_position(t) - object.pose.position;
    }
    
    RotationGizmo(ViewportRaycast rc, Renderable & object, const float3 & axis, const float2 & cursor) : rc(rc), object(object), axis(qrot(object.pose.orientation, axis)), initialOrientation(object.pose.orientation)
    {
        edge1 = compute_edge(cursor);
    }
    
    void on_drag(float2 cursor, bool snap) final
    {
        float3 newEdge = compute_edge(cursor);
        if (snap) object.pose.orientation  = qmul(make_rotation_quat_between_vectors_snapped(edge1, newEdge, to_radians(15)), initialOrientation);
        else object.pose.orientation  = qmul(make_rotation_quat_between_vectors(edge1, newEdge), initialOrientation);
    }
    
    void on_release() final {}
    void on_cancel() final { object.pose.orientation = initialOrientation; }
};

///////////////////
// Scaling Gizmo //
///////////////////

struct ScalingGizmo : public IGizmo
{
    Renderable & object;
    ViewportRaycast rc;
    
    float3 axis, initialScale, scaleDirection;
    float initialFactor;
    
    float compute_scale(float2 cursor)
    {
        Ray first = {object.pose.position, axis};
        Ray second = rc.from(cursor);
        float3 lineDir = second.origin - first.origin;
        const float e1e2 = dot(first.direction, second.direction);
        const float denom = 1 - e1e2 * e1e2;
        const float distance = (dot(lineDir, first.direction) - dot(lineDir, second.direction) * e1e2) / denom;
        return distance;
    }
    
    ScalingGizmo(ViewportRaycast rc, Renderable & object, const float3 & axis, const float2 & cursor) : rc(rc), object(object), axis(qrot(object.pose.orientation, axis)), initialScale(object.scale)
    {
        scaleDirection = axis;
        initialFactor = compute_scale(cursor);
    }
    
    void on_drag(float2 cursor, bool scaleUniformly) final
    {
        float scale = compute_scale(cursor) / initialFactor;
        auto scaled = (scaleUniformly) ? initialScale + (scale - 1) : initialScale + scaleDirection * ((scale - 1) * dot(initialScale, scaleDirection));
        object.scale = clamp(scaled, {0.05f, 0.05f, 0.05f}, {100, 100, 100});
    }
    
    void on_release() final {}
    void on_cancel() final { object.scale = initialScale; }
};

//////////////////
// Gizmo Editor //
//////////////////

GizmoEditor::GizmoEditor(GlCamera & camera) : sceneCamera(camera)
{
    make_gizmo_meshes();
}

GizmoEditor::~GizmoEditor()
{

}

std::shared_ptr<IGizmo> GizmoEditor::make_gizmo(ViewportRaycast & rc, Renderable & object, const float3 & axis, const float2 & cursor)
{
    switch (gizmoMode)
    {
        case GizmoMode::Translate: return std::make_shared<TranslationGizmo>(rc, *selectedObject, axis, cursor);
        case GizmoMode::Rotate: return std::make_shared<RotationGizmo>(rc, *selectedObject, axis, cursor);
        case GizmoMode::Scale: return std::make_shared<ScalingGizmo>(rc, *selectedObject, axis, cursor);
    }
}

Renderable & GizmoEditor::get_gizmo_mesh()
{
    switch (gizmoMode)
    {
        case GizmoMode::Translate: return translationMesh;
        case GizmoMode::Rotate: return rotationMesh;
        case GizmoMode::Scale: return scalingMesh;
    }
}

void GizmoEditor::handle_input(const InputEvent & event, std::vector<Renderable> & sceneModels)
{
    if (event.type == InputEvent::MOUSE && event.action == GLFW_PRESS)
    {
        if (event.value[0] == GLFW_MOUSE_BUTTON_LEFT)
        {
            auto worldRay = sceneCamera.get_world_ray(event.cursor, float2(event.windowSize.x, event.windowSize.y));
            
            // If we've already selected an object, check to see if we can interact with the gizmo:
            if (selectedObject)
            {
                float bestT = std::numeric_limits<float>::max();
                float3 hitAxis;
                for (auto axis : {float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1)})
                {
                    // The handle is composed as a piece of geometry on the X axis:
                    auto p = selectedObject->pose * Pose(make_rotation_quat_between_vectors({1,0,0}, axis), {0,0,0});
                    auto localRay = p.inverse() * worldRay;
                    auto hit = get_gizmo_mesh().check_hit(localRay);
                    if (hit.hit && hit.distance <= bestT)
                    {
                        bestT = hit.distance;
                        hitAxis = axis;
                    }
                }
                if (bestT > 0.f && length(hitAxis) > 0)
                {
                    ViewportRaycast raycast(sceneCamera, float2(event.windowSize.x, event.windowSize.y));
                    activeGizmo = make_gizmo(raycast, *selectedObject, hitAxis, event.cursor);
                }
                else
                {
                    activeGizmo.reset();
                }
            }
            
            // We didn't interact with the gizmo, so check if we should select a new object:
            if (!activeGizmo)
            {
                // Otherwise, select a new object
                for (auto & model : sceneModels)
                {
                    // Todo - correct z order
                    auto worldRay = sceneCamera.get_world_ray(event.cursor, float2(event.windowSize.x, event.windowSize.y));
                    auto hit = model.check_hit(worldRay);
                    if (hit.hit)
                    {
                        selectedObject = &model;
                        break;
                    }
                    selectedObject = nullptr;
                }
            }
        }
        
        if (event.value[0] == GLFW_MOUSE_BUTTON_RIGHT)
            activeGizmo.reset();
    }
    
    else if (event.type == InputEvent::CURSOR && event.drag)
    {
        if (selectedObject && activeGizmo)
            activeGizmo->on_drag(event.cursor, event.using_shift_key());
    }
    
    else if (event.type == InputEvent::KEY && event.action == GLFW_RELEASE)
    {
        if(event.value[0] == hotkey_translate) gizmoMode = GizmoMode::Translate;
        if(event.value[0] == hotkey_rotate) gizmoMode = GizmoMode::Rotate;
        if(event.value[0] == hotkey_scale) gizmoMode = GizmoMode::Scale;
        if (event.drag && selectedObject && activeGizmo)
            if (event.value[0] == GLFW_KEY_ESCAPE) activeGizmo->on_cancel();

    }
}

void GizmoEditor::make_gizmo_meshes()
{
    // Create the translation mesh
    Geometry arrowCap = make_cylinder(0.0f, 1.75f, 2.f, 12, 12, false);
    Geometry axisBoxT = make_cube();
    
    Pose p(make_rotation_quat_between_vectors({0, 1, 0}, {1, 0, 0}), {0,0,0});
    for (auto & v : arrowCap.vertices)
        v = p.transform_coord(v);
    
    for (auto & v : arrowCap.vertices)
        v = float3((v.x * 0.25f) + 4.125f, v.y * 0.25f, v.z * 0.25f);
    
    for (auto & v : axisBoxT.vertices)
        v = float3((v.x * 2) + 2.125f, v.y * 0.125f, v.z * 0.125f);
    
    Geometry tGeometry = concatenate_geometry(axisBoxT, arrowCap);
    tGeometry.compute_normals();
    translationMesh = Renderable(tGeometry);
    
    // Create the rotation mesh
    Geometry rGeometry = make_3d_ring(1.50f, 1.70f, 0.20f);
    for (auto & v : rGeometry.vertices)
        v = float3(v.z * 1.50f, v.x * 1.50f, v.y * 1.50f);
    rotationMesh = Renderable(rGeometry);
    
    // Create the scale mesh
    Geometry endBox = make_cube();
    Geometry axisBoxS = make_cube();
    
    for (auto & v : endBox.vertices)
        v = float3((v.x + 16) * 0.25f, v.y * 0.25f, v.z * 0.25f);
    
    for (auto & v : axisBoxS.vertices)
        v = float3((v.x * 2) + 2.125f, v.y * 0.125f, v.z * 0.125f);
    
    Geometry sGeometry = concatenate_geometry(endBox, axisBoxS);
    sGeometry.compute_normals();
    scalingMesh = Renderable(sGeometry);
}
*/