#ifndef __EMB_HPP__
#define __EMB_HPP__

#include <embree3/rtcore.h>
#include <vector>
#include <limits>
#include <Vec.hpp>

namespace embree
{

struct EmbreeManager
{
    RTCScene scene;
    RTCDevice device;
    EmbreeManager();
    ~EmbreeManager();
};

EmbreeManager::EmbreeManager()
{
    device = rtcNewDevice(NULL);
    scene = rtcNewScene(device);
}

EmbreeManager::~EmbreeManager()
{
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);
}



void LoadScene2Embree(const std::vector<float>& mesh_vertices, 
                      const std::vector<size_t>& mesh_indices,
                      EmbreeManager* emb)
{
    RTCGeometry geom = rtcNewGeometry(emb->device, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vertices = (float*)rtcSetNewGeometryBuffer(geom,
                                                    RTC_BUFFER_TYPE_VERTEX,
                                                    0,
                                                    RTC_FORMAT_FLOAT3,
                                                    3 * sizeof(float),
                                                    mesh_vertices.size() / 3);
    unsigned* indices = (unsigned*)rtcSetNewGeometryBuffer(geom,
                                                        RTC_BUFFER_TYPE_INDEX,
                                                        0,
                                                        RTC_FORMAT_UINT3,
                                                        3 * sizeof(unsigned),
                                                        mesh_indices.size() / 3);
    
    for(int i = 0; i < mesh_vertices.size(); ++i)
    {
        vertices[i] = mesh_vertices[i];
    }

    for(int i = 0; i < mesh_indices.size(); ++i)
    {
        indices[i] = mesh_indices[i];
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(emb->scene, geom);
    rtcReleaseGeometry(geom);

    rtcCommitScene(emb->scene);
}

struct IsectInfo
{
    size_t primID;

    Vec3<float> Ng;
    Vec3<float> hitPos;

    float tfar;
    float u;
    float v;
};


bool intersect(const Vec3<float>& raypos, const Vec3<float>& raydir, EmbreeManager& embman, IsectInfo& isect)
{
    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    struct RTCRayHit rayhit;
    rayhit.ray.org_x = raypos[0];
    rayhit.ray.org_y = raypos[1];
    rayhit.ray.org_z = raypos[2];
    rayhit.ray.dir_x = raydir[0];
    rayhit.ray.dir_y = raydir[1];
    rayhit.ray.dir_z = raydir[2];
    rayhit.ray.tnear = 0.0f;
    rayhit.ray.tfar = std::numeric_limits<float>::infinity();
    rayhit.ray.mask = 0;
    rayhit.ray.flags = 0;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

    rtcIntersect1(embman.scene, &context, &rayhit);

    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
    {
        isect.primID = rayhit.hit.primID;
        isect.tfar = rayhit.ray.tfar;
        isect.u = rayhit.hit.u;
        isect.v = rayhit.hit.v;
        isect.Ng = (Vec3<float>(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z)).normalize();
        isect.hitPos = raypos + isect.tfar * raydir;
        return true;
    }
    return false;
}


}


#endif