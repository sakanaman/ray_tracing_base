#ifndef RAYTRACE_HPP
#define RAYTRACE_HPP

#include "bdpt.hpp"

Vec3<float> Trace_Test(const float* firstRay_dir, const float* firstRay_origin, const Shader<float>& shader,
                       const embree::EmbreeManager& emb, RandomManager& rnd_manager, const SceneData<float>& scenedata)
{
    Vec3<float> throughput(1.0f, 1.0f, 1.0f);
    Vec3<float> I(0.0f, 0.0f, 0.0f);
    Ray comingRay(firstRay_origin, firstRay_dir);
    float Pr = 1.0;
    for(int depth = 0;;++depth)
    {
        embree::IsectInfo<float> isect;
        bool is_hit = embree::intersect<float>(comingRay.origin, comingRay.dir, 
                                                embree::RAYMIN, embree::RAYMAX,
                                                emb, isect);
        
        if(is_hit)
        {
            Vec3<float> wo = -1.0f * comingRay.dir;
            Vec3<float> n = isect.Ng;
            Vec3<float> hitPos = isect.hitPos;
            Vec3<float> orienting_normal = n.dot(comingRay.dir) < 0 ? n : -1.0f * n;

            Vec3<float> e0, e2;
            ONB(orienting_normal, e0, e2);
            Vec3<float> local_wo{e0.dot(wo),
                                 orienting_normal.dot(wo),
                                 e2.dot(wo)};
            

            #ifdef NEE
            float phi, theta;
            float Le[3] = {0.0f, 0.0f, 0.0f};
            float r1 = rnd_manager.GetRND(), r2 = rnd_manager.GetRND();
            float nee_pdf = scenedata.ibl.sample(r1, r2, &phi, &theta, Le);

            Vec3<float> shadowdir = {std::sin(theta) * std::cos(phi),
                                     std::cos(theta),
                                     std::sin(theta) * std::sin(phi)};
            Vec3<float> local_shadowdir = {e0.dot(shadowdir),
                                           orienting_normal.dot(shadowdir),
                                           e2.dot(shadowdir)};
            
            if(shadowdir.dot(orienting_normal) > 0.0f)
            {
                Vec3<float> shadoworigin  = hitPos /*+ 0.001f * orienting_normal*/;
                embree::IsectInfo<float> shadow_hit;
                if(embree::intersect<float>(shadoworigin, shadowdir, 
                                            embree::RAYMIN, embree::RAYMAX,
                                            emb, shadow_hit))
                {
                    //nothing
                }
                else
                {
                    float costerm = local_shadowdir[1];
                    I = I + throughput * eval_bsdf(local_shadowdir, local_wo, isect, shader, scenedata) 
                                       * (1.0f/nee_pdf) * costerm * Vec3<float>(Le);
                }
            }
            #endif

            Vec3<float> local_wi = sample(local_wo, isect, rnd_manager, shader, scenedata);
            Vec3<float> wi = local_wi[0]*e0 + 
                             local_wi[1]*orienting_normal + 
                             local_wi[2]*e2;
            comingRay = Ray(hitPos /*+ 0.001f * orienting_normal*/, wi);

            throughput = throughput * weight(local_wi, local_wo, isect, shader, scenedata);

            Pr *= 0.96;
            if(rnd_manager.GetRND() < Pr)
            {
                throughput  = throughput * (1.0f/Pr);
            }
            else
            {
                break;
            }
        }
        else
        {
            #ifdef NEE
            if(depth == 0)
            {
                float phi = std::atan2(comingRay.dir[2], comingRay.dir[0]);
                float theta = std::acos(comingRay.dir[1]);
                if (phi < 0)
                    phi += 2 * M_PI;
                if (phi > 2 * M_PI)
                    phi -= 2 * M_PI;

                float Le[3];
                scenedata.ibl.GetLe(Le, theta, phi);
                I = I + throughput * Vec3<float>(Le);
            }
            #else
                float phi = std::atan2(comingRay.dir[2], comingRay.dir[0]);
                float theta = std::acos(comingRay.dir[1]);
                if (phi < 0)
                    phi += 2 * M_PI;
                if (phi > 2 * M_PI)
                    phi -= 2 * M_PI;

                float Le[3];
                scenedata.ibl.GetLe(Le, theta, phi);
                I = I + throughput * Vec3<float>(Le);
            #endif
            break;
        }   
    }
    return I;
}


Vec3<float> Trace_debug(const float* firstRay_dir, const float* firstRay_origin,
                        const embree::EmbreeManager& emb, RandomManager& rnd_manager, const SceneData<float>& scenedata)
{
    Vec3<float> light_dir = Vec3<float>(-1.0f, 2.0f, 0.0f).normalize();

    Vec3<float> throughput(1.0f, 1.0f, 1.0f);
    Vec3<float> I(0.0f, 0.0f, 0.0f);
    Ray comingRay(firstRay_origin, firstRay_dir);

    embree::IsectInfo<float> hit;
    bool is_hit = embree::intersect(comingRay.origin, comingRay.dir, 
                                    embree::RAYMIN, embree::RAYMAX,
                                    emb, hit);

    if(is_hit)
    {
        Vec3<float> n = hit.Ng;
        Vec3<float> hitPos = hit.hitPos;
        Vec3<float> orienting_normal = n.dot(comingRay.dir) < 0 ? n : -1.0f * n;
        Vec3<float> norm_color{(n[0] + 1)*0.5f,
                               (n[1] + 1)*0.5f,
                               (n[2] + 1)*0.5f};

        // Vec3<float> shadowdir = light_dir;
        // Vec3<float> shadoworigin  =  hitPos /*+ shadowdir * float(1e-3)*/;
        // embree::IsectInfo<float> shadow_hit;
        // if(embree::intersect<float>(shadoworigin, shadowdir, 
        //                             embree::RAYMIN, embree::RAYMAX,
        //                             emb, shadow_hit))
        // {
        //     if(shadow_hit.primID == hit.primID) //oh..self intersection...
        //     {
        //         return {1.0f, 0.0f, 0.0f};
        //     }
        //     return {1.0, 
        //             1.0,
        //             1.0};
        // }
        // return {1.0, 1.0, 1.0};
        // // return norm_color * norm_color;

        int geomID = hit.primID;
        int matID = scenedata.mat_infos.mat_indices[geomID];
        float u = hit.u;
        float v = hit.v;

        int uv_index0 = scenedata.vertex_infos.uv_indices[3 * geomID + 0]; 
        int uv_index1 = scenedata.vertex_infos.uv_indices[3 * geomID + 1]; 
        int uv_index2 = scenedata.vertex_infos.uv_indices[3 * geomID + 2]; 

        float U0 = scenedata.vertex_infos.uvs[2 * uv_index0 + 0];
        float V0 = scenedata.vertex_infos.uvs[2 * uv_index0 + 1];
        float U1 = scenedata.vertex_infos.uvs[2 * uv_index1 + 0];
        float V1 = scenedata.vertex_infos.uvs[2 * uv_index1 + 1];
        float U2 = scenedata.vertex_infos.uvs[2 * uv_index2 + 0];
        float V2 = scenedata.vertex_infos.uvs[2 * uv_index2 + 1];

        float U = U0 + u * (U1 - U0) + v * (U2 - U0);
        float V = V0 + u * (V1 - V0) + v * (V2 - V0);

        U = std::fmod(U, 1.0); 
        if(U < 0.0) U += 1.0;
        V = std::fmod(V, 1.0); 
        if(V < 0.0) V += 1.0;

        Vec3<float> albedo;
        if(scenedata.mat_infos.albedo_textures[matID].exist())
        {
            scenedata.mat_infos.albedo_textures[matID].getColor(U, V, &(albedo[0]), &(albedo[1]), &(albedo[2]));
        }
        else
        {
            albedo = {scenedata.mat_infos.diffuses[3 * matID + 0],
                      scenedata.mat_infos.diffuses[3 * matID + 1],
                      scenedata.mat_infos.diffuses[3 * matID + 2]};
        }
        return albedo;
        // return {U, V, 0};
    }
    else
    {
        return {0.0f, 0.0f, 0.0f};
    }
}

#endif