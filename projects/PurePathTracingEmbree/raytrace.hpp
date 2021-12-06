#ifndef RAYTRACE_HPP
#define RAYTRACE_HPP

#include "material.hpp"
#include "IBL.hpp"
#include <concurrent.hpp>
#include <wrap_embree.hpp>

class Ray
{
public:
    Ray(){}
    Ray(const Vec3<float>& origin, const Vec3<float>& dir)
    :origin(origin),dir(dir){}
    Ray(const float* _origin, const float* _dir)
    {
        origin = Vec3<float>(_origin);
        dir = Vec3<float>(_dir);
    }

    Vec3<float> origin;
    Vec3<float> dir;
};

Vec3<float> Trace_Test(const float* firstRay_dir, const float* firstRay_origin, const Shader<float>& shader, const IBL& ibl,
                       const embree::EmbreeManager& emb, RandomManager& rnd_manager)
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
            float nee_pdf = ibl.sample(r1, r2, &phi, &theta, Le);

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
                    I = I + throughput * eval_bsdf(local_shadowdir, local_wo, isect, shader) 
                                       * (1.0f/nee_pdf) * costerm * Vec3<float>(Le);
                }
            }
            #endif

            Vec3<float> local_wi = sample(local_wo, isect, rnd_manager, shader);
            Vec3<float> wi = local_wi[0]*e0 + 
                             local_wi[1]*orienting_normal + 
                             local_wi[2]*e2;
            comingRay = Ray(hitPos /*+ 0.001f * orienting_normal*/, wi);

            throughput = throughput * weight(local_wi, local_wo, isect, shader);

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
                ibl.GetLe(Le, theta, phi);
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
                ibl.GetLe(Le, theta, phi);
                I = I + throughput * Vec3<float>(Le);
            #endif
            break;
        }   
    }
    return I;
}


Vec3<float> Trace_debug(const float* firstRay_dir, const float* firstRay_origin, const IBL& ibl,
                        const embree::EmbreeManager& emb, RandomManager& rnd_manager)
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
        return norm_color * norm_color;
    }
    else
    {
        return {0.0f, 0.0f, 0.0f};
    }
}

#endif