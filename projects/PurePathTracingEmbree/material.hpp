#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <variant>
#include "scene.hpp"
#include "bsdf.hpp"
#include "sampling.hpp"
#include <wrap_embree.hpp>
#include "concurrent.hpp"


enum class ShaderType
{
    R, // reflection only
    TR, // Transmit or Reflection
    TRS, // Reflection or Transmit & SSS (Mainly for Uber Shader)
    VOL, // volume only
};

//diffuse only shader(for test)
template<class Real>
class DiffuseShader
{
public:
    DiffuseShader(){};
    //for surface
    Vec3<Real> sample(const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect, RandomManager& rnd_manager,
                      const SceneData<float>& scenedata) const
    {
        Real u = rnd_manager.GetRND();
        Real v = rnd_manager.GetRND();
        return test_randomCosineHemisphere(u, v);
    };
    Real eval_pdf(const Vec3<Real>& wi, const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect,
                  const SceneData<float>& scenedata) const
    {
        return eval_diffuse_pdf(wi, wo);
    };
    Vec3<Real> eval_bsdf(const Vec3<Real>& wi, const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect,
                        const SceneData<float>& scenedata) const
    {
        int geomID = isect.primID;
        Vec3<Real> albedo;
        bool exist_texture = scenedata.mat_infos.albedo_textures[scenedata.mat_infos.mat_indices[geomID]].exist();
        if(exist_texture)
        {
           float u = isect.u;
           float v = isect.v;

           int uv_index0 = scenedata.vertex_infos.uv_indices[3 * geomID + 0]; 
           int uv_index1 = scenedata.vertex_infos.uv_indices[3 * geomID + 1]; 
           int uv_index2 = scenedata.vertex_infos.uv_indices[3 * geomID + 2]; 

           float U0 = scenedata.vertex_infos.uvs[2 * uv_index0 + 0];
           float V0 = scenedata.vertex_infos.uvs[2 * uv_index0 + 1];
           float U1 = scenedata.vertex_infos.uvs[2 * uv_index1 + 0];
           float V1 = scenedata.vertex_infos.uvs[2 * uv_index1 + 1];
           float U2 = scenedata.vertex_infos.uvs[2 * uv_index2 + 0];
           float V2 = scenedata.vertex_infos.uvs[2 * uv_index2 + 1];

           float U = u * (U1 - U0) + v * (U2 - U0);
           float V = u * (V1 - V0) + v * (V2 - V0);

           scenedata.mat_infos.albedo_textures[scenedata.mat_infos.mat_indices[geomID]].getColor(U, V, &(albedo[0]), &(albedo[1]), &(albedo[2]));
        }
        else
        {
            albedo = {scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 0],
                      scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 1],
                      scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 2]};
        }
        return eval_diffuse_bsdf(wi, wo, albedo);
    };
    Vec3<Real> weight(const Vec3<Real>& wi, const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect,
                      const SceneData<float>& scenedata) const
    {
        int geomID = isect.primID;
        Vec3<Real> albedo{scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 0],
                          scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 1],
                          scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 2]};
        return albedo;
    };

    //for subsurface or vpt
    Vec3<Real> sample_phase(const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect){};
    void CalcDensity(const Vec3<Real>& position, Real* sigma_t, Real* sigma_s){}; //this shader is not for vpt
private:
    const ShaderType shader_type = ShaderType::R;
};



template<class T>
using Shader = std::variant<DiffuseShader<T>>;


// wrapper functions 
template<class Real>
Vec3<Real> sample(const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect, RandomManager& rnd_manager
                 ,const Shader<Real>& shader, const SceneData<float>& scenedata)
{
    return std::visit([&](auto& x){return x.sample(wo, isect, rnd_manager, scenedata);}, shader);
}

template<class Real>
Real eval_pdf(const Vec3<Real>& wi, const Vec3<Real>& wo,const embree::IsectInfo<Real>& isect
             ,const Shader<Real>& shader, const SceneData<float>& scenedata)
{
    return std::visit([&](auto& x){return x.eval_pdf(wi, wo, isect, scenedata);},shader);
}

template<class Real>
Vec3<Real> eval_bsdf(const Vec3<Real>& wi, const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect
                    ,const Shader<Real>& shader, const SceneData<float>& scenedata)
{
    return std::visit([&](auto& x){return x.eval_bsdf(wi, wo, isect, scenedata);}, shader);
}

template<class Real>
Vec3<Real> weight(const Vec3<Real>& wi, const Vec3<Real>& wo, const embree::IsectInfo<Real>& isect,
                  const Shader<Real>& shader, const SceneData<float>& scenedata)
{
    return std::visit([&](auto& x){return x.weight(wi, wo, isect, scenedata);}, shader);
}
#endif