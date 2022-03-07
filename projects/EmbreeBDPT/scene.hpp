#ifndef SCENE_HPP
#define SCENE_HPP

#include <memory>
#include <iostream>
#include "texture.hpp"
#include "IBL.hpp"
#include <algorithm>
#include <optional>
#include <sstream>

#define TINYOBJLOADER_IMPLEMENTATION 
#include <tiny_obj_loader.h>

// template<class Real>
// class MaterialData
// {
// public:
//     std::unique_ptr<int[]> mat_indices;
//     std::unique_ptr<Real[]> speculers;
//     std::unique_ptr<Real[]> diffuses;
//     std::unique_ptr<Real[]> transmits;
//     std::unique_ptr<Real[]> emissions;
//     std::unique_ptr<Real[]> roughnesss;
//     std::unique_ptr<Real[]> nis;
// };
template<class Real>
class MaterialData
{
public:
    std::vector<Texture> albedo_textures;
    std::vector<int> mat_indices;
    std::vector<Real> speculers;
    std::vector<Real> diffuses;
    std::vector<Real> transmits;
    std::vector<Real> emissions;
    std::vector<Real> roughnesss;
    std::vector<Real> nis;
};

// template<class Real>
// class VertexData
// {
// public:
//     //uvs
//     std::unique_ptr<int[]> uv_indices;
//     std::unique_ptr<Real[]> uvs;
//     //normals
//     std::unique_ptr<int[]> normal_indices;
//     std::unique_ptr<Real[]> normals;
// };

template<class Real>
class VertexData
{
public:
    //uvs
    std::vector<int> uv_indices;
    std::vector<Real> uvs;
    //normals
    std::vector<int> normal_indices;
    std::vector<Real> normals;
};

template<class Real>
class SceneData
{
public:
    IBL ibl;
    VertexData<Real> vertex_infos;
    MaterialData<Real> mat_infos;
};

class OBJloader
{
public:
    OBJloader(const std::string& filename, const std::string& _mtldir):inputfile(filename), mtldir(_mtldir){}
    std::string inputfile;
    std::string mtldir;
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;
};

template<class Integer>
void LoadObj_Single_Object(OBJloader& loader, std::vector<float>& vertices, std::vector<Integer>& indices, 
                          SceneData<float>& scenedata, float scale)
{

    bool ret = tinyobj::LoadObj(&loader.attrib, &loader.shapes, &loader.materials, 
                                &loader.warn, &loader.err, loader.inputfile.c_str(), loader.mtldir.c_str());

    // error handling
    if (!loader.warn.empty()) {
        std::cout << loader.warn << std::endl;
    }

    if (!loader.err.empty()) {
        std::cerr << loader.err << std::endl;
    }

    if (!ret) {
        exit(1);
    }




    // load vertices and vertex info
    vertices.resize(loader.attrib.vertices.size());
    for(size_t i = 0; i < loader.attrib.vertices.size(); ++i)
    {
        vertices[i] = loader.attrib.vertices[i];
        
        // if(i % 3 == 0) vertices[i] *= -1.0f;
    }

    bool exist_vnormal = false;
    if(loader.attrib.normals.size() > 0)
    {
        exist_vnormal = true;
        scenedata.vertex_infos.normals.resize(loader.attrib.normals.size());
        for(int i = 0; i < loader.attrib.normals.size(); ++i)
        {
            scenedata.vertex_infos.normals[i] = loader.attrib.normals[i];
        }
    }

    bool exist_vcoord = false;
    if(loader.attrib.texcoords.size() > 0)
    {
        exist_vcoord = true;
        scenedata.vertex_infos.uvs.resize(loader.attrib.texcoords.size());
        for(int i = 0; i < loader.attrib.texcoords.size(); ++i)
        {
            scenedata.vertex_infos.uvs[i] = loader.attrib.texcoords[i];
        }
    }

    int num_faces = 0;
    for(size_t i  = 0; i < loader.shapes.size(); ++i)
    {
        num_faces += static_cast<int>(loader.shapes[i].mesh.indices.size()/3);
    }

    printf("ObJ INFO: faces = %d\n", num_faces);
    printf("ObJ INFO: shapes = %lu\n", loader.shapes.size());
    printf("ObJ INFO: materials = %lu\n", loader.materials.size());

    //index setting
    scenedata.vertex_infos.normal_indices.resize(num_faces * 3);
    scenedata.vertex_infos.uv_indices.resize(num_faces * 3);
    indices.resize(num_faces * 3);
    scenedata.mat_infos.mat_indices.resize(num_faces);
    size_t offset = 0;
    for(size_t i = 0; i < loader.shapes.size(); ++i)
    {
        for(size_t f = 0; f < loader.shapes[i].mesh.indices.size()/3; ++f)
        {
            if(exist_vnormal)
            {
                scenedata.vertex_infos.normal_indices[3 * (offset + f) + 0] = 
                    loader.shapes[i].mesh.indices[3 * f + 0].normal_index;
                scenedata.vertex_infos.normal_indices[3 * (offset + f) + 1] = 
                    loader.shapes[i].mesh.indices[3 * f + 1].normal_index;
                scenedata.vertex_infos.normal_indices[3 * (offset + f) + 2] = 
                    loader.shapes[i].mesh.indices[3 * f + 2].normal_index;
            }

            if(exist_vcoord)
            {
                scenedata.vertex_infos.uv_indices[3 * (offset + f) + 0] =
                    loader.shapes[i].mesh.indices[3 * f + 0].texcoord_index;
                scenedata.vertex_infos.uv_indices[3 * (offset + f) + 1] =
                    loader.shapes[i].mesh.indices[3 * f + 1].texcoord_index;
                scenedata.vertex_infos.uv_indices[3 * (offset + f) + 2] =
                    loader.shapes[i].mesh.indices[3 * f + 2].texcoord_index;
            }

            indices[3 * (offset + f) + 0] = 
                loader.shapes[i].mesh.indices[3 * f + 0].vertex_index;
            indices[3 * (offset + f) + 1] = 
                loader.shapes[i].mesh.indices[3 * f + 1].vertex_index;
            indices[3 * (offset + f) + 2] = 
                loader.shapes[i].mesh.indices[3 * f + 2].vertex_index;

            scenedata.mat_infos.mat_indices[offset + f] = 
                loader.shapes[i].mesh.material_ids[f];
        }
        offset += loader.shapes[i].mesh.indices.size()/3;
    }




    //material setting
    scenedata.mat_infos.albedo_textures.resize(loader.materials.size());
    scenedata.mat_infos.diffuses.resize(3 * loader.materials.size());
    scenedata.mat_infos.emissions.resize(3 * loader.materials.size());
    scenedata.mat_infos.speculers.resize(3 * loader.materials.size());
    scenedata.mat_infos.transmits.resize(3 * loader.materials.size());
    scenedata.mat_infos.roughnesss.resize(loader.materials.size());
    //mat_infos.roughnesss.reset(new float[loader.materials.size()]);
    scenedata.mat_infos.nis.resize(loader.materials.size());
    for(int i = 0; i < loader.materials.size(); ++i)
    {
        //diffuse
        scenedata.mat_infos.diffuses[3 * i + 0] = loader.materials[i].diffuse[0];
        scenedata.mat_infos.diffuses[3 * i + 1] = loader.materials[i].diffuse[1];
        scenedata.mat_infos.diffuses[3 * i + 2] = loader.materials[i].diffuse[2];

        //emission
        scenedata.mat_infos.emissions[3 * i + 0] = loader.materials[i].emission[0];
        scenedata.mat_infos.emissions[3 * i + 1] = loader.materials[i].emission[1];
        scenedata.mat_infos.emissions[3 * i + 2] = loader.materials[i].emission[2];

        //speculer
        scenedata.mat_infos.speculers[3 * i + 0] = loader.materials[i].specular[0];
        scenedata.mat_infos.speculers[3 * i + 1] = loader.materials[i].specular[1];
        scenedata.mat_infos.speculers[3 * i + 2] = loader.materials[i].specular[2];

        //transmit
        scenedata.mat_infos.transmits[3 * i + 0] = loader.materials[i].transmittance[0];
        scenedata.mat_infos.transmits[3 * i + 1] = loader.materials[i].transmittance[1];
        scenedata.mat_infos.transmits[3 * i + 2] = loader.materials[i].transmittance[2];

        //roughness
        scenedata.mat_infos.roughnesss[i]= loader.materials[i].roughness;
        // user
        auto it = loader.materials[i].unknown_parameter.find("user-defined-param");
        if(it != loader.materials[i].unknown_parameter.end())
        {
            std::stringstream ss(it->second);
            std::string tmp_s;
            std::vector<float> v;
            while(std::getline(ss, tmp_s, ' '))
            {
                v.push_back(std::stof(tmp_s));
            }
        }

        // user(for filename) 
        auto it_user = loader.materials[i].unknown_parameter.find("user-defined-map");
        if(it_user != loader.materials[i].unknown_parameter.end())
        {
            std::string user_filename = it_user->second;
            printf("find: user(%s)\n", user_filename.c_str());
        } 

        // texture
        if(!loader.materials[i].diffuse_texname.empty())
        {
            std::string albedo_filename = loader.materials[i].diffuse_texname;
            printf("find: map_Kd(%s)\n", albedo_filename.c_str());
            scenedata.mat_infos.albedo_textures[i].LoadTexture(loader.mtldir + "/" + albedo_filename);
        }

        //ior
        scenedata.mat_infos.nis[i] = loader.materials[i].ior;
    }
}

void LoadIBL(const std::string& filename, SceneData<float>& scenedata) 
{
    scenedata.ibl.LoadIBL(filename);
}

#endif