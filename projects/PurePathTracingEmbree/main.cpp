#define NEE
#include "raytrace.hpp"
#include "stb_image_write.h"
#include <camera.hpp>
#include <wrap_embree.hpp>

float gamma(float x)
{
    return std::pow(x, 1/2.2f);
}

//Utility Functions!!
void SaveImage(const float *rgb, int width, int height)  // To PNG
{
   std::unique_ptr<unsigned char[]> pixel_colors = std::make_unique<unsigned char[]>(3 * width * height);
  for (int y = 0; y < height; y++) 
  {
    for (int x = 0; x < width; x++) 
    {
      int index = y * width + x;
      pixel_colors[index * 3 + 0] = static_cast<unsigned char>(std::max(0.0f, std::min(gamma(rgb[index * 3 + 0]) * 255.0f, 255.0f)));
      pixel_colors[index * 3 + 1] = static_cast<unsigned char>(std::max(0.0f, std::min(gamma(rgb[index * 3 + 1]) * 255.0f, 255.0f)));
      pixel_colors[index * 3 + 2] = static_cast<unsigned char>(std::max(0.0f, std::min(gamma(rgb[index * 3 + 2]) * 255.0f, 255.0f)));
    }
  }
  stbi_write_png("output.png", width, height, 3, &(pixel_colors[0]), width * 3);
}



int main(int argc, char** argv)
{
    //shape
    std::vector<float> vertices;
    std::vector<int> indices;
    SceneData<float> scenedata;
    OBJloader load("../../../objects/sponza_crytek/sponza.obj", "../../../objects/sponza_crytek");
    float scale = 1.0;
    LoadObj_Single_Object<int>(load, vertices, indices, scenedata,scale);

    LoadIBL("../../../map_textures/PaperMill_E_3k.hdr", scenedata);

    printf("here\n");

    //Shader
    Shader<float> shader = DiffuseShader<float>();

    //Embree
    embree::EmbreeManager emb;
    embree::LoadScene2Embree<int>(vertices, indices, &emb);


    //Camera
    // float theta = 83.0f * M_PI/180.0f;
    // float phi =   135.0f * M_PI/180.0f;
    // float r = 12.0f;
    // float x = r * std::sin(theta) * std::cos(phi);
    // float y = r * std::cos(theta) + 2.0f;
    // float z = r * std::sin(theta) * std::sin(phi);

    float x = -400.0;
    float y = 100.0;
    float z = -40.0; 


    //Screen, Image
    int width = 812, height = 512;
    float* RGB = new float[width * height * 3];
    for(int i = 0; i < width*height*3; i++)
    {
        RGB[i] = 0.0f;
    }
    float screen_height = 3.0f;
    float screen_width = screen_height * float(width) / height;
    float pixel_size = screen_height/height;

    float cameraPos[3] = {x, y, z};
    // float cameraForward[3] = {0.5f* screen_width * -x/r, 
    //                           0.5f* screen_width * -y/r, 
    //                           0.5f* screen_width * -z/r};
    float cameraForward[3] = {-1.0, 0.0, 0.0};
    PinholeCamera<float> pincam(cameraPos, cameraForward);

    int samples = std::stof(argv[1]);
    printf("spp: %d\n", samples);

    //Rendering
    std::function<void(const int*, const int*, RandomManager&)> render = 
        [&](const int* upper_left, const int* bottom_right, RandomManager& rnd_manager)
        {
            for(int x = upper_left[0]; x <  bottom_right[0]; ++x)
            {
                for(int y = upper_left[1]; y < bottom_right[1]; ++y)
                {
                    for(int i = 0; i < samples; ++i)
                    {
                        float u = x * pixel_size - 0.5f * pixel_size * width + pixel_size * rnd_manager.GetRND();
                        float v = y * pixel_size - 0.5f * pixel_size * height + pixel_size * rnd_manager.GetRND();
                        float ray_dir[3], ray_origin[3];
                        pincam.CreateFirstRay(u, v, ray_origin, ray_dir);
                        Vec3<float> result = Trace_Test(ray_dir, ray_origin, shader, emb, rnd_manager, scenedata);
                        // Vec3<float> result = Trace_debug(ray_dir, ray_origin, emb, rnd_manager, scenedata);
                        RGB[3*(width * y + x) + 0] += result[0]/samples;
                        RGB[3*(width * y + x) + 1] += result[1]/samples;
                        RGB[3*(width * y + x) + 2] += result[2]/samples;
                    }
                }
            }
        };
    ParallelRender paral(render);
    paral.Execute(width, height, 20);
    SaveImage(RGB, width, height);
}