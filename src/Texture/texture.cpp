#include "texture.hpp"
#include <cstdio>


void Texture::LoadTexture(const std::string& filename)
{
    std::printf("Load Texture: %s\n", filename.c_str());
    color = stbi_load(filename.c_str(), &width, &height, &channel, 0);
    if(color)
    {
        std::printf("Success!\n");
        is_exist = true;
        std::printf("----> width: %d\n", width);
        std::printf("----> height: %d\n", height);
        std::printf("----> channel: %d\n", channel);
    }
}

Texture::~Texture() 
{
    if(is_exist)
    stbi_image_free(color);
}

void Texture::getColor(const float u, const float v, float* color_x, float* color_y, float* color_z) const
{
    int i = width * u;
    int j = height * v;

    *color_x = color[3 * (width * j + i) + 0] / 255.0;
    *color_y = color[3 * (width * j + i) + 1] / 255.0;
    *color_z = color[3 * (width * j + i) + 2] / 255.0;
}