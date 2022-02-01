#ifndef __TEXTURE_HPP__
#define __TEXTURE_HPP__
#include <cstddef>
#include <stb_image.h>
#include <string>

class Texture
{
public:
    Texture(const std::string& _filename):filename(_filename){}
    void LoadTexture();
    void getColor(const float u, const float j, float* color_x, float* color_y, float* color_z) const;
private:
    std::string filename; 
    int width;
    int height;
    int channel;
    unsigned char* color;
};

#endif