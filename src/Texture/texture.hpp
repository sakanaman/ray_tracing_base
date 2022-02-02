#ifndef __TEXTURE_HPP__
#define __TEXTURE_HPP__
#include <cstddef>
#include <stb_image.h>
#include <string>

class Texture
{
public:
    Texture(){}
    ~Texture();
    void LoadTexture(const std::string& filename);
    void getColor(const float u, const float j, float* color_x, float* color_y, float* color_z) const;
    bool exist()const{return is_exist;}
private:
    bool is_exist = false;
    int width;
    int height;
    int channel;
    unsigned char* color;
};

#endif