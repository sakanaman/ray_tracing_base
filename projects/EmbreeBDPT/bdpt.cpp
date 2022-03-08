#include "bdpt.hpp"

//============================
//         Utility
//============================

// pdf(p1 --> p2) --> pdf(p2)
// solid angle --> unit area
// NOTE: not consider visibility
float pdfConvert(const float pdf_omega,
                 const Vec3<float>& p1,
                 const Vec3<float>& p2,
                 const Vec3<float>& n2)
{
   float length = (p1 - p2).length();

   float costheta = n2.dot((p1 - p2).normalize());

   return pdf_omega * costheta / (length * length);
}

bool test_visibility(const Vec3<float>& p1, 
                     const Vec3<float>& p2,
                     const embree::EmbreeManager& emb)
{
    embree::IsectInfo<float> isect;
    Vec3<float> origin = p1;
    Vec3<float> dir = (p2 - p1).normalize();
    bool is_hit 
    = embree::intersect<float>(origin, dir, 
                               embree::RAYMIN, embree::RAYMAX,
                               emb, isect);

    float length = (p1 - p2).length();

    return (is_hit && (isect.tfar < length)) || (!is_hit);
}

float calculate_G(const Vec3<float>& p1,
                  const Vec3<float>& p2,
                  const Vec3<float>& n2,
                  const embree::EmbreeManager& emb)
{
    float length = (p1 - p2).length();

    float costheta = n2.dot((p1 - p2).normalize());

    float V = static_cast<float>(test_visibility(p1, p2, emb));

    return V * costheta / (length * length);
}

float Le0_uniform_hemisphere()
{
    return 1.0f / M_PI;
}

Vec3<float> Le1_uniform_hemisphere(const Vec3<float>& intensity = {1.0f, 1.0f, 1.0f})
{
    return intensity * static_cast<float>(M_PI); 
}

float We0_uniform_hemisphere()
{
    return 1.0f / M_PI;
}

Vec3<float> We1_uniform_hemisphere(const Vec3<float>& intensity = {1.0f, 1.0f, 1.0f})
{
    return intensity * static_cast<float>(M_PI); 
}


